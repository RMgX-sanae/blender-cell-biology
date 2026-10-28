"""Diagnose severe inter-protein atomic-center overlap in evaluated frames.

This is a geometric diagnostic, not a molecular dynamics/energy calculation.
The 0.12 nm threshold intentionally flags deep overlaps, not ordinary contact.
"""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

from project_context import ROOT, SKILL, FFMPEG, FONT, LATIN_FONT, OPTIONS

def main():
    p=argparse.ArgumentParser();p.add_argument('--shots',default='a03,a04,a05,a07,c01,c04,c05,f03,f04,f06')
    p.add_argument('--threshold',type=float,default=.12,help='Atomic-center distance in nm; use .20 for the stricter contact screen')
    p.add_argument('--fail-on-overlap',action='store_true')
    p.add_argument('--summary',default='contacts-summary.json',choices=('contacts-summary.json','strict-contacts-summary.json'))
    a=p.parse_args();reports=[]
    if not np.isfinite(a.threshold) or a.threshold<=0:p.error('--threshold must be positive and finite')
    search_radius=max(.5,a.threshold)
    for id in a.shots.split(','):
        folder=ROOT/'shots'/id
        report=json.loads((folder/'evaluated-audit.json').read_text())
        d=np.load(folder/'evaluated-transforms.npz');ms=d['matrices'];sh=d['shapes'];vis=d['visible']
        gate=d['gates'] if 'gates' in d else np.zeros_like(sh)
        assets=[np.load(o['asset']) for o in report['molecules']]
        base=[x['atoms'] if 'atoms' in x else x['f_atoms'] for x in assets]
        trees={};results=[];checks=0
        for i,j in itertools.combinations(range(len(assets)),2):
            # Relative transforms remove unrelated camera composition motion.
            relative=np.linalg.inv(ms[:,j])@ms[:,i]
            sample={}
            for f in range(len(ms)):
                if not vis[f,i] or not vis[f,j]:continue
                key=(np.round(relative[f],6).tobytes(),round(float(sh[f,i]),6),round(float(sh[f,j]),6),round(float(gate[f,i]),6),round(float(gate[f,j]),6))
                if key in sample:continue
                sample[key]=f
            deepest=100.;bad=[];maximum=0
            for _,f in sample.items():
                aa=base[i];bb=base[j]
                if 'g_atoms' in assets[i]:aa=aa+(assets[i]['g_atoms']-aa)*sh[f,i]
                if 'g_atoms' in assets[j]:bb=bb+(assets[j]['g_atoms']-bb)*sh[f,j]
                if 'gate_atoms_delta' in assets[i]:aa=aa+assets[i]['gate_atoms_delta']*gate[f,i]
                if 'gate_atoms_delta' in assets[j]:bb=bb+assets[j]['gate_atoms_delta']*gate[f,j]
                aa=aa@relative[f,:3,:3].T+relative[f,:3,3]
                if np.any(aa.min(0)-bb.max(0)>search_radius) or np.any(bb.min(0)-aa.max(0)>search_radius):continue
                tk=(j,round(float(sh[f,j]),6),round(float(gate[f,j]),6))
                if tk not in trees:trees[tk]=cKDTree(bb)
                dist,_=trees[tk].query(aa,k=1,distance_upper_bound=search_radius)
                minimum=float(dist.min());count=int((dist<a.threshold).sum());checks+=1
                deepest=min(deepest,minimum)
                if count:bad.append(f+1);maximum=max(maximum,count)
            if bad:
                results.append({'a':report['molecules'][i]['name'],'b':report['molecules'][j]['name'],
                                'sampled_overlap_states':len(bad),'first_frame':min(bad),'last_frame':max(bad),
                                'min_center_distance_nm':deepest,'max_atoms_below_threshold':maximum})
        result={'shot':id,'frames':len(ms),'unique_pair_checks':checks,'threshold_nm':a.threshold,'flagged_pairs':results,
                'limits':'Atomic centers only; missing atoms, ligands, environmental proteins excluded. G/F atom interpolation is a source-guided approximation.'}
        (folder/'contacts.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
        reports.append(result);print(id,len(results),'flagged pairs',flush=True)
    (ROOT/'validation'/a.summary).write_text(json.dumps(reports,ensure_ascii=False,indent=2)+'\n')
    if a.fail_on_overlap and any(r['flagged_pairs'] for r in reports):
        raise SystemExit('Contact screen failed; see '+a.summary)

if __name__=='__main__':main()
