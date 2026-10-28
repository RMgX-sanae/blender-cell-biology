"""Couple SOS association, Ras opening and intact nucleotide exchange.

All sampled contacts are screened. This does not infer rates or energetics.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from interaction_motion import matrices, source_contacts
from prepare_phosphotransfer_motion import smooth, covalent_interpolate
from rigid_path_detour import refine


def apply(points, poses):
    return np.einsum('fij,nj->fni',poses[:,:3,:3],points)+poses[:,None,:3,3]


def pair_check(first,second,allowed=None,threshold=.20):
    counts=[];minimum=[]
    for moving,fixed in zip(first,second):
        tree=cKDTree(fixed);distance,_=tree.query(moving)
        bad=set()
        for i in np.flatnonzero(distance<threshold):
            for j in tree.query_ball_point(moving[i],threshold):
                d=float(np.linalg.norm(moving[i]-fixed[j]))
                if allowed is None or (i,j) not in allowed or d<allowed[i,j]-.01:bad.add(int(i))
        counts.append(len(bad));minimum.append(float(distance.min()))
    return np.array(counts),np.array(minimum)


def stage(q,ras0,ras1,sos,nucleotide,spec):
    sos_pose=matrices(sos,spec['sos_direction'],18,q*.5,'loop',spec['sos_turn'],spec['sos_bend'])
    opening=smooth((q-spec['open_start'])/(spec['open_end']-spec['open_start']))
    departure=smooth((q-spec['release_start'])/(spec['release_end']-spec['release_start']))
    nucleotide_pose=matrices(nucleotide,spec['nucleotide_direction'],18,.5+departure*.5,'loop',0,0)
    return {'ras':ras0[None]*(1-opening[:,None,None])+ras1[None]*opening[:,None,None],
            'sos':apply(sos,sos_pose),'nucleotide':apply(nucleotide,nucleotide_pose),
            'sos_pose':sos_pose,'nucleotide_pose':nucleotide_pose,'opening':opening}


def search(ras0,ras1,sos,nucleotide,sos_path,frames):
    q=np.linspace(0,1,49)
    known=source_contacts(nucleotide,ras0,.20)
    sos_known=source_contacts(sos,ras1,.20)
    i=np.arange(384);z=1-2*(i+.5)/len(i);angle=i*2.399963229728653
    directions=np.column_stack([np.sqrt(1-z*z)*np.cos(angle),np.sqrt(1-z*z)*np.sin(angle),z])
    desired=nucleotide.mean(0)-ras0.mean(0);desired/=np.linalg.norm(desired)
    directions=np.vstack([desired,directions[np.argsort(-(directions@desired))]])
    best=None;tested=0
    for begin,end in [(.70,.94),(.75,.97),(.65,.92),(.80,.99),(.60,.88)]:
        for lead in [0,-.035,.035,-.07,.07]:
            spec={'sos_direction':sos_path['direction'],'sos_turn':sos_path['turn_radians'],
                'sos_bend':sos_path['bend_fraction'],'open_start':begin,'open_end':end,
                'release_start':begin+lead,'release_end':min(1,end+lead),'nucleotide_direction':directions[0].tolist()}
            fixed=stage(q,ras0,ras1,sos,nucleotide,spec)
            protein_counts,_=pair_check(fixed['sos'],fixed['ras'],sos_known)
            if protein_counts.max()>0:continue
            trees_r=[cKDTree(x) for x in fixed['ras']];trees_s=[cKDTree(x) for x in fixed['sos']]
            ranking=[]
            departure=smooth((q-spec['release_start'])/(spec['release_end']-spec['release_start']))
            for direction in directions:
                poses=matrices(nucleotide,direction,18,.5+departure*.5,'loop',0,0)
                points_by_frame=apply(nucleotide,poses)
                bad=0;worst=0
                for index,points in enumerate(points_by_frame):
                    for tree in [trees_r[index],trees_s[index]]:
                        d,_=tree.query(points);n=int(np.count_nonzero(d<.20));bad+=n;worst=max(worst,n)
                ranking.append(((worst,bad),direction))
            for _,direction in sorted(ranking,key=lambda item:item[0])[:12]:
                spec['nucleotide_direction']=direction.tolist();sample=stage(np.linspace(0,1,frames+1),ras0,ras1,sos,nucleotide,spec)
                checks={}
                for label,x,y,allow in [('sos_ras',sample['sos'],sample['ras'],sos_known),
                                      ('nucleotide_ras',sample['nucleotide'],sample['ras'],known),
                                      ('nucleotide_sos',sample['nucleotide'],sample['sos'],None)]:
                    counts,d=pair_check(x,y,allow)
                    checks[label]={'maximum_overlap_atoms':int(counts.max()),'flagged_frames':np.flatnonzero(counts>0).tolist(),'minimum_distance_nm':float(d.min())}
                score=(max(c['maximum_overlap_atoms'] for c in checks.values()),sum(len(c['flagged_frames']) for c in checks.values()))
                tested+=1
                if best is None or score<best[0]:
                    best=(score,dict(spec),checks)
                    print('candidate',score,spec,flush=True)
                if score==(0,0):return best[1],best[2],tested
    if best is None:raise ValueError('No protein approach survived the coupled Ras-state screen')
    return best[1],best[2],tested


def refine_entry(points,ras,sos,spec):
    """Add a smooth rigid detour around the constrained, moving pocket rim."""
    counts,_=pair_check(points,ras);bad=np.flatnonzero(counts)
    if not len(bad):return points,spec
    margin=max(8,round(len(points)*.025));lo=max(1,int(bad[0])-margin);hi=min(len(points)-2,int(bad[-1])+margin)
    active=np.arange(lo,hi+1);weight=np.sin(np.pi*(active-lo)/(hi-lo))**2
    ras_trees=[cKDTree(ras[i]) for i in active];sos_trees=[cKDTree(sos[i]) for i in active]
    middle=int(bad[len(bad)//2]);dist,indices=cKDTree(ras[middle]).query(points[middle]);nearest=int(np.argmin(dist))
    desired=points[middle,nearest]-ras[middle,indices[nearest]];desired/=np.linalg.norm(desired)
    i=np.arange(192);z=1-2*(i+.5)/len(i);angle=i*2.399963229728653
    directions=np.column_stack([np.sqrt(1-z*z)*np.cos(angle),np.sqrt(1-z*z)*np.sin(angle),z])
    directions=np.vstack([desired,directions[np.argsort(-(directions@desired))]])
    for distance in [.03,.05,.08,.12,.18,.25,.35]:
        for direction in directions:
            local=points[active]+weight[:,None,None]*direction[None,None,:]*distance
            passed=True
            for coordinates,rt,st in zip(local,ras_trees,sos_trees):
                if rt.query(coordinates)[0].min()<.20 or st.query(coordinates)[0].min()<.20:
                    passed=False;break
            if passed:
                result=points.copy();result[active]=local
                record=dict(spec,local_detour={'first_index':lo,'last_index':hi,'direction':direction.tolist(),
                    'maximum_translation_nm':distance,'weight':'sin(pi*u)^2; zero displacement and slope at both ends'})
                print('constrained entry detour',record['local_detour'],flush=True)
                return result,record
    margin=max(20,round(len(points)*.035))
    result,detour=refine(points,[ras,sos],max(1,int(bad[0])-margin),
        min(len(points)-2,int(bad[-1])+margin))
    if detour['passed']:
        return result,dict(spec,local_detour=detour)
    raise ValueError('No screened local detour; new source-guided pocket motion is required')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--frames',type=int,default=576);p.add_argument('--fps',type=int,default=24)
    p.add_argument('--mode',choices=['teaching','continuous'],default='teaching')
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory')
    if a.frames<240 or a.fps<1:p.error('Use at least 240 frames and a positive fps')
    m=json.loads((a.asset_dir/'manifest.json').read_text())
    d={k:dict(np.load(a.asset_dir/(k+'.npz'))) for k in ['ras','sos','old_gdp','new_gtp']}
    results={};plans={}
    for name,ras_field in [('old_gdp','atoms'),('new_gtp','gtp')]:
        nucleotide=np.vstack([d[name]['atoms'],d[name]['magnesium']])
        plan,checks,tested=search(d['ras'][ras_field],d['ras']['gef'],d['sos']['atoms'],nucleotide,m['paths']['sos'],a.frames//2)
        results[name]={'checks':checks,'candidates_checked':tested};plans[name]=plan
    a.output.mkdir(parents=True)
    passed=all(c['maximum_overlap_atoms']==0 for r in results.values() for c in r['checks'].values())
    report={'kind':'ras-sos-nucleotide-exchange','asset_dir':str(a.asset_dir.resolve()),'frames':a.frames,'fps':a.fps,'mode':a.mode,
        'plans':plans,'path_checks':results,'candidate_path_checks_passed':passed,
        'passed':False,'status':'checking_constrained_motion','threshold_nm':.20}
    (a.output/'manifest.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
    if not passed:raise SystemExit('Coupled motion requires correction; diagnostic manifest retained')
    t=np.linspace(0,1,a.frames+1)
    intervals=[(.06,.43),(.54,.91)] if a.mode=='teaching' else [(0,.5),(.5,1)]
    first=stage(np.clip((t-intervals[0][0])/(intervals[0][1]-intervals[0][0]),0,1),d['ras']['atoms'],d['ras']['gef'],d['sos']['atoms'],
                np.vstack([d['old_gdp']['atoms'],d['old_gdp']['magnesium']]),plans['old_gdp'])
    second=stage(1-np.clip((t-intervals[1][0])/(intervals[1][1]-intervals[1][0]),0,1),d['ras']['gtp'],d['ras']['gef'],d['sos']['atoms'],
                 np.vstack([d['new_gtp']['atoms'],d['new_gtp']['magnesium']]),plans['new_gtp'])
    choose=t>=.5
    ras=first['ras'].copy();ras[choose]=second['ras'][choose]
    for index in range(len(t)):
        part=second if choose[index] else first
        opening=float(part['opening'][index])
        if 0<opening<1:
            start=d['ras']['gtp' if choose[index] else 'atoms']
            ras[index]=covalent_interpolate(start,d['ras']['gef'],d['ras']['bonds'],opening,
                initial=ras[index-1],max_nfev=600,ftol=1e-5)
            if index%12==0:print('constrained Ras frame',index+1,flush=True)
    sos=first['sos_pose'].copy();sos[choose]=second['sos_pose'][choose]
    sos_atoms=apply(d['sos']['atoms'],sos)
    if pair_check(second['nucleotide'],ras,source_contacts(second['nucleotide'][-1],d['ras']['gtp'],.20))[0].max()>0:
        try:
            second['nucleotide'],plans['new_gtp']=refine_entry(second['nucleotide'],ras,sos_atoms,plans['new_gtp'])
        except ValueError as error:
            report.update(status='failed_constrained_entry',failure=str(error))
            (a.output/'manifest.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
            np.savez_compressed(a.output/'unrefined-motion.npz',ras=ras,sos_pose=sos,
                old_nucleotide=first['nucleotide'],new_nucleotide=second['nucleotide'])
            raise
        report['plans']=plans;report['entry_refined_after_protein_constraints']=True
    checks={}
    for label,x,y,allowed in [
        ('sos_ras',sos_atoms,ras,source_contacts(d['sos']['atoms'],d['ras']['gef'],.20)),
        ('old_gdp_ras',first['nucleotide'],ras,source_contacts(first['nucleotide'][0],d['ras']['atoms'],.20)),
        ('new_gtp_ras',second['nucleotide'],ras,source_contacts(second['nucleotide'][-1],d['ras']['gtp'],.20)),
        ('old_gdp_sos',first['nucleotide'],sos_atoms,None),('new_gtp_sos',second['nucleotide'],sos_atoms,None),
        ('old_gdp_new_gtp',first['nucleotide'],second['nucleotide'],None)]:
        counts,distances=pair_check(x,y,allowed)
        checks[label]={'maximum_overlap_atoms':int(counts.max()),'flagged_frames':(np.flatnonzero(counts>0)+1).tolist(),
            'minimum_distance_nm':float(distances.min())}
    bonds=d['ras']['bonds'];bond_lengths=np.linalg.norm(ras[:,bonds[:,0]]-ras[:,bonds[:,1]],axis=2)
    report.update(full_timeline_checks=checks,intervals=intervals,
        maximum_ras_connection_nm=float(bond_lengths.max()),minimum_ras_connection_nm=float(bond_lengths.min()),
        heavy_atoms_in_nucleotides_and_metals=len(d['old_gdp']['atoms'])+len(d['new_gtp']['atoms'])+2,
        chemistry='Intact old GDP leaves; independent new GTP arrives. Neither molecule changes covalent bonds.',
        magnesium_motion='Each displayed Mg ion follows its nucleotide for visual tracking; this is not a hydration/exchange simulation.',
        conformation='Source-guided interpolation with bond and neighboring-angle distance constraints drives an envelope; not a validated atomic dynamics trajectory.',
        passed=all(c['maximum_overlap_atoms']==0 for c in checks.values()) and float(bond_lengths.max())<.22 and float(bond_lengths.min())>.10)
    report['status']='validated_geometry' if report['passed'] else 'failed_full_timeline'
    (a.output/'manifest.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
    np.savez_compressed(a.output/'exchange.npz',ras=ras,sos_pose=sos,old_gdp=first['nucleotide'][:,:-1],
        old_mg=first['nucleotide'][:,-1:],new_gtp=second['nucleotide'][:,:-1],new_mg=second['nucleotide'][:,-1:])
    if not report['passed']:raise SystemExit('Full timeline failed; diagnostic motion retained')
    print(json.dumps({'passed':passed,'output':str(a.output)},ensure_ascii=False))


if __name__=='__main__':main()
