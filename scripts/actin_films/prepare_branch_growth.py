"""Fit complete actin endpoints and coordinate two source-guided branch additions."""
import hashlib
import itertools
import json

import numpy as np
from scipy.interpolate import BSpline
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from prepare_actin_paths import smooth
from prepare_contact_repairs import fit_placement,transform
from project_context import ROOT
from rigid_path_detour import refine


class Growth:
    def __init__(self,assets,targets):
        with np.load(assets/'actin_morph.npz') as data:self.actin={key:data[key] for key in data.files}
        self.core={c:np.load(assets/'7TPT'/(c+'.npz'))['atoms'].astype(float) for c in 'HIJK'}
        plans=json.loads((assets/'approach-plans.json').read_text())
        self.offsets=[np.array(plans['branch_add_'+str(i)]['offset'])*.7 for i in [1,2]]
        self.pivot=self.actin['f_vertices'].mean(0);self.targets=targets;self.repairs={}

    def pose(self,fraction,cycle):
        first=.10+cycle*.31;last=first+.24
        d=1-smooth((fraction-first)/(last-first));far=max(0,(d-.5)*2)
        offset=self.offsets[cycle];side=np.cross(offset,[0,0,1]);side*=1.2/np.linalg.norm(side)
        turn=np.eye(4);turn[:3,:3]=Rotation.from_euler('z',18*smooth(far),degrees=True).as_matrix()
        turn[:3,3]=self.pivot-turn[:3,:3]@self.pivot
        matrix=self.targets[cycle]@turn
        matrix[:3,3]+=offset*d+side*np.sin(np.pi*far)**2
        if cycle in self.repairs:
            record=self.repairs[cycle];a=record['first_fraction'];b=record['last_fraction']
            if a<fraction<b:matrix[:3,3]+=BSpline(record['knots'],record['controls_nm'],3)((fraction-a)/(b-a))
        g=1-smooth((fraction-first-.03)/.08)
        return matrix,g

    def state(self,fraction):
        state=dict(self.core)
        for cycle in [0,1]:
            matrix,g=self.pose(fraction,cycle)
            state['actin'+str(cycle+1)]=transform(self.actin['f_atoms']*(1-g)+self.actin['g_atoms']*g,matrix)
        return state


def main():
    assets=ROOT/'assets';folder=assets/'branch_growth';folder.mkdir(exist_ok=True)
    output=assets/'branch-growth.json'
    record={'passed':False,'limits':'Cross-source rigid actin end fits and CA-guided G/F envelopes. Approach, turning, state interpolation and duration are teaching choices, not observed branch growth or molecular dynamics.',
        'sources':{code:{'url':'https://www.rcsb.org/structure/'+code,
            'sha256':hashlib.sha256((ROOT/'sources'/(code+'.cif')).read_bytes()).hexdigest()} for code in ['7TPT','9B27','2BTF']},
        'actin_asset_sha256':hashlib.sha256((assets/'actin_morph.npz').read_bytes()).hexdigest()}
    output.write_text(json.dumps(record,indent=2)+'\n')
    step=np.array(json.loads((assets/'7TPT/manifest.json').read_text())['step'])
    fit=np.array(json.loads((assets/'actin-target-fits.json').read_text())['branch_K'])
    actin=np.load(assets/'9B27/F.npz')['atoms']
    fixed=np.concatenate([np.load(assets/'7TPT'/(c+'.npz'))['atoms'] for c in 'HIJK'])
    targets=[];placements=[]
    for i in [1,2]:
        initial=np.linalg.matrix_power(step,i)@fit
        correction,report=fit_placement(transform(actin,initial),fixed)
        target=correction@initial;report['matrix']=target.tolist()
        targets.append(target);placements.append(report)
        fixed=np.concatenate([fixed,transform(actin,target)])
    record['placements']=placements
    model=Growth(assets,targets);fractions=np.linspace(0,1,501)
    for cycle in [0,1]:
        states=[model.state(f) for f in fractions];name='actin'+str(cycle+1)
        moving=np.array([s[name] for s in states])
        fixed=np.array([np.concatenate([p for k,p in s.items() if k!=name]) for s in states])
        distances=np.array([cKDTree(q).query(p)[0].min() for p,q in zip(moving,fixed)])
        bad=np.flatnonzero(distances<.215)
        if len(bad):
            first=max(0,int(bad.min())-12);last=min(len(fractions)-1,int(bad.max())+12)
            _,detour=refine(moving,[fixed],first,last,maximum=.65,threshold=.22)
            detour['first_fraction']=float(fractions[first]);detour['last_fraction']=float(fractions[last])
            model.repairs[cycle]=detour
        print('BRANCH_GROWTH_PATH',cycle+1,float(distances.min()),flush=True)
    minima={};worst={};poses=[];shapes=[];radius=.5
    for index,fraction in enumerate(np.linspace(0,1,2001)):
        state=model.state(fraction);bounds={c:(p.min(0),p.max(0)) for c,p in state.items()}
        trees={c:cKDTree(p) for c,p in state.items()}
        for a,b in itertools.combinations(state,2):
            if index and a in 'HIJK' and b in 'HIJK':continue
            low,high=bounds[b];selected=((state[a]>=low-radius)&(state[a]<=high+radius)).all(1)
            value=min(radius,float(trees[b].query(state[a][selected],distance_upper_bound=radius)[0].min())) if selected.any() else radius
            pair=a+'/'+b
            if value<minima.get(pair,float('inf')):minima[pair]=value;worst[pair]=float(fraction)
        pair=[model.pose(fraction,i) for i in [0,1]]
        poses.append([p[0] for p in pair]);shapes.append([p[1] for p in pair])
        if index%500==0:print('BRANCH_GROWTH_SCREEN',index,'of',2000,flush=True)
    record.update({'repairs':model.repairs,'offsets':[v.tolist() for v in model.offsets],
        'screen':{'samples':2001,'threshold_nm':.20,'search_radius_nm':radius,
            'pair_minima_nm_capped_at_radius':minima,'worst_fractions':worst},
        'passed':min(minima.values())>=.20})
    output.write_text(json.dumps(record,indent=2)+'\n')
    if not record['passed']:raise ValueError('Coupled branch growth still overlaps; inspect branch-growth.json')
    np.savez_compressed(folder/'motion.npz',matrices=np.array(poses),g_fractions=np.array(shapes))
    record['motion_sha256']=hashlib.sha256((folder/'motion.npz').read_bytes()).hexdigest()
    output.write_text(json.dumps(record,indent=2)+'\n')
    print('BRANCH_GROWTH_READY',min(minima.values()),flush=True)


if __name__=='__main__':main()
