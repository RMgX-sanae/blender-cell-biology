"""Choose unobstructed illustrative approach directions from source atom sets.

This collision avoidance defines a teaching path, not a reaction pathway.
"""
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from project_context import ROOT, SKILL, FFMPEG, FONT, LATIN_FONT, OPTIONS

def atoms(code,chains):return np.concatenate([np.load(ROOT/'assets'/code/f'{c}.npz')['atoms'] for c in chains])
def trans(a,m):return a@m[:3,:3].T+m[:3,3]
def choose(a,b,desired,length=10):
    tree=cKDTree(b);desired=np.array(desired,dtype=float);desired/=np.linalg.norm(desired)
    j=np.arange(300);z=1-2*(j+.5)/len(j);phi=j*2.399963229728653
    dirs=np.column_stack([np.sqrt(1-z*z)*np.cos(phi),np.sqrt(1-z*z)*np.sin(phi),z])
    dirs=np.vstack([desired,dirs]);order=np.argsort(-(dirs@desired));best=None
    distances=np.r_[np.linspace(.04,.8,9),np.linspace(1,4,9),np.linspace(5,length,6)]
    for k in order:
        v=dirs[k]
        points=(a[::2,None,:]+distances[None,:,None]*v).reshape(-1,3)
        d,_=tree.query(points,k=1,distance_upper_bound=.12)
        count=int(np.isfinite(d).sum());score=(count,1-float(v@desired))
        if best is None or score<best[0]:best=(score,v)
        if count==0:
            d,_=tree.query((a[:,None,:]+distances[None,:,None]*v).reshape(-1,3),k=1,distance_upper_bound=.12)
            count=int(np.isfinite(d).sum())
            if count==0:best=((0,1-float(v@desired)),v);break
    result={'offset':(best[1]*length).tolist(),'sampled_deep_contacts':best[0][0]}
    return result

def main():
    plans={}
    def add(key,a,b,desired,length=10):
        plans[key]=choose(a,b,desired,length);print(key,plans[key],flush=True)
    wrc=atoms('7USE','ABCDE')
    for ch in 'FG':
        moving=atoms('7USE',ch)
        obstacles=np.concatenate([wrc,atoms('7USE','G' if ch=='F' else 'F')])
        add('rac_'+ch,moving,obstacles,moving.mean(0)-wrc.mean(0),10)
    did=atoms('2BAP','AB')
    for ch in 'CD':
        moving=atoms('2BAP',ch)
        obstacles=np.concatenate([did,atoms('2BAP','D' if ch=='C' else 'C')])
        add('dad_'+ch,moving,obstacles,moving.mean(0)-did.mean(0),6)
    arp=atoms('7TPT','ABCDEFG');mother=atoms('7TPT','NOPQRS')
    add('arp',arp,mother,[0,1,.6],11)
    prev=np.concatenate([arp,mother])
    for ch in 'HIJK':
        a=atoms('7TPT',ch);add('branch_'+ch,a,prev,[.3,1,.1],14)
        prev=np.concatenate([prev,a])
    step=np.array(json.loads((ROOT/'assets/7TPT/manifest.json').read_text())['step'])
    fit=np.array(json.loads((ROOT/'assets/actin-target-fits.json').read_text())['branch_K'])
    fa=atoms('9B27','F');prev=np.concatenate([mother,arp,atoms('7TPT','HIJK')])
    for i in [1,2]:
        a=trans(fa,np.linalg.matrix_power(step,i)@fit);add('branch_add_'+str(i),a,prev,[.3,1,-.1],14)
        prev=np.concatenate([prev,a])
    add('debranch',atoms('7TPT','HIJK'),np.concatenate([arp,mother]),[.4,1,.2],12)
    core=atoms('9B27','ABCDEF');heads=atoms('9B27','GH')
    add('formin_release',heads,core,[1,0,0],15)
    prev=heads
    for ch in 'DEF':
        a=atoms('9B27',ch);add('seed_'+ch,a,prev,[0,1,0],14)
        prev=np.concatenate([prev,a])
    step=np.array(json.loads((ROOT/'assets/9B27/manifest.json').read_text())['step'])
    prev=core
    for i in [1,2]:
        m=np.linalg.matrix_power(step,i);a=trans(fa,m)
        add('formin_add_'+str(i),a,np.concatenate([prev,trans(heads,m)]),[.3,1,0],14)
        prev=np.concatenate([prev,a])
    capfit=np.array(json.loads((ROOT/'assets/cap-fit.json').read_text())['matrix'])
    old=ROOT/'assets/filaments/cap'
    cap=np.concatenate([np.load(old/f'chain_{c}.npz')['atoms'] for c in 'GH'])
    add('cap',trans(cap,capfit),core,[1,0,0],12)
    actins=atoms('6VAO','ABCDE')
    for ch in 'FGHIJ':
        a=atoms('6VAO',ch);center=a.mean(0)
        add('cofilin_'+ch,a,actins,[0,center[1],center[2]+1.4],4)
    left=atoms('6VAO','CAFJ');right=atoms('6VAO','DEBIGH')
    add('sever',right,left,[1,0,0],9)
    left_relative=left-np.array(plans['sever']['offset'])
    add('monomer_release',atoms('6VAO','DI'),np.concatenate([atoms('6VAO','EBGH'),left_relative]),[-1,-.5,0],9)
    add('cap_recycle',atoms('6FM2-assembly1',['A']),atoms('6FM2-assembly1',['B','A-2','B-2']),[-1,-1,-.5],10)
    (ROOT/'assets/approach-plans.json').write_text(json.dumps(plans,indent=2)+'\n')

if __name__=='__main__':main()
