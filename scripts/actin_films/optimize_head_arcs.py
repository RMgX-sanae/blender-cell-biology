"""Select illustrative FH2 head arcs against the saved evaluated geometry."""
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from project_context import ROOT, SKILL, FFMPEG, FONT, LATIN_FONT, OPTIONS;folder=ROOT/'shots/f04'
report=json.loads((folder/'evaluated-audit.json').read_text());data=np.load(folder/'evaluated-transforms.npz')
mats=data['matrices'];shape=data['shapes'];gate=data['gates'];assets=[np.load(x['asset']) for x in report['molecules']]
names=[x['name'] for x in report['molecules']];heads=[names.index('9B27 G'),names.index('9B27 H')]
theta=np.deg2rad(-15);rot=np.array([[np.cos(theta),0,np.sin(theta)],[0,1,0],[-np.sin(theta),0,np.cos(theta)]])
def points(i,f):
    a=assets[i];p=a['atoms'] if 'atoms' in a else a['f_atoms']
    if 'g_atoms' in a:p=p+(a['g_atoms']-p)*shape[f,i]
    if 'gate_atoms_delta' in a:p=p+a['gate_atoms_delta']*gate[f,i]
    return p@mats[f,i,:3,:3].T+mats[f,i,:3,3]
results=[]
old_arcs=json.loads((ROOT/'assets/head-arcs.json').read_text()) if (ROOT/'assets/head-arcs.json').exists() else [{'arc':[0,0,-8]},{'arc':[0,0,8]}]
for cycle in [0,1]:
    i=heads[cycle];begin=.08+cycle*.38;end=begin+.14;frames=np.arange(round(begin*(len(mats)-1))+1,round(end*(len(mats)-1)),2)
    old=np.array(old_arcs[cycle]['arc']);states=[]
    for f in frames:
        phase=((f/(len(mats)-1))-begin)/.14;wave=np.sin(np.pi*phase)**2
        obstacle=np.concatenate([points(j,f) for j,n in enumerate(names) if j!=i and (n.startswith('9B27 ') or n.startswith('actin_'))])
        states.append((points(i,f),cKDTree(obstacle),wave))
    candidates=[old,np.array([0,10.,0]),np.array([0,-10.,0]),np.array([0,14.,0]),np.array([0,-14.,0])]
    for k in range(48):
        z=1-2*(k+.5)/48;p=k*2.39996323;v=np.array([np.sqrt(1-z*z)*np.cos(p),np.sqrt(1-z*z)*np.sin(p),z])
        candidates.extend([v*8,v*12])
    best=None
    for v in candidates:
        maximum=0;total=0
        for pa,tree,wave in states:
            query=pa+(rot@(v-old))*wave
            d,_=tree.query(query,k=1,distance_upper_bound=.12);n=int(np.isfinite(d).sum());maximum=max(maximum,n);total+=n
        score=(maximum,total,float(np.linalg.norm(v)))
        if best is None or score<best[0]:
            best=(score,v);print(cycle,score,v.tolist(),flush=True)
        if maximum==0:break
    results.append({'arc':best[1].tolist(),'max_atoms_below_0_12nm':best[0][0],'total_diagnostic_hits':best[0][1]})
(ROOT/'assets/head-arcs.json').write_text(json.dumps(results,indent=2))
