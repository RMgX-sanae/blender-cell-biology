"""Source-guided FH2 conformer exchange for a labelled teaching stepping model."""
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from project_context import ROOT, SKILL, FFMPEG, FONT, LATIN_FONT, OPTIONS
fits=json.loads((ROOT/'assets/formin-head-fit.json').read_text())
for name,other in [('G','H'),('H','G')]:
    a=np.load(ROOT/'assets/9B27'/f'{name}.npz');b=np.load(ROOT/'assets/9B27'/f'{other}.npz')
    m=np.linalg.inv(np.array(fits[f'{name}_to_{other}']))
    ca=dict(zip(a['ca_residues'],a['ca']));cb=dict(zip(b['ca_residues'],b['ca']))
    common=sorted(ca.keys()&cb.keys());x=np.array([ca[n] for n in common])
    y=np.array([cb[n] for n in common])@m[:3,:3].T+m[:3,3]
    delta=y-x
    def warp(v):
        dist,index=cKDTree(x).query(v,k=12)
        weights=1/(dist**2+.02);weights/=weights.sum(1)[:,None]
        return v+np.sum(delta[index]*weights[:,:,None],axis=1)
    verts=warp(a['vertices']);atoms=warp(a['atoms'])
    byatom={(int(n),str(k)):p@m[:3,:3].T+m[:3,3] for n,k,p in zip(b['residue_numbers'],b['atom_names'],b['atoms'])}
    matched=0
    for i,(n,k) in enumerate(zip(a['residue_numbers'],a['atom_names'])):
        if (int(n),str(k)) in byatom:atoms[i]=byatom[(int(n),str(k))];matched+=1
    # An explicitly schematic elastic opening keeps the partner-contact region
    # anchored while the loading side opens. It is not an inferred atomistic
    # intermediate. The intact chain envelope and its topology are retained.
    h=np.load(ROOT/'assets/9B27/H.npz');partner=np.load(ROOT/'assets/9B27/G.npz')['atoms']
    tree=cKDTree(partner);gate_vector=np.array([0.,1.7,1.5])
    def gate_delta(points):
        distance,_=tree.query(points);t=np.clip((distance-.5)/1.4,0,1);weight=t*t*(3-2*t)
        return weight[:,None]*gate_vector
    if name=='G':
        forward=np.array(fits['G_to_H'])
        gv=verts@forward[:3,:3].T+forward[:3,3]
        ga=atoms@forward[:3,:3].T+forward[:3,3]
        gatev=gate_delta(gv)@m[:3,:3].T;gatea=gate_delta(ga)@m[:3,:3].T
    else:
        gatev=gate_delta(a['vertices']);gatea=gate_delta(a['atoms'])
    np.savez_compressed(ROOT/'assets'/f'head_{name}_morph.npz',f_vertices=a['vertices'],g_vertices=verts,
                        faces=a['faces'],f_atoms=a['atoms'],g_atoms=atoms,
                        gate_vertices_delta=gatev,gate_atoms_delta=gatea)
    print(name,matched,'matched atoms of',len(atoms))
