"""Prepare the experimental Rab5–EEA1 interface and endosomal SNARE core in local rigid frames."""
import argparse
import hashlib
import json
from pathlib import Path

import gemmi
import numpy as np
from prepare_molecular_asset import surface


def unit(v):
    return v/np.linalg.norm(v)


def atom_xyz(residues,name=None):
    return np.asarray([list(a.pos) for res in residues for a in res
                       if a.occ>0 and not a.element.is_hydrogen and a.altloc in ('\x00','A')
                       and (name is None or a.name==name)])


def basis(axis,second):
    x=unit(axis);y=unit(second-np.dot(second,x)*x)
    return np.column_stack((x,y,np.cross(x,y)))


def prepare(sources,out):
    scope=json.loads((Path(__file__).resolve().parents[1]/'assets/vesicle-uptake-scope.json').read_text())
    arrays={};evidence={}
    for pdb,config in scope['coordinate_sources'].items():
        path=sources/(pdb+'.cif')
        digest=hashlib.sha256(path.read_bytes()).hexdigest()
        if digest!=config['sha256']:raise ValueError('Requalify changed source '+pdb)
        model=gemmi.read_structure(str(path))[0]
        chains={c:model[c].get_polymer() for c in config['chains']}
        ca={c:atom_xyz(poly,'CA') for c,poly in chains.items()}
        if pdb=='3MJH':
            origin=ca['A'][-1]
            axes=basis(ca['B'].mean(0)-origin,ca['A'][0]-origin)
            for c,poly in chains.items():
                vertices,faces=surface(atom_xyz(poly),1.4,1.8,.8)
                arrays[c+'_verts']=(vertices-origin)@axes/10
                arrays[c+'_faces']=faces
            arrays['finger_tail']=(ca['B'][-1]-origin)@axes/10
            for c,code in [('A','GTP'),('A','MG'),('B','ZN')]:
                residues=[r for r in model[c] if r.name==code]
                if len(residues)!=1:raise ValueError('Missing or ambiguous bound component '+code)
                arrays[code]=(atom_xyz(residues)-origin)@axes/10
        else:
            ends=np.array([p[-4:].mean(0) for p in ca.values()])
            starts=np.array([p[:4].mean(0) for p in ca.values()])
            origin=starts.mean(0);axes=basis(ends.mean(0)-origin,starts[0]-origin)
            for c,p in ca.items():
                # The four experimental C-alpha traces keep their distinct lengths.
                arrays['helix_'+c]=(p-origin)@axes/10
                if np.max(np.linalg.norm(np.diff(p,axis=0),axis=1))>5:
                    raise ValueError('Interrupted SNARE trace '+c)
        evidence[pdb]={'sha256':digest,'chains':config['chains'],'source':config['url'],
                       'display_unit':'nm','origin_angstrom':origin.tolist(),'basis_columns':axes.tolist()}
    out.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(out/'structures.npz',**arrays)
    (out/'structures.json').write_text(json.dumps({'sources':evidence,
        'geometry_sha256':hashlib.sha256((out/'structures.npz').read_bytes()).hexdigest()},indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();prepare(a.sources,a.output)
