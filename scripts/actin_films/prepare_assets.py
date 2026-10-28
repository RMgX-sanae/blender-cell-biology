"""Prepare the actual public structures used in the two complete teaching films."""
import hashlib
import json
from pathlib import Path
import sys

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from project_context import SKILL
sys.path.insert(0,str(SKILL/'scripts'))
from prepare_molecular_asset import atom_table, surface, xyz
from prepare_filament_assets import rigid_fit

from project_context import ROOT, SKILL, FFMPEG, FONT, LATIN_FONT, OPTIONS
CHOICES={'7USC':list('ABCDE'),'7USE':list('ABCDEFG'),'1K8K':list('ABCDEFG'),
         '7TPT':list('ABCDEFGHIJKLMNOPQRSTU'),'9B27':list('ABCDEFGH'),
         '2BTF':['A','P'],'6VAO':list('ABCDEFGHIJ'),
         '6FM2-assembly1':['A','B','A-2','B-2'],'2BAP':list('ABCD'),
         '1FOE':['A','B'],'4JD2':['H']}

def make_frame(table, key):
    all_points=np.concatenate([xyz(table[c]) for c in CHOICES[key]])
    origin=all_points.mean(0)
    basis=np.eye(3)
    if key=='7TPT':
        origin=xyz(table['P']).mean(0)
        x=xyz(table['U']).mean(0)-xyz(table['L']).mean(0)
        z=xyz(table['K']).mean(0)-xyz(table['H']).mean(0)
    elif key=='9B27':
        origin=xyz(table['D']).mean(0)
        x=xyz(table['F']).mean(0)-xyz(table['A']).mean(0)
        z=xyz(table['H']).mean(0)-xyz(table['G']).mean(0)
    elif key=='6VAO':
        origin=xyz(table['D']).mean(0)
        x=xyz(table['B']).mean(0)-xyz(table['C']).mean(0)
        z=xyz(table['J']).mean(0)-origin
    else:
        return origin,basis
    x/=np.linalg.norm(x)
    z-=x*np.dot(z,x);z/=np.linalg.norm(z)
    return origin,np.array([x,np.cross(z,x),z]).T

def main():
    overview={}
    for key,chosen in CHOICES.items():
        out=ROOT/'assets'/key
        if (out/'manifest.json').exists():
            overview[key]=json.loads((out/'manifest.json').read_text())
            continue
        path=ROOT/'sources'/(key+'.cif')
        provenance=json.loads(path.with_suffix('.cif.source.json').read_text())
        assert hashlib.sha256(path.read_bytes()).hexdigest()==provenance['sha256']
        table=atom_table(path)
        origin,basis=make_frame(table,key)
        chains=[];ligands=[]
        model=gemmi.read_structure(str(path))[0]
        out.mkdir(exist_ok=True)
        for c in chosen:
            rows=table[c]
            coords=(xyz(rows)-origin)@basis
            verts,faces=surface(coords,1.15,1.25,.55)
            ca=[r for r in rows if r[2]=='CA']
            np.savez_compressed(out/(c+'.npz'),vertices=verts*.1,faces=faces,
                                atoms=coords.astype(np.float32)*.1,
                                atom_names=np.array([r[2] for r in rows]),
                                residue_numbers=np.array([r[3] for r in rows]),
                                elements=np.array([r[1] for r in rows]),
                                ca=(xyz(ca)-origin)@basis*.1,
                                ca_residues=np.array([r[3] for r in ca]))
            chains.append({'id':c,'atoms':len(rows),'vertices':len(verts),'faces':len(faces),
                           'center_nm':(coords.mean(0)*.1).tolist()})
            for residue in model[c]:
                if residue.entity_type != gemmi.EntityType.NonPolymer:continue
                ats=[a for a in residue if not a.element.is_hydrogen and a.occ>0 and a.altloc in ('\x00','A')]
                if not ats:continue
                ligands.append({'chain':c,'ccd':residue.name,'residue':str(residue.seqid),
                                 'names':[a.name for a in ats], 'elements':[a.element.name for a in ats],
                                 'xyz_nm':((np.array([list(a.pos) for a in ats])-origin)@basis*.1).tolist()})
        info={'source':provenance,'chains':chains,'ligands':ligands,'units':'nm',
              'origin_angstrom':origin.tolist(),'basis_columns':basis.tolist(),
              'representation':'Source-coordinate Gaussian envelope, not experimental density; entire modeled chains retained.',
              'selection_note':'Phalloidin chains omitted from 7TPT display; its mature geometry is a stabilized experimental reference, not a drug-free dynamic observation. Waters/hydrogens omitted.'}
        if key in {'7TPT','9B27'}:
            ca,cb=('J','K') if key=='7TPT' else ('E','F')
            trans,error=rigid_fit(table[ca],table[cb])
            frame=np.eye(4);frame[:3,:3]=basis.T*.1;frame[:3,3]=-basis.T@origin*.1
            info['step']= (frame@trans@np.linalg.inv(frame)).tolist()
            info['step_fit_rmsd_A']=error
        (out/'manifest.json').write_text(json.dumps(info,ensure_ascii=False,indent=2)+'\n')
        overview[key]=info
        print(key,len(chains),'chains prepared',flush=True)
    # Keep the same surface topology for continuous G-to-F teaching morphs.
    # F endpoint is exact source envelope; G endpoint is a CA-guided elastic
    # approximation from bovine profilin-bound 2BTF to rabbit 9B27 F-actin.
    # Neither the interpolated path nor its duration is an experimental result.
    g=np.load(ROOT/'assets/2BTF/A.npz'); f=np.load(ROOT/'assets/9B27/F.npz')
    gi={int(n):p for n,p in zip(g['ca_residues'],g['ca'])}
    fi={int(n):p for n,p in zip(f['ca_residues'],f['ca'])}
    shared=sorted(gi.keys()&fi.keys())
    gp=np.array([gi[n] for n in shared]);fp=np.array([fi[n] for n in shared])
    u,_,vt=np.linalg.svd((gp-gp.mean(0)).T@(fp-fp.mean(0)))
    rot=vt.T@np.diag([1,1,np.linalg.det(vt.T@u.T)])@u.T
    shift=fp.mean(0)-rot@gp.mean(0)
    delta=gp@rot.T+shift-fp
    dist,idx=cKDTree(fp).query(f['vertices'],k=12)
    weights=1/(dist**2+.04);weights/=weights.sum(1)[:,None]
    gverts=f['vertices']+np.sum(delta[idx]*weights[:,:,None],axis=1)
    dist,idx=cKDTree(fp).query(f['atoms'],k=12)
    weights=1/(dist**2+.04);weights/=weights.sum(1)[:,None]
    gatoms=f['atoms']+np.sum(delta[idx]*weights[:,:,None],axis=1)
    np.savez_compressed(ROOT/'assets/actin_morph.npz',f_vertices=f['vertices'],g_vertices=gverts,
                        faces=f['faces'],f_atoms=f['atoms'],g_atoms=gatoms,g_to_f_rotation=rot,g_to_f_shift=shift)
    prof=np.load(ROOT/'assets/2BTF/P.npz')
    np.savez_compressed(ROOT/'assets/profilin_on_morph.npz',vertices=prof['vertices']@rot.T+shift,
                        faces=prof['faces'],atoms=prof['atoms']@rot.T+shift)
    def fit(source,target):
        def cas(d):
            mask=d['atom_names']=='CA'
            return dict(zip(d['residue_numbers'][mask].tolist(),d['atoms'][mask]))
        a,b=cas(source),cas(target);common=sorted(a.keys()&b.keys())
        a=np.array([a[n] for n in common]);b=np.array([b[n] for n in common])
        u,_,vt=np.linalg.svd((a-a.mean(0)).T@(b-b.mean(0)))
        r=vt.T@np.diag([1,1,np.linalg.det(vt.T@u.T)])@u.T
        t=b.mean(0)-r@a.mean(0)
        m=np.eye(4);m[:3,:3]=r;m[:3,3]=t
        return m.tolist(),float(np.sqrt(np.mean(np.sum((a@r.T+t-b)**2,axis=1))))
    cap=np.load(ROOT/'assets/filaments/cap/chain_F.npz')
    matrix,rmsd=fit(cap,f)
    (ROOT/'assets/cap-fit.json').write_text(json.dumps({'matrix':matrix,'CA_RMSD_nm':rmsd,'source':'8F8Q chain F','target':'9B27 chain F'},indent=2))
    matrix,rmsd=fit(f,np.load(ROOT/'assets/7TPT/K.npz'))
    (ROOT/'assets/actin-target-fits.json').write_text(json.dumps({'branch_K':matrix,'CA_RMSD_nm':rmsd,'source':'9B27 chain F','target':'7TPT chain K'},indent=2))
    head_fits={}
    for first,second in [('G','H'),('H','G')]:
        matrix,rmsd=fit(np.load(ROOT/'assets/9B27'/f'{first}.npz'),np.load(ROOT/'assets/9B27'/f'{second}.npz'))
        head_fits[f'{first}_to_{second}']=matrix
        head_fits[f'{first}_to_{second}_CA_RMSD_nm']=rmsd
    (ROOT/'assets/formin-head-fit.json').write_text(json.dumps(head_fits,indent=2))
    (ROOT/'assets/manifest.json').write_text(json.dumps({
      'structures':list(overview),'morph':{'source_G':'2BTF bovine beta-actin, crystallographic Sr/ATP condition',
      'source_F':'9B27 rabbit skeletal alpha-actin ADP state','method':'CA-guided inverse-distance deformation of F surface; schematic state interpolation, not MD',
      'matched_CA':len(shared),'max_CA_displacement_nm':float(np.max(np.linalg.norm(delta,axis=1))),
      'limits':'Chemical states in the growing filament are teaching labels; exact nucleotide/ion states are separately shown using 8A2R/S/T. Cross-isoform morph is an explicit structural analogy.'}},indent=2)+'\n')

if __name__=='__main__':main()
