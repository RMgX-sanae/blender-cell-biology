"""Construct the next U addition after the 5FLM C20 teaching example.

This is constrained geometric modeling, not a force field or catalytic trajectory.
The next template base is T:16 A. An observed A-T pair supplies orientation;
CCD UTP supplies intact chemistry. The A16 nucleotide may move rigidly to close
the RNA link while retaining both DNA backbone connections.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

import gemmi
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation


def fit(x, y):
    cx, cy = x.mean(0), y.mean(0)
    u, _, vt = np.linalg.svd((x-cx).T @ (y-cy))
    if np.linalg.det(vt.T @ u.T) < 0: vt[-1] *= -1
    rotation = vt.T @ u.T
    return rotation, cy-rotation @ cx


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path,required=True)
    p.add_argument('--substrate-dir',type=Path,required=True)
    p.add_argument('--ccd',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Choose a fresh output directory.')
    source=json.loads((a.asset_dir/'manifest.json').read_text())
    first=json.loads((a.substrate_dir/'manifest.json').read_text())
    if source['source_sha256']!=first['source_asset']['source_sha256']:
        p.error('Use the same prepared 5FLM as the first C cycle.')
    data={c:np.load(a.asset_dir/f'chain_{c}.npz') for c in 'PT'}
    def residue(c,r):
        d=data[c];mask=d['residue_numbers']==r
        return {str(n):x.astype(float) for n,x in zip(d['atom_names'][mask],d['atoms'][mask])}
    template=residue('T',16);primer=residue('P',20)
    if set(data['T']['residue_names'][data['T']['residue_numbers']==16])!={'DA'}:
        p.error('The next template position must be T:16 adenine.')
    block=gemmi.cif.read_file(str(a.ccd)).sole_block()
    if block.name!='UTP':p.error('Expected CCD UTP.')
    rows=block.find('_chem_comp_atom.',['atom_id','type_symbol',
        'pdbx_model_Cartn_x_ideal','pdbx_model_Cartn_y_ideal','pdbx_model_Cartn_z_ideal'])
    rows=[r for r in rows if r[1]!='H']
    names=[gemmi.cif.as_string(r[0]) for r in rows]
    elements=np.array([str(r[1]) for r in rows])
    coords=np.array([[float(r[k])*.1 for k in (2,3,4)] for r in rows])
    purine=['N9','C8','N7','C5','C6','N1','C2','N3','C4']
    base=['N1','C2','N3','C4','C5','C6','O2','O4']
    adenine,thymine=residue('P',14),residue('T',23)
    rotation,shift=fit(np.array([adenine[n] for n in purine]),np.array([template[n] for n in purine]))
    target=np.array([thymine[n] for n in base]) @ rotation.T+shift
    rotation,shift=fit(np.array([coords[names.index(n)] for n in base]),target)
    coords=coords @ rotation.T+shift
    extra=['O3A','PB','O1B','O2B','O3B','PG','O1G','O2G','O3G']
    core=[n for n in names if n not in extra]
    ids={n:i for i,n in enumerate(names)}
    core_ids=[ids[n] for n in core]
    axes=[('N1',"C1'",[n for n in names if n not in base]),
          ("C4'","C5'",["C5'","O5'",'PA','O1A','O2A']+extra),
          ("C5'","O5'",["O5'",'PA','O1A','O2A']+extra),
          ("O5'",'PA',['PA','O1A','O2A']+extra)]
    axes=[(ids[x],ids[y],np.array([ids[n] for n in moving])) for x,y,moving in axes]
    origin=template["C1'"]
    def poses(theta):
        rot=Rotation.from_rotvec(theta[4:7]).as_matrix();shift=theta[7:10]
        result=(coords-origin) @ rot.T+origin+shift
        moved={n:rot @ (v-origin)+origin+shift for n,v in template.items()}
        for angle,(x,y,moving) in zip(theta[:4],axes):
            point=result[y];axis=point-result[x];axis/=np.linalg.norm(axis)
            v=result[moving]-point
            result[moving]=point+v*np.cos(angle)+np.cross(axis,v)*np.sin(angle)+np.outer(v @ axis,axis)*(1-np.cos(angle))
        return result,moved
    environment=np.concatenate([data['P']['atoms'],data['T']['atoms']])
    template_ids=np.flatnonzero(data['T']['residue_numbers']==16)+len(data['P']['atoms'])
    mask=np.ones((len(core),len(environment)),bool)
    for i,n in enumerate(core):
        if n in ['PA','O1A','O2A',"O5'"]:
            mask[i,:len(data['P']['atoms'])] &= ~((data['P']['residue_numbers']==20)&
                np.isin(data['P']['atom_names'],["O3'","C3'"]))
    neighbors=np.concatenate([data['T']['atoms'][data['T']['residue_numbers']<15],
        data['T']['atoms'][data['T']['residue_numbers']>17],data['P']['atoms']])
    left,right=residue('T',15)["O3'"],residue('T',17)['P']
    def residual(theta):
        co,moved=poses(theta);pa=co[ids['PA']];bond=pa-primer["O3'"]
        unit=bond/np.linalg.norm(bond)
        tangent=primer["C3'"]-primer["O3'"];tangent/=np.linalg.norm(tangent)
        errors=[(np.linalg.norm(bond)-.16)*30,(tangent @ unit+.5)*.25]
        for n in ["O5'",'O1A','O2A']:
            v=co[ids[n]]-pa;v/=np.linalg.norm(v)
            errors.append((-unit @ v+1/3)*.20)
        for x,y in [(moved['P'],left),(moved["O3'"],right)]:
            errors.append((np.linalg.norm(x-y)-.16)*30)
        errors.extend(theta[4:7]*.03);errors.extend(theta[7:10]*.15)
        ev=environment.copy();ev[template_ids]=list(moved.values())
        errors.extend((np.maximum(0,.24-cdist(ev[template_ids],neighbors))*2).ravel())
        errors.extend((np.maximum(0,.245-cdist(co[core_ids],ev))*mask*2).ravel())
        return errors
    # Seed found by a fixed-seed multistart fit for this specific 5FLM scaffold.
    seed=np.array([-2.683143419,4.896986867,.337571705,1.115742017,
        -.196269884,-.092856189,.122646598,.009085368,-.017813462,-.004470042])
    solution=least_squares(residual,seed,max_nfev=400,ftol=1e-10)
    coords,moved=poses(solution.x)
    length=float(np.linalg.norm(coords[ids['PA']]-primer["O3'"]))
    dna_lengths=[float(np.linalg.norm(moved['P']-left)),float(np.linalg.norm(moved["O3'"]-right))]
    ev=environment.copy();ev[template_ids]=list(moved.values())
    distances=cdist(coords[core_ids],ev);distances[~mask]=np.inf
    if not (.15<length<.17 and all(.15<x<.17 for x in dna_lengths) and distances.min()>.22):
        raise ValueError('Next-cycle modeled geometry failed linkage/clearance checks.')
    ordered=core+extra;order=[ids[n] for n in ordered];index={n:i for i,n in enumerate(ordered)}
    bonds=[];orders=[]
    for x,y,bond_order in block.find('_chem_comp_bond.',['atom_id_1','atom_id_2','value_order']):
        x,y=gemmi.cif.as_string(x),gemmi.cif.as_string(y)
        if x in index and y in index:bonds.append([index[x],index[y]]);orders.append(str(bond_order))
    coords=coords[order];bonds=np.array(bonds,dtype=np.int32)
    lengths=np.linalg.norm(coords[bonds[:,0]]-coords[bonds[:,1]],axis=1)
    if not (.08<lengths.min() and lengths.max()<.20):raise ValueError('CCD bond lengths changed.')
    # Put the newly loaded template A at the former G insertion site. This is
    # a display transform of the entire scaffold, not an experimental state.
    former=residue('T',17)
    trans_r,trans_s=fit(np.array([moved[n] for n in purine]),np.array([former[n] for n in purine]))
    local_r=Rotation.from_rotvec(solution.x[4:7]).as_matrix()
    local_s=origin-local_r @ origin+solution.x[7:10]
    a.output.mkdir(parents=True)
    shutil.copy2(a.substrate_dir/'ctp.npz',a.output/'ctp.npz')
    np.savez_compressed(a.output/'utp.npz',atoms=coords.astype(np.float32),atom_names=np.array(ordered),
        elements=elements[order],bonds=bonds,bond_orders=np.array(orders),product_atom_count=len(core),
        product_names=np.array([{'PA':'P','O1A':'OP1','O2A':'OP2'}.get(n,n) for n in core]),
        product_atoms=coords[:len(core)].astype(np.float32),
        template_rotation=local_r,template_translation=local_s,
        translocation_rotation=trans_r,translocation_translation=trans_s)
    next_info={'component':'UTP','url':'https://files.rcsb.org/ligands/download/UTP.cif',
        'sha256':hashlib.sha256(a.ccd.read_bytes()).hexdigest(),
        'product':{'chain':'P','residue':21,'name':'U','evidence':'modeled, absent from 5FLM'},
        'template_partner':{'chain':'T','residue':16,'name':'DA'},
        'pair_orientation_reference':'5FLM P14 A / T23 T; UTP replaces the pyrimidine in a constructed A-U pair',
        'fit_parameters':solution.x.tolist(),'fit_residual_norm':float(np.linalg.norm(solution.fun)),
        'new_RNA_bond_length_nm':length,'template_backbone_link_lengths_nm':dna_lengths,
        'product_core_nonbonded_clearance_nm':float(distances.min()),
        'CCD_bond_length_range_nm':[float(lengths.min()),float(lengths.max())],
        'template_A16_max_atom_displacement_nm':float(max(np.linalg.norm(moved[n]-v) for n,v in template.items())),
        'limitations':['U21 is an explicitly constructed product, not a deposited experimental residue.',
            'DNA A16 is moved as a complete nucleotide. Two backbone link lengths are constrained; junction angles are not force-field validated.',
            'Global scaffold translocation is an authored one-position display transform; no protein interface or molecular dynamics is reconstructed.',
            'The geometric fit preserves CCD bonds, but does not validate protonation, Mg chemistry or the transition state.']}
    first['next_cycle']=next_info
    (a.output/'manifest.json').write_text(json.dumps(first,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(next_info,ensure_ascii=False))


if __name__=='__main__':main()
