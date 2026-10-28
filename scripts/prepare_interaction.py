"""Prepare selected intact interaction partners and a screened teaching path."""
import argparse
import hashlib
import json
from pathlib import Path

import gemmi
import numpy as np

from prepare_molecular_asset import atom_table, xyz, surface
from interaction_motion import plan


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--examples',type=Path,default=Path(__file__).resolve().parents[1]/'assets/interaction-examples.json')
    p.add_argument('--example',required=True)
    p.add_argument('--sources',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--mode',choices=('loop','teaching'),default='loop')
    p.add_argument('--frames',type=int,default=192)
    p.add_argument('--fps',type=int,default=24)
    a=p.parse_args()
    if a.frames<12 or a.frames%2 or a.fps<=0:p.error('Use an even frame count >=12 and a positive fps')
    if a.output.exists():p.error('Use a fresh output directory')
    e=next((e for e in json.loads(a.examples.read_text())['examples'] if e['id']==a.example),None)
    if e is None:p.error('Unknown example; provide an explicit reviewed specification')
    source=a.sources/(e['pdb']+'.cif');meta=json.loads(source.with_suffix('.cif.source.json').read_text())
    assert hashlib.sha256(source.read_bytes()).hexdigest()==meta['sha256']
    assert meta['id'].upper()==e['pdb'].upper()
    table=atom_table(source);model=gemmi.read_structure(str(source))[0]
    selected=e['anchor']+e['moving']
    if len(set(selected))!=len(selected):raise ValueError('A chain cannot belong to two rigid partners')
    if set(selected)-set(table):raise ValueError('Author chain selection does not match structure')
    origin=np.concatenate([xyz(table[c]) for c in selected]).mean(0)
    a.output.mkdir(parents=True);records=[];sets={}
    for role in ['anchor','moving']:
        sets[role]=[]
        for chain in e[role]:
            rows=table[chain];atoms=(xyz(rows)-origin)*.1
            vertices,faces=surface(atoms*10,.7,1.1,.55)
            name=role+'_'+chain
            np.savez_compressed(a.output/(name+'.npz'),vertices=vertices*.1,faces=faces,atoms=atoms,
                                atom_names=np.array([r[2] for r in rows]),residue_numbers=np.array([r[3] for r in rows]))
            records.append({'id':name,'role':role,'chain':chain,'kind':'polymer','file':name+'.npz','atoms':len(atoms)})
            sets[role].append(atoms)
            for residue in model[chain]:
                if residue.name not in e['ligands'] or residue.entity_type!=gemmi.EntityType.NonPolymer:continue
                atoms_ligand=[atom for atom in residue if not atom.element.is_hydrogen and atom.occ>0 and atom.altloc in ('\x00','A')]
                if not atoms_ligand:continue
                coordinates=(np.array([list(atom.pos) for atom in atoms_ligand])-origin)*.1
                ligand_id=name+'_'+residue.name+'_'+str(residue.seqid)
                np.savez_compressed(a.output/(ligand_id+'.npz'),atoms=coordinates,elements=np.array([atom.element.name for atom in atoms_ligand]))
                sets[role].append(coordinates)
                records.append({'id':ligand_id,'role':role,'chain':chain,'kind':'ligand','ccd':residue.name,'file':ligand_id+'.npz','atoms':len(coordinates)})
    poses,report=plan(np.concatenate(sets['moving']),np.concatenate(sets['anchor']),frames=a.frames,mode=a.mode)
    np.savez_compressed(a.output/'motion.npz',matrices=poses)
    manifest={'example':e,'source':meta,'source_origin_angstrom':origin.tolist(),'units':'nm',
              'actors':records,'fps':a.fps,'frames':a.frames,'mode':a.mode,'motion':report,
              'events':['approach','reorient','source-pose binding','release'],
              'geometry':'Complete modeled selected chains, plus explicitly selected ligands. Unresolved regions are not reconstructed.',
              'limits':'Local interaction cycle only. No catalytic transformation or complete cellular pathway is implied.'}
    (a.output/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
    print(json.dumps({'output':str(a.output),'motion':report},ensure_ascii=False))


if __name__=='__main__':main()
