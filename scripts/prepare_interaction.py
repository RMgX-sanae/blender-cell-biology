"""Prepare selected intact interaction partners and a screened teaching path."""
import argparse
import hashlib
import json
from pathlib import Path
import re

import gemmi
import numpy as np

from prepare_molecular_asset import atom_table, xyz, surface
from interaction_motion import plan


def select_ligands(model, example):
    """Resolve required source instances before writing any prepared output."""
    available = {}
    for chain in model:
        for residue in chain:
            if residue.entity_type != gemmi.EntityType.NonPolymer: continue
            key = (chain.name, residue.seqid.num, residue.seqid.icode.strip(), residue.name)
            if key in available: raise ValueError('Ambiguous nonpolymer source identity: '+str(key))
            atoms = [a for a in residue if not a.element.is_hydrogen and a.occ > 0 and a.altloc in ('\x00','A')]
            available[key] = (chain.name, residue, atoms)
    selected = {'anchor':[], 'moving':[]}
    codes = example.get('ligands', [])
    if not isinstance(codes, list) or any(not isinstance(c,str) or not c for c in codes) or len(codes) != len(set(codes)):
        raise ValueError('ligands must be a list of unique required CCD component ids')
    if 'ligand_instances' in example:
        if codes: raise ValueError('Choose ligand_instances or the legacy ligands list, not both')
        if not example.get('ligand_selection_note'): raise ValueError('Explain ligand selection and rigid-partner ownership')
        used = set()
        for item in example['ligand_instances']:
            role = item['owner_role']; icode = item.get('insertion_code', '')
            if role not in selected or type(item['residue_number']) is not int or not isinstance(icode,str) or len(icode) > 1:
                raise ValueError('Require anchor/moving ownership, integer author residue number and an optional insertion code')
            key = (item['chain'], item['residue_number'], icode.strip(), item['component'])
            if key in used: raise ValueError('A ligand instance cannot be assigned twice: '+str(key))
            if key not in available: raise ValueError('Required nonpolymer instance missing: '+str(key))
            selected[role].append(available[key]); used.add(key)
    else:
        found = set()
        for role in selected:
            for chain in example[role]:
                for key, value in available.items():
                    if key[0] == chain and key[3] in codes:
                        selected[role].append(value); found.add(key[3])
        if set(codes)-found:
            raise ValueError('Required ligands missing from selected author chains: '+', '.join(sorted(set(codes)-found)))
    for entries in selected.values():
        for chain, residue, atoms in entries:
            if not atoms: raise ValueError('Selected ligand has no retained heavy atoms: '+chain+':'+str(residue.seqid))
            if len({a.name for a in atoms}) != len(atoms):
                raise ValueError('Selected ligand has ambiguous atom names after altloc filtering')
    return selected


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
    ligands=select_ligands(model,e)
    names=[role+'_'+chain for role in ligands for chain in e[role]]
    names += [role+'_'+chain+'_'+residue.name+'_'+str(residue.seqid)
              for role,entries in ligands.items() for chain,residue,_ in entries]
    if any(not re.fullmatch(r'[A-Za-z0-9_-]+',name) for name in names) or len({name.casefold() for name in names}) != len(names):
        raise ValueError('Actor identities require distinct filename-safe names, including on case-insensitive filesystems')
    origin=np.concatenate([xyz(table[c]) for c in selected]).mean(0)
    a.output.mkdir(parents=True);records=[];sets={}
    def write_ligand(role,chain,residue,atoms):
        coordinates=(np.array([list(atom.pos) for atom in atoms])-origin)*.1
        ligand_id=role+'_'+chain+'_'+residue.name+'_'+str(residue.seqid)
        np.savez_compressed(a.output/(ligand_id+'.npz'),atoms=coordinates,elements=np.array([atom.element.name for atom in atoms]))
        sets[role].append(coordinates)
        records.append({'id':ligand_id,'role':role,'chain':chain,'kind':'ligand','ccd':residue.name,
            'residue_number':residue.seqid.num,'insertion_code':residue.seqid.icode.strip(),
            'atom_names':[atom.name for atom in atoms],'file':ligand_id+'.npz','atoms':len(coordinates)})
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
            for source_chain,residue,atoms_ligand in ligands[role]:
                if source_chain==chain:write_ligand(role,source_chain,residue,atoms_ligand)
        for source_chain,residue,atoms_ligand in ligands[role]:
            if source_chain not in e[role]:write_ligand(role,source_chain,residue,atoms_ligand)
    poses,report=plan(np.concatenate(sets['moving']),np.concatenate(sets['anchor']),frames=a.frames,mode=a.mode)
    np.savez_compressed(a.output/'motion.npz',matrices=poses)
    manifest={'example':e,'source':meta,'source_origin_angstrom':origin.tolist(),'units':'nm',
              'actors':records,'fps':a.fps,'frames':a.frames,'mode':a.mode,'motion':report,
              'data_sha256':{name:hashlib.sha256((a.output/name).read_bytes()).hexdigest()
                             for name in [*(actor['file'] for actor in records),'motion.npz']},
              'events':['approach','reorient','source-pose binding','release'],
              'geometry':'Complete modeled selected chains, plus explicitly selected ligands. Unresolved regions are not reconstructed.',
              'limits':'Local interaction cycle only. No catalytic transformation or complete cellular pathway is implied.'}
    (a.output/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
    print(json.dumps({'output':str(a.output),'motion':report},ensure_ascii=False))


if __name__=='__main__':main()
