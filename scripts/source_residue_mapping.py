"""Map observed C-alpha identities between reviewed sources using SIFTS.

This maps residues, not conformations or biological equivalence. The caller
selects chains/accessions after reviewing the constructs and experimental state.
"""
import gzip
import hashlib
import xml.etree.ElementTree as ET
from pathlib import Path

import gemmi
import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sifts(path, entry_id=None):
    raw=Path(path).read_bytes()
    root=ET.fromstring(gzip.decompress(raw) if raw[:2]==b'\x1f\x8b' else raw)
    if entry_id is not None and root.attrib.get('dbAccessionId','').lower()!=entry_id.lower():
        raise ValueError('SIFTS accession differs from the coordinate entry')
    result={}
    for residue in root.iter():
        if residue.tag.split('}')[-1]!='residue':continue
        refs=[item.attrib for item in residue if item.tag.split('}')[-1]=='crossRefDb']
        pdb=[r for r in refs if r.get('dbSource')=='PDB']
        uniprot=[r for r in refs if r.get('dbSource')=='UniProt']
        if not pdb:continue
        if len(pdb)!=1:raise ValueError('Ambiguous SIFTS PDB identity')
        p=pdb[0]
        if p['dbResNum'] in ('null','None','.','?'):continue
        key=(p['dbChainId'],p['dbResNum'])
        if key in result:raise ValueError('Duplicate SIFTS author identity: '+str(key))
        result[key]={'name':p['dbResName'],'canonical':[(u['dbAccessionId'],u['dbResNum'],u['dbResName']) for u in uniprot]}
    return result


def observed_ca(path, chains):
    result={}
    for chain in gemmi.read_structure(str(path))[0]:
        if chain.name not in chains:continue
        for residue in chain:
            if residue.entity_type!=gemmi.EntityType.Polymer:continue
            atoms=[a for a in residue if a.name=='CA' and a.element.name=='C' and a.occ>0 and a.altloc in ('\x00','A')]
            if not atoms:continue
            key=(chain.name,str(residue.seqid))
            if len(atoms)!=1 or key in result:raise ValueError('Ambiguous observed C-alpha: '+str(key))
            result[key]={'name':residue.name,'xyz_nm':np.array(atoms[0].pos.tolist())*.1}
    if set(chains)!={k[0] for k in result}:raise ValueError('A selected chain has no observed C-alpha atoms')
    return result


def correspondence(first_cif, second_cif, first_sifts, second_sifts, chains, insertions=()):
    """Return observed-source identities and explicit unmatched regions.

    chains: [{first, second, uniprot}]. insertions: [{first:[chain,resid],
    second:[chain,resid], name, reason}], for explicitly reviewed engineered
    insertions without canonical mapping. Insertion codes remain strings.
    """
    if len({r['first'] for r in chains})!=len(chains) or len({r['second'] for r in chains})!=len(chains):
        raise ValueError('Select one-to-one source chains')
    observed=[observed_ca(p,[r[k] for r in chains]) for p,k in [(first_cif,'first'),(second_cif,'second')]]
    entries=[gemmi.cif.read(str(p)).sole_block().find_value('_entry.id') for p in (first_cif,second_cif)]
    if not all(entries):raise ValueError('Require mmCIF entries with explicit PDB identities')
    mapping=[sifts(p,entry) for p,entry in zip((first_sifts,second_sifts),entries)];indices=[]
    for side,role in enumerate(('first','second')):
        index={}
        for spec in chains:
            count=0
            for key,row in observed[side].items():
                if key[0]!=spec[role]:continue
                if key not in mapping[side] or mapping[side][key]['name']!=row['name']:
                    raise ValueError('SIFTS and coordinate residue disagree: '+str(key))
                refs=[r for r in mapping[side][key]['canonical'] if r[0]==spec['uniprot']]
                if len(refs)>1:raise ValueError('Ambiguous canonical mapping: '+str(key))
                if not refs:continue
                ref=refs[0];canonical=(spec['first'],ref[0],ref[1])
                if canonical in index:raise ValueError('Repeated canonical position within a selected chain')
                if gemmi.find_tabulated_residue(row['name']).one_letter_code.upper()!=ref[2].upper():
                    raise ValueError('Mutation or modified residue needs separate review: '+str(key))
                index[canonical]=key;count+=1
            if not count:raise ValueError('No observed residues for selected accession: '+spec['uniprot'])
        indices.append(index)
    pairs=[]
    for canonical,key in indices[0].items():
        other=indices[1].get(canonical)
        if other is None:continue
        if observed[0][key]['name']!=observed[1][other]['name']:raise ValueError('Mapped residues differ in chemical identity')
        pairs.append({'first':list(key),'second':list(other),'name':observed[0][key]['name'],
                      'uniprot':canonical[1],'canonical_residue':canonical[2],'basis':'SIFTS'})
    used=[{tuple(p[role]) for p in pairs} for role in ('first','second')]
    for entry in insertions:
        if not entry.get('reason'):raise ValueError('An insertion match needs a reviewed reason')
        for side,role in enumerate(('first','second')):
            key=tuple(entry[role])
            if key in used[side] or key not in observed[side] or observed[side][key]['name']!=entry['name']:
                raise ValueError('Invalid reviewed insertion match')
            if mapping[side][key]['canonical']:raise ValueError('Do not override a canonical mapping as an insertion')
            used[side].add(key)
        pairs.append({**entry,'basis':'reviewed engineered insertion','uniprot':None,'canonical_residue':None})
    report={'pairs':pairs,'observed_counts':[len(v) for v in observed],
            'unmatched':{role:[{'author':list(key),'name':v['name']} for key,v in observed[side].items() if key not in used[side]]
                         for side,role in enumerate(('first','second'))},
            'sha256':{role:digest(path) for role,path in [('first_cif',first_cif),('second_cif',second_cif),('first_sifts',first_sifts),('second_sifts',second_sifts)]}}
    return observed,report
