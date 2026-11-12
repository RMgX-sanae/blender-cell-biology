"""Resolve explicit polymer/companion instances from deposited mmCIF assemblies.

Coordinates stay in Angstrom. A source operation is applied once to every atom
of its selected instance; author chain/residue identifiers are never renamed.
Instance ids distinguish repeated copies. Preparation runs outside Blender.
"""
import hashlib
from itertools import product
import json
from pathlib import Path
import re

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from prepare_clathrin_assembly import rigid


def category_rows(block,category):
    columns=block.get_mmcif_category(category)
    return [dict(zip(columns,values)) for values in zip(*columns.values())]


def operation_products(expression,operators):
    """Expand deposited lists/ranges/products; the rightmost matrix acts first.

    Gemmi 0.7.5's assembly reader ignores products, so read the raw category.
    Keep source matrix precision; never fit or normalize deposited operators.
    """
    expression=''.join(expression.split())
    if expression.startswith('('):
        groups=re.findall(r'\(([^()]*)\)',expression)
        if ''.join('('+g+')' for g in groups)!=expression:raise ValueError('Invalid assembly operation expression')
    else:
        groups=[expression]
    expanded=[]
    for group in groups:
        names=[]
        for token in group.split(','):
            if token in operators:names.append(token)
            elif re.fullmatch(r'\d+-\d+',token):
                first,last=map(int,token.split('-'))
                if last<first:raise ValueError('Descending assembly operation range')
                names.extend(str(i) for i in range(first,last+1))
            else:raise ValueError('Missing or invalid assembly operation: '+token)
        if not names or len(set(names))!=len(names) or set(names)-operators.keys():
            raise ValueError('Missing or duplicate assembly operations')
        expanded.append(names)
    for factors in product(*expanded):
        matrix=np.eye(4)
        for name in factors:matrix=matrix@operators[name]
        rigid(matrix)
        name=factors[0] if len(factors)==1 else ''.join('('+f+')' for f in factors)
        yield name,matrix


def retained_rows(residues):
    rows=[]
    for residue in residues:
        for atom in residue:
            if atom.element.is_hydrogen or atom.occ<=0 or atom.altloc not in ('\x00','A'):continue
            point=np.asarray(list(atom.pos),float)
            if not np.isfinite(point).all():raise ValueError('Nonfinite source atom position')
            rows.append((point.tolist(),atom.element.name,atom.name,residue.seqid.num,
                         residue.name,residue.seqid.icode.strip()))
    keys=[(r[3],r[5],r[2]) for r in rows]
    if not rows or len(keys)!=len(set(keys)):
        raise ValueError('Empty or ambiguous retained source atom identities')
    return rows


def fragments(model):
    """Keep label/auth chains separate; nonpolymers remain individual residues."""
    groups={}
    for chain in model:
        for residue in chain:
            label=residue.subchain
            if not label:raise ValueError('Assembly sources require explicit mmCIF label asym ids')
            if label not in groups:groups[label]={'chain':chain.name,'residues':[]}
            if groups[label]['chain']!=chain.name:raise ValueError('Ambiguous label/auth chain correspondence')
            groups[label]['residues'].append(residue)
    return groups


class SourceAssembly:
    def __init__(self,path,assembly_id):
        self.path=Path(path)
        self.structure=gemmi.read_structure(str(path))
        if len(self.structure)!=1:raise ValueError('Select a single reviewed coordinate model before assembly preparation')
        self.groups=fragments(self.structure[0])
        block=gemmi.cif.read(str(path)).sole_block()
        self.entities={r['id']:r for r in category_rows(block,'_entity.')}
        self.polymer_types={r['entity_id']:r['type'] for r in category_rows(block,'_entity_poly.')}
        matches=[a for a in category_rows(block,'_pdbx_struct_assembly.') if a['id']==assembly_id]
        if len(matches)!=1:raise ValueError('Missing or ambiguous deposited assembly')
        self.metadata=matches[0]
        self.assembly_id=assembly_id;self.operations={}
        operators={}
        for record in category_rows(block,'_pdbx_struct_oper_list.'):
            matrix=np.eye(4)
            for i in range(3):
                matrix[i,3]=float(record[f'vector[{i+1}]'])
                for j in range(3):matrix[i,j]=float(record[f'matrix[{i+1}][{j+1}]'])
            rigid(matrix)
            if record['id'] in operators:raise ValueError('Duplicate source operation id')
            operators[record['id']]=matrix
        generators=[g for g in category_rows(block,'_pdbx_struct_assembly_gen.') if g['assembly_id']==assembly_id]
        if not generators:raise ValueError('Missing assembly generators')
        for generator in generators:
            labels=[label.strip() for label in generator['asym_id_list'].split(',')]
            if not labels or len(labels)!=len(set(labels)) or set(labels)-self.groups.keys():
                raise ValueError('Assembly generator has missing or ambiguous source labels')
            for name,matrix in operation_products(generator['oper_expression'],operators):
                for label in labels:
                    key=(label,name)
                    if key in self.operations:raise ValueError('Duplicate assembly instance: '+str(key))
                    self.operations[key]=matrix.copy()

    def member(self,spec,kind):
        instance=spec.get('instance_id')
        if not isinstance(instance,str) or not re.fullmatch(r'[A-Za-z0-9_-]+',instance):
            raise ValueError('Require a filename-safe explicit assembly instance id')
        key=(spec['label_chain'],spec['operator'])
        if key not in self.operations:raise ValueError('Selected instance is outside the deposited assembly: '+str(key))
        group=self.groups[spec['label_chain']]
        if spec['chain']!=group['chain']:raise ValueError('Author chain does not match the selected label chain')
        residues=group['residues']
        if kind=='polymer':
            if any(r.entity_type!=gemmi.EntityType.Polymer for r in residues):
                raise ValueError('Selected label is not a polymer')
        elif kind=='companion':
            number=spec['residue_number'];icode=spec.get('insertion_code','')
            if type(number) is not int or not isinstance(icode,str) or len(icode)>1:
                raise ValueError('Companions require integer author numbers and an optional insertion code')
            residues=[r for r in residues if r.entity_type in (gemmi.EntityType.NonPolymer,gemmi.EntityType.Water)
                      and r.seqid.num==number and r.seqid.icode.strip()==icode and r.name==spec['component']]
            if len(residues)!=1:raise ValueError('Missing or ambiguous required companion')
        else:raise ValueError('Unknown assembly member kind')
        entities={r.entity_id for r in residues}
        if len(entities)!=1 or not entities<=self.entities.keys():raise ValueError('Ambiguous source entity identity')
        entity_id=next(iter(entities));entity=self.entities[entity_id]
        if spec.get('description') and spec['description']!=entity.get('pdbx_description'):
            raise ValueError('Source entity differs from the reviewed description')
        raw=retained_rows(residues);matrix=self.operations[key]
        coordinates=np.asarray([r[0] for r in raw])@matrix[:3,:3].T+matrix[:3,3]
        transformed=[(point.tolist(),*row[1:]) for point,row in zip(coordinates,raw)]
        return {**spec,'kind':kind,'source_rows':raw,'rows':transformed,'entity_id':entity_id,
                'source_description':entity.get('pdbx_description'),'source_mutation_field':entity.get('pdbx_mutation'),
                **({'polymer_type':self.polymer_types[entity_id]} if kind=='polymer' else {}),
                'source_operation_angstrom':matrix.tolist()}

    def select(self,polymers,companions=()):
        if not polymers:raise ValueError('Select at least one complete polymer instance')
        selected=[self.member(spec,'polymer') for spec in polymers]
        owners={m['instance_id']:m for m in selected}
        for spec in companions:
            owner=owners.get(spec.get('owner_instance'))
            if owner is None or spec['operator']!=owner['operator']:
                raise ValueError('Companion must follow a declared owner in the same source operation')
            selected.append(self.member(spec,'companion'))
        ids=[m['instance_id'].casefold() for m in selected]
        keys=[(m['kind'],m['label_chain'],m['operator'],
               m.get('residue_number'),m.get('insertion_code',''),m.get('component')) for m in selected]
        if len(ids)!=len(set(ids)) or len(keys)!=len(set(keys)):
            raise ValueError('Assembly instances may not be duplicated or aliased')
        return selected


def verify_reference(members,path,expected_sha256,*,complete=False):
    """Match full atom identities/positions to an independently downloaded assembly.

    The 0.001 Angstrom bound is the existing 0.0001 nm assembly check, not a
    fitting tolerance. No coordinate refitting or source-matrix normalization.
    """
    path=Path(path)
    ids=[m['instance_id'].casefold() for m in members]
    if len(ids)!=len(set(ids)):raise ValueError('Duplicate comparison instance identities')
    if hashlib.sha256(path.read_bytes()).hexdigest()!=expected_sha256:
        raise ValueError('Independent assembly source hash changed')
    structure=gemmi.read_structure(str(path))
    if len(structure)!=1:raise ValueError('Independent assembly must contain one coordinate model')
    reference=[]
    for label,group in fragments(structure[0]).items():
        polymers=[r for r in group['residues'] if r.entity_type==gemmi.EntityType.Polymer]
        if polymers:reference.append((label,group['chain'],'polymer',retained_rows(polymers)))
        for residue in group['residues']:
            if residue.entity_type in (gemmi.EntityType.NonPolymer,gemmi.EntityType.Water):
                reference.append((label,group['chain'],'companion',retained_rows([residue])))
    if not reference or not members:raise ValueError('Empty assembly comparison')
    coordinates=[np.array([r[0] for r in entry[3]]) for entry in reference]
    identities=[{tuple(row[1:]):i for i,row in enumerate(entry[3])} for entry in reference]
    tree=cKDTree([points.mean(0) for points in coordinates]);used=set();matched={};errors=[]
    for member in members:
        rows=member['rows'];points=np.array([r[0] for r in rows]);candidates=[]
        for index in tree.query_ball_point(points.mean(0),np.sqrt(3)*.001):
            entry=reference[index]
            keys=[tuple(r[1:]) for r in rows]
            if entry[2]!=member['kind'] or len(keys)!=len(identities[index]) or set(keys)!=identities[index].keys():continue
            target=coordinates[index][[identities[index][key] for key in keys]]
            error=float(np.abs(points-target).max())
            if error<=.001:candidates.append((index,error))
        if len(candidates)!=1 or candidates[0][0] in used:
            raise ValueError('Assembly identity/coordinates do not match one independent instance: '+member['instance_id'])
        index,error=candidates[0];used.add(index);errors.append(error)
        matched[member['instance_id']]={'label_chain':reference[index][0],'chain':reference[index][1]}
    if complete and len(used)!=len(reference):raise ValueError('Assembly selection omits independent reference instances')
    return {'passed':True,'matched_instances':len(used),'reference_instances':len(reference),
            'maximum_coordinate_error_angstrom':max(errors),'tolerance_angstrom':.001,
            'complete_reference_inventory':len(used)==len(reference),'instance_map':matched,
            'reference_sha256':expected_sha256}


def reviewed_members(path,spec):
    """Load explicit instances and require an independent assembly source check."""
    assembly=SourceAssembly(path,spec['assembly_id'])
    members=assembly.select(spec['chains'],spec['ligand_instances'])
    reference=Path(path).with_name(spec['pdb']+'-assembly'+spec['assembly_id']+'.cif')
    checksum=spec['assembly_reference']['sha256']
    meta=json.loads(reference.with_suffix('.cif.source.json').read_text())
    if (meta['id'].upper()!=spec['pdb'] or str(meta['assembly'])!=spec['assembly_id']
            or meta['sha256']!=checksum):raise ValueError('Independent assembly provenance differs')
    report=verify_reference(members,reference,checksum,
                            complete=spec['assembly_reference'].get('complete',False))
    return assembly,members,report
