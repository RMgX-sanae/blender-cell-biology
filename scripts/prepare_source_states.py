"""Prepare reviewed, separate experimental states; never morph or merge them.

Uses the existing source, polymer, ligand, surface and rigid-fit helpers.
Supports legacy identity selections and explicit, independently checked
assembly instances, with reviewed CA maps between retained polymer instances.
"""
import argparse
import json
from pathlib import Path

import gemmi
import numpy as np

from prepare_binding_candidates import checked_source, fit_template, residue_id, valid_residue_id, sha, transform
from prepare_interaction import select_ligands
from prepare_molecular_asset import atom_table, surface, xyz
from source_assembly import reviewed_members


def selection_key(spec):
    return spec.get('instance_id', spec['chain'])


def rows(block, category):
    columns = block.get_mmcif_category(category)
    return [dict(zip(columns, values)) for values in zip(*columns.values())]


def alignment_maps(config):
    """Resolve reviewed maps without inferring chain or residue equivalence."""
    fit = config['alignment']; states = {s['pdb']: s for s in config['states']}
    if config['reference'] not in states: raise ValueError('Reference state must be selected')
    if 'maps' in fit:
        if any(key in fit for key in ('chain', 'author_residues', 'minimum_CA_coverage', 'maximum_rmsd_angstrom')):
            raise ValueError('Choose explicit per-state maps or the same-numbered alignment, not both')
        reference_chain = fit.get('reference_instance_id', fit['reference_chain']); mappings = fit['maps']
        if set(mappings) != set(states) or any(not m.get('note') for m in mappings.values()):
            raise ValueError('Supply one reviewed mapping and its rationale for every state')
    else:
        reference_chain = fit['chain']
        mapping = {'chain': reference_chain, 'residue_pairs': [[n, n] for n in fit['author_residues']],
                   'expected_substitutions': [], 'min_CA_coverage': fit['minimum_CA_coverage'],
                   'max_CA_rmsd_angstrom': fit['maximum_rmsd_angstrom']}
        mappings = {code: dict(mapping) for code in states}
    if reference_chain not in {selection_key(c) for c in states[config['reference']]['chains']}:
        raise ValueError('Reference fitting chain must be retained')
    reference_spec = next(c for c in states[config['reference']]['chains'] if selection_key(c) == reference_chain)
    if reference_spec['chain'] != fit.get('reference_chain',fit.get('chain')):
        raise ValueError('Reference fitting author chain differs from the selected instance')
    for code, mapping in mappings.items():
        selected = [c for c in states[code]['chains'] if selection_key(c) == selection_key(mapping)]
        if len(selected) != 1 or selected[0]['chain'] != mapping['chain']:
            raise ValueError('State fitting chain must be retained: '+code)
    reference = mappings[config['reference']]
    if selection_key(reference) != reference_chain or any(a != b for a, b in reference['residue_pairs']):
        raise ValueError('The reference state must retain its own chain and residue identities')
    return reference_chain, mappings


def qualify(sources, state):
    path = checked_source(sources, state)
    provenance = json.loads(path.with_suffix('.cif.source.json').read_text())
    block = gemmi.cif.read(str(path)).sole_block()
    model = gemmi.make_structure_from_block(block)[0]
    entities = {r['id']: r for r in rows(block, '_entity.')}
    polymers = rows(block, '_entity_poly.')
    refs = rows(block, '_struct_ref.')
    diffs = rows(block, '_struct_ref_seq_dif.')
    organisms = []
    for category, field in [('_entity_src_gen.', 'pdbx_gene_src_scientific_name'),
                            ('_entity_src_nat.', 'pdbx_organism_scientific'),
                            ('_pdbx_entity_src_syn.', 'organism_scientific')]:
        organisms.extend({'entity_id': r['entity_id'], 'category': category, 'organism': r.get(field)}
                         for r in rows(block, category))
    assemblies = rows(block, '_pdbx_struct_assembly.')
    assembly = next((r for r in assemblies if r['id'] == state['assembly_id']), None)
    generators = [r for r in rows(block, '_pdbx_struct_assembly_gen.')
                  if r['assembly_id'] == state['assembly_id']]
    operators = {r['id']: r for r in rows(block, '_pdbx_struct_oper_list.')}
    included = set()
    if assembly is None or not generators:
        raise ValueError('Missing reviewed assembly')
    explicit = 'assembly_reference' in state or any('instance_id' in c or 'operator' in c or 'label_chain' in c for c in state['chains'])
    members = {}
    if explicit:
        _, selected_members, proof = reviewed_members(path, state)
        members = {m['instance_id']: m for m in selected_members if m['kind'] == 'polymer'}
        table = {key: m['rows'] for key,m in members.items()}
        assembly = {**assembly, 'instance_check': proof}
        ligand_rows = {'assembly_members': [m for m in selected_members if m['kind'] == 'companion']}
    for gen in ([] if explicit else generators):
        # A restricted, explicit boundary: do not silently ignore symmetry.
        op = operators.get(gen['oper_expression'])
        if op is None:
            raise ValueError('Prepare non-identity/multiple-operation assemblies separately')
        matrix = [[float(op[f'matrix[{i}][{j}]']) for j in (1, 2, 3)] for i in (1, 2, 3)]
        shift = [float(op[f'vector[{i}]']) for i in (1, 2, 3)]
        if not np.allclose(matrix, np.eye(3), atol=1e-8, rtol=0) or not np.allclose(shift, 0, atol=1e-8, rtol=0):
            raise ValueError('Prepare non-identity/multiple-operation assemblies separately')
        included.update(gen['asym_id_list'].split(','))
    if not explicit: table = atom_table(path, include_insertion_codes=True)
    selected, info = {}, []
    for spec in state['chains']:
        chain = spec['chain']; key = selection_key(spec)
        if key in selected or key not in table:
            raise ValueError('Missing or duplicate selected polymer chain: '+chain)
        if explicit:
            labels = [spec['label_chain']]
            entity_id = next(r['entity_id'] for r in rows(block, '_struct_asym.') if r['id'] == labels[0])
            matches = [r for r in polymers if r['entity_id'] == entity_id]
        else:
            matches = [r for r in polymers if chain in r['pdbx_strand_id'].split(',')]
            labels = sorted({r.subchain for c in model if c.name == chain for r in c
                             if r.entity_type == gemmi.EntityType.Polymer})
        if len(matches) != 1:
            raise ValueError('Ambiguous author-chain entity: '+chain)
        entity = entities[matches[0]['entity_id']]
        if not explicit and not set(labels) <= included:
            raise ValueError('Selected polymer is outside the reviewed assembly')
        if entity['pdbx_description'] != spec['description']:
            raise ValueError('Source chain identity differs from the reviewed description: '+chain)
        accessions = [r['pdbx_db_accession'] for r in refs
                      if r['entity_id'] == entity['id'] and r['db_name'] == 'UNP']
        if spec.get('uniprot') and spec['uniprot'] not in accessions:
            raise ValueError('Source chain UniProt identity differs: '+chain)
        source_organisms = [r for r in organisms if r['entity_id'] == entity['id']]
        if spec.get('source_organism') and spec['source_organism'] not in {r['organism'] for r in source_organisms}:
            raise ValueError('Source organism differs from the reviewed protein: '+chain)
        selected[key] = table[key]
        actual = {}
        for row in table[key]:
            atom_key = (residue_id(row), row[2])
            if atom_key in actual:
                raise ValueError('Ambiguous retained polymer atom identity: '+str(atom_key))
            actual[atom_key] = row
        for site in spec.get('sites', []):
            if not valid_residue_id(site['residue']):
                raise ValueError('Reviewed sites require integer author numbers or explicit insertion codes')
            site_rows = [r for r in table[key] if residue_id(r) == site['residue']]
            if not site_rows or {r[4] for r in site_rows} != {site['name']}:
                raise ValueError('Reviewed site residue differs or is not retained: '+chain+str(site))
            if not set(site.get('atoms', [])) <= {r[2] for r in site_rows}:
                raise ValueError('Reviewed site atom is not retained: '+chain+str(site))
        operation = {'source_operation_angstrom': members[key]['source_operation_angstrom']} if explicit else {}
        info.append({**spec, **operation, 'entity_id': entity['id'], 'label_chains': labels,
                     'polymer_type': matches[0]['type'],
                     'uniprot_accessions': accessions,
                     'source_organisms': source_organisms,
                     'deposited_mutation_field': entity.get('pdbx_mutation'),
                     'sequence_differences': [r for r in diffs if r.get('pdbx_pdb_strand_id') == chain]})
    if not explicit:
        ligand_spec = {'ligand_instances': state['ligand_instances'],
                       'ligand_selection_note': state['selection_note']}
        ligand_rows = select_ligands(model, ligand_spec)
        for _, residue, _ in ligand_rows['anchor'] + ligand_rows['moving']:
            if residue.subchain not in included:
                raise ValueError('Selected ligand is outside the reviewed assembly')
    return path, provenance, selected, info, ligand_rows, assembly


def prepare(sources, config_path, output):
    if output.exists():
        raise ValueError('Choose a fresh output directory')
    config = json.loads(config_path.read_text())
    states = config['states']
    ids = [s['pdb'] for s in states]
    if len(ids) != len(set(ids)) or len(states) < 2:
        raise ValueError('Require at least two distinct source states')
    if not config.get('limitations') or any(not s.get('selection_note') for s in states):
        raise ValueError('Document selection and state-comparison limits')
    qualified = {s['pdb']: qualify(sources, s) for s in states}
    reference_chain, mappings = alignment_maps(config)
    reference = qualified[config['reference']][2][reference_chain]
    origin = np.asarray(config['origin_angstrom'], dtype=float)
    basis = np.asarray(config['basis_columns'], dtype=float)
    if origin.shape != (3,) or not np.isfinite(origin).all() or basis.shape != (3, 3):
        raise ValueError('Invalid common viewing frame')
    if not np.isfinite(basis).all() or not np.allclose(basis.T @ basis, np.eye(3), atol=1e-7, rtol=0) or not np.isclose(np.linalg.det(basis), 1, atol=1e-7, rtol=0):
        raise ValueError('Viewing basis must be a proper rigid rotation')
    records, arrays = [], {}
    for state in states:
        code = state['pdb']
        path, provenance, selected, info, ligands, assembly = qualified[code]
        mapping = mappings[code]
        matrix, fit_report = fit_template(selected[selection_key(mapping)], reference, mapping)
        def local(points):
            return (transform(points, matrix) - origin) @ basis * .1
        actors = []
        for i, spec in enumerate(info):
            source_rows = selected[selection_key(spec)]
            atoms = local(xyz(source_rows))
            vertices, faces = surface(atoms * 10, .75, 1.15, .55)
            name = f'{code}_polymer_{i}.npz'
            arrays[name] = dict(vertices=vertices * .1, faces=faces, atoms=atoms,
                elements=np.array([r[1] for r in source_rows]), atom_names=np.array([r[2] for r in source_rows]),
                residue_numbers=np.array([r[3] for r in source_rows]),
                insertion_codes=np.array([r[5] for r in source_rows]), residue_names=np.array([r[4] for r in source_rows]))
            actors.append({**spec, 'file': name, 'atoms': len(atoms), 'kind': 'polymer'})
        companions = []
        for member in ligands.get('assembly_members', []):
            source_rows = member['rows']
            identity = {k:v for k,v in member.items() if k not in ('source_rows','rows','kind')}
            companions.append({**identity, 'residue': str(member['residue_number'])+member.get('insertion_code',''),
                'atom_names': [r[2] for r in source_rows], 'elements': [r[1] for r in source_rows],
                'positions_nm': local(xyz(source_rows)).tolist()})
        for chain, residue, atoms in ligands.get('anchor', []) + ligands.get('moving', []):
            companions.append({'chain': chain, 'residue': str(residue.seqid), 'component': residue.name,
                'atom_names': [a.name for a in atoms], 'elements': [a.element.name for a in atoms],
                'positions_nm': local(np.array([list(a.pos) for a in atoms])).tolist()})
        records.append({**state, 'actors': actors, 'ligands': companions,
            'source': provenance, 'assembly': assembly, 'alignment_matrix_angstrom': matrix.tolist(),
            'alignment_check': fit_report})
    output.mkdir(parents=True)
    for name, data in arrays.items():
        np.savez_compressed(output / name, **data)
    manifest = {'config': config, 'config_sha256': sha(config_path), 'units': 'nm', 'states': records,
                'data_sha256': {name: sha(output / name) for name in arrays},
                'interpretation': 'Independent experimental state snapshots; no trajectory, ligand exchange, or reaction modeled.'}
    (output / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'output': str(output), 'states': ids, 'polymer_atoms': [sum(a['atoms'] for a in r['actors']) for r in records]}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.sources, args.config, args.output)
