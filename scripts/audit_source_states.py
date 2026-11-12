"""Check selected raw-source atoms, shared transforms, and saved static panels.

Outside Blender: --sources --config --asset-dir --output.
Inside a reopened Blender scene: --asset-dir --output (geometry/identity check).
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition: raise ValueError(message)


def assets(directory):
    path = directory/'manifest.json'
    m = json.loads(path.read_text())
    require(m['units'] == 'nm', 'Prepared state units must remain nm')
    reviewed_states = m['config']['states']
    require([s['pdb'] for s in m['states']] == [s['pdb'] for s in reviewed_states], 'State inventory mismatch')
    for state, reviewed in zip(m['states'], reviewed_states):
        require(all(state.get(key) == value for key, value in reviewed.items()), 'Prepared state declaration differs from reviewed config')
        require(state['source']['id'].upper() == reviewed['pdb'] and state['source']['sha256'] == reviewed['sha256']
                and state['assembly']['id'] == reviewed['assembly_id'], 'Prepared source/assembly identity differs')
        require(len(state['actors']) == len(reviewed['chains']), 'Polymer chain inventory mismatch')
        for actor, spec in zip(state['actors'], reviewed['chains']):
            require(all(actor.get(key) == value for key, value in spec.items()), 'Prepared polymer declaration differs from reviewed config')
        if any('instance_id' in spec for spec in reviewed['chains']):
            proof=state['assembly'].get('instance_check',{})
            ids={s['instance_id'] for s in reviewed['chains']+reviewed['ligand_instances']}
            require(proof.get('passed') and proof.get('reference_sha256')==reviewed['assembly_reference']['sha256']
                    and proof.get('matched_instances')==len(ids) and set(proof.get('instance_map',{}))==ids,
                    'Missing or inconsistent independent assembly evidence')
        expected = [(s.get('instance_id'), s['chain'], str(s['residue_number'])+s.get('insertion_code','').strip(), s['component'])
                    for s in reviewed['ligand_instances']]
        actual = [(s.get('instance_id'), s['chain'], s['residue'], s['component']) for s in state['ligands']]
        require(len(expected) == len(set(expected)) and len(actual) == len(set(actual))
                and set(actual) == set(expected), 'Companion inventory/identity mismatch')
        for ligand in state['ligands']:
            if 'instance_id' in ligand:
                spec = next(s for s in reviewed['ligand_instances'] if s['instance_id'] == ligand['instance_id'])
                require(all(ligand.get(k) == v for k,v in spec.items()), 'Prepared companion instance differs from reviewed config')
            names = ligand['atom_names']; positions = np.asarray(ligand['positions_nm'], dtype=float)
            require(bool(names) and len(names) == len(set(names)) and len(ligand['elements']) == len(names)
                    and positions.shape == (len(names),3) and np.isfinite(positions).all(), 'Invalid companion atom inventory')
    files = [actor['file'] for s in m['states'] for actor in s['actors']]
    require(len(files) == len(set(files)) and set(files) == set(m['data_sha256']), 'Actor file inventory mismatch')
    for name in files:
        require(Path(name).name == name, 'Prepared actor files must be directly inside the asset directory')
        require(digest(directory/name) == m['data_sha256'][name], 'Asset hash mismatch: '+name)
    return m


def check_sources(directory, sources, config_path):
    import gemmi
    from prepare_source_states import alignment_maps, selection_key
    from prepare_binding_candidates import valid_residue_id
    from source_assembly import reviewed_members
    m = assets(directory); config = json.loads(config_path.read_text())
    require(m['config'] == config and m['config_sha256'] == digest(config_path), 'Reviewed config mismatch')
    require([s['pdb'] for s in m['states']] == [s['pdb'] for s in config['states']], 'State inventory mismatch')
    basis, origin = np.array(config['basis_columns']), np.array(config['origin_angstrom'])
    require(np.allclose(basis.T@basis, np.eye(3), rtol=0, atol=1e-7) and
            np.isclose(np.linalg.det(basis), 1, rtol=0, atol=1e-7), 'Nonrigid viewing frame')
    _, mappings = alignment_maps(config)
    counts, source_ca, aligned_ca, reports = [], {}, {}, {}
    for state, reviewed in zip(m['states'], config['states']):
        path = sources/(state['pdb']+'.cif')
        require(digest(path) == reviewed['sha256'], 'Raw source hash mismatch')
        block = gemmi.cif.read(str(path)).sole_block()
        model = gemmi.make_structure_from_block(block)[0]
        polymer_table = block.get_mmcif_category('_entity_poly.')
        polymer_types = {chain: kind for chains, kind in zip(polymer_table['pdbx_strand_id'], polymer_table['type'])
                         for chain in chains.split(',')}
        entity_types = dict(zip(polymer_table['entity_id'],polymer_table['type']))
        asym = block.get_mmcif_category('_struct_asym.')
        label_entities = dict(zip(asym['id'],asym['entity_id']))
        members = {}
        if 'instance_check' in state['assembly']:
            _, selected, proof = reviewed_members(path, reviewed)
            members = {member['instance_id']: member for member in selected}
            require(proof == state['assembly']['instance_check'], 'Independent assembly check differs')
        matrix = np.array(state['alignment_matrix_angstrom'])
        require(np.allclose(matrix[3], [0, 0, 0, 1]) and
                np.allclose(matrix[:3,:3].T@matrix[:3,:3], np.eye(3), rtol=0, atol=1e-7) and
                np.isclose(np.linalg.det(matrix[:3,:3]), 1, rtol=0, atol=1e-7), 'Nonrigid state alignment')
        def local(points):
            return (np.array(points)@matrix[:3,:3].T+matrix[:3,3]-origin)@basis*.1
        require([a['chain'] for a in state['actors']] == [a['chain'] for a in reviewed['chains']], 'Polymer chain inventory mismatch')
        state_count, fit_positions = 0, {}
        for actor in state['actors']:
            if 'polymer_type' in actor:
                kind = entity_types[label_entities[actor['label_chain']]] if 'instance_id' in actor else polymer_types[actor['chain']]
                require(actor['polymer_type'] == kind, 'Raw polymer type mismatch')
            raw = []
            operation = np.eye(4)
            if 'instance_id' in actor:
                operation = np.array(members[actor['instance_id']]['source_operation_angstrom'])
                require(np.array_equal(operation, actor['source_operation_angstrom']), 'Source assembly operation differs')
            for chain in model:
                if chain.name != actor['chain']: continue
                for res in chain:
                    if res.entity_type != gemmi.EntityType.Polymer: continue
                    if 'instance_id' in actor and res.subchain != actor['label_chain']: continue
                    for atom in res:
                        if atom.element.is_hydrogen or atom.occ <= 0 or atom.altloc not in ('\x00', 'A'): continue
                        point = np.array(list(atom.pos))@operation[:3,:3].T+operation[:3,3]
                        raw.append((point, atom.name, res.seqid.num, res.seqid.icode.strip(), res.name, atom.element.name))
                        if selection_key(actor) == selection_key(mappings[state['pdb']]) and atom.name == 'CA':
                            identity = str(res.seqid) if res.seqid.icode.strip() else res.seqid.num
                            require(identity not in fit_positions and atom.element.name == 'C', 'Ambiguous raw C-alpha identity')
                            fit_positions[identity] = (point, res.name)
            d = np.load(directory/actor['file'], allow_pickle=False)
            require(len(raw) == actor['atoms'] == len(d['atoms']), 'Atom count mismatch')
            for key, col in [('atom_names',1), ('residue_numbers',2), ('insertion_codes',3), ('residue_names',4), ('elements',5)]:
                require(np.array_equal(d[key], np.array([r[col] for r in raw])), 'Raw atom identity mismatch: '+key)
            require(np.allclose(d['atoms'], local([r[0] for r in raw]), rtol=0, atol=1e-8), 'Raw atom coordinate mismatch')
            state_count += len(raw)
        source_ca[state['pdb']] = fit_positions
        aligned_ca[state['pdb']] = {identity: np.array(point)@matrix[:3,:3].T+matrix[:3,3]
                                  for identity, (point, _) in fit_positions.items()}
        reports[state['pdb']] = state['alignment_check']
        require(np.allclose(reports[state['pdb']]['matrix_angstrom'],matrix,rtol=0,atol=1e-10), 'Prepared fit matrix differs')
        required_ligands = reviewed['ligand_instances']
        require(len(state['ligands']) == len(required_ligands), 'Companion inventory mismatch')
        for ligand in state['ligands']:
            key = (ligand['chain'], ligand['residue'], ligand['component'])
            matched = [r for c in model if c.name == key[0] for r in c
                       if str(r.seqid) == key[1] and r.name == key[2]
                       and r.entity_type in (gemmi.EntityType.NonPolymer, gemmi.EntityType.Water)
                       and ('instance_id' not in ligand or r.subchain == ligand['label_chain'])]
            require(len(matched) == 1, 'Ambiguous raw companion')
            atoms = [a for a in matched[0] if not a.element.is_hydrogen and a.occ > 0 and a.altloc in ('\x00', 'A')]
            require(ligand['atom_names'] == [a.name for a in atoms] and ligand['elements'] == [a.element.name for a in atoms], 'Companion atom identity mismatch')
            points = np.array([list(a.pos) for a in atoms])
            if 'instance_id' in ligand:
                operation = np.array(members[ligand['instance_id']]['source_operation_angstrom'])
                require(np.array_equal(operation, ligand['source_operation_angstrom']), 'Companion assembly operation differs')
                points = points@operation[:3,:3].T+operation[:3,3]
            require(np.allclose(ligand['positions_nm'], local(points), rtol=0, atol=1e-8), 'Companion coordinate mismatch')
        counts.append({'state': state['pdb'], 'polymer_atoms': state_count, 'ligand_atoms': sum(len(l['atom_names']) for l in state['ligands'])})
    reference = source_ca[config['reference']]; rmsds = {}
    for code, mapping in mappings.items():
        pairs = mapping['residue_pairs']; current = source_ca[code]
        require(len(pairs) >= 3 and all(len(pair) == 2 and all(valid_residue_id(v) for v in pair) for pair in pairs)
                and len({a for a,b in pairs}) == len(pairs) == len({b for a,b in pairs})
                and all(a in current and b in reference for a,b in pairs), 'Invalid reviewed CA mapping')
        coverage_gate, rmsd_gate = mapping['min_CA_coverage'], mapping['max_CA_rmsd_angstrom']
        require(all(not isinstance(v,bool) and isinstance(v,(int,float)) and np.isfinite(v) for v in (coverage_gate,rmsd_gate))
                and 0 <= coverage_gate <= 1 and rmsd_gate >= 0, 'Invalid CA mapping gates')
        substitutions = [{'source_residue':a, 'template_residue':b, 'source_name':current[a][1], 'template_name':reference[b][1]}
                         for a,b in pairs if current[a][1] != reference[b][1]]
        require(substitutions == mapping['expected_substitutions'], 'Raw mapped residue substitutions differ')
        points = np.array([aligned_ca[code][a] for a,b in pairs]); target = np.array([reference[b][0] for a,b in pairs])
        rmsd = float(np.sqrt(np.mean(np.sum((points-target)**2, axis=1)))); coverage = len(pairs)/len(current)
        require(coverage >= coverage_gate and rmsd <= rmsd_gate, 'Common-frame alignment mismatch')
        report = reports[code]
        require(report['CA_pairs'] == pairs and report['substitutions'] == substitutions
                and np.isclose(report['source_CA_coverage'],coverage,rtol=0,atol=1e-12)
                and np.isclose(report['CA_rmsd_angstrom'],rmsd,rtol=0,atol=1e-8), 'Prepared alignment evidence differs from raw coordinates')
        rmsds[code] = rmsd
    return {'passed': True, 'source_counts': counts, 'alignment_rmsd_angstrom': rmsds,
            'scope': 'Selected raw atom identities/coordinates and shared transforms; not a trajectory or biological mechanism validation.'}


def check_scene(directory):
    import bpy
    from bpy_extras.object_utils import world_to_camera_view
    from mathutils import Vector
    m = assets(directory); scene = bpy.context.scene
    require(scene['source_states_manifest_sha256'] == digest(directory/'manifest.json'), 'Saved manifest binding mismatch')
    offsets = json.loads(scene['source_state_panel_offsets'])
    require(len(offsets) == len(m['states']), 'Panel count mismatch')
    polymer_objects = [o for o in scene.objects if 'source_state' in o]
    ligand_objects = [o for o in scene.objects if 'source_ligand_state' in o]
    require(len(polymer_objects) == sum(len(s['actors']) for s in m['states']), 'Saved polymer inventory mismatch')
    require(len(ligand_objects) == sum(len(l['atom_names']) for s in m['states'] for l in s['ligands']), 'Saved companion inventory mismatch')
    screen_bounds = []
    for state, offset in zip(m['states'], offsets):
        expected = np.eye(4); expected[:3,3] = offset
        for actor in state['actors']:
            objs = [o for o in polymer_objects if o['source_state'] == state['pdb'] and o['source_file'] == actor['file']]
            require(len(objs) == 1, 'Missing/duplicate saved polymer')
            obj = objs[0]; d = np.load(directory/actor['file'], allow_pickle=False)
            require(obj['source_author_chain'] == actor['chain'], 'Saved author-chain identity mismatch')
            for key in ('instance_id', 'label_chain', 'operator'):
                if key in actor: require(obj.get('source_'+key) == actor[key], 'Saved assembly instance identity mismatch')
            vertices = np.empty(len(obj.data.vertices)*3, np.float32)
            obj.data.vertices.foreach_get('co', vertices)
            indices = np.empty(len(obj.data.loops), np.int32)
            obj.data.loops.foreach_get('vertex_index', indices)
            require(np.array_equal(vertices.reshape(-1,3), d['vertices']) and np.array_equal(indices, d['faces'].ravel()), 'Saved mesh changed')
            require(not obj.modifiers and not obj.hide_render and not obj.hide_viewport and
                    np.allclose(obj.matrix_world, expected, rtol=0, atol=1e-6), 'Saved molecule moved, hidden or deformed')
            screen = np.array([world_to_camera_view(scene, scene.camera, obj.matrix_world@Vector(corner))
                               for corner in obj.bound_box])
            require(np.all(screen[:, :2] >= .015) and np.all(screen[:, :2] <= .985) and
                    np.all(screen[:, 2] > 0), 'Saved molecule clipped by camera')
            screen_bounds.append({'state': state['pdb'], 'chain': actor['chain'],
                                  'min_xy': screen[:, :2].min(0).tolist(), 'max_xy': screen[:, :2].max(0).tolist()})
        for j, ligand in enumerate(state['ligands']):
            for k, point in enumerate(ligand['positions_nm']):
                objs = [o for o in ligand_objects if o['source_ligand_state'] == state['pdb'] and o['ligand_index'] == j and o['atom_index'] == k]
                require(len(objs) == 1, 'Missing/duplicate saved companion atom')
                obj = objs[0]
                for key in ('instance_id', 'label_chain', 'operator', 'owner_instance'):
                    if key in ligand: require(obj.get('source_'+key) == ligand[key], 'Saved companion instance identity mismatch')
                require(not obj.hide_render and not obj.hide_viewport and not obj.modifiers and
                        np.allclose(obj.matrix_world.translation, np.array(point)+offset, rtol=0, atol=2e-6) and
                        np.allclose(np.array(obj.matrix_world)[:3,:3], np.eye(3)*.13, rtol=0, atol=1e-6), 'Saved companion transform changed')
    return {'passed': True, 'polymer_objects': len(polymer_objects), 'companion_atom_objects': len(ligand_objects),
            'projected_molecule_bounds': screen_bounds,
            'scope': 'Raw saved mesh topology, transforms and visible inventory; visual legibility checked separately.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir', type=Path, required=True)
    p.add_argument('--sources', type=Path)
    p.add_argument('--config', type=Path)
    p.add_argument('--output', type=Path, required=True)
    argv = sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else sys.argv[1:]
    a = p.parse_args(argv)
    try:
        report = check_sources(a.asset_dir, a.sources, a.config) if a.sources and a.config else check_scene(a.asset_dir)
    except (ValueError, KeyError) as error:
        a.output.write_text(json.dumps({'passed': False, 'error': str(error)}, indent=2)+'\n')
        raise
    a.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))
