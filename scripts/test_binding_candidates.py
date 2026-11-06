"""Audit monomer preservation and reject invalid maps, poses and evidence claims."""
import argparse
import copy
import json
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

from interaction_assets import verify_assets
from prepare_binding_candidates import fit_template, load_chain, load_companions, residue_id, restraint_indices, rigid


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources', type=Path, required=True)
    p.add_argument('--prepared', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    report = json.loads((a.prepared/'candidates.json').read_text())
    config = report['config']; folder = a.prepared/'interaction'
    manifest = json.loads((folder/'manifest.json').read_text())
    verify_assets(folder, manifest)
    rows = {role: load_chain(a.sources, config[role]) for role in ['anchor', 'moving']}
    selected = next(c for c in report['candidates'] if c['id'] == manifest['binding_candidate']['id'])
    checks = []; coordinates = {}; errors = {}
    for role in rows:
        source = np.array([r[0] for r in rows[role]])
        matrix = rigid(report['monomer_to_base_matrices_angstrom'][role])
        expected = (source@matrix[:3, :3].T+matrix[:3, 3])*.1-report['shared_origin_nm']
        if role == 'moving':
            matrix = rigid(selected['matrix_nm'])
            expected = expected@matrix[:3, :3].T+matrix[:3, 3]
        data = np.load(folder/(role+'.npz'))
        for key, index in [('atom_names', 2), ('residue_numbers', 3), ('residue_names', 4), ('elements', 1)]:
            if not np.array_equal(data[key], np.array([r[index] for r in rows[role]])):
                raise AssertionError('Source atom identity changed: '+role+' '+key)
        codes = data['insertion_codes'] if 'insertion_codes' in data else np.full(len(rows[role]), '')
        if not np.array_equal(codes, [r[5] for r in rows[role]]):
            raise AssertionError('Source insertion codes changed: '+role)
        error = float(np.abs(data['atoms']-expected).max())
        if error > 1e-10: raise AssertionError('Source coordinates changed beyond rigid placement')
        errors[role] = error; coordinates[role] = expected
    companion_count = 0
    for role in rows:
        for chain, residue, atoms in load_companions(a.sources, config[role], role):
            name = role+'_'+chain+'_'+residue.name+'_'+str(residue.seqid)
            data = np.load(folder/(name+'.npz'))
            source = np.array([list(atom.pos) for atom in atoms])
            matrix = rigid(report['monomer_to_base_matrices_angstrom'][role])
            expected = (source@matrix[:3, :3].T+matrix[:3, 3])*.1-report['shared_origin_nm']
            if role == 'moving':
                matrix = rigid(selected['matrix_nm']); expected = expected@matrix[:3, :3].T+matrix[:3, 3]
            if not np.array_equal(data['atom_names'], [atom.name for atom in atoms]) or not np.array_equal(data['elements'], [atom.element.name for atom in atoms]):
                raise AssertionError('Source companion identity changed')
            error = float(np.abs(data['atoms']-expected).max())
            if error > 1e-10: raise AssertionError('Source companion no longer follows its monomer')
            errors[name] = error; coordinates[role] = np.concatenate([coordinates[role], expected]); companion_count += len(atoms)
    # Independent dense all-pair distances: no KD tree or source-contact exemptions.
    poses = np.load(folder/'motion.npz')['matrices']; minimum = float('inf')
    for pose in poses:
        rigid(pose)
        moving = coordinates['moving']@pose[:3, :3].T+pose[:3, 3]
        minimum = min(minimum, float(cdist(coordinates['anchor'], moving).min()))
    if minimum < config['search']['atom_guard_nm']: raise AssertionError('Deep-overlap path passed preparation')
    for item in config['restraints']:
        selected_points = []
        for role in ['anchor', 'moving']:
            mask = [residue_id(r) == item[role] and (not item.get(role+'_atom') or r[2] == item[role+'_atom']) for r in rows[role]]
            selected_points.append(coordinates[role][:len(rows[role])][mask])
        distance = float(cdist(*selected_points).min()); lo, hi = item['range_nm']
        if not lo-1e-4 <= distance <= hi+1e-4: raise AssertionError('Unmet selected restraint')
    def reject(name, call):
        try: call()
        except ValueError as error: checks.append({'case': name, 'rejected': True, 'diagnostic': str(error)})
        else: raise AssertionError('Accepted '+name)
    for name, value in [('reflection', -1.), ('stretch', 1.1), ('nonfinite', float('nan'))]:
        matrix = np.eye(4); matrix[0, 0] = value
        reject(name, lambda matrix=matrix: rigid(matrix))
    changed = copy.deepcopy(config['restraints']); changed[0]['moving'] = 9999
    reject('missing-residue', lambda: restraint_indices(rows, changed))
    changed = copy.deepcopy(config['restraints']); changed[0]['anchor_atom'] = 'NOT_AN_ATOM'
    reject('missing-atom', lambda: restraint_indices(rows, changed))
    changed = copy.deepcopy(config['restraints']); changed[0]['range_nm'] = [float('nan'), 1.]
    reject('nonfinite-restraint', lambda: restraint_indices(rows, changed))
    changed = copy.deepcopy(config['restraints']); changed[0]['evidence'] = ''
    reject('unattributed-restraint', lambda: restraint_indices(rows, changed))
    changed = copy.deepcopy(config['restraints']); changed[1] = changed[0]
    reject('duplicate-restraint', lambda: restraint_indices(rows, changed))
    if report['method'] == 'template':
        template = config['template']
        target = load_chain(a.sources, {**template, 'chain': template['moving']['chain']})
        mapping = copy.deepcopy(template['moving']); mapping['expected_substitutions'].append({'unreviewed': True})
        reject('silently-matched-mutant', lambda: fit_template(rows['moving'], target, mapping))
        mapping = copy.deepcopy(template['moving']); mapping['residue_pairs'][1] = mapping['residue_pairs'][0]
        reject('duplicated-fit-identity', lambda: fit_template(rows['moving'], target, mapping))
        mapping = copy.deepcopy(template['moving']); mapping['residue_pairs'] = mapping['residue_pairs'][:-1]; mapping['min_CA_coverage'] = 1.
        reject('insufficient-declared-coverage', lambda: fit_template(rows['moving'], target, mapping))
    for role in rows:
        if not config[role].get('ligand_instances'): continue
        spec = copy.deepcopy(config[role]); spec['ligand_instances'][0]['residue_number'] = 9999
        reject('missing-companion-'+role, lambda: load_companions(a.sources, spec, role))
        spec = copy.deepcopy(config[role]); spec['ligand_instances'][0]['owner_role'] = 'moving' if role == 'anchor' else 'anchor'
        reject('wrong-companion-owner-'+role, lambda: load_companions(a.sources, spec, role))
        spec = copy.deepcopy(config[role]); spec['ligand_instances'].append(spec['ligand_instances'][0])
        reject('duplicate-companion-'+role, lambda: load_companions(a.sources, spec, role))
    if companion_count:
        changed = copy.deepcopy(manifest); removed = next(actor for actor in changed['actors'] if actor['kind']=='ligand')
        changed['actors'].remove(removed); del changed['data_sha256'][removed['file']]
        reject('companion-silently-omitted', lambda: verify_assets(folder, changed))
    for role in rows:
        inserted = next((r for r in rows[role] if r[5]), None)
        if inserted:
            item = copy.deepcopy(config['restraints'][0]); item[role] = str(inserted[3])+('B' if inserted[5] != 'B' else 'C')
            changed = [item, *copy.deepcopy(config['restraints'][1:])]
            reject('wrong-insertion-code-'+role, lambda: restraint_indices(rows, changed))
    for name, mutate in [
        ('candidate-as-measured', lambda m: m['binding_candidate'].update(evidence_class='measured')),
        ('overlapping-endpoint', lambda m: m['binding_candidate']['endpoint_screen'].update(deep_overlap_atoms=1)),
        ('exempted-predicted-contact', lambda m: m['motion'].update(endpoint_source_overlap_pairs=1)),
        ('wrong-monomer-source', lambda m: m['actors'][0].update(source_pdb='1BRS')),
        ('candidate-label-removed', lambda m: m.update(display_title='Observed binding')),
        ('path-has-overlap', lambda m: m['motion'].update(strict_overlap_atoms_max=1)),
    ]:
        changed = copy.deepcopy(manifest); mutate(changed)
        reject(name, lambda changed=changed: verify_assets(folder, changed))
    result = {'passed': True, 'method': report['method'], 'source_atom_counts': {role: len(r) for role, r in rows.items()},
              'source_companion_atoms': companion_count,
              'source_coordinate_error_max_nm': errors, 'independent_dense_path_minimum_nm': minimum,
              'checked_poses': len(poses), 'negative_cases': len(checks), 'checks': checks,
              'limits': 'Geometric/source identity audit, not a validation of biological binding or affinity.'}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['passed', 'method', 'source_atom_counts', 'checked_poses', 'negative_cases']}))


if __name__ == '__main__':
    main()
