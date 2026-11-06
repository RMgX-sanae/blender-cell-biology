"""Generate constrained rigid binding hypotheses from two separate polymer sources.

Geometric screening only: no energy, affinity, kinetics, or blind docking claim.
Template transfer and residue-restraint initialization share the existing motion
and Blender display workflow after a strict, non-exempt endpoint check.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

from interaction_motion import contacts, plan
from prepare_molecular_asset import atom_table, surface, xyz


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rigid(matrix):
    matrix = np.asarray(matrix, float)
    if (matrix.shape != (4, 4) or not np.isfinite(matrix).all()
            or not np.allclose(matrix[3], [0, 0, 0, 1], atol=1e-9, rtol=0)
            or not np.allclose(matrix[:3, :3].T @ matrix[:3, :3], np.eye(3), atol=1e-8, rtol=0)
            or not np.isclose(np.linalg.det(matrix[:3, :3]), 1, atol=1e-8, rtol=0)):
        raise ValueError('Require a finite proper rigid transform, without reflection or scaling')
    return matrix


def transform(points, matrix):
    matrix = rigid(matrix)
    return points @ matrix[:3, :3].T + matrix[:3, 3]


def load_chain(sources, spec):
    path = sources / (spec['pdb'] + '.cif')
    meta = json.loads(path.with_suffix('.cif.source.json').read_text())
    if sha(path) != spec['sha256'] or meta['sha256'] != spec['sha256'] or meta['id'].upper() != spec['pdb']:
        raise ValueError('Requalify the source file and metadata: ' + spec['pdb'])
    rows = atom_table(path)[spec['chain']]
    if len({(r[3], r[2]) for r in rows}) != len(rows):
        raise ValueError('Ambiguous author residue/atom identities; insertion-code mapping is not supported here')
    if not spec.get('selection_note'):
        raise ValueError('Explain the construct, missing regions and exclusion of nonpolymers')
    return rows


def fit_template(first, second, mapping):
    a = {r[3]: r for r in first if r[2] == 'CA'}
    b = {r[3]: r for r in second if r[2] == 'CA'}
    pairs = mapping['residue_pairs']
    if (len(pairs) < 3 or len({p[0] for p in pairs}) != len(pairs)
            or len({p[1] for p in pairs}) != len(pairs)
            or any(p[0] not in a or p[1] not in b for p in pairs)):
        raise ValueError('Require unique, present source/template CA residue pairs')
    differences = [{'source_residue': i, 'template_residue': j, 'source_name': a[i][4], 'template_name': b[j][4]}
                   for i, j in pairs if a[i][4] != b[j][4]]
    if differences != mapping['expected_substitutions']:
        raise ValueError('Sequence substitutions differ from the reviewed mapping')
    first_points = np.array([a[i][0] for i, _ in pairs])
    second_points = np.array([b[j][0] for _, j in pairs])
    ca, cb = first_points.mean(0), second_points.mean(0)
    if min(np.linalg.matrix_rank(first_points-ca), np.linalg.matrix_rank(second_points-cb)) < 2:
        raise ValueError('Mapped CA sites must define a noncollinear fit')
    u, _, vt = np.linalg.svd((first_points-ca).T @ (second_points-cb))
    rotation = vt.T @ np.diag([1, 1, np.linalg.det(vt.T @ u.T)]) @ u.T
    matrix = np.eye(4); matrix[:3, :3] = rotation; matrix[:3, 3] = cb-rotation@ca
    rmsd = float(np.sqrt(np.mean(np.sum((transform(first_points, matrix)-second_points)**2, axis=1))))
    coverage = len(pairs)/len(a)
    if coverage < mapping['min_CA_coverage'] or rmsd > mapping['max_CA_rmsd_angstrom']:
        raise ValueError('Template mapping fails the declared coverage/RMSD gate')
    return matrix, {'CA_pairs': pairs, 'substitutions': differences, 'source_CA_coverage': coverage,
                    'CA_rmsd_angstrom': rmsd, 'matrix_angstrom': matrix.tolist()}


def restraint_indices(rows, restraints):
    result = []
    if len(restraints) < 3:
        raise ValueError('Supply at least three reviewed residue-pair restraints; these still may not define a unique pose')
    seen = set()
    for item in restraints:
        pair = (item['anchor'], item['moving'])
        key = (*pair, item.get('anchor_atom'), item.get('moving_atom'))
        if key in seen or any(type(value) is not int for value in pair):
            raise ValueError('Require distinct integer author residue pairs')
        seen.add(key)
        lo, hi = item['range_nm']
        if not np.isfinite([lo, hi]).all() or not 0 < lo < hi or not item.get('evidence'):
            raise ValueError('Each restraint requires positive distance bounds and its evidence/assumption')
        indices = [np.flatnonzero([r[3] == item[role] and (not item.get(role+'_atom') or r[2] == item[role+'_atom'])
                                  for r in rows[role]]) for role in ['anchor', 'moving']]
        if any(len(idx) == 0 for idx in indices):
            raise ValueError('Restraint selects a missing residue')
        result.append((*indices, lo, hi))
    return result


def distances(anchor, moving, restraints):
    return np.array([cdist(anchor[a], moving[b]).min() for a, b, _, _ in restraints])


def check_pose(anchor, moving, restraints, guard):
    nearest = cKDTree(anchor).query(moving)[0]
    measured = distances(anchor, moving, restraints)
    lower = np.array([r[2] for r in restraints]); upper = np.array([r[3] for r in restraints])
    violation = np.maximum(lower-measured, 0) + np.maximum(measured-upper, 0)
    overlap = int(np.count_nonzero(nearest < guard))
    return {'minimum_atom_center_nm': float(nearest.min()), 'deep_overlap_atoms': overlap,
            'restraint_distances_nm': measured.tolist(), 'max_restraint_violation_nm': float(violation.max()),
            'interface_moving_atoms_within_0.5nm': int(np.count_nonzero(nearest < .5)),
            'passed': bool(overlap == 0 and violation.max() <= 1e-4)}


def patch_frame(points, indices):
    center = points[np.unique(np.concatenate(indices))].mean(0)
    normal = center-points.mean(0)
    if np.linalg.norm(normal) < 1e-6:
        raise ValueError('An internal/central patch needs a manually reviewed initial pose')
    normal /= np.linalg.norm(normal)
    tangent = np.eye(3)[np.argmin(abs(normal))]
    tangent -= normal*np.dot(normal, tangent); tangent /= np.linalg.norm(tangent)
    return center, np.column_stack([tangent, np.cross(normal, tangent), normal])


def starts(anchor, moving, restraints, method, seed):
    if method == 'template':
        rng = np.random.default_rng(seed)
        yield np.eye(4)
        center = moving.mean(0)
        for _ in range(15):
            rotation = Rotation.from_rotvec(rng.normal(0, .12, 3)).as_matrix()
            matrix = np.eye(4); matrix[:3, :3] = rotation
            matrix[:3, 3] = center-rotation@center+rng.normal(0, .12, 3)
            yield matrix
    else:
        ca, fa = patch_frame(anchor, [r[0] for r in restraints])
        cm, fm = patch_frame(moving, [r[1] for r in restraints])
        for roll in np.linspace(0, 2*np.pi, 16, endpoint=False):
            # Oppose the two outward patch normals, then sample the remaining roll.
            rotation = fa @ Rotation.from_rotvec([0, 0, roll]).as_matrix() @ np.diag([1, -1, -1]) @ fm.T
            matrix = np.eye(4); matrix[:3, :3] = rotation
            matrix[:3, 3] = ca+fa[:, 2]*.35-rotation@cm
            yield matrix


def search(anchor, moving, restraints, method, config):
    tree = cKDTree(anchor); guard = config['atom_guard_nm']
    lower = np.array([r[2] for r in restraints]); upper = np.array([r[3] for r in restraints])
    attempts = []; accepted = []
    for index, initial in enumerate(starts(anchor, moving, restraints, method, config['seed'])):
        placed = transform(moving, initial); center = placed.mean(0)
        def delta_matrix(x):
            matrix = np.eye(4); matrix[:3, :3] = Rotation.from_rotvec(x[:3]).as_matrix()
            matrix[:3, 3] = center-matrix[:3, :3]@center+x[3:]
            return matrix
        def residual(x):
            points = transform(placed, delta_matrix(x))
            values = distances(anchor, points, restraints)
            # Padding favors strict clearance; scores remain geometry penalties, not energies.
            overlap = np.maximum(guard+.02-tree.query(points)[0], 0)
            violation = np.maximum(lower+.005-values, 0) + np.maximum(values-(upper-.005), 0)
            return np.concatenate([overlap*4, violation*4, x*.004])
        fit = least_squares(residual, np.zeros(6), max_nfev=250, diff_step=1e-4)
        matrix = rigid(delta_matrix(fit.x) @ initial); points = transform(moving, matrix)
        report = check_pose(anchor, points, restraints, guard)
        displacement = float(np.sqrt(np.mean(np.sum((points-moving)**2, axis=1))))
        if method == 'template' and displacement > config['max_template_displacement_nm']:
            report['passed'] = False
            report['rejection'] = 'Outside the declared local template-transfer neighborhood'
        record = {'start': index, 'matrix_nm': matrix.tolist(), 'screen': report,
                  'initial_frame_displacement_rms_nm': displacement,
                  'geometric_penalty': float(np.sum(fit.fun**2)), 'optimizer_converged': bool(fit.success)}
        attempts.append(record)
    for record in sorted(attempts, key=lambda r: r['geometric_penalty']):
        if not record['screen']['passed']: continue
        points = transform(moving, record['matrix_nm'])
        if any(np.sqrt(np.mean(np.sum((points-transform(moving, r['matrix_nm']))**2, axis=1))) < config['cluster_rms_nm'] for r in accepted):
            record['clustered_with_retained_candidate'] = True
            continue
        record['id'] = 'candidate_' + str(len(accepted)+1).zfill(2)
        accepted.append(record)
    return attempts, accepted


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--sources', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--method', choices=('template', 'restraints'), required=True)
    parser.add_argument('--candidate', type=int, default=1, help='1-based geometric candidate rank to animate')
    parser.add_argument('--mode', choices=('loop', 'teaching'), default='loop')
    parser.add_argument('--frames', type=int, default=144)
    parser.add_argument('--fps', type=int, default=24)
    args = parser.parse_args()
    if args.output.exists(): parser.error('Use a fresh output directory')
    if args.frames < 12 or args.frames % 2 or args.fps <= 0 or args.candidate < 1:
        parser.error('Require even frames >=12, positive fps and a 1-based candidate number')
    config = json.loads(args.config.read_text()); settings = config['search']
    for key in ['atom_guard_nm', 'max_template_displacement_nm', 'cluster_rms_nm']:
        if not isinstance(settings[key], (int, float)) or isinstance(settings[key], bool) or not np.isfinite(settings[key]) or settings[key] <= 0:
            raise ValueError('Require a positive finite search setting: ' + key)
    rows = {role: load_chain(args.sources, config[role]) for role in ['anchor', 'moving']}
    restraints = restraint_indices(rows, config['restraints'])
    points = {}; fits = {}; matrices = {}
    for role in rows:
        if args.method == 'template':
            template = config['template']
            target = load_chain(args.sources, {**template, 'chain': template[role]['chain']})
            matrix, fits[role] = fit_template(rows[role], target, template[role])
        else:
            matrix = np.eye(4); matrix[:3, 3] = -xyz(rows[role]).mean(0)
        matrices[role] = matrix
        points[role] = transform(xyz(rows[role]), matrix)*.1
    # The template origin is shared; independent monomer coordinates have no shared frame.
    origin = points['anchor'].mean(0) if args.method == 'template' else np.zeros(3)
    points = {role: p-origin for role, p in points.items()}
    direct = check_pose(points['anchor'], points['moving'], restraints, settings['atom_guard_nm']) if args.method == 'template' else None
    attempts, candidates = search(points['anchor'], points['moving'], restraints, args.method, settings)
    args.output.mkdir(parents=True)
    report = {'schema': 1, 'evidence_class': 'inferred', 'pose_kind': args.method+'_guided_candidate',
              'config': config, 'config_sha256': sha(args.config), 'method': args.method,
              'template_coordinates_used': args.method == 'template', 'template_fits': fits,
              'monomer_to_base_matrices_angstrom': {role: m.tolist() for role, m in matrices.items()},
              'shared_origin_nm': origin.tolist(), 'direct_transfer_screen': direct,
              'attempts': attempts, 'candidates': candidates,
              'limits': 'Candidates describe geometric hypotheses on one specified interface, not demonstrated alternative biological binding modes. No force field, side-chain relaxation, water, affinity, confidence probability or physical trajectory. Restraint-derived examples are not independent blind predictions.'}
    report_path = args.output/'candidates.json'
    report_path.write_text(json.dumps(report, indent=2)+'\n')
    if args.candidate > len(candidates):
        raise ValueError('Requested candidate is unavailable; preserve the failed search and revise evidence or model, not the clearance gate')
    chosen = candidates[args.candidate-1]
    points['moving'] = transform(points['moving'], chosen['matrix_nm'])
    endpoint = check_pose(points['anchor'], points['moving'], restraints, settings['atom_guard_nm'])
    if not endpoint['passed']: raise ValueError('Unqualified candidate endpoint')
    if config.get('reference_contacts'):
        diagnostic = restraint_indices(rows, config['reference_contacts'])
        endpoint['reference_contact_diagnostic'] = check_pose(points['anchor'], points['moving'], diagnostic, settings['atom_guard_nm'])
        endpoint['reference_contact_note'] = 'Diagnostic only. Failing distances remain unresolved contacts; candidate acceptance does not assert their formation.'
    poses, motion = plan(points['moving'], points['anchor'], frames=args.frames, mode=args.mode, threshold=settings['atom_guard_nm'])
    counts, minimum = contacts(points['moving'], points['anchor'], poses, settings['atom_guard_nm'])
    motion['endpoint_evidence'] = report['pose_kind']
    motion['source_contact_policy'] = 'No endpoint or path overlap exemptions for predicted candidates'
    motion['strict_overlap_atoms_max'] = int(counts.max())
    motion['collision_screen_passed'] = bool(motion['collision_screen_passed'] and counts.max() == 0)
    (args.output/'motion-screen.json').write_text(json.dumps(motion, indent=2)+'\n')
    if not motion['collision_screen_passed']: raise ValueError('Candidate passed endpoint screening but its approach path failed')
    folder = args.output/'interaction'; folder.mkdir(); actors = []
    for role, source_rows in rows.items():
        vertices, faces = surface(points[role]*10, .7, 1.1, .55)
        file = role+'.npz'
        np.savez_compressed(folder/file, atoms=points[role], vertices=vertices*.1, faces=faces,
                            atom_names=np.array([r[2] for r in source_rows]), elements=np.array([r[1] for r in source_rows]),
                            residue_numbers=np.array([r[3] for r in source_rows]), residue_names=np.array([r[4] for r in source_rows]))
        actors.append({'id': role+'_monomer', 'role': role, 'kind': 'polymer', 'file': file,
                       'chain': config[role]['chain'], 'source_pdb': config[role]['pdb'], 'atoms': len(source_rows)})
    np.savez_compressed(folder/'motion.npz', matrices=poses)
    manifest = {'input_kind': 'separate_monomers', 'example': {'id': config['id'], 'pdb': config['anchor']['pdb']+'+'+config['moving']['pdb'],
                           'anchor_label': config['anchor']['label'], 'moving_label': config['moving']['label'], 'limits': report['limits']},
                'actors': actors, 'units': 'nm', 'fps': args.fps, 'frames': args.frames, 'mode': args.mode,
                'motion': motion, 'binding_candidate': {'id': chosen['id'], 'evidence_class': 'inferred', 'pose_kind': report['pose_kind'],
                    'endpoint_screen': endpoint, 'report_file': str(report_path.resolve()), 'report_sha256': sha(report_path)},
                'display_title': 'Binding hypothesis: '+config['anchor']['label']+' / '+config['moving']['label'],
                'animation_evidence': 'Illustrative rigid motion toward a constrained binding hypothesis; neither a measured complex nor kinetics',
                'source_config': config, 'data_sha256': {name: sha(folder/name) for name in ['anchor.npz', 'moving.npz', 'motion.npz']}}
    (folder/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({'method': args.method, 'candidates': len(candidates), 'selected': chosen['id'], 'endpoint': endpoint,
                      'path_passed': motion['collision_screen_passed'], 'minimum_path_center_nm': float(minimum.min())}))


if __name__ == '__main__':
    main()
