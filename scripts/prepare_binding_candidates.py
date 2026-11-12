"""Generate constrained rigid binding hypotheses from two separate polymer sources.

Geometric screening only: no energy, affinity, kinetics, or blind docking claim.
Template transfer and residue-restraint initialization share the existing motion
and Blender display workflow after a strict, non-exempt endpoint check.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re

import gemmi
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

from interaction_motion import contacts, plan
from prepare_molecular_asset import atom_table, surface, xyz
from prepare_interaction import select_ligands
from source_assembly import reviewed_members


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


def checked_source(sources, spec):
    path = sources / (spec['pdb'] + '.cif')
    meta = json.loads(path.with_suffix('.cif.source.json').read_text())
    if sha(path) != spec['sha256'] or meta['sha256'] != spec['sha256'] or meta['id'].upper() != spec['pdb']:
        raise ValueError('Requalify the source file and metadata: ' + spec['pdb'])
    return path


def residue_id(row):
    return str(row[3])+row[5] if len(row) > 5 and row[5] else row[3]


def valid_residue_id(value):
    return type(value) is int or isinstance(value, str) and re.fullmatch(r'-?\d+[A-Za-z]', value) is not None


def load_chain(sources, spec):
    rows = atom_table(checked_source(sources, spec), include_insertion_codes=True)[spec['chain']]
    if len({(residue_id(r), r[2]) for r in rows}) != len(rows):
        raise ValueError('Ambiguous author residue/insertion-code/atom identities')
    if not spec.get('selection_note'):
        raise ValueError('Explain the construct, missing regions and exclusion of nonpolymers')
    return rows


def load_companions(sources, spec, role):
    if 'ligand_instances' not in spec: return []
    selected = select_ligands(gemmi.read_structure(str(checked_source(sources, spec)))[0], spec)
    if any(entries for owner, entries in selected.items() if owner != role):
        raise ValueError('A monomer companion must belong to its declared source partner')
    for chain, residue, atoms in selected[role]:
        name = role+'_'+chain+'_'+residue.name+'_'+str(residue.seqid)
        if not re.fullmatch(r'[A-Za-z0-9_-]+', name):
            raise ValueError('Require a filename-safe companion identity')
    return selected[role]


def load_partner(sources, spec, role):
    """Keep complete polymer instances separate, with one shared partner frame."""
    if 'chains' in spec:
        if 'chain' in spec or not spec.get('selection_note'):
            raise ValueError('Declare reviewed assembly instances instead of an ambiguous single chain')
        _, members, proof = reviewed_members(checked_source(sources, spec), spec)
        result = []
        for member in members:
            if member.get('owner_role', role) != role: raise ValueError('Companion belongs to another source partner')
            name = role+'_'+member['instance_id']
            identity = {k:v for k,v in member.items() if k not in ('source_rows','kind','component')}
            item = {**identity, 'id': name, 'file': name+'.npz', 'role': role,
                    'source_pdb': spec['pdb'], 'assembly_id': spec['assembly_id'],
                    'kind': 'polymer' if member['kind'] == 'polymer' else 'ligand'}
            if item['kind'] == 'ligand': item['ccd'] = member['component']
            result.append(item)
        return result, proof
    result = [{'id': role+'_monomer', 'role': role, 'kind': 'polymer', 'file': role+'.npz',
               'chain': spec['chain'], 'source_pdb': spec['pdb'], 'rows': load_chain(sources, spec)}]
    for chain, residue, atoms in load_companions(sources, spec, role):
        name = role+'_'+chain+'_'+residue.name+'_'+str(residue.seqid)
        result.append({'id': name, 'role': role, 'kind': 'ligand', 'file': name+'.npz', 'chain': chain,
            'source_pdb': spec['pdb'], 'ccd': residue.name, 'residue_number': residue.seqid.num,
            'insertion_code': residue.seqid.icode.strip(),
            'rows': [(list(a.pos),a.element.name,a.name,residue.seqid.num,residue.name,residue.seqid.icode.strip()) for a in atoms]})
    return result, None


def rigid_fit_points(first_points, second_points):
    ca, cb = first_points.mean(0), second_points.mean(0)
    if min(np.linalg.matrix_rank(first_points-ca), np.linalg.matrix_rank(second_points-cb)) < 2:
        raise ValueError('Mapped CA sites must define a noncollinear fit')
    u, _, vt = np.linalg.svd((first_points-ca).T @ (second_points-cb))
    rotation = vt.T @ np.diag([1, 1, np.linalg.det(vt.T @ u.T)]) @ u.T
    matrix = np.eye(4); matrix[:3, :3] = rotation; matrix[:3, 3] = cb-rotation@ca
    rmsd = float(np.sqrt(np.mean(np.sum((transform(first_points, matrix)-second_points)**2, axis=1))))
    return matrix, rmsd


def fit_template(first, second, mapping):
    def ca_sites(rows):
        sites = {}
        for row in rows:
            if row[2] != 'CA': continue
            key = residue_id(row); point = np.asarray(row[0], dtype=float)
            if (not valid_residue_id(key) or key in sites or row[1] != 'C'
                    or point.shape != (3,) or not np.isfinite(point).all()):
                raise ValueError('Require unique, finite carbon-alpha identities in each source chain')
            sites[key] = row
        return sites
    a, b = ca_sites(first), ca_sites(second)
    coverage_gate, rmsd_gate = mapping['min_CA_coverage'], mapping['max_CA_rmsd_angstrom']
    if (any(isinstance(v, bool) or not isinstance(v, (int, float)) or not np.isfinite(v)
            for v in (coverage_gate, rmsd_gate)) or not 0 <= coverage_gate <= 1 or rmsd_gate < 0):
        raise ValueError('Require finite CA coverage in [0,1] and a finite nonnegative RMSD gate')
    pairs = mapping['residue_pairs']
    if any(not isinstance(p, (list, tuple)) or len(p) != 2 or
           any(not valid_residue_id(value) for value in p) for p in pairs):
        raise ValueError('CA pairs require integer author numbers or explicit insertion-code identities')
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
    matrix, rmsd = rigid_fit_points(first_points, second_points)
    coverage = len(pairs)/len(a)
    if coverage < coverage_gate or rmsd > rmsd_gate:
        raise ValueError('Template mapping fails the declared coverage/RMSD gate')
    return matrix, {'CA_pairs': pairs, 'substitutions': differences, 'source_CA_coverage': coverage,
                    'CA_rmsd_angstrom': rmsd, 'matrix_angstrom': matrix.tolist()}


def fit_partner_template(first, second, mappings):
    source = {m['instance_id']: m['rows'] for m in first if m['kind'] == 'polymer'}
    target = {m['instance_id']: m['rows'] for m in second if m['kind'] == 'polymer'}
    if (len(mappings) != len(source) or {m['instance_id'] for m in mappings} != set(source)
            or len(mappings) != len(target) or {m['template_instance_id'] for m in mappings} != set(target)):
        raise ValueError('Map every retained partner/template polymer instance exactly once')
    points, reports = [], []
    for mapping in mappings:
        a, b = source[mapping['instance_id']], target[mapping['template_instance_id']]
        _, report = fit_template(a, b, mapping)
        sites = [{residue_id(r): r[0] for r in rows if r[2] == 'CA'} for rows in (a,b)]
        pair = [np.array([sites[k][p[k]] for p in mapping['residue_pairs']]) for k in (0,1)]
        points.append(pair);reports.append(report)
    matrix, rmsd = rigid_fit_points(*(np.concatenate([p[k] for p in points]) for k in (0,1)))
    for mapping, pair, report in zip(mappings,points,reports):
        error = float(np.sqrt(np.mean(np.sum((transform(pair[0],matrix)-pair[1])**2,axis=1))))
        if error > mapping['max_CA_rmsd_angstrom']:
            raise ValueError('One shared partner transform fails a member CA RMSD gate')
        report.update(instance_id=mapping['instance_id'],template_instance_id=mapping['template_instance_id'],
                      CA_rmsd_angstrom=error,matrix_angstrom=matrix.tolist())
    return matrix, {'matrix_angstrom':matrix.tolist(),'CA_rmsd_angstrom':rmsd,'members':reports}


def restraint_indices(rows, restraints, instance_ids=None):
    result = []
    if len(restraints) < 3:
        raise ValueError('Supply at least three reviewed residue-pair restraints; these still may not define a unique pose')
    seen = set()
    for item in restraints:
        pair = (item['anchor'], item['moving'])
        key = (*pair, item.get('anchor_atom'), item.get('moving_atom'),item.get('anchor_instance'),item.get('moving_instance'))
        if key in seen or any(not valid_residue_id(value) for value in pair):
            raise ValueError('Require distinct author residue pairs, with explicit insertion codes when present')
        seen.add(key)
        lo, hi = item['range_nm']
        if not np.isfinite([lo, hi]).all() or not 0 < lo < hi or not item.get('evidence'):
            raise ValueError('Each restraint requires positive distance bounds and its evidence/assumption')
        indices = []
        for role in ['anchor','moving']:
            ids = instance_ids.get(role) if instance_ids else None
            member = item.get(role+'_instance')
            if (ids is not None and (len(ids) != len(rows[role]) or member not in ids)) or (ids is None and member is not None):
                raise ValueError('Restraint must identify a retained polymer instance in its declared partner')
            indices.append(np.flatnonzero([residue_id(r) == item[role]
                and (ids is None or ids[i] == member) and (not item.get(role+'_atom') or r[2] == item[role+'_atom'])
                for i,r in enumerate(rows[role])]))
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
    loaded = {role: load_partner(args.sources, config[role], role) for role in ['anchor', 'moving']}
    partners = {role: value[0] for role,value in loaded.items()}
    rows = {role: [row for member in members if member['kind'] == 'polymer' for row in member['rows']]
            for role,members in partners.items()}
    instance_ids = {role: [m['instance_id'] for m in members if m['kind'] == 'polymer' for _ in m['rows']]
                    if 'chains' in config[role] else None for role,members in partners.items()}
    restraints = restraint_indices(rows, config['restraints'], instance_ids)
    points = {}; fits = {}; matrices = {}
    for role in rows:
        if args.method == 'template':
            template = config['template']
            if 'chains' in config[role]:
                target, proof = load_partner(args.sources, {**template, **template[role]}, role)
                matrix, fits[role] = fit_partner_template(partners[role], target, template[role]['maps'])
                fits[role]['template_assembly_check'] = proof
            else:
                target = load_chain(args.sources, {**template, 'chain': template[role]['chain']})
                matrix, fits[role] = fit_template(rows[role], target, template[role])
        else:
            matrix = np.eye(4); matrix[:3, 3] = -xyz(rows[role]).mean(0)
        matrices[role] = matrix
        points[role] = transform(xyz(rows[role]), matrix)*.1
    # The template origin is shared; independent monomer coordinates have no shared frame.
    origin = points['anchor'].mean(0) if args.method == 'template' else np.zeros(3)
    points = {role: p-origin for role, p in points.items()}
    companion_points = {role: [transform(xyz(m['rows']), matrices[role])*.1-origin
                              for m in members if m['kind'] == 'ligand'] for role, members in partners.items()}
    # Polymer atoms stay first so residue restraint indices remain unchanged.
    points = {role: np.concatenate([p, *companion_points[role]]) for role, p in points.items()}
    direct = check_pose(points['anchor'], points['moving'], restraints, settings['atom_guard_nm']) if args.method == 'template' else None
    attempts, candidates = search(points['anchor'], points['moving'], restraints, args.method, settings)
    args.output.mkdir(parents=True)
    report = {'schema': 1, 'evidence_class': 'inferred', 'pose_kind': args.method+'_guided_candidate',
              'config': config, 'config_sha256': sha(args.config), 'method': args.method,
              'template_coordinates_used': args.method == 'template', 'template_fits': fits,
              'monomer_to_base_matrices_angstrom': {role: m.tolist() for role, m in matrices.items()},
              'shared_origin_nm': origin.tolist(), 'direct_transfer_screen': direct,
              'attempts': attempts, 'candidates': candidates,
              'limits': 'Candidates describe geometric hypotheses on one specified interface, not demonstrated alternative biological binding modes. Explicit source-bound nonpolymers follow their owner rigidly and enter cross-partner screening. No force field, side-chain relaxation, bulk solvent, affinity, confidence probability or physical trajectory. Restraint-derived examples are not independent blind predictions.'}
    if any(proof for _,proof in loaded.values()):
        report['partner_assembly_checks'] = {role: proof for role,(_,proof) in loaded.items() if proof}
    report_path = args.output/'candidates.json'
    report_path.write_text(json.dumps(report, indent=2)+'\n')
    if args.candidate > len(candidates):
        raise ValueError('Requested candidate is unavailable; preserve the failed search and revise evidence or model, not the clearance gate')
    chosen = candidates[args.candidate-1]
    points['moving'] = transform(points['moving'], chosen['matrix_nm'])
    endpoint = check_pose(points['anchor'], points['moving'], restraints, settings['atom_guard_nm'])
    if not endpoint['passed']: raise ValueError('Unqualified candidate endpoint')
    if config.get('reference_contacts'):
        diagnostic = restraint_indices(rows, config['reference_contacts'], instance_ids)
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
    for role, members in partners.items():
        offset = 0
        for member in members:
            source_rows = member['rows']; atoms = points[role][offset:offset+len(source_rows)]
            data = {'atoms':atoms,'elements':np.array([r[1] for r in source_rows]),
                    'atom_names':np.array([r[2] for r in source_rows])}
            actor = {k:v for k,v in member.items() if k != 'rows'}
            if member['kind'] == 'polymer':
                vertices, faces = surface(atoms*10, .7, 1.1, .55)
                data.update(vertices=vertices*.1,faces=faces,residue_numbers=np.array([r[3] for r in source_rows]),
                            residue_names=np.array([r[4] for r in source_rows]))
                if any(r[5] for r in source_rows): data['insertion_codes'] = np.array([r[5] for r in source_rows])
            else: actor['atom_names'] = [r[2] for r in source_rows]
            np.savez_compressed(folder/member['file'], **data)
            actors.append({**actor,'atoms':len(source_rows)})
            offset += len(source_rows)
    np.savez_compressed(folder/'motion.npz', matrices=poses)
    manifest = {'input_kind': 'separate_monomers', 'example': {'id': config['id'], 'pdb': config['anchor']['pdb']+'+'+config['moving']['pdb'],
                           'anchor_label': config['anchor']['label'], 'moving_label': config['moving']['label'], 'limits': report['limits']},
                'actors': actors, 'units': 'nm', 'fps': args.fps, 'frames': args.frames, 'mode': args.mode,
                'motion': motion, 'binding_candidate': {'id': chosen['id'], 'evidence_class': 'inferred', 'pose_kind': report['pose_kind'],
                    'endpoint_screen': endpoint, 'report_file': str(report_path.resolve()), 'report_sha256': sha(report_path)},
                'display_title': 'Binding hypothesis: '+config['anchor']['label']+' / '+config['moving']['label'],
                'animation_evidence': 'Illustrative rigid motion toward a constrained binding hypothesis; neither a measured complex nor kinetics',
                'source_config': config, 'data_sha256': {name: sha(folder/name) for name in [*(a['file'] for a in actors), 'motion.npz']}}
    if report.get('partner_assembly_checks'):
        manifest['input_kind'] = 'separate_partners'
        manifest['partner_assembly_checks'] = report['partner_assembly_checks']
    (folder/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({'method': args.method, 'candidates': len(candidates), 'selected': chosen['id'], 'endpoint': endpoint,
                      'path_passed': motion['collision_screen_passed'], 'minimum_path_center_nm': float(minimum.min())}))


if __name__ == '__main__':
    main()
