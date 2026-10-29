"""Rebuild a source-based RNF4/Ubc13/UEV K63 net-transfer teaching sequence.

Starts with precharged E2 and ends at the bound product by default. --release
appends a qualified inferred product departure. No E1 charging, elementary
proton relay, proteasome targeting, or biochemical time inference is implemented.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

from interaction_motion import matrices
from prepare_molecular_asset import surface
from ubiquitin_geometry import UbiquitinGeometry, PROTEINS, SOURCE_HASHES, digest
from reaction_ledger import validate as validate_reaction_ledger


def smooth(value):
    q = np.clip(value, 0, 1)
    return q*q*q*(10+q*(-15+6*q))


def sample(model, config, times):
    curve = PchipInterpolator(config['torsion_times'], config['torsion_parameters'])
    reactant = model.pose(curve(0)); start, end = config['transfer_interval']
    arrays = {k: [] for k in model.rows}; states = []; fractions = []
    time_scale = config['release']['transfer_fraction'] if config.get('include_release') else 1.
    for time in times:
        t = min(time/time_scale, 1.)
        fraction = float(smooth((t-start)/(end-start))); points = model.pose(curve(fraction))
        for spec in config['paths'].values():
            begin, finish = spec['interval']; q = np.clip((t-begin)/(finish-begin), 0, 1)
            moving = np.concatenate([reactant[k] for k in spec['moving']])*.1
            pose = matrices(moving, spec['direction'], spec['distance_nm'], np.array([q/2]), 'loop', spec['turn_radians'], spec['bend_fraction'])[0]
            for k in spec['moving']: points[k] = points[k]@pose[:3, :3].T+pose[:3, 3]*10
        if config.get('include_release') and time >= time_scale:
            points = model.release_pose(points, config['release'], time)
        for k in arrays: arrays[k].append(points[k])
        fractions.append(fraction); states.append('reactant' if fraction == 0 else 'product' if fraction == 1 else 'transition')
    return {k: np.array(v) for k, v in arrays.items()}, states, np.array(fractions)


def audit(model, arrays, states, fractions):
    cross_minimum = np.inf; self_minimum = np.inf; failures = []; bond_error = 0.; endpoint = {}
    bonds = {}
    for k, edges in model.bonds.items():
        reference = np.linalg.norm(model.source[k][edges[:, 0]]-model.source[k][edges[:, 1]], axis=1)
        observed = np.linalg.norm(arrays[k][:, edges[:, 0]]-arrays[k][:, edges[:, 1]], axis=2)
        error = float(np.max(abs(observed-reference))); bond_error = max(bond_error, error)
        bonds[k] = {'edges': len(edges), 'maximum_length_change_A': error}
    for frame, state in enumerate(states):
        points = {k: v[frame] for k, v in arrays.items()}; check = model.check(points, state)
        cross = min(check['nonbonded_minima_A'].values()); internal = min(check['mobile_self_nonbonded_minima_A'].values())
        cross_minimum = min(cross_minimum, cross); self_minimum = min(self_minimum, internal)
        if cross < 2.38 or internal < 2.38: failures.append({'frame': frame+1, 'cross_min_A': cross, 'mobile_self_min_A': internal})
        if state == 'reactant' and abs(check['C_S_A']-1.82) > .03: failures.append({'frame': frame+1, 'thioester_length_A': check['C_S_A']})
        if state == 'product' and (abs(check['C_N_A']-1.33) > .025 or check['C_S_A'] < 3): failures.append({'frame': frame+1, 'product_geometry': check})
        angles = check['reaction_angles']
        if state == 'reactant' and (abs(angles['O_C_S_degrees']-123) > 10 or abs(angles['CB_S_C_degrees']-100) > 10 or angles['thioester_carbonyl_plane_error_degrees'] > 10):
            failures.append({'frame': frame+1, 'thioester_angles': angles})
        if state == 'product' and (abs(angles['O_C_N_degrees']-123) > 10 or abs(angles['C_N_CE_degrees']-122) > 10 or abs(angles['O_C_N_CE_degrees']) > 10 or abs(abs(angles['CA_C_N_CE_degrees'])-180) > 10):
            failures.append({'frame': frame+1, 'amide_angles': angles})
        if state != 'transition': endpoint[state] = check
    # Whole-chain torsions and rigid transforms preserve adjacent bond angles;
    # independently check the corresponding two-bond endpoint distances too.
    angle_error = 0.
    for k, neighbors in model.neighbors.items():
        pairs = sorted({tuple(sorted([a, b])) for group in neighbors for a in group for b in group if a != b})
        pairs = np.array(pairs)
        ref = np.linalg.norm(model.source[k][pairs[:, 0]]-model.source[k][pairs[:, 1]], axis=1)
        distances = np.linalg.norm(arrays[k][:, pairs[:, 0]]-arrays[k][:, pairs[:, 1]], axis=2)
        angle_error = max(angle_error, float(np.max(abs(distances-ref))))
    zinc_relative_error = float(np.max(abs(arrays['zinc']-model.source['zinc'])))
    passed = not failures and bond_error < 1e-6 and angle_error < 1e-6 and zinc_relative_error < 1e-6
    release = model.release_checks(arrays, states) if hasattr(model, 'release_checks') else None
    if release is not None: passed = passed and release['passed']
    return {'passed': passed, 'frames_checked': len(states), 'threshold_A': 2.38,
            'minimum_cross_nonbonded_A': cross_minimum, 'minimum_mobile_self_nonbonded_A': self_minimum,
            'covalent_geometry': bonds, 'max_adjacent_bond_pair_distance_change_A': angle_error,
            'structural_zinc_coordinate_error_A': zinc_relative_error, 'endpoints': endpoint, 'failures': failures, 'release': release,
            'contact_policy': 'Check all protein pairs and non-RING contacts with all four Zn. Exclude only graph-derived 1-2 and 1-3 pairs for the declared reaction bonds; transition allows both endpoint bond neighborhoods. Static RING/Zn coordination remains unchanged.',
            'limits': 'Heavy-atom center overlap screen and geometric invariants; not a force-field calculation, energy minimization, molecular dynamics, proton pathway, or observed trajectory.'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--config', type=Path, default=Path(__file__).resolve().parents[1]/'assets/ubiquitin-transfer.json')
    p.add_argument('--frames', type=int); p.add_argument('--fps', type=int, default=24)
    p.add_argument('--release', action='store_true', help='Append inferred E2/UEV departure and intact diubiquitin release; default 36 seconds at 24 fps')
    a = p.parse_args()
    if a.frames is None: a.frames = 864 if a.release else 576
    if a.output.exists(): p.error('Use a fresh output directory')
    if a.frames < 96 or a.fps < 1: p.error('Use at least 96 frames and a positive frame rate')
    config = json.loads(a.config.read_text()); config['include_release'] = a.release
    if a.release:
        from ubiquitin_release import UbiquitinReleaseGeometry
        model = UbiquitinReleaseGeometry(a.sources)
    else: model = UbiquitinGeometry(a.sources)
    selected_atoms = [[f'{k}:{r[3]}:{r[2]}', r[1]] for k, rows in model.rows.items() for r in rows]
    ledger = {'version': 1,
        'scope': 'All selected source heavy atoms, including four Zn. The same source indices persist in both endpoint poses; hydrogens are implicit.',
        'bond_scope': 'Net-transfer reactive bonds and donor terminal C-O adjacency only; bond order and Zn coordination are not represented.',
        'reactant_atoms': selected_atoms, 'product_atoms': selected_atoms,
        'atom_map': [[identity, identity] for identity, _ in selected_atoms],
        'reactant_bonds': [['e2:87:SG', 'donor:76:C'], ['donor:76:C', 'donor:76:O']],
        'product_bonds': [['donor:76:C', 'acceptor:63:NZ'], ['donor:76:C', 'donor:76:O']],
        'break': [['e2:87:SG', 'donor:76:C']], 'form': [['donor:76:C', 'acceptor:63:NZ']]}
    ledger_check = validate_reaction_ledger(ledger)
    arrays, states, fractions = sample(model, config, np.arange(a.frames+1)/a.frames)
    checks = audit(model, arrays, states, fractions)
    if not checks['passed']:
        print(json.dumps(checks, indent=2)); raise ValueError('Fix motion geometry before writing a renderable asset')
    a.output.mkdir(parents=True)
    (a.output/'reaction-ledger.json').write_text(json.dumps(ledger, indent=2)+'\n')
    origin = np.concatenate([model.source[k] for k in PROTEINS]).mean(0)
    np.savez_compressed(a.output/'motion.npz', **{k: (v-origin)*.1 for k, v in arrays.items()}, transfer_fraction=fractions)
    identities = {}
    for k, rows in model.rows.items():
        ids = [f'{k}:{r[3]}:{r[2]}' for r in rows]
        if len(set(ids)) != len(ids): raise ValueError('Repeated atom identity: '+k)
        identities[k] = ids
        np.savez_compressed(a.output/(k+'.npz'), atoms=(model.source[k]-origin)*.1,
            bonds=model.bonds.get(k, np.empty((0, 2), dtype=int)), atom_names=[r[2] for r in rows],
            residue_numbers=[r[3] for r in rows], elements=[r[1] for r in rows], residue_names=[r[4] for r in rows])
    # Background shapes are complete source proteins from this reconstituted
    # molecular example. Their placements are artistic, not cellular abundance.
    for k in ('donor', 'e2', 'uev'):
        points = model.source[k]-model.source[k].mean(0)
        vertices, faces = surface(points, .8, 1.15, .55)
        np.savez_compressed(a.output/('background_'+k+'.npz'), vertices=vertices*.1, faces=faces)
    count = sum(map(len, model.rows.values())); assert count == 4476
    scale = config['release']['transfer_fraction'] if a.release else 1.
    captions = [[begin*scale, end*scale, body] for begin, end, body in config['captions']]
    if a.release: captions += config['release']['captions']
    manifest = {'kind': config['kind'], 'passed': True, 'frames': a.frames, 'fps': a.fps,
        'reaction_ledger': ledger_check,
        'source_directory': str(a.sources.resolve()), 'source_sha256': SOURCE_HASHES, 'source_changes': model.source_record,
        'origin_A': origin.tolist(), 'config': config, 'config_sha256': digest(a.config), 'checks': checks,
        'atom_ids': identities, 'atom_counts': {k: len(v) for k, v in model.rows.items()}, 'reaction_atom_indices': model.ids,
        'heavy_atoms_initial': count, 'heavy_atoms_final': count, 'conserved_structural_zinc': 4,
        'release_included': a.release, 'captions': captions,
        'timing': {'acceptor_arrival': config['paths']['acceptor_in']['interval'][1]*scale, 'product_formed': config['transfer_interval'][1]*scale},
        'reaction': {'break': ['e2:87:SG', 'donor:76:C'], 'form': ['donor:76:C', 'acceptor:63:NZ'],
            'donor_terminal_OXT': False, 'acceptor_terminal_OXT': True, 'heavy_atom_byproducts': [],
            'hydrogens': 'Implicit. No definitive general-base or proton-transfer sequence is asserted.'},
        'scope': 'Precharged E2/RING/UEV/acceptor assembly and inferred net K63 ubiquitin transfer, '+('followed by E2/UEV departure and intact diubiquitin release.' if a.release else 'ending at bound diubiquitin.'),
        'limits': ['5AIT is a 3.4 A engineered crosslink construct; native restoration and catalytic poses are inferred.',
            'Whole-source proteins retain resolved residues only: Ub 2-76, Ubc13 4-150, UEV 6-145 and linked RING 131-259. Other unresolved termini/atoms are not reconstructed.',
            'Rat linked RNF4 RING dimer, human Ubc13/UEV, bovine ubiquitin: a mixed experimental construct, not one native cellular assembly.',
            'Relative recruitment/release order and all timing are illustrative; movie duration is not a biological turnover time.',
            'Source support does not resolve elementary reaction intermediates or the lysine-deprotonation route.',
            'No E1 charging, ubiquitin recycling, proteasome pathway, or general ubiquitination solver is implemented here.',
            'Release uses inferred interface side-chain torsions and rigid departures, keeping the E2/UEV pair and product bond intact.' if a.release else 'This selected sequence ends at bound product; use --release to append the qualified release example.'],
        'data_sha256': {f.name: digest(f) for f in sorted(a.output.glob('*.npz'))}}
    (a.output/'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'passed': True, 'frames': a.frames, 'atoms': count, 'cross_min_A': checks['minimum_cross_nonbonded_A'], 'output': str(a.output)}))


if __name__ == '__main__': main()
