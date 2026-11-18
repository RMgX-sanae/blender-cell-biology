"""Measure explicitly reviewed covalent bonds in prepared molecular display coordinates."""
import numpy as np


def _measure(rule, graph, locations, motion, sample):
    atoms = rule['atoms']
    if not isinstance(atoms, list) or len(atoms) != 2 or not all(isinstance(a, str) for a in atoms):
        raise ValueError('A reviewed bond needs two state-local atom IDs')
    edge = tuple(sorted(atoms))
    order = rule['order']
    if not isinstance(order, str) or graph.get(edge) != order.lower():
        raise ValueError(f'Bond {edge} with order {order} is absent from the selected state')
    bounds = rule['length_A']
    if (not isinstance(bounds, list) or len(bounds) != 2 or any(type(x) not in (int, float) for x in bounds)
            or not np.isfinite(bounds).all() or not 0 < bounds[0] < bounds[1]):
        raise ValueError('Supply a finite, positive bond-length interval in Angstrom')
    basis = rule['basis']
    if not isinstance(basis, str) or not basis.strip():
        raise ValueError('Identify the source or modeling basis for the bond-length interval')
    # Motion uses initial physical IDs and nm throughout, even when product atom names change.
    points = []
    physical_atoms = []
    for atom in atoms:
        component, index, original = locations[atom]
        physical_atoms.append(original)
        points.append(np.asarray(motion[component][sample, index], dtype=float))
    distance = 10 * float(np.linalg.norm(points[0] - points[1]))
    if not bounds[0] <= distance <= bounds[1]:
        raise ValueError(f'Sample {sample}, {edge} ({order}): {distance:.5f} Å outside {bounds} Å')
    return dict(atoms=atoms, physical_atoms=physical_atoms, order=order.lower(),
                length_A=round(distance, 8), interval_A=bounds, basis=basis)


def check_pose_distances(review, ledger, atom_ids, motion, frames):
    """Run after the ledger, complete initial atom set and finite nm arrays have been verified."""
    if type(review['version']) is not int or review['version'] != 1 or ledger['version'] != 2:
        raise ValueError('Reaction geometry requires a version-1 review and a version-2 reaction ledger')
    stages = ledger['stages'] if 'stages' in ledger else [ledger]
    by_state = {}
    timeline = {}
    if not review['poses']:
        raise ValueError('Select at least one pose for the bond-length check')
    for pose in review['poses']:
        state, sample = pose['state'], pose['sample']
        if type(state) is not int or not 0 <= state <= len(stages):
            raise ValueError('Select state 0 (initial) or a stage product numbered from 1')
        if type(sample) is not int or not 0 <= sample <= frames or sample in timeline:
            raise ValueError('Each reviewed pose needs a unique in-range source sample')
        if not pose['bonds']:
            raise ValueError('Select at least one bond in each reviewed pose')
        timeline[sample] = state
        by_state.setdefault(state, []).append(pose)
    states = [timeline[sample] for sample in sorted(timeline)]
    if states != sorted(states):
        raise ValueError('Reaction states run backwards in source sample time')

    locations = {atom: (component, i, atom) for component, atoms in atom_ids.items() for i, atom in enumerate(atoms)}
    result = []
    endpoints = [(stages[0], 'reactant'), *((stage, 'product') for stage in stages)]
    for state, (stage, side) in enumerate(endpoints):
        if side == 'product':
            locations = {target: locations[start] for start, target in stage['atom_map']}
        if state not in by_state:
            continue
        graph = {tuple(sorted((a, b))): order.lower() for a, b, order in stage[side + '_bonds']}
        for pose in by_state[state]:
            measured = []
            seen = set()
            for rule in pose['bonds']:
                row = _measure(rule, graph, locations, motion, pose['sample'])
                edge = tuple(sorted(row['atoms']))
                if edge in seen:
                    raise ValueError(f'Repeated bond in one reviewed pose: {edge}')
                seen.add(edge)
                measured.append(row)
            result.append(dict(state=state, sample=pose['sample'], bonds=measured))
    return sorted(result, key=lambda pose: pose['sample'])
