"""Measure local covalent geometry at source samples and saved Blender poses."""
from typing import NamedTuple

import numpy as np


class Measure(NamedTuple):
    size: int
    field: str
    interval: str
    unit: str
    saved_field: str


MEASURES = {
    'bonds': Measure(2, 'length_A', 'interval_A', 'Å', 'saved_bond_lengths'),
    'angles': Measure(3, 'angle_deg', 'interval_deg', '°', 'saved_valence_angles'),
    'signed_heights': Measure(4, 'height_A', 'interval_A', 'Å', 'saved_signed_heights'),
}


def signed_height(points):
    """Height of the fourth point above the ordered first three, in the input length unit."""
    origin, a, b, site = points
    u, v = a - origin, b - origin
    normal = np.cross(u, v)
    length = np.linalg.norm(normal)
    if length <= 1e-6 * np.linalg.norm(u) * np.linalg.norm(v):
        raise ValueError("Cannot orient a plane from coincident or collinear anchors")
    return float((site - origin) @ (normal / length))


def measure_nm(kind, points):
    points = np.asarray(points, dtype=float)
    if kind == 'bonds':
        return 10 * float(np.linalg.norm(points[0] - points[1]))
    if kind == 'angles':
        u, v = points[0] - points[1], points[2] - points[1]
        if not np.linalg.norm(u) or not np.linalg.norm(v):
            raise ValueError('Cannot measure a valence angle with coincident sites')
        return float(np.degrees(np.arctan2(np.linalg.norm(np.cross(u, v)), u @ v)))
    if kind == 'signed_heights':
        return 10 * signed_height(points)
    raise ValueError(f'Unknown pose measurement: {kind}')


def _measure(kind, rule, graph, locations, motion, sample):
    measure = MEASURES[kind]
    atoms = rule['atoms']
    if (not isinstance(atoms, list) or len(atoms) != measure.size or not all(isinstance(a, str) for a in atoms)
            or len(set(atoms)) != measure.size):
        raise ValueError(f'{kind} requires {measure.size} distinct state-local atom IDs')
    if kind == 'bonds':
        identity = tuple(sorted(atoms))
        order = rule['order']
        if not isinstance(order, str) or graph.get(identity) != order.lower():
            raise ValueError(f'Bond {identity} with order {order} is absent from the selected state')
        chemistry = {'order': order.lower()}
    else:
        # Angle: [arm, center, arm]. Height: [neighbor 1, neighbor 2, neighbor 3, center].
        center, neighbors = (atoms[1], [atoms[0], atoms[2]]) if kind == 'angles' else (atoms[3], atoms[:3])
        edges = [tuple(sorted((center, neighbor))) for neighbor in neighbors]
        if any(edge not in graph for edge in edges):
            raise ValueError(f'{kind}: covalent neighbor of {center} absent from the selected state: {neighbors}')
        identity = (center, *sorted(neighbors))
        chemistry = {'bond_orders': [graph[edge] for edge in edges]}
    bounds = rule[measure.field]
    if (not isinstance(bounds, list) or len(bounds) != 2 or any(type(x) not in (int, float) for x in bounds)
            or not np.isfinite(bounds).all() or not bounds[0] < bounds[1]
            or (kind == 'bonds' and bounds[0] <= 0) or (kind == 'angles' and not 0 <= bounds[0] < bounds[1] <= 180)):
        raise ValueError(f'Supply a finite, ordered {measure.field} interval in {measure.unit}')
    basis = rule['basis']
    if not isinstance(basis, str) or not basis.strip():
        raise ValueError(f'Identify the source or modeling basis for the {measure.field} interval')
    # Motion uses initial physical IDs and nm throughout, even when product atom names change.
    points = []
    physical_atoms = []
    for atom in atoms:
        component, index, original = locations[atom]
        physical_atoms.append(original)
        points.append(np.asarray(motion[component][sample, index], dtype=float))
    value = measure_nm(kind, points)
    if not bounds[0] <= value <= bounds[1]:
        raise ValueError(f'Sample {sample}, {kind} {atoms}: {value:.8f} {measure.unit} outside {bounds} {measure.unit}')
    row = {'atoms': atoms, 'physical_atoms': physical_atoms, **chemistry,
           measure.field: round(value, 8), measure.interval: bounds, 'basis': basis}
    return identity, row


def check_pose_geometry(review, ledger, atom_ids, motion, frames):
    """Run after the ledger, complete initial atom set and finite nm arrays have been verified."""
    version = review['version']
    if type(version) is not int or version not in (1, 2) or ledger['version'] != 2:
        raise ValueError('Reaction geometry requires a version-1 or version-2 review and a version-2 reaction ledger')
    stages = ledger['stages'] if 'stages' in ledger else [ledger]
    by_state = {}
    timeline = {}
    if not review['poses']:
        raise ValueError('Select at least one pose for the geometry check')
    for pose in review['poses']:
        state, sample = pose['state'], pose['sample']
        if type(state) is not int or not 0 <= state <= len(stages):
            raise ValueError('Select state 0 (initial) or a stage product numbered from 1')
        if type(sample) is not int or not 0 <= sample <= frames or sample in timeline:
            raise ValueError('Each reviewed pose needs a unique in-range source sample')
        if version == 1 and ('angles' in pose or 'signed_heights' in pose):
            raise ValueError('Valence angles and signed heights require a version-2 geometry review')
        if version == 1 and not pose['bonds']:
            raise ValueError('Select at least one bond in each reviewed pose')
        if version == 2:
            if pose.keys() - {'state', 'sample'} - MEASURES.keys():
                raise ValueError('Unknown pose measurement in version-2 review')
            if any(not isinstance(pose[kind], list) for kind in MEASURES if kind in pose):
                raise ValueError('Supply pose measurements as lists of rules')
            if not any(pose.get(kind) for kind in MEASURES):
                raise ValueError('Select at least one geometry measurement in each pose')
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
            measured = dict(state=state, sample=pose['sample'])
            for kind in MEASURES:
                if kind not in pose:
                    continue
                seen = set()
                measured[kind] = []
                for rule in pose[kind]:
                    identity, row = _measure(kind, rule, graph, locations, motion, pose['sample'])
                    if identity in seen:
                        raise ValueError(f'Repeated {kind} measurement in one pose: {identity}')
                    seen.add(identity)
                    measured[kind].append(row)
            result.append(measured)
    return sorted(result, key=lambda pose: pose['sample'])


# Original Python entry point for bond-only reviews.
check_pose_distances = check_pose_geometry
