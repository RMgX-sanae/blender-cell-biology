"""Rebuild the explicitly illustrative FH2 seed assembly from 9B27.

A functional lobe contains one chain's body and the partner's lasso. Separate
only pre-existing disconnected mesh components at unresolved linker gaps: no
face, vertex or modeled atom is deleted, cut or elastically opened.
"""
import hashlib
import itertools
import json

import numpy as np
from scipy.interpolate import BSpline
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from prepare_contact_repairs import repair_path, repair_rigid_release
from project_context import ROOT


def smooth(value):
    value = np.clip(value, 0, 1)
    return value**3 * (10 + value * (-15 + 6 * value))


def motion_matrix(fraction, record):
    offset = np.array(record['offset']);side = np.cross(offset, [0, 0, 1])
    side *= 1.2 / np.linalg.norm(side);far = max(0, (fraction - .5) * 2)
    shift = offset * fraction + side * np.sin(np.pi * far)**2
    correction = np.zeros(6);pivot = np.zeros(3);detour = record['detour']
    if detour and detour['first_fraction'] < fraction < detour['last_fraction']:
        u = (fraction - detour['first_fraction']) / (detour['last_fraction'] - detour['first_fraction'])
        controls = detour.get('controls_nm_radians', detour.get('controls_nm'))
        value = BSpline(detour['knots'], controls, 3)(u)
        correction[:len(value)] = value;pivot = np.array(detour.get('pivot_nm', pivot))
    matrix = np.eye(4);matrix[:3, :3] = Rotation.from_rotvec(correction[3:]).as_matrix()
    matrix[:3, 3] = pivot - matrix[:3, :3] @ pivot + shift + correction[:3]
    return matrix


def main():
    assets = ROOT / 'assets';folder = assets / 'fh2_seed';folder.mkdir(exist_ok=True)
    data = {c:np.load(assets / '9B27' / (c + '.npz')) for c in 'DEFGH'}
    atoms = {c:data[c]['atoms'].astype(float) for c in 'DEF'};lobes = {};ownership = []
    output = {'schema_version':1, 'source':'9B27',
        'source_sha256':hashlib.sha256((ROOT / 'sources/9B27.cif').read_bytes()).hexdigest(),
        'source_url':'https://www.rcsb.org/structure/9B27',
        'primary_reference':'https://doi.org/10.1038/s41586-024-07637-0',
        'limits':'9B27 is a mature mDia1 barbed-end reference, not an observed nucleation trajectory. Sequential assembly, lobe opening and missing-linker curves are authored teaching motions; no MD or nucleation rate is inferred.',
        'lobes':{}, 'motions':{}, 'intervals':{'D':[.07,.25], 'E':[.29,.48], 'F':[.58,.78]},
        'opening':[.48,.57], 'closing':[.80,.91], 'linkers':[]}
    for body, partner in [('G','H'), ('H','G')]:
        vertices = [];faces = [];points = [];materials = [];vertex_ids = [];atom_ids = []
        for material_index, (chain, lasso) in enumerate([(body,False), (partner,True)]):
            source = data[chain];mask = (source['residue_numbers'] <= 804) == lasso
            nearest = cKDTree(source['atoms']).query(source['vertices'])[1];selected = mask[nearest]
            membership = selected[source['faces']]
            if np.any(np.ptp(membership.astype(int), axis=1)):
                raise ValueError('FH2 partition would cut a real mesh face; review the source')
            keep = membership.all(1);remap = np.full(len(selected), -1)
            remap[selected] = np.arange(selected.sum()) + sum(len(v) for v in vertices)
            vertices.append(source['vertices'][selected]);faces.append(remap[source['faces'][keep]])
            points.append(source['atoms'][mask]);materials.extend([material_index] * int(keep.sum()))
            vertex_ids.extend((chain,int(i)) for i in np.flatnonzero(selected))
            atom_ids.extend((chain,int(i)) for i in np.flatnonzero(mask))
        lobes[body] = np.concatenate(points).astype(float);ownership.extend(atom_ids)
        np.savez_compressed(folder / (body + '.npz'), vertices=np.concatenate(vertices),
            faces=np.concatenate(faces), atoms=lobes[body], material_indices=np.array(materials))
        output['lobes'][body] = {'body_chain':body, 'partner_lasso_chain':partner,
            'source_atom_ids':atom_ids, 'source_vertex_ids':vertex_ids,
            'source_faces':sum(len(f) for f in faces)}
    expected = {(c,i) for c in 'GH' for i in range(len(data[c]['atoms']))}
    if len(ownership) != len(expected) or set(ownership) != expected:
        raise ValueError('FH2 atom ownership is incomplete or duplicated')
    all_vertices = [tuple(v) for lobe in output['lobes'].values() for v in lobe['source_vertex_ids']]
    if len(all_vertices) != len(set(all_vertices)) or len(all_vertices) != sum(len(data[c]['vertices']) for c in 'GH'):
        raise ValueError('FH2 source vertices were not conserved')
    if sum(l['source_faces'] for l in output['lobes'].values()) != sum(len(data[c]['faces']) for c in 'GH'):
        raise ValueError('FH2 source faces were not conserved')
    for chain, first, last, start_lobe, end_lobe in [('G',803,824,'H','G'), ('H',804,836,'G','H')]:
        ca = dict(zip(data[chain]['ca_residues'], data[chain]['ca']))
        if any(first<n<last for n in ca):raise ValueError('Expected unresolved FH2 linker is modeled in this source')
        output['linkers'].append({'chain':chain, 'residue_interval':[first,last],
            'start_lobe':start_lobe, 'end_lobe':end_lobe, 'start':ca[first].tolist(), 'end':ca[last].tolist(),
            'contour_upper_bound_nm':.38 * (last - first), 'representation':'unresolved linker, schematic curve'})
    # Authored approach directions in the prepared 9B27 frame. Recalculate all
    # local corrections and clearance checks from current source coordinates.
    offsets = {'D':[-9.085450920998348,4.141063122478634,.5533333333333335],
               'E':[-8.841442599388511,1.8272637362236843,-4.3],
               'F':[7.882706520505707,6.075817651460504,.9733333333333338]}
    head_offset = np.array([2.6165610346469768,-1.050318660915411,2.07375])
    fixed = np.concatenate([lobes['G'],lobes['H']])
    for chain in 'DE':
        output['motions'][chain] = repair_rigid_release(atoms[chain], fixed, atoms[chain].mean(0),
            offsets[chain], samples=151, clearance=.23, iterations=700)
        fixed = np.concatenate([fixed, atoms[chain]])
        print('FH2_SEED_PATH', chain, output['motions'][chain]['after_minimum_nm'], flush=True)
    output['motions']['H'] = repair_path(lobes['H'], np.concatenate([atoms[c] for c in 'DEF'] + [lobes['G']]), head_offset)
    output['motions']['F'] = repair_rigid_release(atoms['F'],
        np.concatenate([atoms['D'], atoms['E'], lobes['G'], lobes['H'] + head_offset]),
        atoms['F'].mean(0), offsets['F'], samples=101, clearance=.23, iterations=600)
    # Include parked future monomers and all lobe motion, not just each isolated
    # approach. Dense parameter sampling is independent of the chosen film fps.
    base = {**atoms, **lobes};minima = {a+b:float('inf') for a,b in itertools.combinations(base,2)}
    for u in np.linspace(0,1,1001):
        posed = {}
        for chain, points in base.items():
            if chain == 'G':matrix = np.eye(4)
            else:
                if chain == 'H':fraction = smooth((u-.48)/.09) * (1-smooth((u-.80)/.11))
                else:
                    first,last = output['intervals'][chain];fraction = 1-smooth((u-first)/(last-first))
                matrix = motion_matrix(fraction, output['motions'][chain])
            posed[chain] = points @ matrix[:3,:3].T + matrix[:3,3]
        for a,b in itertools.combinations(base,2):
            minima[a+b] = min(minima[a+b],float(cKDTree(posed[b]).query(posed[a])[0].min()))
    output['sequence_screen'] = {'samples':1001, 'threshold_nm':.20, 'pair_minima_nm':minima}
    if min(minima.values()) < .20:raise ValueError('FH2 seed scene still has an overlapping pair: '+str(minima))
    (assets / 'fh2-seed-repairs.json').write_text(json.dumps(output,indent=2)+'\n')
    print('FH2_SEED_READY', min(minima.values()), flush=True)


if __name__ == '__main__':main()
