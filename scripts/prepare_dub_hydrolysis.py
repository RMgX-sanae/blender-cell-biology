"""Transfer the established atom/motion workflow to AMSH-LP K63 net hydrolysis.

This source-specific adapter starts bound and ends with separated products.
It does not implement association, a proton relay, solvent reload or a full cycle.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.spatial import cKDTree

from dub_geometry import DubGeometry, SOURCE_HASHES, oxygen_transport, angle
from interaction_motion import matrices
from prepare_molecular_asset import surface
from prepare_ubiquitin_transfer import smooth
from reaction_ledger import validate
from ubiquitin_geometry import digest


def sample(model, config, frames):
    curves = []
    for name in ('reaction', 'departure'):
        times = np.asarray(config[name+'_times']); values = np.asarray(config[name+'_parameters'])
        if times[0] != 0 or times[-1] != 1 or values.shape != (len(times), 24) or not np.isfinite(values).all():
            raise ValueError('Expected qualified 24-parameter control poses: '+name)
        curves.append(PchipInterpolator(times, values))
    reaction, departure = curves
    if not np.allclose(reaction(1), departure(0), atol=1e-9):
        raise ValueError('Reaction and release do not share the same product pose')
    endpoint = model.pose(reaction(1)); oxygen = oxygen_transport(model, endpoint)
    grid = np.linspace(0, 1, 2049)
    positions = []
    for q in grid:
        points = model.pose(departure(q))['distal']
        positions.append(np.vstack([points, oxygen(points)]))
    lengths = np.r_[0, np.cumsum(np.linalg.norm(np.diff(positions, axis=0), axis=2).max(axis=1))]
    lengths /= lengths[-1]
    time = np.arange(frames+1)/frames
    def progress(name):
        start, end = config['timing'][name]
        return smooth((time-start)/(end-start))
    fractions = progress('reaction'); release_q = np.interp(progress('distal_release'), lengths, grid)
    path = config['proximal_release']
    poses = matrices(endpoint['proximal']*.1, path['direction'], path['distance_nm'],
                     .5+.5*progress('proximal_release'), 'loop', path['turn_radians'], path['bend_fraction'])
    result = {k: [] for k in (*model.rows, 'zinc', 'water')}
    for i, q in enumerate(fractions):
        points = model.pose(reaction(q)); water = model.water
        if time[i] >= config['timing']['distal_release'][0]:
            points['distal'] = model.pose(departure(release_q[i]))['distal']; water = oxygen(points['distal'])
        points['proximal'] = points['proximal']@poses[i, :3, :3].T+poses[i, :3, 3]*10
        points.update(zinc=model.metals, water=water[None])
        for key, coordinates in points.items(): result[key].append(coordinates)
    return {k: np.asarray(v) for k, v in result.items()}, fractions


def audit(model, arrays, fractions):
    failures = []; minima = {}; endpoint = {}; bond_error = 0.; angle_error = 0.
    for key, edges in model.graphs.items():
        points = arrays[key]; source = model.source[key]
        reference = np.linalg.norm(source[edges[:, 0]]-source[edges[:, 1]], axis=1)
        observed = np.linalg.norm(points[:, edges[:, 0]]-points[:, edges[:, 1]], axis=2)
        bond_error = max(bond_error, float(np.max(abs(observed-reference))))
        neighbors = [set() for _ in source]
        for i, j in edges: neighbors[i].add(j); neighbors[j].add(i)
        pairs = np.array(sorted({tuple(sorted((i,j))) for group in neighbors for i in group for j in group if i != j}))
        reference = np.linalg.norm(source[pairs[:,0]]-source[pairs[:,1]], axis=1)
        observed = np.linalg.norm(points[:,pairs[:,0]]-points[:,pairs[:,1]], axis=2)
        angle_error = max(angle_error, float(np.max(abs(observed-reference))))
    for i, q in enumerate(fractions):
        p = {k: v[i] for k, v in arrays.items()}; w = p['water'][0]; metrics = {}
        for a, b in [('enzyme','distal'), ('enzyme','proximal'), ('distal','proximal')]:
            distances, ids = cKDTree(p[a]).query(p[b], k=5)
            if (a,b) == ('distal','proximal') and q < 1:
                for first, second in model.excluded: distances[second, ids[second] == first] = np.inf
            metrics[a+'/'+b] = float(distances.min())
        for key, edges in model.self_pairs.items():
            metrics[key+'_mobile_self'] = float(np.linalg.norm(p[key][edges[:,0]]-p[key][edges[:,1]], axis=1).min())
        for key in model.rows:
            distances = np.linalg.norm(p[key]-w, axis=1)
            if key == 'distal' and q > 0: distances[[model.ids[n] for n in ('C','CA','O')]] = np.inf
            metrics[key+'/water'] = float(distances.min())
        for key in ('distal','proximal'): metrics[key+'/zinc'] = float(cKDTree(p['zinc']).query(p[key])[0].min())
        for metal, excluded in enumerate(model.metal_ligands):
            distances = np.linalg.norm(p['enzyme']-p['zinc'][metal], axis=1)
            distances[excluded] = np.inf; metrics['enzyme/zinc_'+str(metal)] = float(distances.min())
        metrics['water/structural_zinc'] = float(np.linalg.norm(w-p['zinc'][0]))
        for key, value in metrics.items():
            minima[key] = min(minima.get(key, float('inf')), value)
            if value < 2.38: failures.append([i+1, key, value])
        if np.linalg.norm(w-p['zinc'][1]) < np.linalg.norm(model.water-model.zinc)-.01:
            failures.append([i+1, 'water/catalytic_zinc_overlap'])
        c,o,ca,n,ce,normal = model.values(p)
        cn = float(np.linalg.norm(c-n)); co = float(np.linalg.norm(c-w))
        if q == 0:
            geometry = {'C_N_A':cn, 'C_water_A':co, 'N_plane_A':float(abs((n-c)@normal)),
                        'CE_plane_A':float(abs((ce-c)@normal)), 'O_C_N_deg':float(np.degrees(angle(o,c,n))),
                        'C_N_CE_deg':float(np.degrees(angle(c,n,ce)))}
            if abs(cn-1.33) > .03 or co < 2.38 or geometry['N_plane_A'] > .1 or geometry['CE_plane_A'] > .1 or abs(geometry['O_C_N_deg']-123) > 10 or abs(geometry['C_N_CE_deg']-122) > 10:
                failures.append([i+1, 'reactant_geometry', geometry])
            endpoint['reactant'] = geometry
        elif q == 1:
            geometry = {'C_N_A':cn, 'C_water_A':co, 'water_plane_A':float(abs((w-c)@normal)),
                        'O_C_water_deg':float(np.degrees(angle(o,c,w)))}
            if abs(co-1.25) > .03 or cn < 2.38 or geometry['water_plane_A'] > .1 or abs(geometry['O_C_water_deg']-123) > 10:
                failures.append([i+1, 'product_geometry', geometry])
            endpoint['product'] = geometry
    metal_error = float(np.max(abs(arrays['zinc']-model.metals)))
    coordination_error = 0.
    for metal, ligands in enumerate(model.metal_ligands):
        expected = np.linalg.norm(model.source['enzyme'][ligands]-model.metals[metal], axis=1)
        observed = np.linalg.norm(arrays['enzyme'][:,ligands]-arrays['zinc'][:,metal,None], axis=2)
        coordination_error = max(coordination_error, float(np.max(abs(observed-expected))))
    finite = all(np.isfinite(v).all() for v in arrays.values())
    passed = finite and not failures and max(bond_error, angle_error, metal_error, coordination_error) < 1e-6
    return {'passed':passed, 'frames_checked':len(fractions), 'minimum_nonbonded_A':minima,
            'max_bond_length_change_A':bond_error, 'max_adjacent_bond_pair_distance_change_A':angle_error,
            'static_zinc_error_A':metal_error, 'metal_coordination_length_error_A':coordination_error,
            'endpoints':endpoint, 'failures':failures,
            'contact_policy':'2.38 A heavy-atom screen. Only named metal ligands and declared forming/breaking bond 1-2/1-3 neighborhoods are exempt. Restored Glu292 and mobile ubiquitin self contacts checked; other unchanged protein self contacts inherit the source.',
            'limits':'Geometric display constraints, not force-field energy, kinetics, proton relay or observed motion.'}


def ledger(model):
    atoms = [[f'{k}:{r[3]}:{r[2]}', r[1]] for k, rows in model.rows.items() for r in rows]
    atoms += [['zinc:structural','Zn'], ['zinc:catalytic','Zn'], ['water:517:O','O']]
    product = [[identity if identity != 'water:517:O' else 'distal:76:OXT', element] for identity, element in atoms]
    return {'version':1, 'scope':'Selected resolved experimental protein constructs, two zinc ions and one water oxygen; no hydrogens or bulk solvent.',
            'bond_scope':'Scissile isopeptide and distal terminal C-O adjacency only; no bond orders or proton transfer.',
            'reactant_atoms':atoms, 'product_atoms':product, 'atom_map':[[a[0],b[0]] for a,b in zip(atoms,product)],
            'reactant_bonds':[['distal:76:C','proximal:63:NZ'], ['distal:76:C','distal:76:O']],
            'product_bonds':[['distal:76:C','distal:76:OXT'], ['distal:76:C','distal:76:O']],
            'break':[['distal:76:C','proximal:63:NZ']], 'form':[['distal:76:C','water:517:O']]}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--config', type=Path, default=Path(__file__).resolve().parents[1]/'assets/dub-hydrolysis.json')
    p.add_argument('--frames', type=int, default=720); p.add_argument('--fps', type=int, default=24)
    a = p.parse_args()
    if a.output.exists(): p.error('Use a fresh preparation directory')
    if a.frames < 96 or a.fps < 1: p.error('Use at least 96 frames and positive fps')
    config = json.loads(a.config.read_text()); model = DubGeometry(a.sources)
    book = ledger(model); receipt = validate(book)
    arrays, fractions = sample(model, config, a.frames); checks = audit(model, arrays, fractions)
    if not checks['passed']: raise ValueError(checks)
    origin = model.source['enzyme'].mean(0)
    a.output.mkdir(parents=True)
    np.savez_compressed(a.output/'motion.npz', **{k:(v-origin)*.1 for k,v in arrays.items()}, chemical_fraction=fractions)
    for key, rows in model.rows.items():
        np.savez_compressed(a.output/(key+'.npz'), atoms=(model.source[key]-origin)*.1, bonds=model.graphs[key],
                            atom_names=[r[2] for r in rows], residue_numbers=[r[3] for r in rows])
        vertices, faces = surface(model.source[key]-model.source[key].mean(0), .8, 1.15, .55)
        np.savez_compressed(a.output/('background_'+key+'.npz'), vertices=vertices*.1, faces=faces)
    (a.output/'reaction-ledger.json').write_text(json.dumps(book, indent=2)+'\n')
    identities = {k:[f'{k}:{r[3]}:{r[2]}' for r in rows] for k,rows in model.rows.items()}
    identities.update(zinc=['zinc:structural','zinc:catalytic'], water=['water:517:O'])
    manifest = {'kind':'amsh-k63-net-hydrolysis', 'passed':True, 'frames':a.frames, 'fps':a.fps, 'atom_ids':identities,
                'reaction_ledger':receipt, 'source_reconstruction':model.record, 'source_sha256':SOURCE_HASHES,
                'checks':checks, 'origin_A':origin.tolist(), 'reaction_atom_indices':model.ids, 'config':config,
                'scope':config['scope'], 'retiming':'Max atom-travel arc-length and quintic easing; presentation time.',
                'maximum_frame_step_nm':float(max(np.linalg.norm(np.diff(v,axis=0),axis=2).max() for v in arrays.values()))*.1,
                'data_sha256':{f.name:digest(f) for f in a.output.iterdir()}}
    (a.output/'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'passed':True, 'atoms':receipt['mapped_atoms'], 'checked_poses':len(fractions), 'output':str(a.output)}))


if __name__ == '__main__': main()
