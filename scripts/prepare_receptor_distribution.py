"""Place intact oriented receptor assets on a world-space membrane (outside Blender).

Fixed-count, area-weighted rejection sampling with a relative intensity field.
Hard-core exclusion changes the realized field; this is not equilibrium sampling,
measured abundance, a diffusion simulation, or a molecular force-field calculation.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def surface_table(vertices, faces):
    vertices = np.asarray(vertices, dtype=float)
    raw_faces = np.asarray(faces)
    if vertices.ndim != 2 or vertices.shape[1:] != (3,) or not np.isfinite(vertices).all():
        raise ValueError('Expected finite world-space vertices in nm')
    if (raw_faces.ndim != 2 or raw_faces.shape[1:] != (3,) or
            not np.issubdtype(raw_faces.dtype, np.integer) or
            np.any(raw_faces < 0) or np.any(raw_faces >= len(vertices))):
        raise ValueError('Expected valid triangular vertex indices')
    triangles = vertices[raw_faces]
    cross = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
    length = np.linalg.norm(cross, axis=1)
    keep = length > 1e-12
    if not np.any(keep):
        raise ValueError('No nondegenerate membrane surface')
    return triangles[keep], cross[keep]/length[keep, None], length[keep]*.5, np.flatnonzero(keep)


def draw_surface(table, count, rng):
    triangles, normals, area, original = table
    ids = rng.choice(len(area), count, p=area/area.sum())
    u, v = np.sqrt(rng.random(count)), rng.random(count)
    bary = np.column_stack((1-u, u*(1-v), u*v))
    points = np.einsum('ni,nij->nj', bary, triangles[ids])
    return points, normals[ids], original[ids]


def ellipsoid(axes, latitudes=64, longitudes=128):
    axes = np.asarray(axes, float)
    if axes.shape != (3,) or not np.isfinite(axes).all() or np.min(axes) <= 0:
        raise ValueError('Ellipsoid semiaxes must be three positive lengths in nm')
    vertices = [[0., 0., axes[2]]]
    for i in range(1, latitudes):
        theta = np.pi*i/latitudes
        for j in range(longitudes):
            phi = 2*np.pi*j/longitudes
            vertices.append(axes*np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]))
    vertices.append([0., 0., -axes[2]])
    def ring(i, j): return 1+i*longitudes+j % longitudes
    faces = []
    for j in range(longitudes):
        faces.append([0, ring(0, j), ring(0, j+1)])
        faces.append([len(vertices)-1, ring(latitudes-2, j+1), ring(latitudes-2, j)])
    for i in range(latitudes-2):
        for j in range(longitudes):
            a, b, c, d = ring(i, j), ring(i, j+1), ring(i+1, j), ring(i+1, j+1)
            faces.extend([[a, c, b], [b, c, d]])
    vertices, faces = np.asarray(vertices), np.asarray(faces)
    tri = vertices[faces]
    inward = np.einsum('ij,ij->i', np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0]), tri.mean(1)) < 0
    faces[inward] = faces[inward, ::-1]
    return vertices, faces


def intensity(config, point):
    for region in config.get('excluded_spheres', []):
        if np.linalg.norm(point-np.asarray(region['center_nm'])) < region['radius_nm']:
            return 0.
    field = config['field']
    if field['kind'] == 'uniform': return 1.
    if field['kind'] == 'clustered':
        return field['background']+sum(c['amplitude']*np.exp(-np.sum((point-c['center_nm'])**2)/(2*c['sigma_nm']**2)) for c in field['clusters'])
    if field['kind'] == 'polarized':
        normal = np.asarray(field['plane_normal'], float)
        return field['enrichment'] if np.dot(point, normal/np.linalg.norm(normal)) >= field['plane_offset_nm'] else 1.
    raise ValueError('Unknown intensity field')


def field_bound(config):
    f = config['field']
    if f['kind'] == 'uniform': return 1.
    if f['kind'] == 'polarized':
        if f['enrichment'] <= 0 or np.linalg.norm(f['plane_normal']) == 0:
            raise ValueError('Invalid polar enrichment/normal')
        return max(1., f['enrichment'])
    if f['kind'] == 'clustered':
        if f['background'] < 0 or any(c['amplitude'] < 0 or c['sigma_nm'] <= 0 for c in f['clusters']):
            raise ValueError('Invalid cluster field')
        return f['background']+sum(c['amplitude'] for c in f['clusters'])
    raise ValueError('Unknown intensity field')


def sample(table, count, atoms, weight, max_weight, *, min_anchor_distance,
           atom_clearance=.20, normal_sign=1, seed=23, max_attempts=100000):
    atoms = np.asarray(atoms, float)
    if atoms.ndim != 2 or atoms.shape[1:] != (3,) or not len(atoms) or not np.isfinite(atoms).all():
        raise ValueError('Expected nonempty finite prototype heavy atoms in nm')
    if (count < 0 or int(count) != count or min_anchor_distance < 0 or atom_clearance <= 0 or
            not np.isfinite([min_anchor_distance, atom_clearance, max_weight]).all() or
            max_weight <= 0 or normal_sign not in (-1, 1) or max_attempts < 1):
        raise ValueError('Invalid count, distance, weight bound, normal sign or attempt limit')
    rng = np.random.default_rng(seed)
    positions, rotations, face_ids, trees = [], [], [], []
    radius = np.linalg.norm(atoms, axis=1).max()
    attempts = 0
    while len(positions) < count and attempts < max_attempts:
        attempts += 1
        points, normals, ids = draw_surface(table, 1, rng)
        point, normal = points[0], normals[0]*normal_sign
        value = float(weight(point))
        if not np.isfinite(value) or value < 0 or value > max_weight+1e-12:
            raise ValueError('Intensity field outside its declared [0,max_weight] bound')
        if rng.uniform(0, max_weight) >= value: continue
        if positions and np.linalg.norm(np.array(positions)-point, axis=1).min() < min_anchor_distance: continue
        axis = np.array([1., 0, 0]) if abs(normal[0]) < .9 else np.array([0., 1, 0])
        x = axis-normal*np.dot(axis, normal); x /= np.linalg.norm(x); y = np.cross(normal, x)
        angle = rng.uniform(0, 2*np.pi)
        basis = np.column_stack((x*np.cos(angle)+y*np.sin(angle), -x*np.sin(angle)+y*np.cos(angle), normal))
        transformed = atoms@basis.T+point
        nearby = [i for i, p in enumerate(positions) if np.linalg.norm(p-point) < 2*radius+atom_clearance]
        if any(np.min(trees[i].query(transformed, k=1, distance_upper_bound=atom_clearance)[0]) < atom_clearance for i in nearby): continue
        positions.append(point); rotations.append(basis); face_ids.append(int(ids[0])); trees.append(cKDTree(transformed))
    if len(positions) != count:
        raise ValueError(f'Packing not reached: {len(positions)}/{count}; do not silently lower the count')
    return {'positions': np.asarray(positions).reshape(-1, 3), 'rotations': np.asarray(rotations).reshape(-1, 3, 3),
            'face_ids': np.asarray(face_ids, dtype=int), 'attempts': attempts, 'surface_area_nm2': float(table[2].sum())}


def lipid_points(vertices, faces, protein_atoms, seed):
    # Complete schematic lipids: one head and two tails. Surface-area proposals
    # and hard-core thinning avoid a UV lattice. These are not MD equilibrated.
    tri = vertices[faces]; cross = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
    normals = np.zeros_like(vertices)
    for corner in range(3): np.add.at(normals, faces[:, corner], cross)
    lengths = np.linalg.norm(normals, axis=1)
    if np.any(lengths == 0): raise ValueError('Membrane has unused vertices or undefined vertex normals')
    normals /= lengths[:, None]
    rng = np.random.default_rng(seed); protein = cKDTree(protein_atoms)
    heads, directions = [], []
    for side in (-1, 1):
        table = surface_table(vertices+side*1.8*normals, faces)
        points, outward, _ = draw_surface(table, int(table[2].sum()*6), rng)
        into = -side*outward
        keep = np.ones(len(points), dtype=bool)
        # .48 nm bounds head/tail radii plus molecular display-envelope margin.
        for distance in np.linspace(0, 1.6, 9):
            keep &= protein.query(points+distance*into, distance_upper_bound=.48)[0] >= .48
        points, into = points[keep], into[keep]
        neighbors = cKDTree(points); available = np.ones(len(points), dtype=bool); accepted = []
        for i, point in enumerate(points):
            if not available[i]: continue
            accepted.append(i); available[neighbors.query_ball_point(point, .72)] = False
        heads.append(points[accepted]); directions.append(into[accepted])
    return np.concatenate(heads).astype(np.float32), np.concatenate(directions).astype(np.float32)


def membrane_shell(vertices, faces, protein_atoms):
    # Refine only the continuum display mesh so conservative protein holes do
    # not become large triangular craters on a coarsely tessellated membrane.
    for _ in range(2):
        edges = np.sort(faces[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), axis=1)
        unique, inverse = np.unique(edges, axis=0, return_inverse=True)
        ab, bc, ca = (inverse.reshape(-1, 3)+len(vertices)).T
        vertices = np.vstack([vertices, vertices[unique].mean(1)])
        a, b, c = faces.T
        faces = np.vstack([np.column_stack([a, ab, ca]), np.column_stack([ab, b, bc]),
                           np.column_stack([ca, bc, c]), np.column_stack([ab, bc, ca])])
    tri = vertices[faces]; cross = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
    normals = np.zeros_like(vertices)
    for corner in range(3): np.add.at(normals, faces[:, corner], cross)
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    protein = cKDTree(protein_atoms); shell_vertices, shell_faces = [], []
    for side in (-1, 1):
        points = vertices+side*1.52*normals; triangles = points[faces]
        centers = triangles.mean(1); radius = np.linalg.norm(triangles-centers[:, None], axis=2).max(1)
        # Conservative whole-triangle exclusion. The continuum core fills gaps
        # between illustrative lipids; it never trims receptor or lipid meshes.
        keep = protein.query(centers)[0]-radius >= .30
        shell_faces.extend(faces[keep]+len(shell_vertices)); shell_vertices.extend(points)
    return np.asarray(shell_vertices), np.asarray(shell_faces)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir', type=Path, required=True)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--distribution', choices=['uniform', 'clustered', 'polarized'])
    a = p.parse_args()
    if a.output.exists(): p.error('Use a fresh distribution output directory')
    config = json.loads(a.config.read_text())
    if config['count'] < 1 or len(set(config['chains'])) != len(config['chains']):
        p.error('Require a positive receptor count and unique selected chains')
    mode = a.distribution or config['distribution']
    config['distribution'] = mode; config['field'] = config['fields'][mode]
    assets = a.asset_dir.resolve(); source = json.loads((assets/'manifest.json').read_text())
    if source['units'] != 'nm' or not config.get('orientation_evidence'):
        raise ValueError('Require nm assets and explicit local +Z orientation evidence')
    if config.get('source_preset') != source['preset'] or config.get('source_id') != Path(source['source_file']).stem.upper():
        raise ValueError('Configuration/source identity mismatch; requalify another receptor')
    chains = config['chains']
    if set(item['chain'] for item in config['display']['legend']) != set(chains):
        p.error('Legend must cover exactly the selected chains')
    atom_arrays = [np.load(assets/f'chain_{chain}.npz')['atoms'] for chain in chains]
    atoms = np.concatenate(atom_arrays)
    if config['surface']['kind'] == 'ellipsoid':
        vertices, faces = ellipsoid(config['surface']['semiaxes_nm'])
    elif config['surface']['kind'] == 'mesh':
        surface_file = (a.config.parent/config['surface']['file']).resolve()
        data = np.load(surface_file); vertices, faces = data['vertices'], data['faces']
        config['surface']['sha256'] = digest(surface_file)
    else: raise ValueError('Use an ellipsoid or an explicitly oriented world-space mesh')
    table = surface_table(vertices, faces)
    sign = {'outward': 1, 'inward': -1}[config['local_positive_z']]
    placed = sample(table, config['count'], atoms, lambda point: intensity(config, point), field_bound(config),
        min_anchor_distance=config['minimum_anchor_nm'], atom_clearance=config['atom_clearance_nm'],
        normal_sign=sign, seed=config['seed'], max_attempts=config.get('max_attempts', 100000))
    world = [atoms@rotation.T+position for position, rotation in zip(placed['positions'], placed['rotations'])]
    minimum = None
    for i, points in enumerate(world):
        if i: minimum = min(minimum if minimum is not None else float('inf'), float(cKDTree(np.concatenate(world[:i])).query(points)[0].min()))
    world_atoms = np.concatenate(world)
    heads, directions = lipid_points(vertices, faces, world_atoms, config['seed']+19)
    shell_vertices, shell_faces = membrane_shell(vertices, faces, world_atoms)
    a.output.mkdir(parents=True)
    np.savez_compressed(a.output/'distribution.npz', vertices=vertices, faces=faces,
        positions=placed['positions'], rotations=placed['rotations'], face_ids=placed['face_ids'],
        lipid_heads=heads, lipid_directions=directions, shell_vertices=shell_vertices, shell_faces=shell_faces)
    record = {'kind': 'source-receptor-distribution', 'config': config, 'asset_dir': str(assets),
        'source': source, 'source_manifest_sha256': digest(assets/'manifest.json'),
        'chain_hashes': {c: digest(assets/f'chain_{c}.npz') for c in chains},
        'atoms_per_instance': len(atoms), 'atoms_by_chain': dict(zip(chains, map(len, atom_arrays))),
        'surface_area_nm2': placed['surface_area_nm2'], 'attempts': placed['attempts'],
        'minimum_intercomplex_atom_distance_nm': minimum, 'lipids': len(heads),
        'data_sha256': digest(a.output/'distribution.npz'),
        'limits': ['Fixed count and relative intensity field are authored hypotheses, not measured density or equilibrium.',
            'Anchor exclusions apply to anchors; intercomplex clearance screens all selected resolved heavy atoms.',
            'Missing atoms, glycans, hydrogens, ligands and membrane molecular energetics are outside this geometry screen.',
            'Two intact lipid leaflets are a procedural illustration, not a lipidomics or MD model.',
            'A continuum bilayer core fills display gaps between lipids; triangles near protein atoms are omitted conservatively.',
            'No receptor diffusion or activation kinetics are inferred from this static placement.']}
    (a.output/'manifest.json').write_text(json.dumps(record, ensure_ascii=False, indent=2))
    print(json.dumps({k: record[k] for k in ['atoms_per_instance', 'surface_area_nm2', 'minimum_intercomplex_atom_distance_nm', 'lipids']}))


if __name__ == '__main__': main()
