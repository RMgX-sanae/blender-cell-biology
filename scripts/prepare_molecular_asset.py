"""Prepare source-derived Gaussian molecular surfaces for Blender (outside Blender).

Requires numpy, scipy, scikit-image and gemmi. This is an illustrative molecular
envelope, NOT an experimental density map or an exact solvent-excluded surface.
The three presets encode selections for the named PDB entries only.
"""
import argparse
import hashlib
import json
from pathlib import Path

import gemmi
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
from skimage.measure import marching_cubes


def atom_table(path):
    structure = gemmi.read_structure(str(path))
    result = {}
    for chain in structure[0]:
        rows = []
        for residue in chain:
            # Entity-aware selection keeps modified polymer residues in mmCIF,
            # while excluding free ligands, water and crystallization aids.
            if residue.entity_type != gemmi.EntityType.Polymer:
                continue
            for atom in residue:
                if atom.element.is_hydrogen or atom.occ <= 0:
                    continue
                if atom.altloc not in ('\x00', 'A'):
                    continue
                rows.append((list(atom.pos), atom.element.name, atom.name,
                             residue.seqid.num, residue.name))
        if rows:
            result[chain.name] = rows
    return result


def xyz(rows):
    return np.array([r[0] for r in rows], dtype=float)


def orientation(table, preset):
    if preset == 'beta2':
        # UniProt P07550 canonical TM spans; PDB 3SN6 author R numbering.
        spans = [(30, 56), (71, 96), (104, 129), (150, 171),
                 (198, 220), (275, 298), (305, 326)]
        ends = []
        for i, (first, last) in enumerate(spans):
            rows = sorted((r for r in table['R'] if first <= r[3] <= last
                           and r[2] == 'CA'), key=lambda r: r[3])
            a, b = xyz(rows[:4]).mean(0), xyz(rows[-4:]).mean(0)
            ends.append((a, b) if i % 2 == 0 else (b, a))
        exterior = np.mean([v[0] for v in ends], axis=0)
        interior = np.mean([v[1] for v in ends], axis=0)
        origin = (exterior + interior) / 2
        z = exterior - interior
        z /= np.linalg.norm(z)
        x = np.array([1., 0., 0.])
        x -= z * np.dot(x, z)
        x /= np.linalg.norm(x)
        y = np.cross(z, x)
        return origin, np.array([x, y, z]).T
    if preset == 'motor':
        a, b, k = (xyz(table[c]).mean(0) for c in ('A', 'B', 'K'))
        origin = (a + b) / 2
        x = b - a  # Alpha to beta; illustrative protofilament + direction.
        x /= np.linalg.norm(x)
        z = k - origin
        z -= x * np.dot(z, x)
        z /= np.linalg.norm(z)
        y = np.cross(z, x)
        return origin, np.array([x, y, z]).T
    points = np.concatenate([xyz(r) for r in table.values()])
    return points.mean(0), np.eye(3)


def surface(points, spacing, sigma, level):
    padding = 4 * sigma
    origin = points.min(0) - padding
    shape = np.ceil((points.max(0) - origin + padding) / spacing).astype(int) + 1
    if np.prod(shape) > 100_000_000:
        raise ValueError('Grid exceeds 100 million voxels; increase --spacing or split chains.')
    grid = np.zeros(tuple(shape), dtype=np.float32)
    indices = np.rint((points - origin) / spacing).astype(int)
    np.add.at(grid, tuple(indices.T), 1.)
    pixels = sigma / spacing
    grid = gaussian_filter(grid, pixels) * (np.sqrt(2 * np.pi) * pixels) ** 3
    vertices, faces, _, _ = marching_cubes(grid, level=level,
                                          spacing=(spacing,) * 3)
    return (vertices + origin).astype(np.float32), faces.astype(np.int32)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('structure', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--preset', choices=['generic', 'beta2', 'transcription', 'motor'],
                        default='generic')
    parser.add_argument('--spacing', type=float, default=.75, help='Angstrom')
    parser.add_argument('--sigma', type=float, default=1.15, help='Gaussian sigma, Angstrom')
    parser.add_argument('--level', type=float, default=.55)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit('Output already exists; choose a fresh output directory.')
    expected = {'beta2': '3SN6', 'transcription': '5FLM', 'motor': '4HNA'}
    if args.preset in expected and args.structure.stem.upper() != expected[args.preset]:
        raise SystemExit(f'{args.preset} selection is validated only for {expected[args.preset]}')
    if min(args.spacing, args.sigma, args.level) <= 0:
        raise SystemExit('Surface parameters must be positive.')
    table = atom_table(args.structure)
    if args.preset == 'beta2':
        table = {c: table[c] for c in ('R', 'A', 'B', 'G')}
        table['R'] = [r for r in table['R'] if r[3] < 400]
    elif args.preset == 'motor':
        table = {c: table[c] for c in ('A', 'B', 'K')}
    origin, basis = orientation(table, args.preset)
    manifest = {
        'source_file': args.structure.name,
        'source_sha256': hashlib.sha256(args.structure.read_bytes()).hexdigest(),
        'entry_url': 'https://www.rcsb.org/structure/' + args.structure.stem.upper(),
        'preset': args.preset, 'units': 'nm', 'model_index': 0,
        'source_origin_angstrom': origin.tolist(), 'source_basis_columns': basis.tolist(),
        'surface': {'method': 'sum of equal atom-centred Gaussians; illustrative envelope',
                    'spacing_angstrom': args.spacing, 'sigma_angstrom': args.sigma,
                    'isovalue': args.level},
        'excluded': 'Hydrogens, non-polymer entities, nonpositive occupancy, alternative locations except A. Modified polymer residues retained.',
        'limitations': [
            'Coordinates are experimental structure models, with construct and missing-region limits.',
            'Surface is not experimental EM density, an exact SES, or a live-cell conformation.',
            'Any later replication, placement or animation is a schematic unless separately evidenced.'],
        'chains': []}
    if args.preset == 'beta2':
        manifest['selection_note'] = '3SN6 R residues <400 plus Gs A/B/G; exclude engineered T4 fusion and stabilizing nanobody N; omit ligand/ions.'
        manifest['orientation_note'] = 'Inferred from UniProt P07550 TM1–7 endpoint averages; +Z extracellular. Not a measured membrane plane or OPM result.'
        manifest['orientation_source'] = 'https://www.uniprot.org/uniprotkb/P07550/entry'
    if args.preset == 'motor':
        manifest['selection_note'] = '4HNA A alpha-tubulin, B beta-tubulin, K kinesin motor; exclude DARPin D. Does not supply an intact two-headed motor or microtubule lattice.'
    args.output.mkdir(parents=True)
    for chain, rows in table.items():
        coords = (xyz(rows) - origin) @ basis
        vertices, faces = surface(coords, args.spacing, args.sigma, args.level)
        names = np.array([r[2] for r in rows])
        residues = np.array([r[4] for r in rows])
        nucleic = bool(np.isin(residues, ['A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT']).mean() > .5)
        bonds = np.empty((0, 2), dtype=np.int32)
        if nucleic:
            pairs = cKDTree(coords).query_pairs(1.85, output_type='ndarray')
            distances = np.linalg.norm(coords[pairs[:, 0]] - coords[pairs[:, 1]], axis=1)
            bonds = pairs[distances > .4].astype(np.int32)
        np.savez_compressed(args.output / f'chain_{chain}.npz',
                            vertices=vertices * .1, faces=faces, atoms=coords.astype(np.float32) * .1,
                            elements=np.array([r[1] for r in rows]), atom_names=names,
                            residue_numbers=np.array([r[3] for r in rows]),
                            residue_names=residues, bonds=bonds)
        info = {'chain': chain, 'atoms': len(rows), 'vertices': len(vertices),
                'faces': len(faces), 'nucleic_acid': nucleic,
                'centroid_nm': (coords.mean(0) * .1).tolist(),
                'bounds_nm': [(coords.min(0) * .1).tolist(), (coords.max(0) * .1).tolist()]}
        manifest['chains'].append(info)
        print(json.dumps(info), flush=True)
    (args.output / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
