"""Prepare reusable actin, branch, cap and tubulin assets from four public mmCIFs.

Uses the existing Gemmi / NumPy / SciPy / scikit-image environment outside
Blender. Input: 8F8R, 8F8Q, 7AQK, 3J6F plus fetch_structure provenance files.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from prepare_molecular_asset import atom_table, surface, xyz


def rigid_fit(first, second):
    a = {r[3]: r[0] for r in first if r[2] == 'CA'}
    b = {r[3]: r[0] for r in second if r[2] == 'CA'}
    shared = sorted(a.keys() & b.keys())
    x, y = np.array([a[i] for i in shared]), np.array([b[i] for i in shared])
    u, _, vt = np.linalg.svd((x-x.mean(0)).T @ (y-y.mean(0)))
    rotation = vt.T @ np.diag([1, 1, np.linalg.det(vt.T @ u.T)]) @ u.T
    shift = y.mean(0) - rotation @ x.mean(0)
    rmsd = np.sqrt(np.mean(np.sum((x @ rotation.T + shift - y)**2, axis=1)))
    matrix = np.eye(4)
    matrix[:3, :3], matrix[:3, 3] = rotation, shift
    return matrix, float(rmsd)


def basis_for(x, z):
    x = x / np.linalg.norm(x)
    z = z - x * np.dot(x, z)
    z /= np.linalg.norm(z)
    return np.array([x, np.cross(z, x), z]).T


def transform_frame(origin, basis):
    matrix = np.eye(4)
    matrix[:3, :3] = basis.T * .1
    matrix[:3, 3] = -basis.T @ origin * .1
    return matrix


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        p.error('Use a fresh output directory; keep cached assets with their hashes')
    tables, sources = {}, {}
    for code in ('8F8R', '8F8Q', '7AQK', '3J6F'):
        file = args.source_dir / (code + '.cif')
        record = json.loads(file.with_suffix('.cif.source.json').read_text())
        digest = hashlib.sha256(file.read_bytes()).hexdigest()
        if digest != record['sha256'] or record['id'] != code:
            raise ValueError(f'{code}: source identity/hash does not match provenance')
        sources[code] = {k: record[k] for k in ('id', 'url', 'entry_url', 'sha256', 'license', 'license_url')}
        sources[code]['title'] = record['metadata']['struct']['title']
        sources[code]['resolution_angstrom'] = record['metadata']['rcsb_entry_info'].get('resolution_combined')
        tables[code] = atom_table(file)

    a = tables['8F8R']
    step, rmsd = rigid_fit(a['D'], a['E'])
    rot, shift = step[:3, :3], step[:3, 3]
    _, _, vt = np.linalg.svd(rot-np.eye(3))
    axis = vt[-1]
    if np.dot(axis, shift) < 0:
        axis *= -1
    rise = float(np.dot(axis, shift))
    axis_point = np.linalg.pinv(np.eye(3)-rot) @ (shift-rise*axis)
    center = xyz(a['D']).mean(0)
    origin = axis_point + axis*np.dot(center-axis_point, axis)
    basis = basis_for(axis, center-origin)
    frame = transform_frame(origin, basis)
    step_nm = frame @ step @ np.linalg.inv(frame)
    if not 2.6 < rise*.1 < 3.0 or rmsd > .8:
        raise ValueError('Unexpected actin repeat; inspect the source before extending it')
    specs = {'actin': {'source': '8F8R', 'chains': ['D'], 'origin': origin, 'basis': basis,
                       'step_matrix_nm': step_nm.tolist(), 'fit_rmsd_angstrom': rmsd,
                       'rise_nm': rise*.1,
                       'polarity': '+X is toward the barbed end, D→E→F→G in 8F8R; cap-side cross-check with 8F8Q',
                       'limits': 'One internal F-actin ADP conformation is repeated. Incoming units are rigid F-actin shape proxies, not a measured G→F transition; nucleotide chemistry is omitted.'}}
    cap = tables['8F8Q']
    origin = np.mean([xyz(cap[c]).mean(0) for c in 'ABCDEF'], axis=0)
    axis = xyz(cap['F']).mean(0)-xyz(cap['A']).mean(0)
    basis = basis_for(axis, xyz(cap['E']).mean(0)-xyz(cap['F']).mean(0))
    specs['cap'] = {'source': '8F8Q', 'chains': list('ABCDEFGH'), 'origin': origin, 'basis': basis,
                    'actin_chains': list('ABCDEF'), 'cap_chains': ['G', 'H'],
                    'limits': 'Rabbit actin and human CapZ are the published mixed-species construct. Bound endpoint retained; approach is authored; no CapZ conformational cycle.'}
    branch = tables['7AQK']
    origin = np.mean([xyz(branch[c]).mean(0) for c in 'abcdefg'], axis=0)
    basis = basis_for(xyz(branch['r']).mean(0)-xyz(branch['k']).mean(0),
                      xyz(branch['j']).mean(0)-xyz(branch['h']).mean(0))
    step, error = rigid_fit(branch['i'], branch['j'])
    frame = transform_frame(origin, basis)
    specs['branch'] = {'source': '7AQK', 'chains': list('abcdefghijklmnopqr'), 'origin': origin, 'basis': basis,
                       'mother_chains': list('klmnopqr'), 'daughter_chains': list('hij'),
                       'arp_chains': list('abcdefg'), 'step_matrix_nm': (frame@step@np.linalg.inv(frame)).tolist(),
                       'fit_rmsd_angstrom': error,
                       'limits': '9 Å cellular branch model; backbone-derived display envelope. Existing seven-subunit junction is retained; daughter extension is authored, without NPF activation or de novo nucleation.'}
    tub = tables['3J6F']
    shift = xyz(tub['I']).mean(0)-xyz(tub['A']).mean(0)
    origin = np.mean([xyz(tub[c]).mean(0) for c in 'ABCDEF'], axis=0)
    basis = basis_for(shift, xyz(tub['E']).mean(0)-xyz(tub['C']).mean(0))
    specs['tubulin'] = {'source': '3J6F', 'chains': list('ABCDEF'), 'origin': origin, 'basis': basis,
                        'dimers': [['A', 'B'], ['C', 'D'], ['E', 'F']],
                        'repeat_nm': (shift@basis*.1).tolist(),
                        'polarity': '+X follows alpha A→beta B→alpha I; beta is the terminal subunit',
                        'limits': 'Three complete protofilament strips from a 14-protofilament reference. Other protofilaments are omitted whole. GDP-lattice shape proxy; growth, peeling and rescue are scripted, with no GTP cap or hydrolysis simulation.'}
    args.output.mkdir(parents=True)
    for name, spec in specs.items():
        folder = args.output / name
        folder.mkdir()
        table = tables[spec['source']]
        origin, basis = spec.pop('origin'), spec.pop('basis')
        spec.update(units='nm', source_origin_angstrom=origin.tolist(), source_basis_columns=basis.tolist(),
                    provenance=sources[spec['source']], geometry=[])
        for chain in spec['chains']:
            rows = table[chain]
            coordinates = (xyz(rows)-origin)@basis
            sigma = 1.65 if name == 'branch' else 1.15
            vertices, faces = surface(coordinates, .85, sigma, .55)
            np.savez_compressed(folder/f'chain_{chain}.npz', vertices=vertices*.1, faces=faces,
                                atoms=coordinates.astype(np.float32)*.1,
                                elements=np.array([r[1] for r in rows]), atom_names=np.array([r[2] for r in rows]),
                                residue_numbers=np.array([r[3] for r in rows]))
            spec['geometry'].append({'chain': chain, 'atoms': len(rows), 'vertices': len(vertices), 'faces': len(faces)})
        spec['surface'] = 'Atom-centred Gaussian display envelope, not experimental density or exact SES'
        (folder/'manifest.json').write_text(json.dumps(spec, ensure_ascii=False, indent=2)+'\n')
    (args.output/'manifest.json').write_text(json.dumps({'sources': sources, 'presets': list(specs),
            'units': 'nm', 'asset_pipeline': 'Reuse public coordinates → prepare once → link meshes per instance'}, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'output': str(args.output), 'presets': list(specs), 'actin_rise_nm': rise*.1,
                       'actin_repeat_fit_rmsd_angstrom': rmsd}, ensure_ascii=False))


if __name__ == '__main__':
    main()
