"""Preserve polymer AND nucleotide/ion coordinates for 8A2R/S/T state comparisons.

Run outside Blender with Gemmi, NumPy, SciPy and scikit-image. These three
structures are rabbit skeletal alpha F-actin, not free G-actin or human beta-actin.
"""
import argparse
import hashlib
import json
from pathlib import Path

import gemmi
import numpy as np

from prepare_molecular_asset import atom_table, surface, xyz
from prepare_filament_assets import rigid_fit


STATES = {
    '8A2R': ('Mg-ADP-BeF3', {'ADP': 27, 'MG': 1, 'BEF': 4}),
    '8A2S': ('Mg-ADP-Pi', {'ADP': 27, 'MG': 1, 'PO4': 5}),
    '8A2T': ('Mg-ADP', {'ADP': 27, 'MG': 1}),
}


def prepare(source, output):
    if output.exists():
        raise ValueError('Choose a fresh asset directory')
    tables = {code: atom_table(source / (code + '.cif')) for code in STATES}
    reference = xyz(tables['8A2S']['C'])
    origin = reference.mean(0)
    _, _, axes = np.linalg.svd(reference-origin, full_matrices=False)
    z = axes[0]
    if z[2] < 0:
        z = -z
    # A common viewing frame; no per-state shape scaling or morphing.
    y = np.array([0., 1., 0.])
    y -= z * np.dot(y, z)
    y /= np.linalg.norm(y)
    basis = np.array([np.cross(y, z), y, z]).T
    output.mkdir(parents=True)
    records = []
    for code, (state, expected) in STATES.items():
        file = source / (code + '.cif')
        provenance = json.loads(file.with_suffix('.cif.source.json').read_text())
        digest = hashlib.sha256(file.read_bytes()).hexdigest()
        if provenance['id'] != code or provenance['sha256'] != digest:
            raise ValueError(f'{code}: provenance mismatch')
        transform, rmsd = rigid_fit(tables[code]['C'], tables['8A2S']['C'])
        def aligned(points):
            return ((points @ transform[:3, :3].T + transform[:3, 3])-origin) @ basis * .1
        protein = aligned(xyz(tables[code]['C']))
        vertices, faces = surface(protein * 10, .75, 1.15, .55)
        ligands = []
        chain = gemmi.read_structure(str(file))[0]['C']
        for residue in chain:
            if residue.name not in expected or residue.entity_type == gemmi.EntityType.Polymer:
                continue
            atoms = [a for a in residue if not a.element.is_hydrogen and a.occ > 0
                     and a.altloc in ('\x00', 'A')]
            ligands.append({'ccd': residue.name, 'author_chain': 'C',
                            'author_residue': str(residue.seqid),
                            'atom_names': [a.name for a in atoms],
                            'elements': [a.element.name for a in atoms],
                            'positions_nm': aligned(np.array([list(a.pos) for a in atoms])).tolist()})
        actual = {r['ccd']: len(r['atom_names']) for r in ligands}
        if actual != expected or len(ligands) != len(expected):
            raise ValueError(f'{code}: expected ligands {expected}; found {actual}')
        np.savez_compressed(output / (code + '.npz'), vertices=vertices*.1, faces=faces,
                            atoms=protein.astype(np.float32))
        records.append({'id': code, 'state': state, 'author_chain': 'C',
                        'species': 'Oryctolagus cuniculus', 'isoform': 'skeletal alpha-actin',
                        'assembly_view': 'One complete modeled F-actin subunit; other four chains omitted whole',
                        'ligands': ligands, 'protein_heavy_atoms': len(protein),
                        'alignment_to_8A2S_C_rmsd_angstrom': rmsd,
                        'alignment_matrix_angstrom': transform.tolist(),
                        'source_sha256': digest, 'entry_url': provenance['entry_url'],
                        'source_url': provenance['url'], 'license': provenance['license'],
                        'resolution_angstrom': provenance['metadata']['rcsb_entry_info']['resolution_combined']})
    manifest = {'units': 'nm', 'states': records,
                'source_origin_angstrom': origin.tolist(), 'source_basis_columns': basis.tolist(),
                'paper': 'https://www.nature.com/articles/s41586-022-05241-8',
                'limitations': [
                    '8A2R contains ADP-BeF3, an ATP-state mimic; it does not contain ATP.',
                    'Static experimental coordinate models, not a reaction trajectory or free G-actin.',
                    'Hydrogens and waters omitted; no catalytic water or proton-transfer claim.',
                    'Protein surface is a Gaussian display envelope, not measured EM density.',
                    'No single molecule is cut. Ligands keep the same rigid transform as their protein.']}
    (output / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'output': str(output), 'states': [r['id'] for r in records],
                      'ligand_atom_counts': [{l['ccd']: len(l['atom_names']) for l in r['ligands']} for r in records]}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    prepare(args.source_dir, args.output)
