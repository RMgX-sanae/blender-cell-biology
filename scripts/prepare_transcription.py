"""Prepare an explicitly reconstructed one-C addition example from 5FLM + CCD CTP.

Run outside Blender with numpy and gemmi. The product RNA geometry is observed;
the incoming substrate pose, phosphate appendage and motion are constructions.
"""
import argparse
import hashlib
import json
from pathlib import Path

import gemmi
import numpy as np


def fit_torsions(coords, names, product, extras):
    """Rotate four explicit acyclic bonds; preserve the CCD bond lengths/angles."""
    base=['N1','C2','N3','C4','C5','C6','O2','N4']
    axes=[('N1',"C1'",[n for n in names if n not in base]),
          ("C4'","C5'",["C5'","O5'",'PA','O1A','O2A']+extras),
          ("C5'","O5'",["O5'",'PA','O1A','O2A']+extras),
          ("O5'",'PA',['PA','O1A','O2A']+extras)]
    rotations=[]
    for _ in range(3):
        for a,b,moving in axes:
            ids=[names.index(n) for n in moving]
            origin=coords[names.index(b)]
            axis=origin-coords[names.index(a)];axis/=np.linalg.norm(axis)
            v=coords[ids]-origin
            best=(float(np.sum((coords[:len(product)]-product)**2)),0,coords)
            for degrees in range(-180,180):
                theta=np.deg2rad(degrees)
                candidate=coords.copy()
                candidate[ids]=origin+v*np.cos(theta)+np.cross(axis,v)*np.sin(theta)+np.outer(v@axis,axis)*(1-np.cos(theta))
                error=float(np.sum((candidate[:len(product)]-product)**2))
                if error<best[0]:best=(error,degrees,candidate)
            coords=best[2]
            rotations.append({'axis':[a,b],'degrees':best[1]})
    return coords,rotations


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir', type=Path, required=True)
    p.add_argument('--ccd', type=Path, required=True, help='RCSB Chemical Component Dictionary CTP.cif')
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.output.exists(): p.error('Choose a fresh output directory.')
    source = json.loads((args.asset_dir/'manifest.json').read_text())
    if not source['source_file'].upper().startswith('5FLM'):
        p.error('The residue map is specific to 5FLM.')
    rna = np.load(args.asset_dir/'chain_P.npz')
    mask = rna['residue_numbers'] == 20
    product = {str(n): np.array(x, dtype=float) for n, x in
               zip(rna['atom_names'][mask], rna['atoms'][mask])}
    if set(rna['residue_names'][mask]) != {'C'} or "O3'" not in product:
        p.error('Expected a modeled cytidine with a 3-prime oxygen at P:20.')
    block = gemmi.cif.read_file(str(args.ccd)).sole_block()
    if block.name != 'CTP': p.error('Expected CCD CTP.')
    rows = block.find('_chem_comp_atom.', ['atom_id', 'type_symbol',
        'pdbx_model_Cartn_x_ideal', 'pdbx_model_Cartn_y_ideal', 'pdbx_model_Cartn_z_ideal'])
    unquote = gemmi.cif.as_string
    ideal = {unquote(r[0]): (unquote(r[1]), np.array([float(r[k]) for k in (2,3,4)])*.1)
             for r in rows if r[1] != 'H'}
    mapping = {'PA': 'P', 'O1A': 'OP1', 'O2A': 'OP2'}
    common = [n for n in ideal if mapping.get(n, n) in product]
    extras = [n for n in ideal if n not in common]
    # Align the complementary base, keeping the ideal CTP molecule rigid.
    # Fitting the flexible triphosphate together with the base misorients the
    # base and can create collisions with the template.
    fitting = ['N1','C2','N3','C4','C5','C6','O2','N4']
    x = np.array([ideal[n][1] for n in fitting])
    y = np.array([product[mapping.get(n, n)] for n in fitting])
    cx, cy = x.mean(0), y.mean(0)
    u, _, vt = np.linalg.svd((x-cx).T @ (y-cy))
    if np.linalg.det(vt.T @ u.T) < 0: vt[-1] *= -1
    rot = vt.T @ u.T
    shift = cy - rot @ cx
    rmsd = float(np.sqrt(np.mean(np.sum((x @ rot.T + shift-y)**2, axis=1))))
    names = common + extras
    coords = np.array([rot @ ideal[n][1]+shift for n in names])
    product_coords=np.array([product[mapping.get(n,n)] for n in common])
    coords,torsions=fit_torsions(coords,names,product_coords,extras)
    index = {n: i for i, n in enumerate(names)}
    bond_rows = block.find('_chem_comp_bond.', ['atom_id_1', 'atom_id_2', 'value_order'])
    bonds, orders = [], []
    for a, b, order in bond_rows:
        a, b = unquote(a), unquote(b)
        if a in index and b in index:
            bonds.append([index[a],index[b]])
            orders.append(order)
    bonds = np.array(bonds, dtype=np.int32)
    lengths = np.linalg.norm(coords[bonds[:,0]]-coords[bonds[:,1]],axis=1)
    if not (.08 < lengths.min() and lengths.max() < .2):
        raise ValueError('Reconstructed CTP contains an excessive bond-length distortion.')
    args.output.mkdir(parents=True)
    np.savez_compressed(args.output/'ctp.npz', atoms=coords.astype(np.float32),
        atom_names=np.array(names), elements=np.array([ideal[n][0] for n in names]),
        bonds=bonds, bond_orders=np.array(orders),
        product_atom_count=len(common), product_names=np.array([mapping.get(n,n) for n in common]),
        product_atoms=product_coords.astype(np.float32))
    manifest = {'source_asset': source,
        'ccd': {'component': 'CTP', 'url': 'https://files.rcsb.org/ligands/download/CTP.cif',
                'sha256': hashlib.sha256(args.ccd.read_bytes()).hexdigest(),
                'coordinates': 'ideal CTP with four acyclic bond torsions fitted to the product; not an experimentally bound substrate'},
        'product': {'chain':'P','residue':20,'name':'C','source_atoms_retained':len(common)},
        'primer': {'chain':'P','last_residue':19,'name':'G', 'nucleophile':"O3'"},
        'template_partner': {'chain':'T','residue':17,'name':'DG'},
        'constructed_atoms':extras, 'cleaved_ccd_bond':['PA','O3A'],
        'new_bond':["P:19:O3'",'P:20:P'],
        'base_fit_atoms':fitting,'base_fit_rmsd_nm':rmsd,
        'fitted_torsions':torsions,
        'full_core_endpoint_rmsd_nm':float(np.sqrt(np.mean(np.sum((coords[:len(common)]-product_coords)**2,axis=1)))),
        'bond_length_range_nm':[float(lengths.min()),float(lengths.max())],
        'interpretation':'Product-coordinate teaching reconstruction, NOT an observed CTP-bound Pol II state.',
        'limitations':['5FLM supplies a static bovine Pol II complex with synthetic nucleic-acid scaffold.',
            'Ideal CTP is fitted to the product and offset before joining; interpolation to the product is illustrative, not catalytic geometry.',
            'No hydrogen, protonation, metal chemistry, trigger-loop closure or translocation is simulated.',
            'Missing source residues remain absent. Do not connect across the N-chain gap.']}
    (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:manifest[k] for k in ['product','constructed_atoms','base_fit_rmsd_nm','full_core_endpoint_rmsd_nm','bond_length_range_nm']}))


if __name__ == '__main__': main()
