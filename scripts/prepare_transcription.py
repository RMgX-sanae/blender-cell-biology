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
from molecular_asset_io import chain_file
from reaction_ledger import validate


def addition_ledger(rna, primer_residue, component, names, elements, bonds,
                    product_names, product_residue, *, product_observed):
    """Local RNA-end + NTP heavy-atom graph; not a complete catalytic mechanism."""
    names = list(map(str, names)); product_names = list(map(str, product_names))
    core = len(product_names)
    leaving = {'O3A','PB','O1B','O2B','O3B','PG','O1G','O2G','O3G'}
    rename = {'PA':'P','O1A':'OP1','O2A':'OP2'}
    if (len(set(names)) != len(names) or set(names[core:]) != leaving or
            product_names != [rename.get(n,n) for n in names[:core]]):
        raise ValueError('Qualify the RNA core / PPi partition and CCD-to-RNA atom names')
    indices = np.flatnonzero(rna['residue_numbers'] == primer_residue)
    primer_ids = {int(i):f'RNA:P:{primer_residue}:'+str(rna['atom_names'][i]) for i in indices}
    primer_atoms = [[primer_ids[i],str(rna['elements'][i])] for i in indices]
    nucleophile = f"RNA:P:{primer_residue}:O3'"
    if dict(primer_atoms).get(nucleophile) != 'O':
        raise ValueError('The selected RNA end must contain a resolved 3-prime oxygen')
    primer_bonds = [[primer_ids[i],primer_ids[j]] for i,j in rna['bonds'] if i in primer_ids and j in primer_ids]
    incoming = [f'{component}_for_P{product_residue}:'+n for n in names]
    products = [f'RNA:P:{product_residue}:'+n for n in product_names]
    products += [f'PPi_from_P{product_residue}:'+n for n in names[core:]]
    product_elements = list(map(str,elements))
    if product_observed:
        ids = np.flatnonzero(rna['residue_numbers'] == product_residue)
        native = {str(rna['atom_names'][i]):int(i) for i in ids}
        if set(native) != set(product_names) or len(ids) != core:
            raise ValueError('Observed RNA product does not match the selected nucleotide core')
        product_elements[:core] = [str(rna['elements'][native[n]]) for n in product_names]
        native_ids = {i:f'RNA:P:{product_residue}:'+str(rna['atom_names'][i]) for i in ids}
        core_bonds = [[native_ids[i],native_ids[j]] for i,j in rna['bonds'] if i in native_ids and j in native_ids]
        all_ids = {**primer_ids, **native_ids}
        joining = [[all_ids[i],all_ids[j]] for i,j in rna['bonds'] if i in all_ids and j in all_ids and
                   ((i in primer_ids) != (j in primer_ids))]
    else:
        core_bonds = [[products[i],products[j]] for i,j in bonds if i < core and j < core]
        joining = [[nucleophile,f'RNA:P:{product_residue}:P']]
    ppi_bonds = [[products[i],products[j]] for i,j in bonds if i >= core and j >= core]
    book = {'version':1,
        'scope':f'Local terminal RNA P:{primer_residue} plus {component}; new RNA P:{product_residue} and PPi. Other RNA residues, polymerase, DNA, solvent, metals and hydrogens omitted.',
        'bond_scope':'Resolved primer intranucleotide graph; CCD NTP graph; product RNA core and PPi graph plus the 3-prime-O/alpha-P link. Upstream RNA linkage is outside this local subgraph. Bond orders/protonation are not checked.',
        'reactant_atoms':primer_atoms+[[identity,str(e)] for identity,e in zip(incoming,elements)],
        'product_atoms':primer_atoms+[[identity,e] for identity,e in zip(products,product_elements)],
        'atom_map':[[identity,identity] for identity,_ in primer_atoms]+list(map(list,zip(incoming,products))),
        'reactant_bonds':primer_bonds+[[incoming[i],incoming[j]] for i,j in bonds],
        'product_bonds':primer_bonds+core_bonds+ppi_bonds+joining,
        'break':[[f'{component}_for_P{product_residue}:PA',f'{component}_for_P{product_residue}:O3A']],
        'form':[[nucleophile,f'{component}_for_P{product_residue}:PA']],
        'product_core_evidence':'Observed 5FLM RNA coordinates and inferred covalent graph' if product_observed else 'Modeled CCD-derived RNA core; not a deposited product',
        'sources':['https://www.rcsb.org/structure/5FLM','https://www.rcsb.org/ligand/'+component,
                   'https://www-ssrl.slac.stanford.edu/research/highlights_archive/pol-II.pdf']}
    return book, validate(book)


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
    rna = np.load(chain_file(args.asset_dir, 'P'))
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
    book, receipt = addition_ledger(rna, 19, 'CTP', names, [ideal[n][0] for n in names], bonds,
        [mapping.get(n,n) for n in common], 20, product_observed=True)
    args.output.mkdir(parents=True)
    ledger_file = args.output/'reaction-ledger-ctp.json'
    ledger_file.write_text(json.dumps(book,indent=2)+'\n')
    np.savez_compressed(args.output/'ctp.npz', atoms=coords.astype(np.float32),
        atom_names=np.array(names), elements=np.array([ideal[n][0] for n in names]),
        bonds=bonds, bond_orders=np.array(orders),
        product_atom_count=len(common), product_names=np.array([mapping.get(n,n) for n in common]),
        product_atoms=product_coords.astype(np.float32))
    manifest = {'source_asset': source,
        'reaction_ledger':{'file':ledger_file.name,'sha256':hashlib.sha256(ledger_file.read_bytes()).hexdigest(),'audit':receipt},
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
