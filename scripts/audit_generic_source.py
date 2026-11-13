"""Independently compare a generic mmCIF cache with raw atom_site rows.

This bypasses the preparation atom_table/Structure parser and checks every
selected chain, atom row and coordinate. Named presets/project selections need
their own selection checks; this entry deliberately accepts generic full sets.
"""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import gemmi
import numpy as np
from molecular_asset_io import chain_file


def audit(source, directory):
    m = json.loads((directory/'manifest.json').read_text())
    if m['preset'] != 'generic' or 'project_config' in m:
        raise ValueError('This raw-row audit requires the generic full-polymer selection')
    if hashlib.sha256(source.read_bytes()).hexdigest() != m['source_sha256']:
        raise ValueError('Raw source hash differs')
    block = gemmi.cif.read_file(str(source)).sole_block()
    entity = block.get_mmcif_category('_entity.')
    polymers = {key for key, kind in zip(entity['id'], entity['type']) if kind == 'polymer'}
    table = block.get_mmcif_category('_atom_site.'); selected = defaultdict(list)
    first_model = table['pdbx_PDB_model_num'][0]
    for i, chain in enumerate(table['auth_asym_id']):
        if (table['label_entity_id'][i] not in polymers or table['pdbx_PDB_model_num'][i] != first_model
                or table['type_symbol'][i].upper() in ('H', 'D') or float(table['occupancy'][i]) <= 0
                or table['label_alt_id'][i] not in (False, None, 'A')):
            continue
        selected[chain].append(i)
    entries = {row['chain']: row for row in m['chains']}
    if set(entries) != set(selected) or len(entries) != len(m['chains']):
        raise ValueError('Prepared chain inventory differs from raw polymers')
    origin = np.asarray(m['source_origin_angstrom']); basis = np.asarray(m['source_basis_columns'])
    points = np.column_stack([np.asarray(table['Cartn_'+axis], float) for axis in 'xyz'])
    rows = []
    for chain, indices in selected.items():
        with np.load(chain_file(directory, chain), allow_pickle=False) as data:
            expected = ((points[indices]-origin)@basis).astype(np.float32)*.1
            if not np.array_equal(expected, data['atoms']):
                raise ValueError('Raw atom coordinates differ: '+chain)
            for field, raw, transform in [('atom_names', 'auth_atom_id', str),
                                           ('residue_names', 'auth_comp_id', str),
                                           ('residue_numbers', 'auth_seq_id', int),
                                           ('elements', 'type_symbol', lambda s: gemmi.Element(s).name)]:
                values = np.array([transform(table[raw][i]) for i in indices])
                if not np.array_equal(values, data[field]):
                    raise ValueError('Raw atom row identity differs: '+chain+'/'+field)
            rows.append({'chain': chain, 'file': chain_file(directory, chain).name,
                         'selected_atoms': len(indices), 'exact_coordinates_and_rows': True})
    return {'passed': True, 'source': str(source.resolve()), 'source_sha256': m['source_sha256'],
            'manifest_sha256': hashlib.sha256((directory/'manifest.json').read_bytes()).hexdigest(),
            'chains': rows, 'chain_count': len(rows), 'atoms': sum(r['selected_atoms'] for r in rows),
            'scope': 'All generic selected polymer atom rows/coordinates compared with raw mmCIF categories, independent of the preparation atom_table/Structure path. Not a check of missing residues, complete chemical assembly or biological dynamics.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True); p.add_argument('--asset-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); result = audit(a.source, a.asset_dir)
    a.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'chains'}))
