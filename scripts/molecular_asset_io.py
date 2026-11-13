"""Resolve source chain identities without relying on filesystem case sensitivity."""
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
import numpy as np


def chain_filenames(chains):
    names = list(chains)
    if len(names) != len(set(names)):
        raise ValueError('Duplicate source chain identity')
    counts = Counter(name.casefold() for name in names)
    return {name: (f'chain_{name}.npz' if counts[name.casefold()] == 1 and
                   re.fullmatch(r'[A-Za-z0-9_-]+', name) else
                   'chain-id-'+name.encode('utf-8').hex()+'.npz') for name in names}


def chain_file(directory, chain):
    root = Path(directory)
    manifest = json.loads((root/'manifest.json').read_text())
    rows = manifest.get('geometry') if manifest.get('chains') and isinstance(manifest['chains'][0], str) else manifest['chains']
    entries = {row['chain']: row for row in rows}
    if len(entries) != len(rows) or chain not in entries:
        raise ValueError('Missing or duplicate source chain identity: '+chain)
    filenames = {name: row.get('file', f'chain_{name}.npz') for name, row in entries.items()}
    if any(Path(name).name != name or not name.endswith('.npz') for name in filenames.values()):
        raise ValueError('Chain assets require simple NPZ filenames')
    if len({name.casefold() for name in filenames.values()}) != len(filenames):
        raise ValueError('Ambiguous chain filenames; rebuild from raw source with case-safe chain identities')
    row = entries[chain]; path = root/filenames[chain]
    if manifest.get('chain_file_schema') == 1:
        if hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']:
            raise ValueError('Source chain file hash differs: '+chain)
    with np.load(path, allow_pickle=False) as data:
        if manifest.get('chain_file_schema') == 1 and (data['source_chain'].shape != () or data['source_chain'].item() != chain):
            raise ValueError('Source chain payload identity differs: '+chain)
        for key in ('atoms', 'vertices', 'faces'):
            if key in row and len(data[key]) != row[key]:
                raise ValueError('Source chain inventory differs: '+chain+'/'+key)
        if 'centroid_nm' in row and not np.allclose(data['atoms'].mean(0), row['centroid_nm'], atol=2e-5, rtol=1e-5):
            raise ValueError('Source chain coordinates differ from declared centroid: '+chain)
    return path
