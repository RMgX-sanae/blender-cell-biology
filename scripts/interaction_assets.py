"""Bind prepared interaction geometry and motion to their preparation record."""
import hashlib
from pathlib import Path


def verify_assets(folder, manifest):
    folder = Path(folder)
    required = {'motion.npz', *(actor['file'] for actor in manifest['actors'])}
    hashes = manifest.get('data_sha256', {})
    if set(hashes) != required:
        raise ValueError('Prepared file hashes must cover every actor and motion.npz; rerun prepare_interaction.py')
    for name, expected in hashes.items():
        if Path(name).name != name:
            raise ValueError('Prepared actor files must be directly inside the asset directory')
        path = folder/name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError('Prepared interaction data missing or changed: '+name)
    return {'verified': True, 'files': len(required),
            'manifest_sha256': hashlib.sha256((folder/'manifest.json').read_bytes()).hexdigest()}
