"""Bind prepared interaction geometry and motion to their preparation record."""
import hashlib
from pathlib import Path


def verify_assets(folder, manifest):
    folder = Path(folder)
    if 'binding_candidate' in manifest or manifest.get('input_kind') == 'separate_monomers':
        candidate = manifest.get('binding_candidate', {})
        if (candidate.get('evidence_class') != 'inferred'
                or candidate.get('pose_kind') not in {'template_guided_candidate', 'restraints_guided_candidate'}
                or not candidate.get('endpoint_screen', {}).get('passed')
                or candidate.get('endpoint_screen', {}).get('deep_overlap_atoms') != 0
                or manifest['motion'].get('strict_overlap_atoms_max') != 0
                or manifest['motion'].get('endpoint_source_overlap_pairs') != 0):
            raise ValueError('Binding hypotheses require explicit candidate evidence and strict non-exempt endpoint/path checks')
        if not manifest.get('display_title', '').startswith('Binding hypothesis:') or not manifest.get('animation_evidence'):
            raise ValueError('Keep a readable candidate label and its motion evidence')
        for actor in manifest['actors']:
            if actor.get('source_pdb') != manifest['source_config'][actor['role']]['pdb']:
                raise ValueError('Monomer actor lost its individual source identity')
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
