"""Bind prepared interaction geometry and motion to their preparation record."""
import hashlib
from pathlib import Path


def verify_assets(folder, manifest):
    folder = Path(folder)
    if 'binding_candidate' in manifest or manifest.get('input_kind') in ('separate_monomers','separate_partners'):
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
        config = manifest['source_config']
        def instance(spec):
            return (spec.get('instance_id',''),spec.get('label_chain',''),spec.get('operator',''))
        polymers = sorted((a['role'], a['chain'], *instance(a)) for a in manifest['actors'] if a['kind'] == 'polymer')
        expected = sorted((role, item['chain'], *instance(item)) for role in ['anchor','moving']
                          for item in config[role].get('chains',[config[role]]))
        if polymers != expected or len(polymers) != len(set(polymers)):
            raise ValueError('Retain every declared partner polymer instance exactly once')
        for role in ['anchor','moving']:
            if 'chains' not in config[role]: continue
            proof = manifest.get('partner_assembly_checks',{}).get(role,{})
            if (not proof.get('passed') or proof.get('reference_sha256') != config[role]['assembly_reference']['sha256']
                    or proof.get('matched_instances') != len(config[role]['chains'])+len(config[role]['ligand_instances'])
                    or set(proof.get('instance_map',{})) != {s['instance_id'] for s in config[role]['chains']+config[role]['ligand_instances']}):
                raise ValueError('Partner instances require their independent assembly check')
            if any(a.get('assembly_id') != config[role]['assembly_id'] for a in manifest['actors'] if a['role'] == role):
                raise ValueError('Partner lost its reviewed assembly identity')
        expected = sorted((role, item['chain'], item['component'], item['residue_number'], item.get('insertion_code', ''),
                           *instance(item),item.get('owner_instance',''))
                          for role in ['anchor', 'moving'] for item in config[role].get('ligand_instances', []))
        actual = sorted((a['role'], a['chain'], a['ccd'], a['residue_number'], a.get('insertion_code', ''),
                         *instance(a),a.get('owner_instance',''))
                        for a in manifest['actors'] if a['kind'] == 'ligand')
        if actual != expected:
            raise ValueError('Retain each declared monomer companion with its original partner ownership')
    ids = [a['id'].casefold() for a in manifest['actors']]
    files = [a['file'].casefold() for a in manifest['actors']]
    if len(ids) != len(set(ids)) or len(files) != len(set(files)):
        raise ValueError('Prepared actor ids and files must be unique')
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
