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


def snapshot_handoff(manifest, poses, selector, visible, manifest_sha256):
    """Describe one frame of a verified rigid-partner interaction in its prepared nm coordinates."""
    import numpy as np

    if manifest.get('units') != 'nm':
        raise ValueError('Interaction handoffs require prepared coordinates in nm')
    count = manifest['frames']
    if poses.shape not in {(count, 4, 4), (count + 1, 4, 4)}:
        raise ValueError('Expected a partner transform for each playback frame and an optional virtual seam')
    identity = np.eye(4)
    # Coordinates in the actor files already describe the prepared endpoint.
    # Only numerical roundoff is allowed when recognizing that placement.
    at_endpoint = np.all(np.isclose(poses, identity, rtol=0, atol=1e-10), axis=(1, 2))
    if selector == 'endpoint':
        matches = np.flatnonzero(at_endpoint[:count])
        if not len(matches):
            raise ValueError('The prepared endpoint is absent from the playback frames')
        frame = int(matches[np.argmin(np.abs(matches - count / 2))]) + 1
    else:
        try:
            frame = int(selector)
        except (TypeError, ValueError) as exc:
            raise ValueError('Choose a 1-based frame number or endpoint') from exc
        if not 1 <= frame <= len(poses):
            raise ValueError('Snapshot frame is outside the prepared poses')
    matrix = poses[frame - 1]
    rotation = matrix[:3, :3]
    if (not np.isfinite(matrix).all() or not np.allclose(matrix[3], identity[3], atol=1e-7, rtol=0)
            or not np.allclose(rotation.T @ rotation, identity[:3, :3], atol=1e-7, rtol=0)
            or np.linalg.det(rotation) <= 0):
        raise ValueError('The selected partner pose must preserve lengths and handedness')

    example = manifest['example']
    candidate = manifest.get('binding_candidate')
    source = manifest.get('source') or {}
    if candidate:
        endpoint = {'evidence': 'inferred', 'kind': candidate['pose_kind'], 'candidate_id': candidate.get('id'),
                    'candidate_report_sha256': candidate.get('report_sha256')}
    elif (manifest.get('input_kind') in (None, 'deposited_complex')
          and source.get('evidence') == 'experimental_structure_model' and source.get('sha256')
          and str(source.get('id', '')).casefold() == example['pdb'].casefold()
          and all(a.get('source_pdb', source['id']).casefold() == source['id'].casefold()
                  for a in manifest['actors'])):
        endpoint = {'evidence': 'deposited', 'kind': 'selected_coordinates'}
    else:
        endpoint = {'evidence': 'unclassified', 'kind': 'prepared_coordinates'}

    config = manifest.get('source_config', {})
    partners = {}
    for role in ('anchor', 'moving'):
        spec = config[role] if candidate else source
        partners[role] = {'label': example[role + '_label'], 'source_id': spec.get('pdb', spec.get('id')),
                          'source_sha256': spec.get('sha256'), 'assembly': spec.get('assembly_id', spec.get('assembly')),
                          'selection_note': spec.get('selection_note'),
                          'matrix_from_prepared_nm': (matrix if role == 'moving' else identity).tolist()}
    if candidate and candidate['pose_kind'] == 'template_guided_candidate' and config.get('template'):
        template = config['template']
        endpoint['template_sources'] = {}
        for role in partners:
            overrides = template[role] if 'chains' in config[role] else {}
            endpoint['template_sources'][role] = {'source_id': overrides.get('pdb', template.get('pdb')),
                                                  'source_sha256': overrides.get('sha256', template.get('sha256'))}
    identity_fields = ('id', 'role', 'kind', 'chain', 'file', 'source_pdb', 'instance_id', 'label_chain',
                       'operator', 'assembly_id', 'owner_instance', 'ccd', 'residue_number', 'insertion_code')
    actors = [{**{k: a[k] for k in identity_fields if k in a}, 'visible': a['id'] in visible}
              for a in manifest['actors']]
    matched = bool(at_endpoint[frame - 1])
    handoff = {'prepared_manifest_sha256': manifest_sha256, 'frame': frame, 'virtual_seam': frame > count,
               'at_prepared_endpoint': matched, 'endpoint': endpoint,
               'relative_pose_evidence': endpoint['evidence'] if matched else 'schematic',
               'motion_evidence': 'schematic', 'partners': partners, 'actors': actors,
               'interaction_scope': config.get('scope', example.get('interaction')), 'limits': example['limits']}
    if manifest.get('membrane'):
        membrane = manifest['membrane']
        fields = ('kind', 'origin_nm', 'normal', 'map_url', 'map_sha256', 'cargo_tether_atom', 'cargo_lipid_index', 'limits')
        handoff['membrane'] = {key: membrane[key] for key in fields if key in membrane}
    return handoff
