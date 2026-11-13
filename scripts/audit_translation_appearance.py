"""Inspect effective role materials, context composition and animated typography."""
import argparse
import json
import sys
from pathlib import Path

import bpy
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_surface_materials import check_materials
from detail_scene import linear
from translation_appearance import read_profile
from audit_layered_appearance import equal, check_request, check_lighting, check_composition, check_typography


def audit(scene):
    if 'translation_appearance' not in scene:
        return None
    record = json.loads(scene['translation_appearance'])
    check_request(scene, record, read_profile)
    manifest = json.loads(Path(bpy.data.filepath).with_name('manifest.json').read_text())
    p = record['profile']; studio = record['studio']; palette = studio['palette']
    source = {r['chain']: r for r in manifest['source_asset']['chains']}
    roles = set(); bodies = []
    for obj in scene.objects:
        if obj.type not in ('MESH', 'CURVE') or obj.name.startswith('Polarity guide '):
            continue
        chain = obj.get('source_chain')
        if chain:
            role = 'nucleic_acid' if source[chain]['nucleic_acid'] else 'protein'
            key = ('alpha' if chain.startswith('L') else 'beta') if role == 'nucleic_acid' else 'heads'
            mods = list(obj.modifiers)
            if (len(mods) != 1 or mods[0].type != 'SMOOTH' or mods[0].iterations != 2
                    or not mods[0].show_render or not mods[0].show_viewport):
                raise ValueError('Context source envelope modifier differs')
            equal(mods[0].factor, .35, 'Context smoothing differs')
            if bool(obj.hide_render) != (chain in manifest['presentation']['hidden_complete_chains']):
                raise ValueError('Whole subunit visibility differs')
        else:
            if obj.modifiers:
                raise ValueError('Substrate display modifier changed')
            role = 'nucleic_acid'
            if obj.get('source_nucleic_chain') in ('At', 'Pt', 'Et'):
                key = {'At': 'protein', 'Pt': 'rna', 'Et': 'accent'}[obj['source_nucleic_chain']]
            elif obj.get('source_nucleic_chain') == 'Mr':
                key = 'dna1' if obj['residue_number'] in (5, 6, 7) else 'dna2'
            elif obj.name.startswith('mRNA phosphodiester '):
                key = 'dna1'
            else:
                role = 'peptide_symbol'
                key = 'amino_acid' if obj.name in ('Peptide residue Ile', 'A-tRNA ester connection') else 'peptide'
        mats = list(obj.data.materials)
        if len(mats) != 1 or mats[0] is None or mats[0].get('surface_role') != role:
            raise ValueError('Object role assignment differs: '+obj.name)
        if obj.type == 'MESH' and any(face.material_index != 0 for face in obj.data.polygons):
            raise ValueError('Face material assignment differs')
        socket = mats[0].node_tree.nodes['Principled BSDF'].inputs['Base Color']
        if socket.is_linked:
            raise ValueError('Semantic color overridden')
        color = palette.get(key, {'peptide': 'F8E8CD', 'amino_acid': 'FE9E56'}.get(key))
        equal(socket.default_value, (*linear(color), 1), 'Semantic color differs: '+obj.name)
        roles.add(role); bodies.append(obj.name)
    if roles & p.get('surface_roles', {}).keys():
        check_materials(p, roles & p['surface_roles'].keys())
    camera = scene.camera
    target = np.array([-2.4, -.3, -2.2]); front = np.array([-.254887, -.951251, -.173648])
    position = target+50*front; front /= np.linalg.norm(front)
    right = np.cross([0, 0, 1], front); right /= np.linalg.norm(right); up = np.cross(front, right)
    settings = p.get('camera', {}); pan = settings.get('center_offset_nm', [0, 0]); shift = right*pan[0]+up*pan[1]
    if camera.data.type != 'ORTHO' or camera.data.dof.use_dof or camera.parent or camera.constraints or camera.animation_data or camera.data.animation_data:
        raise ValueError('Authored translation camera changed')
    equal(camera.location, position+shift, 'Requested camera pan ignored')
    equal(camera.data.ortho_scale, 26*settings.get('scale_multiplier', 1), 'Requested camera scale ignored')
    equal(np.array(camera.matrix_world), record['camera_matrix'], 'Camera matrix differs')
    equal(np.array(camera.matrix_world.to_3x3()), np.column_stack([right, up, front]), 'Camera aim differs')
    check_lighting(scene, studio, target+shift)
    check_composition(scene, record, 'Whole ribosomal context')
    typography = check_typography(scene, record, 12)
    return {'passed': True, 'roles': sorted(roles), 'objects': len(bodies), **typography,
            'scope': 'Actual appearance, isolated whole context and protected bilingual captions. Source geometry, motion and actual pixels are separately checked.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(sys.argv[sys.argv.index('--')+1:]); result = audit(bpy.context.scene)
    args.output.write_text(json.dumps(result, indent=2)+'\n'); print(json.dumps(result))
