"""Inspect effective two-cycle appearance without claiming catalytic dynamics."""
import argparse
import json
import sys
from pathlib import Path
import bpy
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_layered_appearance import equal, check_request, check_lighting, check_composition, check_typography
from audit_surface_materials import check_materials
from detail_scene import linear
from transcription_appearance import read_profile


def audit(scene):
    if 'transcription_appearance' not in scene:
        return None
    record = json.loads(scene['transcription_appearance'])
    check_request(scene, record, read_profile)
    manifest = json.loads(Path(bpy.data.filepath).with_name('manifest.json').read_text())
    cycles = 2 if 'next_cycle_geometry' in manifest else 1
    if record['cycles'] != cycles:
        raise ValueError('Appearance cycle count differs')
    p = record['profile']; studio = record['studio']; palette = studio['palette']
    mapping = {
        'Template DNA — all modeled residues': ('nucleic_acid', 'dna2'),
        'Template DNA except complete A16': ('nucleic_acid', 'dna2'),
        'Complete template A16 — modeled loading pose': ('nucleic_acid', 'dna2'),
        'Template T15 O3′ to T16 P': ('nucleic_acid', 'dna2'),
        'Template T16 O3′ to T17 P': ('nucleic_acid', 'dna2'),
        'Non-template DNA — hidden whole chain': ('nucleic_acid', 'dna1'),
        'RNA primer P residues 7–19': ('nucleic_acid', 'rna'),
        'Incoming CTP core → RNA C20': ('nucleic_acid', 'protein'),
        'Incoming UTP core → modeled RNA U21': ('nucleic_acid', 'protein'),
        'CTP beta/gamma phosphates → PPi': ('phosphate', 'accent'),
        'UTP beta/gamma phosphates → second PPi': ('phosphate', 'accent'),
        'Breaking alpha–beta phosphate connection': ('phosphate', 'accent'),
        'Second alpha–beta phosphate connection': ('phosphate', 'accent'),
        'New RNA O3′–P bond': ('nucleic_acid', 'rna'),
        'Second RNA O3′–P bond': ('nucleic_acid', 'rna'),
    }
    roles = set(); bodies = []
    for obj in scene.objects:
        if obj.type not in ('MESH', 'CURVE') or obj.name.startswith('Leader '):
            continue
        if obj.get('source_chain'):
            role, key = 'protein', 'heads'
            mods = list(obj.modifiers)
            if len(mods) != 1 or mods[0].type != 'SMOOTH' or mods[0].iterations != 2 or not mods[0].show_render:
                raise ValueError('Context source envelope modifier differs')
            equal(mods[0].factor, .35, 'Context smoothing differs')
            if obj.hide_render:
                raise ValueError('Whole Pol II chain unexpectedly hidden')
        else:
            if obj.name not in mapping:
                raise ValueError('Unidentified reaction display object: '+obj.name)
            role, key = mapping[obj.name]
            if obj.modifiers:
                raise ValueError('Substrate display modifier changed')
        mats = list(obj.data.materials)
        if len(mats) != 1 or mats[0] is None or mats[0].get('surface_role') != role:
            raise ValueError('Object role assignment differs: '+obj.name)
        if obj.type == 'MESH' and any(face.material_index != 0 for face in obj.data.polygons):
            raise ValueError('Face material assignment differs')
        socket = mats[0].node_tree.nodes['Principled BSDF'].inputs['Base Color']
        if socket.is_linked:
            raise ValueError('Semantic color overridden')
        if cycles == 1 or obj.name != 'Incoming CTP core → RNA C20':
            equal(socket.default_value, (*linear(palette[key]), 1), 'Semantic color differs: '+obj.name)
        roles.add(role); bodies.append(obj.name)
    check_materials(p, roles & p.get('surface_roles', {}).keys())
    target = np.array([-.1, -2.3, 1.2]); front = np.array([.9396926, .3420201, 0]); front /= np.linalg.norm(front)
    right = np.cross([0, 0, 1], front); up = np.cross(front, right)
    settings = p.get('camera', {}); pan = settings.get('center_offset_nm', [0, 0])
    shift = right*pan[0]+up*pan[1]
    camera = scene.camera
    if camera.data.type != 'ORTHO' or camera.data.dof.use_dof or camera.parent or camera.constraints:
        raise ValueError('Authored transcription camera changed')
    frame = scene.frame_current
    try:
        for f in range(scene.frame_start, scene.frame_end+1):
            scene.frame_set(f); bpy.context.view_layer.update()
            u = np.clip(((f-1)/24-11)/2.5, 0, 1) if cycles == 2 else 0
            advance = u*u*(3-2*u)
            equal(camera.location, target+shift+50*front+up*.30*advance, 'Requested camera path/pan differs')
            equal(camera.data.ortho_scale, 14.5*settings.get('scale_multiplier', 1)+1.5*advance, 'Requested camera scale differs')
            equal(np.array(camera.matrix_world.to_3x3()), np.column_stack([right, up, front]), 'Camera aim differs')
            if cycles == 2:
                core = scene.objects['Incoming CTP core → RNA C20']
                color = core.data.materials[0].node_tree.nodes['Principled BSDF'].inputs['Base Color']
                expected = (1-advance)*np.array([*linear(palette['protein']), 1])+advance*np.array([*linear(palette['rna']), 1])
                equal(color.default_value, expected, 'Incorporated C20 color handoff differs')
        scene.frame_set(manifest['preview_frame']); bpy.context.view_layer.update()
        equal(np.array(camera.matrix_world), record['camera_matrix'], 'Saved preview camera matrix differs')
    finally:
        scene.frame_set(frame); bpy.context.view_layer.update()
    check_lighting(scene, studio, target+shift)
    check_composition(scene, record, 'Whole Pol II context')
    def check_caption_frame(f):
        # Evaluate caption actions in their actual render layer, not the beauty
        # graph where they are deliberately excluded from artistic filtering.
        prefix = 'Cycle stage ' if cycles == 2 else 'Stage '
        time = (f-1)/24 if cycles == 2 else (f-1)/(scene.frame_end-1)
        start, end = ('start_seconds', 'end_seconds') if cycles == 2 else ('start', 'end')
        expected = {prefix+str(i+1) for i, row in enumerate(manifest['stages']) if row[start] <= time < row[end]}
        graph = bpy.context.evaluated_depsgraph_get()
        actual = {o.name for o in scene.objects if o.name.startswith(prefix) and not o.evaluated_get(graph).hide_render}
        if actual != expected:
            raise ValueError('Caption event timing differs at '+str(f))
    typography = check_typography(scene, record, 19 if cycles == 2 else 12, check_caption_frame)
    has_chinese = any(0x4e00 <= ord(c) <= 0x9fff for o in scene.objects if o.type == 'FONT' for c in o.data.body)
    if has_chinese != (record['language'] == 'zh'):
        raise ValueError('Caption language differs')
    return {'passed': True, 'roles': sorted(roles), 'objects': len(bodies), 'cycles': cycles, **typography,
            'scope': 'Requested shaders, lights, camera, two-cycle color handoff, context layers and protected bilingual captions. Geometry and actual pixels require their separate checks.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(sys.argv[sys.argv.index('--')+1:])
    result = audit(bpy.context.scene)
    args.output.write_text(json.dumps(result, indent=2)+'\n'); print(json.dumps(result))
