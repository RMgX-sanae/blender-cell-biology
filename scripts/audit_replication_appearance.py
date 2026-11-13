"""Check the effective materials of linked fork units and their layered display."""
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
from replication_appearance import read_profile


def audit(scene):
    if 'replication_appearance' not in scene:
        return None
    record = json.loads(scene['replication_appearance'])
    check_request(scene, record, read_profile)
    p = record['profile']; studio = record['studio']; palette = studio['palette']
    roles = set(); count = 0; object_slots = 0; guides = 0
    for obj in scene.objects:
        if obj.type not in ('MESH', 'CURVE') or obj.name.startswith(('Prototype ', 'Primer symbol prototype ')):
            continue
        if obj.get('source_chain'):
            role, key = 'protein', 'heads'
            mods = list(obj.modifiers)
            if len(mods) != 1 or mods[0].type != 'SMOOTH' or mods[0].iterations != 2 or not mods[0].show_render:
                raise ValueError('Context source envelope modifier differs')
            equal(mods[0].factor, .35, 'Context smoothing differs')
            if obj.hide_render: raise ValueError('Whole CMG chain unexpectedly hidden')
        else:
            if obj.get('strand_role'):
                if len(obj.material_slots) != 1 or obj.material_slots[0].link != 'OBJECT':
                    raise ValueError('Object material slot link differs: '+obj.name)
                object_slots += 1
                strand = obj['strand_role']
                key = 'accent' if strand == 'lagging' and int(obj['index']) % 12 >= 9 else {
                    'parent_A':'dna2', 'parent_B':'dna2', 'leading':'protein', 'lagging':'dna1'}[strand]
            elif obj.name.startswith(('parent_', 'Parent color key')): key = 'dna2'
            elif obj.name.startswith(('leading ', 'Leading ')): key = 'protein'
            elif obj.name.startswith(('lagging ', 'Lagging fragment ')): key = 'dna1'
            elif obj.name == 'Primer color key': key = 'accent'
            else: raise ValueError('Unidentified fork display object: '+obj.name)
            role = 'primer_symbol' if key == 'accent' else 'nucleic_acid'
            if obj.modifiers: raise ValueError('Fork unit display modifier changed')
        # Linked nucleotide meshes have one DATA material but different OBJECT
        # overrides. Only effective object slots describe the rendered strand.
        mats = [slot.material for slot in obj.material_slots]
        if obj.type == 'CURVE' and obj.name.startswith(('Parent color key', 'Primer color key', 'Leading ', 'Lagging fragment ')):
            if len(mats) != 1 or mats[0] is None: raise ValueError('Missing annotation line material')
            output = next(n for n in mats[0].node_tree.nodes if n.type == 'OUTPUT_MATERIAL' and n.is_active_output)
            links = output.inputs['Surface'].links
            if len(links) != 1 or links[0].from_node.type != 'EMISSION':
                raise ValueError('Direction/color guide must be unlit: '+obj.name)
            emit = links[0].from_node
            equal(emit.inputs['Color'].default_value, (*linear(palette[key]), 1), 'Guide semantic color differs')
            equal(emit.inputs['Strength'].default_value, 1, 'Guide strength differs')
            guides += 1; continue
        if len(mats) != 1 or mats[0] is None or mats[0].get('surface_role') != role:
            raise ValueError('Object role assignment differs: '+obj.name)
        if obj.type == 'MESH' and any(face.material_index != 0 for face in obj.data.polygons):
            raise ValueError('Face material assignment differs')
        socket = mats[0].node_tree.nodes['Principled BSDF'].inputs['Base Color']
        if socket.is_linked: raise ValueError('Semantic color overridden')
        equal(socket.default_value, (*linear(palette[key]), 1), 'Semantic color differs: '+obj.name)
        roles.add(role); count += 1
    if object_slots != 288: raise ValueError('Fork nucleotide identity inventory differs')
    if guides != 16: raise ValueError('Direction/color guide inventory differs')
    check_materials(p, roles & p.get('surface_roles', {}).keys())
    front = np.array([0, -.35, 1.]); front /= np.linalg.norm(front)
    right = np.array([1., 0, 0]); up = np.cross(front, right)
    settings = p.get('camera', {}); pan = settings.get('center_offset_nm', [0, 0])
    shift = right*pan[0]+up*pan[1]; camera = scene.camera
    if (camera.data.type != 'ORTHO' or camera.data.dof.use_dof or camera.parent or camera.constraints
            or camera.animation_data or camera.data.animation_data):
        raise ValueError('Authored replication camera changed')
    equal(camera.location, shift+65*front, 'Requested camera pan differs')
    equal(camera.data.ortho_scale, 33*settings.get('scale_multiplier', 1), 'Requested camera scale differs')
    equal(np.array(camera.matrix_world.to_3x3()), np.column_stack([right, up, front]), 'Camera aim differs')
    equal(np.array(camera.matrix_world), record['camera_matrix'], 'Saved preview camera matrix differs')
    check_lighting(scene, studio, shift)
    check_composition(scene, record, 'Whole CMG context')
    def check_caption_frame(f):
        t = (f-1)/(scene.frame_end-1)
        u = np.clip((t-.06)/(.88-.06), 0, 1); fork = 16+35*u*u*(3-2*u)
        graph = bpy.context.evaluated_depsgraph_get()
        for block in range(6):
            for end in ('5′', '3′'):
                obj = scene.objects[f'Lagging fragment {block+1} {end}'].evaluated_get(graph)
                if obj.hide_render != (fork < block*12+16):
                    raise ValueError('Fragment polarity caption timing differs at '+str(f))
    typography = check_typography(scene, record, 23, check_caption_frame)
    has_chinese = any(0x4e00 <= ord(c) <= 0x9fff for o in scene.objects if o.type == 'FONT' for c in o.data.body)
    if has_chinese != (record['language'] == 'zh'): raise ValueError('Caption language differs')
    return {'passed':True, 'roles':sorted(roles), 'objects':count, 'object_material_slots':object_slots, 'unlit_guides':guides,
            **typography, 'scope':'Effective object material overrides, whole CMG layer, camera, lights, protected bilingual text and fragment caption timing. Source geometry, fork growth and actual pixels require separate checks.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(sys.argv[sys.argv.index('--')+1:]); result = audit(bpy.context.scene)
    args.output.write_text(json.dumps(result, indent=2)+'\n'); print(json.dumps(result))
