"""Saved-state checks shared by intact layered translation/transcription views."""
import hashlib
import json
from pathlib import Path
import bpy
import numpy as np
from mathutils import Vector
from bpy_extras.object_utils import world_to_camera_view
from detail_scene import linear
from export_composite_layers import published_sources
from layered_molecular_appearance import digest


def equal(a, b, message):
    if np.shape(a) != np.shape(b) or not np.allclose(a, b, atol=2e-6, rtol=1e-5):
        raise ValueError(message)


def check_request(scene, record, read_profile):
    manifest = json.loads(Path(bpy.data.filepath).with_name('manifest.json').read_text())
    if (manifest.get('appearance') != record or digest(record['profile_path']) != record['profile_sha256']
            or read_profile(record['profile_path']) != (record['profile'], True)):
        raise ValueError('Saved layered molecular appearance differs from immutable request')
    p = record['profile']; studio = record['studio']
    defaults = json.loads((Path(__file__).resolve().parents[1]/'assets/style-molecular-reference.json').read_text())
    expected = {**defaults, **p, 'palette': {**defaults['palette'], **p.get('palette', {})}, 'camera': {}}
    if studio != expected:
        raise ValueError('Requested layered molecular appearance ignored')


def check_lighting(scene, studio, target):
    lights = {o.name: o for o in scene.objects if o.type == 'LIGHT'}
    if set(lights) != {r['name'] for r in studio['lighting']}:
        raise ValueError('Light inventory differs')
    for row in studio['lighting']:
        obj = lights[row['name']]
        if obj.data.type != 'AREA' or obj.hide_render or obj.parent or obj.constraints or obj.animation_data or obj.data.animation_data:
            raise ValueError('Light display differs')
        equal(obj.location, row['position'], 'Light position differs'); equal(obj.data.energy, row['power'], 'Light power differs')
        equal(obj.data.size, row['size'], 'Light size differs'); equal(obj.data.color, linear(row['color']), 'Light color differs')
        direction = target-np.array(obj.location); direction /= np.linalg.norm(direction)
        equal(-np.array(obj.matrix_world.to_3x3())[:, 2], direction, 'Light aim differs')
    world = scene.world.node_tree.nodes['Background']
    equal(world.inputs['Color'].default_value, (*linear(studio['world']), 1), 'World color differs')
    equal(world.inputs['Strength'].default_value, .35, 'World strength differs')
    if (scene.view_settings.view_transform != studio.get('view_transform', 'AgX') or scene.view_settings.look !=
            studio.get('look', 'AgX - Medium High Contrast' if scene.view_settings.view_transform == 'AgX' else 'None')):
        raise ValueError('Display transform differs')
    equal(scene.view_settings.exposure, studio.get('exposure', 0), 'Exposure differs')


def check_composition(scene, record, context_name):
    p = record['profile']; studio = record['studio']
    tree = scene.compositing_node_group; ports = published_sources(tree)
    def node(label):
        found = [n for n in tree.nodes if n.label == label]
        if len(found) != 1: raise ValueError('Missing compositor node: '+label)
        return found[0]
    def linked(socket, expected, message):
        if len(socket.links) != 1 or socket.links[0].from_socket != expected: raise ValueError(message)
    beauty = scene.view_layers[scene['compositor_beauty_layer']]
    text = scene.view_layers[scene['compositor_label_layer']]
    context = scene.view_layers[record['context_layer']]
    if (not beauty.layer_collection.children[record['context_collection']].exclude or
            not beauty.layer_collection.children[scene['compositor_label_collection']].exclude):
        raise ValueError('Context or text leaks into sharp beauty')
    if {c.name for c in context.layer_collection.children if not c.exclude} != {record['context_collection'], '90 Camera and lighting'}:
        raise ValueError('Substrates or labels leak into context blur')
    if {c.name for c in text.layer_collection.children if not c.exclude} != {scene['compositor_label_collection']}:
        raise ValueError('Foreground leaks into protected typography')
    base = node('Straight background plus premultiplied rendered foreground')
    blur = node(context_name+' defocus'); over = node(context_name+' over background')
    equal(blur.inputs['Size'].default_value, [p.get('context_blur_pixels_1280', 0)*scene.render.resolution_x/1280]*2, 'Context blur differs')
    equal(over.inputs['Background'].default_value, (*linear(studio['background']), 1), 'Background differs')
    equal(over.inputs['Factor'].default_value, 1, 'Context opacity differs')
    if len(blur.inputs['Image'].links) != 1 or blur.inputs['Image'].links[0].from_node.type != 'R_LAYERS' or blur.inputs['Image'].links[0].from_node.layer != context.name:
        raise ValueError('Context render layer differs')
    linked(over.inputs['Foreground'], blur.outputs['Image'], 'Context blur bypassed')
    linked(base.inputs['Background'], over.outputs['Image'], 'Context composition bypassed')
    if ports['Labels'].node.type != 'R_LAYERS' or ports['Labels'].node.layer != text.name:
        raise ValueError('Text routed through artistic filters')
    if scene['selected_compositor'] != p.get('render_recipe', 'neutral'):
        raise ValueError('Render recipe differs')
    if scene['selected_compositor'] == 'neutral' and ports['Beauty'] != base.outputs['Image']:
        raise ValueError('Neutral foreground bypassed')
    output = next(n for n in tree.nodes if n.type == 'GROUP_OUTPUT' and n.is_active_output)
    final = output.inputs['Image'].links[0].from_node
    if final.type != 'ALPHAOVER': raise ValueError('Protected text composite bypassed')
    linked(final.inputs['Background'], ports['Beauty'], 'Beauty output bypassed')
    linked(final.inputs['Foreground'], ports['Labels'], 'Text output bypassed')


def check_typography(scene, record, expected_count, check_frame=None):
    studio = record['studio']; palette = studio['palette']
    text = scene.view_layers[scene['compositor_label_layer']]
    group = bpy.data.collections[scene['compositor_label_collection']]
    camera = scene.camera
    fonts = [o for o in scene.objects if o.type == 'FONT']
    if {o.name: o.data.body for o in fonts} != {r['object']: r['body'] for r in record['fonts']} or len(fonts) != expected_count:
        raise ValueError('Caption inventory differs')
    for obj in fonts:
        if set(obj.users_collection) != {group}:
            raise ValueError('Caption is outside protected typography: '+obj.name)
        if not obj.get('compositor_preserve_text_color'): raise ValueError('Text protection removed')
        for font, sha in zip([obj.data.font, obj.data.font_bold], record['font_sha256']):
            if not font.packed_file or hashlib.sha256(bytes(font.packed_file.data)).hexdigest() != sha:
                raise ValueError('Packed bilingual font differs')
        if any(fmt.use_bold != (ord(c) < 0x400 or 0x2010 <= ord(c) <= 0x22ff) for c, fmt in zip(obj.data.body, obj.data.body_format)):
            raise ValueError('Mixed font routing differs')
        mat = obj.data.materials[0]
        out = next(n for n in mat.node_tree.nodes if n.type == 'OUTPUT_MATERIAL' and n.is_active_output)
        if len(out.inputs['Surface'].links) != 1 or out.inputs['Surface'].links[0].from_node.type != 'EMISSION':
            raise ValueError('Protected text is not unlit')
        emit = out.inputs['Surface'].links[0].from_node
        light = np.dot(linear(studio['background']), [.2126, .7152, .0722]) > .45
        equal(emit.inputs['Color'].default_value, (*linear(palette.get('ink', '24323A' if light else 'E5EFF2')), 1), 'Text ink differs')
        equal(emit.inputs['Strength'].default_value, 1, 'Text strength differs')
    original = bpy.context.window.view_layer; frame = scene.frame_current; evaluated = 0
    try:
        bpy.context.window.view_layer = text
        for f in range(scene.frame_start, scene.frame_end+1):
            scene.frame_set(f); bpy.context.view_layer.update(); graph = bpy.context.evaluated_depsgraph_get()
            if check_frame:
                check_frame(f)
            for obj in fonts:
                eo = obj.evaluated_get(graph)
                if eo.hide_render: continue
                points = [world_to_camera_view(scene, camera, eo.matrix_world@Vector(p)) for p in eo.bound_box]
                if any(v.z <= 0 or min(v.x, v.y) < 0 or max(v.x, v.y) > 1 for v in points):
                    raise ValueError('Caption outside frame: '+obj.name+' at '+str(f))
                evaluated += 1
    finally:
        bpy.context.window.view_layer = original; scene.frame_set(frame); bpy.context.view_layer.update()
    return {'fonts': len(fonts), 'caption_evaluations': evaluated, 'frames': scene.frame_end}
