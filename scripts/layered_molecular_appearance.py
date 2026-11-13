"""Shared appearance for intact layered molecular teaching scenes."""
import hashlib
import json
import math
import re
from pathlib import Path

import bpy
import numpy as np
from detail_scene import material
from interaction_appearance import build_compositor


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_profile(path, *, roles, palette, mechanism):
    p = json.loads(Path(path).read_text())
    if not ({'surface_roles', 'render_recipe', 'context_blur_pixels_1280'} & p.keys()):
        return p, False
    allowed = {'name', 'reference', 'reference_page', 'usage', 'observations',
               'implementation_assumptions', 'palette', 'surface_roles', 'roughness',
               'bump', 'background', 'world', 'lighting', 'view_transform', 'look',
               'exposure', 'camera', 'render_recipe', 'context_blur_pixels_1280'}
    if set(p)-allowed:
        raise ValueError('Unsupported '+mechanism+' appearance fields: '+str(set(p)-allowed)+
                         '; geometry is fixed to intact surface context and atomic substrates')
    if set(p.get('palette', {}))-palette or set(p.get('surface_roles', {}))-roles:
        raise ValueError('Unknown '+mechanism+' semantic palette key or surface role')
    colors = list(p.get('palette', {}).values())+[p[k] for k in ('background', 'world') if k in p]
    if any(not isinstance(c, str) or not re.fullmatch(r'#?[0-9a-fA-F]{6}', c) for c in colors):
        raise ValueError('Appearance colors must be sRGB hex')
    camera = p.get('camera', {})
    if set(camera)-{'scale_multiplier', 'center_offset_nm'}:
        raise ValueError('Authored '+mechanism+' view supports only orthographic scale/pan; other views require path and label reauthoring')
    scale = camera.get('scale_multiplier', 1)
    pan = camera.get('center_offset_nm', [0, 0])
    if np.shape(scale) != () or not math.isfinite(scale) or scale <= 0:
        raise ValueError('Camera scale must be positive and finite')
    if np.shape(pan) != (2,) or not np.isfinite(pan).all():
        raise ValueError('Camera pan must be a finite right/up pair')
    blur = p.get('context_blur_pixels_1280', 0)
    if np.shape(blur) != () or not math.isfinite(blur) or blur < 0:
        raise ValueError('Context blur must be finite and nonnegative')
    if 'exposure' in p and not math.isfinite(p['exposure']):
        raise ValueError('Exposure must be finite')
    if 'lighting' in p:
        if not p['lighting'] or len({v['name'] for v in p['lighting']}) != len(p['lighting']):
            raise ValueError('Provide uniquely named area lights')
        for light in p['lighting']:
            if set(light) != {'name', 'position', 'power', 'color', 'size'}:
                raise ValueError('Area lights require name/position/power/color/size')
            if np.shape(light['position']) != (3,) or not np.isfinite(light['position']).all():
                raise ValueError('Light position must be finite')
            if any(not math.isfinite(light[k]) or light[k] <= 0 for k in ('power', 'size')):
                raise ValueError('Light power and size must be positive')
            if not re.fullmatch(r'#?[0-9a-fA-F]{6}', light['color']):
                raise ValueError('Invalid light color')
    return p, True


def shade(profile, name, color, role, roughness, bump):
    mat = material(name, color, roughness, bump,
                   settings=profile.get('surface_roles', {}).get(role))
    mat['surface_role'] = role
    return mat


def apply_fonts(scene, font_cjk, font_latin):
    latin = font_latin if isinstance(font_latin, bpy.types.VectorFont) else bpy.data.fonts.load(str(font_latin), check_existing=True)
    for obj in scene.objects:
        if obj.type != 'FONT':
            continue
        if font_cjk:
            obj.data.font = font_cjk if isinstance(font_cjk, bpy.types.VectorFont) else bpy.data.fonts.load(str(font_cjk), check_existing=True)
        else:
            obj.data.font = latin
        obj.data.font_bold = latin
        for char, fmt in zip(obj.data.body, obj.data.body_format):
            fmt.use_bold = ord(char) < 0x400 or 0x2010 <= ord(char) <= 0x22ff
        obj['compositor_preserve_text_color'] = True


def finish(scene, args, request, profile, context, *, record_key, context_name, scope):
    scene.frame_set(args.frame)
    bpy.context.view_layer.update()
    apply_fonts(scene, args.font, args.font_latin)
    compositor = build_compositor(scene, request, profile['background'], args.width)
    beauty = scene.view_layers[scene['compositor_beauty_layer']]
    beauty.layer_collection.children[context.name].exclude = True
    layer = scene.view_layers.new(context_name)
    for child in layer.layer_collection.children:
        child.exclude = child.name not in {context.name, '90 Camera and lighting'}
    tree = scene.compositing_node_group
    base = next(n for n in tree.nodes if n.label == 'Straight background plus premultiplied rendered foreground')
    render = tree.nodes.new('CompositorNodeRLayers'); render.layer = layer.name
    blur = tree.nodes.new('CompositorNodeBlur'); blur.label = context_name+' defocus'
    blur.inputs['Size'].default_value = [request.get('context_blur_pixels_1280', 0)*args.width/1280]*2
    tree.links.new(render.outputs['Image'], blur.inputs['Image'])
    over = tree.nodes.new('CompositorNodeAlphaOver'); over.label = context_name+' over background'
    over.inputs['Background'].default_value = base.inputs['Background'].default_value
    tree.links.new(blur.outputs['Image'], over.inputs['Foreground'])
    tree.links.new(over.outputs['Image'], base.inputs['Background'])
    snapshot = args.output.resolve()/'appearance.json'
    fonts = [{'object': o.name, 'body': o.data.body} for o in scene.objects if o.type == 'FONT']
    record = {'profile': request, 'profile_path': str(snapshot),
              'profile_sha256': digest(args.style_profile), 'studio': profile,
              'context_collection': context.name, 'context_layer': layer.name,
              'camera_matrix': np.array(scene.camera.matrix_world).tolist(),
              'compositor': compositor, 'fonts': fonts,
              'font_sha256': [digest(args.font or args.font_latin), digest(args.font_latin)],
              'scope': scope}
    scene[record_key] = json.dumps(record)
    return record
