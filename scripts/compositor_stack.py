"""Editable, ordered artistic filters with explicit render-pass requirements.

These nodes change presentation, not scientific data or molecular geometry.
Validated combinations must be recorded per project and Blender version.
"""
import bpy
from mathutils import Vector

from detail_scene import linear
from style_nodes import node, math_node, ramp, grain_image


def mix(tree, first, second, factor, operation='MIX'):
    n = node(tree, 'ShaderNodeMixRGB')
    n.blend_type = operation
    for value, socket in [(factor, n.inputs[0]), (first, n.inputs[1]), (second, n.inputs[2])]:
        if hasattr(value, 'node'):
            tree.links.new(value, socket)
        else:
            socket.default_value = value
    return n.outputs[0]


def depth_bounds(scene):
    """Fixed artistic depth window from the current camera and visible meshes."""
    camera = scene.camera.matrix_world.inverted()
    depths = []
    for obj in scene.objects:
        if obj.hide_render or obj.type not in {'MESH', 'CURVE', 'VOLUME'}:
            continue
        for point in obj.bound_box:
            z = -(camera @ obj.matrix_world @ Vector(point)).z
            if z > scene.camera.data.clip_start:
                depths.append(z)
    if not depths or max(depths) - min(depths) < 1e-5:
        raise ValueError('Depth fog needs visible geometry with a nonzero camera depth range')
    return [min(depths), max(depths)]


def isolate_labels(scene, background, minimum_pixels=26):
    """Keep typography out of artistic filters in a separately rendered layer."""
    fonts = [o for o in scene.objects if o.type == 'FONT']
    if not fonts: return None
    name = scene.get('compositor_label_layer')
    if name and name in scene.view_layers:
        label_layer = scene.view_layers[name]
        group = bpy.data.collections[scene['compositor_label_collection']]
        beauty = scene['compositor_beauty_layer']
    else:
        beauty = bpy.context.view_layer.name
        group = bpy.data.collections.new('Compositor typography')
        scene.collection.children.link(group)
        for obj in fonts:
            for previous in list(obj.users_collection): previous.objects.unlink(obj)
            group.objects.link(obj)
        group.objects.link(scene.camera)
        label_layer = scene.view_layers.new('Typography')
        for child in label_layer.layer_collection.children:
            child.exclude = child.collection != group
        scene.view_layers[beauty].layer_collection.children[group.name].exclude = True
        scene['compositor_label_layer'] = label_layer.name
        scene['compositor_label_collection'] = group.name
        scene['compositor_beauty_layer'] = beauty
    color = (.001, .001, .001, 1) if sum(linear(background))/3 > .4 else (.8, .85, .87, 1)
    mat = bpy.data.materials.new('Protected typography contrast')
    mat.use_nodes = True; mat.node_tree.nodes.clear()
    emission = mat.node_tree.nodes.new('ShaderNodeEmission'); emission.inputs['Color'].default_value = color
    output = mat.node_tree.nodes.new('ShaderNodeOutputMaterial')
    mat.node_tree.links.new(emission.outputs[0], output.inputs['Surface'])
    for obj in fonts:
        obj.data.materials.clear(); obj.data.materials.append(mat)
        camera = scene.camera
        depth = -(camera.matrix_world.inverted() @ obj.matrix_world.translation).z
        span = camera.data.ortho_scale if camera.data.type == 'ORTHO' else depth*camera.data.sensor_width/camera.data.lens
        pixels_per_unit = scene.render.resolution_x/span
        scale = max(obj.matrix_world.to_scale().x, 1e-6)
        minimum = minimum_pixels/(pixels_per_unit*scale)
        estimated_em = max(1, sum(1 if ord(c)>0x300 else .55 for c in obj.data.body))
        maximum = scene.render.resolution_x*.82/(estimated_em*pixels_per_unit*scale)
        obj.data.size = min(max(obj.data.size, minimum), maximum)
    scene.render.film_transparent = True
    return {'layer': label_layer.name, 'beauty': beauty, 'objects': len(fonts), 'color': color,
            'requested_minimum_font_pixels': minimum_pixels, 'long_labels_limited_to_frame_width': True}


def build_stack(scene, recipe, name=None):
    if not hasattr(scene, 'compositing_node_group'):
        raise RuntimeError('This stack uses Blender 5.x compositor sockets; check the installed API before adapting')
    filters = recipe.get('filters', [])
    kinds = {f['type'] for f in filters}
    supported = {'exposure', 'color_balance', 'saturation', 'contrast', 'bloom', 'outline',
                 'depth_fog', 'ambient_occlusion', 'vignette', 'grain', 'duotone', 'subject_focus'}
    if kinds - supported:
        raise ValueError('Unknown filters: ' + ', '.join(sorted(kinds - supported)))
    if 'subject_focus' in kinds and scene.render.engine != 'CYCLES':
        raise ValueError('Object-index subject isolation is supported here only with Cycles; choose a compatible recipe')
    labels = isolate_labels(scene, recipe.get('background', '#071f2b'), recipe.get('label_min_pixels',26)) if recipe.get('protect_labels', True) else None
    for layer in scene.view_layers:
        if 'outline' in kinds: layer.use_pass_normal = True
        if 'depth_fog' in kinds: layer.use_pass_z = True
        if 'ambient_occlusion' in kinds: layer.use_pass_ambient_occlusion = True
        if 'subject_focus' in kinds: layer.use_pass_object_index = True
    t = bpy.data.node_groups.new('COMP_STACK_' + (name or recipe['id']), 'CompositorNodeTree')
    t.use_fake_user = True
    t.interface.new_socket(name='Image', in_out='OUTPUT', socket_type='NodeSocketColor')
    layers = node(t, 'CompositorNodeRLayers')
    if labels: layers.layer = labels['beauty']
    over = node(t, 'CompositorNodeAlphaOver', 'Straight background plus premultiplied rendered foreground')
    over.inputs['Background'].default_value = (*linear(recipe.get('background', '#071f2b')), 1)
    t.links.new(layers.outputs['Image'], over.inputs['Foreground'])
    image = over.outputs['Image']
    actual = []
    for settings in filters:
        f = dict(settings)
        kind = f['type']
        if kind == 'exposure':
            n = node(t, 'CompositorNodeExposure', 'Artistic exposure in stops')
            n.inputs['Exposure'].default_value = f['stops']
            t.links.new(image, n.inputs['Image']); image = n.outputs[0]
        elif kind == 'color_balance':
            n = node(t, 'CompositorNodeColorBalance', 'Lift / gamma / gain; preserve identity colors when required')
            n.inputs['Type'].default_value = 'Lift/Gamma/Gain'
            for key in ['Lift', 'Gamma', 'Gain']:
                socket = next(s for s in n.inputs if s.name == key and s.bl_idname == 'NodeSocketColor')
                socket.default_value = (*f.get(key.lower(), [1, 1, 1]), 1)
            t.links.new(image, n.inputs['Image']); image = n.outputs[0]
        elif kind == 'saturation':
            n = node(t, 'CompositorNodeHueSat')
            n.inputs['Saturation'].default_value = f['amount']
            t.links.new(image, n.inputs['Image']); image = n.outputs[0]
        elif kind == 'contrast':
            n = node(t, 'CompositorNodeBrightContrast')
            n.inputs['Contrast'].default_value = f['amount']
            t.links.new(image, n.inputs['Image']); image = n.outputs[0]
        elif kind == 'bloom':
            n = node(t, 'CompositorNodeGlare', 'Small highlights only; no inferred fluorescent intensity')
            n.inputs['Type'].default_value = 'Fog Glow'
            n.inputs['Threshold'].default_value = f.get('threshold', 1.4)
            n.inputs['Strength'].default_value = f.get('strength', .08)
            n.inputs['Size'].default_value = f.get('size', .25)
            t.links.new(image, n.inputs['Image']); image = n.outputs['Image']
        elif kind == 'outline':
            normal = node(t, 'CompositorNodeFilter', 'Normals show structural creases without cutting the mesh')
            normal.inputs['Type'].default_value = 'Sobel'
            t.links.new(layers.outputs['Normal'], normal.inputs['Image'])
            gray = node(t, 'CompositorNodeRGBToBW'); t.links.new(normal.outputs[0], gray.inputs[0])
            edge = math_node(t, 'MULTIPLY', gray.outputs[0], f.get('strength', .18))
            silhouette = node(t, 'CompositorNodeFilter'); silhouette.inputs['Type'].default_value = 'Sobel'
            t.links.new(layers.outputs['Alpha'], silhouette.inputs['Image'])
            edge = math_node(t, 'ADD', edge, math_node(t, 'MULTIPLY', silhouette.outputs[0], f.get('silhouette', .3)))
            tone = math_node(t, 'MAXIMUM', math_node(t, 'SUBTRACT', 1, edge), .1)
            image = mix(t, image, tone, 1, 'MULTIPLY')
        elif kind == 'depth_fog':
            near, far = f.get('range', depth_bounds(scene))
            if far <= near: raise ValueError('depth_fog range must increase')
            f['range'] = [near, far]
            mapped = node(t, 'ShaderNodeMapRange', 'Fixed scene-unit range; not optical attenuation')
            mapped.inputs['From Min'].default_value = near
            mapped.inputs['From Max'].default_value = far
            mapped.clamp = True
            t.links.new(layers.outputs['Depth'], mapped.inputs['Value'])
            factor = math_node(t, 'MULTIPLY', mapped.outputs[0], layers.outputs['Alpha'])
            factor = math_node(t, 'MULTIPLY', factor, f.get('strength', .2))
            image = mix(t, image, (*linear(f.get('color', '#597c8b')), 1), factor)
        elif kind == 'ambient_occlusion':
            ao = mix(t, (1, 1, 1, 1), layers.outputs['Ambient Occlusion'], f.get('strength', .2))
            image = mix(t, image, ao, 1, 'MULTIPLY')
        elif kind == 'vignette':
            ellipse = node(t, 'CompositorNodeEllipseMask', 'Optional frame-edge grading; inspect peripheral events')
            ellipse.inputs['Position'].default_value = (.5, .5)
            ellipse.inputs['Size'].default_value = tuple(f.get('size', [.9, .85]))
            blur = node(t, 'CompositorNodeBlur')
            size = min(scene.render.resolution_x, scene.render.resolution_y) * .2
            blur.inputs['Size'].default_value = (size, size)
            t.links.new(ellipse.outputs[0], blur.inputs['Image'])
            tone = math_node(t, 'ADD', math_node(t, 'MULTIPLY', blur.outputs[0], f.get('strength', .12)), 1-f.get('strength', .12))
            image = mix(t, image, tone, 1, 'MULTIPLY')
        elif kind == 'grain':
            n = node(t, 'CompositorNodeImage', 'Packed fixed grain; no frame-to-frame random flicker')
            n.image = grain_image(scene.render.resolution_x, scene.render.resolution_y, f.get('amount', .02))
            image = mix(t, image, n.outputs['Image'], 1, 'MULTIPLY')
        elif kind == 'duotone':
            gray = node(t, 'CompositorNodeRGBToBW'); t.links.new(image, gray.inputs[0])
            colored = ramp(t, gray.outputs[0], [(0, f['shadows']), (1, f['highlights'])])
            image = mix(t, image, colored, f.get('strength', 1))
        elif kind == 'subject_focus':
            index = int(f.get('object_index', 1))
            if not any(o.pass_index == index and not o.hide_render for o in scene.objects):
                raise ValueError(f'No visible objects assigned pass_index={index}')
            mask = node(t, 'CompositorNodeIDMask', 'Whole-object selection; no geometry slicing')
            mask.inputs['Index'].default_value = index
            mask.inputs['Anti-Alias'].default_value = True
            t.links.new(layers.outputs['Object Index'], mask.inputs['ID value'])
            grade = node(t, 'CompositorNodeHueSat')
            grade.inputs['Saturation'].default_value = f.get('background_saturation', .15)
            grade.inputs['Value'].default_value = f.get('background_value', .65)
            t.links.new(image, grade.inputs['Image'])
            image = mix(t, grade.outputs[0], image, mask.outputs[0])
        actual.append(f)
    if labels:
        text_layer = node(t, 'CompositorNodeRLayers', 'Typography bypasses filters; contrast follows the chosen background')
        text_layer.layer = labels['layer']
        over = node(t, 'CompositorNodeAlphaOver')
        t.links.new(image, over.inputs['Background'])
        t.links.new(text_layer.outputs['Image'], over.inputs['Foreground'])
        image = over.outputs[0]
    out = node(t, 'NodeGroupOutput'); t.links.new(image, out.inputs['Image'])
    scene.compositing_node_group = t
    scene['selected_compositor'] = name or recipe['id']
    return {'group': t.name, 'filters': actual, 'engine': scene.render.engine, 'protected_labels': labels,
            'requires_color_legend_review': bool(kinds & {'color_balance','saturation','duotone','subject_focus'}),
            'limits': 'Appearance operations only. Inspect real final pixels, encoded motion and scientific color labels.'}
