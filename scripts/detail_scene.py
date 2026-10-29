"""Shared editable geometry and image-reference style implementation for Blender.

Run through build_molecular_scene.py. Colors in profiles are sRGB hex, converted
to linear shader values; material settings are artistic, not measured biology.
"""
import json
import hashlib
import math
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Vector


def linear(hex_color):
    values = [int(hex_color.lstrip('#')[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    return tuple(v / 12.92 if v <= .04045 else ((v + .055) / 1.055) ** 2.4 for v in values)


def collection(name):
    c = bpy.data.collections.new(name)
    bpy.context.scene.collection.children.link(c)
    return c


def mesh(name, vertices, faces, materials, group, indices=None):
    data = bpy.data.meshes.new(name)
    data.from_pydata(vertices, [], faces)
    data.update()
    obj = bpy.data.objects.new(name, data)
    group.objects.link(obj)
    for material in materials:
        data.materials.append(material)
    for polygon in data.polygons:
        polygon.use_smooth = True
    if indices is not None:
        data.polygons.foreach_set('material_index', indices)
    return obj


def material(name, color, roughness=.48, bump=.012):
    m = bpy.data.materials.new(name)
    m.diffuse_color = (*linear(color), 1)
    m.use_nodes = True
    nodes, links = m.node_tree.nodes, m.node_tree.links
    shader = nodes.get('Principled BSDF')
    shader.inputs['Base Color'].default_value = m.diffuse_color
    shader.inputs['Roughness'].default_value = roughness
    shader.inputs['Specular IOR Level'].default_value = .32
    shader.inputs['Subsurface Weight'].default_value = .025
    shader.inputs['Subsurface Radius'].default_value = (.08, .05, .025)
    if bump:
        noise = nodes.new('ShaderNodeTexNoise')
        noise.label = 'Subtle surface texture; not atom coordinates'
        noise.inputs['Scale'].default_value = 16
        noise.inputs['Detail'].default_value = 2
        texture = nodes.new('ShaderNodeTexCoord')
        links.new(texture.outputs['Object'], noise.inputs['Vector'])
        node = nodes.new('ShaderNodeBump')
        node.inputs['Strength'].default_value = .22
        node.inputs['Distance'].default_value = bump
        links.new(noise.outputs['Fac'], node.inputs['Height'])
        links.new(node.outputs['Normal'], shader.inputs['Normal'])
    return m


class Parts:
    """Combine small components; avoids thousands of object datablocks."""
    def __init__(self):
        self.v, self.f, self.mi = [], [], []

    def ball(self, position, radius, mat=0, seg=12, rings=8):
        p = np.asarray(position)
        scale = np.asarray(radius if hasattr(radius, '__len__') else [radius] * 3)
        start = len(self.v)
        self.v.append(tuple(p + (0, 0, scale[2])))
        for j in range(1, rings):
            a = math.pi * j / rings
            for i in range(seg):
                t = math.tau * i / seg
                self.v.append(tuple(p + scale * (math.sin(a) * math.cos(t),
                                                math.sin(a) * math.sin(t), math.cos(a))))
        bottom = len(self.v)
        self.v.append(tuple(p - (0, 0, scale[2])))
        faces = []
        for i in range(seg):
            faces.append((start, start + 1 + i, start + 1 + (i + 1) % seg))
        for j in range(rings - 2):
            for i in range(seg):
                a = start + 1 + j * seg + i
                b = start + 1 + j * seg + (i + 1) % seg
                faces.append((a, a + seg, b + seg, b))
        for i in range(seg):
            faces.append((bottom, bottom - seg + (i + 1) % seg, bottom - seg + i))
        self.f.extend(faces)
        self.mi.extend([mat] * len(faces))

    def tube(self, points, radius, mat=0, sides=7):
        pts = [Vector(p) for p in points]
        start = len(self.v)
        # Parallel-transporting the basis prevents twist around sharp curves.
        last_u = Vector((0, 0, 1))
        for j, p in enumerate(pts):
            axis = (pts[min(j + 1, len(pts) - 1)] - pts[max(j - 1, 0)]).normalized()
            u = last_u - axis * last_u.dot(axis)
            if u.length < .01:
                u = axis.cross(Vector((0, 1, 0)))
            u.normalize()
            v = axis.cross(u)
            last_u = u
            for i in range(sides):
                self.v.append(tuple(p + radius * (u * math.cos(math.tau * i / sides)
                                                   + v * math.sin(math.tau * i / sides))))
        for j in range(len(pts) - 1):
            for i in range(sides):
                a = start + j * sides + i
                b = start + j * sides + (i + 1) % sides
                self.f.append((a, b, b + sides, a + sides))
                self.mi.append(mat)
        self.f.extend([tuple(start + i for i in reversed(range(sides))),
                       tuple(start + (len(pts) - 1) * sides + i for i in range(sides))])
        self.mi.extend([mat, mat])

    def object(self, name, materials, group):
        return mesh(name, self.v, self.f, materials, group, self.mi)


def path(name, points, radius, mat, group, cyclic=False):
    data = bpy.data.curves.new(name, 'CURVE')
    data.dimensions = '3D'
    data.bevel_depth = radius
    data.bevel_resolution = 3
    spline = data.splines.new('POLY')
    spline.points.add(len(points) - 1)
    for p, co in zip(spline.points, points):
        p.co = (*co, 1)
    spline.use_cyclic_u = cyclic
    data.materials.append(mat)
    obj = bpy.data.objects.new(name, data)
    group.objects.link(obj)
    return obj


def load_chain(directory, chain, mat, group, name=None, representation='surface'):
    data = np.load(Path(directory) / f'chain_{chain}.npz', allow_pickle=False)
    if representation == 'spacefill':
        parts = Parts()
        radii = {'C': .17, 'N': .155, 'O': .152, 'S': .18, 'P': .18}
        for pos, element in zip(data['atoms'], data['elements']):
            parts.ball(pos, radii.get(element, .17), seg=10, rings=7)
        obj = parts.object(name or f'Experimental chain {chain}', [mat], group)
        obj['evidence'] = 'experimental heavy-atom coordinates; conventional VDW radii; H omitted'
    else:
        obj = mesh(name or f'Experimental chain {chain}', data['vertices'].tolist(),
                   data['faces'].tolist(), [mat], group)
        smooth = obj.modifiers.new('Envelope smoothing below structural detail', 'SMOOTH')
        smooth.factor, smooth.iterations = .35, 2
        obj['evidence'] = 'experimental-coordinate-derived Gaussian envelope'
    obj['source_chain'] = chain
    return obj, data


def duplicate(obj, name, group):
    copy = obj.copy()
    copy.data = obj.data
    copy.name = name
    group.objects.link(copy)
    return copy


def aim(obj, target):
    obj.rotation_euler = (Vector(target) - obj.location).to_track_quat('-Z', 'Y').to_euler()


def fit_orthographic(camera, objects, margin=.08):
    """Fit selected subject bounds, preserving viewing direction and clear margins."""
    if camera.data.type != 'ORTHO':
        raise ValueError('Subject fitting here is for orthographic cameras only.')
    bpy.context.view_layer.update()
    inverse = camera.matrix_world.inverted()
    points = [inverse @ (obj.matrix_world @ Vector(p))
              for obj in objects if obj.type in {'MESH', 'CURVE'} for p in obj.bound_box]
    lo = np.min(points, axis=0)
    hi = np.max(points, axis=0)
    center = (hi+lo)/2
    camera.location += camera.rotation_euler.to_matrix() @ Vector((center[0],center[1],0))
    scene = bpy.context.scene
    aspect = scene.render.resolution_x / scene.render.resolution_y
    camera.data.ortho_scale = max(hi[0]-lo[0],(hi[1]-lo[1])*aspect)/(1-2*margin)


def compositors(scene, profile, selected=None):
    """Save independent named branches, including a neutral diagnostic branch."""
    selected = selected or profile['compositor']['preset']
    if not hasattr(scene, 'compositing_node_group'):
        raise RuntimeError('Detailed compositor builder requires Blender 5.x; use the documented legacy adaptation for 4.x.')
    groups = {}
    for mode in ('clean', 'reference', 'monochrome'):
        tree = bpy.data.node_groups.new('COMP_' + mode, 'CompositorNodeTree')
        tree.use_fake_user = True  # Preserve unselected alternatives after reopening.
        tree.interface.new_socket(name='Image', in_out='OUTPUT', socket_type='NodeSocketColor')
        render = tree.nodes.new('CompositorNodeRLayers')
        render.location = (-600, 80)
        output = tree.nodes.new('NodeGroupOutput')
        output.location = (400, 80)
        over = tree.nodes.new('CompositorNodeAlphaOver')
        over.location = (-360, 80)
        over.inputs['Background'].default_value = (*linear(profile['background']), 1)
        over.inputs['Factor'].default_value = 1
        tree.links.new(render.outputs['Image'], over.inputs['Foreground'])
        socket = over.outputs['Image']
        if mode == 'reference':
            grade = tree.nodes.new('CompositorNodeBrightContrast')
            grade.location = (80, 80)
            grade.inputs['Contrast'].default_value = profile['compositor'].get('contrast', 2)
            tree.links.new(socket, grade.inputs['Image'])
            socket = grade.outputs['Image']
        if mode == 'monochrome':
            gray = tree.nodes.new('CompositorNodeRGBToBW')
            gray.location = (80, 80)
            tree.links.new(socket, gray.inputs['Image'])
            socket = gray.outputs[0]
        # Grade after compositing: a contrast offset must not create negative
        # premultiplied RGB around an otherwise transparent render background.
        tree.links.new(socket, output.inputs['Image'])
        groups[mode] = tree
    if selected not in groups:
        raise ValueError(f'Unknown compositor {selected}')
    scene.compositing_node_group = groups[selected]
    scene['selected_compositor'] = selected


def studio(profile, args, camera_defaults):
    scene = bpy.context.scene
    scene.render.engine = 'CYCLES' if args.engine == 'cycles' else 'BLENDER_EEVEE'
    if scene.render.engine == 'BLENDER_EEVEE': scene.eevee.taa_render_samples = args.samples
    scene.render.resolution_x = args.width
    scene.render.resolution_y = round(args.width * 5 / 8)
    scene.render.resolution_percentage = 100
    scene.render.image_settings.file_format = 'PNG'
    scene.render.image_settings.color_mode = 'RGBA'
    scene.render.film_transparent = True
    scene.render.fps = args.fps
    scene.frame_start, scene.frame_end = 1, args.frames
    scene.view_settings.view_transform = profile.get('view_transform', 'AgX')
    scene.view_settings.look = 'AgX - Medium High Contrast'
    scene.view_settings.exposure = profile.get('exposure', 0)
    scene.world = bpy.data.worlds.new('Environment')
    scene.world.use_nodes = True
    scene.world.node_tree.nodes['Background'].inputs['Color'].default_value = (*linear(profile['world']), 1)
    scene.world.node_tree.nodes['Background'].inputs['Strength'].default_value = .35
    if args.engine == 'cycles':
        scene.cycles.samples = args.samples
        scene.cycles.use_denoising = True
        scene.cycles.adaptive_threshold = .025
        scene.cycles.max_bounces = 8
        prefs = bpy.context.preferences.addons['cycles'].preferences
        try:
            prefs.compute_device_type = 'METAL'
            prefs.get_devices()
            metal = [d for d in prefs.devices if d.type == 'METAL']
            for device in prefs.devices:
                device.use = device in metal
            if metal:
                scene.cycles.device = 'GPU'
        except (TypeError, RuntimeError):
            pass
    group = collection('90 Camera and lighting')
    light_scale = camera_defaults.get('light_scale', 1)
    for item in profile['lighting']:
        data = bpy.data.lights.new(item['name'], 'AREA')
        data.energy = item['power'] * light_scale ** 2
        data.color = linear(item['color'])
        data.shape = 'DISK'
        data.size = item['size'] * light_scale
        obj = bpy.data.objects.new(item['name'], data)
        group.objects.link(obj)
        obj.location = tuple(x * light_scale for x in item['position'])
        aim(obj, camera_defaults.get('target', (0, 0, 0)))
    settings = {**camera_defaults, **profile.get('camera', {})}
    data = bpy.data.cameras.new('Reference-controlled camera')
    camera = bpy.data.objects.new('Reference-controlled camera', data)
    group.objects.link(camera)
    camera.location = settings['position']
    aim(camera, settings.get('target', (0, 0, 0)))
    data.type = settings.get('type', 'ORTHO')
    data.ortho_scale = settings['scale']
    data.lens = settings.get('lens', 52)
    data.clip_end = 1000
    if settings.get('fstop'):
        data.dof.use_dof = True
        data.dof.aperture_fstop = settings['fstop']
        data.dof.focus_distance = (camera.location - Vector(settings.get('target', (0, 0, 0)))).length
    scene.camera = camera
    compositors(scene, profile, args.compositor)
    scene['style_profile'] = json.dumps(profile, ensure_ascii=False)
    return camera


def save_and_render(args, manifest):
    scene = bpy.context.scene
    scene.frame_set(args.frame)
    args.output.mkdir(parents=True, exist_ok=False)
    scene.render.filepath = str(args.output / 'preview.png')
    manifest.update({'blender': bpy.app.version_string, 'engine': scene.render.engine,
                     'samples': args.samples, 'width': args.width,
                     'height': scene.render.resolution_y, 'fps': args.fps,
                     'frames': args.frames, 'preview_frame': args.frame,
                     'compositor': scene['selected_compositor'],
                     'style': json.loads(scene['style_profile'])})
    manifest['build_parameters'] = {k:str(v) if isinstance(v,Path) else v
                                    for k,v in vars(args).items()}
    manifest['command_arguments'] = sys.argv
    manifest['camera_resolved'] = {'location':list(scene.camera.location),
                                   'rotation_euler':list(scene.camera.rotation_euler),
                                   'type':scene.camera.data.type,
                                   'ortho_scale':scene.camera.data.ortho_scale,
                                   'lens':scene.camera.data.lens}
    manifest['generator_hashes'] = {p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                                    for p in Path(__file__).resolve().parent.glob('*.py')}
    (args.output / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    bpy.ops.wm.save_as_mainfile(filepath=str(args.output / 'scene.blend'))
    if args.animate:
        (args.output / 'frames').mkdir()
        scene.render.filepath = str(args.output / 'frames' / 'frame_')
        bpy.ops.render.render(animation=True)
    elif args.render:
        bpy.ops.render.render(write_still=True)
