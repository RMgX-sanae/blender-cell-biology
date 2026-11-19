"""Build an editable marked overview and matching detail from a saved, single-camera Blender scene."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile

import bpy
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Vector

sys.path.insert(0, str(Path(__file__).resolve().parent))


def frame_bounds(data, scene):
    corners = data.view_frame(scene=scene)
    if data.type == 'PERSP':
        corners = [corner / -corner.z for corner in corners]
    return (min(p.x for p in corners), max(p.x for p in corners),
            min(p.y for p in corners), max(p.y for p in corners))


def frame_point(bounds, u, v):
    left, right, bottom, top = bounds
    return Vector((left + u*(right-left), bottom + v*(top-bottom)))


def render_objects(layer):
    if layer.exclude or layer.collection.hide_render:
        return
    yield from layer.collection.objects
    for child in layer.children:
        yield from render_objects(child)


def pick_surface(scene, camera, graph, bounds, eligible, pick):
    # Input coordinates are measured from the upper left of the overview image.
    u, top = pick['at']
    if not all(type(n) in (int, float) and 0 <= n <= 1 for n in (u, top)):
        raise ValueError('Landmark coordinates must lie inside the overview')
    name, group = pick.get('object'), pick.get('collection')
    if group:
        collection = bpy.data.collections.get(group)
        if collection is None or (scene.collection != collection and collection not in scene.collection.children_recursive[:]):
            raise ValueError('Collection is absent from the source scene: ' + group)
        objects = [o for o in collection.all_objects if o.type == 'MESH']
    elif name and (obj := scene.objects.get(name)) is not None and obj.type == 'MESH':
        objects = [obj]
    else:
        raise ValueError('Choose a concrete mesh object or a collection of mesh objects')

    xy = frame_point(bounds, u, 1-top)
    basis = camera.matrix_world.normalized()
    if camera.data.type == 'ORTHO':
        origin = basis @ Vector((*xy, 0))
        direction = basis.to_quaternion() @ Vector((0, 0, -1))
    else:
        origin = basis.translation
        direction = basis.to_quaternion() @ Vector((*xy, -1)).normalized()
    hits = []
    for original in objects:
        if original not in eligible or original.hide_render or not original.visible_get(view_layer=bpy.context.view_layer):
            continue
        mesh = original.evaluated_get(graph)
        inverse = mesh.matrix_world.inverted()
        hit, local, normal, face = mesh.ray_cast(inverse @ origin, inverse.to_3x3() @ direction)
        if hit:
            world = mesh.matrix_world @ local
            screen = world_to_camera_view(scene, camera, world)
            if camera.data.clip_start <= screen.z <= camera.data.clip_end:
                normal_world = (inverse.to_3x3().transposed() @ normal).normalized()
                hits.append((screen.z, original, mesh.matrix_world.copy(), face, local, world, normal_world))
    if not hits:
        raise ValueError('No selected source mesh at ' + str(pick['at']))
    _, obj, matrix, face, local, world, normal = min(hits, key=lambda item: item[0])
    if name and obj.name != name:
        raise ValueError(f'Expected {name}, but the nearer selected surface is {obj.name}')
    if any(m.show_viewport != m.show_render or
           (m.type == 'SUBSURF' and m.show_render and m.levels != m.render_levels) for m in obj.modifiers):
        raise ValueError('Match viewport and render geometry before picking the landmark on ' + obj.name)
    return {'object':obj.name, 'mesh':obj.data.name, 'ray_face_index':face,
            'overview_at':list(pick['at']), 'local':list(local), 'world':list(world),
            'world_normal':list(normal), 'materials':[m.name for m in obj.data.materials if m],
            'matrix_world':[list(row) for row in matrix],
            'source':{k:obj[k] for k in ('source_role','source_sha256','source_chain','instance_id') if k in obj}}


def detail_camera(scene, original, bounds, focus, zoom):
    data = original.data.copy()
    data.animation_data_clear()
    if data.type == 'ORTHO':
        data.ortho_scale /= zoom
    else:
        data.lens *= zoom
    # Native camera shifts work in sensor units, including vertical sensor fit.
    data.shift_x = data.shift_y = 0
    origin = frame_point(frame_bounds(data, scene), .5, .5)
    data.shift_x = data.shift_y = 1
    unit_shift = frame_point(frame_bounds(data, scene), .5, .5) - origin
    u, top = focus['overview_at']
    target = frame_point(bounds, u, 1-top) - origin
    data.shift_x = target.x / unit_shift.x
    data.shift_y = target.y / unit_shift.y
    camera = bpy.data.objects.new('Detail ' + focus['object'], data)
    scene.collection.objects.link(camera)
    camera.matrix_world = original.matrix_world.normalized()
    bpy.context.view_layer.update()
    return camera


def add_locator(scene, at, radius_fraction):
    original = scene.compositing_node_group if scene.render.use_compositing else None
    if original and original.get('linked_view_overlay'):
        original = original.get('linked_view_base')
    if original:
        tree = original.copy()
        tree['linked_view_base'] = original
        output = next(n for n in tree.nodes if n.type == 'GROUP_OUTPUT' and n.is_active_output)
        image = output.inputs['Image'].links[0].from_socket
    else:
        tree = bpy.data.node_groups.new('Overview compositor', 'CompositorNodeTree')
        tree.interface.new_socket(name='Image', in_out='OUTPUT', socket_type='NodeSocketColor')
        output = tree.nodes.new('NodeGroupOutput')
        image = tree.nodes.new('CompositorNodeRLayers').outputs['Image']
    tree.name = 'Overview locator ' + tree.name
    tree['linked_view_overlay'] = True
    scene.compositing_node_group = tree
    scene.render.use_compositing = True
    width, height = scene.render.resolution_x, scene.render.resolution_y
    radius = radius_fraction * min(width, height) / width
    stroke = .0024 * min(width, height) / width
    parent = tree.nodes.new('NodeFrame')
    parent.label = 'Overview landmark'
    parent.location = (output.location.x, output.location.y - 250)

    def stroke_mask(label, outer, inner, row):
        previous = None
        gap = 1.15 * radius
        shapes = [('Add',2*outer,2*outer), ('Subtract',2*inner,2*inner),
                  ('Subtract',gap,2.2*outer), ('Subtract',2.2*outer,gap)]
        for i, (operation, w, h) in enumerate(shapes):
            box = tree.nodes.new('CompositorNodeBoxMask')
            box.parent, box.location = parent, (210*i, -270*row)
            box.label = label
            box.inputs['Position'].default_value = (at[0], 1-at[1])
            box.inputs['Size'].default_value = (w,h)
            box.inputs['Operation'].default_value = operation
            if previous is not None:
                tree.links.new(previous, box.inputs['Mask'])
            previous = box.outputs['Mask']
        return previous

    for i, (label, outer, inner, color) in enumerate([
        ('Dark edge',radius+1.3*stroke,radius-1.3*stroke,(.006,.009,.012,1)),
        ('Light stroke',radius+.48*stroke,radius-.48*stroke,(1.5,1.48,1.36,1)),
    ]):
        mask = stroke_mask(label, outer, inner, i)
        over = tree.nodes.new('CompositorNodeAlphaOver')
        over.parent, over.location = parent, (865+i*200,-270*i)
        over.inputs['Foreground'].default_value = color
        tree.links.new(mask, over.inputs['Factor'])
        tree.links.new(image, over.inputs['Background'])
        image = over.outputs['Image']
    tree.links.new(image, output.inputs['Image'])
    output.location.x += 1330
    scene['linked_view_locator'] = json.dumps({'at':at,'radius':radius_fraction})


def save_view(scene, camera, role, directory):
    scene.camera = camera
    scene['linked_view_role'] = role
    scene.render.filepath = '//' + role
    target = directory / (role + '.blend')
    result = bpy.ops.wm.save_as_mainfile(filepath=str(target), copy=True, relative_remap=True)
    if 'FINISHED' not in result:
        raise RuntimeError('Could not save ' + role)


def build(brief, output):
    scene = bpy.context.scene
    source = Path(bpy.data.filepath)
    if not source.is_file() or output.exists():
        raise ValueError('Open a saved source .blend and select a fresh output directory')
    recipe_bytes = brief.read_bytes()
    recipe = json.loads(recipe_bytes)['view_handoff']
    if not all(isinstance(recipe.get(key), str) and recipe[key].strip() for key in ('entity','scale_basis')):
        raise ValueError('Name the focused entity and the basis of the scale in the brief')
    if scene.get('linked_view_locator'):
        previous = scene.compositing_node_group
        if previous and previous.get('linked_view_overlay'):
            base = previous.get('linked_view_base')
            scene.compositing_node_group = base
            scene.render.use_compositing = base is not None
        del scene['linked_view_locator']
    frame = recipe.get('frame','saved')
    if frame != 'saved':
        if type(frame) is not int or not scene.frame_start <= frame <= scene.frame_end:
            raise ValueError('Frame is outside the saved timeline')
        scene.frame_set(frame)
    graph = bpy.context.evaluated_depsgraph_get()
    overview = scene.camera
    if overview is None or overview.data.type not in {'ORTHO','PERSP'} or scene.render.use_border or scene.render.use_multiview:
        raise ValueError('Use a full-frame, single orthographic or perspective camera')
    camera = overview.evaluated_get(graph)
    bounds = frame_bounds(camera.data, scene)
    zoom, radius = recipe['zoom'], recipe.get('locator_radius',.035)
    if (type(zoom) not in (float,int) or not math.isfinite(zoom) or zoom <= 1 or
            type(radius) not in (float,int) or not .008 <= radius <= .45):
        raise ValueError('Use a finite zoom greater than one and a locator radius in [.008, .45]')
    eligible = set(render_objects(bpy.context.view_layer.layer_collection))
    focus = pick_surface(scene,camera,graph,bounds,eligible,recipe['focus'])
    neighbors = [pick_surface(scene,camera,graph,bounds,eligible,p) for p in recipe.get('neighbors',[])]
    if any(p['overview_at'] == focus['overview_at'] for p in neighbors):
        raise ValueError('Choose a neighboring landmark distinct from the focus')
    u, top = focus['overview_at']
    half = .5 / zoom
    crop = [u-half, top-half, u+half, top+half]
    sx, sy = scene.render.resolution_x, scene.render.resolution_y
    r = (radius+.0032)*min(sx,sy)
    if min(crop) < 0 or max(crop) > 1 or min(u*sx,(1-u)*sx,top*sy,(1-top)*sy) <= r:
        raise ValueError('The detail field or locator extends beyond the overview')
    detail = detail_camera(scene,camera,bounds,focus,zoom)
    actual = frame_bounds(detail.data,scene)
    lo = frame_point(bounds,crop[0],1-crop[3])
    hi = frame_point(bounds,crop[2],1-crop[1])
    expected = (lo.x,hi.x,lo.y,hi.y)
    tolerance = max(bounds[1]-bounds[0],bounds[3]-bounds[2])*2e-5
    if max(abs(a-b) for a,b in zip(actual,expected)) > tolerance:
        raise ValueError('The native camera cannot represent the requested field of view')
    for point in [focus,*neighbors]:
        projected = world_to_camera_view(scene,detail,Vector(point['world']))
        uv = [projected.x,1-projected.y]
        reference = [.5+zoom*(a-b) for a,b in zip(point['overview_at'],focus['overview_at'])]
        if min(uv)<0 or max(uv)>1 or max(abs(a-b) for a,b in zip(uv,reference))>1e-4:
            raise ValueError('Landmark falls outside the detail or fails its projection: '+point['object'])
        point['detail_at'] = uv
    with source.open('rb') as stream:
        source_hash = hashlib.file_digest(stream,'sha256').hexdigest()
    record = {'entity':recipe['entity'],'scale_basis':recipe['scale_basis'],'source_file':str(source),
              'brief_sha256':hashlib.sha256(recipe_bytes).hexdigest(),
              'source_sha256':source_hash,'scene':scene.name,'frame':scene.frame_current,'subframe':scene.frame_subframe,
              'cameras':[overview.name,detail.name],'camera_type':camera.data.type,'zoom':zoom,
              'camera_to_world':[list(row) for row in detail.matrix_world],
              'detail_rect_in_overview':crop,'locator_radius':radius,'focus':focus,'neighbors':neighbors}
    # Each copy is an editable still at the same instant. Keep the original file untouched.
    for marker in scene.timeline_markers:
        if marker.camera:
            marker.camera = None
    scene.frame_start = scene.frame_end = scene.frame_current
    scene['linked_view_source'] = source_hash
    scene['linked_view_subframe'] = scene.frame_subframe
    output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as temporary:
        staging = Path(temporary)
        original_tree = scene.compositing_node_group
        compositing = scene.render.use_compositing
        add_locator(scene, focus['overview_at'], radius)
        save_view(scene, overview, 'overview', staging)
        scene.compositing_node_group = original_tree
        scene.render.use_compositing = compositing
        del scene['linked_view_locator']
        if recipe.get('detail_spacefill') is not None:
            from linked_spacefill import switch_to_spacefill
            record['detail_representation'] = switch_to_spacefill(
                scene, detail, focus, recipe['detail_spacefill'], brief.parent)
        save_view(scene, detail, 'detail', staging)
        (staging/'handoff.json').write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n')
        if output.exists():
            raise FileExistsError(output)
        staging.rename(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--brief',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args = parser.parse_args(sys.argv[sys.argv.index('--')+1:])
    build(args.brief.resolve(),args.output.resolve())
