"""Editable whole-protein K63 net-transfer views with protected typography."""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Matrix, Euler

sys.path.insert(0, str(Path(__file__).resolve().parent))
import detail_scene as ds
import style_nodes as sn
import cellular_environment as environment
from build_phosphotransfer_scene import visibility
from compositor_stack import build_stack


def digest(path):
    import hashlib
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def keyed_points(name, coordinates, group, edges=()):
    mesh = bpy.data.meshes.new(name); mesh.from_pydata(coordinates[0].tolist(), list(edges), []); mesh.update()
    obj = bpy.data.objects.new(name, mesh); group.objects.link(obj)
    for i, points in enumerate(coordinates):
        key = obj.shape_key_add(name=f'frame_{i+1:04d}')
        key.data.foreach_set('co', points.astype(np.float32).ravel()); key.interpolation = 'KEY_LINEAR'
    keys = mesh.shape_keys; keys.use_relative = False
    for i, key in enumerate(keys.key_blocks):
        keys.eval_time = key.frame; keys.keyframe_insert('eval_time', frame=i+1)
    return obj


def envelope(obj, material):
    tree = sn.density(obj, material, {'volume_voxel_nm': .065, 'volume_radius_nm': .20}, isosurface=True)
    target = next(n for n in tree.nodes if n.bl_idname == 'NodeGroupOutput')
    source = target.inputs['Geometry'].links[0].from_socket
    smooth = sn.node(tree, 'GeometryNodeSetShadeSmooth')
    tree.links.new(source, smooth.inputs['Geometry']); tree.links.new(smooth.outputs['Geometry'], target.inputs['Geometry'])


def annotation_line(obj, material, radius):
    tree, src, target = sn.geometry_tree(obj, 'GN_net_bond_annotation')
    curve = sn.node(tree, 'GeometryNodeMeshToCurve'); tree.links.new(src, curve.inputs['Mesh'])
    circle = sn.node(tree, 'GeometryNodeCurvePrimitiveCircle'); circle.inputs['Resolution'].default_value = 12; circle.inputs['Radius'].default_value = radius
    solid = sn.node(tree, 'GeometryNodeCurveToMesh'); tree.links.new(curve.outputs['Curve'], solid.inputs['Curve']); tree.links.new(circle.outputs['Curve'], solid.inputs['Profile Curve'])
    paint = sn.node(tree, 'GeometryNodeSetMaterial'); paint.inputs['Material'].default_value = material
    tree.links.new(solid.outputs['Mesh'], paint.inputs['Geometry']); tree.links.new(paint.outputs['Geometry'], target)


def background(asset_dir, scale, components=('donor', 'e2', 'uev')):
    groups = [ds.collection('Molecular background '+n) for n in ('far', 'middle', 'near')]
    colors = dict(zip(components, ['#778f9a', '#789994', '#8f839e'])); meshes = {}; records = []
    for k, color in colors.items():
        d = np.load(asset_dir/('background_'+k+'.npz'))
        mesh = bpy.data.meshes.new('Whole background '+k); mesh.from_pydata(d['vertices'].tolist(), [], d['faces'].tolist()); mesh.update()
        mesh.materials.append(ds.material('Background '+k, color, .65, 0))
        for polygon in mesh.polygons: polygon.use_smooth = True
        meshes[k] = mesh
    rng = np.random.default_rng(924)
    # Uneven depths, orientations and spacings; source proteins supply shape.
    placements = [(-11, 3.4, 0, 'e2'), (9.5, 2.6, 0, 'uev'), (3.4, -4.4, 0, 'e2'),
                  (-7.7, -3.1, 1, 'uev'), (9.4, -2.8, 1, 'e2'), (-1.1, 3.9, 1, 'donor'),
                  (-12.3, -1.3, 2, 'uev'), (12.3, 2.2, 2, 'e2')]
    for i, (x, y, layer, k) in enumerate(placements):
        k = dict(zip(('donor', 'e2', 'uev'), components))[k]
        obj = bpy.data.objects.new(f'Background {k} {i}', meshes[k]); groups[layer].objects.link(obj)
        angles = rng.uniform(0, 6.28, 3); projection = [.55, .9, 1.4][layer]
        matrix = Matrix.Translation((x, y, [-20, -12, -6][layer]))@Euler(angles).to_matrix().to_4x4()@Matrix.Scale(scale*projection, 4)
        obj.matrix_world = matrix; obj['environment_only'] = True; obj['evidence'] = 'Artistic placement of complete source protein shape'
        records.append({'object': obj.name, 'component': k, 'matrix': np.array(matrix).tolist(), 'depth_layer': layer})
    return groups, records


def paired_composite(scene, close, annotations, bg_layers, background_color, width):
    """Compose synchronized whole-molecule panels after build_stack initializes labels."""
    tree = scene.compositing_node_group
    # The label collection is created by build_stack after our view layers.
    for layer in [close, annotations, *bg_layers]: layer.layer_collection.children[scene['compositor_label_collection']].exclude = True
    base = next(n for n in tree.nodes if n.label == 'Straight background plus premultiplied rendered foreground')
    left_render = base.inputs['Foreground'].links[0].from_socket
    right_render = sn.node(tree, 'CompositorNodeRLayers'); right_render.layer = close.name
    domain = sn.math_node(tree, 'MULTIPLY', left_render.node.outputs['Alpha'], 0)
    mask = sn.node(tree, 'CompositorNodeBoxMask'); mask.inputs['Position'].default_value = (.25, .5); mask.inputs['Size'].default_value = (.5, 1)
    tree.links.new(domain, mask.inputs['Mask'])
    mix = sn.node(tree, 'ShaderNodeMixRGB'); tree.links.new(mask.outputs['Mask'], mix.inputs[0]); tree.links.new(right_render.outputs['Image'], mix.inputs[1]); tree.links.new(left_render, mix.inputs[2])
    # MixRGB does not select both premultiplied alpha channels with the RGB.
    # Explicitly select alpha too, otherwise the left panel adds onto a light
    # background instead of occluding it and looks incorrectly translucent.
    inverse = sn.math_node(tree, 'SUBTRACT', 1, mask.outputs['Mask'])
    left_alpha = sn.math_node(tree, 'MULTIPLY', left_render.node.outputs['Alpha'], mask.outputs['Mask'])
    right_alpha = sn.math_node(tree, 'MULTIPLY', right_render.outputs['Alpha'], inverse)
    alpha = sn.math_node(tree, 'ADD', left_alpha, right_alpha)
    panels = sn.node(tree, 'CompositorNodeSetAlpha', 'Select matching premultiplied panel alpha')
    panels.inputs['Type'].default_value = 'Replace Alpha'
    tree.links.new(mix.outputs[0], panels.inputs['Image']); tree.links.new(alpha, panels.inputs['Alpha'])
    window = sn.node(tree, 'CompositorNodeBoxMask', 'Framing window, no molecular mesh section')
    # Box-mask height is relative to image width in this Blender compositor.
    window.inputs['Position'].default_value = (.5, .5); window.inputs['Size'].default_value = (1, .32)
    tree.links.new(domain, window.inputs['Mask'])
    crop = sn.node(tree, 'CompositorNodeSetAlpha', 'Apply framing to premultiplied color and alpha')
    crop.inputs['Type'].default_value = 'Apply Mask'
    tree.links.new(panels.outputs[0], crop.inputs['Image']); tree.links.new(window.outputs['Mask'], crop.inputs['Alpha'])
    annotate = sn.node(tree, 'CompositorNodeRLayers'); annotate.layer = annotations.name
    over = sn.node(tree, 'CompositorNodeAlphaOver'); tree.links.new(crop.outputs[0], over.inputs['Background']); tree.links.new(annotate.outputs['Image'], over.inputs['Foreground'])
    tree.links.new(over.outputs['Image'], base.inputs['Foreground'])
    tree.links.new(environment.composite(tree, bg_layers, background_color, width), base.inputs['Background'])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--font-zh', type=Path, required=True); p.add_argument('--font-en', type=Path, required=True)
    p.add_argument('--style', choices=['dark', 'light'], default='dark')
    p.add_argument('--recipe', choices=['neutral', 'soft-cinema', 'fluorescent-soft', 'warm-print', 'cool-duotone', 'crisp-neutral'], default='neutral')
    p.add_argument('--engine', choices=['eevee', 'cycles'], default='eevee')
    p.add_argument('--width', type=int, default=1280); p.add_argument('--samples', type=int, default=24)
    p.add_argument('--render', action='store_true')
    a = p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.output.exists(): p.error('Use a fresh scene directory')
    if not a.font_zh.is_file() or not a.font_en.is_file(): p.error('Supply Songti and Times New Roman fonts')
    source = json.loads((a.asset_dir/'manifest.json').read_text())
    if source['kind'] != 'rnf4-ubc13-uev-k63-net-transfer' or not source['passed']: p.error('Prepare and check the source transfer first')
    for name, expected in source['data_sha256'].items():
        if digest(a.asset_dir/name) != expected: p.error('Prepared data changed: '+name)
    data = dict(np.load(a.asset_dir/'motion.npz')); count = source['frames']; config = source['config']
    bpy.ops.object.select_all(action='SELECT'); bpy.ops.object.delete(use_global=False)
    scene = bpy.context.scene; overview = ds.collection('Whole reaction complex'); detail = ds.collection('Whole ubiquitin pair')
    overlay = ds.collection('Reaction annotations'); labels = ds.collection('Figure labels')
    root = Path(__file__).resolve().parents[1]
    profile = json.loads((root/'assets'/('style-molecular-reference.json' if a.style == 'dark' else 'style-studio.json')).read_text())
    profile.update(camera={}, bump=0, background='#071f2b' if a.style == 'dark' else '#edf0e9')
    a.fps, a.frames, a.compositor = source['fps'], count, 'clean'
    camera = ds.studio(profile, a, {'position': [0, 0, 40], 'target': [0, 0, 0], 'scale': 22, 'light_scale': 1.2})
    camera.matrix_world = Matrix.Translation((0, 0, 40)); camera.data.dof.use_dof = False
    scene.render.resolution_y = 2*round(a.width*9/32)
    if hasattr(scene, 'eevee'): scene.eevee.taa_render_samples = a.samples
    front = np.array(config['display']['view_direction'], float); front /= np.linalg.norm(front)
    right = np.cross([0, 0, 1], front); right /= np.linalg.norm(right); up = np.cross(front, right)
    basis = np.column_stack([right, up, front]); transforms = {}; parents = {}
    framing = ds.collection('View framing transforms')
    for view, keys, offset in [('overview', list(source['atom_ids']), [-5.5, 0, 0]), ('detail', ['donor', 'acceptor'], [5.5, 0, 0])]:
        coordinates = np.concatenate([data[k] for k in keys], axis=1)@basis
        low, high = coordinates.min(1), coordinates.max(1)
        center = (low+high)/2
        kernel = np.ones(13)/13
        center = np.column_stack([np.convolve(np.pad(center[:, axis], 6, mode='edge'), kernel, mode='valid') for axis in range(3)])
        extent = np.maximum(high-center, center-low)
        allowed = np.minimum(8.6/(2*extent[:, 0]), 6.2/(2*extent[:, 1]))
        # A conservative moving minimum followed by smoothing avoids cropping
        # or abrupt zoom when a new atom defines the projected bounding box.
        padded = np.pad(allowed, 12, mode='edge')
        scale = np.array([padded[i:i+25].min() for i in range(count+1)])
        scale = np.convolve(np.pad(scale, 6, mode='edge'), kernel, mode='valid')*.98
        assert np.all(scale <= allowed+1e-10)
        matrices = np.repeat(np.eye(4)[None], count+1, axis=0)
        matrices[:, :3, :3] = basis.T[None]*scale[:, None, None]
        matrices[:, :3, 3] = offset-center*scale[:, None]
        transforms[view] = matrices
        parent = bpy.data.objects.new(view+' smooth framing', None); framing.objects.link(parent); parents[view] = parent
        parent.rotation_mode = 'QUATERNION'; parent.rotation_quaternion = Matrix(basis.T).to_quaternion()
        for i, matrix in enumerate(matrices):
            parent.location = matrix[:3, 3]; parent.scale = [scale[i]]*3
            parent.keyframe_insert('location', frame=i+1); parent.keyframe_insert('scale', frame=i+1)
    materials = {k: ds.material(k, color, .52, 0) for k, color in config['display']['colors'].items()}
    actors = []
    for k in source['atom_ids']:
        obj = keyed_points(k+' overview', data[k], overview); obj.parent = parents['overview']
        if k == 'zinc': sn.instance_spheres(obj, materials[k], radius=.13)
        else: envelope(obj, materials[k])
        obj['component'] = k; obj['view'] = 'overview'; obj['atom_ids'] = json.dumps(source['atom_ids'][k]); actors.append(obj.name)
        if k in ('donor', 'acceptor'):
            copy = ds.duplicate(obj, k+' detail', detail); copy.parent = parents['detail']
            copy['view'] = 'detail'; actors.append(copy.name)
    # Explicit graphic overlays identify the net bond exchange. Full molecular
    # shapes remain present; no residue fragments are cut out of a protein.
    ids = source['reaction_atom_indices']; fraction = data['transfer_fraction']; marks = []
    timing = source.get('timing', {'acceptor_arrival': .39, 'product_formed': .76})
    ink = sn.material('High contrast net-bond annotation', {'material': 'paint'}, '#fff1b1' if a.style == 'dark' else '#3c2454', 'paint')
    for name, first, second, view, visible in [
        ('E2 thioester annotation', data['e2'][:, ids['sg']], data['donor'][:, ids['c']], 'overview', (np.arange(count+1)/count >= timing['acceptor_arrival'])&(fraction < .05)),
        ('K63 product bond annotation', data['donor'][:, ids['c']], data['acceptor'][:, ids['nz']], 'detail', fraction >= .98)]:
        points = np.stack([first, second], axis=1); obj = keyed_points(name, points, overlay, [(0, 1)])
        annotation_line(obj, ink, .065)
        obj.parent = parents[view]; obj['display_annotation'] = True
        for i, show in enumerate(visible):
            obj.hide_render = not bool(show); obj.keyframe_insert('hide_render', frame=i+1)
        marks.append({'object': obj.name, 'view': view, 'kind': 'net covalent bond highlight, enlarged stroke'})
    font_zh = bpy.data.fonts.load(str(a.font_zh)); font_en = bpy.data.fonts.load(str(a.font_en)); font_zh.pack(); font_en.pack()
    def text(body, x, y, size):
        curve = bpy.data.curves.new(body, 'FONT'); curve.body = body; curve.font = font_zh; curve.font_bold = font_en
        for char, form in zip(body, curve.body_format): form.use_bold = ord(char) < 128 or char in '·—'
        curve.size = size; curve.align_x = 'CENTER'
        if a.style == 'light': curve.offset = size*.01
        obj = bpy.data.objects.new(body, curve); labels.objects.link(obj); obj.location = (x, y, 12)
        return obj
    text('K63 双泛素连接', 0, 5.1, .85)
    text('组装、转移与释放' if source.get('release_included') else '复合物组装与转移', -5.5, 4.1, .64); text('供体与受体泛素', 5.5, 4.1, .64)
    for x, k, title in [(-8.5, 'ring', 'RNF4'), (-4.3, 'e2', 'Ubc13 / E2'), (0, 'uev', 'Ube2V2'), (4.4, 'donor', '供体泛素'), (8.7, 'acceptor', '受体泛素')]:
        text(title, x, -4.3, .49)
        dot = sn.point_object(title+' legend dot', [[x-1.45, -4.12, 12]], overlay)
        sn.instance_spheres(dot, sn.material(title+' legend color', {'material': 'paint'}, config['display']['colors'][k], 'paint'), radius=.09)
    for begin, end, caption in source.get('captions', config['captions']):
        obj = text(caption, 0, -5.4, .66); visibility(obj, round(begin*count)+1, round(end*count))
    product_label = text('Gly76—Lys63', 7.9, -2.2, .49)
    visibility(product_label, round(timing['product_formed']*count)+1, count)
    contact = (data['donor'][:, ids['c']]+data['acceptor'][:, ids['nz']])/2
    contact = np.einsum('nij,nj->ni', transforms['detail'][:, :3, :3], contact)+transforms['detail'][:, :3, 3]
    line = np.broadcast_to([[0., 0., 12.], [7.2, -1.6, 12.], [8.7, -1.6, 12.]], (count+1, 3, 3)).copy()
    line[:, 0, :2] = contact[:, :2]
    leader = keyed_points('Product bond leader', line, overlay, [(0, 1), (1, 2)])
    annotation_line(leader, ink, .012)
    visibility(leader, round(timing['product_formed']*count)+1, count)
    divider = sn.material('High contrast divider', {'material': 'paint'}, '#bad6d9' if a.style == 'dark' else '#365a64', 'paint')
    ds.path('Panel separator', [(0, -3.5, 12), (0, 3.55, 12)], .016, divider, overlay)
    bg_groups, bg_record = background(a.asset_dir, float(np.linalg.norm(transforms['overview'][-1, :3, 0])))
    all_groups = [overview, detail, overlay, labels, *bg_groups]
    def isolate(layer, allowed):
        for g in all_groups: layer.layer_collection.children[g.name].exclude = g not in allowed
    beauty = bpy.context.view_layer; beauty.name = 'Overview'; isolate(beauty, [overview])
    close = scene.view_layers.new('Ubiquitin detail'); isolate(close, [detail])
    annotations = scene.view_layers.new('Reaction annotation render'); isolate(annotations, [overlay])
    bg_layers = []
    for g in bg_groups:
        layer = scene.view_layers.new(g.name+' render'); isolate(layer, [g]); bg_layers.append(layer)
    recipe = next(r.copy() for r in json.loads((root/'assets/render-recipes.json').read_text())['recipes'] if r['id'] == a.recipe)
    recipe.update(background=profile['background'], label_min_pixels=25*a.width/1280)
    report = build_stack(scene, recipe)
    paired_composite(scene, close, annotations, bg_layers, profile['background'], a.width)
    for action in bpy.data.actions:
        for layer in action.layers:
            for strip in layer.strips:
                for bag in strip.channelbags:
                    for curve in bag.fcurves:
                        for key in curve.keyframe_points: key.interpolation = 'CONSTANT' if curve.data_path in ('hide_render', 'hide_viewport') else 'LINEAR'
    scene['asset_directory'] = str(a.asset_dir.resolve()); scene['mechanism_scope'] = source['scope']
    scene['section_method'] = 'Whole source proteins. Detail omits the entire E2, UEV and RING. Smooth framing keeps all displayed atom centers inside their panels.'
    a.output.mkdir(parents=True); scene.frame_set(round(count*.82)+1); bpy.context.view_layer.update()
    bpy.ops.wm.save_as_mainfile(filepath=str(a.output/'scene.blend'))
    result = {'kind': 'ubiquitin-net-transfer-scene', 'fps': a.fps, 'frames': count, 'preparation': source, 'asset_directory': str(a.asset_dir.resolve()),
        'transforms': {k: v.tolist() for k, v in transforms.items()}, 'actors': actors, 'annotations': marks, 'product_leader_tracks_atoms': True,
        'environment': {'objects': bg_record, 'evidence': 'Source-shaped molecular crowding backdrop; no cell type, measured abundance or common scene scale asserted.'},
        'font_zh_sha256': digest(a.font_zh), 'font_en_sha256': digest(a.font_en), 'style': a.style, 'compositor_report': report,
        'limits': source['limits']+['The two views use different magnifications. Covalent annotation strokes are enlarged for readability.', 'Full atom-derived Gaussian envelopes are rebuilt per frame; no protein surface stretching.']}
    (a.output/'manifest.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    if a.render:
        for t in ((.12, .32, .43, .60, .73, .83, .89, .95) if source.get('release_included') else (.12, .32, .43, .60, .82)):
            scene.frame_set(round(count*t)+1); scene.render.filepath = str(a.output/f'key_{scene.frame_current:04d}.png'); bpy.ops.render.render(write_still=True)
    print(json.dumps({'saved': str(a.output/'scene.blend'), 'frames': count, 'style': a.style, 'recipe': a.recipe}))


if __name__ == '__main__': main()
