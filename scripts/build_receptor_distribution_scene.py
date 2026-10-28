"""Render whole source receptors on an oriented bilayer, with selectable looks.

Run in Blender. Preparation and all generated assets stay outside the skill.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Matrix, Vector, Euler

sys.path.insert(0, str(Path(__file__).resolve().parent))
import detail_scene as ds
import style_nodes as sn
import cellular_environment as environment
from compositor_stack import build_stack


def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def lipids(data, group, head_material, tail_material):
    obj = sn.point_object('Intact two-tail lipid instances', data['lipid_heads'], group)
    rotations = [Vector(direction).to_track_quat('Z', 'Y').to_euler() for direction in data['lipid_directions']]
    attr = obj.data.attributes.new('lipid_rotation', 'FLOAT_VECTOR', 'POINT')
    attr.data.foreach_set('vector', np.asarray(rotations, dtype=np.float32).reshape(-1))
    tree, source, output = sn.geometry_tree(obj, 'GN_Intact_lipid_leaflets')
    sphere = sn.node(tree, 'GeometryNodeMeshIcoSphere'); sphere.inputs['Radius'].default_value = .28
    sphere.inputs['Subdivisions'].default_value = 2
    head = sn.node(tree, 'GeometryNodeSetMaterial'); head.inputs['Material'].default_value = head_material
    tree.links.new(sphere.outputs['Mesh'], head.inputs['Geometry'])
    cylinder = sn.node(tree, 'GeometryNodeMeshCylinder'); cylinder.inputs['Vertices'].default_value = 8
    cylinder.inputs['Radius'].default_value = .065; cylinder.inputs['Depth'].default_value = 1.44
    tails = sn.node(tree, 'GeometryNodeSetMaterial'); tails.inputs['Material'].default_value = tail_material
    tree.links.new(cylinder.outputs['Mesh'], tails.inputs['Geometry'])
    joined = sn.node(tree, 'GeometryNodeJoinGeometry'); tree.links.new(head.outputs['Geometry'], joined.inputs['Geometry'])
    for offset in (-.11, .11):
        move = sn.node(tree, 'GeometryNodeTransform'); move.inputs['Translation'].default_value = (offset, 0, .88)
        tree.links.new(tails.outputs['Geometry'], move.inputs['Geometry'])
        tree.links.new(move.outputs['Geometry'], joined.inputs['Geometry'])
    smooth = sn.node(tree, 'GeometryNodeSetShadeSmooth'); tree.links.new(joined.outputs['Geometry'], smooth.inputs['Geometry'])
    points = sn.node(tree, 'GeometryNodeMeshToPoints'); points.mode = 'VERTICES'
    tree.links.new(source, points.inputs['Mesh'])
    rotation = sn.node(tree, 'GeometryNodeInputNamedAttribute'); rotation.data_type = 'FLOAT_VECTOR'
    rotation.inputs['Name'].default_value = 'lipid_rotation'
    instances = sn.node(tree, 'GeometryNodeInstanceOnPoints')
    tree.links.new(points.outputs['Points'], instances.inputs['Points'])
    tree.links.new(smooth.outputs['Geometry'], instances.inputs['Instance'])
    tree.links.new(rotation.outputs['Attribute'], instances.inputs['Rotation'])
    tree.links.new(instances.outputs['Instances'], output)
    obj['evidence'] = 'Intact schematic lipid head plus two tails; no lipid species or dynamics inferred'
    return obj


def backdrop(asset_dir, basis):
    source = json.loads((asset_dir/'manifest.json').read_text())
    if source['kind'] != 'schematic-cytosolic-environment': raise ValueError('Prepare cellular environment assets first')
    groups = [ds.collection('Cytosol '+name) for name in ('far', 'middle', 'near')]
    records, data, meshes = [], {}, {}
    for name, color in [('ribosome', '#638e8d'), ('actin', '#8b829a'), ('ubiquitin', '#728d9b')]:
        data[name] = dict(np.load(asset_dir/(name+'.npz')))
        mesh = bpy.data.meshes.new('Intact environment '+name)
        mesh.from_pydata(data[name]['vertices'].tolist(), [], data[name]['faces'].tolist()); mesh.update()
        mesh.materials.append(ds.material('Background '+name, color, .65, 0))
        for polygon in mesh.polygons: polygon.use_smooth = True
        meshes[name] = mesh
    def place(name, layer, center, angles, local=None):
        matrix = np.eye(4); matrix[:3, :3] = basis@np.array(Euler(angles).to_matrix()); matrix[:3, 3] = basis@center
        if local is not None: matrix = matrix@local
        obj = bpy.data.objects.new(f'Context {name} {len(records):03d}', meshes[name]); groups[layer].objects.link(obj)
        obj.matrix_world = Matrix(matrix); obj['environment_only'] = True
        obj['source_structure'] = source['sources'][name]['id']
        records.append({'object': obj.name, 'asset': name, 'matrix_nm': matrix.tolist(), 'depth_layer': layer})
    for center, angles in [([-72, 12, -70], [.5, .7, .6]), ([58, 32, -95], [1.8, .1, -.3]), ([76, -37, -65], [.8, 1.5, .2])]:
        place('ribosome', 0, np.array(center), angles)
    for center, angles in [([-67, -34, -40], [.2, 1.3, -.45]), ([59, 38, -46], [1.3, .8, .35])]:
        for repeat in range(-6, 7): place('actin', 1, np.array(center), angles, np.linalg.matrix_power(data['actin']['repeat_nm'], repeat))
    rng = np.random.default_rng(741)
    for index in range(22):
        side = -1 if index % 2 else 1
        center = np.array([side*rng.uniform(66, 92), rng.uniform(-49, 44), rng.uniform(-70, -32)])
        place('ubiquitin', 2 if index > 18 else 1, center, rng.uniform(0, 6.28, 3))
    return groups, {'source_manifest': source, 'objects': records,
        'layout': 'Source meshes at native nm scale, illustrative cytosolic placements; not measured concentration.',
        'defocus': 'Separate artistic depth-layer blur; not an optical PSF or blur of the subject.'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--distribution-dir', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--font-zh', type=Path, required=True); p.add_argument('--font-en', type=Path, required=True)
    p.add_argument('--environment-dir', type=Path)
    p.add_argument('--style', choices=['dark', 'light'], default='dark')
    p.add_argument('--recipe', default='neutral')
    p.add_argument('--engine', choices=['eevee', 'cycles'], default='eevee')
    p.add_argument('--width', type=int, default=1280); p.add_argument('--samples', type=int, default=32)
    p.add_argument('--render', action='store_true')
    a = p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.output.exists(): p.error('Use a fresh scene output directory')
    if not a.font_zh.is_file() or not a.font_en.is_file(): p.error('Supply Songti and Times New Roman fonts')
    root = Path(__file__).resolve().parents[1]
    record = json.loads((a.distribution_dir/'manifest.json').read_text())
    if record['kind'] != 'source-receptor-distribution': p.error('Prepare a source receptor distribution first')
    if digest(a.distribution_dir/'distribution.npz') != record['data_sha256']: p.error('Distribution data hash mismatch')
    assets = Path(record['asset_dir'])
    if digest(assets/'manifest.json') != record['source_manifest_sha256']: p.error('Source manifest changed after preparation')
    for chain, expected in record['chain_hashes'].items():
        if digest(assets/f'chain_{chain}.npz') != expected: p.error('Source chain changed: '+chain)
    data = dict(np.load(a.distribution_dir/'distribution.npz'))
    bpy.ops.object.select_all(action='SELECT'); bpy.ops.object.delete(use_global=False)
    scene = bpy.context.scene
    receptor_group = ds.collection('Source receptor complexes'); membrane_group = ds.collection('Intact membrane lipids')
    labels = ds.collection('Distribution labels')
    profile_file = 'style-molecular-reference.json' if a.style == 'dark' else 'style-studio.json'
    profile = json.loads((root/'assets'/profile_file).read_text()); profile['camera'] = {}; profile['bump'] = 0
    profile['background'] = '#071f2b' if a.style == 'dark' else '#edf0e9'
    a.fps, a.frames, a.frame, a.animate, a.compositor = 24, 1, 1, False, 'clean'
    display = record['config']['display']; span = display['camera_scale_nm']
    camera = ds.studio(profile, a, {'position': display['camera_position_nm'], 'target': display['camera_target_nm'], 'scale': span, 'light_scale': span/29.375})
    camera.data.dof.use_dof = False
    if hasattr(scene, 'eevee'): scene.eevee.taa_render_samples = a.samples
    bpy.context.view_layer.update(); basis = np.array(camera.matrix_world.to_3x3())
    colors = {c: profile['palette'][key] for c, key in [('R', 'receptor'), ('A', 'alpha'), ('B', 'beta'), ('G', 'gamma')]}
    # Other explicitly qualified assets can supply chain colors in their config.
    colors.update(record['config'].get('chain_colors', {}))
    objects = []
    for chain in record['config']['chains']:
        mat = ds.material('Source chain '+chain, colors[chain], .54, 0)
        template, source = ds.load_chain(assets, chain, mat, receptor_group, name=f'Receptor_000_chain_{chain}')
        for i, (position, rotation) in enumerate(zip(data['positions'], data['rotations'])):
            obj = template if i == 0 else ds.duplicate(template, f'Receptor_{i:03d}_chain_{chain}', receptor_group)
            matrix = np.eye(4); matrix[:3, :3] = rotation; matrix[:3, 3] = position
            obj.matrix_world = Matrix(matrix); obj['instance_id'] = i; obj.pass_index = 1
            objects.append({'object': obj.name, 'chain': chain, 'instance': i, 'matrix': matrix.tolist()})
    membrane_palette = display.get('membrane_palette', {}).get(a.style, {})
    lipids(data, membrane_group, ds.material('Lipid heads', membrane_palette.get('heads', profile['palette']['heads']), .65, 0),
           ds.material('Lipid tails', profile['palette']['tails'], .65, 0))
    shell = ds.mesh('Continuum bilayer core', data['shell_vertices'].tolist(), data['shell_faces'].tolist(),
        [ds.material('Membrane core fill', membrane_palette.get('core', profile['palette']['heads']), .68, 0)], membrane_group)
    shell['evidence'] = 'Schematic continuum core with receptor exclusion; complete lipids and proteins remain intact'
    font_zh = bpy.data.fonts.load(str(a.font_zh)); font_en = bpy.data.fonts.load(str(a.font_en))
    font_zh.pack(); font_en.pack()
    def label(body, x, y, size):
        curve = bpy.data.curves.new(body, 'FONT'); curve.body = body; curve.font = font_zh; curve.font_bold = font_en
        for char, form in zip(body, curve.body_format): form.use_bold = ord(char) < 128 or char in 'αβγ'
        curve.size = size; curve.align_x = 'CENTER'
        if a.style == 'light': curve.offset = size*.01
        obj = bpy.data.objects.new(body, curve); labels.objects.link(obj)
        matrix = np.eye(4); matrix[:3, :3] = basis; matrix[:3, 3] = np.array(display['camera_target_nm'])+basis@np.array([x*span, y*span, .48*span])
        obj.matrix_world = Matrix(matrix)
    mode = record['config']['distribution']
    label({'uniform': '均匀分布（示意）', 'clustered': '成簇分布（示意）', 'polarized': '极性分布（示意）'}[mode], 0, .272, .03*span)
    label(display['context_label'], 0, -.263, .0225*span)
    for x, item in zip(np.linspace(-.3, .3, len(display['legend'])), display['legend']):
        chain, text = item['chain'], item['label']
        label(text, x, -.294, .021*span)
        dot_position = np.array(display['camera_target_nm'])+basis@np.array([(x-.045)*span, -.2865*span, .48*span])
        dot = sn.point_object(text+' legend', [dot_position], labels)
        sn.instance_spheres(dot, sn.material(text+' identity', {'material': 'paint'}, colors[chain], 'paint'), radius=.0045*span)
    env_groups, env_record = [], None
    if a.environment_dir: env_groups, env_record = backdrop(a.environment_dir, basis)
    recipes = json.loads((root/'assets/render-recipes.json').read_text())['recipes']
    recipe = next((r.copy() for r in recipes if r['id'] == a.recipe), None)
    if recipe is None: p.error('Unknown compositor recipe')
    recipe['background'] = profile['background']
    recipe['label_min_pixels'] = 26*a.width/1280
    result = build_stack(scene, recipe)
    # The existing filter chain consumes the composed depth backdrop; the main
    # molecules and protected text never pass through the background blur.
    if env_groups:
        beauty = scene.view_layers[scene['compositor_beauty_layer']]
        layers = []
        for group in env_groups:
            beauty.layer_collection.children[group.name].exclude = True
            layer = scene.view_layers.new(group.name+' render')
            for child in layer.layer_collection.children:
                child.exclude = child.collection not in [group, camera.users_collection[0]]
            layers.append(layer)
        tree = scene.compositing_node_group
        image = environment.composite(tree, layers, profile['background'], a.width)
        base = next(n for n in tree.nodes if n.label == 'Straight background plus premultiplied rendered foreground')
        tree.links.new(image, base.inputs['Background'])
    ds.save_and_render(a, {'kind': 'source-receptor-distribution-scene', 'preparation': record,
        'objects': objects, 'environment': env_record, 'compositor_report': result,
        'font_zh_sha256': digest(a.font_zh), 'font_en_sha256': digest(a.font_en),
        'display': 'Static arrangement of complete source complexes and intact procedural lipids. No protein slicing or anisotropic scaling.'})


if __name__ == '__main__': main()
