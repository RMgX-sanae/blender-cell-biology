"""Apply layer-dependent filters before panels are framed, then grade the layout.

Each published panel binds an actual Render Layers image to a specific layout
input. Depth/normal/object passes never cross between these bindings. Geometry,
camera, annotations, background layers and existing panel masks stay intact.
"""
import copy
import json
from pathlib import Path
import sys

import bpy

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compositor_stack import build_stack, filter_image
from style_nodes import node


PASSES = {'outline': ['Normal', 'Alpha'], 'depth_fog': ['Depth', 'Alpha'],
          'ambient_occlusion': ['Ambient Occlusion'], 'subject_focus': ['Object Index']}
FLAGS = {'Normal': 'use_pass_normal', 'Depth': 'use_pass_z',
         'Ambient Occlusion': 'use_pass_ambient_occlusion', 'Object Index': 'use_pass_object_index'}


def socket_ref(socket, output):
    return [socket.node.name, list(socket.node.outputs if output else socket.node.inputs).index(socket)]


def publish_panel_sources(tree, panels):
    """Register named (Render Layers node, downstream image input) pairs."""
    bindings = {}
    for name, (render, target) in panels.items():
        if render.bl_idname != 'CompositorNodeRLayers' or target.node.id_data != tree:
            raise ValueError('Publish a Render Layers node and an input in this compositor')
        if len(target.links) != 1 or target.links[0].from_socket != render.outputs['Image']:
            raise ValueError('Publish the direct image link before modifying a panel')
        bindings[name] = {'render_node': render.name, 'view_layer': render.layer,
                          'target': socket_ref(target, False), 'image': socket_ref(render.outputs['Image'], True)}
    if not bindings: raise ValueError('Publish at least one explicit panel binding')
    tree['display_panel_bindings'] = json.dumps(bindings)
    tree['display_panel_contract_version'] = 1


def bindings_for(scene, tree):
    if tree.get('display_panel_contract_version') != 1:
        raise ValueError('Scene has no qualified per-panel source bindings; rebuild with a supporting generator')
    bindings = json.loads(tree['display_panel_bindings']); targets = set()
    for name, item in bindings.items():
        render = tree.nodes.get(item['render_node'])
        if (render is None or render.bl_idname != 'CompositorNodeRLayers'
                or render.layer != item['view_layer'] or render.layer not in scene.view_layers
                or render.scene not in (None, scene)):
            raise ValueError('Panel render scene/layer changed: '+name)
        target_key = tuple(item['target'])
        if target_key in targets: raise ValueError('Two panels cannot own the same layout input')
        targets.add(target_key)
        try:
            target = tree.nodes[item['target'][0]].inputs[item['target'][1]]
            image = tree.nodes[item['image'][0]].outputs[item['image'][1]]
        except (KeyError, IndexError):
            raise ValueError('Panel image/target socket missing: '+name) from None
        if len(target.links) != 1 or target.links[0].from_socket != image:
            raise ValueError('Panel image no longer feeds its declared layout input: '+name)
    return bindings


def audit_panel_stack(scene, tree):
    bindings = bindings_for(scene, tree)
    history = json.loads(tree.get('display_panel_filter_history', '[]'))
    for item in history:
        render = tree.nodes[bindings[item['panel']]['render_node']]
        for pass_name in item['required_passes']:
            if pass_name in FLAGS and not getattr(scene.view_layers[render.layer], FLAGS[pass_name]):
                raise ValueError('Panel render pass disabled: '+pass_name)
        alpha = item['alpha_nodes']
        for key, kind, value in [('straight', 'CompositorNodePremulKey', 'To Straight'),
                                 ('opaque', 'CompositorNodeSetAlpha', 'Replace Alpha'),
                                 ('restore', 'CompositorNodeSetAlpha', 'Replace Alpha'),
                                 ('premultiplied', 'CompositorNodePremulKey', 'To Premultiplied')]:
            step = tree.nodes.get(alpha[key])
            if step is None or step.bl_idname != kind or step.inputs['Type'].default_value != value:
                raise ValueError('Panel alpha conversion changed: '+item['panel']+' '+key)
        opaque = tree.nodes[alpha['opaque']]
        if opaque.inputs['Alpha'].is_linked or opaque.inputs['Alpha'].default_value != 1:
            raise ValueError('Panel color is not opaque during grading')
        for target_ref, image_ref in item['alpha_image_links']:
            target = tree.nodes[target_ref[0]].inputs[target_ref[1]]
            image = tree.nodes[image_ref[0]].outputs[image_ref[1]]
            if len(target.links) != 1 or target.links[0].from_socket != image:
                raise ValueError('Panel alpha image routing changed: '+item['panel'])
        for pass_name, target_ref in item['pass_links']:
            target = tree.nodes[target_ref[0]].inputs[target_ref[1]]
            if pass_name not in render.outputs or len(target.links) != 1 or target.links[0].from_socket != render.outputs[pass_name]:
                raise ValueError('Mismatched render pass in panel '+item['panel']+': '+pass_name)
        for name in item['nodes']:
            if name not in tree.nodes: raise ValueError('Panel filter node missing: '+name)
    return {'passed': True, 'panels': {name: item['view_layer'] for name, item in bindings.items()},
            'filter_applications': len(history), 'pass_links_checked': sum(len(i['pass_links']) for i in history),
            'limits': 'Graph/source alignment only; inspect rendered alpha edges, color identity, labels and animation separately.'}


def build_panel_stack(scene, recipe, source, selected=None):
    """Apply the leading pass-dependent operators per panel, then global filters.

    Fog requires explicit per-panel ranges, in each rendered layer's scene units.
    These artistic ranges are fixed for the shot; they are not optical constants.
    """
    if source is None: raise ValueError('Per-panel filtering requires an existing qualified compositor')
    bindings = bindings_for(scene, source)
    selected = list(bindings) if selected is None else selected
    if not selected or len(selected) != len(set(selected)) or set(selected)-set(bindings):
        raise ValueError('Select distinct published panel names')
    prefix = []; suffix = []; after_image_filter = False
    for settings in recipe.get('filters', []):
        if settings['type'] in PASSES:
            if after_image_filter:
                raise ValueError('Per-panel pass filters must precede image-only filters; do not silently reorder a recipe')
            prefix.append(settings)
        else:
            after_image_filter = True; suffix.append(settings)
    if any(f['type'] == 'subject_focus' for f in prefix) and scene.render.engine != 'CYCLES':
        raise ValueError('Per-panel object-index isolation requires Cycles')
    if prefix and (scene.camera.data.dof.use_dof or getattr(scene.render, 'use_motion_blur', False)):
        raise ValueError('Per-panel data passes require unblurred beauty; use separate composited background defocus instead of camera DOF/motion blur')
    for f in prefix:
        if f['type'] != 'subject_focus': continue
        index = int(f.get('object_index', 1))
        for name in selected:
            layer = scene.view_layers[bindings[name]['view_layer']]
            if not any(obj.pass_index == index and not obj.hide_render and obj.visible_get(view_layer=layer)
                       and obj.type in {'MESH', 'CURVE', 'SURFACE', 'META', 'FONT'} for obj in layer.objects):
                raise ValueError('No visible indexed subject in panel '+name)
    for f in prefix:
        if f['type'] != 'depth_fog': continue
        ranges = f.get('panel_ranges', {})
        for name in selected:
            values = ranges.get(name)
            if (not isinstance(values, list) or len(values) != 2
                    or any(isinstance(v, bool) or not isinstance(v, (float, int)) for v in values)
                    or not all(float('-inf') < v < float('inf') for v in values)
                    or not 0 <= values[0] < values[1]):
                raise ValueError('Provide a finite increasing depth_fog panel_ranges entry for '+name)
    tree = source.copy(); tree.name = 'COMP_PANELS_'+recipe['id']; tree.use_fake_user = True
    history = json.loads(tree.get('display_panel_filter_history', '[]')); records = []
    for name in selected:
        if not prefix: continue
        item = bindings[name]; render = tree.nodes[item['render_node']]
        required = {p for f in prefix for p in PASSES[f['type']]}
        for pass_name in required:
            if pass_name in FLAGS: setattr(scene.view_layers[render.layer], FLAGS[pass_name], True)
        filters = copy.deepcopy(prefix)
        for f in filters:
            if f['type'] == 'depth_fog': f['range'] = f.pop('panel_ranges')[name]
        old_nodes = {n.name for n in tree.nodes}
        image = tree.nodes[item['image'][0]].outputs[item['image'][1]]
        straight = node(tree, 'CompositorNodePremulKey', name+' unpremultiply')
        straight.inputs['Type'].default_value = 'To Straight'; tree.links.new(image, straight.inputs['Image'])
        opaque = node(tree, 'CompositorNodeSetAlpha', name+' opaque color for filters')
        opaque.inputs['Type'].default_value = 'Replace Alpha'; opaque.inputs['Alpha'].default_value = 1
        tree.links.new(straight.outputs[0], opaque.inputs['Image'])
        filtered, actual = filter_image(scene, tree, opaque.outputs[0], render.outputs, filters, panel=True)
        restore = node(tree, 'CompositorNodeSetAlpha', name+' restore source coverage')
        restore.inputs['Type'].default_value = 'Replace Alpha'
        tree.links.new(filtered, restore.inputs['Image']); tree.links.new(render.outputs['Alpha'], restore.inputs['Alpha'])
        premultiplied = node(tree, 'CompositorNodePremulKey', name+' premultiply for layout')
        premultiplied.inputs['Type'].default_value = 'To Premultiplied'; tree.links.new(restore.outputs[0], premultiplied.inputs['Image'])
        target = tree.nodes[item['target'][0]].inputs[item['target'][1]]
        tree.links.new(premultiplied.outputs[0], target)
        item['image'] = socket_ref(premultiplied.outputs[0], True)
        added = [n.name for n in tree.nodes if n.name not in old_nodes]
        pass_links = [[link.from_socket.name, socket_ref(link.to_socket, False)] for link in tree.links
                      if link.from_node == render and link.to_node.name in added]
        record = {'panel': name, 'view_layer': render.layer, 'filters': actual, 'nodes': added,
                  'required_passes': sorted(required), 'pass_links': pass_links,
                  'alpha_nodes': {key: step.name for key, step in [('straight', straight), ('opaque', opaque),
                                   ('restore', restore), ('premultiplied', premultiplied)]},
                  'alpha_image_links': [[socket_ref(target, False), socket_ref(source, True)] for target, source in
                      [(straight.inputs['Image'], image), (opaque.inputs['Image'], straight.outputs[0]),
                       (restore.inputs['Image'], filtered), (premultiplied.inputs['Image'], restore.outputs[0])]]}
        history.append(record); records.append(record)
    tree['display_panel_bindings'] = json.dumps(bindings)
    tree['display_panel_filter_history'] = json.dumps(history)
    report = build_stack(scene, {**recipe, 'filters': suffix}, source=tree)
    report.update(scope='panels_then_layout', panels=records,
                  panel_check=audit_panel_stack(scene, scene.compositing_node_group),
                  requires_color_legend_review=report['requires_color_legend_review'] or any(f['type'] in {'subject_focus', 'depth_fog'} for f in prefix),
                  panel_alpha='Unpremultiply, filter opaque color, restore original coverage, premultiply before unchanged layout')
    return report


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Check saved panel-pass routing and alpha handling')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--')+1:])
    if args.output.exists(): parser.error('Use a fresh audit output')
    result = audit_panel_stack(bpy.context.scene, bpy.context.scene.compositing_node_group)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))
