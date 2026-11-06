"""Reject cross-panel passes and invalid alpha/engine/range configurations."""
import argparse
import copy
import json
from pathlib import Path
import sys

import bpy

sys.path.insert(0, str(Path(__file__).resolve().parent))
from panel_compositor import audit_panel_stack, bindings_for, build_panel_stack


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(sys.argv[sys.argv.index('--')+1:])
    scene = bpy.context.scene; base = scene.compositing_node_group
    bindings = bindings_for(scene, base); names = list(bindings)
    if len(names) < 2: raise ValueError('This test requires two independently published panels')
    scene.camera.data.dof.use_dof = False; scene.render.use_motion_blur = False
    depth = {'id': 'test-depth', 'filters': [{'type': 'depth_fog', 'panel_ranges': {n: [20., 60.] for n in names}}]}
    build_panel_stack(scene, depth, base); valid = scene.compositing_node_group
    audit_panel_stack(scene, valid); results = []
    def reject(name, expected, call):
        try: call()
        except ValueError as error:
            if expected not in str(error): raise AssertionError((name, str(error))) from error
            results.append({'case': name, 'rejected': True, 'diagnostic': str(error)})
        else: raise AssertionError('Accepted '+name)
    bad = base.copy(); del bad['display_panel_contract_version']
    reject('unqualified-layout', 'no qualified', lambda: bindings_for(scene, bad))
    bad = base.copy(); bad.nodes[bindings[names[0]]['render_node']].layer = bindings[names[1]]['view_layer']
    reject('wrong-render-layer', 'scene/layer changed', lambda: bindings_for(scene, bad))
    bad = base.copy(); target = bindings[names[0]]['target']
    bad.links.new(bad.nodes[bindings[names[1]]['render_node']].outputs['Image'], bad.nodes[target[0]].inputs[target[1]])
    reject('wrong-panel-image', 'no longer feeds', lambda: bindings_for(scene, bad))
    bad = base.copy(); changed = copy.deepcopy(bindings); changed[names[1]]['target'] = changed[names[0]]['target']
    bad['display_panel_bindings'] = json.dumps(changed)
    reject('duplicate-layout-owner', 'same layout input', lambda: bindings_for(scene, bad))
    reject('unknown-panel', 'distinct published', lambda: build_panel_stack(scene, depth, base, ['absent']))
    bad_recipe = copy.deepcopy(depth); bad_recipe['filters'][0]['panel_ranges'][names[0]] = [float('nan'), 10]
    reject('nonfinite-depth', 'finite increasing', lambda: build_panel_stack(scene, bad_recipe, base))
    bad_recipe = copy.deepcopy(depth); bad_recipe['filters'][0]['panel_ranges'][names[0]] = [40, 20]
    reject('reversed-depth', 'finite increasing', lambda: build_panel_stack(scene, bad_recipe, base))
    bad_recipe = copy.deepcopy(depth); bad_recipe['filters'][0].pop('panel_ranges')
    reject('missing-panel-depths', 'finite increasing', lambda: build_panel_stack(scene, bad_recipe, base))
    recipe = {'id': 'bad-order', 'filters': [{'type': 'contrast', 'amount': 2}, {'type': 'outline'}]}
    reject('silently-reordered-filters', 'must precede', lambda: build_panel_stack(scene, recipe, base))
    scene.render.engine = 'BLENDER_EEVEE'
    focus = {'id': 'focus', 'filters': [{'type': 'subject_focus', 'object_index': 137}]}
    reject('incompatible-engine', 'requires Cycles', lambda: build_panel_stack(scene, focus, base))
    scene.render.engine = 'CYCLES'
    reject('absent-subject', 'No visible indexed subject', lambda: build_panel_stack(scene, focus, base))
    layer = scene.view_layers[bindings[names[0]]['view_layer']]
    subject = next(obj for obj in layer.objects if obj.type == 'MESH' and not obj.hide_render and obj.visible_get(view_layer=layer))
    previous_index = subject.pass_index; subject.pass_index = 137
    reject('subject-only-in-other-panel', 'No visible indexed subject in panel '+names[1], lambda: build_panel_stack(scene, focus, base))
    subject.pass_index = previous_index
    scene.camera.data.dof.use_dof = True
    reject('blurred-beauty-unblurred-passes', 'unblurred beauty', lambda: build_panel_stack(scene, depth, base))
    scene.camera.data.dof.use_dof = False
    bad = valid.copy(); history = json.loads(bad['display_panel_filter_history'])
    link = next(item for item in history[0]['pass_links'] if item[0] == 'Depth'); target = bad.nodes[link[1][0]].inputs[link[1][1]]
    bad.links.new(bad.nodes[bindings[names[1]]['render_node']].outputs['Depth'], target)
    reject('cross-panel-depth', 'Mismatched render pass', lambda: audit_panel_stack(scene, bad))
    bad = valid.copy(); alpha = history[0]['alpha_nodes']
    bad.links.new(bad.nodes[bindings[names[1]]['render_node']].outputs['Alpha'], bad.nodes[alpha['restore']].inputs['Alpha'])
    reject('cross-panel-alpha', 'Mismatched render pass', lambda: audit_panel_stack(scene, bad))
    bad = valid.copy(); alpha = history[0]['alpha_nodes']; bad.nodes[alpha['straight']].inputs['Type'].default_value = 'To Premultiplied'
    reject('wrong-alpha-conversion', 'alpha conversion changed', lambda: audit_panel_stack(scene, bad))
    bad = valid.copy(); bad.nodes[alpha['opaque']].inputs['Alpha'].default_value = .4
    reject('nonopaque-grading-color', 'not opaque', lambda: audit_panel_stack(scene, bad))
    layer = scene.view_layers[bindings[names[0]]['view_layer']]; layer.use_pass_z = False
    reject('disabled-depth-pass', 'pass disabled', lambda: audit_panel_stack(scene, valid)); layer.use_pass_z = True
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps({'passed': True, 'negative_cases': len(results), 'checks': results}, indent=2)+'\n')
    print(json.dumps({'passed': True, 'negative_cases': len(results)}))


if __name__ == '__main__': main()
