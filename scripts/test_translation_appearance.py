"""Reject broken requests and actual saved-scene mutations without saving them."""
import argparse
import hashlib
import json
import sys
import tempfile
from pathlib import Path
import bpy
sys.path.insert(0, str(Path(__file__).resolve().parent))
from translation_appearance import read_profile
from audit_translation_appearance import audit


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(sys.argv[sys.argv.index('--')+1:])
    path = Path(bpy.data.filepath); before = hashlib.sha256(path.read_bytes()).hexdigest()
    baseline = audit(bpy.context.scene)
    invalid = [{'camera': {'type': 'PERSP'}}, {'camera': {'scale_multiplier': 0}},
               {'camera': {'center_offset_nm': [0, float('nan')]}}, {'surface_roles': {'unknown': {}}},
               {'palette': {'unknown': 'FFFFFF'}}, {'palette': {'rna': 'wrong'}},
               {'lighting': []}, {'exposure': float('inf')}, {'context_blur_pixels_1280': -1},
               {'representation': 'spacefill'}]
    with tempfile.TemporaryDirectory() as folder:
        target = Path(folder)/'profile.json'
        for case in invalid:
            target.write_text(json.dumps({'render_recipe': 'neutral', **case}))
            try: read_profile(target)
            except ValueError: pass
            else: raise AssertionError('Accepted invalid profile '+str(case))
    cases = {'role': 'role assignment', 'color': 'Semantic color', 'shader': 'Saved shader',
             'smoothing': 'smoothing', 'camera': 'camera pan', 'scale': 'camera scale',
             'light': 'Light power', 'font': 'Packed bilingual font', 'routing': 'Mixed font routing',
             'ink': 'Text ink', 'caption': 'Caption inventory', 'blur': 'Context blur differs',
             'context': 'leak into context blur', 'typography': 'Foreground leaks',
             'crop': 'Caption outside frame', 'text_output': 'Labels-over-Beauty'}
    results = []
    for case, diagnostic in cases.items():
        bpy.ops.wm.open_mainfile(filepath=str(path)); scene = bpy.context.scene
        record = json.loads(scene['translation_appearance']); tree = scene.compositing_node_group
        obj = next(o for o in scene.objects if o.get('source_chain') == 'LA')
        mat = obj.data.materials[0]; caption = scene.objects['Stage 1']
        if case == 'role': mat['surface_role'] = 'peptide_symbol'
        if case == 'color': scene.objects['tRNA identity At'].data.materials[0].node_tree.nodes['Principled BSDF'].inputs['Base Color'].default_value[0] += .1
        if case == 'shader': mat.node_tree.nodes['Principled BSDF'].inputs['Roughness'].default_value = .9
        if case == 'smoothing': obj.modifiers[0].factor = .8
        if case == 'camera': scene.camera.location.x += 1
        if case == 'scale': scene.camera.data.ortho_scale *= .8
        if case == 'light': scene.objects[record['studio']['lighting'][0]['name']].data.energy += 1000
        if case == 'font': caption.data.font = caption.data.font_bold
        if case == 'routing':
            index = caption.data.body.index('t'); caption.data.body_format[index].use_bold = False
        if case == 'ink': next(n for n in caption.data.materials[0].node_tree.nodes if n.type == 'EMISSION').inputs['Color'].default_value[0] += .1
        if case == 'caption': caption.data.body += '?'
        if case == 'blur': next(n for n in tree.nodes if n.label == 'Whole ribosomal context defocus').inputs['Size'].default_value = (50, 50)
        if case == 'context': scene.view_layers[record['context_layer']].layer_collection.children['02 Intact tRNAs — identity retained throughout'].exclude = False
        if case == 'typography': scene.view_layers[scene['compositor_label_layer']].layer_collection.children[record['context_collection']].exclude = False
        if case == 'crop': caption.location.x += 100
        if case == 'text_output':
            out = next(n for n in tree.nodes if n.type == 'GROUP_OUTPUT' and n.is_active_output)
            render = next(n for n in tree.nodes if n.type == 'R_LAYERS' and n.layer == scene['compositor_beauty_layer'])
            tree.links.new(render.outputs['Image'], out.inputs['Image'])
        bpy.context.view_layer.update()
        try: audit(scene)
        except ValueError as error:
            assert diagnostic in str(error), (case, str(error))
            results.append({'case': case, 'rejected': True, 'diagnostic': str(error)})
        else: raise AssertionError('Accepted scene mutation '+case)
    bpy.ops.wm.open_mainfile(filepath=str(path))
    assert audit(bpy.context.scene) == baseline
    assert hashlib.sha256(path.read_bytes()).hexdigest() == before
    result = {'passed': True, 'invalid_profiles': len(invalid), 'native_cases': results,
              'restored_positive': True, 'saved_file_unchanged': True}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'passed': True, 'invalid_profiles': len(invalid), 'native_cases': len(results)}))


if __name__ == '__main__': main()
