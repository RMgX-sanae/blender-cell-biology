"""Reject actual two-cycle appearance corruptions; never save mutated scenes."""
import argparse
import hashlib
import json
import sys
import tempfile
from pathlib import Path
import bpy
sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_transcription_appearance import audit
from transcription_appearance import read_profile


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(sys.argv[sys.argv.index('--')+1:])
    path = Path(bpy.data.filepath); before = hashlib.sha256(path.read_bytes()).hexdigest()
    baseline = audit(bpy.context.scene)
    if baseline['cycles'] != 2:
        raise ValueError('Run this extension regression against a two-cycle appearance scene')
    invalid = [{'surface_roles': {'peptide_symbol': {}}}, {'palette': {'peptide': 'FFFFFF'}},
               {'camera': {'center_offset_nm': [0, float('nan')]}}, {'camera': {'scale_multiplier': 0}},
               {'context_blur_pixels_1280': -1}, {'representation': 'ribbon'}]
    with tempfile.TemporaryDirectory() as folder:
        profile = Path(folder)/'profile.json'
        for case in invalid:
            profile.write_text(json.dumps({'render_recipe':'neutral', **case}))
            try: read_profile(profile)
            except ValueError: pass
            else: raise AssertionError('Invalid transcription request accepted')
    cases = {'phosphate_role':'role assignment', 'phosphate_color':'Semantic color',
             'shader':'Saved shader', 'color_handoff':'color handoff', 'camera_path':'camera path',
             'new_caption_layer':'outside protected typography', 'context_layer':'leak into context blur',
             'caption_timing':'Caption event timing', 'font':'Packed bilingual font',
             'context_blur':'Context blur differs'}
    results = []
    for case, diagnostic in cases.items():
        bpy.ops.wm.open_mainfile(filepath=str(path)); scene = bpy.context.scene
        record = json.loads(scene['transcription_appearance'])
        core = scene.objects['Incoming CTP core → RNA C20']
        phosphate = scene.objects['UTP beta/gamma phosphates → second PPi'].data.materials[0]
        caption = scene.objects['Second UTP label']
        if case == 'phosphate_role': phosphate['surface_role'] = 'protein'
        if case == 'phosphate_color': phosphate.node_tree.nodes['Principled BSDF'].inputs['Base Color'].default_value[0] += .1
        if case == 'shader': phosphate.node_tree.nodes['Principled BSDF'].inputs['Roughness'].default_value = .95
        if case == 'color_handoff': core.data.materials[0].node_tree.animation_data_clear()
        if case == 'camera_path': scene.camera.animation_data_clear()
        if case == 'new_caption_layer':
            for group in list(caption.users_collection): group.objects.unlink(caption)
            bpy.data.collections['02 Native modeled DNA and RNA'].objects.link(caption)
        if case == 'context_layer':
            scene.view_layers[record['context_layer']].layer_collection.children['02 Native modeled DNA and RNA'].exclude = False
        if case == 'caption_timing':
            stage = scene.objects['Cycle stage 1']; stage.animation_data_clear(); stage.hide_render = False
        if case == 'font': caption.data.font = bpy.data.curves.new('Unpacked font mutation', 'FONT').font
        if case == 'context_blur':
            next(n for n in scene.compositing_node_group.nodes if n.label == 'Whole Pol II context defocus').inputs['Size'].default_value = (30, 30)
        bpy.context.view_layer.update()
        try: audit(scene)
        except ValueError as error:
            assert diagnostic in str(error), (case, str(error))
            results.append({'case':case, 'rejected':True, 'diagnostic':str(error)})
        else: raise AssertionError('Accepted mutation '+case)
    bpy.ops.wm.open_mainfile(filepath=str(path)); assert audit(bpy.context.scene) == baseline
    assert hashlib.sha256(path.read_bytes()).hexdigest() == before
    report = {'passed':True, 'invalid_profiles':len(invalid), 'native_cases':results,
              'restored_positive':True, 'saved_file_unchanged':True}
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'passed':True, 'invalid_profiles':len(invalid), 'native_cases':len(results)}))


if __name__ == '__main__':main()
