"""Reject actual fork appearance corruptions without saving mutated scenes."""
import argparse
import hashlib
import json
import sys
import tempfile
from pathlib import Path
import bpy
sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_replication_appearance import audit
from replication_appearance import read_profile


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(sys.argv[sys.argv.index('--')+1:])
    path = Path(bpy.data.filepath); before = hashlib.sha256(path.read_bytes()).hexdigest()
    baseline = audit(bpy.context.scene)
    invalid = [{'surface_roles':{'phosphate':{}}}, {'palette':{'rna':'FFFFFF'}},
               {'camera':{'center_offset_nm':[0, float('nan')]}}, {'camera':{'scale_multiplier':0}},
               {'context_blur_pixels_1280':-1}, {'representation':'ribbon'}]
    with tempfile.TemporaryDirectory() as folder:
        profile = Path(folder)/'profile.json'
        for case in invalid:
            profile.write_text(json.dumps({'render_recipe':'neutral', **case}))
            try: read_profile(profile)
            except ValueError: pass
            else: raise AssertionError('Invalid replication request accepted')
    cases = {'effective_slot':'Semantic color', 'mesh_slot':'slot link', 'primer_role':'role assignment',
             'shader':'Saved shader', 'guide_light':'must be unlit', 'camera':'camera pan',
             'context_layer':'leak into context blur', 'caption_timing':'caption timing',
             'font':'Packed bilingual font', 'context_blur':'Context blur differs'}
    results = []
    for case, diagnostic in cases.items():
        bpy.ops.wm.open_mainfile(filepath=str(path)); scene = bpy.context.scene
        record = json.loads(scene['replication_appearance'])
        leading = scene.objects['leading nucleotide 00']
        if case == 'effective_slot': leading.material_slots[0].material = scene.objects['parent_A nucleotide 00'].material_slots[0].material
        if case == 'mesh_slot': leading.material_slots[0].link = 'DATA'
        if case == 'primer_role': scene.objects['lagging nucleotide 09'].material_slots[0].material['surface_role'] = 'protein'
        if case == 'shader': scene.objects['CMG complete chain 2'].material_slots[0].material.node_tree.nodes['Principled BSDF'].inputs['Roughness'].default_value = .95
        if case == 'guide_light': scene.objects['Leading synthesis direction'].data.materials[0] = leading.material_slots[0].material
        if case == 'camera': scene.camera.location.x += .2
        if case == 'context_layer': scene.view_layers[record['context_layer']].layer_collection.children['02 Newly synthesized DNA and primer symbols'].exclude = False
        if case == 'caption_timing':
            obj = scene.objects['Lagging fragment 2 3′']; obj.animation_data_clear(); obj.hide_render = False
        if case == 'font': scene.objects['Primer legend'].data.font = bpy.data.curves.new('Unpacked font mutation', 'FONT').font
        if case == 'context_blur': next(n for n in scene.compositing_node_group.nodes if n.label == 'Whole CMG context defocus').inputs['Size'].default_value = (30, 30)
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
