"""Check saved AMSH-LP identities, all displayed poses, oxygen tracking and labels."""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_phosphotransfer_scene import vertices
from build_ubiquitin_transfer_scene import digest
from reaction_ledger import validate


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.output.exists(): p.error('Use a new audit output')
    scene = bpy.context.scene; root = Path(bpy.data.filepath).parent
    record = json.loads((root/'manifest.json').read_text()); source = record['preparation']
    assets = Path(record['asset_directory']); frames = source['frames']
    assert source['passed'] and source['checks']['passed']
    assert scene.frame_end == frames and scene.render.fps == source['fps']
    for name, expected in source['data_sha256'].items(): assert digest(assets/name) == expected, name
    book = validate(json.loads((assets/'reaction-ledger.json').read_text()))
    assert book['mapped_atoms'] == 2569
    assert book['formed_bonds_in_reactant_ids'] == [('distal:76:C','water:517:O')]
    data = np.load(assets/'motion.npz'); groups = {'overview':[], 'detail':[]}
    for name in record['actors']:
        obj = scene.objects[name]; groups[obj['view']].append(obj)
        assert not any(mod.type == 'BOOLEAN' for mod in obj.modifiers), name
        assert json.loads(obj['atom_ids']) == source['atom_ids'][obj['component']], name
        for mod in obj.modifiers: mod.show_viewport = False
    identities = [identity for obj in groups['overview'] for identity in json.loads(obj['atom_ids'])]
    assert len(identities) == len(set(identities)) == 2569
    assert set(obj['component'] for obj in groups['overview']) == set(source['atom_ids'])
    assert set(obj['component'] for obj in groups['detail']) == {'distal','proximal','water'}
    maximum = 0.; checked = 0; final = {}
    for view, layer in [('overview','Overview'), ('detail','Hydrolysis detail')]:
        bpy.context.window.view_layer = scene.view_layers[layer]
        for frame in range(1, frames+1):
            scene.frame_set(frame)
            for obj in groups[view]:
                observed = vertices(obj); expected = data[obj['component']][frame-1]
                evaluated = obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
                matrix = np.asarray(evaluated.matrix_world)
                assert not evaluated.hide_render, (frame, obj.name)
                assert np.allclose(matrix, record['transforms'][view][frame-1], atol=1e-5)
                assert np.ptp(np.linalg.svd(matrix[:3,:3], compute_uv=False)) < 1e-5, 'Stretched protein'
                assert observed.shape == expected.shape and np.isfinite(observed).all()
                error = float(np.max(abs(observed-expected))); maximum = max(maximum, error)
                assert error < 2e-4, (frame, obj.name, error)
                projected = observed@matrix[:3,:3].T+matrix[:3,3]
                center = -5.5 if view == 'overview' else 5.5
                assert abs(projected[:,0]-center).max() < 4.4 and abs(projected[:,1]).max() < 3.2, 'Clipped molecule'
                if view == 'overview' and frame == frames: final[obj['component']] = projected[:,:2]
                checked += 1
    gaps = {}
    for first, second in [('enzyme','distal'), ('enzyme','proximal'), ('distal','proximal')]:
        gap = float(np.maximum(final[first].min(0)-final[second].max(0), final[second].min(0)-final[first].max(0)).max())
        assert gap > .1, (first, second, 'Final projected products overlap'); gaps[first+'/'+second] = gap
    bpy.context.window.view_layer = scene.view_layers['Reaction annotations']
    leader = scene.objects['Water identity leader']
    for obj in [leader, *[scene.objects['Water oxygen highlight '+view] for view in groups]]:
        for mod in obj.modifiers: mod.show_viewport = False
    for frame in range(1, frames+1):
        scene.frame_set(frame)
        for view in groups:
            marker = scene.objects['Water oxygen highlight '+view]
            assert np.max(abs(vertices(marker)-data['water'][frame-1])) < 2e-4
            matrix = np.asarray(marker.evaluated_get(bpy.context.evaluated_depsgraph_get()).matrix_world)
            assert np.allclose(matrix, record['transforms'][view][frame-1], atol=1e-5)
        matrix = np.asarray(record['transforms']['detail'][frame-1]); expected = matrix[:3,:3]@data['water'][frame-1,0]+matrix[:3,3]
        assert np.max(abs(vertices(leader)[0,:2]-expected[:2])) < 2e-4, 'Detached water identity leader'
    bpy.context.window.view_layer = scene.view_layers[scene['compositor_label_layer']]
    fonts = [obj for obj in scene.objects if obj.type == 'FONT']; captions = {row[2] for row in source['config']['captions']}
    for obj in fonts:
        assert obj.data.font.packed_file and obj.data.font_bold.packed_file
        for char, style in zip(obj.data.body, obj.data.body_format): assert style.use_bold == (ord(char) < 128 or char in '·—')
    for frame in range(1, frames+1):
        scene.frame_set(frame); graph = bpy.context.evaluated_depsgraph_get()
        active = [obj for obj in fonts if obj.data.body in captions and not obj.evaluated_get(graph).hide_render]
        assert len(active) == 1, (frame, 'Caption missing or overlaps')
    tree = scene.compositing_node_group
    mask = next(n for n in tree.nodes if n.label == 'Select matching premultiplied panel alpha')
    assert mask.inputs['Type'].default_value == 'Replace Alpha' and mask.inputs['Alpha'].is_linked
    result = {'passed':True, 'frames_checked':frames, 'actor_frame_checks':checked, 'heavy_atoms':len(identities),
              'maximum_coordinate_error_nm':maximum, 'two_views_and_water_annotation_checked':True,
              'packed_fonts_and_caption_timing_checked':True, 'final_projected_gaps':gaps,
              'source_blend':bpy.data.filepath, 'scope':source['scope'],
              'limits':'Saved coordinates, identity and display setup; rendered visibility and motion review remain separate.'}
    a.output.parent.mkdir(parents=True, exist_ok=True); a.output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__': main()
