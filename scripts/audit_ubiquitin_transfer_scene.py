"""Check reopened Blender atom motion, full-protein views and caption timing."""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_phosphotransfer_scene import vertices


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.output.exists(): p.error('Use a new audit output file')
    scene = bpy.context.scene; root = Path(bpy.data.filepath).parent
    record = json.loads((root/'manifest.json').read_text()); source = record['preparation']
    data = np.load(Path(record['asset_directory'])/'motion.npz'); frames = source['frames']
    assert source['passed'] and source['checks']['passed']
    assert scene.frame_end == frames and scene.render.fps == source['fps']
    objects = [scene.objects[name] for name in record['actors']]; groups = {'overview': [], 'detail': []}
    for obj in objects:
        groups[obj['view']].append(obj)
        assert not any(mod.type == 'BOOLEAN' for mod in obj.modifiers), obj.name
        assert json.loads(obj['atom_ids']) == source['atom_ids'][obj['component']]
        if obj['view'] == 'detail': assert obj.data == scene.objects[obj['component']+' overview'].data
        if obj['component'] != 'zinc':
            assert any(n.bl_idname == 'GeometryNodePointsToVolume' for mod in obj.modifiers if mod.type == 'NODES' for n in mod.node_group.nodes)
        for mod in obj.modifiers: mod.show_viewport = False
    identities = [i for obj in groups['overview'] for i in json.loads(obj['atom_ids'])]
    assert len(identities) == source['heavy_atoms_initial'] == source['heavy_atoms_final'] == 4476
    assert len(set(identities)) == 4476 and len(source['atom_ids']['zinc']) == 4
    assert set(obj['component'] for obj in groups['detail']) == {'donor', 'acceptor'}
    assert 'donor:76:OXT' not in identities and 'acceptor:76:OXT' in identities
    maximum = 0.; checked = 0
    for view, layer in [('overview', 'Overview'), ('detail', 'Ubiquitin detail')]:
        bpy.context.window.view_layer = scene.view_layers[layer]
        for frame in range(1, frames+1):
            scene.frame_set(frame)
            for obj in groups[view]:
                expected = data[obj['component']][frame-1]; observed = vertices(obj)
                evaluated = obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
                matrix = np.asarray(evaluated.matrix_world)
                assert np.allclose(matrix, record['transforms'][view][frame-1], atol=1e-5), (frame, obj.name, 'Framing transform')
                assert np.ptp(np.linalg.svd(matrix[:3, :3], compute_uv=False)) < 1e-5, 'Nonuniform protein scale'
                projected = observed@matrix[:3, :3].T+matrix[:3, 3]
                center_x = -5.5 if view == 'overview' else 5.5
                assert abs(projected[:, 0]-center_x).max() < 4.4 and abs(projected[:, 1]).max() < 3.2, 'Protein cropped by panel'
                assert observed.shape == expected.shape, (frame, obj.name)
                error = float(np.max(abs(observed-expected))); maximum = max(maximum, error)
                assert error < 2e-4, (frame, obj.name, error)
                assert np.isfinite(observed).all(); checked += 1
    bpy.context.window.view_layer = scene.view_layers[scene['compositor_label_layer']]
    captions = {row[2] for row in source['config']['captions']}
    fonts = [obj for obj in scene.objects if obj.type == 'FONT']
    for obj in fonts:
        assert obj.data.font.packed_file and obj.data.font_bold.packed_file
        for char, style in zip(obj.data.body, obj.data.body_format):
            assert style.use_bold == (ord(char) < 128 or char in '·—')
    for frame in range(1, frames+1):
        scene.frame_set(frame); graph = bpy.context.evaluated_depsgraph_get()
        active = [obj for obj in fonts if obj.data.body in captions and not obj.evaluated_get(graph).hide_render]
        assert len(active) == 1, (frame, 'Absent or overlapping step caption')
    windows = [n for n in scene.compositing_node_group.nodes if n.bl_idname == 'CompositorNodeBoxMask']
    assert len(windows) == 2 and all(n.inputs['Mask'].is_linked for n in windows), 'Masks require an image domain'
    assert any(n.bl_idname == 'CompositorNodeSetAlpha' for n in scene.compositing_node_group.nodes), 'Framing must mask color and alpha'
    panel_alpha = next(n for n in scene.compositing_node_group.nodes if n.label == 'Select matching premultiplied panel alpha')
    assert panel_alpha.inputs['Type'].default_value == 'Replace Alpha' and panel_alpha.inputs['Alpha'].is_linked
    result = {'passed': True, 'source_blend': bpy.data.filepath, 'frames_checked': frames,
        'actor_frame_checks': checked, 'maximum_coordinate_error_nm': maximum, 'coordinate_tolerance_nm': 2e-4,
        'overview_heavy_atoms': 4476, 'detail_heavy_atoms': len(source['atom_ids']['donor'])+len(source['atom_ids']['acceptor']),
        'four_source_zinc_preserved': True, 'source_atom_identity_and_terminal_chemistry_checked': True,
        'shared_coordinates_between_views': True, 'whole_atom_envelopes': True, 'caption_timing_checked': True,
        'packed_fonts_checked': True, 'image_domain_masks_checked': True,
        'limits': 'Checks transferred coordinates, identities and display setup; final visual/encoded-video review is separate. Product release is outside this example.'}
    a.output.parent.mkdir(parents=True, exist_ok=True); a.output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__': main()
