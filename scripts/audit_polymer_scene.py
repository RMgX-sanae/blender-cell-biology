"""Evaluate saved polymer templates: rigid subunits, identity, motion and endpoints.

Run in Blender with the saved scene loaded. This is geometric validation, not
an assessment of biological kinetics or a substitute for inspecting the video.
"""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Vector


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--poses', type=Path, help='Optional actual world transforms for external contact checks')
    args = p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    folder = Path(bpy.data.filepath).parent
    manifest = json.loads((folder/'manifest.json').read_text())
    scene = bpy.context.scene
    units = [o for o in scene.objects if o.get('role') == 'polymer_unit']
    errors = []
    if sorted(o.name for o in units) != sorted(manifest['unit_ids']):
        errors.append('Stable unit identity/count mismatch')
    if len({o['polymer_unit_id'] for o in units}) != len(units):
        errors.append('Duplicate unit IDs')
    previous, maximum_step, maximum_scale_error, endpoint_error = {}, 0., 0., 0.
    branches = [g.name for g in bpy.data.node_groups if g.name.startswith('COMP_')]
    if any('COMP_'+m not in branches for m in ('clean','reference','monochrome')):
        errors.append('Saved compositor alternatives missing')
    baseline_sizes = {o.name: len(o.data.vertices) for o in units}
    tracks = []
    for frame in range(scene.frame_start, scene.frame_end+1):
        scene.frame_set(frame)
        bpy.context.view_layer.update()
        if args.poses:
            tracks.append([np.array(o.matrix_world) for o in units])
        for obj in units:
            matrix = np.array(obj.matrix_world)
            if not np.isfinite(matrix).all():
                errors.append(f'{obj.name} non-finite transform at {frame}')
                continue
            # Singular values detect shears and reflected/nonuniform scale.
            scale_error = float(np.max(np.abs(np.linalg.svd(matrix[:3,:3], compute_uv=False)-1)))
            maximum_scale_error = max(maximum_scale_error, scale_error)
            if scale_error > 2e-5 or np.linalg.det(matrix[:3,:3]) < 0:
                errors.append(f'{obj.name} changed molecular scale/handedness at {frame}')
            if len(obj.data.vertices) != baseline_sizes[obj.name] or any(m.type=='BOOLEAN' for m in obj.modifiers):
                errors.append(f'{obj.name} no longer an intact source mesh')
            if obj.hide_render:
                errors.append(f'{obj.name} hidden instead of continuous motion at {frame}')
            if obj.name in previous:
                maximum_step = max(maximum_step, float(np.linalg.norm(matrix[:3,3]-previous[obj.name])))
            previous[obj.name] = matrix[:3,3].copy()
        for event in manifest['events']:
            if frame in (event['first_frame'], event['last_frame']):
                key = 'start_local_matrix' if frame==event['first_frame'] else 'end_local_matrix'
                expected = np.array(event[key])
                actual = np.array(bpy.data.objects[event['object']].matrix_basis)
                endpoint_error = max(endpoint_error, float(np.max(np.abs(expected-actual))))
        for obj in scene.objects:
            if obj.get('role') == 'annotation':
                corners = [world_to_camera_view(scene, scene.camera, obj.matrix_world @ Vector(p)) for p in obj.bound_box]
                if any(p.x < .01 or p.x > .99 or p.y < .01 or p.y > .99 for p in corners):
                    errors.append(f'Annotation outside safe frame: {obj.name} at {frame}')
    if endpoint_error > 2e-4:
        errors.append('Animation endpoint mismatch')
    # Detect severe jumps, not a claimed physical speed limit.
    if maximum_step > 2.5:
        errors.append('Frame-to-frame origin displacement exceeds 2.5 nm; inspect timing')
    final_error = max(float(np.max(np.abs(np.array(o.matrix_world)-np.array(manifest['endpoints_world'][o.name])))) for o in units)
    if final_error > 2e-4:
        errors.append('Saved final assembly differs from manifest')
    report = {'source_blend': bpy.data.filepath, 'frames_evaluated': scene.frame_end-scene.frame_start+1,
              'units': len(units), 'distinct_meshes': len({o.data.name for o in units}),
              'maximum_scale_error': maximum_scale_error, 'maximum_origin_step_nm': maximum_step,
              'endpoint_matrix_max_error': endpoint_error, 'final_matrix_max_error': final_error,
              'compositors': branches, 'failures': sorted(set(errors)),
              'limits': ['No molecular force/energy or chemical kinetics tested.',
                         'Does not prove collision-free paths or visual readability; inspect actual rendered frames.']}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    if args.poses:
        args.poses.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(args.poses, matrices=np.array(tracks, dtype=np.float32),
                            ids=np.array([o.name for o in units]),
                            kinds=np.array([o['asset_kind'] for o in units]),
                            chains=np.array([o['source_chain'] for o in units]))
    print(json.dumps(report, ensure_ascii=False))
    if errors:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
