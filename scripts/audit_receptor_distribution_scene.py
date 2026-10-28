"""Check a reopened receptor-distribution .blend against its prepared sources.

This verifies saved geometry/transforms/labels. It does not replace pixel review
or the outside-Blender heavy-atom distance screen in the preparation step.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Vector, Euler


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project', type=Path, required=True)
    a = p.parse_args(sys.argv[sys.argv.index('--')+1:])
    manifest = json.loads((a.project/'manifest.json').read_text()); prepared = manifest['preparation']
    source_dir = Path(prepared['asset_dir']); config = prepared['config']
    distribution_dir = Path(manifest['build_parameters']['distribution_dir'])
    data = np.load(distribution_dir/'distribution.npz')
    if hashlib.sha256((distribution_dir/'distribution.npz').read_bytes()).hexdigest() != prepared['data_sha256']:
        raise ValueError('Prepared coordinates changed')
    scene = bpy.context.scene; issues = []; checked = 0; bounds = []; meshes_checked = set()
    expected_names = {row['object'] for row in manifest['objects']}
    present_names = {obj.name for obj in scene.objects if 'instance_id' in obj}
    if expected_names != present_names: issues.append('Missing or unexpected receptor instances')
    for chain in config['chains']:
        source = np.load(source_dir/f'chain_{chain}.npz')
        if hashlib.sha256((source_dir/f'chain_{chain}.npz').read_bytes()).hexdigest() != prepared['chain_hashes'][chain]:
            issues.append('Changed source chain '+chain)
        for row in [r for r in manifest['objects'] if r['chain'] == chain]:
            obj = scene.objects.get(row['object'])
            if obj is None: continue
            checked += 1
            mesh_key = (chain, obj.data.as_pointer())
            if mesh_key not in meshes_checked:
                coordinates = np.empty(len(obj.data.vertices)*3, dtype=np.float32); obj.data.vertices.foreach_get('co', coordinates)
                if coordinates.size != source['vertices'].size or not np.allclose(coordinates.reshape(-1, 3), source['vertices'], atol=1e-6, rtol=0):
                    issues.append('Source mesh changed: '+obj.name)
                faces = np.array([tuple(face.vertices) for face in obj.data.polygons])
                if not np.array_equal(faces, source['faces']): issues.append('Source topology changed: '+obj.name)
                meshes_checked.add(mesh_key)
            i = row['instance']; matrix = np.array(obj.matrix_world)
            if not np.allclose(matrix[:3, :3], data['rotations'][i], atol=2e-6, rtol=0) or not np.allclose(matrix[:3, 3], data['positions'][i], atol=1e-5, rtol=0):
                issues.append('Instance transform differs: '+obj.name)
            if obj.hide_render or [m.type for m in obj.modifiers] != ['SMOOTH']:
                issues.append('Unexpected hiding or geometry modifier: '+obj.name)
            for corner in obj.bound_box:
                point = world_to_camera_view(scene, scene.camera, obj.matrix_world@Vector(corner))
                bounds.append([point.x, point.y])
    if checked != config['count']*len(config['chains']): issues.append('Wrong complex/chain count')
    bounds = np.asarray(bounds)
    if len(bounds) and (bounds.min() < .025 or bounds.max() > .975): issues.append('Subject framing outside safe image bounds')
    lipid = scene.objects['Intact two-tail lipid instances']
    heads = np.empty(len(lipid.data.vertices)*3, dtype=np.float32); lipid.data.vertices.foreach_get('co', heads)
    if not np.array_equal(heads.reshape(-1, 3), data['lipid_heads']): issues.append('Changed lipid head positions')
    # Blender FLOAT_VECTOR foreach_get requires a float32 buffer here; float64
    # returned raw misinterpreted bytes in Blender 5.2, despite scalar reads being correct.
    angles = np.empty(len(lipid.data.vertices)*3, dtype=np.float32)
    lipid.data.attributes['lipid_rotation'].data.foreach_get('vector', angles)
    directions = np.array([Euler(angle).to_matrix()@Vector((0, 0, 1)) for angle in angles.reshape(-1, 3)])
    direction_error = np.linalg.norm(directions-data['lipid_directions'], axis=1)
    if direction_error.max() > 5e-6: issues.append('Changed leaflet direction')
    fonts = [obj for obj in scene.objects if obj.type == 'FONT']
    for obj in fonts:
        if not obj.data.font.packed_file or not obj.data.font_bold.packed_file: issues.append('Unpacked font '+obj.name)
        for corner in obj.bound_box:
            point = world_to_camera_view(scene, scene.camera, obj.matrix_world@Vector(corner))
            if not (0 <= point.x <= 1 and 0 <= point.y <= 1): issues.append('Text outside frame: '+obj.name); break
    if not scene.get('compositor_label_layer'): issues.append('Typography is not isolated from artistic filters')
    report = {'passed': not issues, 'issues': issues, 'receptor_chain_instances': checked, 'complexes': config['count'],
        'source_atoms_per_complex': prepared['atoms_per_instance'], 'lipids': len(heads)//3,
        'font_objects': len(fonts), 'compositor': scene['selected_compositor'],
        'maximum_lipid_direction_vector_error': float(direction_error.max()),
        'subject_bounds_normalized': [bounds.min(0).tolist(), bounds.max(0).tolist()] if len(bounds) else None,
        'minimum_intercomplex_atom_distance_nm': prepared['minimum_intercomplex_atom_distance_nm'],
        'scope': 'Saved source meshes, topology, rigid placements, leaflet direction, packed fonts and framing. Distance value is from preparation; final visual review is separate.'}
    (a.project/'scene-audit.json').write_text(json.dumps(report, ensure_ascii=False, indent=2))
    print(json.dumps(report))
    if issues: raise RuntimeError('Saved distribution scene audit failed')


if __name__ == '__main__': main()
