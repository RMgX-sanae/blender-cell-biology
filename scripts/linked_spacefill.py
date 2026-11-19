"""Substitute source heavy-atom centers for one registered chain in a native detail view."""
import hashlib
import json
import numpy as np
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Vector

from detail_scene import SPACEFILL_RADII_NM
from molecular_asset_io import chain_file
from style_nodes import instance_spheres, point_mesh


def check_surface(obj, data):
    mesh = obj.data
    if mesh.shape_keys or any(m.show_render and m.type != 'SMOOTH' for m in obj.modifiers):
        raise ValueError('The visible surface must use the source coordinates with at most display smoothing')
    vertices = np.empty((len(mesh.vertices), 3), dtype=np.float32)
    mesh.vertices.foreach_get('co', vertices.ravel())
    indices = np.empty(len(mesh.loops), dtype=np.int32)
    mesh.loops.foreach_get('vertex_index', indices)
    faces = data['faces']
    if (not np.array_equal(vertices, data['vertices']) or faces.shape != (len(mesh.polygons), 3)
            or not np.array_equal(indices, faces.ravel())):
        raise ValueError('Prepared vertices or triangles do not match the saved surface of ' + obj.name)
    basis = np.asarray(obj.matrix_world.to_3x3())
    gram = basis.T @ basis
    scale_sq = np.trace(gram) / 3
    if scale_sq <= 0 or np.linalg.det(basis) <= 0 or not np.allclose(gram, scale_sq*np.eye(3), atol=2e-5*scale_sq, rtol=0):
        raise ValueError('The focused chain must have a proper rotation and uniform scene scale')
    return float(np.sqrt(scale_sq))


def switch_to_spacefill(scene, camera, focus, spec, brief_dir):
    obj = scene.objects[focus['object']]
    chain = obj.get('source_chain')
    materials = [slot.material for slot in obj.material_slots]
    if not chain or len(materials) != 1 or materials[0] is None:
        raise ValueError('Space-fill requires a source chain and its single identity material')
    material = materials[0]
    directory = (brief_dir / spec['asset_dir']).resolve()
    raw_manifest = (directory/'manifest.json').read_bytes()
    manifest_hash = hashlib.sha256(raw_manifest).hexdigest()
    if manifest_hash != spec['manifest_sha256']:
        raise ValueError('Molecular asset manifest differs from the reference in the scene brief')
    manifest = json.loads(raw_manifest)
    if manifest.get('units') != 'nm' or not manifest.get('source_sha256'):
        raise ValueError('The asset manifest must identify its source and give coordinates in nm')
    path = chain_file(directory, chain)
    with path.open('rb') as stream:
        chain_hash = hashlib.file_digest(stream,'sha256').hexdigest()
    expected = spec.get('chain_sha256')
    if (manifest.get('chain_file_schema') != 1 or expected is not None) and expected != chain_hash:
        raise ValueError('The chain file must match its qualified fingerprint')

    with np.load(path, allow_pickle=False) as data:
        scale = check_surface(obj, data)
        atoms = data['atoms']
        fields = {k:data[k] for k in ('elements','atom_names','residue_numbers','residue_names')}
        if (atoms.ndim != 2 or atoms.shape[1] != 3 or not len(atoms) or not np.isfinite(atoms).all()
                or any(values.shape != (len(atoms),) for values in fields.values())):
            raise ValueError('The molecular asset has incomplete atom coordinates or identities')
        if 'insertion_codes' in data and data['insertion_codes'].shape != (len(atoms),):
            raise ValueError('Insertion codes must correspond to the prepared atom rows')
        elements = fields['elements']
        element_radii = {str(e):SPACEFILL_RADII_NM.get(str(e), .17) for e in np.unique(elements)}
        radii = np.array([element_radii[e] for e in elements], dtype=np.float32)
        local = np.array(focus['local'])
        distances = np.linalg.norm(atoms-local, axis=1)
        row = int(distances.argmin())
        center = obj.matrix_world @ Vector(atoms[row])
        screen = world_to_camera_view(scene, camera, center)
        landmark = {'row':row, 'element':str(elements[row]), 'atom':str(fields['atom_names'][row]),
                    'residue_name':str(fields['residue_names'][row]),
                    'residue_number':int(fields['residue_numbers'][row]),
                    'insertion_code':str(data['insertion_codes'][row]) if 'insertion_codes' in data else None,
                    'local_nm':atoms[row].tolist(), 'world':list(center), 'detail_at':[screen.x,1-screen.y],
                    'distance_to_surface_pick_nm':float(distances[row])}
        # Give only this object a new datablock. Other copies may share the original mesh.
        mesh = point_mesh(obj.name+' source atoms', atoms, radii)
        mesh.materials.append(material)
        for name, values in [('source_atom_row',np.arange(len(atoms),dtype=np.int32)),
                             ('source_residue_number',fields['residue_numbers'].astype(np.int32))]:
            attribute = mesh.attributes.new(name, 'INT', 'POINT')
            attribute.data.foreach_set('value', values)
        obj.data = mesh
        obj.modifiers.clear()
        instance_spheres(obj, material, radius=.17, subdivisions=2)
        obj['evidence'] = 'Source heavy-atom centers; space-filling radii are display conventions.'
        obj['source_atom_file_sha256'] = chain_hash

    return {'kind':'source_spacefill','object':obj.name,'source_chain':chain,'source_instance':obj.get('instance_id'),
            'source_sha256':manifest['source_sha256'],'entry_url':manifest.get('entry_url'),
            'model_index':manifest.get('model_index'),'asset_manifest_sha256':manifest_hash,
            'chain_file':path.name,'chain_sha256':chain_hash,'atom_count':len(atoms),'material':material.name,
            'scene_units_per_nm':scale,'radii_nm':element_radii,
            'elements_with_illustrative_radius':sorted(element_radii.keys()-SPACEFILL_RADII_NM.keys()),
            'nearest_prepared_atom':landmark}
