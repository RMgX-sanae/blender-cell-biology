"""Check actual saved interaction geometry, transforms, framing and loop seam."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Vector

sys.path.insert(0,str(Path(__file__).resolve().parent))
from interaction_assets import verify_assets
from interaction_membrane import audit_display,plane_clearance
from detail_scene import linear


def displayed_materials(obj):
    """Read effective surface slots or reachable atom-instance material setters."""
    if not obj.modifiers:
        slots=obj.material_slots
        return {slots[p.material_index].material if p.material_index<len(slots) else None for p in obj.data.polygons}
    assigned=set()
    for modifier in obj.modifiers:
        if modifier.type!='NODES' or not modifier.show_render or modifier.node_group is None:return {None}
        pending=[n for n in modifier.node_group.nodes if n.type=='GROUP_OUTPUT' and n.is_active_output]
        reached=set()
        while pending:
            node=pending.pop()
            if node in reached:continue
            reached.add(node)
            pending.extend(link.from_node for socket in node.inputs for link in socket.links)
        for node in reached:
            if node.type=='SET_MATERIAL':
                socket=node.inputs['Material']
                assigned.add(None if socket.is_linked else socket.default_value)
    return assigned


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    scene=bpy.context.scene;folder=Path(scene['interaction_asset_dir'])
    manifest_path=folder/'manifest.json'
    manifest=json.loads(manifest_path.read_text())
    fingerprint=scene.get('interaction_manifest_sha256')
    if fingerprint and hashlib.sha256(manifest_path.read_bytes()).hexdigest()!=fingerprint:
        raise ValueError('Prepared interaction manifest changed since this scene was built')
    integrity=verify_assets(folder,manifest) if fingerprint or 'data_sha256' in manifest else {
        'verified':False,'reason':'Legacy preparation has no file hashes; only saved geometry and poses are compared'}
    integrity['scene_manifest_bound']=bool(fingerprint)
    expected=np.load(folder/'motion.npz')['matrices']
    snapshot=scene.get('interaction_snapshot_frame')
    if snapshot is not None:
        if type(snapshot) is not int or not 1<=snapshot<=len(expected):raise ValueError('Invalid snapshot source frame')
        expected=expected[snapshot-1:snapshot]
    visible=set(json.loads(scene.get('interaction_visible_actors',json.dumps([actor['id'] for actor in manifest['actors']]))))
    if visible-set(actor['id'] for actor in manifest['actors']) or any(not any(actor['id'] in visible and actor['role']==role for actor in manifest['actors']) for role in ['moving','anchor']):
        raise ValueError('Require known visible actors from both partners')
    errors=[];meshes=[];actor_data={};membrane_report=None;membrane_minimum=None
    colors=json.loads(scene.get('interaction_resolved_colors','{}'))
    if colors and set(colors)!={actor['id'] for actor in manifest['actors']}:errors.append('Resolved actor color inventory differs')
    if manifest.get('binding_candidate'):
        candidate=manifest['binding_candidate']
        if scene.get('binding_pose_kind')!=candidate['pose_kind'] or scene.get('binding_evidence_class')!='inferred':
            errors.append('Binding hypothesis evidence changed')
        evidence='Static review of one prepared pose; no motion validation' if snapshot is not None else manifest['animation_evidence']
        if scene.get('animation_evidence')!=evidence:errors.append('Candidate motion evidence changed')
        titles=[obj for obj in scene.objects if obj.type=='FONT' and obj.data.body==manifest['display_title'] and not obj.hide_render]
        if not titles:errors.append('Readable binding-hypothesis label missing')
    if manifest['example'].get('requires_membrane') and not (manifest.get('membrane') and manifest.get('membrane_screen_passed')):
        raise ValueError('Required membrane preparation is missing')
    if manifest.get('membrane'):membrane_report=audit_display(scene,manifest['membrane'])
    for actor in manifest['actors']:
        obj=scene.objects.get(actor['id'])
        if obj is None:errors.append('Missing actor '+actor['id']);continue
        data=np.load(folder/actor['file'])
        actor_data[actor['id']]={key:data[key] for key in ['atoms','vertices','atom_names','residue_numbers'] if key in data}
        reference=data['atoms'] if obj.modifiers else data['vertices']
        actual=np.empty(len(obj.data.vertices)*3);obj.data.vertices.foreach_get('co',actual);actual=actual.reshape(-1,3)
        intact=actual.shape==reference.shape and np.allclose(actual,reference,atol=2e-6)
        if not obj.modifiers and 'faces' in data:
            faces=np.array([list(face.vertices) for face in obj.data.polygons])
            intact=bool(intact and np.array_equal(faces,data['faces']))
        meshes.append({'id':actor['id'],'intact_source_geometry':bool(intact),'source_chain':actor['chain']})
        if not intact:errors.append('Geometry changed '+actor['id'])
        if obj.get('biological_id')!=actor['id']:errors.append('Identity changed '+actor['id'])
        if obj.get('source_structure')!=actor.get('source_pdb',manifest['example']['pdb']):errors.append('Source identity changed '+actor['id'])
        if 'instance_id' in actor and obj.get('source_author_chain')!=actor['chain']:errors.append('Author chain changed '+actor['id'])
        if colors:
            mat=bpy.data.materials.get(actor['id']);expected_color=(*linear(colors[actor['id']]),1.)
            if displayed_materials(obj)!={mat}:errors.append('Actor material assignment changed '+actor['id'])
            if mat is None or not np.allclose(mat.diffuse_color,expected_color,rtol=0,atol=1e-6):
                errors.append('Actor display color changed '+actor['id'])
            elif mat.use_nodes:
                shader=mat.node_tree.nodes.get('Principled BSDF')
                if shader and (shader.inputs['Base Color'].is_linked or not np.allclose(shader.inputs['Base Color'].default_value,expected_color,rtol=0,atol=1e-6)):
                    errors.append('Actor surface color changed '+actor['id'])
        for key in ('instance_id','label_chain','operator','assembly_id','owner_instance'):
            if key in actor and obj.get('source_'+key)!=actor[key]:errors.append('Source instance changed '+actor['id']+' '+key)
        if obj.hide_render!=(actor['id'] not in visible):errors.append('Whole-actor visibility changed '+actor['id'])
        if obj.parent!=scene.objects.get(actor['role']):errors.append('Rigid partner ownership changed '+actor['id'])
        if any(m.type=='BOOLEAN' for m in obj.modifiers):errors.append('Cutting modifier on '+actor['id'])
    if manifest.get('membrane'):
        atom=manifest['membrane']['cargo_tether_atom'];data=actor_data[atom['actor']]
        index=np.flatnonzero((data['atom_names']==atom['atom_name'])&(data['residue_numbers']==atom['residue_number']))
        if len(index)!=1 or not np.allclose(data['atoms'][index[0]],manifest['membrane']['cargo_tether_nm'][0],atol=1e-9,rtol=0):
            raise ValueError('Cargo tether lost its source peptide endpoint')
    maximum_error=0.;minimum_margin=1.;clipped=[];actor_errors={a['id']:0. for a in manifest['actors']}
    for frame,pose in enumerate(expected,1):
        scene.frame_set(frame);bpy.context.view_layer.update()
        actual=np.array(scene.objects['moving'].matrix_world)
        maximum_error=max(maximum_error,float(abs(actual-pose).max()))
        for actor in manifest['actors']:
            obj=scene.objects[actor['id']]
            target=pose if actor['role']=='moving' else np.eye(4)
            actor_errors[actor['id']]=max(actor_errors[actor['id']],float(abs(np.array(obj.matrix_world)-target).max()))
            if manifest.get('membrane'):
                matrix=np.asarray(obj.matrix_world)[None]
                data=actor_data[actor['id']]
                gap=min(float(plane_clearance(data['atoms'],matrix,manifest['membrane'],manifest['membrane']['atom_guard_nm']).min()),
                        float(plane_clearance(data.get('vertices',data['atoms']),matrix,manifest['membrane']).min()))
                membrane_minimum=gap if membrane_minimum is None else min(membrane_minimum,gap)
            for point in (obj.bound_box if actor['id'] in visible else []):
                v=world_to_camera_view(scene,scene.camera,obj.matrix_world@Vector(point))
                margin=min(v.x,v.y,1-v.x,1-v.y);minimum_margin=min(minimum_margin,margin)
                if margin<0 or v.z<=0:clipped.append({'frame':frame,'actor':actor['id']});break
    if maximum_error>2e-5:errors.append('Saved animation differs from prepared poses')
    for name,error in actor_errors.items():
        if error>2e-5:errors.append('Actor no longer follows its rigid partner '+name)
    if clipped:errors.append('Main actor bounding boxes leave the camera frame')
    if membrane_minimum is not None and membrane_minimum<0:errors.append('Saved actor crosses the membrane core')
    if not manifest['motion']['collision_screen_passed']:errors.append('Prepared path contact screen failed')
    if scene.frame_end!=(1 if snapshot is not None else manifest['frames']) or scene.render.fps!=manifest['fps']:errors.append('Timing differs from manifest')
    result={'scene':bpy.data.filepath,'actors':meshes,'pose_error_max':maximum_error,
            'prepared_file_integrity':integrity,
            'visible_actors':sorted(visible),'omitted_whole_actors':sorted(set(actor['id'] for actor in manifest['actors'])-visible),
            'membrane_display':membrane_report,'membrane_minimum_clearance_nm':membrane_minimum,
            'actor_transform_error_max':actor_errors,
            'camera_minimum_fractional_margin':minimum_margin,'clipped_examples':clipped[:15],
            'checked_pose_count':len(expected),'virtual_seam_checked':snapshot is None,
            'checked_frames_including_virtual_seam':len(expected) if snapshot is None else None,
            'prepared_contact_screen':manifest['motion'],
            'snapshot_source_frame':snapshot,
            'errors':errors,'passed':not errors,
            'limits':'Technical and geometric checks only. Inspect readability, source suitability and biological scope separately; inspect timing when delivering animation. A snapshot checks one saved pose only.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2))
    print(json.dumps({'passed':not errors,'errors':errors,'minimum_margin':minimum_margin}))
    if errors:raise SystemExit(1)


if __name__=='__main__':main()
