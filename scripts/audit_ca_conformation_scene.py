"""Check saved CA animation, active density/material links and representative bounds."""
import argparse,hashlib,json,sys
from pathlib import Path
import bpy,numpy as np
from mathutils import Vector
from bpy_extras.object_utils import world_to_camera_view

sys.path.insert(0,str(Path(__file__).resolve().parent))
from ca_conformation_assets import load_assets,digest
from coordinate_animation import verify_point_animation
from audit_surface_materials import check_materials
from audit_interaction_scene import displayed_materials
from detail_scene import linear


def audit(asset_dir,display_path,style_path):
    scene=bpy.context.scene;record=json.loads(scene['ca_conformation'])
    source,data=load_assets(asset_dir);display=json.loads(display_path.read_text());profile=json.loads(style_path.read_text())
    for field,path in [('source_manifest_sha256',asset_dir/'manifest.json'),('display_sha256',display_path),('style_sha256',style_path)]:
        if record[field]!=digest(path):raise ValueError('Saved input binding differs: '+field)
    if record['display']!=display or record['style']!=profile:raise ValueError('Saved display/profile record differs')
    origin=data['points'][0].mean(0)
    def equal(actual,expected,label):
        if not np.allclose(actual,expected,atol=1e-6,rtol=1e-6):raise ValueError(label)
    equal(record['origin_nm'],origin,'Source origin differs')
    shaders=check_materials(profile,set(profile['component_roles'].values()))
    bodies=[];probes=[]
    wanted=['NodeGroupInput','GeometryNodeMeshToPoints','GeometryNodePointsToVolume','GeometryNodeVolumeToMesh',
            'GeometryNodeSetMaterial','GeometryNodeSetShadeSmooth','NodeGroupOutput']
    for chain in source['chains']:
        obj=scene.objects['Chain '+chain];probe=scene.objects['CA probe '+chain];ids=np.flatnonzero(data['keys'][:,0]==chain)
        if obj.data!=probe.data or probe.modifiers or not probe.hide_render or obj.hide_render or len(obj.modifiers)!=1:
            raise ValueError('Displayed body no longer shares its retained coordinate probe')
        equal(obj.matrix_world,np.eye(4),'Source body transform differs');equal(probe.matrix_world,np.eye(4),'Probe transform differs')
        if json.loads(obj['source_indices'])!=ids.tolist() or json.loads(obj['source_ids'])!=[source['point_ids'][i] for i in ids]:
            raise ValueError('Saved source identity differs')
        verify_point_animation(obj,data['points'][:,ids]-origin,data['frames'])
        mod=obj.modifiers[0]
        if mod.type!='NODES' or not mod.show_render or not mod.show_viewport:raise ValueError('Density display disabled')
        tree=mod.node_group
        if sorted(n.bl_idname for n in tree.nodes)!=sorted(wanted):raise ValueError('Unreviewed density modifier nodes')
        nodes={n.bl_idname:n for n in tree.nodes}
        for first,second in zip(wanted,wanted[1:]):
            if not any(link.from_node==nodes[first] and link.to_node==nodes[second] for link in tree.links):
                raise ValueError('Disconnected source-to-surface geometry')
        volume=nodes['GeometryNodePointsToVolume'];iso=nodes['GeometryNodeVolumeToMesh']
        equal(volume.inputs['Radius'].default_value,display['density']['volume_radius_nm'],'Density radius differs')
        equal(volume.inputs['Voxel Size'].default_value,display['density']['volume_voxel_nm'],'Density voxel size differs')
        equal(iso.inputs['Threshold'].default_value,.35,'Density threshold differs');equal(iso.inputs['Adaptivity'].default_value,.08,'Surface adaptation differs')
        bodies.append(obj);probes.append((probe,ids));mod.show_viewport=False
    expected_components={**{'Chain '+c:c for c in source['chains']},**({'Membrane heads':'heads','Membrane tails':'tails'} if display.get('membrane') else {})}
    for name,component in expected_components.items():
        obj=scene.objects[name];material=bpy.data.materials[name]
        if displayed_materials(obj)!={material} or material['surface_role']!=profile['component_roles'][component]:
            raise ValueError('Actual role material assignment differs: '+name)
        color=(*linear(profile['semantic_colors'][component]),1);equal(material.diffuse_color,color,'Display color differs')
        socket=material.node_tree.nodes['Principled BSDF'].inputs['Base Color']
        if socket.is_linked:raise ValueError('Unreviewed color override')
        equal(socket.default_value,color,'Shader color differs')
    for name,expected in record['environment_mesh_sha256'].items():
        obj=scene.objects[name];equal(obj.matrix_world,np.eye(4),'Membrane transform differs')
        if obj.modifiers or obj.hide_render:raise ValueError('Membrane geometry display differs')
        xyz=np.empty(len(obj.data.vertices)*3,np.float32);obj.data.vertices.foreach_get('co',xyz)
        if hashlib.sha256(xyz.tobytes()).hexdigest()!=expected:raise ValueError('Saved membrane geometry changed')
    timing=source['config']['presentation']
    if scene.frame_start!=1 or scene.frame_end!=timing['frames']+1 or scene.render.fps!=timing['fps']:
        raise ValueError('Saved finite presentation differs')
    frames=data['frames'];samples=np.unique(np.r_[frames,(frames[:-1]+frames[1:])/2,np.arange(1,scene.frame_end+1)])
    maximum=0.
    for frame in samples:
        scene.frame_set(int(frame),subframe=float(frame-int(frame)));graph=bpy.context.evaluated_depsgraph_get()
        actual=np.empty_like(data['points'][0])
        for probe,ids in probes:
            evaluated=probe.evaluated_get(graph);mesh=evaluated.to_mesh();xyz=np.empty(len(mesh.vertices)*3)
            mesh.vertices.foreach_get('co',xyz);evaluated.to_mesh_clear();actual[ids]=xyz.reshape(-1,3)+origin
        index=min(len(frames)-2,max(0,int(np.searchsorted(frames,frame,side='right')-1)))
        fraction=np.clip((frame-frames[index])/(frames[index+1]-frames[index]),0,1)
        expected=data['points'][index]*(1-fraction)+data['points'][index+1]*fraction
        maximum=max(maximum,float(np.linalg.norm(actual-expected,axis=1).max()))
        if maximum>3e-5:raise ValueError('Saved evaluated coordinates differ from the checked linear path')
    for obj in bodies:obj.modifiers[0].show_viewport=True
    representative=[1,1+timing['frames']//2,1+timing['frames']];bounds={}
    for frame in representative:
        scene.frame_set(frame);bpy.context.view_layer.update();graph=bpy.context.evaluated_depsgraph_get();bounds[str(frame)]={}
        for obj in bodies+[scene.objects[n] for n in record['environment_mesh_sha256']]:
            evaluated=obj.evaluated_get(graph)
            coordinates=np.array([world_to_camera_view(scene,scene.camera,evaluated.matrix_world@Vector(c))[:2] for c in evaluated.bound_box])
            if np.any(coordinates<.015) or np.any(coordinates>.985):raise ValueError('Representative complete surface cropped: '+obj.name)
            bounds[str(frame)][obj.name]=[coordinates.min(0).tolist(),coordinates.max(0).tolist()]
    camera=scene.camera
    if camera.data.type!='ORTHO' or camera.data.dof.use_dof!=(display.get('fstop') is not None):raise ValueError('Camera mode differs')
    if camera.data.dof.use_dof:equal(camera.data.dof.aperture_fstop,display['fstop'],'Depth of field differs')
    equal(scene.world.node_tree.nodes['Background'].inputs['Color'].default_value,(*linear(profile['world']),1),'World color differs')
    if scene.view_settings.view_transform!=profile['view_transform']:raise ValueError('View transform differs')
    equal(scene.view_settings.exposure,profile['exposure'],'Exposure differs')
    for item in profile['lighting']:
        light=scene.objects[item['name']];equal(light.location,item['position'],'Light position differs')
        equal(light.data.energy,item['power'],'Light power differs');equal(light.data.size,item['size'],'Light size differs')
        equal(light.data.color,linear(item['color']),'Light color differs')
    labels=[obj for obj in scene.objects if obj.type=='FONT']
    if {obj.data.body for obj in labels}!={display['title'],display['subtitle']}:raise ValueError('Main labels differ')
    if scene.get('compositor_label_layer') not in scene.view_layers:raise ValueError('Protected typography layer missing')
    for obj in labels:
        if obj.data.font.packed_file is None or obj.data.font_bold.packed_file is None:raise ValueError('Pack both supplied fonts')
        if any(bool(fmt.use_bold)!=(ord(char)<128 or 0x370<=ord(char)<=0x3ff) for char,fmt in zip(obj.data.body,obj.data.body_format)):
            raise ValueError('Chinese and English/Greek font assignment differs')
        for corner in obj.bound_box:
            xy=world_to_camera_view(scene,camera,obj.matrix_world@Vector(corner))
            if not (0<=xy.x<=1 and 0<=xy.y<=1):raise ValueError('Main label cropped')
    return {'passed':True,'points':len(source['point_ids']),'stored_samples':len(frames),'evaluated_times':len(samples),
            'maximum_coordinate_error_nm':maximum,'materials':len(shaders['materials']),'representative_surface_bounds':bounds,
            'limits':'Stored linear path and evaluated probes, active density/material links, selected surface bounds; not atomic dynamics or continuous membrane-envelope clearance'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('asset-dir','display','style-profile','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);result=audit(a.asset_dir,a.display,a.style_profile)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='representative_surface_bounds'}))


if __name__=='__main__':main()
