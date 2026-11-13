"""Build or sample an editable whole-triskelion patch motion, without auto-run code."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import bpy
import numpy as np
from mathutils import Matrix, Quaternion, Vector, kdtree

sys.path.insert(0,str(Path(__file__).resolve().parent))
from molecular_asset_io import chain_file
import blender_clathrin_assembly as assembly
import detail_scene as ds
from compositor_stack import build_stack


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(folder):
    record=json.loads((folder/'manifest.json').read_text())
    if record['kind']!='clathrin-articulated-patch' or not record['checks']['passed']:
        raise ValueError('Require qualified articulated source assets')
    for name,expected in record['data_sha256'].items():
        if digest(folder/name)!=expected:raise ValueError('Changed motion file: '+name)
    assets=Path(record['assembly_directory']);source,_=assembly.load_assets(assets)
    if digest(assets/'assembly.json')!=record['assembly_sha256']:raise ValueError('Changed assembly')
    return record,source,assets,dict(np.load(folder/'motion.npz'))


def attribute(mesh,numbers):
    item=mesh.attributes.get('source_residue') or mesh.attributes.new('source_residue','INT','POINT')
    item.data.foreach_set('value',np.asarray(numbers,dtype=np.int32))


def hinge_group(chain,joints,angles):
    tree=bpy.data.node_groups.new('Whole-chain hinges '+chain,'GeometryNodeTree')
    for direction in ('INPUT','OUTPUT'):tree.interface.new_socket(name='Geometry',in_out=direction,socket_type='NodeSocketGeometry')
    src=tree.nodes.new('NodeGroupInput');target=tree.nodes.new('NodeGroupOutput')
    number=tree.nodes.new('GeometryNodeInputNamedAttribute');number.data_type='INT';number.inputs['Name'].default_value='source_residue'
    position=tree.nodes.new('GeometryNodeInputPosition').outputs['Position']
    for index,joint in joints:
        turn=tree.nodes.new('ShaderNodeVectorRotate');turn.rotation_type='EULER_XYZ'
        turn.name=f'Pivot {joint["pivot_residue"]}';turn.inputs['Center'].default_value=joint['pivot_nm']
        tree.links.new(position,turn.inputs['Vector'])
        select=tree.nodes.new('ShaderNodeMath');select.operation='LESS_THAN';select.inputs[1].default_value=joint['pivot_residue']
        tree.links.new(number.outputs['Attribute'],select.inputs[0])
        switch=tree.nodes.new('GeometryNodeSwitch');switch.input_type='VECTOR'
        tree.links.new(select.outputs[0],switch.inputs['Switch']);tree.links.new(position,switch.inputs['False'])
        tree.links.new(turn.outputs['Vector'],switch.inputs['True']);position=switch.outputs['Output']
        for frame,vector in enumerate(angles[:,index*3:index*3+3],1):
            angle=float(np.linalg.norm(vector));axis=Vector(vector/angle) if angle>1e-12 else Vector((1,0,0))
            turn.inputs['Rotation'].default_value=Quaternion(axis,angle).to_euler('XYZ')
            turn.inputs['Rotation'].keyframe_insert('default_value',frame=frame)
    setter=tree.nodes.new('GeometryNodeSetPosition');tree.links.new(src.outputs['Geometry'],setter.inputs['Geometry'])
    tree.links.new(position,setter.inputs['Position']);tree.links.new(setter.outputs['Geometry'],target.inputs['Geometry'])
    return tree


def build(args):
    record,source,assets,data=load(args.motion)
    if not args.font_zh.is_file() or not args.font_en.is_file():raise ValueError('Supply Songti and Times New Roman fonts')
    moving=record['config']['moving_unit'];visible=[moving]+record['config']['fixed_units']
    bpy.context.preferences.filepaths.save_version=0
    assembly.build(SimpleNamespace(assets=assets,output=args.output,show=visible,highlight=moving,
        width=args.width,samples=args.samples,engine=args.engine,style_profile=args.style_profile,compositor='clean',save=False))
    scene=bpy.context.scene;scene.frame_end=record['frames']+1;scene.render.fps=record['fps']
    scene.render.resolution_y=round(args.width*9/16)
    group=ds.collection('Source CA probes');units={u['id']:u for u in source['instances']}
    for unit_id in visible:
        for chain in units[unit_id]['chains']:
            raw=dict(np.load(chain_file(assets/'prototypes', chain)));surface=scene.objects[unit_id+'/'+chain]
            mesh=bpy.data.meshes.new('CA '+unit_id+'/'+chain);mesh.from_pydata(raw['atoms'].tolist(),[],[]);mesh.update()
            probe=bpy.data.objects.new('CA '+unit_id+'/'+chain,mesh);group.objects.link(probe)
            probe.parent=surface.parent;probe.hide_render=True;probe['source_chain']=chain;probe['assembly_unit_id']=unit_id
            attribute(mesh,raw['residue_numbers'])
            joints=[(i,j) for i,j in enumerate(record['joints']) if j['chain']==chain]
            if unit_id!=moving or not joints:continue
            search=kdtree.KDTree(len(raw['atoms']))
            for i,p in enumerate(raw['atoms']):search.insert(p,i)
            search.balance()
            nearest=[search.find(v.co)[1] for v in surface.data.vertices]
            attribute(surface.data,raw['residue_numbers'][nearest])
            tree=hinge_group(chain,joints,data['angles'])
            for obj in (surface,probe):obj.modifiers.new('Source-preserving hinges','NODES').node_group=tree
    origin=np.array(record['origin_nm']);rig=scene.objects[moving];rig.rotation_mode='QUATERNION';previous=None
    for frame,(center,rotation) in enumerate(zip(data['centers'],data['rotations']),1):
        q=Matrix(rotation).to_quaternion()
        if previous is not None and previous.dot(q)<0:q.negate()
        previous=q.copy();rig.rotation_quaternion=q;rig.location=center-rotation@origin
        rig.keyframe_insert('rotation_quaternion',frame=frame);rig.keyframe_insert('location',frame=frame)
    for action in bpy.data.actions:
        for layer in action.layers:
            for strip in layer.strips:
                for bag in strip.channelbags:
                    for curve in bag.fcurves:
                        for key in curve.keyframe_points:key.interpolation='LINEAR'
    # Fit every prepared pose; an oblique angle makes radial travel readable.
    points=np.concatenate([data['moving_ca'].reshape(-1,3),data['obstacle_ca']])
    target=points.mean(0);front=np.array(units[moving]['hub_nm']);front/=np.linalg.norm(front)
    side=np.cross(front,[0,0,1]);side/=np.linalg.norm(side)
    direction=.78*front+.63*side;direction/=np.linalg.norm(direction)
    camera=scene.camera;camera.location=target+direction*180;ds.aim(camera,target);camera.data.dof.use_dof=False
    bpy.context.view_layer.update();basis=np.array(camera.rotation_euler.to_matrix())
    projected=(points-np.array(camera.location))@basis;low=projected.min(0);high=projected.max(0)
    camera.location+=Vector(basis@np.array([(low[0]+high[0])/2,(low[1]+high[1])/2,0]))
    camera.data.ortho_scale=max(high[0]-low[0]+3,(high[1]-low[1]+3)*16/9)/.72
    labels=ds.collection('Figure labels');zh=bpy.data.fonts.load(str(args.font_zh));en=bpy.data.fonts.load(str(args.font_en));zh.pack();en.pack()
    def label(body,x,y,pixels):
        text=bpy.data.curves.new(body,'FONT');text.body=body;text.font=zh;text.font_bold=en
        for char,formatting in zip(body,text.body_format):formatting.use_bold=ord(char)<128
        text.size=camera.data.ortho_scale*pixels/args.width
        obj=bpy.data.objects.new(body,text);labels.objects.link(obj);obj.parent=camera
        obj.location=(x*camera.data.ortho_scale,y*camera.data.ortho_scale*9/16,-20);return obj
    label('网格蛋白装配示意',-.45,.40,34*args.width/1280)
    label('金色：运动单体    蓝色：相邻单体',-.45,-.43,27*args.width/1280)
    recipe=json.loads(args.recipe.read_text()) if args.recipe else {'id':'clathrin-clean','filters':[]}
    recipe.update(background=json.loads(scene['style_profile'])['background'],label_min_pixels=27*args.width/1280)
    report=build_stack(scene,recipe)
    scene['clathrin_motion_directory']=str(args.motion.resolve());scene['clathrin_motion_sha256']=digest(args.motion/'manifest.json')
    scene['clathrin_motion_scope']=record['scope'];scene['clathrin_recipe']=json.dumps(report)
    scene['clathrin_surface_mapping']='Nearest source CA residue, ordered rigid-prefix rotations; coarse surface interpolation at joints'
    scene.frame_set(record['frames']//2+1);args.output.parent.mkdir(parents=True,exist_ok=True)
    bpy.ops.wm.save_as_mainfile(filepath=str(args.output))
    print(json.dumps({'saved':str(args.output),'frames':scene.frame_end,'visible_whole_units':visible}))


def export(output):
    """Sample actual saved GN and parent animation; external audit checks the result."""
    if output.exists():raise ValueError('Use a fresh sample directory')
    scene=bpy.context.scene;folder=Path(scene['clathrin_motion_directory']);record,source,assets,data=load(folder)
    if digest(folder/'manifest.json')!=scene['clathrin_motion_sha256']:raise ValueError('Motion manifest changed since build')
    moving=record['config']['moving_unit'];visible=[moving]+record['config']['fixed_units'];units={u['id']:u for u in source['instances']}
    if scene.frame_end!=record['frames']+1 or scene.render.fps!=record['fps']:raise ValueError('Saved timing differs')
    if set(json.loads(scene['assembly_visible_units']))!=set(visible):raise ValueError('Visible patch differs')
    geometry={};weights={};trees={}
    expected={u['id']+'/'+c for u in source['instances'] for c in u['chains']}
    actual={o.name for o in scene.objects if o.get('source_reference_chain') is not None}
    if actual!=expected:raise ValueError('Whole source chain inventory changed')
    probe_names={'CA '+key+'/'+c for key in visible for c in units[key]['chains']}
    if {o.name for o in scene.objects if o.type=='MESH'}!=expected|probe_names:raise ValueError('Unexpected/missing display or probe mesh')
    for unit in source['instances']:
        for chain in unit['chains']:
            obj=scene.objects[unit['id']+'/'+chain];raw=dict(np.load(chain_file(assets/'prototypes', chain)))
            coords=np.empty(len(obj.data.vertices)*3);obj.data.vertices.foreach_get('co',coords)
            faces=np.array([list(p.vertices) for p in obj.data.polygons])
            if obj.parent!=scene.objects[unit['id']] or not np.allclose(obj.matrix_basis,np.eye(4),atol=1e-7):raise ValueError('Whole chain moved within unit')
            if coords.size!=raw['vertices'].size or not np.allclose(coords.reshape(-1,3),raw['vertices'],atol=2e-6,rtol=0) or not np.array_equal(faces,raw['faces']):raise ValueError('Source surface changed/cut')
            if obj.hide_render!=(unit['id'] not in visible):raise ValueError('Partial unit visibility')
            changing=unit['id']==moving and any(j['chain']==chain for j in record['joints'])
            if len(obj.modifiers)!=int(changing):raise ValueError('Unexpected surface modifier')
            if changing:
                probe=scene.objects['CA '+unit['id']+'/'+chain]
                if len(probe.modifiers)!=1 or probe.modifiers[0].node_group!=obj.modifiers[0].node_group:raise ValueError('Probe differs from displayed surface deformation')
                numbers=np.empty(len(obj.data.vertices),np.int32);obj.data.attributes['source_residue'].data.foreach_get('value',numbers)
                weights[chain]=numbers;trees[chain]=obj.modifiers[0].node_group
                obj.modifiers[0].show_viewport=False
            if unit['id']!=moving and not np.allclose(obj.parent.matrix_world,unit['matrix_nm'],atol=2e-6,rtol=0):raise ValueError('Fixed source unit moved')
            if chain not in geometry:geometry[chain]=coords.reshape(-1,3)
    probes={key:[scene.objects['CA '+key+'/'+c] for c in units[key]['chains']] for key in visible}
    def coordinates(objects):
        result=[];graph=bpy.context.evaluated_depsgraph_get()
        for obj in objects:
            evaluated=obj.evaluated_get(graph);mesh=evaluated.to_mesh();values=np.empty(len(mesh.vertices)*3)
            mesh.vertices.foreach_get('co',values);matrix=np.array(evaluated.matrix_world)
            result.append(values.reshape(-1,3)@matrix[:3,:3].T+matrix[:3,3]);evaluated.to_mesh_clear()
        return np.concatenate(result)
    values=[];fixed=None;fixed_error=0.;rotations=[];angles=[];times=np.arange(record['frames']*4+1)/4+1
    for time in times:
        frame=int(time);scene.frame_set(frame,subframe=float(time-frame));values.append(coordinates(probes[moving]))
        observed=coordinates([obj for key in record['config']['fixed_units'] for obj in probes[key]])
        if fixed is None:fixed=observed
        fixed_error=max(fixed_error,float(np.max(abs(observed-fixed))))
        rotations.append(np.array(scene.objects[moving].matrix_world))
        angles.append([list(trees[j['chain']].nodes[f'Pivot {j["pivot_residue"]}'].inputs['Rotation'].default_value) for j in record['joints']])
    output.mkdir(parents=True)
    if fixed_error>2e-5:raise ValueError('Fixed probe geometry changed during motion')
    np.savez_compressed(output/'samples.npz',moving=np.array(values,dtype=np.float32),fixed=fixed.astype(np.float32),
        matrices=rotations,eulers=angles,times=times,**{'surface_'+k:v for k,v in geometry.items()},**{'weights_'+k:v for k,v in weights.items()})
    result={'scene_sha256':digest(Path(bpy.data.filepath)),'motion_sha256':digest(folder/'manifest.json'),'samples':len(times),
        'all_source_chains':len(expected),'maximum_fixed_pose_change_nm':fixed_error,
        'shared_display_and_probe_nodes':True,'camera_matrix':np.array(scene.camera.matrix_world).tolist(),
        'ortho_scale':scene.camera.data.ortho_scale,'resolution':[scene.render.resolution_x,scene.render.resolution_y],
        'samples_sha256':digest(output/'samples.npz'),'limits':'Quarter-frame sampling, not a continuous-time or atomic steric proof.'}
    (output/'manifest.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'samples':len(times),'output':str(output)}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['build','export'])
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--motion',type=Path)
    parser.add_argument('--font-zh',type=Path);parser.add_argument('--font-en',type=Path);parser.add_argument('--style-profile',type=Path)
    parser.add_argument('--recipe',type=Path);parser.add_argument('--width',type=int,default=1280);parser.add_argument('--samples',type=int,default=24)
    parser.add_argument('--engine',choices=['eevee','cycles'],default='eevee');args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
    if args.action=='build':build(args)
    else:export(args.output)
