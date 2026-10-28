"""Editable synchronized Ras/SOS and intact nucleotide-exchange views."""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Matrix
sys.path.insert(0,str(Path(__file__).resolve().parent))
import detail_scene as ds
import style_nodes as sn
from build_phosphotransfer_scene import animated_mesh, visibility, view_layer
import cellular_environment as environment


COLORS={'ras':'#70b7cb','sos':'#d68ca2','old_gdp':'#82c9e6','new_gtp':'#80dfb0',
        'gamma':'#f0c45e','old_mg':'#c0a2e1','new_mg':'#c0a2e1'}


def captions(motion):
    (a,b),(c,d)=motion['intervals']
    old=motion['plans']['old_gdp'];new=motion['plans']['new_gtp']
    release=a+(b-a)*old['release_start']
    loaded=c+(d-c)*(1-new['release_start'])
    intervals=[(0,a,'Ras 与 GDP 结合'),(a,release,'SOS 接近 Ras'),(release,b,'GDP 与 Mg2+ 释放'),
        (b,c,'无核苷酸的 Ras–SOS 复合物'),(c,loaded,'新的 GTP 与 Mg2+ 进入'),
        (loaded,d,'SOS 释放'),(d,1,'形成 Ras–GTP 状态')]
    return [row for row in intervals if row[1]>row[0]]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--motion-dir',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--font-zh',type=Path,required=True);p.add_argument('--font-en',type=Path,required=True)
    p.add_argument('--width',type=int,default=1280);p.add_argument('--samples',type=int,default=24)
    p.add_argument('--engine',choices=['eevee','cycles'],default='eevee');p.add_argument('--render',action='store_true')
    p.add_argument('--environment-dir',type=Path,help='Prepared cytosolic microstructures for the preferred depth backdrop')
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.output.exists():p.error('Use a fresh scene directory')
    if not a.font_zh.is_file() or not a.font_en.is_file():p.error('Supply Songti Regular and Times New Roman Regular')
    motion=json.loads((a.motion_dir/'manifest.json').read_text())
    if not motion['passed'] or 'full_timeline_checks' not in motion:p.error('Complete full motion checks before scene generation')
    assets=Path(motion['asset_dir']);sources=json.loads((assets/'manifest.json').read_text())
    data={k:dict(np.load(assets/(k+'.npz'))) for k in ['ras','sos','old_gdp','new_gtp']}
    samples=dict(np.load(a.motion_dir/'exchange.npz'));count=motion['frames']
    if len(samples['ras'])!=count+1:p.error('Motion frame count mismatch')
    bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
    scene=bpy.context.scene;left=ds.collection('Overview');right=ds.collection('Nucleotide detail')
    labels=ds.collection('Typography');envelope=ds.collection('Intact protein envelopes')
    a.fps=motion['fps'];a.frames=count;a.compositor='clean'
    profile=json.loads((Path(__file__).resolve().parents[1]/'assets/style-molecular-reference.json').read_text())
    profile.update(background='#071f2b',bump=0,exposure=0,camera={})
    ds.studio(profile,a,{'position':[0,0,40],'target':[0,0,0],'scale':22,'light_scale':1.4})
    camera=scene.camera;camera.matrix_world=Matrix.Translation((0,0,40));camera.data.type='ORTHO'
    camera.data.ortho_scale=22;camera.data.dof.use_dof=False
    scene.render.resolution_y=2*round(a.width*9/32)
    if hasattr(scene,'eevee'):scene.eevee.taa_render_samples=a.samples
    ras_center=data['ras']['atoms'].mean(0);chemical_center=data['old_gdp']['atoms'].mean(0)
    horizontal=data['sos']['atoms'].mean(0)-ras_center
    horizontal/=np.linalg.norm(horizontal)
    up=np.array(motion['plans']['old_gdp']['nucleotide_direction'])
    up-=horizontal*np.dot(up,horizontal);up/=np.linalg.norm(up)
    front=np.cross(horizontal,up)
    if np.dot(front,chemical_center-ras_center)<0:front=-front;up=-up
    basis=np.column_stack([horizontal,up,front])
    overview=np.vstack([data['ras']['atoms'],data['ras']['gef'],data['sos']['vertices']])@basis
    detail=np.vstack([data['old_gdp']['atoms'],data['new_gtp']['atoms'],data['old_gdp']['magnesium'],data['new_gtp']['magnesium']])@basis
    transforms={}
    for name,coordinates,offset in [('overview',overview,[-5.5,0,0]),('detail',detail,[5.5,0,0])]:
        lo=coordinates.min(0);hi=coordinates.max(0);center=(lo+hi)/2
        if name=='detail':lo-=.10;hi+=.10
        scale=min(8.8/(hi[0]-lo[0]),6.3/(hi[1]-lo[1]))
        transforms[name]={'basis':basis.tolist(),'center':center.tolist(),'scale':float(scale),'offset':offset}
    def placement(view):
        t=transforms[view];matrix=np.eye(4);matrix[:3,:3]=basis.T*t['scale']
        matrix[:3,3]=np.array(t['offset'])-np.array(t['center'])*t['scale']
        return Matrix(matrix)
    material_profile={'material':'pbr','roughness':.44,'bump':0}
    materials={k:sn.material(k,material_profile,color) for k,color in COLORS.items()}
    ras=animated_mesh('Ras G-domain',samples['ras'],np.empty((0,2),dtype=int),envelope,materials['ras'],0,1,count)
    ras.modifiers.clear()
    density=sn.density(ras,materials['ras'],{'volume_voxel_nm':.055,'volume_radius_nm':.24},isosurface=True)
    output=next(n for n in density.nodes if n.bl_idname=='NodeGroupOutput')
    surface=output.inputs['Geometry'].links[0].from_socket;smooth=sn.node(density,'GeometryNodeSetShadeSmooth')
    density.links.new(surface,smooth.inputs['Geometry']);density.links.new(smooth.outputs['Geometry'],output.inputs['Geometry'])
    ras.matrix_world=placement('overview');ras['biological_id']='H-Ras G-domain 1-166';ras['atom_count']=len(data['ras']['atoms'])
    sos=ds.mesh('SOS catalytic region',data['sos']['vertices'].tolist(),data['sos']['faces'].tolist(),[materials['sos']],envelope)
    sos.rotation_mode='QUATERNION'
    for i,pose in enumerate(samples['sos_pose']):
        sos.matrix_world=placement('overview')@Matrix(pose)
        for field in ['location','rotation_quaternion','scale']:sos.keyframe_insert(field,frame=i+1)
    sos['biological_id']='SOS1 modeled catalytic region';sos['source_structure']='1BKD author S'
    actors=[]
    for key in ['old_gdp','new_gtp','old_mg','new_mg']:
        d=data.get(key);edges=d['bonds'] if d is not None else np.empty((0,2),dtype=int)
        ids=[f'{key}:{n}' for n in d['atom_names']] if d is not None else [key+':Mg']
        obj=animated_mesh(key,samples[key],edges,left,materials[key],.085 if key.endswith('mg') else .055,1,count)
        obj.matrix_world=placement('overview');obj['atom_ids']=json.dumps(ids);obj['component']=key;obj['view']='overview'
        copy=ds.duplicate(obj,key+'_detail',right);copy.matrix_world=placement('detail');copy['view']='detail'
        actors.extend([obj.name,copy.name])
    # A separately colored display of the four existing gamma atoms; no new
    # chemical atoms or phosphate product are introduced by this overlay.
    indices=[i for i,n in enumerate(data['new_gtp']['atom_names']) if n in ['PG','O1G','O2G','O3G']]
    gamma=animated_mesh('GTP gamma phosphate highlight',samples['new_gtp'][:,indices],np.empty((0,2),dtype=int),left,materials['gamma'],.060,1,count)
    gamma.matrix_world=placement('overview');gamma['display_only']=True
    copy=ds.duplicate(gamma,gamma.name+'_detail',right);copy.matrix_world=placement('detail')
    font_zh=bpy.data.fonts.load(str(a.font_zh));font_en=bpy.data.fonts.load(str(a.font_en));font_zh.pack();font_en.pack()
    ink=sn.material('Songti and Times labels',material_profile,'#e8efed','paint')
    def text(body,x,y,size=.6):
        curve=bpy.data.curves.new(body,'FONT');curve.body=body;curve.font=font_zh;curve.font_bold=font_en
        for char,style in zip(body,curve.body_format):style.use_bold=ord(char)<128 or char=='γ'
        curve.size=size;curve.align_x='CENTER';curve.materials.append(ink)
        obj=bpy.data.objects.new(body,curve);labels.objects.link(obj);obj.location=(x,y,8)
        return obj
    text('Ras 核苷酸交换',0,5.2,.86);text('Ras 与 SOS',-5.5,4.25,.68);text('完整核苷酸',5.5,4.25,.68)
    for x,label,key in [(-9,'Ras','ras'),(-5.4,'SOS','sos'),(-1.8,'GDP','old_gdp'),(1.8,'GTP','new_gtp'),(5.4,'γ-P','gamma'),(9,'Mg2+','old_mg')]:
        text(label,x,-4.2,.60)
        dot=sn.point_object(label+' color',np.array([[x-1.1,-3.94,8]]),labels)
        sn.instance_spheres(dot,sn.material(label+' legend',material_profile,COLORS[key],'paint'),radius=.10)
    for begin,end,caption in captions(motion):
        obj=text(caption,0,-5.4,.76);visibility(obj,round(begin*count)+1,round(end*count))
    ds.path('Panel separator',[(0,-3.55,8),(0,3.65,8)],22/a.width,
        sn.material('Divider',material_profile,'#dde6e9','paint'),labels)
    environment_groups=[];environment_record=None
    if a.environment_dir:
        environment_groups,environment_record=environment.build(a.environment_dir,transforms['overview']['scale'],count)
    groups=[left,right,labels,envelope,*environment_groups]
    layers=[view_layer(scene,name,allowed,groups) for name,allowed in [
        ('Overview render',[left]),('Detail render',[right]),('Typography render',[labels]),('Protein envelopes render',[envelope])]]
    environment_layers=[view_layer(scene,group.name+' render',[group],groups) for group in environment_groups]
    for layer in list(scene.view_layers):
        if layer not in layers+environment_layers:scene.view_layers.remove(layer)
    tree=bpy.data.node_groups.new('COMP_Ras_exchange_synchronized_views','CompositorNodeTree')
    tree.interface.new_socket(name='Image',in_out='OUTPUT',socket_type='NodeSocketColor');images=[]
    backdrop=environment.composite(tree,environment_layers,profile['background'],a.width) if environment_layers else None
    for index,layer in enumerate(layers[:2]):
        render=sn.node(tree,'CompositorNodeRLayers');render.layer=layer.name
        background=sn.node(tree,'CompositorNodeAlphaOver');background.inputs['Background'].default_value=(*ds.linear(profile['background']),1)
        if backdrop is not None:tree.links.new(backdrop,background.inputs['Background'])
        if index==0:
            whole=sn.node(tree,'CompositorNodeRLayers');whole.layer=layers[3].name
            base=sn.node(tree,'CompositorNodeAlphaOver');base.inputs['Background'].default_value=(*ds.linear(profile['background']),1)
            if backdrop is not None:tree.links.new(backdrop,base.inputs['Background'])
            tree.links.new(whole.outputs['Image'],base.inputs['Foreground'])
            faint=sn.node(tree,'ShaderNodeMixRGB','Intact envelopes with visibility overlay')
            faint.inputs[0].default_value=.65;faint.inputs[1].default_value=(*ds.linear(profile['background']),1)
            if backdrop is not None:tree.links.new(backdrop,faint.inputs[1])
            tree.links.new(base.outputs['Image'],faint.inputs[2]);tree.links.new(faint.outputs[0],background.inputs['Background'])
        tree.links.new(render.outputs['Image'],background.inputs['Foreground']);images.append(background.outputs['Image'])
    mask=sn.node(tree,'CompositorNodeBoxMask');mask.inputs['Position'].default_value=(.25,.5);mask.inputs['Size'].default_value=(.5,1)
    mix=sn.node(tree,'ShaderNodeMixRGB');tree.links.new(mask.outputs['Mask'],mix.inputs[0]);tree.links.new(images[1],mix.inputs[1]);tree.links.new(images[0],mix.inputs[2])
    window=sn.node(tree,'CompositorNodeBoxMask');window.inputs['Position'].default_value=(.5,.5);window.inputs['Size'].default_value=(1,.32)
    crop=sn.node(tree,'ShaderNodeMixRGB');crop.inputs[1].default_value=(*ds.linear(profile['background']),1)
    if backdrop is not None:tree.links.new(backdrop,crop.inputs[1])
    tree.links.new(window.outputs['Mask'],crop.inputs[0]);tree.links.new(mix.outputs[0],crop.inputs[2])
    typography=sn.node(tree,'CompositorNodeRLayers');typography.layer=layers[2].name
    over=sn.node(tree,'CompositorNodeAlphaOver');tree.links.new(crop.outputs[0],over.inputs['Background']);tree.links.new(typography.outputs['Image'],over.inputs['Foreground'])
    output=sn.node(tree,'NodeGroupOutput');tree.links.new(over.outputs['Image'],output.inputs['Image']);scene.compositing_node_group=tree
    for obj in scene.objects:
        for modifier in obj.modifiers:
            if modifier.type!='NODES':continue
            for node in modifier.node_group.nodes:
                if node.bl_idname=='GeometryNodeMeshIcoSphere':node.inputs['Subdivisions'].default_value=3
                if node.bl_idname=='GeometryNodeSetShadeSmooth':node.inputs['Shade Smooth'].default_value=True
    for action in bpy.data.actions:
        for layer in action.layers:
            for strip in layer.strips:
                for bag in strip.channelbags:
                    for curve in bag.fcurves:
                        for key in curve.keyframe_points:key.interpolation='CONSTANT' if curve.data_path in ['hide_render','hide_viewport'] else 'LINEAR'
    scene['motion_directory']=str(a.motion_dir.resolve());scene['animation_mode']=motion['mode']
    scene['section_method']='Intact protein envelopes in overview; whole proteins omitted in detail. No molecular mesh cutting.'
    a.output.mkdir(parents=True);scene.frame_set(1);bpy.ops.wm.save_as_mainfile(filepath=str(a.output/'scene.blend'))
    result={'frames':count,'fps':motion['fps'],'motion':motion,'sources':sources,'transforms':transforms,'actors':actors,
        'font_zh':str(a.font_zh),'font_en':str(a.font_en),'captions':captions(motion),
        'environment':environment_record,
        'limits':['Ras G-domain and SOS catalytic constructs only; no membrane anchor or upstream receptor.',
            'The independently enlarged detail shows intact nucleotides. The overview overlays them on intact protein envelopes.',
            'Nucleotide exchange only: no hydrolysis, GAP catalysis, or full Ras signaling pathway is claimed.']}
    (a.output/'manifest.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
    if a.render:
        for fraction in [.15,.33,.37,.48,.57,.61,.68,.94]:
            frame=round(fraction*count)+1;scene.frame_set(frame);scene.render.filepath=str(a.output/f'key_{frame:04d}.png');bpy.ops.render.render(write_still=True)


if __name__=='__main__':main()
