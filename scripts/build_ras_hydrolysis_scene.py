"""Build synchronized Ras hydrolysis or integrated exchange/hydrolysis views."""
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
import cellular_environment as environment
from build_phosphotransfer_scene import animated_mesh, visibility, view_layer


COLORS={'ras':'#70b7cb','gap':'#d68ca2','nucleotide_core':'#80dfb0',
        'phosphate':'#f0c45e','water_oxygen':'#ef967c','magnesium':'#c0a2e1',
        'sos':'#aaa0d5','old_gdp':'#82c9e6','old_mg':'#c0a2e1'}


def captions(motion):
    binding=motion['intervals']['binding'];hydrolysis=motion['intervals']['hydrolysis']
    release=motion['intervals']['release'];depart=release[0]+.22*(release[1]-release[0])
    rows=[(0,binding[0],'Ras 与 GTP 结合'),(binding[0],binding[1],'GAP 接近并稳定反应构象'),
          (binding[1],hydrolysis[0],'催化位点对齐'),
          (hydrolysis[0],hydrolysis[1],'GTP 水解：水中的氧进入 Pi'),
          (hydrolysis[1],depart,'GDP 与 Pi 暂时留在复合物中'),
          (depart,release[1],'Pi 释放，GAP 离开，Ras 复位'),
          (release[1],1,'Ras 保持 GDP 状态')]
    return [row for row in rows if row[1]>row[0]]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--motion-dir',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--font-zh',type=Path,required=True);p.add_argument('--font-en',type=Path,required=True)
    p.add_argument('--width',type=int,default=1280);p.add_argument('--samples',type=int,default=24)
    p.add_argument('--engine',choices=['eevee','cycles'],default='eevee');p.add_argument('--render',action='store_true')
    p.add_argument('--environment-dir',type=Path)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.output.exists():p.error('Use a fresh scene directory')
    if not a.font_zh.is_file() or not a.font_en.is_file():p.error('Supply Songti Regular and Times New Roman Regular')
    motion=json.loads((a.motion_dir/'manifest.json').read_text())
    if not motion['passed'] or motion['kind'] not in ['ras-gap-hydrolysis','ras-state-cycle']:p.error('Complete motion geometry checks first')
    cycle=motion['kind']=='ras-state-cycle';caption_rows=motion['captions'] if cycle else captions(motion)
    assets=Path(motion['asset_dir']);sources=json.loads((assets/'manifest.json').read_text())
    samples=dict(np.load(a.motion_dir/('state-cycle.npz' if cycle else 'hydrolysis.npz')));count=motion['frames']
    data={k:dict(np.load(assets/(k+'.npz'))) for k in ['ras','gap','new_gtp']+(['sos','old_gdp'] if cycle else [])}
    if len(samples['ras'])!=count+1:p.error('Motion frame count mismatch')
    bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
    scene=bpy.context.scene
    left=ds.collection('Overview');right=ds.collection('Reaction detail');labels=ds.collection('Typography')
    envelope=ds.collection('Intact protein envelopes')
    a.fps=motion['fps'];a.frames=count;a.compositor='clean'
    profile=json.loads((Path(__file__).resolve().parents[1]/'assets/style-molecular-reference.json').read_text())
    profile.update(background='#071f2b',bump=0,exposure=0,camera={})
    ds.studio(profile,a,{'position':[0,0,40],'target':[0,0,0],'scale':22,'light_scale':1.4})
    camera=scene.camera;camera.matrix_world=Matrix.Translation((0,0,40));camera.data.type='ORTHO'
    camera.data.ortho_scale=22;camera.data.dof.use_dof=False
    scene.render.resolution_y=2*round(a.width*9/32)
    if hasattr(scene,'eevee'):scene.eevee.taa_render_samples=a.samples
    ras_center=data['ras']['gap'].mean(0)
    horizontal=data['gap']['atoms'].mean(0)-ras_center;horizontal/=np.linalg.norm(horizontal)
    up=samples['water_oxygen'][motion.get('handoff_index',0),0]-ras_center
    up-=horizontal*np.dot(up,horizontal);up/=np.linalg.norm(up)
    front=np.cross(horizontal,up);basis=np.column_stack([horizontal,up,front])
    bound=round(motion['intervals']['hydrolysis'][1]*count)
    overview=np.vstack([data['ras']['atoms'],data['ras']['gap'],data['gap']['atoms']]+
                       ([data['ras']['gef'],data['sos']['vertices']] if cycle else []))@basis
    # A separate camera orientation exposes both reacting P-O vectors. Only the
    # view changes; the detail shares the overview's chemical coordinates.
    detail_basis=np.array([[-.0205326086,.6971455765,-.7166355120],
                          [-.8842978830,-.3470944096,-.3123183071],
                          [-.4664715062,.6273065566,.6236111111]])
    detail=np.vstack([samples[k][f] for k in ['nucleotide_core','phosphate','water_oxygen','magnesium']
                      for f in [motion.get('handoff_index',0),bound]]+
                     ([data['old_gdp']['atoms'],data['old_gdp']['magnesium']] if cycle else []))@detail_basis
    transforms={}
    for name,coordinates,view_basis,offset in [('overview',overview,basis,[-5.5,0,0]),('detail',detail,detail_basis,[5.5,0,0])]:
        lo=coordinates.min(0)-.12;hi=coordinates.max(0)+.12;center=(lo+hi)/2
        scale=min(8.6/(hi[0]-lo[0]),6.2/(hi[1]-lo[1]))
        transforms[name]={'basis':view_basis.tolist(),'center':center.tolist(),'scale':float(scale),'offset':offset}
    def placement(view):
        t=transforms[view];matrix=np.eye(4);matrix[:3,:3]=np.array(t['basis']).T*t['scale']
        matrix[:3,3]=np.array(t['offset'])-np.array(t['center'])*t['scale']
        return Matrix(matrix)
    material_profile={'material':'pbr','roughness':.44,'bump':0}
    materials={k:sn.material(k,material_profile,color) for k,color in COLORS.items()}
    proteins={}
    for key,label in [('ras','Ras G-domain'),('gap','p120GAP catalytic region')]:
        obj=animated_mesh(label,samples[key],np.empty((0,2),dtype=int),envelope,materials[key],0,1,count)
        obj.modifiers.clear()
        density=sn.density(obj,materials[key],{'volume_voxel_nm':.055,'volume_radius_nm':.24},isosurface=True)
        output=next(n for n in density.nodes if n.bl_idname=='NodeGroupOutput')
        surface=output.inputs['Geometry'].links[0].from_socket;smooth=sn.node(density,'GeometryNodeSetShadeSmooth')
        density.links.new(surface,smooth.inputs['Geometry']);density.links.new(smooth.outputs['Geometry'],output.inputs['Geometry'])
        obj.matrix_world=placement('overview');obj['biological_id']=label;obj['atom_count']=len(samples[key][0])
        proteins[key]=obj.name
    sos=None
    if cycle:
        sos=ds.mesh('SOS catalytic region',data['sos']['vertices'].tolist(),data['sos']['faces'].tolist(),[materials['sos']],envelope)
        sos.rotation_mode='QUATERNION';sos['biological_id']='SOS1 modeled catalytic region';sos['source_structure']='1BKD author S'
        for i,pose in enumerate(samples['sos_pose']):
            sos.matrix_world=placement('overview')@Matrix(pose)
            for field in ['location','rotation_quaternion','scale']:sos.keyframe_insert(field,frame=i+1)
    prefix='GTP_B:' if cycle else 'GTP:'
    identities={
        'nucleotide_core':[prefix+n for n in motion['core_atom_names']],
        'phosphate':[prefix+n for n in motion['gamma_atom_names']],
        'water_oxygen':['water:O'],'magnesium':['Mg_B' if cycle else 'Mg:1']}
    if cycle:identities.update(old_gdp=['GDP_A:'+n for n in motion['old_gdp_atom_names']],old_mg=['Mg_A'])
    actors=[]
    for key in identities:
        edges=samples['core_bonds'] if key=='nucleotide_core' else np.array([[0,1],[0,2],[0,3]]) if key=='phosphate' else np.empty((0,2),dtype=int)
        if key=='old_gdp':edges=samples['old_gdp_bonds']
        obj=animated_mesh(key,samples[key],edges,left,materials[key],.085 if key in ['magnesium','old_mg'] else .055,1,count)
        obj.matrix_world=placement('overview');obj['atom_ids']=json.dumps(identities[key]);obj['component']=key;obj['view']='overview'
        copy=ds.duplicate(obj,key+'_detail',right);copy.matrix_world=placement('detail');copy['view']='detail'
        actors.extend([obj.name,copy.name])
    bonds={}
    bridge=motion['core_atom_names'].index('O3B')
    for name,endpoints,weight in [
            ('Leaving P-O3B',np.stack([samples['phosphate'][:,0],samples['nucleotide_core'][:,bridge]],axis=1),1-samples['transfer']),
            ('Forming P-waterO',np.stack([samples['phosphate'][:,0],samples['water_oxygen'][:,0]],axis=1),samples['transfer'])]:
        obj=animated_mesh(name,endpoints,np.array([[0,1]]),left,materials['phosphate'],0,1,count)
        obj.matrix_world=placement('overview');obj['display_only']=True
        circle=next(n for n in obj.modifiers[0].node_group.nodes if n.bl_idname=='GeometryNodeCurvePrimitiveCircle')
        for i,value in enumerate(weight):
            circle.inputs['Radius'].default_value=.026*float(value);circle.inputs['Radius'].keyframe_insert('default_value',frame=i+1)
        copy=ds.duplicate(obj,name+'_detail',right);copy.matrix_world=placement('detail')
        bonds[name]={'object':obj.name,'detail':copy.name,'kind':'leaving' if name.startswith('Leaving') else 'forming'}
    font_zh=bpy.data.fonts.load(str(a.font_zh));font_en=bpy.data.fonts.load(str(a.font_en));font_zh.pack();font_en.pack()
    ink=sn.material('Songti and Times labels',material_profile,'#e8efed','paint')
    def text(body,x,y,size=.6):
        curve=bpy.data.curves.new(body,'FONT');curve.body=body;curve.font=font_zh;curve.font_bold=font_en
        for char,style in zip(body,curve.body_format):style.use_bold=ord(char)<128 or char=='γ'
        curve.size=size;curve.align_x='CENTER';curve.materials.append(ink)
        obj=bpy.data.objects.new(body,curve);labels.objects.link(obj);obj.location=(x,y,8)
        return obj
    text('Ras 开关：交换 → 水解' if cycle else 'GAP 促进 Ras–GTP 水解',0,5.2,.86)
    text('Ras / SOS / GAP' if cycle else 'Ras 与 GAP',-5.5,4.25,.68);text('核苷酸与水中氧',5.5,4.25,.68)
    legend=[(-9,'Ras','ras'),(-5.4,'GAP','gap'),(-1.8,'GDP 骨架','nucleotide_core'),
            (1.8,'γ-P','phosphate'),(5.4,'水中氧','water_oxygen'),(9,'Mg2+','magnesium')]
    if cycle:
        legend=[(x,label,key) for x,(label,key) in zip(np.linspace(-9.6,9.9,8),
                [('Ras','ras'),('SOS','sos'),('GAP','gap'),('旧 GDP','old_gdp'),
                 ('新核苷酸','nucleotide_core'),('γ-P','phosphate'),('水中氧','water_oxygen'),('Mg2+','magnesium')])]
    for x,label,key in legend:
        text(label,x,-4.2,.53 if cycle else .58)
        dot=sn.point_object(label+' color',np.array([[x-(1.21 if cycle else 1.32),-3.96,8]]),labels)
        sn.instance_spheres(dot,sn.material(label+' legend',material_profile,COLORS[key],'paint'),radius=.10)
    for begin,end,caption in caption_rows:
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
    tree=bpy.data.node_groups.new('COMP_Ras_'+('state_cycle' if cycle else 'hydrolysis')+'_synchronized_views','CompositorNodeTree')
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
            faint.inputs[0].default_value=.68;faint.inputs[1].default_value=(*ds.linear(profile['background']),1)
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
    for action in bpy.data.actions:
        for layer in action.layers:
            for strip in layer.strips:
                for bag in strip.channelbags:
                    for curve in bag.fcurves:
                        for key in curve.keyframe_points:key.interpolation='CONSTANT' if curve.data_path in ['hide_render','hide_viewport'] else 'LINEAR'
    scene['motion_directory']=str(a.motion_dir.resolve());scene['animation_mode']=motion['mode']
    scene['section_method']='Intact protein envelopes in overview; whole proteins omitted in detail. No molecular mesh cutting.'
    a.output.mkdir(parents=True);scene.frame_set(1);bpy.ops.wm.save_as_mainfile(filepath=str(a.output/'scene.blend'))
    result={'frames':count,'fps':motion['fps'],'motion':motion,'sources':sources,'transforms':transforms,
        'actors':actors,'proteins':proteins,'reaction_bonds':bonds,'font_zh':str(a.font_zh),'font_en':str(a.font_en),
        'captions':caption_rows,'sos_object':sos.name if sos else None,'environment':environment_record,
        'limits':motion['chemistry_limits']+['The detail omits whole proteins and enlarges the chemical actors independently.',
            'Source-modeled catalytic constructs; no membrane anchor or full upstream pathway.',
            'Continuous motion is a finite sequence with open substrate/product pools, not a closed chemical loop.']}
    (a.output/'manifest.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
    if a.render:
        for fraction in ([0,.18,.25,.31,.44,.5,.68,.72,.76,.86,.94,.99] if cycle else [0,.2,.36,.44,.52,.7,.8,.99]):
            frame=min(count,round(fraction*count)+1);scene.frame_set(frame);scene.render.filepath=str(a.output/f'key_{frame:04d}.png');bpy.ops.render.render(write_still=True)


if __name__=='__main__':main()
