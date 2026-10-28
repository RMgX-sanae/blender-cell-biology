"""Build an editable, synchronized PKA overview and complete-substrate detail."""
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


COLORS={'peptide':'#de829c','nucleotide_core':'#75c6de','phosphoryl':'#f0c45e',
        'M1':'#94d696','M2':'#94d696'}
CAPTIONS=[(0,.04,'新一轮催化循环'),(.04,.23,'ATP 与 Mg2+ 结合'),
          (.23,.35,'第二个 Mg2+ 到位'),(.35,.50,'SP20 底物进入'),
          (.50,.54,'反应位点对齐'),(.54,.60,'磷酸基团转移'),
          (.60,.69,'磷酸化产物形成'),(.69,.82,'磷酸化肽释放'),
          (.82,.89,'Mg2+ 释放'),(.89,1,'ADP 与 Mg2+ 离开')]


def visibility(obj,start,end):
    for frame,hidden in [(1,True),(start,False),(end+1,True)]:
        obj.hide_render=hidden;obj.hide_viewport=hidden
        obj.keyframe_insert('hide_render',frame=frame);obj.keyframe_insert('hide_viewport',frame=frame)


def animated_mesh(name,points,edges,group,material,radius,start,end):
    mesh=bpy.data.meshes.new(name);mesh.from_pydata(points[0].tolist(),edges.tolist(),[]);mesh.update()
    obj=bpy.data.objects.new(name,mesh);group.objects.link(obj)
    for index,coordinates in enumerate(points):
        key=obj.shape_key_add(name=f't_{index:04d}')
        key.data.foreach_set('co',coordinates.astype(np.float32).ravel())
        key.interpolation='KEY_LINEAR'
    keys=mesh.shape_keys;keys.use_relative=False
    # Absolute key times accumulate small float32 offsets in Blender. Key each
    # sampled time explicitly so the two views reproduce the prepared frames.
    for index,key in enumerate(keys.key_blocks):
        keys.eval_time=key.frame;keys.keyframe_insert('eval_time',frame=start+index)
    tree,source,target=sn.geometry_tree(obj,'GN_atoms_and_covalent_bonds')
    join=sn.node(tree,'GeometryNodeJoinGeometry')
    if radius:
        sphere=sn.node(tree,'GeometryNodeMeshIcoSphere');sphere.inputs['Radius'].default_value=radius
        sphere.inputs['Subdivisions'].default_value=2
        shade=sn.node(tree,'GeometryNodeSetShadeSmooth');tree.links.new(sphere.outputs['Mesh'],shade.inputs['Geometry'])
        paint=sn.node(tree,'GeometryNodeSetMaterial');paint.inputs['Material'].default_value=material
        tree.links.new(shade.outputs['Geometry'],paint.inputs['Geometry'])
        instance=sn.node(tree,'GeometryNodeInstanceOnPoints')
        tree.links.new(source,instance.inputs['Points']);tree.links.new(paint.outputs['Geometry'],instance.inputs['Instance'])
        tree.links.new(instance.outputs['Instances'],join.inputs['Geometry'])
    if len(edges):
        curve=sn.node(tree,'GeometryNodeMeshToCurve');tree.links.new(source,curve.inputs['Mesh'])
        circle=sn.node(tree,'GeometryNodeCurvePrimitiveCircle');circle.inputs['Resolution'].default_value=8
        circle.inputs['Radius'].default_value=.026
        solid=sn.node(tree,'GeometryNodeCurveToMesh');tree.links.new(curve.outputs['Curve'],solid.inputs['Curve'])
        tree.links.new(circle.outputs['Curve'],solid.inputs['Profile Curve'])
        paint=sn.node(tree,'GeometryNodeSetMaterial');paint.inputs['Material'].default_value=material
        tree.links.new(solid.outputs['Mesh'],paint.inputs['Geometry']);tree.links.new(paint.outputs['Geometry'],join.inputs['Geometry'])
    tree.links.new(join.outputs['Geometry'],target)
    visibility(obj,start,end)
    return obj


def view_layer(scene,name,allowed,groups):
    layer=scene.view_layers.new(name)
    for group in groups:
        layer.layer_collection.children[group.name].exclude=group not in allowed
    return layer


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--motion-dir',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--font-zh',type=Path,required=True);p.add_argument('--font-en',type=Path,required=True)
    p.add_argument('--width',type=int,default=1280);p.add_argument('--samples',type=int,default=24)
    p.add_argument('--engine',choices=['eevee','cycles'],default='eevee');p.add_argument('--render',action='store_true')
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.output.exists():p.error('Use a fresh scene directory')
    if not a.font_zh.is_file() or not a.font_en.is_file():p.error('Supply Songti Regular and Times New Roman Regular font files')
    motion=json.loads((a.motion_dir/'manifest.json').read_text())
    if not motion['passed']:p.error('Fix failed motion/geometry checks before building the scene')
    assets=Path(motion['asset_dir']);source=json.loads((assets/'manifest.json').read_text())
    data={k:dict(np.load(assets/(k+'.npz'))) for k in ['kinase','peptide','nucleotide_core','phosphoryl','metals']}
    samples=dict(np.load(a.motion_dir/'cycle.npz'));count=motion['frames_per_cycle']
    if len(samples['peptide'])!=count+1:p.error('Motion frame count mismatch')
    bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
    scene=bpy.context.scene;left=ds.collection('Overview');right=ds.collection('Substrate detail');labels=ds.collection('Typography')
    envelope=ds.collection('Intact enzyme envelope')
    a.fps=motion['fps'];a.frames=motion['frames'];a.compositor='clean'
    profile=json.loads((Path(__file__).resolve().parents[1]/'assets/style-molecular-reference.json').read_text())
    profile.update(background='#071f2b',bump=0,exposure=0,camera={})
    ds.studio(profile,a,{'position':[0,0,40],'target':[0,0,0],'scale':22,'light_scale':1.4})
    camera=scene.camera;camera.matrix_world=Matrix.Translation((0,0,40));camera.data.type='ORTHO'
    camera.data.ortho_scale=22;camera.data.dof.use_dof=False
    scene.render.resolution_y=2*round(a.width*9/32)
    if hasattr(scene,'eevee'):scene.eevee.taa_render_samples=a.samples
    chemical_center=data['phosphoryl']['atoms'][0]
    front=chemical_center-data['kinase']['atoms'].mean(0);front/=np.linalg.norm(front)
    horizontal=data['peptide']['atoms'].mean(0)-data['nucleotide_core']['atoms'].mean(0)
    horizontal-=front*np.dot(horizontal,front);horizontal/=np.linalg.norm(horizontal)
    up=np.cross(front,horizontal);basis=np.column_stack([horizontal,up,front])
    detail=np.concatenate([data[k][field] for k in ['peptide','nucleotide_core'] for field in ['atoms','product_atoms']])@basis
    overview=data['kinase']['vertices']@basis
    transforms={}
    for name,coordinates,offset in [('overview',overview,[-5.5,0,0]),('detail',detail,[5.5,0,0])]:
        lo=coordinates.min(0);hi=coordinates.max(0);center=(lo+hi)/2
        scale=min(8.8/(hi[0]-lo[0]),6.3/(hi[1]-lo[1]))
        transforms[name]={'basis':basis.tolist(),'center':center.tolist(),'scale':float(scale),'offset':offset}
    def placement(view):
        t=transforms[view];matrix=np.eye(4);matrix[:3,:3]=basis.T*t['scale']
        matrix[:3,3]=np.array(t['offset'])-np.array(t['center'])*t['scale']
        return Matrix(matrix)
    mat_profile={'material':'pbr','roughness':.44,'bump':0}
    materials={k:sn.material(k,mat_profile,color) for k,color in COLORS.items()}
    kinase=data['kinase'];shell=sn.material('Intact PKAc envelope',mat_profile,'#8faebc','pbr')
    enzyme=sn.point_object('PKAc intact overview',kinase['atoms'],envelope)
    density=sn.density(enzyme,shell,{'volume_voxel_nm':.055,'volume_radius_nm':.24},isosurface=True)
    output=next(n for n in density.nodes if n.bl_idname=='NodeGroupOutput')
    surface=output.inputs['Geometry'].links[0].from_socket
    smooth=sn.node(density,'GeometryNodeSetShadeSmooth')
    density.links.new(surface,smooth.inputs['Geometry']);density.links.new(smooth.outputs['Geometry'],output.inputs['Geometry'])
    enzyme.matrix_world=placement('overview');enzyme['biological_id']='PKAc';enzyme['source_structure']='4IAC; 4IAD; 1J3H'
    enzyme.shape_key_add(name='closed')
    for name,field in [('product','product_atoms'),('apo','open_atoms'),('bound_open','bound_open_atoms'),('product_bound_open','product_bound_open_atoms')]:
        key=enzyme.shape_key_add(name=name);key.data.foreach_set('co',kinase[field].astype(np.float32).ravel())
    for cycle in range(motion['cycles']):
        for index,weights in enumerate(samples['enzyme_weights']):
            for j,name in enumerate(motion['enzyme_states'][1:],1):
                key=enzyme.data.shape_keys.key_blocks[name];key.value=float(weights[j]);key.keyframe_insert('value',frame=cycle*count+index+1)
    actors=[]
    for cycle in range(motion['cycles']):
        start=cycle*count+1;end=(cycle+1)*count
        for key,color in COLORS.items():
            d=data.get(key);edges=d['bonds'] if d is not None else np.empty((0,2),dtype=int)
            if key=='peptide':ids=[f'SP20_{cycle+1}:{r}:{n}' for r,n in zip(d['residue_numbers'],d['atom_names'])]
            elif key in ['nucleotide_core','phosphoryl']:ids=[f'ATP_{cycle+1}:{n}' for n in d['atom_names']]
            else:ids=[f'{key}_{cycle+1}']
            obj=animated_mesh(f'packet_{cycle+1}_{key}',samples[key],edges,left,materials[key],.085 if key.startswith('M') else .072,start,end)
            obj.matrix_world=placement('overview');obj['atom_ids']=json.dumps(ids);obj['component']=key;obj['packet']=cycle+1;obj['view']='overview'
            copy=ds.duplicate(obj,obj.name+'_detail',right);copy.matrix_world=placement('detail');copy['view']='detail'
            actors.extend([obj.name,copy.name])
        for name,component,atom in [('ATP donor bond','nucleotide_core','O3B'),('Product acceptor bond','peptide','OG')]:
            d=data[component]
            index=next(i for i,n in enumerate(d['atom_names']) if n==atom and (component!='peptide' or d['residue_numbers'][i]==621))
            coordinates=np.stack([samples[component][:,index],samples['phosphoryl'][:,0]],axis=1)
            obj=animated_mesh(f'{cycle+1} {name}',coordinates,np.array([[0,1]]),left,materials['phosphoryl'],0,start,end)
            obj.matrix_world=placement('overview')
            for i,points in enumerate(coordinates[:-1]):
                hidden=bool(np.linalg.norm(points[0]-points[1])>.19)
                obj.hide_render=hidden;obj.hide_viewport=hidden
                obj.keyframe_insert('hide_render',frame=start+i);obj.keyframe_insert('hide_viewport',frame=start+i)
            copy=ds.duplicate(obj,obj.name+'_detail',right);copy.matrix_world=placement('detail')
    font_zh=bpy.data.fonts.load(str(a.font_zh));font_en=bpy.data.fonts.load(str(a.font_en));font_zh.pack();font_en.pack()
    ink=sn.material('Songti and Times labels',mat_profile,'#e8efed','paint')
    def text(body,x,y,size=.6):
        curve=bpy.data.curves.new(body,'FONT');curve.body=body;curve.font=font_zh;curve.font_bold=font_en
        for char,style in zip(body,curve.body_format):style.use_bold=ord(char)<128 or char in '²⁺·'
        curve.size=size;curve.align_x='CENTER';curve.materials.append(ink)
        obj=bpy.data.objects.new(body,curve);labels.objects.link(obj);obj.location=(x,y,8)
        return obj
    text('PKA 催化循环',0,5.2,.86);text('完整激酶',-5.5,4.25,.68);text('底物与磷酸转移',5.5,4.25,.68)
    for x,label,key in [(-8,'SP20 底物','peptide'),(-2.9,'ATP / ADP','nucleotide_core'),(2.5,'磷酸基团','phosphoryl'),(7.5,'Mg2+','M1')]:
        text(label,x,-4.2,.60)
        dot=sn.point_object(label+' color',np.array([[x-1.65,-3.94,8]]),labels)
        sn.instance_spheres(dot,sn.material(label+' ink',mat_profile,COLORS[key],'paint'),radius=.10)
    for cycle in range(motion['cycles']):
        for begin,end,caption in CAPTIONS:
            obj=text(caption,0,-5.4,.76)
            visibility(obj,cycle*count+round(begin*count)+1,cycle*count+round(end*count))
    ds.path('Panel separator',[(0,-3.55,8),(0,3.65,8)],.009,sn.material('Divider',mat_profile,'#37515c','paint'),labels)
    groups=[left,right,labels,envelope]
    overview_layer=view_layer(scene,'Overview render',[left],groups)
    detail_layer=view_layer(scene,'Detail render',[right],groups)
    typography_layer=view_layer(scene,'Typography render',[labels],groups)
    envelope_layer=view_layer(scene,'Enzyme envelope render',[envelope],groups)
    for layer in list(scene.view_layers):
        if layer not in [overview_layer,detail_layer,typography_layer,envelope_layer]:scene.view_layers.remove(layer)
    tree=bpy.data.node_groups.new('COMP_PKA_synchronized_views','CompositorNodeTree')
    tree.interface.new_socket(name='Image',in_out='OUTPUT',socket_type='NodeSocketColor')
    images=[]
    for layer in [overview_layer,detail_layer]:
        render=sn.node(tree,'CompositorNodeRLayers');render.layer=layer.name
        background=sn.node(tree,'CompositorNodeAlphaOver');background.inputs['Background'].default_value=(*ds.linear(profile['background']),1)
        if layer==overview_layer:
            whole=sn.node(tree,'CompositorNodeRLayers');whole.layer=envelope_layer.name
            base=sn.node(tree,'CompositorNodeAlphaOver');base.inputs['Background'].default_value=(*ds.linear(profile['background']),1)
            tree.links.new(whole.outputs['Image'],base.inputs['Foreground'])
            faint=sn.node(tree,'ShaderNodeMixRGB','Low-contrast intact envelope; illustrative visibility overlay')
            faint.inputs[0].default_value=.30;faint.inputs[1].default_value=(*ds.linear(profile['background']),1)
            tree.links.new(base.outputs['Image'],faint.inputs[2]);tree.links.new(faint.outputs[0],background.inputs['Background'])
        tree.links.new(render.outputs['Image'],background.inputs['Foreground']);images.append(background.outputs['Image'])
    mask=sn.node(tree,'CompositorNodeBoxMask');mask.inputs['Position'].default_value=(.25,.5);mask.inputs['Size'].default_value=(.5,1)
    mix=sn.node(tree,'ShaderNodeMixRGB');tree.links.new(mask.outputs['Mask'],mix.inputs[0]);tree.links.new(images[1],mix.inputs[1]);tree.links.new(images[0],mix.inputs[2])
    window=sn.node(tree,'CompositorNodeBoxMask','Keep entering/leaving molecules within the view windows')
    window.inputs['Position'].default_value=(.5,.5);window.inputs['Size'].default_value=(1,.32)
    crop=sn.node(tree,'ShaderNodeMixRGB');crop.inputs[1].default_value=(*ds.linear(profile['background']),1)
    tree.links.new(window.outputs['Mask'],crop.inputs[0]);tree.links.new(mix.outputs[0],crop.inputs[2])
    typography=sn.node(tree,'CompositorNodeRLayers');typography.layer=typography_layer.name
    over=sn.node(tree,'CompositorNodeAlphaOver');tree.links.new(crop.outputs[0],over.inputs['Background']);tree.links.new(typography.outputs['Image'],over.inputs['Foreground'])
    output=sn.node(tree,'NodeGroupOutput');tree.links.new(over.outputs['Image'],output.inputs['Image']);scene.compositing_node_group=tree
    for action in bpy.data.actions:
        for layer in action.layers:
            for strip in layer.strips:
                for bag in strip.channelbags:
                    for curve in bag.fcurves:
                        for key in curve.keyframe_points:key.interpolation='CONSTANT' if curve.data_path in ['hide_render','hide_viewport'] else 'LINEAR'
    scene['motion_directory']=str(a.motion_dir.resolve());scene['animation_mode']=motion['mode'];scene['source_evidence']=source['source_model']
    scene['section_method']='Whole kinase omitted in detail; no molecule mesh is cut. Views use different magnifications.'
    scene['envelope_representation']='Isosurface rebuilt from all source atom positions each frame; no displaced surface vertices.'
    a.output.mkdir(parents=True);scene.frame_set(1);bpy.context.view_layer.update()
    bpy.ops.wm.save_as_mainfile(filepath=str(a.output/'scene.blend'))
    if a.render:
        for fraction in [.18,.45,.57,.65,.76,.94]:
            frame=round(fraction*count)+1;scene.frame_set(frame);scene.render.filepath=str(a.output/f'key_{frame:04d}.png');bpy.ops.render.render(write_still=True)
    result={'fps':a.fps,'frames':a.frames,'motion':motion,'sources':source,'transforms':transforms,'actors':actors,
            'font_zh':str(a.font_zh),'font_en':str(a.font_en),'rendered_keyframes':a.render,
            'limits':['The detail omits the whole kinase and is enlarged independently.',
                      'Overview substrates are composited over a low-contrast intact envelope for visibility; this is not optical transparency or a cut surface.',
                      'A complete 20-residue engineered SP20 substrate is shown, not a native full-length target protein.',
                      'The absence of a rendered bond during transfer is a teaching convention, not an observed free-metaphosphate intermediate.']}
    (a.output/'manifest.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
