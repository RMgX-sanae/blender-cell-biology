"""Render qualified molecular coordinate sequences from a project display specification.

The renderer does not infer chemistry or motion. Full selected component atoms
remain in the scene; detail framing is a camera window, not a mesh section.
"""
import argparse,json,sys
from pathlib import Path
import bpy,numpy as np
from mathutils import Matrix

sys.path.insert(0,str(Path(__file__).resolve().parent))
import detail_scene as ds
import style_nodes as sn
from event_timing import frame_window
from build_ubiquitin_transfer_scene import keyed_points,envelope,annotation_line,background,paired_composite,digest
from build_phosphotransfer_scene import visibility,animated_mesh
from compositor_stack import build_stack
from molecular_sequence_assets import load_assets

SKILL=Path(ds.__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--asset-dir',type=Path,required=True);parser.add_argument('--display',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--font-zh',type=Path,required=True);parser.add_argument('--font-en',type=Path,required=True)
    parser.add_argument('--style',choices=['dark','light'],default='dark');parser.add_argument('--recipe',default='neutral')
    parser.add_argument('--engine',choices=['eevee','cycles'],default='eevee');parser.add_argument('--width',type=int,default=1280)
    parser.add_argument('--samples',type=int,default=24);parser.add_argument('--render',action='store_true')
    parser.add_argument('--style-profile',type=Path,help='Project role materials, semantic colors, lighting and depth-background settings')
    args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
    if args.output.exists():parser.error('Use a fresh scene directory')
    assets=args.asset_dir.resolve();display=json.loads(args.display.read_text())
    source,data,input_checks=load_assets(assets,display);count=source['frames']
    appearance=json.loads(args.style_profile.read_text()) if args.style_profile else None
    if appearance is not None:
        components={key for view in display['views'] for key in view['actors']}
        roles=appearance.get('component_roles',{});surfaces=appearance.get('surface_roles',{})
        if set(roles)!=components or set(roles.values())-set(surfaces):
            parser.error('Map every displayed component to an explicit configured surface role')
        if set(appearance.get('semantic_colors',{}))-set(display['colors'][args.style]):
            parser.error('Semantic colors must use existing display color identities')
        if any(appearance.get(key) for key in ('camera','representation','palette','compositor')):
            parser.error('Use display views/representations and --recipe; camera, representation, palette and compositor profile controls are not consumed here')
        if 'atom_subdivisions' in appearance and (type(appearance['atom_subdivisions']) is not int or not 1<=appearance['atom_subdivisions']<=5):
            parser.error('Atom display subdivisions must be an integer from 1 to 5')
    if not args.font_zh.is_file() or not args.font_en.is_file():parser.error('Supply Songti and Times New Roman font files')
    recipes=json.loads((SKILL/'assets/render-recipes.json').read_text())['recipes']
    recipe=next((r.copy() for r in recipes if r['id']==args.recipe),None)
    if recipe is None:parser.error('Unknown compositor recipe')
    bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
    scene=bpy.context.scene;groups=[ds.collection('Sequence overview'),ds.collection('Sequence detail')]
    overlays=ds.collection('Sequence annotations');labels=ds.collection('Figure labels');framing=ds.collection('Source view frames')
    profile=json.loads((SKILL/'assets'/('style-molecular-reference.json' if args.style=='dark' else 'style-studio.json')).read_text())
    profile.update(camera={},bump=0,background='#071f2b' if args.style=='dark' else '#edf0e9')
    if appearance is not None:
        profile.update({key:appearance[key] for key in ('lighting','world','view_transform','look','exposure','background') if key in appearance})
    args.fps,args.frames,args.compositor=source['fps'],count,'clean'
    camera=ds.studio(profile,args,{'position':[0,0,40],'target':[0,0,0],'scale':22,'light_scale':1.2})
    camera.matrix_world=Matrix.Translation((0,0,40));camera.data.dof.use_dof=False
    scene.render.resolution_y=2*round(args.width*9/32)
    colors=dict(display['colors'][args.style]);colors.update((appearance or {}).get('semantic_colors',{}))
    parents={};transforms={};actors=[];markers=[]
    for view,group,offset in zip(display['views'],groups,[[-5.5,0,0],[5.5,0,0]]):
        name=view['id'];basis=np.array(view['basis_columns'],float);center=np.array(view['center_nm'],float)
        if basis.shape!=(3,3) or not np.allclose(basis.T@basis,np.eye(3),atol=1e-7) or not np.isclose(np.linalg.det(basis),1):parser.error('Require a proper rigid viewing basis')
        if view['span_nm']<=0:parser.error('Require a positive horizontal view span')
        scale=8.6/view['span_nm'];matrix=np.eye(4);matrix[:3,:3]=basis.T*scale;matrix[:3,3]=offset-center@basis*scale
        parent=bpy.data.objects.new(name+' source frame',None);framing.objects.link(parent);parent.matrix_world=Matrix(matrix)
        parents[name]=parent;transforms[name]=matrix.tolist()
        def annotation_visibility(obj,coordinates,radius,interval=None):
            projected=coordinates@matrix[:3,:3].T+matrix[:3,3]
            visible=(abs(projected[:,:,0]-offset[0]).max(1)+radius*scale<4.3)&(abs(projected[:,:,1]).max(1)+radius*scale<3.45)
            if interval is not None:
                time=np.arange(count+1)/count;visible&=(time>=interval[0])&(time<interval[1])
            for i,show in enumerate(visible):obj.hide_render=not bool(show);obj.keyframe_insert('hide_render',frame=i+1)
            return visible.tolist()
        for key,settings in view['actors'].items():
            if key not in source['atom_ids']:parser.error('Display refers to an absent component: '+key)
            role=appearance['component_roles'][key] if appearance is not None else None
            material=ds.material(name+' '+key,settings.get('color',colors[settings.get('color_key',key)]),.52,0,
                                 settings=appearance['surface_roles'][role] if role else None)
            if role:material['surface_role']=role
            representation=settings['representation']
            if representation=='surface':
                obj=keyed_points(key+' '+name,data[key],group);envelope(obj,material)
            elif representation in ['sticks','spheres']:
                edges=settings.get('bonds',[])
                if representation=='sticks' and 'bond_file' in settings:edges=np.load(assets/settings['bond_file'])['bonds']
                edges=np.asarray(edges,dtype=np.int32).reshape(-1,2)
                obj=animated_mesh(key+' '+name,data[key],edges,group,material,settings.get('radius_nm',.065),1,count+1)
                if 'bond_radius_nm' in settings:
                    for modifier in obj.modifiers:
                        if modifier.type=='NODES':
                            for node in modifier.node_group.nodes:
                                if node.bl_idname=='GeometryNodeCurvePrimitiveCircle':node.inputs['Radius'].default_value=settings['bond_radius_nm']
            else:parser.error('Use surface, sticks or spheres for prepared component display')
            if appearance is not None and 'atom_subdivisions' in appearance:
                for modifier in obj.modifiers:
                    for node in modifier.node_group.nodes:
                        if node.bl_idname=='GeometryNodeMeshIcoSphere':
                            node.inputs['Subdivisions'].default_value=appearance['atom_subdivisions']
            obj.parent=parent;obj['component']=key;obj['view']=name;obj['atom_ids']=json.dumps(source['atom_ids'][key]);actors.append(obj.name)
        for highlight in display.get('highlights',[]):
            if name not in highlight['views']:continue
            key=highlight['component'];indices=highlight['indices']
            coordinates=data[key][:,indices];obj=keyed_points(highlight['id']+' '+name,coordinates,overlays)
            obj.parent=parent;obj['display_annotation']=True
            color=highlight.get('color',colors[highlight.get('color_key',key)]);sn.instance_spheres(obj,sn.material(obj.name,{'material':'paint'},color,'paint'),radius=highlight['radius_nm'],
                subdivisions=(appearance or {}).get('atom_subdivisions',2))
            visible=annotation_visibility(obj,coordinates,highlight['radius_nm'])
            markers.append({'object':obj.name,'component':key,'indices':indices,'view':name,'visible':visible})
        for bond in display.get('bond_annotations',[]):
            if name not in bond['views']:continue
            points=np.stack([data[key][:,index] for key,index in bond['atoms']],axis=1)
            obj=keyed_points(bond['id']+' '+name,points,overlays,[(0,1)]);obj.parent=parent;obj['display_annotation']=True
            ink=sn.material(obj.name+' ink',{'material':'paint'},'#fff1b1' if args.style=='dark' else '#422554','paint')
            radius=bond.get('radius_nm',.024);annotation_line(obj,ink,radius)
            visible=annotation_visibility(obj,points,radius,bond['interval'])
            markers.append({'object':obj.name,'atoms':bond['atoms'],'view':name,'interval':bond['interval'],'visible':visible})
    zh=bpy.data.fonts.load(str(args.font_zh));en=bpy.data.fonts.load(str(args.font_en));zh.pack();en.pack()
    def text(body,x,y,size):
        curve=bpy.data.curves.new(body,'FONT');curve.body=body;curve.font=zh;curve.font_bold=en
        for char,form in zip(body,curve.body_format):form.use_bold=ord(char)<128 or char in '·—'
        curve.align_x='CENTER';curve.size=size
        if args.style=='light':curve.offset=size*.018
        obj=bpy.data.objects.new(body,curve);labels.objects.link(obj);obj.location=(x,y,12);return obj
    text(display['title'],0,5.1,.78)
    for view,x in zip(display['views'],[-5.5,5.5]):text(view['title'],x,4.1,.64)
    for x,(key,title) in zip(np.linspace(-7.8,7.8,len(display['legend'])),display['legend']):
        text(title,float(x),-4.3,.52)
        dot=sn.point_object(title+' color',[[x-1.4,-4.12,12]],overlays)
        sn.instance_spheres(dot,sn.material(title,{'material':'paint'},colors[key],'paint'),radius=.09,
                            subdivisions=(appearance or {}).get('atom_subdivisions',2))
    for begin,end,body in display['captions']:visibility(text(body,0,-5.4,.64),*frame_window(begin,end,count))
    ds.path('Panel separator',[(0,-3.5,12),(0,3.55,12)],.016,sn.material('Divider',{'material':'paint'},'#bad6d9' if args.style=='dark' else '#365a64','paint'),overlays)
    bg,bg_record=background(assets,8.6/display['views'][0]['span_nm'],tuple(display['background_components']),profile=appearance)
    if args.style=='dark' and not (appearance or {}).get('environment_palette'):
        seen=set()
        for group in bg:
            for obj in group.objects:
                for material in obj.data.materials:
                    if material.name in seen:continue
                    seen.add(material.name)
                    shader=material.node_tree.nodes.get('Principled BSDF');color=list(shader.inputs['Base Color'].default_value)
                    color[:3]=[value*.55 for value in color[:3]];shader.inputs['Base Color'].default_value=color;material.diffuse_color=color
    all_groups=[*groups,overlays,labels,*bg]
    def isolate(layer,allowed):
        for group in all_groups:layer.layer_collection.children[group.name].exclude=group not in allowed
    beauty=bpy.context.view_layer;beauty.name='Sequence overview';isolate(beauty,[groups[0]])
    close=scene.view_layers.new('Sequence detail');isolate(close,[groups[1]])
    annotations=scene.view_layers.new('Sequence annotations');isolate(annotations,[overlays])
    bg_layers=[]
    for group in bg:
        layer=scene.view_layers.new(group.name+' render');isolate(layer,[group]);bg_layers.append(layer)
    recipe.update(background=profile['background'],label_min_pixels=25*args.width/1280)
    environment_settings=(appearance or {}).get('environment_compositor',{}).get(args.style,{})
    compositor=build_stack(scene,recipe);paired_composite(scene,close,annotations,bg_layers,profile['background'],args.width,
                                                      environment_settings=environment_settings)
    for action in bpy.data.actions:
        for layer in action.layers:
            for strip in layer.strips:
                for bag in strip.channelbags:
                    for curve in bag.fcurves:
                        for key in curve.keyframe_points:key.interpolation='CONSTANT' if curve.data_path in ('hide_render','hide_viewport') else 'LINEAR'
    args.output.mkdir(parents=True)
    scene['sequence_source_sha256']=digest(assets/'manifest.json');scene['sequence_display_sha256']=digest(args.display)
    (args.output/'display.json').write_bytes(args.display.read_bytes())
    if appearance is not None:
        (args.output/'style-profile.json').write_bytes(args.style_profile.read_bytes())
        scene['sequence_style_profile_sha256']=digest(args.style_profile)
    record={'kind':'molecular-sequence-scene','frames':count,'fps':source['fps'],'preparation':source,'asset_directory':str(assets),'display':display,
            'source_sha256':scene['sequence_source_sha256'],'display_sha256':scene['sequence_display_sha256'],'input_checks':input_checks,
            'actors':actors,'markers':markers,'transforms':transforms,'style':args.style,'recipe':args.recipe,'engine':args.engine,
            'fonts':{'zh':str(args.font_zh),'en':str(args.font_en)},'compositor':compositor,'environment':bg_record,
            'appearance':appearance,'appearance_sha256':digest(args.style_profile) if appearance is not None else None,
            'resolved_colors':colors,'environment_compositor':environment_settings,
            'scope':source['scope'],'limits':'Geometry/chemistry belong to the preparation; window framing and enlarged foreground markers are display choices.'}
    (args.output/'manifest.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
    scene.frame_set(round(count*.25));bpy.ops.wm.save_as_mainfile(filepath=str(args.output/'scene.blend'))
    if args.render:
        for fraction in display.get('key_fractions',[0.,.3,.6,.9]):
            frame=min(count,round(count*fraction)+1);scene.frame_set(frame);scene.render.filepath=str(args.output/f'key_{frame:04d}.png');bpy.ops.render.render(write_still=True)
    print(json.dumps({'saved':str(args.output/'scene.blend'),'actors':len(actors)}),flush=True)


if __name__=='__main__':main()
