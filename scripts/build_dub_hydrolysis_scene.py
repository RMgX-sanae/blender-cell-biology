"""Build the scoped AMSH-LP transfer example with existing whole-molecule scene tools."""
import argparse,json,sys
from pathlib import Path
import bpy,numpy as np
from mathutils import Matrix

SKILL=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(SKILL/'scripts'))
import detail_scene as ds
import style_nodes as sn
from build_ubiquitin_transfer_scene import keyed_points,envelope,annotation_line,background,paired_composite,digest
from build_phosphotransfer_scene import visibility,animated_mesh
from compositor_stack import build_stack

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--font-zh',type=Path,required=True);p.add_argument('--font-en',type=Path,required=True)
    p.add_argument('--style',choices=['dark','light'],default='dark')
    p.add_argument('--recipe',choices=['neutral','soft-cinema','fluorescent-soft','warm-print','cool-duotone','crisp-neutral'],default='neutral')
    p.add_argument('--engine',choices=['eevee','cycles'],default='eevee')
    p.add_argument('--width',type=int,default=1280);p.add_argument('--samples',type=int,default=24);p.add_argument('--render',action='store_true')
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);a.compositor='clean'
    if a.output.exists():p.error('Use a fresh scene directory')
    if not a.font_zh.is_file() or not a.font_en.is_file():p.error('Supply Songti and Times New Roman fonts')
    assets=a.asset_dir.resolve();source=json.loads((assets/'manifest.json').read_text())
    if source['kind']!='amsh-k63-net-hydrolysis' or not source['passed']:p.error('Prepare and check the source hydrolysis first')
    for name,expected in source['data_sha256'].items():
        if digest(assets/name)!=expected:p.error('Prepared data changed: '+name)
    data=dict(np.load(assets/'motion.npz'));a.frames=source['frames'];a.fps=source['fps'];count=a.frames
    bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
    scene=bpy.context.scene;left=ds.collection('AMSH overview');right=ds.collection('Whole ubiquitins');overlay=ds.collection('Reaction annotations');labels=ds.collection('Figure labels')
    profile=json.loads((SKILL/'assets'/('style-molecular-reference.json' if a.style=='dark' else 'style-studio.json')).read_text())
    profile.update(camera={},bump=0,background='#071f2b' if a.style=='dark' else '#edf0e9')
    camera=ds.studio(profile,a,{'position':[0,0,40],'target':[0,0,0],'scale':22,'light_scale':1.2});camera.matrix_world=Matrix.Translation((0,0,40));camera.data.dof.use_dof=False
    scene.render.resolution_y=2*round(a.width*9/32)
    ids=source['reaction_atom_indices'];c=data['distal'][0,ids['C']];n=data['proximal'][0,ids['NZ']];w=data['water'][0,0]
    front=np.cross(n-c,w-c);front/=np.linalg.norm(front)
    if front@(c-data['enzyme'][0].mean(0))<0:front=-front
    horizontal=data['proximal'][0].mean(0)-data['distal'][0].mean(0);horizontal-=front*(horizontal@front);horizontal/=np.linalg.norm(horizontal)
    up=np.cross(front,horizontal);basis=np.column_stack([horizontal,up,front]);transforms={};parents={}
    framing=ds.collection('Shared source framing')
    for view,keys,offset in [('overview',list(source['atom_ids']),[-5.5,0,0]),('detail',['distal','proximal','water'],[5.5,0,0])]:
        points=np.concatenate([data[k] for k in keys],axis=1)@basis
        lo,hi=points.min(1),points.max(1);center=(lo+hi)/2;kernel=np.ones(13)/13
        center=np.column_stack([np.convolve(np.pad(center[:,i],6,mode='edge'),kernel,mode='valid') for i in range(3)])
        extent=np.maximum(hi-center,center-lo);allowed=np.minimum(8.6/(2*extent[:,0]+.4),6.2/(2*extent[:,1]+.4))
        pad=np.pad(allowed,12,mode='edge');scale=np.array([pad[i:i+25].min() for i in range(count+1)])
        scale=np.convolve(np.pad(scale,6,mode='edge'),kernel,mode='valid')*.98
        matrices=np.repeat(np.eye(4)[None],count+1,axis=0);matrices[:,:3,:3]=basis.T[None]*scale[:,None,None];matrices[:,:3,3]=offset-center*scale[:,None]
        transforms[view]=matrices
        obj=bpy.data.objects.new(view+' frame',None);framing.objects.link(obj);parents[view]=obj
        obj.rotation_mode='QUATERNION';obj.rotation_quaternion=Matrix(basis.T).to_quaternion()
        for i in range(count+1):
            obj.location=matrices[i,:3,3];obj.scale=[scale[i]]*3;obj.keyframe_insert('location',frame=i+1);obj.keyframe_insert('scale',frame=i+1)
    colors={'enzyme':'#8baab8','distal':'#5fbfc3','proximal':'#d78da8','water':'#f4a05d','zinc':'#baabd9'}
    if a.style=='light':colors.update(enzyme='#5a7d90',distal='#006671',proximal='#9a365f',water='#b7541b',zinc='#6d4a9a')
    materials={k:ds.material(k,color,.52,0) for k,color in colors.items()};actors=[]
    for k in source['atom_ids']:
        for view,group in [('overview',left),('detail',right)]:
            if view=='detail' and k in ['enzyme','zinc']:continue
            if view=='detail' and k in ['distal','proximal']:
                edges=np.load(assets/(k+'.npz'))['bonds']
                obj=animated_mesh(k+' '+view,data[k],edges,group,materials[k],.045,1,count)
            else:
                obj=keyed_points(k+' '+view,data[k],group)
                if k in ['water','zinc']:sn.instance_spheres(obj,materials[k],radius=.075 if k=='water' else .13)
                else:envelope(obj,materials[k])
            obj.parent=parents[view]
            obj['component']=k;obj['view']=view;obj['atom_ids']=json.dumps(source['atom_ids'][k]);actors.append(obj.name)
    ink=sn.material('Covalent change highlight',{'material':'paint'},'#fff1b1' if a.style=='dark' else '#422554','paint')
    for view in parents:
        marker=keyed_points('Water oxygen highlight '+view,data['water'],overlay);marker.parent=parents[view]
        marker['display_annotation']=True;marker['component']='water';marker['view']=view
        sn.instance_spheres(marker,sn.material('Water tracking '+view,{'material':'paint'},colors['water'],'paint'),radius=.075)
        for title,first,second,visible in [
            ('K63 isopeptide',data['distal'][:,ids['C']],data['proximal'][:,ids['NZ']],data['chemical_fraction']<.05),
            ('Water-derived terminal oxygen',data['distal'][:,ids['C']],data['water'][:,0],data['chemical_fraction']>.98)]:
            obj=keyed_points(title+' '+view,np.stack([first,second],axis=1),overlay,[(0,1)]);annotation_line(obj,ink,.045);obj.parent=parents[view]
            for i,show in enumerate(visible):obj.hide_render=not bool(show);obj.keyframe_insert('hide_render',frame=i+1)
    zh=bpy.data.fonts.load(str(a.font_zh));en=bpy.data.fonts.load(str(a.font_en));zh.pack();en.pack()
    def text(body,x,y,size):
        curve=bpy.data.curves.new(body,'FONT');curve.body=body;curve.font=zh;curve.font_bold=en
        for c,f in zip(body,curve.body_format):f.use_bold=ord(c)<128 or c in '·—'
        curve.align_x='CENTER';curve.size=size
        if a.style=='light':curve.offset=size*.018
        obj=bpy.data.objects.new(body,curve);labels.objects.link(obj);obj.location=(x,y,12);return obj
    text('AMSH-LP · K63 链水解',0,5.1,.80);text('整体复合物',-5.5,4.1,.64);text('完整泛素与水氧',5.5,4.1,.64)
    text('水氧',7.9,-2.2,.49)
    water_position=np.einsum('nij,nj->ni',transforms['detail'][:,:3,:3],data['water'][:,0])+transforms['detail'][:,:3,3]
    line=np.broadcast_to([[0.,0.,12.],[7.2,-1.6,12.],[8.7,-1.6,12.]],(count+1,3,3)).copy();line[:,0,:2]=water_position[:,:2]
    leader=keyed_points('Water identity leader',line,overlay,[(0,1),(1,2)]);annotation_line(leader,ink,.012)
    for x,k,title in [(-8.4,'enzyme','AMSH-LP'),(-4.4,'distal','远端泛素'),(-.3,'proximal','近端泛素'),(4.2,'water','水中的氧'),(8.4,'zinc','Zn2+')]:
        text(title,x,-4.3,.49);dot=sn.point_object(title+' dot',[[x-1.3,-4.12,12]],overlay);sn.instance_spheres(dot,sn.material(title+' ink',{'material':'paint'},colors[k],'paint'),radius=.09)
    for start,end,body in source['config']['captions']:
        visibility(text(body,0,-5.4,.64),round(start*count)+1,round(end*count))
    ds.path('Panel separator',[(0,-3.5,12),(0,3.55,12)],.016,sn.material('Divider',{'material':'paint'},'#bad6d9' if a.style=='dark' else '#365a64','paint'),overlay)
    bg,bg_record=background(assets,float(np.linalg.norm(transforms['overview'][0,:3,0])),('distal','enzyme','proximal'))
    groups=[left,right,overlay,labels,*bg]
    def isolate(layer,allowed):
        for group in groups:layer.layer_collection.children[group.name].exclude=group not in allowed
    beauty=bpy.context.view_layer;beauty.name='Overview';isolate(beauty,[left]);close=scene.view_layers.new('Hydrolysis detail');isolate(close,[right]);annotations=scene.view_layers.new('Reaction annotations');isolate(annotations,[overlay])
    bg_layers=[]
    for group in bg:
        layer=scene.view_layers.new(group.name+' render');isolate(layer,[group]);bg_layers.append(layer)
    recipe=next(r.copy() for r in json.loads((SKILL/'assets/render-recipes.json').read_text())['recipes'] if r['id']==a.recipe);recipe.update(background=profile['background'],label_min_pixels=25*a.width/1280)
    compositor=build_stack(scene,recipe);paired_composite(scene,close,annotations,bg_layers,profile['background'],a.width)
    for action in bpy.data.actions:
        for layer in action.layers:
            for strip in layer.strips:
                for bag in strip.channelbags:
                    for curve in bag.fcurves:
                        for key in curve.keyframe_points:key.interpolation='CONSTANT' if curve.data_path in ('hide_render','hide_viewport') else 'LINEAR'
    a.output.mkdir(parents=True);scene.frame_set(min(count,round(count*.35)));bpy.ops.wm.save_as_mainfile(filepath=str(a.output/'scene.blend'))
    manifest={'kind':'amsh-k63-net-hydrolysis-scene','preparation':source,'frames':count,'fps':a.fps,'asset_directory':str(assets),'actors':actors,'transforms':{k:v.tolist() for k,v in transforms.items()},'style':a.style,'recipe':a.recipe,'engine':a.engine,'fonts':{'zh':str(a.font_zh),'en':str(a.font_en)},'display_annotations':'Enlarged bond and water-oxygen strokes; same coordinates in both panels.','compositor_report':compositor,'environment':bg_record,'scope':source['scope']}
    (a.output/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    if a.render:
        for frame in [1,round(count*.35),round(count*2/3),round(count*.9)]:
            scene.frame_set(frame);scene.render.filepath=str(a.output/f'key_{frame:04d}.png');bpy.ops.render.render(write_still=True)
    print(json.dumps({'saved':str(a.output/'scene.blend'),'style':a.style}),flush=True)


if __name__=='__main__':main()
