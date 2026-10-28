"""Make an editable local protein-interaction cycle from reviewed source actors."""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Matrix, Vector

sys.path.insert(0,str(Path(__file__).resolve().parent))
import detail_scene as ds
import style_nodes as sn
from compositor_stack import isolate_labels


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--style',choices=('clay','cinematic','teaching','illustrative','fluorescence'),default='cinematic')
    p.add_argument('--engine',choices=('eevee','cycles'),default='eevee')
    p.add_argument('--font-en',type=Path,required=True)
    p.add_argument('--crowder-dir',type=Path)
    p.add_argument('--width',type=int,default=1280)
    p.add_argument('--samples',type=int,default=24)
    p.add_argument('--render',action='store_true')
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.output.exists():p.error('Use a fresh output directory')
    if a.style=='illustrative' and a.engine!='eevee':p.error('The selected toon material requires Eevee')
    if not a.font_en.is_file():p.error('Supply an installed Times New Roman Regular font file')
    m=json.loads((a.asset_dir/'manifest.json').read_text())
    if not m['motion']['collision_screen_passed']:
        p.error('The prepared path has unresolved deep-overlap frames; inspect and refine it first')
    poses=np.load(a.asset_dir/'motion.npz')['matrices'];n=m['frames']
    library=json.loads((Path(__file__).resolve().parents[1]/'assets/style-library.json').read_text())
    profile=dict(library['profiles'][a.style])
    if a.style=='illustrative':
        profile['palette']=['#78a8e8','#d69a68','#5d87c0']
    if a.style=='fluorescence':
        profile.update(emission_strength=2.4,emission_depth=.85)
    bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
    molecules=ds.collection('Interaction partners');environment=ds.collection('Environment');labels=ds.collection('Labels')
    rigs={}
    for role in ['anchor','moving']:
        rig=bpy.data.objects.new(role,None);molecules.objects.link(rig);rigs[role]=rig
    points={'anchor':[],'moving':[]};bounds={'anchor':[],'moving':[]};objects=[]
    for actor in m['actors']:
        d=np.load(a.asset_dir/actor['file']);points[actor['role']].append(d['atoms'])
        color=profile['palette'][0 if actor['role']=='anchor' else 1]
        if actor['kind']=='polymer':
            mat=sn.material(actor['id'],profile,color)
            if profile['representation']=='atoms':
                obj=sn.point_object(actor['id'],d['atoms'],molecules)
                sn.instance_spheres(obj,mat,radius=.13,subdivisions=2)
            else:
                obj=ds.mesh(actor['id'],d['vertices'],d['faces'],[mat],molecules)
        else:
            mat=sn.material(actor['id'],profile,'#e9bd65','pbr')
            obj=sn.point_object(actor['id'],d['atoms'],molecules)
            sn.instance_spheres(obj,mat,radius=.11,subdivisions=2)
        obj.parent=rigs[actor['role']];obj['biological_id']=actor['id']
        obj['source_structure']=m['example']['pdb'];obj['source_asset']=str((a.asset_dir/actor['file']).resolve())
        objects.append(obj)
        extent=d['atoms'] if obj.modifiers else d['vertices']
        padding=.14 if obj.modifiers else 0
        lo=extent.min(0)-padding;hi=extent.max(0)+padding
        bounds[actor['role']].append(np.array([[x,y,z] for x in [lo[0],hi[0]] for y in [lo[1],hi[1]] for z in [lo[2],hi[2]]]))
    moving=rigs['moving'];moving.rotation_mode='QUATERNION'
    for f,pose in enumerate(poses,1):
        position,rotation,scale=Matrix(pose).decompose()
        moving.location=position;moving.rotation_quaternion=rotation
        moving.keyframe_insert('location',frame=f);moving.keyframe_insert('rotation_quaternion',frame=f)
    for action in bpy.data.actions:
        for layer in action.layers:
            for strip in layer.strips:
                for bag in strip.channelbags:
                    for curve in bag.fcurves:
                        for key in curve.keyframe_points:key.interpolation='LINEAR'
    anchor=np.concatenate(points['anchor']);mobile=np.concatenate(points['moving'])
    right=mobile.mean(0)-anchor.mean(0);right/=np.linalg.norm(right)
    up=np.eye(3)[np.argmin(abs(right))];up-=right*np.dot(up,right);up/=np.linalg.norm(up)
    front=np.cross(right,up)
    projected=[]
    basis=np.column_stack([right,up,front])
    anchor_bounds=np.concatenate(bounds['anchor']);mobile_bounds=np.concatenate(bounds['moving'])
    for pose in poses:
        xyz=np.concatenate([anchor_bounds,mobile_bounds@pose[:3,:3].T+pose[:3,3]])
        projected.append(xyz@basis)
    projected=np.concatenate(projected);lo=projected.min(0);hi=projected.max(0)
    center=((lo+hi)/2)@basis.T
    width=max((hi[0]-lo[0])*1.25,(hi[1]-lo[1])*16/9*1.5)
    a.fps=m['fps'];a.frames=n;a.compositor='clean'
    studio=json.loads((Path(__file__).resolve().parents[1]/'assets/style-molecular-reference.json').read_text())
    studio['background']=profile['background'];studio['bump']=0;studio['exposure']=0
    ds.studio(studio,a,{'position':(center+front*width*2).tolist(),'target':center.tolist(),'scale':width,'light_scale':1})
    scene=bpy.context.scene
    # Preserve the source coordinate frame; rotate the camera roll to the selected right/up axes.
    camera_matrix=np.eye(4);camera_matrix[:3,:3]=np.column_stack([right,up,front]);camera_matrix[:3,3]=center+front*width*2
    scene.camera.matrix_world=Matrix(camera_matrix)
    scene.camera.data.type='ORTHO';scene.camera.data.ortho_scale=width
    scene.camera.data.dof.use_dof=False
    if a.style=='cinematic':
        scene.camera.data.type='PERSP';scene.camera.data.lens=72
        scene.camera.data.sensor_fit='HORIZONTAL';scene.camera.data.sensor_width=36
        scene.camera.data.dof.use_dof=True
        scene.camera.data.dof.focus_distance=width*2
        scene.camera.data.dof.aperture_fstop=.3
    scene.render.resolution_y=2*round(a.width*9/32)
    scene.render.film_transparent=True
    if hasattr(scene,'eevee'):scene.eevee.taa_render_samples=a.samples
    if profile['compositor']=='ink':
        for layer in scene.view_layers:layer.use_pass_normal=True
    sn.compositor(scene,profile)
    font=bpy.data.fonts.load(str(a.font_en));font.pack()
    ink_color='#20272a' if np.mean(ds.linear(profile['background']))>.4 else '#e8eeed'
    ink=sn.material('Readable labels',profile,ink_color,'paint')
    def text(body,x,y,size,color=None):
        data=bpy.data.curves.new(body,'FONT');data.body=body;data.font=font;data.size=size;data.align_x='CENTER'
        obj=bpy.data.objects.new(body,data);labels.objects.link(obj)
        obj.location=scene.camera.matrix_world@Vector((x,y,-width*2))
        obj.rotation_euler=scene.camera.rotation_euler
        data.materials.append(ink if color is None else sn.material(body,profile,color,'paint'))
        return obj
    text(m['example']['anchor_label']+' / '+m['example']['moving_label'],0,width*.24,width*.035)
    text(m['example']['anchor_label'],-width*.22,-width*.235,width*.028)
    text(m['example']['moving_label'],width*.22,-width*.235,width*.028)
    for role,x,color in [('anchor',-.40,profile['palette'][0]),('moving',.04,profile['palette'][1])]:
        bpy.ops.mesh.primitive_uv_sphere_add(segments=16,ring_count=8,radius=width*.009,
            location=scene.camera.matrix_world@Vector((width*x,-width*.225,-width*2)))
        dot=bpy.context.object;dot.name=role+' legend color'
        dot['legend_for_objects']=json.dumps([actor['id'] for actor in m['actors'] if actor['role']==role])
        dot.data.materials.append(sn.material(dot.name,profile,color,'paint'))
    typography=isolate_labels(scene,profile['background'])
    tree=scene.compositing_node_group
    for node in tree.nodes:
        if node.bl_idname=='CompositorNodeRLayers':node.layer=typography['beauty']
    output=next(n for n in tree.nodes if n.bl_idname=='NodeGroupOutput')
    previous=output.inputs['Image'].links[0].from_socket
    layer=tree.nodes.new('CompositorNodeRLayers');layer.layer=typography['layer']
    over=tree.nodes.new('CompositorNodeAlphaOver')
    tree.links.new(previous,over.inputs['Background']);tree.links.new(layer.outputs['Image'],over.inputs['Foreground'])
    tree.links.new(over.outputs['Image'],output.inputs['Image'])
    if a.crowder_dir:
        d=np.load(a.crowder_dir/'chain_A.npz')
        mat=sn.material('Distant ubiquitin',{**profile,'roughness':.85},'#263f49','pbr')
        for i,(x,y) in enumerate([(-.78,.20),(.85,-.18),(-.85,-.24),(.78,.24)]):
            obj=ds.mesh('Context ubiquitin '+str(i),d['vertices'],d['faces'],[mat],environment)
            obj.location=Vector(center+right*x*width+up*y*width-front*width*(1.1+.15*i))
            obj.rotation_euler=(i*.5,.3*i,.7*i)
            obj['source_structure']='1UBQ';obj['distribution_evidence']='illustrative, not measured concentration'
    scene['interaction_asset_dir']=str(a.asset_dir.resolve())
    scene['animation_evidence']='Illustrative rigid-partner motion with a measured bound-pose anchor; no kinetics'
    scene['animation_mode']=m['mode'];scene['source_limits']=m['example']['limits']
    scene.frame_set(1);bpy.context.view_layer.update()
    a.output.mkdir(parents=True)
    bpy.ops.wm.save_as_mainfile(filepath=str(a.output/'scene.blend'))
    if a.render:
        for frame in [1,n//4+1,n//2+1,3*n//4+1,n]:
            scene.frame_set(frame);scene.render.filepath=str(a.output/f'key_{frame:04d}.png')
            bpy.ops.render.render(write_still=True)
    (a.output/'manifest.json').write_text(json.dumps({'asset':m,'style':a.style,'engine':a.engine,
        'source_objects':[o.name for o in objects],'camera_scale_nm':width,'fps':m['fps'],'frames':n,
        'partner_colors':profile['palette'][:2],'protected_typography':typography,
        'rendered_keyframes':a.render,'font_en':str(a.font_en),'environment_illustrative':bool(a.crowder_dir),
        'camera_dof':'Art-directed for readability; camera optics are not a microscope simulation'},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
