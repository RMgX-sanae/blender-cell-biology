"""Build an editable source-qualified C-alpha conformation view.

Consumes a separately reviewed path; it does not solve protein deformation.
The optional two-leaflet membrane is a schematic source-oriented environment.
"""
import argparse,json,sys
from pathlib import Path
from types import SimpleNamespace
import bpy,numpy as np
from mathutils import Vector
from mathutils.kdtree import KDTree

sys.path.insert(0,str(Path(__file__).resolve().parent))
import detail_scene as ds
import style_nodes as sn
from ca_conformation_assets import load_assets,digest
from coordinate_animation import keyed_points
from interaction_membrane import tail_paths
from compositor_stack import build_stack


def role_material(name,component,profile):
    role=profile['component_roles'][component]
    mat=ds.material(name,profile['semantic_colors'][component],settings=profile['surface_roles'][role])
    mat['surface_role']=role
    return mat


def membrane_geometry(settings,points,origin,profile):
    heads=ds.Parts();tails=ds.Parts();rng=np.random.default_rng(settings['seed'])
    mid=settings['midplane_z_nm']-origin[2];half=settings['half_thickness_nm']
    span=np.array(settings['span_xy_nm'])/2
    tree=KDTree(points.size//3)
    for i,p in enumerate(points.reshape(-1,3)):tree.insert(p,i)
    tree.balance();lipids=[]
    for side in [-1,1]:
        accepted=[]
        for xy in rng.uniform(-span,span,(settings['candidates_per_leaflet'],2)):
            if accepted and np.min(np.linalg.norm(np.array(accepted)-xy,axis=1))<settings['spacing_nm']:continue
            head=np.r_[xy,mid+side*half+rng.uniform(-.06,.06)]
            angle=float(rng.uniform(0,2*np.pi));paths=tail_paths(head,side,angle,half)
            # Conservative display clearance from all input knots. Continuous
            # CA/surface-lipid clearance is not certified by this placement.
            if tree.find(head)[2]<settings['head_CA_clearance_nm'] or any(
                    tree.find(p)[2]<settings['tail_CA_clearance_nm'] for path in paths for p in path):continue
            accepted.append(xy);lipids.append([*head,side,angle]);heads.ball(head,.26,seg=12,rings=8)
            for path in paths:tails.tube(path,.05,sides=6)
    group=ds.collection('Source-oriented schematic membrane')
    h=heads.object('Membrane heads',[role_material('Membrane heads','heads',profile)],group)
    t=tails.object('Membrane tails',[role_material('Membrane tails','tails',profile)],group)
    return [h,t],lipids,np.concatenate([np.array(heads.v),np.array(tails.v)])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('asset-dir','display','style-profile','output','font-zh','font-en'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--engine',choices=['eevee','cycles'],default='cycles');p.add_argument('--width',type=int,default=1600)
    p.add_argument('--samples',type=int,default=32);p.add_argument('--render',action='store_true')
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.output.exists():p.error('Use a fresh scene directory')
    source,data=load_assets(a.asset_dir);display=json.loads(a.display.read_text());profile=json.loads(a.style_profile.read_text())
    components=set(source['chains'])|({'heads','tails'} if display.get('membrane') else set())
    if set(profile['component_roles'])!=components or set(profile['semantic_colors'])!=components:
        p.error('Map every retained chain and environment component to a role and color')
    if set(profile['component_roles'].values())-set(profile['surface_roles']):p.error('Missing surface role')
    if any(profile.get(k) for k in ('camera','representation','palette','compositor')):
        p.error('Use the display specification for camera/density and the neutral protected-label compositor')
    if not a.font_zh.is_file() or not a.font_en.is_file():p.error('Supply Songti and Times New Roman fonts')
    if display['version']!=1 or display['density']['volume_voxel_nm']<=0 or display['density']['volume_radius_nm']<=0:
        p.error('Require a version-1 display and positive density lengths in nm')
    bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
    bpy.context.preferences.filepaths.save_version=0
    origin=data['points'][0].mean(0);points=data['points']-origin
    group=ds.collection('Complete selected modeled chains');probes=ds.collection('Coordinate probes');actors=[]
    for chain in source['chains']:
        ids=np.flatnonzero(data['keys'][:,0]==chain)
        obj=keyed_points('Chain '+chain,points[:,ids],group,frames=data['frames']);obj['source_indices']=json.dumps(ids.tolist())
        obj['source_ids']=json.dumps([source['point_ids'][i] for i in ids]);actors.append(obj.name)
        probe=bpy.data.objects.new('CA probe '+chain,obj.data);probes.objects.link(probe);probe.hide_render=True
        material=role_material('Chain '+chain,chain,profile)
        tree=sn.density(obj,material,display['density'],isosurface=True)
        setter=next(n for n in tree.nodes if n.bl_idname=='GeometryNodeSetMaterial')
        output=next(n for n in tree.nodes if n.bl_idname=='NodeGroupOutput')
        smooth=tree.nodes.new('GeometryNodeSetShadeSmooth')
        tree.links.new(setter.outputs['Geometry'],smooth.inputs['Geometry']);tree.links.new(smooth.outputs['Geometry'],output.inputs['Geometry'])
    environment=[];lipids=[];cloud=points.reshape(-1,3)
    if display.get('membrane'):
        environment,lipids,context=membrane_geometry(display['membrane'],points,origin,profile);cloud=np.concatenate([cloud,context])
    timing=source['config']['presentation'];settings=SimpleNamespace(engine=a.engine,samples=a.samples,width=a.width,
                 frames=timing['frames']+1,fps=timing['fps'],compositor='clean')
    direction=np.array(display['view_direction'],float)
    if direction.shape!=(3,) or not np.isfinite(direction).all() or np.linalg.norm(direction)<1e-6:p.error('Invalid viewing direction')
    direction/=np.linalg.norm(direction);target=(cloud.min(0)+cloud.max(0))/2
    camera=ds.studio({**profile,'compositor':{}},settings,{'position':(target+direction*55).tolist(),'target':target.tolist(),'scale':32})
    scene=bpy.context.scene;scene.render.resolution_y=2*round(a.width*9/32)
    bpy.context.view_layer.update();basis=np.array(camera.rotation_euler.to_matrix());projected=(cloud-np.array(camera.location))@basis
    low,high=projected.min(0),projected.max(0)
    camera.location+=Vector(basis@np.array([(low[0]+high[0])/2,(low[1]+high[1])/2,0]))
    camera.data.ortho_scale=max((high[0]-low[0]+1.2)/.90,(high[1]-low[1]+1.2)*16/9/.72)
    focus=points.mean((0,1));camera.data.dof.use_dof=display.get('fstop') is not None
    if camera.data.dof.use_dof:
        camera.data.dof.aperture_fstop=display['fstop'];camera.data.dof.focus_distance=float(np.dot(np.array(camera.location)-focus,direction))
    labels=ds.collection('Figure labels');zh=bpy.data.fonts.load(str(a.font_zh));en=bpy.data.fonts.load(str(a.font_en));zh.pack();en.pack()
    for body,x,y,pixels in [(display['title'],-.45,.41,36),(display['subtitle'],-.45,-.44,29)]:
        curve=bpy.data.curves.new(body,'FONT');curve.body=body;curve.font=zh;curve.font_bold=en
        for char,item in zip(body,curve.body_format):item.use_bold=ord(char)<128 or 0x370<=ord(char)<=0x3ff
        size=camera.data.ortho_scale;curve.size=size*pixels/1600
        obj=bpy.data.objects.new(body,curve);labels.objects.link(obj);obj.parent=camera;obj.location=(x*size,y*size*9/16,-20)
    build_stack(scene,{'id':'ca-conformation-clean','filters':[],'background':profile['background'],'label_min_pixels':29*a.width/1600})
    scene.frame_set(1+timing['frames']//2)
    record={'source_manifest_sha256':digest(a.asset_dir/'manifest.json'),'display_sha256':digest(a.display),
            'style_sha256':digest(a.style_profile),'origin_nm':origin.tolist(),'actors':actors,'lipids':lipids,
            'environment_mesh_sha256':{},'asset_dir':str(a.asset_dir.resolve()),'display':display,'style':profile,
            'source_scope':source['scope'],'full_movie_rendered':False,'command':sys.argv}
    import hashlib
    for obj in environment:
        xyz=np.empty(len(obj.data.vertices)*3,dtype=np.float32);obj.data.vertices.foreach_get('co',xyz)
        record['environment_mesh_sha256'][obj.name]=hashlib.sha256(xyz.tobytes()).hexdigest()
    scene['ca_conformation']=json.dumps(record,ensure_ascii=False)
    a.output.mkdir(parents=True);(a.output/'manifest.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
    bpy.ops.wm.save_as_mainfile(filepath=str(a.output/'scene.blend'))
    if a.render:
        scene.render.filepath=str(a.output/'preview.png');bpy.ops.render.render(write_still=True)


if __name__=='__main__':main()
