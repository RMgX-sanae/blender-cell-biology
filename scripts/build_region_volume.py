"""Build an editable measured or regionally modeled cell view in the same physical coordinate frame."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Matrix

sys.path.insert(0,str(Path(__file__).resolve().parent))
import detail_scene as ds
from blender_volume import place,audit
from scene_chapters import empty,caption_set,font_pair
from apply_scale_environment import surface


def make_scene(opts):
    p=opts.region;manifest=json.loads((p/'manifest.json').read_text())
    meta=json.loads((p/'surfaces.json').read_text())
    if (meta['source_masks_sha256']!=manifest['data_sha256']['source-masks.npz']
            or meta['source_masks_sha256']!=sha256((p/'source-masks.npz').read_bytes()).hexdigest()
            or sha256((p/'surfaces.npz').read_bytes()).hexdigest()!=meta['surfaces_sha256']):
        raise ValueError('Use context surfaces from the same checked source-mask bundle')
    modeled=opts.layout is not None
    if modeled:
        layout=json.loads((opts.layout/'manifest.json').read_text())
        if (meta['source_masks_sha256']!=layout['files']['source-masks.npz']
                or layout.get('placement_policy')!='replace-measured-structure'):
            raise ValueError('Use a measured nuclear-shell layout from this same cell')
    bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
    scene=bpy.context.scene;group=ds.collection('Registered segmentation')
    anchor=empty('Source voxel coordinate frame',group)
    cellmat=surface('Cell boundary','#418B9A',opacity=.09,bump=0,emission=.09)
    tree=cellmat.node_tree;fresnel=tree.nodes.new('ShaderNodeFresnel');fresnel.inputs['IOR'].default_value=1.38
    optic=tree.nodes.new('ShaderNodeMath');optic.operation='MULTIPLY_ADD'
    optic.inputs[1].default_value=.29;optic.inputs[2].default_value=.022
    tree.links.new(fresnel.outputs[0],optic.inputs[0])
    tree.links.new(optic.outputs[0],next(n for n in tree.nodes if n.bl_idname=='ShaderNodeMixShader').inputs[0])
    nucleus=surface('Nuclear mask','#655A8D',opacity=.25,bump=0,emission=.07)
    marker=surface('TOMM20 marker' if manifest['config']['expected_structure']=='TOMM20' else 'Structure marker','#F7B25F',bump=0,emission=.30)
    with np.load(p/'surfaces.npz') as z:
        points=z['cell_vertices']
        center=(points.max(0)+points.min(0))/2
        anchor.location=-center
        for role,mat in [('cell',cellmat),('nucleus',nucleus),('measured',marker)]:
            obj=ds.mesh('Source '+role,z[role+'_vertices'],z[role+'_faces'],[mat],group)
            obj.data.normals_split_custom_set_from_vertices(z[role+'_normals'])
            obj.parent=anchor;obj['source_role']=role
            obj.hide_render=obj.hide_viewport=bool(modeled and role=='measured')
        v=z['body_vertices'];origin=(v.max(0)+v.min(0))/2
        prototype=ds.mesh('Complete marker component',v-origin,z['body_faces'],[marker],group)
        normals=z['body_normals']
        prototype.data.normals_split_custom_set_from_vertices(normals)
        if modeled:
            with np.load(opts.layout/'body.npz') as obj:
                if not np.array_equal(v-origin,obj['vertices']) or not np.array_equal(z['body_faces'],obj['faces']):
                    raise ValueError('The source-linked marker view requires the selected source component as its exported body')
        prototype['source_center_xyz_um']=origin.tolist()
        prototype.hide_render=True;prototype.hide_viewport=True
    opts.fps,opts.frames,opts.engine,opts.compositor=24,1,'cycles','clean'
    profile=json.loads((Path(__file__).resolve().parents[1]/'assets/style-organelle.json').read_text())
    profile.update(background='071A24',camera={},exposure=.12)
    centered=points-center
    _,vectors=np.linalg.eigh(np.cov(centered[:,:2].T));right=np.r_[vectors[:,-1],0.]
    if right[0]<0:right*=-1
    minor=np.cross([0.,0.,1.],right);front=np.array([0.,0.,.84])-.54*minor;front/=np.linalg.norm(front)
    up=np.cross(front,right)
    xx,yy=centered@right,centered@up
    target=.5*(xx.min()+xx.max())*right+.5*(yy.min()+yy.max())*up
    scale=float(max(np.ptp(xx)/.86,np.ptp(yy)*1.6/.76))
    camera=ds.studio(profile,opts,{'position':target+front*67,'target':target,'scale':scale,'light_scale':1.25})
    camera.rotation_euler=Matrix(np.column_stack((right,up,front)).tolist()).to_euler()
    scene.cycles.transparent_max_bounces=22
    if modeled:
        objects=place(opts.layout,anchor)
        objects[0].data.materials.append(marker)
        for polygon in objects[0].data.polygons:polygon.use_smooth=True
        objects[0].data.normals_split_custom_set_from_vertices(normals)
        audit(opts.layout)
    label=manifest['config']['expected_structure']
    cn,en=(('按核距加权的 '+label+' 标记片段','Nuclear-distance model of '+label+' segments') if modeled else
           (label+' 标记的实测分布','Measured '+label+' distribution'))
    caption_set(scene,camera,font_pair(opts,opts.language),opts.language,[(1,cn,en)],[(1,scale)],outline=True)
    scene['measured_volume_context']=json.dumps({'surfaces':meta,'marker':label,'view':'weighted placements' if modeled else 'measured segmentation'})
    scene.frame_set(1);opts.output.parent.mkdir(parents=True,exist_ok=True)
    bpy.ops.wm.save_as_mainfile(filepath=str(opts.output.resolve()),compress=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--region',type=Path,required=True);p.add_argument('--layout',type=Path)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--width',type=int,default=1600)
    p.add_argument('--samples',type=int,default=24);p.add_argument('--language',choices=['zh','en'],default='zh')
    p.add_argument('--font-zh',type=Path);p.add_argument('--font-latin',type=Path)
    opts=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if opts.output.exists():p.error('Use a fresh output scene')
    make_scene(opts)
