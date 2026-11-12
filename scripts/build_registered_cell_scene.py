"""Show one Allen cell's registered surfaces without moving its organelles.

Run in Blender. This preserves source display resolution, not membrane topology
or atomic detail. Source OBJ, generated meshes, fonts and renders stay external.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Euler, Matrix, Vector

sys.path.insert(0,str(Path(__file__).resolve().parent))
import detail_scene as ds
from apply_scale_environment import surface
from compositor_stack import build_stack
from prepare_allen_environment import load_cell_bundle


def main():
    root=Path(__file__).resolve().parents[1]
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cell',choices=['er','mito','golgi'],required=True)
    p.add_argument('--display',type=Path,default=root/'assets/registered-cell-display.json')
    p.add_argument('--style',choices=['dark','light'],default='dark')
    p.add_argument('--style-profile',type=Path,help='Optional surface roles, lighting, camera and registered_palette overrides; orthographic display only.')
    p.add_argument('--recipe',default='neutral');p.add_argument('--engine',choices=['cycles','eevee'],default='cycles')
    p.add_argument('--width',type=int,default=1280);p.add_argument('--samples',type=int,default=32)
    p.add_argument('--font-zh',type=Path,required=True);p.add_argument('--font-en',type=Path,required=True)
    p.add_argument('--render',action='store_true')
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.output.exists():p.error('Use a fresh output directory')
    if not a.font_zh.is_file() or not a.font_en.is_file():p.error('Supply Songti and Times New Roman fonts')
    data,sources=load_cell_bundle(a.asset_dir,a.cell)
    display=json.loads(a.display.read_text());colors=display['styles'][a.style].copy()
    profile=json.loads((root/'assets/style-organelle.json').read_text())
    profile.update(background=colors['background'],camera={},exposure=0)
    if a.style_profile:
        profile.update(json.loads(a.style_profile.read_text()))
        colors.update(profile.get('registered_palette',{}))
        colors['background']=profile['background']
    if profile.get('camera',{}).get('type','ORTHO')!='ORTHO':
        raise ValueError('Registered-cell labels and framing require an orthographic camera')
    angles=np.asarray(display['rotation_degrees'],float);span=float(display['display_span'])
    if display['version']!=1 or angles.shape!=(3,) or not np.isfinite(angles).all() or not np.isfinite(span) or span<=0:
        raise ValueError('Expected a finite rigid view rotation and positive common display span')
    recipe=next((r.copy() for r in json.loads((root/'assets/render-recipes.json').read_text())['recipes'] if r['id']==a.recipe),None)
    if recipe is None:p.error('Unknown compositor recipe')
    if profile.get('composited_beauty_denoise',False):
        recipe['filters']=[{'type':'denoise'},*recipe.get('filters',[])]
    bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
    scene=bpy.context.scene;group=ds.collection('Registered source cell');labels=ds.collection('Cell labels')
    rotation=Euler(np.radians(angles)).to_matrix().to_4x4()
    points=np.concatenate([d['vertices'] for d in data.values()])@np.asarray(rotation)[:3,:3].T
    center=(points.min(0)+points.max(0))/2;factor=span/np.ptp(points,axis=0).max()
    transform=Matrix.Scale(factor,4)@Matrix.Translation(-Vector(center))@rotation
    objects=[]
    for role,d in data.items():
        surface_role={'cell':'cell_boundary','dna':'chromatin','structure':'organelle_marker'}[role]
        settings=profile.get('surface_roles',{}).get(surface_role)
        mat=surface('Registered '+role,colors[role],{'cell':.12,'dna':.35,'structure':1}[role],
                    bump=0,settings=settings,emission=0 if settings is not None else .18)
        mat['surface_role']=surface_role
        if role=='cell':
            tree=mat.node_tree;fresnel=tree.nodes.new('ShaderNodeFresnel');fresnel.inputs['IOR'].default_value=1.38
            opacity=tree.nodes.new('ShaderNodeMath');opacity.operation='MULTIPLY_ADD'
            opacity.inputs[1].default_value=.26 if a.style=='dark' else .65
            opacity.inputs[2].default_value=.04 if a.style=='dark' else .08
            tree.links.new(fresnel.outputs[0],opacity.inputs[0])
            tree.links.new(opacity.outputs[0],next(n for n in tree.nodes if n.bl_idname=='ShaderNodeMixShader').inputs[0])
        obj=ds.mesh('Source '+role,d['vertices'].tolist(),d['faces'].tolist(),[mat],group)
        obj.matrix_world=transform;obj['source_role']=role;obj['source_sha256']=sources[role]['sha256']
        obj['evidence']='Source display surface preserved; no remeshing, organelle relocation or membrane topology inference.'
        objects.append(obj)
    a.fps,a.frames,a.frame,a.animate,a.compositor=24,1,1,False,'clean'
    camera=ds.studio(profile,a,{'position':(0,0,30),'target':(0,0,0),'scale':span*1.5,'light_scale':.55})
    # Fit actual rotated source vertices; transformed axis-aligned boxes are
    # excessively loose for these flattened, obliquely viewed cells.
    bpy.context.view_layer.update()
    camera_points=np.asarray([camera.matrix_world.inverted()@Vector(v) for v in (points-center)*factor])
    if np.any(camera_points[:,2]>=-camera.data.clip_start):
        raise ValueError('Registered source geometry must remain in front of the camera')
    extent=2*np.abs(camera_points[:,:2]).max(axis=0)
    if not {'position','target'}.intersection(profile.get('camera',{})):
        extent=np.ptp(points,axis=0)[:2]*factor
    fit=max(extent[0]/.78,extent[1]*(scene.render.resolution_x/scene.render.resolution_y)/.68)
    requested=profile.get('camera',{}).get('scale',fit)
    if not np.isfinite(requested) or requested<fit-1e-6:
        raise ValueError('Requested camera scale clips the registered cell; use at least '+str(fit))
    camera.data.ortho_scale=requested
    scene.cycles.transparent_max_bounces=24
    zh=bpy.data.fonts.load(str(a.font_zh));en=bpy.data.fonts.load(str(a.font_en));zh.pack();en.pack()
    bpy.context.view_layer.update();width=camera.data.ortho_scale
    def label(body,x,y):
        curve=bpy.data.curves.new(body,'FONT');curve.body=body;curve.font=zh;curve.font_bold=en
        for char,form in zip(body,curve.body_format):form.use_bold=ord(char)<128
        curve.size=width*.024;curve.align_x='CENTER'
        if a.style=='light':curve.offset=width*.00024
        obj=bpy.data.objects.new(body,curve);labels.objects.link(obj)
        obj.matrix_world=camera.matrix_world@Matrix.Translation((x*width,y*width,-18))
    label(display['titles'][a.cell],0,.267)
    for x,role in zip([-.29,0,.29],data):
        body=display['labels'][a.cell if role=='structure' else role]
        label(body,x,-.278)
        # Separate colored swatches, not tiny annotation text.
        pos=camera.matrix_world@Vector(((x-.09)*width,-.270*width,-18))
        parts=ds.Parts();parts.ball(pos,width*.006,seg=16,rings=8)
        parts.object(body+' color',[ds.material(body+' color',colors[role],.6,0)],labels)
    recipe['background']=colors['background'];recipe['label_min_pixels']=26*a.width/1280
    result=build_stack(scene,recipe)
    record={'kind':'registered-source-cell-scene','cell':a.cell,'source_assets':sources,'display':display,
        'style_name':a.style,'resolved_colors':colors,
        'asset_manifest_sha256':hashlib.sha256((a.asset_dir/'manifest.json').read_bytes()).hexdigest(),
        'objects':{o['source_role']:o.name for o in objects},'common_matrix':list(map(list,transform)),
        'font_sha256':{k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in [('zh',a.font_zh),('en',a.font_en)]},
        'compositor_report':result,'coordinate_units':'Uncalibrated display units; common rigid transform and uniform scale only.',
        'limits':['One source cell and one marked structure per model; the three examples are different cells.',
                  'Source confocal display surfaces; no cristae, ER lumen, Golgi cisternae, molecular structure or complete-cell inventory inferred.',
                  'Transparency and colors are display choices; no measured abundance, density or concentration.',
                  'No generated placement, topology/contact test, dynamics or quantitative scale bar.']}
    scene['registered_source_cell']=json.dumps(record,ensure_ascii=False)
    ds.save_and_render(a,record)


if __name__=='__main__':main()
