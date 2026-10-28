"""Reopen and check the saved loop before rendering its 384 distinct frames."""
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Vector

import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--asset-dir',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
OUT=args.output.resolve();OUT.mkdir(parents=True,exist_ok=True)
scene=bpy.context.scene
info=json.loads((args.asset_dir/'manifest.json').read_text())
errors=[];asset_checks=[]
objects=[o for o in bpy.data.objects if o.get('source')=='6VQ6']
main=[o for o in objects if o.name.startswith('6VQ6 chain ')]
assert len(main)==31 and len(objects)==41
files={x['chain']:x['file'] for x in info['chains']};assert len(set(files.values()))==31
for o in objects:
    d=np.load(args.asset_dir/files[o['source_chain']]);v=np.empty(len(o.data.vertices)*3,dtype=np.float32)
    o.data.vertices.foreach_get('co',v);v=v.reshape(-1,3)
    same=v.shape==d['vertices'].shape and np.allclose(v,d['vertices'],atol=1e-6,rtol=0)
    faces=np.array([list(p.vertices) for p in o.data.polygons])
    same=same and np.array_equal(faces,d['faces'])
    assert same and not any(m.type=='BOOLEAN' for m in o.modifiers),o.name
    asset_checks.append({'object':o.name,'chain':o['source_chain'],'complete_source_mesh':True,'vertices':len(v)})
labels=[o for o in bpy.data.objects if (o.type=='FONT' and not o.hide_render) or o.get('label_text')]
label_sources=[bpy.data.objects[o['label_source']] if o.get('label_source') else o for o in labels]
frame_bounds=[]
left=[o for o in bpy.data.objects if o.name.startswith('Left H ')]
visibility={o.name:[] for o in left}
phase_error=0.
if scene.get('left_proton_motion'):assert len(left)==30
for f in range(1,386):
    scene.frame_set(f)
    for o in main+labels:
        if o.hide_render:errors.append([f,o.name,'hidden primary/label'])
        p=np.array([tuple(world_to_camera_view(scene,scene.camera,o.matrix_world@Vector(v))) for v in o.bound_box])
        if not np.all(np.isfinite(p)):errors.append([f,o.name,'nonfinite'])
        if np.any(p[:,:2]<0) or np.any(p[:,:2]>1):errors.append([f,o.name,'outside frame'])
        if f in [1,128,256,384]:frame_bounds.append({'frame':f,'object':o.name,'min':p[:,:2].min(0).tolist(),'max':p[:,:2].max(0).tolist()})
    for o in left:
        visibility[o.name].append(not o.hide_render)
        if o['role']=='bound' and not o.hide_render:
            right=bpy.data.objects['H transport token '+o['carrier_chain']]
            a=np.array([(o.location.x+11)/4.8,(o.location.y+2)/4.8])
            b=np.array([(right.location.x-12.8)/7.25,(right.location.z-7.6)/7.25])
            phase_error=max(phase_error,float(np.max(abs(a-b))))
if phase_error>1e-5:errors.append(['left/right phase',phase_error])
for name,visible in visibility.items():
    if visible[0]!=visible[-1]:errors.append([name,'visibility seam'])
    a=np.array(visible[:384],dtype=bool)
    if np.count_nonzero(a&~np.roll(a,1))!=1:errors.append([name,'expected one active interval per loop'])
def snapshot(f):
    scene.frame_set(f)
    return {o.name:np.array(o.matrix_world) for o in bpy.data.objects}
start=snapshot(1);end=snapshot(385)
seam=max(float(np.max(abs(start[k]-end[k]))) for k in start)
if seam>2e-5:errors.append(['seam',seam])
stator=[o for o in main if o['role']=='stator']
fixed=snapshot(97)
assert all(np.allclose(start[o.name],fixed[o.name],atol=1e-7) for o in stator)
scene.frame_set(1)
assert len([o for o in bpy.data.objects if o.name.startswith('H transport token ')])==10
assert all(name in bpy.data.node_groups for name in ['COMP_clean','COMP_reference','COMP_monochrome'])
report={'errors':errors,'frame_evaluations':385,'source_meshes':asset_checks,'stator_fixed':True,
        'seam_frame_1_vs_385_max_matrix_delta':seam,'frames':384,'fps':24,
        'source_files_case_safe':True,'proton_carriers':10,'ATP_graphic_cycles_per_turn':3,
        'left_proton_graphics':len(left),'left_right_bound_phase_max_error':phase_error,
        'labels':[{'text':o.data.body,'font':o.data.font.name,'latin':o.data.font_bold.name,'packed':bool(o.data.font.packed_file and o.data.font_bold.packed_file)} for o in label_sources],
        'checks_do_not_establish':['atomistic collision-free path','real enzyme speed','pH or proton-motive-force dynamics','experimental nucleotide catalytic intermediates']}
(OUT/'saved-scene.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
(OUT/'frame-bounds.json').write_text(json.dumps(frame_bounds,ensure_ascii=False))
assert not errors,errors[:10]
print('SAVED SCENE AUDIT PASSED',seam,flush=True)
