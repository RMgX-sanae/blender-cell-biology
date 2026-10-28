"""Reopen saved chapters, audit every evaluated frame and render unique states.

Repeated frames are exact holds, not interpolated or AI-generated frames.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import bpy
import numpy as np
from mathutils import Vector

sys.path.insert(0,str(Path(__file__).resolve().parent))
from project_context import ROOT, SKILL, FFMPEG, FONT, LATIN_FONT, OPTIONS

def scan(folder):
    scene=bpy.context.scene
    obs=sorted(scene.objects,key=lambda o:o.name)
    molecules=[o for o in obs if o.get('atoms_file')]
    transforms=[];shapes=[];gates=[];visible=[];keys=[];framing=[]
    W=scene.camera.data.ortho_scale
    for f in range(scene.frame_start,scene.frame_end+1):
        scene.frame_set(f);bpy.context.view_layer.update()
        deps=bpy.context.evaluated_depsgraph_get()
        h=hashlib.sha256()
        for o in obs:
            h.update(o.name.encode());h.update(bytes([o.hide_render]))
            if o.hide_render:continue
            h.update(np.round(np.array(o.matrix_world),7).tobytes())
            if o.type=='MESH' and o.data.shape_keys:
                h.update(np.array([k.value for k in o.data.shape_keys.key_blocks],dtype=np.float32).tobytes())
        for node in scene.compositing_node_group.nodes:
            if node.bl_idname=='CompositorNodeAlphaOver':
                h.update(np.array([node.inputs['Factor'].default_value],dtype=np.float32).tobytes())
        keys.append(h.hexdigest())
        transforms.append([np.array(o.matrix_world) for o in molecules])
        shapes.append([o.data.shape_keys.key_blocks[1].value if o.type=='MESH' and o.data.shape_keys else 0 for o in molecules])
        gates.append([o.data.shape_keys.key_blocks[2].value if o.type=='MESH' and o.data.shape_keys and len(o.data.shape_keys.key_blocks)>2 else 0 for o in molecules])
        visible.append([not o.hide_render for o in molecules])
        for o in molecules:
            if o.hide_render:continue
            eo=o.evaluated_get(deps)
            bounds=np.array([o.matrix_world@Vector(p) for p in eo.bound_box])
            # Reserve the fixed title/beat and subtitle/footer bands.
            if max(abs(bounds[:,0]))>W*.49 or max(abs(bounds[:,2]))>W*.178:
                framing.append({'frame':f,'object':o.name,'x':[float(min(bounds[:,0])),float(max(bounds[:,0]))],
                                'z':[float(min(bounds[:,2])),float(max(bounds[:,2]))]})
    matrices=np.array(transforms,dtype=np.float32).reshape(len(keys),len(molecules),4,4)
    determinants=np.linalg.det(matrices[:,:,:3,:3]) if molecules else np.ones(1)
    error=float(np.max(abs(np.matmul(matrices[:,:,:3,:3].transpose(0,1,3,2),matrices[:,:,:3,:3])-np.eye(3)))) if molecules else 0
    np.savez_compressed(folder/'evaluated-transforms.npz',matrices=matrices,shapes=np.array(shapes),gates=np.array(gates),visible=np.array(visible))
    report={'frames':len(keys),'unique_visual_states':len(set(keys)),
            'positive_orientation':bool(np.all(determinants>0)),'max_orthogonality_error':error,
            'finite':bool(np.isfinite(matrices).all()),'framing_warning_frames':len(set(r['frame'] for r in framing)),
            'framing_examples':framing[::max(1,len(framing)//30)][:30],
            'molecules':[{'name':o.name,'asset':o['atoms_file'],'source':o.get('source_structure'),
                          'biological_id':o.get('biological_id')} for o in molecules]}
    (folder/'evaluated-audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    return keys,report

def main():
    p=argparse.ArgumentParser();p.add_argument('--shots',default='all');p.add_argument('--audit-only',action='store_true')
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);chosen=a.shots.split(',')
    storyboard=json.loads((ROOT/'storyboard.json').read_text())
    for c in storyboard['shots']:
        if chosen!=['all'] and c['id'] not in chosen:continue
        folder=ROOT/'shots'/c['id'];blend=folder/'scene.blend'
        bpy.ops.wm.open_mainfile(filepath=str(blend))
        keys,report=scan(folder)
        print('AUDITED',c['id'],report['frames'],report['unique_visual_states'],report['framing_warning_frames'],flush=True)
        if a.audit_only:continue
        cache=folder/'frames';cache.mkdir(exist_ok=True)
        scene=bpy.context.scene;scene.render.image_settings.compression=15
        signature=hashlib.sha256(blend.read_bytes()).hexdigest()
        previous=cache/'scene-sha256.txt'
        # A changed scene uses a new namespace; stale frames never count as hits.
        cache=cache/signature[:12];cache.mkdir(exist_ok=True)
        start=time.monotonic();done=set()
        for f,key in enumerate(keys,1):
            output=cache/(key+'.png')
            if key in done:continue
            if not output.exists():
                scene.frame_set(f);scene.render.filepath=str(output)
                bpy.ops.render.render(write_still=True)
            done.add(key)
            if len(done)%25==0:
                progress={'shot':c['id'],'unique_complete':len(done),'unique_total':len(set(keys)),
                          'elapsed_seconds':round(time.monotonic()-start,1)}
                (ROOT/'render-progress.json').write_text(json.dumps(progress))
                print('RENDER_PROGRESS',json.dumps(progress),flush=True)
        index={'scene_sha256':signature,'fps':scene.render.fps,'frames':len(keys),
               'paths':[str(cache/(key+'.png')) for key in keys],
               'unique_states':len(done),'elapsed_seconds':time.monotonic()-start}
        (folder/'render-index.json').write_text(json.dumps(index,indent=2)+'\n')
        print('RENDER_COMPLETE',c['id'],len(keys),len(done),index['elapsed_seconds'],flush=True)

if __name__=='__main__':main()
