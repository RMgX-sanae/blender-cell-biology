"""Audit the saved fine-scene motor/division rigs at every rendered frame.

This checks evaluated scene behavior, not aesthetic quality or biological kinetics.
"""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Vector
from bpy_extras.object_utils import world_to_camera_view


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kind',choices=['motor','division'],required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    scene=bpy.context.scene
    failures=[]
    records=[]
    if args.kind=='motor':
        heads=[bpy.data.objects['Trailing motor — schematic duplicated head'],
               bpy.data.objects['Leading kinesin motor — K coordinate-derived shape']]
        cargo=bpy.data.objects['Transport assembly']
    else:
        membrane=bpy.data.collections['01 Continuous-volume shape cache — one frame visible at a time']
        chromosomes=bpy.data.collections['02 Four duplicated chromosomes — eight sister chromatids']
    for frame in range(scene.frame_start,scene.frame_end+1):
        scene.frame_set(frame)
        if args.kind=='motor':
            positions=[list(obj.matrix_world.translation) for obj in heads]
            center=list(cargo.matrix_world.translation)
            record={'frame':frame,'head_locations_nm':positions,'cargo_location_nm':center}
            if abs(center[0]-(positions[0][0]+positions[1][0])/2)>1e-5:
                failures.append(f'frame {frame}: cargo disconnected from midpoint rig')
            subjects=heads
        else:
            visible=[obj for obj in membrane.objects if not obj.hide_render]
            if len(visible)!=1:failures.append(f'frame {frame}: {len(visible)} membrane caches visible')
            if len(chromosomes.objects)!=8:failures.append('Chromatid identity count changed')
            record={'frame':frame,'visible_membranes':[obj.name for obj in visible]}
            subjects=list(chromosomes.objects)
        bounds=[]
        for obj in subjects:
            for point in obj.bound_box:
                q=world_to_camera_view(scene,scene.camera,obj.matrix_world@Vector(point))
                bounds.append(tuple(q))
        bounds=np.asarray(bounds)
        if not np.isfinite(bounds).all():failures.append(f'frame {frame}: nonfinite projection')
        if np.any(bounds[:,:2]<0) or np.any(bounds[:,:2]>1) or np.any(bounds[:,2]<=0):
            failures.append(f'frame {frame}: primary subject outside camera')
        record['primary_projected_bounds_xy']=[bounds[:,:2].min(0).tolist(),bounds[:,:2].max(0).tolist()]
        records.append(record)
    if args.kind=='motor':
        midpoint=(scene.frame_end-scene.frame_start)/2
        for item in records:
            phase=(item['frame']-scene.frame_start)/midpoint
            pos=np.asarray(item['head_locations_nm'])
            fixed=pos[1] if phase<=1 else pos[0]
            expected=(0,0,12) if phase<=1 else (8,0,12)
            if not np.allclose(fixed,expected,atol=1e-5):
                failures.append(f'frame {item["frame"]}: stance head slid during swing')
    report={'blend':bpy.data.filepath,'blender':bpy.app.version_string,
            'kind':args.kind,'fps':scene.render.fps,'frames':records,'failures':failures,
            'limits':'Does not certify biological kinetics, visual quality or all environment collisions.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2))
    print(json.dumps({'frames_checked':len(records),'failures':failures}))
    if failures:raise RuntimeError('Animation audit failed; inspect report.')


if __name__=='__main__':main()
