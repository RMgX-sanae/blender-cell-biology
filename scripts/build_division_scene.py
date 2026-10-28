"""Metaphase-to-cytokinesis teaching animation with a continuous shape sequence.

No mechanical or molecular kinetics claim. Four duplicated chromosomes (eight
sister chromatids) are tracked; a fixed cutaway exposes the schematic interior.
"""
import argparse
import json
import math
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Vector

sys.path.insert(0, str(Path(__file__).resolve().parent))
from detail_scene import Parts, collection, material, mesh, path, save_and_render, studio


def ease(t):
    t = min(1,max(0,t))
    return t*t*(3-2*t)


def chromatid_points(sign, row, segregation):
    separation=.14+4.8*ease((segregation-.13)/.56)
    points=[]
    for u in np.linspace(-1,1,61):
        trail = (.32-.86*ease(segregation/.6))*abs(u)
        points.append((sign*(separation+trail),row+.13*math.sin(u*2.8),u*.88))
    return points


def key_points(curve, coords, frame):
    for point,co in zip(curve.data.splines[0].points,coords):
        point.co=(*co,1)
        point.keyframe_insert(data_path='co',frame=frame)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--style-profile',type=Path,
                   default=Path(__file__).resolve().parent.parent/'assets/style-organelle.json')
    p.add_argument('--engine',choices=['cycles','eevee'],default='cycles')
    p.add_argument('--compositor',choices=['clean','reference','monochrome'])
    p.add_argument('--width',type=int,default=1600)
    p.add_argument('--samples',type=int,default=48)
    p.add_argument('--fps',type=int,default=24)
    p.add_argument('--frame',type=int,default=62)
    p.add_argument('--render',action='store_true')
    p.add_argument('--animate',action='store_true')
    args=p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    if args.output.exists():raise SystemExit('Choose a fresh output directory.')
    source=json.loads((args.asset_dir/'manifest.json').read_text())
    args.frames=source['frames']
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)
    profile=json.loads(args.style_profile.read_text())
    shellmat=material('Cell membrane — schematic cutaway','73A6A3',.5,.008)
    chrom_colors=['B7A1CB','CEAA7E','8FAEBB','BA8F9E']
    chrom_mats=[material(f'Chromosome identity {i+1}',color,.52,.006) for i,color in enumerate(chrom_colors)]
    spindle_mat=material('Spindle microtubules','D5C6A0',.5,0)
    ringmat=material('Contractile ring','D4B58E',.52,0)
    membrane=collection('01 Continuous-volume shape cache — one frame visible at a time')
    for info in source['frame_data']:
        frame=info['frame']
        data=np.load(args.asset_dir/f'frame_{frame:04}.npz',allow_pickle=False)
        obj=mesh(f'Cell membrane at frame {frame:04}',data['vertices'].tolist(),data['faces'].tolist(),[shellmat],membrane)
        thick=obj.modifiers.new('Readable membrane cut edge','SOLIDIFY')
        thick.thickness=.035
        for f,hidden in [(0,True),(frame-1,True),(frame,False),(frame+1,True)]:
            obj.hide_render=hidden
            obj.hide_viewport=hidden
            obj.keyframe_insert(data_path='hide_render',frame=f)
            obj.keyframe_insert(data_path='hide_viewport',frame=f)
    chromosomes=collection('02 Four duplicated chromosomes — eight sister chromatids')
    spindle=collection('03 Bipolar spindle and contractile ring')
    curves=[]
    for index,row in enumerate([-1.22,-.42,.42,1.22]):
        for sign in (-1,1):
            obj=path(f'Chromosome {index+1} sister {sign}',chromatid_points(sign,row,0),.14,chrom_mats[index],chromosomes)
            for j,point in enumerate(obj.data.splines[0].points):
                point.radius=.94+.07*math.sin(j*1.8)
            fibers=[]
            for k in range(4):
                fiber=path(f'Kinetochore microtubule {index} {sign} {k}',[(0,0,0)]*25,.014,spindle_mat,spindle)
                fibers.append(fiber)
            curves.append((sign,row,obj,fibers))
    ring=path('Equatorial contractile ring',[(0,0,0)]*100,.045,ringmat,spindle,cyclic=True)
    schedule=sorted(set([1,args.frames]+list(range(1,args.frames+1,3))))
    for frame in schedule:
        t=(frame-1)/(args.frames-1)
        info=source['frame_data'][frame-1]
        pole=2.9+3.8*ease(t)
        for sign,row,obj,fibers in curves:
            coords=chromatid_points(sign,row,t)
            key_points(obj,coords,frame)
            centroid=Vector(coords[len(coords)//2])
            for k,fiber in enumerate(fibers):
                a=Vector((sign*pole,.0,0))
                b=centroid+Vector((0,(k-1.5)*.025,0))
                points=[]
                for u in np.linspace(0,1,25):
                    pnt=a.lerp(b,float(u))
                    pnt.z+=(k-1.5)*.11*math.sin(math.pi*u)
                    points.append(tuple(pnt))
                key_points(fiber,points,frame)
                fiber.data.bevel_depth=.014*(1-ease((t-.72)/.16))
                fiber.data.keyframe_insert(data_path='bevel_depth',frame=frame)
        radius=info['waist_radius']
        key_points(ring,[(0,radius*math.cos(a),radius*math.sin(a)) for a in np.linspace(0,math.tau,100,endpoint=False)],frame)
        ring.data.bevel_depth=.045*(1-ease((t-.86)/.03))
        ring.data.keyframe_insert(data_path='bevel_depth',frame=frame)
    studio(profile,args,{'position':(2,-20,16),'target':(0,0,0),'scale':21,'light_scale':.5})
    save_and_render(args,{'scene':'Metaphase to cytokinesis — schematic',
                          'cell_shape':source,
                          'duplicated_chromosomes_initial':4,'sister_chromatids_initial':8,
                          'daughter_chromosomes_per_cell':4,
                          'limitations':['Illustrative chromosome count, not a human karyotype.',
                                         'Nuclear-envelope reassembly, chromatin decondensation, midbody and ESCRT abscission details are omitted.',
                                         'Gross cell separation is a smooth volume-conserving geometric sequence, not measured mechanics.',
                                         'Cell membranes have fixed cutaways; filament thickness and timing are chosen for readability.']})


if __name__=='__main__':main()
