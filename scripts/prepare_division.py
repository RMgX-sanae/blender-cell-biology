"""Continuous-volume Cassini-ovoid cytokinesis mesh sequence (schematic).

The level set changes topology continuously at the waist. This geometric model
does not solve cell mechanics or the molecular mechanism of abscission.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq


def radius_squared(x, a, b):
    return np.maximum(0, np.sqrt(b**4 + 4*a*a*x*x) - x*x - a*a)


def volume(a, b):
    extent = np.sqrt(a*a+b*b)
    x = np.linspace(-extent, extent, 8001)
    return float(np.pi*np.trapezoid(radius_squared(x,a,b),x))


def build(a,b):
    extent = math.sqrt(a*a+b*b)
    inner = math.sqrt(max(0,a*a-b*b))
    intervals = [(-extent,extent)] if a <= b else [(-extent,-inner),(inner,extent)]
    vertices, faces = [], []
    # Cutaway sector: y<0,z>0. This is a viewing strategy, not a biological hole.
    theta = np.linspace(math.pi, 2.5*math.pi, 73)
    for left,right in intervals:
        start = len(vertices)
        x_values = np.linspace(left,right,151)
        for x in x_values:
            radius = math.sqrt(float(radius_squared(x,a,b)))
            vertices.extend([(x,radius*math.cos(t),radius*math.sin(t)) for t in theta])
        for i in range(150):
            for j in range(72):
                k = start+i*73+j
                faces.append((k,k+73,k+74,k+1))
    return np.array(vertices,dtype=np.float32),np.array(faces,dtype=np.int32),len(intervals)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--frames',type=int,default=120)
    args=p.parse_args()
    if args.output.exists():raise SystemExit('Choose a fresh output directory.')
    args.output.mkdir(parents=True)
    target=4*math.pi*3.8**3/3
    frames=[]
    for frame in range(1,args.frames+1):
        t=(frame-1)/(args.frames-1)
        # Prescribe the waist/gap itself so a square-root singularity cannot
        # turn a small parameter increment into a visible one-frame neck jump.
        phase=min(1,t/.88)
        ease=phase*phase*(3-2*phase)
        waist=3.8*(1-ease)
        gap_phase=max(0,(t-.88)/.12)
        gap=2.0*gap_phase*gap_phase*(3-2*gap_phase)
        delta=waist*waist-gap*gap
        lower=math.sqrt(max(0,-delta))
        a=brentq(lambda a:volume(a,math.sqrt(max(0,a*a+delta)))-target,
                 lower,12,xtol=1e-9)
        b=math.sqrt(a*a+delta)
        v,f,components=build(a,b)
        np.savez_compressed(args.output/f'frame_{frame:04}.npz',vertices=v,faces=f)
        frames.append({'frame':frame,'a':a,'b':b,'components':components,
                       'volume':volume(a,b),'extent':math.sqrt(a*a+b*b),
                       'waist_radius':math.sqrt(max(0,b*b-a*a))})
    manifest={'model':'Volume-conserving Cassini-ovoid mesh sequence',
              'frames':args.frames,'target_volume':target,'units':'schematic um',
              'frame_data':frames,
              'limitations':['Shape evolution is authored geometry, not a mechanobiological simulation.',
                             'The topology change illustrates gross separation; ESCRT, midbody and abscission machinery are not modeled.',
                             'A fixed cutaway sector reveals the internal schematic.']}
    (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps({'frames':args.frames,'volume_relative_error':max(abs(f['volume']/target-1) for f in frames),
                      'first_separated_frame':next(f['frame'] for f in frames if f['components']==2)}))


if __name__=='__main__':main()
