"""Trace reference-image XY branches into an illustrative 3D ER background.

Depth and tube radii are authored display choices, not microscopy reconstruction.
Requires numpy, Pillow, scipy and scikit-image; writes assets outside the skill.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import distance_transform_edt, gaussian_filter, gaussian_filter1d
from skimage.morphology import remove_small_objects, skeletonize


def trace_paths(skeleton):
    pixels=set(map(tuple,np.argwhere(skeleton)))
    adjacency={}
    for y,x in pixels:
        neighbors=[]
        for dy in [-1,0,1]:
            for dx in [-1,0,1]:
                if not(dx or dy) or (y+dy,x+dx) not in pixels:continue
                if dx and dy and ((y,x+dx) in pixels or (y+dy,x) in pixels):continue
                neighbors.append((y+dy,x+dx))
        adjacency[(y,x)]=neighbors
    vertices={p for p,n in adjacency.items() if len(n)!=2}
    visited=set();paths=[]
    def edge(a,b):return tuple(sorted([a,b]))
    for start in sorted(vertices):
        for nxt in adjacency[start]:
            if edge(start,nxt) in visited:continue
            path=[start];previous=start;current=nxt
            visited.add(edge(previous,current))
            while True:
                path.append(current)
                if current in vertices:break
                following=next(n for n in adjacency[current] if n!=previous)
                visited.add(edge(current,following));previous,current=current,following
            # Keep junction links; discard only short terminal spurs.
            if len(path)>=6 or (len(adjacency[path[0]])>2 and len(adjacency[path[-1]])>2):paths.append(path)
    for start in sorted(pixels):
        for nxt in adjacency[start]:
            if edge(start,nxt) in visited:continue
            path=[start];current=nxt;visited.add(edge(start,nxt))
            while current!=start:
                path.append(current)
                choices=[n for n in adjacency[current] if edge(current,n) not in visited]
                if not choices:break
                following=choices[0];visited.add(edge(current,following));current=following
            if current==start:path.append(start)
            if len(path)>=8:paths.append(path)
    return paths


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--image',type=Path,required=True)
    p.add_argument('--profile',type=Path,default=Path(__file__).resolve().parents[1]/'assets/environment-reference-er.json')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Choose a fresh output directory.')
    profile=json.loads(a.profile.read_text());settings=profile['network']
    if hashlib.sha256(a.image.read_bytes()).hexdigest()!=profile['sha256']:
        raise ValueError('Reference image does not match its source profile.')
    image=Image.open(a.image).convert('RGB')
    x,y,w,h=profile['crop_bottom_left_xywh'];image_height=image.height
    if min(x,y)<0 or min(w,h)<=0 or x+w>image.width or y+h>image.height:
        raise ValueError('Reference crop is outside the image.')
    image=image.crop((x,image_height-y-h,x+w,image_height-y))
    width=settings['trace_width'];height=round(width*h/w)
    rgb=np.asarray(image.resize((width,height),Image.Resampling.LANCZOS),dtype=float)/255
    # Max-channel brightness avoids treating the source depth colors as intensity measurements.
    strength=gaussian_filter(rgb.max(2),.75)
    mask=remove_small_objects(strength>settings['threshold'],min_size=45)
    paths=trace_paths(skeletonize(mask));distance=distance_transform_edt(mask)
    points=[];radii=[];offsets=[0]
    rng=np.random.default_rng(settings['seed'])
    frequencies=rng.normal(size=(14,2))*np.linspace(.65,2,14)[:,None]
    phases=rng.uniform(0,2*np.pi,14)
    for path in paths:
        q=np.asarray(path,dtype=float)
        smooth=gaussian_filter1d(q,1.4,axis=0,mode='nearest');smooth[0]=q[0];smooth[-1]=q[-1]
        selected=np.unique(np.r_[np.arange(0,len(q),3),len(q)-1])
        smooth=smooth[selected];source=q[selected].astype(int)
        xy=np.column_stack((smooth[:,1]/width-.5,.5*height/width-smooth[:,0]/width))
        z=(np.sin(xy @ frequencies.T*2*np.pi+phases)*np.linspace(1,.2,14)).sum(1)/14*settings['depth_amplitude']
        points.extend(np.column_stack((xy,z)).tolist())
        radii.extend((np.clip(distance[source[:,0],source[:,1]],.9,2.4)/width*1.35).tolist())
        offsets.append(len(points))
    if not points:raise ValueError('No usable branch paths; inspect the crop and threshold.')
    a.output.mkdir(parents=True)
    np.savez_compressed(a.output/'network.npz',points=np.asarray(points,np.float32),
                        radii=np.asarray(radii,np.float32),offsets=np.asarray(offsets,np.int32))
    record={'kind':'reference-guided-er-network','source':profile,'paths':len(paths),'points':len(points),
            'geometry_sha256':hashlib.sha256((a.output/'network.npz').read_bytes()).hexdigest(),
            'evidence':'Schematic: XY morphology guided by a thresholded reference image; depth, smoothing and tube thickness are authored. Not tomography or quantitative segmentation.',
            'topology':'Retains image-guided branches and loops after noise/spur pruning. Projected crossings can be fused; this does not establish native ER connectivity.',
            'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (a.output/'manifest.json').write_text(json.dumps(record,indent=2))
    print(json.dumps({'paths':len(paths),'points':len(points),'output':str(a.output)}))


if __name__=='__main__':main()
