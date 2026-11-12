"""Monotone presentation schedules for an already qualified coordinate path.

Retain every source sample. Repeated source anchors insert observation holds;
they do not create a reaction cycle or provide kinetic evidence.
"""
import numpy as np


def presentation(source, profile=None, source_sha256=None):
    count=source['frames']
    if profile is None:
        return {'mode':'source','frames':count,'fps':source['fps'],'captions':'events',
                'source_samples':list(range(count+1)),'sample_frames':list(range(1,count+2))}
    if set(profile)-{'version','mode','source_sha256','fps','frames','anchors','captions'}:
        raise ValueError('Unknown presentation settings')
    if type(profile.get('version')) is not int or profile['version']!=1:
        raise ValueError('Require a version-1 presentation profile')
    if profile.get('mode') not in {'teaching','continuous'}:
        raise ValueError('Choose teaching or continuous; a finite source path is not a closed loop')
    if not source_sha256 or profile.get('source_sha256')!=source_sha256:
        raise ValueError('Presentation profile does not match the qualified source manifest')
    if any(type(profile.get(key)) is not int or profile[key]<1 for key in ('fps','frames')):
        raise ValueError('Require positive integer presentation fps and frames')
    anchors=profile.get('anchors')
    if (not isinstance(anchors,list) or len(anchors)<2 or
            any(not isinstance(pair,list) or len(pair)!=2 or any(type(v) is not int for v in pair) for pair in anchors)):
        raise ValueError('Anchors must be integer [source sample, presentation frame] pairs')
    if anchors[0]!=[0,1] or anchors[-1]!=[count,profile['frames']+1]:
        raise ValueError('Anchors must retain both source endpoints and the complete presentation interval')
    samples=[0];frames=[1.]
    for (a,start),(b,end) in zip(anchors,anchors[1:]):
        if not 0<=a<=b<=count or end<=start:
            raise ValueError('Presentation must move forward in source and frame order')
        if a==b:
            if profile['mode']!='teaching':raise ValueError('Explicit observation holds require teaching mode')
            samples.append(a);frames.append(float(end))
        else:
            samples.extend(range(a+1,b+1))
            frames.extend(np.linspace(start,end,b-a+1)[1:].tolist())
    captions=profile.get('captions','events' if profile['mode']=='teaching' else 'none')
    if captions not in {'events','none'}:raise ValueError('Presentation captions must be events or none')
    return {'mode':profile['mode'],'frames':profile['frames'],'fps':profile['fps'],'captions':captions,
            'source_samples':samples,'sample_frames':frames}


def source_positions(schedule,frames):
    frames=np.asarray(frames,float)
    if not np.isfinite(frames).all() or np.any(frames<1) or np.any(frames>schedule['frames']+1):
        raise ValueError('Query frames must lie inside the presentation interval')
    return np.interp(frames,schedule['sample_frames'],schedule['source_samples'])


def coordinates_at(coordinates,positions):
    positions=np.asarray(positions,float)
    if not np.isfinite(positions).all() or np.any(positions<0) or np.any(positions>len(coordinates)-1):
        raise ValueError('Coordinate queries must stay inside the qualified source path')
    lower=np.floor(positions).astype(int);upper=np.minimum(lower+1,len(coordinates)-1)
    weight=(positions-lower)[:,None,None]
    return coordinates[lower]*(1-weight)+coordinates[upper]*weight


def caption_windows(captions,schedule,source_frames):
    if schedule['captions']=='none':return []
    time=source_positions(schedule,np.arange(1,schedule['frames']+1))/source_frames
    windows=[];coverage=np.zeros(len(time),int)
    for start,end,text in captions:
        active=(time>=start)&(time<end);where=np.flatnonzero(active);coverage+=active
        if not len(where):raise ValueError('Presentation removes every displayed frame of a caption; lengthen it or select captions none')
        windows.append([int(where[0])+1,int(where[-1])+1,text])
    if not np.all(coverage==1):raise ValueError('Event captions must cover each presentation frame exactly once')
    return windows
