"""Virtual CA lengths/angles/dihedrals for source-reviewed terminal fragments."""
import numpy as np


def basis(a,b,c):
    x=c-b
    if np.linalg.norm(x)<1e-9:raise ValueError('Coincident virtual backbone anchors')
    x=x/np.linalg.norm(x);z=np.cross(b-a,x)
    if np.linalg.norm(z)<1e-9:raise ValueError('Collinear virtual backbone anchors')
    z/=np.linalg.norm(z)
    return x,np.cross(z,x),z


def parameters(trace):
    trace=np.asarray(trace,float)
    if trace.ndim!=2 or trace.shape[1]!=3 or len(trace)<4 or not np.isfinite(trace).all():
        raise ValueError('Require at least four finite ordered CA points')
    result=[]
    for i in range(3,len(trace)):
        a,b,c,d=trace[i-3:i+1];x,y,z=basis(a,b,c);v=d-c;length=np.linalg.norm(v)
        if length<1e-9:raise ValueError('Coincident virtual backbone points')
        result.append([length,np.arccos(np.clip(-v@x/length,-1,1)),np.arctan2(v@z,v@y)])
    return np.array(result)


def reconstruct(anchors,values):
    trace=list(np.asarray(anchors,float).copy())
    if len(trace)!=3:raise ValueError('Use three ordered anchors for a virtual backbone')
    for length,angle,dihedral in values:
        x,y,z=basis(*trace[-3:])
        trace.append(trace[-1]+length*(-np.cos(angle)*x+np.sin(angle)*(np.cos(dihedral)*y+np.sin(dihedral)*z)))
    return np.array(trace)


def interpolate(anchors,first,second,fraction):
    delta=second-first;delta[:,2]=(delta[:,2]+np.pi)%(2*np.pi)-np.pi
    return reconstruct(anchors,first+fraction*delta)
