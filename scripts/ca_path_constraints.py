"""Continuous geometry checks for piecewise-linear C-alpha coordinates in nm.

This is not a motion solver or an atomic/surface collision test. Rotating rigid
bodies and hinge interpolants need their own checks; their chords are different
paths. Connectivity and limits are supplied by the source-qualified caller.
"""
import numpy as np
from scipy.spatial import cKDTree


def pairs_array(values, count):
    pairs=np.asarray(values)
    if not pairs.size:return np.empty((0,2),dtype=int)
    if (pairs.ndim!=2 or pairs.shape[1]!=2 or not np.issubdtype(pairs.dtype,np.integer)
            or np.any(pairs<0) or np.any(pairs>=count) or np.any(pairs[:,0]==pairs[:,1])):
        raise ValueError('Require pairs of distinct retained point indices')
    pairs=np.sort(pairs,axis=1)
    if len(np.unique(pairs,axis=0))!=len(pairs):raise ValueError('Repeated connectivity pair')
    return pairs


def neighborhoods(edges, count, steps=3):
    """Local exclusions follow modeled connectivity, never author-number arithmetic."""
    neighbors=[set() for _ in range(count)]
    for i,j in edges:neighbors[i].add(j);neighbors[j].add(i)
    result=set()
    for i in range(count):
        seen={i};frontier={i}
        for _ in range(steps):
            frontier={j for node in frontier for j in neighbors[node]}-seen;seen.update(frontier)
        result.update((i,j) for j in seen if i<j)
    return result


def distance_bounds(first, second, pairs):
    """Exact minimum and maximum pair distances along one linear segment."""
    i,j=pairs.T;v=first[i]-first[j];w=second[i]-second[j];delta=w-v
    denom=np.einsum('ij,ij->i',delta,delta)
    t=np.clip(-np.einsum('ij,ij->i',v,delta)/np.maximum(denom,1e-30),0,1)
    low=np.linalg.norm(v+t[:,None]*delta,axis=1)
    high=np.maximum(np.linalg.norm(v,axis=1),np.linalg.norm(w,axis=1))
    return low,high


def inspect(points, edges, policy, observed_endpoints):
    points=np.asarray(points,float)
    if (points.ndim!=3 or points.shape[0]<2 or points.shape[1]<1 or points.shape[2]!=3
            or not np.isfinite(points).all()):raise ValueError('Require finite (samples, points, 3) coordinates')
    count=points.shape[1];edges=pairs_array(edges,count)
    guard=float(policy['nonlocal_guard_nm']);allowance=float(policy['source_contact_allowance_nm'])
    bond_low,bond_high=map(float,policy['adjacent_distance_nm']);step_limit=float(policy['maximum_step_nm'])
    if not (np.isfinite([guard,allowance,bond_low,bond_high,step_limit]).all()
            and 0<=allowance<guard and 0<bond_low<bond_high and step_limit>0):
        raise ValueError('Invalid reviewed CA geometry limits')
    excluded=neighborhoods(edges,count)
    def near(x,radius):
        return np.array([pair for pair in cKDTree(x).query_pairs(radius) if pair not in excluded],int).reshape(-1,2)
    floors={}
    for xyz,mask in observed_endpoints:
        xyz=np.asarray(xyz,float);mask=np.asarray(mask)
        if xyz.shape!=(count,3) or not np.isfinite(xyz).all() or mask.shape!=(count,) or mask.dtype!=bool:
            raise ValueError('Observed endpoint coordinates and boolean coverage must retain point identity')
        for i,j in near(xyz,guard):
            if mask[i] and mask[j]:
                value=float(np.linalg.norm(xyz[i]-xyz[j]))-allowance
                if value<=0:raise ValueError('Source contact cannot authorize overlapping CA centers')
                floors[i,j]=min(floors.get((i,j),guard),value)
    gap_specs=policy.get('unmodeled_gaps',[])
    gaps=pairs_array([g['indices'] for g in gap_specs],count)
    gap_limits=np.array([g['maximum_nm'] for g in gap_specs],float)
    if np.any(~np.isfinite(gap_limits)) or np.any(gap_limits<=0):raise ValueError('Invalid unmodeled-gap limit')
    margin=0.;minimum=float('inf');maximum=0.;step_max=0.;gap_max=np.zeros(len(gaps))
    worst=None
    for index,(first,second) in enumerate(zip(points[:-1],points[1:])):
        travel=float(np.linalg.norm(second-first,axis=1).max());step_max=max(step_max,travel)
        # A pair can close by at most twice the largest point displacement.
        pairs=near(first,guard+2*travel+1e-8)
        if len(pairs):
            low,_=distance_bounds(first,second,pairs)
            slack=low-np.array([floors.get(tuple(pair),guard) for pair in pairs]);k=int(slack.argmin())
            if slack[k]<margin:margin=float(slack[k]);worst={'segment':index,'indices':pairs[k].tolist()}
        if len(edges):
            low,high=distance_bounds(first,second,edges);minimum=min(minimum,float(low.min()));maximum=max(maximum,float(high.max()))
        if len(gaps):gap_max=np.maximum(gap_max,distance_bounds(first,second,gaps)[1])
    plane_min=None
    if 'halfspace' in policy:
        plane=policy['halfspace'];normal=np.asarray(plane['normal'],float);offset=float(plane['offset_nm'])
        if (normal.shape!=(3,) or not np.isfinite(normal).all() or not np.isfinite(offset)
                or not np.isclose(np.linalg.norm(normal),1.,atol=1e-7)):
            raise ValueError('Halfspace needs a finite unit normal and offset')
        plane_min=float((points@normal-offset).min())
    result={'linear_segments':len(points)-1,'points':count,'modeled_edges':len(edges),
            'local_excluded_pairs':len(excluded),'source_contact_pairs':len(floors),
            'nonlocal_CA_margin_nm':margin,'worst_contact':worst,
            'adjacent_CA_min_nm':minimum if len(edges) else None,'adjacent_CA_max_nm':maximum if len(edges) else None,
            'maximum_step_nm':step_max,'gap_max_nm':gap_max.tolist(),'halfspace_min_margin_nm':plane_min,
            'scope':'Exact linear-segment CA distances only; not atomic geometry, chirality, surface clearance or physical dynamics'}
    result['passed']=(margin>=-1e-5 and (not len(edges) or (minimum>=bond_low-1e-5 and maximum<=bond_high+1e-5))
                      and step_max<=step_limit+1e-5 and np.all(gap_max<=gap_limits+1e-5)
                      and (plane_min is None or plane_min>=-1e-5))
    result['passed']=bool(result['passed'])
    return result
