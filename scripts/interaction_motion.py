"""Source-pose-anchored approach, orientation, binding and release in nm.

The endpoint is supplied by a selected experimental complex. This does not
predict docking poses, affinities, reaction rates or a molecular trajectory.
"""
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation


def unit(v):
    v = np.asarray(v, dtype=float)
    if np.linalg.norm(v) < 1e-9:
        raise ValueError('A nonzero approach direction is required')
    return v/np.linalg.norm(v)


def progress(t, mode):
    t = np.asarray(t, dtype=float)
    if mode == 'loop':
        return (1-np.cos(2*np.pi*t))/2, np.where(t <= .5, 1, -1)
    if mode != 'teaching':
        raise ValueError('mode must be loop or teaching')
    arrive = np.clip((t-.12)/.30, 0, 1)
    leave = np.clip((t-.60)/.30, 0, 1)
    smooth = lambda u: u*u*u*(10+u*(-15+6*u))
    return smooth(arrive)-smooth(leave), np.where(t <= .6, 1, -1)


def matrices(points, direction, distance, fractions, mode='loop', turn=.4, bend=.18):
    direction = unit(direction)
    center = points.mean(0)
    cross_axis = np.eye(3)[np.argmin(abs(direction))]
    tangent = unit(np.cross(direction, cross_axis))
    u, side = progress(fractions, mode)
    # Finish most turning before entering the binding pocket.
    rotation = Rotation.from_rotvec(((1-u)**3)[:, None]*turn*tangent).as_matrix()
    offset = (1-u)[:, None]*direction*distance
    offset += side[:, None]*(np.sin(np.pi*u)*(1-u))[:, None]*tangent*distance*bend
    result = np.repeat(np.eye(4)[None], len(u), axis=0)
    result[:, :3, :3] = rotation
    result[:, :3, 3] = center-(rotation@center)+offset
    return result


def source_contacts(moving, obstacle, threshold):
    tree=cKDTree(obstacle)
    return {(i,j):float(np.linalg.norm(moving[i]-obstacle[j]))
            for i,neighbors in enumerate(tree.query_ball_point(moving,threshold)) for j in neighbors}


def contacts(moving, obstacle, poses, threshold=.12, allowed=None, trees=None):
    if trees is None:trees=[cKDTree(obstacle)] if obstacle.ndim==2 else [cKDTree(p) for p in obstacle]
    counts, minimum = [], []
    for frame,pose in enumerate(poses):
        fixed=obstacle if obstacle.ndim==2 else obstacle[frame]
        tree=trees[0] if len(trees)==1 else trees[frame]
        points = moving@pose[:3, :3].T+pose[:3, 3]
        distances, _ = tree.query(points, k=1)
        if allowed is None:
            count=int(np.count_nonzero(distances < threshold))
        else:
            bad=set()
            for i in np.flatnonzero(distances<threshold):
                for j in tree.query_ball_point(points[i],threshold):
                    distance=float(np.linalg.norm(points[i]-fixed[j]))
                    if (i,j) not in allowed or distance < allowed[i,j]-.01:bad.add(int(i))
            count=len(bad)
        counts.append(count)
        minimum.append(float(distances.min()))
    return np.array(counts), np.array(minimum)


def plan(moving, obstacle, frames=192, mode='loop', distance=None, threshold=.20, obstacle_open=None, closure_start=.8):
    """Choose an illustrative arc, then screen every output pose against atoms."""
    desired = unit(moving.mean(0)-obstacle.mean(0))
    distance = distance or max(3., np.linalg.norm(np.ptp(moving, axis=0))*1.1)
    indices = np.arange(192)
    z = 1-2*(indices+.5)/len(indices)
    angle = indices*2.399963229728653
    directions = np.column_stack([np.sqrt(1-z*z)*np.cos(angle), np.sqrt(1-z*z)*np.sin(angle), z])
    directions = np.vstack([desired, directions[np.argsort(-(directions@desired))]])
    coarse = np.linspace(0, 1, 25)
    fine=np.arange(frames+1)/frames
    def environment(fractions):
        if obstacle_open is None:return obstacle,[cKDTree(obstacle)]
        if obstacle_open.shape!=obstacle.shape:raise ValueError('Open/closed obstacle atoms need identical identities and order')
        u,_=progress(fractions,mode);opening=np.clip((1-u)/(1-closure_start),0,1)
        dynamic=obstacle[None]+opening[:,None,None]*(obstacle_open-obstacle)[None]
        return dynamic,[cKDTree(points) for points in dynamic]
    coarse_fixed,coarse_trees=environment(coarse);fine_fixed,fine_trees=environment(fine)
    baseline, _ = contacts(moving, obstacle, np.eye(4)[None], threshold)
    allowed=source_contacts(moving,obstacle,threshold)
    best = None
    for turn,bend in [(.4,.18),(.2,.08),(0.,0.)]:
        scores=[]
        for direction in directions:
            poses=matrices(moving,direction,distance,coarse,mode,turn,bend)
            count,_=contacts(moving[::3],coarse_fixed,poses,threshold,trees=coarse_trees)
            scores.append(((int(count.max()),int(count.sum()),1-float(direction@desired)),direction))
        # Sparse screening ranks candidates; it never grants the final pass.
        for _,direction in sorted(scores,key=lambda x:x[0])[:10]:
            poses=matrices(moving,direction,distance,fine,mode,turn,bend)
            count,minimum=contacts(moving,fine_fixed,poses,threshold,allowed,fine_trees)
            score=(int(count.max()),int(count.sum()),1-float(direction@desired))
            if best is None or score<best[0]:
                best=(score,direction,poses,count,minimum,turn,bend)
            if count.max()==0:break
        if best[0][0]==0:break
    _,direction,poses,count,minimum,turn,bend=best
    report = {'frames': frames, 'mode': mode, 'units': 'nm', 'direction': direction.tolist(),
              'turn_radians':turn,'bend_fraction':bend,
              'distance_nm': float(distance), 'threshold_nm': threshold,
              'endpoint_source_overlaps': int(baseline[0]),
              'endpoint_source_overlap_pairs':len(allowed),
              'max_new_or_worsened_overlap_atoms': int(count.max()),
              'flagged_frames': (np.flatnonzero(count > 0)+1).tolist(),
              'minimum_atom_center_distance_nm': float(minimum.min()),
              'seam_matrix_error': float(np.abs(poses[-1]-poses[0]).max()),
              'bound_pose_error': float(np.abs(poses[frames//2]-np.eye(4)).max()),
              'collision_screen_passed': bool(not np.any(count > 0)),
              'source_contact_policy':'Only identical source atom pairs are exempt, with at most 0.01 nm additional overlap; equal counts at other pairs cannot pass.',
              'obstacle_opening':None if obstacle_open is None else {'closure_starts_at_progress':closure_start,
                  'interpolation':'Source-guided coordinate interpolation, not molecular dynamics'},
              'limits': 'Atomic-center deep-overlap screening only; no van der Waals, solvent, ligand dynamics, energy or kinetics. Static bound-state conformations remain rigid.'}
    return poses, report
