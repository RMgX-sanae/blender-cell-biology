"""Illustrative rigid-body release around atomic obstacles; not physical kinetics.

Use when a direct arc fails for a wrapping interface. Local constrained steps
translate and rotate the complete body, then check the actual nonlinear path.
The caller must still screen all simultaneous partners on its final timeline.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation, Slerp

from interaction_motion import unit


def escape(points, obstacles, direction, guard=.24, distance=3., clearance=.8, max_steps=120):
    points = np.asarray(points, float); obstacles = np.asarray(obstacles, float)
    if (points.ndim != 2 or points.shape[1] != 3 or len(points) < 3
            or obstacles.ndim != 2 or obstacles.shape[1] != 3 or not len(obstacles)
            or not np.isfinite(points).all() or not np.isfinite(obstacles).all()):
        raise ValueError('Require finite complete-body and obstacle coordinates in nm')
    if not np.isfinite([guard, distance, clearance]).all() or not 0 < guard < clearance or distance <= 0:
        raise ValueError('Require positive travel and a parking clearance greater than the atomic guard')
    direction = unit(direction); tree = cKDTree(obstacles)
    if tree.query(points)[0].min() < guard:
        raise ValueError('Bound hypothesis overlaps an obstacle; do not exempt source contacts')
    center = points.mean(0); original_center = center.copy(); local = points-center
    rotation = np.eye(3); centers = [center.copy()]; rotations = [rotation.copy()]
    sampled_minimum = float(tree.query(points)[0].min())
    for step in range(max_steps+1):
        current = local@rotation.T+center
        parked = np.linalg.norm(center-original_center) >= distance and tree.query(current)[0].min() >= clearance
        if parked: break
        if step == max_steps: raise ValueError('Rigid release search exhausted; retain the failed trial and revise the path/model')
        pairs = [(i,j) for i, neighbors in enumerate(tree.query_ball_point(current, guard+.41)) for j in neighbors]
        i = np.array([a for a,b in pairs], dtype=int); j = np.array([b for a,b in pairs], dtype=int)
        delta = current[i]-obstacles[j]; lengths = np.linalg.norm(delta, axis=1)
        normal = delta/lengths[:,None]
        jacobian = np.column_stack([normal, np.cross(current[i]-center, normal)])
        constraint = {'type':'ineq', 'fun':lambda x: lengths-(guard+.005)+jacobian@x,
                      'jac':lambda x: jacobian}
        fit = minimize(lambda x: np.sum((x[:3]-.10*direction)**2)+.025*np.sum(x[3:]**2),
            np.zeros(6), method='SLSQP', constraints=constraint if pairs else (),
            bounds=[(-.12,.12)]*3+[(-.15,.15)]*3, options={'ftol':1e-10,'maxiter':100})
        if not fit.success: raise ValueError('Rigid release local step failed: '+fit.message)
        accepted = False
        for scale in (1.,.5,.25,.125,.0625):
            minima = []
            for fraction in np.linspace(0,1,17):
                r = Rotation.from_rotvec(fit.x[3:]*scale*fraction).as_matrix()@rotation
                c = center+fit.x[:3]*scale*fraction
                minima.append(float(tree.query(local@r.T+c)[0].min()))
            if min(minima) >= guard:
                center += fit.x[:3]*scale
                rotation = Rotation.from_rotvec(fit.x[3:]*scale).as_matrix()@rotation
                sampled_minimum = min(sampled_minimum, min(minima)); accepted = True; break
        if not accepted or np.linalg.norm(fit.x[:3]*scale) < .00001:
            raise ValueError('Rigid release stalled; a different route or qualified flexibility is needed')
        centers.append(center.copy()); rotations.append(rotation.copy())
    centers = np.asarray(centers); rotations = np.asarray(rotations)
    angular = Rotation.from_matrix(np.einsum('fij,fkj->fik',rotations[1:],rotations[:-1])).magnitude()
    segments = np.linalg.norm(np.diff(centers,axis=0),axis=1)+np.linalg.norm(local,axis=1).max()*angular
    knots = np.r_[0., np.cumsum(segments)]; knots /= knots[-1]
    path = {'centers':centers, 'rotations':rotations, 'knots':knots}
    report = {'method':'constrained rigid steps with nonlinear substep screening', 'direction':direction.tolist(),
        'waypoints':len(centers), 'guard_nm':guard, 'subsamples_per_step':17,
        'sampled_minimum_nm':sampled_minimum, 'park_clearance_nm':float(tree.query(local@rotation.T+center)[0].min()),
        'center_travel_nm':float(np.linalg.norm(center-original_center)),
        'limits':'Geometric illustrative path; no deformation, energy, affinity, rate or guaranteed global/continuous collision-free solution.'}
    return path, report


def sample(points, path, fractions):
    fractions = np.asarray(fractions, float)
    if fractions.ndim != 1 or not np.isfinite(fractions).all() or np.any((fractions < 0)|(fractions > 1)):
        raise ValueError('Path fractions must lie in [0,1]')
    centers = np.column_stack([np.interp(fractions,path['knots'],path['centers'][:,axis]) for axis in range(3)])
    rotations = Slerp(path['knots'],Rotation.from_matrix(path['rotations']))(fractions).as_matrix()
    return np.einsum('fij,nj->fni',rotations,points-points.mean(0))+centers[:,None,:]
