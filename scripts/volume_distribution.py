"""Conservative rigid-body placement inside a declared binary voxel region.

Arrays use ZYX; coordinates, spacing and origin use XYZ in one caller-declared
length unit. Origin is voxel (0,0,0)'s center. True voxels denote available
closed boxes; exterior space is forbidden. Placements are schematic, not a
measured distribution or equilibrium simulation.
"""
import numpy as np
from scipy.ndimage import distance_transform_edt
from scipy.spatial.transform import Rotation


def grid(mask, spacing_xyz, origin_xyz):
    mask=np.asarray(mask);spacing=np.asarray(spacing_xyz,float);origin=np.asarray(origin_xyz,float)
    if mask.dtype!=bool or mask.ndim!=3 or not mask.size or not mask.any():
        raise ValueError('Supply a nonempty 3D boolean available-region mask')
    if spacing.shape!=(3,) or not np.isfinite(spacing).all() or np.any(spacing<=0):
        raise ValueError('Supply positive finite XYZ voxel spacing')
    if origin.shape!=(3,) or not np.isfinite(origin).all():raise ValueError('Supply finite XYZ voxel-center origin')
    return mask,spacing,origin


def enclosing_sphere(points):
    points=np.asarray(points,float)
    if points.ndim!=2 or points.shape[1:]!=(3,) or len(points)<2 or not np.isfinite(points).all():
        raise ValueError('Supply finite XYZ body vertices in the same units as the mask')
    center=(points.min(0)+points.max(0))/2
    radius=float(np.linalg.norm(points-center,axis=1).max())
    if radius<=0:raise ValueError('Body must have positive spatial extent')
    return center,radius


def sample(mask, spacing_xyz, origin_xyz, points, count, *, seed=23,
           boundary_gap=0., pair_gap=0., max_attempts=100000):
    """Place copies of one body; all of its enclosing ball must fit.

    EDT measures distances between voxel centers. Subtracting the voxel half
    diagonal and proposal jitter norm gives a lower bound on distance to any
    forbidden voxel box. This can reject feasible placements; it never claims
    exact continuous clearance. Pair exclusion uses complete enclosing balls.
    """
    mask,spacing,origin=grid(mask,spacing_xyz,origin_xyz);local_center,radius=enclosing_sphere(points)
    if type(count)!=int or count<1 or type(max_attempts)!=int or max_attempts<1:
        raise ValueError('Count and attempt limit must be positive integers')
    if not np.isfinite([boundary_gap,pair_gap]).all() or min(boundary_gap,pair_gap)<0:
        raise ValueError('Gaps must be finite and nonnegative')
    distance=distance_transform_edt(np.pad(mask,1),sampling=spacing[::-1])[1:-1,1:-1,1:-1]
    lower=distance-np.linalg.norm(spacing)/2
    tolerance=1e-12*max(1.,np.linalg.norm(np.asarray(mask.shape)[::-1]*spacing))
    sites=np.argwhere(mask & (lower>radius+boundary_gap+tolerance))
    if not len(sites):raise ValueError('No conservative whole-body placement sites')
    rng=np.random.default_rng(seed);centers=[];transforms=[];bounds=[]
    for attempt in range(1,max_attempts+1):
        index=sites[rng.integers(len(sites))]
        jitter=rng.uniform(-.5,.5,3)*spacing
        bound=float(lower[tuple(index)]-np.linalg.norm(jitter))
        if bound<=radius+boundary_gap+tolerance:continue
        center=origin+index[::-1]*spacing+jitter
        if centers and np.any(np.linalg.norm(np.asarray(centers)-center,axis=1)<=2*radius+pair_gap+tolerance):continue
        rotation=Rotation.random(random_state=rng).as_matrix()
        matrix=np.eye(4);matrix[:3,:3]=rotation;matrix[:3,3]=center-rotation@local_center
        centers.append(center);transforms.append(matrix);bounds.append(bound)
        if len(centers)==count:
            return {'transforms':np.asarray(transforms),'centers_xyz':np.asarray(centers),
                    'body_center_xyz':local_center,'body_radius':radius,
                    'clearance_lower_bounds':np.asarray(bounds),'attempts':attempt}
    raise ValueError(f'Packing not reached: {len(centers)}/{count}; change the declared count or geometry explicitly')


def ball_inside(mask, spacing_xyz, origin_xyz, center_xyz, radius):
    """Independent box-distance check, without EDT or center-only assumptions.

    Tests the ball against the domain boundary and every forbidden voxel box
    in its bounding neighborhood. Touching a forbidden box is rejected.
    """
    mask,spacing,origin=grid(mask,spacing_xyz,origin_xyz)
    center=np.asarray(center_xyz,float);shape=np.asarray(mask.shape)[::-1]
    if center.shape!=(3,) or not np.isfinite(center).all() or not np.isfinite(radius) or radius<0:
        raise ValueError('Supply a finite center and nonnegative radius')
    if np.any(center-radius<=origin-spacing/2) or np.any(center+radius>=origin+(shape-.5)*spacing):return False
    low=np.maximum(0,np.floor((center-radius-origin)/spacing-.5).astype(int))
    high=np.minimum(shape-1,np.ceil((center+radius-origin)/spacing+.5).astype(int))
    patch=mask[low[2]:high[2]+1,low[1]:high[1]+1,low[0]:high[0]+1]
    forbidden=np.argwhere(~patch)[:,::-1]+low
    if not len(forbidden):return True
    nearest=np.maximum(np.abs(center-(origin+forbidden*spacing))-spacing/2,0)
    return bool(np.all(np.linalg.norm(nearest,axis=1)>radius))


def validate(mask, spacing_xyz, origin_xyz, points, transforms, *, boundary_gap=0., pair_gap=0.):
    """Check saved rigid transforms with independent voxel-box containment."""
    mask,spacing,origin=grid(mask,spacing_xyz,origin_xyz);local_center,radius=enclosing_sphere(points)
    transforms=np.asarray(transforms,float)
    if transforms.ndim!=3 or transforms.shape[1:]!=(4,4) or not len(transforms) or not np.isfinite(transforms).all():
        raise ValueError('Supply finite rigid 4x4 transforms')
    if not np.isfinite([boundary_gap,pair_gap]).all() or min(boundary_gap,pair_gap)<0:raise ValueError('Invalid gaps')
    centers=[]
    for matrix in transforms:
        rotation=matrix[:3,:3]
        if (not np.allclose(matrix[3],[0,0,0,1],atol=1e-10,rtol=0) or
            not np.allclose(rotation.T@rotation,np.eye(3),atol=1e-10,rtol=0) or
            not np.isclose(np.linalg.det(rotation),1,atol=1e-10,rtol=0)):
            raise ValueError('Body transform stretches or reflects source geometry')
        center=rotation@local_center+matrix[:3,3]
        if not ball_inside(mask,spacing,origin,center,radius+boundary_gap):
            raise ValueError('Whole-body bound intersects forbidden space')
        if centers and np.any(np.linalg.norm(np.asarray(centers)-center,axis=1)<=2*radius+pair_gap):
            raise ValueError('Whole-body bounds violate pair separation')
        centers.append(center)
    return {'passed':True,'count':len(transforms),'body_radius':radius,
            'independent_voxel_box_checks':len(transforms),
            'limits':'Containment of enclosing balls in the declared voxel union; no source segmentation accuracy, actual organelle abundance or motion/path validation.'}
