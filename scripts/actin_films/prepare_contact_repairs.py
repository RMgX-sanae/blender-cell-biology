"""Fit intact CapZ placements and screen local teaching-path corrections.

These are bounded geometric adjustments to cross-structure illustrations, not
experimental interfaces, energy minimization, or molecular dynamics.
"""
import hashlib
import json
import sys

import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from project_context import ROOT, SKILL
sys.path.insert(0, str(SKILL / 'scripts'))
from rigid_path_detour import refine


def transform(points, matrix):
    return points @ matrix[:3, :3].T + matrix[:3, 3]


def fit_placement(points, obstacles):
    center = points.mean(0)
    centered = points - center
    tree = cKDTree(obstacles)
    target = .225

    def objective(q):
        rotation = Rotation.from_rotvec(q[3:]).as_matrix()
        posed = centered @ rotation.T + center + q[:3]
        distance, neighbor = tree.query(posed, distance_upper_bound=target)
        selected = distance < target
        penetration = target - distance[selected]
        force = -20000 * penetration[:, None] * (
            posed[selected] - obstacles[neighbor[selected]]) / np.maximum(distance[selected, None], 1e-9)
        gradient = np.r_[.06 * q[:3], 2 * q[3:]]
        gradient[:3] += force.sum(0)
        for axis in range(3):
            h = np.zeros(3);h[axis] = 1e-6
            derivative = (Rotation.from_rotvec(q[3:] + h).as_matrix()
                          - Rotation.from_rotvec(q[3:] - h).as_matrix()) / 2e-6
            gradient[3 + axis] += np.sum(force * (centered[selected] @ derivative.T))
        loss = .03 * (q[:3] @ q[:3]) + q[3:] @ q[3:] + 10000 * (penetration @ penetration)
        return loss, gradient

    fit = minimize(objective, np.zeros(6), method='L-BFGS-B', jac=True,
        bounds=[(-.4, .4)] * 3 + [(-.09, .09)] * 3,
        options={'maxiter':350, 'ftol':1e-13, 'gtol':1e-7, 'maxls':40})
    correction = np.eye(4)
    correction[:3, :3] = Rotation.from_rotvec(fit.x[3:]).as_matrix()
    correction[:3, 3] = center + fit.x[:3] - correction[:3, :3] @ center
    posed = transform(points, correction)
    report = {'correction':correction.tolist(),
        'before_minimum_nm':float(tree.query(points)[0].min()),
        'after_minimum_nm':float(tree.query(posed)[0].min()),
        'maximum_atom_displacement_nm':float(np.linalg.norm(posed - points, axis=1).max()),
        'rotation_degrees':float(np.linalg.norm(fit.x[3:]) * 180 / np.pi),
        'optimizer_success':bool(fit.success), 'optimizer_message':str(fit.message)}
    if report['after_minimum_nm'] < .21 or report['maximum_atom_displacement_nm'] > .5:
        raise ValueError('CapZ placement requires further source/interface review: ' + json.dumps(report))
    return correction, report


def repair_path(points, obstacles, offset, far_shift=(0, 0, 0)):
    # Parameterize by travel fraction, independent of film fps and duration.
    fraction = np.linspace(0, 1, 241)
    offset = np.array(offset, dtype=float)
    side = np.cross(offset, [0, 0, 1]);side *= 1.2 / np.linalg.norm(side)
    far = np.maximum(0, (fraction - .5) * 2)
    eased = far**3 * (10 + far * (-15 + 6 * far))
    travel = (fraction[:, None] * offset + np.sin(np.pi * far[:, None])**2 * side
              + eased[:, None] * np.array(far_shift))
    moving = points[None] + travel[:, None]
    fixed = np.broadcast_to(obstacles, (len(fraction), len(obstacles), 3))
    tree = cKDTree(obstacles)
    minima = np.array([tree.query(p)[0].min() for p in moving])
    result = {'offset':offset.tolist(), 'far_shift':list(far_shift),
              'parameter':'travel_fraction', 'sample_count':len(fraction),
              'before_minimum_nm':float(minima.min()), 'detour':None}
    bad = np.flatnonzero(minima < .215)
    if len(bad):
        first = max(0, int(bad.min()) - 8)
        last = min(len(fraction) - 1, int(bad.max()) + 8)
        moving, detour = refine(moving, [fixed], first, last, maximum=.65, threshold=.21)
        detour['first_fraction'] = float(fraction[first])
        detour['last_fraction'] = float(fraction[last])
        result['detour'] = detour
    result['after_minimum_nm'] = float(min(tree.query(p)[0].min() for p in moving))
    if result['after_minimum_nm'] < .20:
        raise ValueError('Teaching path still has overlapping atomic centers: ' + json.dumps(result))
    return result


def repair_rigid_release(points, obstacles, pivot, offset, far_shift=(0,0,0)):
    """Withdraw an intact daughter filament or FH2 dimer with a rigid turn."""
    centered = points-pivot
    tree = cKDTree(obstacles)
    offset = np.array(offset, dtype=float)
    fraction = np.linspace(0, .25, 61)
    travel = fraction[:, None]*offset
    knots = np.r_[np.zeros(4), np.linspace(0,1,7)[1:-1], np.ones(4)]
    basis = BSpline(knots, np.eye(len(knots)-4), 3)(fraction/.25)[:,2:-2]
    count = basis.shape[1];target = .218

    def objective(flat):
        controls = flat.reshape(count,6);poses = basis@controls
        gradient = np.zeros_like(poses)
        loss = .03*np.sum(controls[:,:3]**2)+np.sum(controls[:,3:]**2)
        for i,q in enumerate(poses):
            rotation = Rotation.from_rotvec(q[3:]).as_matrix()
            moved = centered@rotation.T+pivot+travel[i]+q[:3]
            distance,neighbor = tree.query(moved,distance_upper_bound=target)
            selected = distance<target;penetration = target-distance[selected]
            loss += 10000*(penetration@penetration)
            force = -20000*penetration[:,None]*(moved[selected]-obstacles[neighbor[selected]])/np.maximum(distance[selected,None],1e-9)
            gradient[i,:3] = force.sum(0)
            for axis in range(3):
                h = np.zeros(3);h[axis] = 1e-6
                derivative = (Rotation.from_rotvec(q[3:]+h).as_matrix()-Rotation.from_rotvec(q[3:]-h).as_matrix())/2e-6
                gradient[i,axis+3] = np.sum(force*(centered[selected]@derivative.T))
        result = basis.T@gradient
        result[:,:3] += .06*controls[:,:3];result[:,3:] += 2*controls[:,3:]
        return loss,result.ravel()

    fit = minimize(objective,np.zeros(count*6),jac=True,method='L-BFGS-B',
        bounds=([(-1.,1.)]*3+[(-.3,.3)]*3)*count,
        options={'maxiter':500,'ftol':1e-12,'gtol':1e-7,'maxls':40})
    controls = np.vstack([np.zeros((2,6)),fit.x.reshape(count,6),np.zeros((2,6))])
    spline = BSpline(knots,controls,3)
    side = np.cross(offset,[0,0,1]);side *= 1.2/np.linalg.norm(side)
    minima = [];angles = [];shifts = []
    for value in np.linspace(0,1,1001):
        q = spline(value/.25) if value<.25 else np.zeros(6)
        far = max(0,(value-.5)*2)
        eased = far**3*(10+far*(-15+6*far))
        moved = centered@Rotation.from_rotvec(q[3:]).as_matrix().T+pivot+offset*value+side*np.sin(np.pi*far)**2+q[:3]+np.array(far_shift)*eased
        minima.append(float(tree.query(moved)[0].min()))
        angles.append(float(np.linalg.norm(q[3:])*180/np.pi));shifts.append(float(np.linalg.norm(q[:3])))
    result = {'offset':offset.tolist(),'far_shift':list(far_shift), 'parameter':'travel_fraction',
        'sample_count':1001,'after_minimum_nm':min(minima),
        'detour':{'kind':'rigid pose spline about the recorded pivot; authored teaching release',
            'controls_nm_radians':controls.tolist(),'pivot_nm':pivot.tolist(),'knots':knots.tolist(),
            'first_fraction':0.,'last_fraction':.25,'maximum_rotation_degrees':max(angles),
            'maximum_translation_nm':max(shifts),'optimizer_success':bool(fit.success),
            'optimizer_message':str(fit.message)}}
    if min(minima)<.20:raise ValueError('Rigid release needs further review: '+json.dumps(result))
    return result


def main():
    assets = ROOT / 'assets'
    def atoms(code, chains):
        return np.concatenate([np.load(assets/code/(c+'.npz'))['atoms'] for c in chains]).astype(float)
    plans = json.loads((assets/'approach-plans.json').read_text())
    cap_fit = np.array(json.loads((assets/'cap-fit.json').read_text())['matrix'])
    branch_fit = np.array(json.loads((assets/'actin-target-fits.json').read_text())['branch_K'])
    cap = np.concatenate([np.load(assets/'filaments/cap'/f'chain_{c}.npz')['atoms'] for c in 'GH'])
    core = atoms('9B27', 'ABCDEF')
    branch = atoms('7TPT', 'ABCDEFGHIJKNOPQRS')
    output = {'schema_version':1, 'threshold_nm':.20, 'placements':{}, 'motions':{},
        'limits':'Intact rigid CapZ adjustments and whole-group spline detours. Static source fits and authored paths are not validated docking, conformational dynamics or reaction kinetics.',
        'sources_sha256':{name:hashlib.sha256((ROOT/'sources'/(name+'.cif')).read_bytes()).hexdigest()
                          for name in ['8F8Q', '9B27', '7TPT', '6FM2-assembly1']}}
    for site, initial, fixed in [('formin', cap_fit, core), ('branch', branch_fit @ cap_fit, branch)]:
        correction, record = fit_placement(transform(cap, initial), fixed)
        record['matrix'] = (correction @ initial).tolist()
        record['source'] = '8F8Q G/H; original heterodimer interface retained'
        output['placements'][site] = record
        offset = np.array(plans['cap']['offset'])
        if site == 'branch':offset = branch_fit[:3, :3] @ offset * .6
        output['motions']['cap_'+site] = repair_path(transform(cap, correction @ initial), fixed, offset)
        print('CAP_CONTACT_REPAIR', site, output['motions']['cap_'+site]['after_minimum_nm'], flush=True)
    output['motions']['debranch'] = repair_rigid_release(atoms('7TPT','HIJK'), atoms('7TPT','ABCDEFGNOPQRS'),
                                                   atoms('7TPT','H').mean(0),plans['debranch']['offset'])
    print('PATH_CONTACT_REPAIR', 'debranch', output['motions']['debranch']['after_minimum_nm'], flush=True)
    jobs = [
        ('cap_recycle', atoms('6FM2-assembly1', ['A']), atoms('6FM2-assembly1', ['B','A-2','B-2']), np.array(plans['cap_recycle']['offset'])*.6, (0,0,0))]
    for name, moving, fixed, offset, far_shift in jobs:
        output['motions'][name] = repair_path(moving, fixed, offset, far_shift)
        print('PATH_CONTACT_REPAIR', name, output['motions'][name]['after_minimum_nm'], flush=True)
    heads = atoms('9B27','GH')
    output['motions']['formin_release'] = repair_rigid_release(heads,core,heads.mean(0),plans['formin_release']['offset'],(12,16,0))
    print('PATH_CONTACT_REPAIR', 'formin_release', output['motions']['formin_release']['after_minimum_nm'], flush=True)
    (assets/'contact-repairs.json').write_text(json.dumps(output, indent=2)+'\n')


if __name__ == '__main__':main()
