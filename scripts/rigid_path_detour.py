"""Locally translate an intact group along a smooth screened spline detour."""
import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import minimize
from scipy.spatial import cKDTree


def refine(points, obstacles, first, last, maximum=.45, threshold=.20):
    """maximum bounds each control component; report the actual offset norm."""
    frames=np.arange(first,last+1)
    u=(frames-first)/(last-first)
    # Two zero endpoint controls give zero displacement and velocity there.
    knots=np.r_[np.zeros(4),np.linspace(0,1,7)[1:-1],np.ones(4)]
    count=len(knots)-4
    basis=BSpline(knots,np.eye(count),3)(u)[:,2:-2]
    trees=[[cKDTree(obstacle[f]) for f in frames] for obstacle in obstacles]
    target=threshold+.005
    def objective(flat):
        controls=flat.reshape(-1,3);offset=basis@controls
        loss=float(np.sum(controls**2))*.05
        gradient=controls*.1
        differences=np.diff(controls,n=2,axis=0)
        loss+=float(np.sum(differences**2))*.2
        gradient[:-2]+=.4*differences;gradient[1:-1]-=.8*differences;gradient[2:]+=.4*differences
        forces=np.zeros_like(offset)
        for i,f in enumerate(frames):
            mobile=points[f]+offset[i]
            for group,group_trees in zip(obstacles,trees):
                distance,neighbor=group_trees[i].query(mobile,distance_upper_bound=target)
                penetration=np.maximum(target-distance,0)
                loss+=float(penetration@penetration)*10000
                selected=penetration>0
                if np.any(selected):
                    delta=mobile[selected]-group[f,neighbor[selected]]
                    forces[i]-=20000*np.sum(penetration[selected,None]*delta/np.maximum(distance[selected,None],1e-8),axis=0)
        gradient+=basis.T@forces
        return loss,gradient.ravel()
    fit=minimize(objective,np.zeros((count-4)*3),method='L-BFGS-B',jac=True,
        bounds=[(-maximum,maximum)]*((count-4)*3),options={'maxiter':500,'ftol':1e-12,'gtol':1e-7,'maxls':40})
    controls=fit.x.reshape(-1,3);offset=basis@controls
    output=points.copy();output[frames]+=offset[:,None]
    minimum=min(float(tree.query(output[f])[0].min()) for group in trees for f,tree in zip(frames,group))
    return output,{'first_index':int(first),'last_index':int(last),'knots':knots.tolist(),
        'controls_nm':np.vstack([np.zeros((2,3)),controls,np.zeros((2,3))]).tolist(),
        'maximum_translation_nm':float(np.linalg.norm(offset,axis=1).max()),
        'minimum_screened_distance_nm':minimum,'passed':minimum>=threshold,
        'optimizer_success':bool(fit.success),'optimizer_message':str(fit.message),
        'kind':'cubic spline translation; zero endpoint displacement and velocity; not molecular dynamics'}
