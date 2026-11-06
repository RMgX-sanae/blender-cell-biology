"""Coarse CA hinge kinematics and screened escape; not atomistic dynamics.

Joint selections/pivots and angle bounds must be source-reviewed by the caller.
Rigid-segment labels identify pairs whose internal distances cannot change.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation, Slerp


def shape(local, joints, angles):
    points=local.copy()
    for joint,vector in zip(joints,angles.reshape(-1,3)):
        ids=joint['indices'];pivot=joint['pivot']
        points[ids]=(points[ids]-pivot)@Rotation.from_rotvec(vector).as_matrix().T+pivot
    return points


class Contacts:
    def __init__(self, source, obstacle, segments, chains, residues, guard, source_tolerance):
        self.obstacle=obstacle;self.tree=cKDTree(obstacle);self.segments=segments
        self.chains=chains;self.residues=residues;self.guard=guard
        if not 0 <= source_tolerance < guard:raise ValueError('Invalid reviewed source-contact tolerance')
        self.external_floor={(i,j):max(1e-6,float(np.linalg.norm(source[i]-obstacle[j]))-source_tolerance)
            for i,neighbors in enumerate(self.tree.query_ball_point(source,guard)) for j in neighbors}
        self.internal_floor={(int(i),int(j)):max(1e-6,float(np.linalg.norm(source[i]-source[j]))-source_tolerance)
                             for i,j in self.internal_pairs(source,guard)}

    def internal_pairs(self, points, radius):
        pairs=np.array(list(cKDTree(points).query_pairs(radius)),dtype=int).reshape(-1,2)
        if not len(pairs):return pairs
        i,j=pairs.T
        return pairs[(self.segments[i]!=self.segments[j])&
                     ~((self.chains[i]==self.chains[j])&(abs(self.residues[i]-self.residues[j])<=3))]

    def pairs(self, points, radius):
        external=np.array([(i,j) for i,neighbors in enumerate(self.tree.query_ball_point(points,radius))
                           for j in neighbors],dtype=int).reshape(-1,2)
        internal=self.internal_pairs(points,radius)
        floors=np.r_[[self.external_floor.get((int(i),int(j)),self.guard) for i,j in external],
                     [self.internal_floor.get((int(i),int(j)),self.guard) for i,j in internal]]
        return external,internal,floors

    def safe(self, points):
        external,internal,floors=self.pairs(points,self.guard)
        distances=np.r_[np.linalg.norm(points[external[:,0]]-self.obstacle[external[:,1]],axis=1),
                        np.linalg.norm(points[internal[:,0]]-points[internal[:,1]],axis=1)]
        return bool(np.all(distances>=floors-1e-10))


def plan(points, joints, contacts, direction, distance=20., clearance=1., max_steps=350):
    origin=points.mean(0);local=points-origin;center=origin.copy();rotation=np.eye(3)
    # Pivots are in the same origin-relative frame as local coordinates.
    caps=np.array([j['maximum_radians'] for j in joints]);size=3*len(joints);angles=np.zeros(size)
    direction=np.asarray(direction,float);direction/=np.linalg.norm(direction)
    centers=[center.copy()];rotations=[rotation.copy()];values=[angles.copy()]
    if not contacts.safe(points):raise ValueError('Initial source pose fails the declared CA contact policy')
    for step in range(max_steps+1):
        body=shape(local,joints,angles);current=body@rotation.T+center
        if np.linalg.norm(center-origin)>=distance and contacts.tree.query(current)[0].min()>=clearance:break
        if step==max_steps:raise ValueError('Articulated release search exhausted')
        external,internal,floors=contacts.pairs(current,contacts.guard+1.3)
        ie,je=external.T;ii,ji=internal.T
        delta=np.concatenate([current[ie]-contacts.obstacle[je],current[ii]-current[ji]])
        lengths=np.linalg.norm(delta,axis=1);normal=delta/lengths[:,None]
        jacobian=np.zeros((len(lengths),6+size))
        jacobian[:len(external),:6]=np.column_stack([normal[:len(external)],np.cross(current[ie]-center,normal[:len(external)])])
        for k in range(size):
            shifted=angles.copy();shifted[k]+=1e-5
            derivative=(shape(local,joints,shifted)-body)@rotation.T/1e-5
            jacobian[:,6+k]=np.sum(np.concatenate([derivative[ie],derivative[ii]-derivative[ji]])*normal,axis=1)
        fit=minimize(lambda q:np.sum((q[:3]-.2*direction)**2)+.05*np.sum(q[3:6]**2)+
                     .02*np.sum(q[6:]**2)+.0001*np.sum((angles+q[6:])**2),np.zeros(6+size),method='SLSQP',
                     constraints=[{'type':'ineq','fun':lambda q:lengths-floors-.005+jacobian@q,'jac':lambda q:jacobian},
                                  {'type':'ineq','fun':lambda q:caps*caps-np.sum((angles+q[6:]).reshape(-1,3)**2,axis=1)}],
                     bounds=[(-.24,.24)]*3+[(-.06,.06)]*3+[(-.04,.04)]*size,options={'ftol':1e-9,'maxiter':150})
        if not fit.success:raise ValueError('Articulated local step failed: '+fit.message)
        accepted=False
        for scale in (1.,.5,.25,.125,.0625):
            if all(contacts.safe(shape(local,joints,angles+fit.x[6:]*scale*f)@
                   (Rotation.from_rotvec(fit.x[3:6]*scale*f).as_matrix()@rotation).T+center+fit.x[:3]*scale*f)
                   for f in np.linspace(0,1,17)):
                center+=fit.x[:3]*scale;rotation=Rotation.from_rotvec(fit.x[3:6]*scale).as_matrix()@rotation
                angles+=fit.x[6:]*scale;accepted=True;break
        if not accepted or np.linalg.norm(fit.x*scale)<1e-5:raise ValueError('Articulated release stalled')
        centers.append(center.copy());rotations.append(rotation.copy());values.append(angles.copy())
    centers=np.array(centers);rotations=np.array(rotations);values=np.array(values)
    poses=[shape(local,joints,a)@r.T+c for a,r,c in zip(values,rotations,centers)]
    lengths=[float(np.linalg.norm(b-a,axis=1).max()) for a,b in zip(poses,poses[1:])]
    knots=np.r_[0.,np.cumsum(lengths)];knots/=knots[-1]
    return {'centers':centers,'rotations':rotations,'angles':values,'knots':knots}


def sample(path, fractions):
    fractions=np.asarray(fractions,float)
    if not np.isfinite(fractions).all() or np.any((fractions<0)|(fractions>1)):raise ValueError('Path fractions must be in [0,1]')
    center=np.column_stack([np.interp(fractions,path['knots'],path['centers'][:,k]) for k in range(3)])
    angles=np.column_stack([np.interp(fractions,path['knots'],path['angles'][:,k]) for k in range(path['angles'].shape[1])])
    rotation=Slerp(path['knots'],Rotation.from_matrix(path['rotations']))(fractions).as_matrix()
    return center,rotation,angles
