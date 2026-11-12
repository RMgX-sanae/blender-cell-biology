"""Bounded optimization toward source-guided CA goals; not molecular dynamics."""
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from ca_path_constraints import neighborhoods,inspect


def traces(edges,count):
    neighbors=[set() for _ in range(count)]
    for i,j in edges:neighbors[i].add(j);neighbors[j].add(i)
    if any(len(row)>2 for row in neighbors):raise ValueError('Candidate solver supports open unbranched CA traces')
    visited=set();result=[]
    for start in range(count):
        if start in visited or len(neighbors[start])>1:continue
        trace=[];current=start
        while current is not None:
            trace.append(current);visited.add(current);next_nodes=neighbors[current]-visited
            current=next(iter(next_nodes)) if next_nodes else None
        result.append(trace)
    if len(visited)!=count:raise ValueError('Cyclic CA traces need a separately reviewed solver')
    return result


class Constraints:
    def __init__(self,model):
        self.m=model;self.source=model.source;self.target=model.target;self.policy=model.policy
        self.settings=model.config['solver'];self.guard=self.policy['nonlocal_guard_nm']
        s=self.settings;weights=np.asarray(s['local_spring_weights'],float)
        positive=['step_bound_nm','fold_radius_nm','ftol','gtol']
        nonnegative=['fold_length_change_nm','fold_weight','volume_weight','contact_weight','boundary_weight','contact_buffer_nm']
        if (weights.shape!=(3,) or not np.isfinite(weights).all() or np.any(weights<=0)
                or any(not np.isfinite(s[k]) or s[k]<=0 for k in positive)
                or any(not np.isfinite(s[k]) or s[k]<0 for k in nonnegative)
                or any(type(s[k]) is not int or s[k]<1 for k in ('max_iterations','max_line_search'))
                or np.sqrt(3)*s['step_bound_nm']>self.policy['maximum_step_nm']):
            raise ValueError('Invalid bounded solver weights, tolerances or step limit')
        self.excluded=neighborhoods(model.edges,len(self.source));self.floors={}
        for xyz,mask in [(self.source,np.ones(len(self.source),bool)),(self.target,model.mapped)]:
            for i,j in self.pairs(xyz,self.guard):
                if mask[i] and mask[j]:
                    value=np.linalg.norm(xyz[i]-xyz[j])-self.policy['source_contact_allowance_nm']
                    if value<=0:raise ValueError('Source CA centers overlap within the declared allowance')
                    self.floors[i,j]=min(self.floors.get((i,j),self.guard),float(value))
        paths=traces(model.edges,len(self.source));springs=[];quads=[]
        for path in paths:
            springs.extend((min(i,j),max(i,j),step) for at,i in enumerate(path) for step,j in enumerate(path[at+1:at+4],1))
            quads.extend(path[at:at+4] for at in range(len(path)-3))
        springs=np.array(sorted(springs),int).reshape(-1,3);self.ei,self.ej,step=springs.T
        self.weight=np.array(self.settings['local_spring_weights'])[step-1]
        self.length0=np.linalg.norm(self.source[self.ei]-self.source[self.ej],axis=1)
        self.length1=np.linalg.norm(self.target[self.ei]-self.target[self.ej],axis=1)
        pairs=set(map(tuple,self.pairs(self.source,self.settings['fold_radius_nm'])))|set(map(tuple,self.pairs(self.target,self.settings['fold_radius_nm'])))
        pairs=np.array([pair for pair in sorted(pairs) if model.groups[pair[0]]==model.groups[pair[1]]],int).reshape(-1,2)
        ai,aj=pairs.T;el0=np.linalg.norm(self.source[ai]-self.source[aj],axis=1);el1=np.linalg.norm(self.target[ai]-self.target[aj],axis=1)
        keep=abs(el1-el0)<self.settings['fold_length_change_nm'];self.ai,self.aj=ai[keep],aj[keep];self.el0,self.el1=el0[keep],el1[keep]
        self.terminal_edges=model.terminal_mask[self.ei]|model.terminal_mask[self.ej]
        self.quads=np.array(quads,int).reshape(-1,4);self.volume0=self.volumes(self.source);self.volume1=self.volumes(self.target)
        self.terminal_quads=np.any(model.terminal_mask[self.quads],axis=1)

    def volumes(self,points):
        a,b,c,d=points[self.quads].transpose(1,0,2)
        return np.einsum('ij,ij->i',np.cross(b-a,c-b),d-c)

    def pairs(self,points,radius):
        return np.array([pair for pair in sorted(cKDTree(points).query_pairs(radius)) if pair not in self.excluded],int).reshape(-1,2)

    def objective(self,previous,goal,phase,goal_weight):
        s=self.settings;bound=s['step_bound_nm'];pairs=self.pairs(previous,self.guard+2*np.sqrt(3)*bound+.01);ci,cj=pairs.T
        minimum=np.array([self.floors.get((int(i),int(j)),self.guard) for i,j in pairs])+s['contact_buffer_nm']
        vector=previous[ci]-previous[cj];normal=vector/np.linalg.norm(vector,axis=1)[:,None]
        edge_phase=(phase[self.m.groups[self.ei]]+phase[self.m.groups[self.ej]])/2
        lengths=self.length0*(1-edge_phase)+self.length1*edge_phase
        lengths[self.terminal_edges]=np.linalg.norm(goal[self.ei[self.terminal_edges]]-goal[self.ej[self.terminal_edges]],axis=1)
        fold_phase=phase[self.m.groups[self.ai]];elastic=self.el0*(1-fold_phase)+self.el1*fold_phase
        volume_phase=phase[self.m.groups[self.quads]].mean(1);volumes=self.volume0*(1-volume_phase)+self.volume1*volume_phase
        volumes[self.terminal_quads]=self.volumes(goal)[self.terminal_quads]
        def loss(flat):
            points=flat.reshape(-1,3);diff=points-goal;value=goal_weight*float(np.sum(diff*diff));grad=2*goal_weight*diff
            def spring(i,j,desired,weight):
                nonlocal value
                vector=points[i]-points[j];distance=np.linalg.norm(vector,axis=1);delta=distance-desired
                value+=float(np.sum(weight*delta*delta));force=(2*weight*delta/np.maximum(distance,1e-9))[:,None]*vector
                np.add.at(grad,i,force);np.add.at(grad,j,-force)
            spring(self.ei,self.ej,lengths,self.weight);spring(self.ai,self.aj,elastic,s['fold_weight'])
            a,b,c,d=points[self.quads].transpose(1,0,2);v=b-a;w=c-b;x=d-c
            vw=np.cross(v,w);wx=np.cross(w,x);xv=np.cross(x,v);error=np.einsum('ij,ij->i',vw,x)-volumes
            value+=s['volume_weight']*float(error@error)
            for column,derivative in enumerate([-wx,wx-xv,xv-vw,vw]):
                np.add.at(grad,self.quads[:,column],2*s['volume_weight']*error[:,None]*derivative)
            # A bounded step keeps currently nearby pairs in their separating
            # halfspaces. An independent exact linear-segment check follows.
            deficit=np.minimum(np.sum((points[ci]-points[cj])*normal,axis=1)-minimum,0)
            value+=s['contact_weight']*float(deficit@deficit);force=2*s['contact_weight']*deficit[:,None]*normal
            np.add.at(grad,ci,force);np.add.at(grad,cj,-force)
            if 'halfspace' in self.policy:
                plane=self.policy['halfspace'];axis=np.asarray(plane['normal']);below=np.minimum(points@axis-plane['offset_nm'],0)
                value+=s['boundary_weight']*float(below@below);grad+=2*s['boundary_weight']*below[:,None]*axis
            for gap in self.policy.get('unmodeled_gaps',[]):
                i,j=gap['indices'];vector=points[i]-points[j];length=np.linalg.norm(vector);excess=max(0,length-gap['maximum_nm'])
                value+=s['boundary_weight']*excess**2;force=2*s['boundary_weight']*excess*vector/max(length,1e-9)
                grad[i]+=force;grad[j]-=force
            return value,grad.ravel()
        return loss

    def step(self,previous,goal,phase,goal_weight=1.):
        s=self.settings;bound=s['step_bound_nm']
        result=minimize(self.objective(previous,goal,phase,goal_weight),previous.ravel(),jac=True,method='L-BFGS-B',
                        bounds=list(zip((previous-bound).ravel(),(previous+bound).ravel())),
                        options={'ftol':s['ftol'],'gtol':s['gtol'],'maxiter':s['max_iterations'],'maxls':s['max_line_search']})
        points=result.x.reshape(-1,3)
        checks=inspect(np.array([previous,points]),self.m.edges,self.policy,
                       [(self.source,np.ones(len(self.source),bool)),(self.target,self.m.mapped)])
        if not checks['passed']:raise ValueError('Constrained candidate step fails independent geometry: '+str(checks))
        return points,{'margin_nm':checks['nonlocal_CA_margin_nm'],'step_nm':checks['maximum_step_nm'],
                       'rms_to_goal_nm':float(np.sqrt(np.mean(np.sum((points-goal)**2,axis=1)))),
                       'optimizer_converged':bool(result.success),'optimizer_status':int(result.status),'iterations':int(result.nit)}
