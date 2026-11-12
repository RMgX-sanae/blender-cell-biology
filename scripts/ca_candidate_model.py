"""Source-qualified fitting partitions and authored CA targets, not MD domains."""
import gemmi
import numpy as np
from scipy.spatial.transform import Rotation
from prepare_ca_conformation import qualified_endpoints
import ca_internal_coordinates as internal


def fit(first,second):
    a,b=first.mean(0),second.mean(0)
    if len(first)<3 or np.linalg.matrix_rank(first-a)<2 or np.linalg.matrix_rank(second-b)<2:
        raise ValueError('Fitting partition needs three noncollinear matched points')
    u,_,v=np.linalg.svd((first-a).T@(second-b))
    rotation=v.T@np.diag([1,1,np.linalg.det(v.T@u.T)])@u.T
    return rotation,b-rotation@a


class CandidateModel:
    def __init__(self,sources,source_spec,config):
        if config.get('version')!=1 or not config.get('rationale'):
            raise ValueError('Review the candidate assumptions and version before solving')
        self.manifest,data=qualified_endpoints(sources,source_spec)
        self.source=data['first'];self.target=data['second'].copy();self.mapped=data['second_observed']
        self.keys=data['keys'];self.edges=data['edges'];self.policy=self.manifest['geometry'];self.config=config
        lookup={tuple(k):i for i,k in enumerate(self.keys)}
        raw=gemmi.read_structure(str(sources/source_spec['sources']['first']['cif']['file']))[0]
        labels={(c.name,str(r.seqid)):r.label_seq for c in raw for r in c}
        sequence=np.array([labels[tuple(key)] for key in self.keys],int);chains=self.keys[:,0]
        def select(rules):
            mask=np.zeros(len(self.source),bool)
            for rule in rules:
                if rule['chain'] not in chains:raise ValueError('Partition refers to an unselected source chain')
                for start,end in rule['label_ranges']:
                    if type(start) is not int or type(end) is not int or start>end:raise ValueError('Use increasing source label sequence ranges')
                    mask|=(chains==rule['chain'])&(sequence>=start)&(sequence<=end)
            return np.flatnonzero(mask)
        self.names=[row['id'] for row in config['partitions']]
        if len(self.names)!=len(set(self.names)):raise ValueError('Fitting partition IDs must be unique')
        self.groups=np.full(len(self.source),-1,int);self.centers=[];rotations=[];shifts=[];self.fits={}
        fit_config=config['robust_fit']
        if (type(fit_config['iterations']) is not int or fit_config['iterations']<1
                or not 0<fit_config['retained_quantile']<=1 or not np.isfinite(fit_config['distance_floor_nm'])
                or fit_config['distance_floor_nm']<0):raise ValueError('Invalid robust fitting settings')
        for group,row in enumerate(config['partitions']):
            ids=select(row['select'])
            if not len(ids) or np.any(self.groups[ids]>=0):raise ValueError('Partitions must be nonempty and disjoint')
            self.groups[ids]=group;valid=ids[self.mapped[ids]];used=valid.copy()
            for iteration in range(fit_config['iterations']):
                r,t=fit(self.source[used],self.target[used]);errors=np.linalg.norm(self.source[valid]@r.T+t-self.target[valid],axis=1)
                nxt=valid[errors<max(fit_config['distance_floor_nm'],float(np.quantile(errors,fit_config['retained_quantile'])))]
                if len(nxt)<3 or np.array_equal(nxt,used) or iteration==fit_config['iterations']-1:break
                used=nxt
            if row.get('fixed_reference'):r=np.eye(3);t=np.zeros(3)
            center=self.source[ids].mean(0);self.centers.append(center);rotations.append(r);shifts.append(r@center+t-center)
            self.fits[row['id']]={'mapped':len(valid),'fit_used':len(used),'rms_nm':float(np.sqrt(np.mean(errors**2))),
                                  'rotation_degrees':float(np.linalg.norm(Rotation.from_matrix(r).as_rotvec())*180/np.pi)}
        if np.any(self.groups<0):raise ValueError('Every retained CA must belong to an explicit fitting partition')
        self.centers=np.array(self.centers);self.rotations=np.array(rotations);self.shifts=np.array(shifts)
        self.rotvec=Rotation.from_matrix(self.rotations).as_rotvec();assigned=set()
        for row in config['inferred_targets']:
            if not row.get('reason'):raise ValueError('An inferred target needs a reviewed rationale')
            ids=select(row['select']);anchor=lookup[tuple(row['anchor'])];group=self.names.index(row['partition'])
            if (not len(ids) or np.any(self.mapped[ids]) or not self.mapped[anchor] or set(ids)&assigned
                    or np.any(self.groups[ids]!=group) or self.groups[anchor]!=group):
                raise ValueError('Infer only unmatched target points from an observed anchor in the same partition')
            self.target[ids]=(self.source[ids]-self.source[anchor])@self.rotations[group].T+self.target[anchor];assigned.update(ids)
        if assigned!=set(np.flatnonzero(~self.mapped)) or not np.isfinite(self.target).all():
            raise ValueError('All unobserved target positions need explicit inference')
        self.local0=self.source-self.centers[self.groups]
        self.local1=np.einsum('nji,nj->ni',self.rotations[self.groups],self.target-(self.centers+self.shifts)[self.groups])
        self.terminals=[];self.terminal_mask=np.zeros(len(self.source),bool);edge_set={tuple(sorted(edge)) for edge in self.edges}
        for row in config['terminals']:
            start,end=row['label_range'];step=1 if end>start else -1
            positions=list(range(start,end+step,step));ids=[]
            for position in positions:
                found=np.flatnonzero((chains==row['chain'])&(sequence==position))
                if len(found)!=1:raise ValueError('Terminal fragment must retain every source label position')
                ids.append(int(found[0]))
            ids=np.array(ids,int)
            if any(tuple(sorted(edge)) not in edge_set for edge in zip(ids,ids[1:])):
                raise ValueError('Terminal interpolation cannot bridge a missing source region')
            first,second=internal.parameters(self.source[ids]),internal.parameters(self.target[ids])
            for points,values in [(self.source[ids],first),(self.target[ids],second)]:
                if np.max(abs(internal.reconstruct(points[:3],values)-points))>1e-9:raise ValueError('Terminal endpoint reconstruction differs')
            self.terminals.append((ids,first,second,row['phase']));self.terminal_mask[ids[3:]]=True
        self.rates=np.array([r['phase_rate'] for r in config['partitions']],float)
        self.detours=np.array([r.get('detour_nm',[0,0,0]) for r in config['partitions']],float)
        if not np.isfinite(self.rates).all() or np.any(self.rates<1) or self.detours.shape!=(len(self.names),3) or not np.isfinite(self.detours).all():
            raise ValueError('Partition phase rates must reach the endpoint and detours must be finite nm vectors')

    def goal(self,q):
        if not 0<=q<=1:raise ValueError('Candidate progress must be in [0,1]')
        phase=np.minimum(1.,q*self.rates);per_point=phase[self.groups,None]
        rotation=Rotation.from_rotvec(self.rotvec*phase[:,None]).as_matrix()
        local=self.local0*(1-per_point)+self.local1*per_point
        result=np.einsum('nij,nj->ni',rotation[self.groups],local)+self.centers[self.groups]+per_point*self.shifts[self.groups]
        result+=self.detours[self.groups]*np.sin(np.pi*q)
        for ids,first,second,(start,end) in self.terminals:
            if not 0<=start<end<=1:raise ValueError('Terminal phase must increase within [0,1]')
            fraction=np.clip((q-start)/(end-start),0,1);fraction=fraction*fraction*(3-2*fraction)
            result[ids[3:]]=internal.interpolate(result[ids[:3]],first,second,fraction)[3:]
        return result,phase
