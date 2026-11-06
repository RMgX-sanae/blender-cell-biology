"""Authored lumen topology and conserved cargo; not a membrane dynamics model."""
import json
import argparse
import hashlib
from pathlib import Path
import sys
import time
import numpy as np
from scipy.ndimage import label,generate_binary_structure
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from skimage.measure import marching_cubes

sys.path.insert(0,str(Path(__file__).resolve().parent))
from event_timing import phase_progress
from volume_distribution import ball_inside

INTERVALS={'capture':[0,.06],'invagination':[.08,.24],'neck':[.24,.34],
           'uncoat':[.36,.43],'rab':[.43,.50],'transport':[.43,.62],
           'tether':[.62,.69],'dock':[.69,.77],'pore_expand':[.80,.88],
           'mix':[.89,.97],'recycle':[.94,1.]}
SCISSION=.35;FUSION=.80;HALF=.10;TARGET=np.array([3.35,0.,-6.])
ORIGIN=np.array([-5.6,-3.6,-10.]);STOP=np.array([6.,3.6,3.])
CARGO_OFFSETS=np.array([[-.18,0,0],[.10,.10,0],[.10,-.10,0]])
CARGO_RADIUS=.06
CARGO_IDS=['cargo-1','cargo-2','cargo-3']


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def smooth(x):return x*x*(3-2*x)


def state(t):
    q={k:float(v[0]) for k,v in phase_progress([t],INTERVALS).items()}
    inv=smooth(q['invagination']);transport=smooth(q['transport'])
    center=np.array([-1.1*transport+.22*np.sin(np.pi*transport),
                     .18*np.sin(np.pi*transport),.5-1.85*inv-4.65*transport])
    center[0]+=.7*smooth(q['tether'])+.4*smooth(q['dock'])
    radius=inv;neck=min(radius,.75-.55*smooth(q['neck']))
    receptor=center+np.array([0,0,-radius]);receptor[2]=min(receptor[2],0)
    normal=np.array([0.,0.,-1.])
    ligand=receptor-normal*(.28+.6*(1-smooth(q['capture'])))
    cargo=center+CARGO_OFFSETS+np.array([0.,0.,.75*(1-inv)])+(TARGET-center)*smooth(q['mix'])
    return {'t':t,'center':center,'radius':radius,'neck':neck,'cargo':cargo,'cargo_ids':CARGO_IDS.copy(),
            'receptor':receptor,'normal':normal,'ligand':ligand,
            'pore_radius':(.16+.54*smooth(q['pore_expand'])) if t>=FUSION else 0.}


def capsule(x,y,z,start,end,radius):
    direction=end-start;fraction=np.clip(((x-start[0])*direction[0]+(y-start[1])*direction[1]+(z-start[2])*direction[2])/np.dot(direction,direction),0,1)
    return np.sqrt((x-start[0]-fraction*direction[0])**2+(y-start[1]-fraction*direction[1])**2+(z-start[2]-fraction*direction[2])**2)-radius


def field(s,x,y,z):
    c=s['center'];r=s['radius']
    phi=np.minimum(-z,np.sqrt((x-c[0])**2+(y-c[1])**2+(z-c[2])**2)-r)
    phi=np.minimum(phi,np.sqrt((x-TARGET[0])**2+(y-TARGET[1])**2+(z-TARGET[2])**2)-2.)
    if s['t']<SCISSION:
        phi=np.minimum(phi,capsule(x,y,z,c,np.array([0.,0.,1.]),s['neck']))
    if s['t']>=FUSION:
        phi=np.minimum(phi,capsule(x,y,z,c,TARGET,s['pore_radius']))
    return phi


def prepare_grid(spacing):
    axes=[np.arange(lo,hi+spacing*.1,spacing) for lo,hi in zip(ORIGIN,STOP)]
    x=axes[0][None,None,:];y=axes[1][None,:,None];z=axes[2][:,None,None]
    return x,y,z


def inspect(s,phi,spacing):
    if s['cargo_ids']!=CARGO_IDS or np.asarray(s['cargo']).shape!=(len(CARGO_IDS),3):
        raise ValueError('Persistent cargo inventory changed')
    water=phi<-HALF;cytosol=phi>HALF;membrane=np.abs(phi)<=HALF
    positions={'outside':np.array([-3.,0.,1.]),'carrier':s['center'],'target':TARGET}
    expected={('outside','carrier'):s['t']<SCISSION,('carrier','target'):s['t']>=FUSION,('outside','target'):False}
    reports={}
    def at(labels,point):
        point=np.asarray(point,float)
        if point.shape!=(3,) or not np.isfinite(point).all():raise ValueError('Invalid tracked position')
        index=np.rint((point-ORIGIN)/spacing).astype(int)[::-1]
        if np.any(index<0) or np.any(index>=labels.shape):raise ValueError('Tracked point outside the grid domain')
        return int(labels[tuple(index)])
    for connectivity in [1,3]:
        labels,count=label(water,generate_binary_structure(3,connectivity))
        ids={name:at(labels,p) for name,p in positions.items()}
        if not all(ids.values()):raise ValueError('Named lumen seed left the aqueous space')
        for pair,connected in expected.items():
            if (ids[pair[0]]==ids[pair[1]])!=connected:raise ValueError('Wrong lumen connectivity: '+str(pair))
        if count!=(3 if SCISSION<=s['t']<FUSION else 2):raise ValueError('Unexpected disconnected water pocket')
        membrane_count=label(membrane,generate_binary_structure(3,connectivity))[1]
        if membrane_count!=(3 if SCISSION<=s['t']<FUSION else 2):raise ValueError('Membrane components touch or fragment unexpectedly')
        for point in [*s['cargo'],s['ligand']]:
            if at(labels,point)!=ids['carrier']:raise ValueError('Cargo left its permitted connected compartment')
        reports[str(connectivity)]={'water_components':int(count),'membrane_components':int(membrane_count),'seed_labels':ids}
    for index,p in enumerate(s['cargo']):
        if field(s,*p)+CARGO_RADIUS>=-HALF or not ball_inside(water,[spacing]*3,ORIGIN,p,CARGO_RADIUS):
            raise ValueError('Whole cargo crosses membrane: '+str(index))
    # SDF min-unions are 1-Lipschitz; this analytic check covers each cargo ball.
    for name,p,radius,side in [('ligand',s['ligand'],.08,-1),
                               ('cytosolic-tail',s['receptor']+s['normal']*.28,.05,1)]:
        value=field(s,*p)
        if side==-1 and value+radius>=-HALF:raise ValueError('Lumen/extracellular domain crossed membrane')
        if side==1 and value-radius<=HALF:raise ValueError('Cytosolic tail crossed membrane')
    if abs(field(s,*s['receptor']))>1e-9:raise ValueError('Receptor anchor left the membrane mid-surface')
    for i,p in enumerate(s['cargo']):
        if np.any(np.linalg.norm(s['cargo'][:i]-p,axis=1)<=2*CARGO_RADIUS):raise ValueError('Cargo overlaps another cargo')
    cyto_components=label(cytosol,generate_binary_structure(3,1))[1]
    if cyto_components!=1:raise ValueError('Unexpected cytosol disconnection')
    return {'time':s['t'],'connectivity':reports,'cargo_count':len(s['cargo']),
            'cargo_centers':s['cargo'].tolist(),'carrier_center':s['center'].tolist(),
            'membrane_voxels':int(membrane.sum()),'cytosol_components':int(cyto_components)}


def inspect_mesh(s,vertices,faces,level,spacing):
    """Check exported surfaces independently of the voxel component labels."""
    if not np.isfinite(vertices).all() or not len(faces):raise ValueError('Invalid surface mesh')
    oriented=np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]])
    edges,inverse,counts=np.unique(np.sort(oriented,axis=1),axis=0,return_inverse=True,return_counts=True)
    if np.any(counts>2):raise ValueError('Non-manifold mesh edge')
    sums=np.bincount(inverse,weights=np.where(oriented[:,0]<oriented[:,1],1,-1))
    if np.any(sums[counts==2]!=0):raise ValueError('Inconsistent adjacent face winding')
    graph=coo_matrix((np.ones(len(edges)),(edges[:,0],edges[:,1])),shape=(len(vertices),len(vertices))).tocsr()
    number,labels=connected_components(graph,directed=False)
    expected=3 if SCISSION<=s['t']<FUSION else 2
    if number!=expected:raise ValueError('Wrong mesh component count')
    euler=np.bincount(labels,minlength=number)-np.bincount(labels[edges[:,0]],minlength=number)+np.bincount(labels[faces[:,0]],minlength=number)
    if sorted(euler)!=[1]+[2]*(expected-1):raise ValueError('Wrong component Euler characteristics')
    ends=vertices[edges[counts==1]];grid=prepare_grid(spacing)
    stop=np.array([grid[0].max(),grid[1].max(),grid[2].max()])
    on_crop=np.any(np.all(np.isclose(ends,ORIGIN,atol=1e-5),axis=1)|np.all(np.isclose(ends,stop,atol=1e-5),axis=1),axis=1)
    if not np.all(on_crop):raise ValueError('Open mesh edge away from the crop boundary')
    error=float(np.max(np.abs(field(s,*vertices.T)-level)))
    if error>spacing/2+1e-5:raise ValueError('Surface departs from the declared field level')
    centers=vertices[faces].mean(1)
    normals=np.cross(vertices[faces[:,1]]-vertices[faces[:,0]],vertices[faces[:,2]]-vertices[faces[:,0]])
    if np.any(np.linalg.norm(normals,axis=1)==0):raise ValueError('Zero-area face')
    gradient=np.column_stack([(field(s,*(centers+np.eye(3)[axis]*1e-4).T)-field(s,*(centers-np.eye(3)[axis]*1e-4).T))/(2e-4) for axis in range(3)])
    dot=np.einsum('ij,ij->i',normals,gradient)*np.sign(level)
    if np.any(np.bincount(labels[faces[:,0]],weights=dot,minlength=number)<=0):raise ValueError('Reversed component orientation')
    return {'components':int(number),'vertices':len(vertices),'triangles':len(faces),
            'euler_by_component':euler.tolist(),'boundary_edges_only_on_crop':True,
            'max_field_level_error':error,'component_orientation_checked':True,
            'positive_point_gradient_fraction':float(np.mean(dot>0))}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    output=p.parse_args().output;output.mkdir(parents=True,exist_ok=False)
    times=np.unique(np.r_[np.linspace(0,1,33),.24,.34,.349,.35,.36,.43,.62,.69,.77,.799,.8,.801,.88,.89,.97])
    keys={0.,.24,.34,.35,.43,.62,.77,.8,.88,.97,1.};reports=[];snapshots=[];started=time.monotonic()
    for spacing,subset in [(.08,times),(.04,[.34,.35,.77,.799,.8,.801,.97])]:
        grid=prepare_grid(spacing)
        for t in subset:
            s=state(float(t));phi=field(s,*grid).astype(np.float32)
            report=inspect(s,phi,spacing);report['spacing']=spacing;reports.append(report)
            if spacing==.08 and t in keys:
                meshes={};checks={}
                for name,level in [('lumen_surface',-HALF),('cytosolic_surface',HALF)]:
                    vertices,faces,_,_=marching_cubes(phi,level=level,spacing=(spacing,)*3,allow_degenerate=False)
                    meshes[name+'_vertices']=(vertices[:,::-1]+ORIGIN).astype(np.float32)
                    meshes[name+'_faces']=(faces if name=='lumen_surface' else faces[:,::-1]).astype(np.int32)
                    checks[name]=inspect_mesh(s,meshes[name+'_vertices'],meshes[name+'_faces'],level,spacing)
                path=output/f'state_{t:.3f}.npz'
                np.savez_compressed(path,**meshes,cargo=s['cargo'],cargo_ids=np.array(CARGO_IDS),receptor=s['receptor'],ligand=s['ligand'],normal=s['normal'])
                snapshots.append({'file':path.name,'sha256':digest(path),'time':float(t),'spacing':spacing,'mesh_checks':checks})
            print(json.dumps({'t':float(t),'spacing':spacing,'passed':True}),flush=True)
    result={'passed':True,'states':len(reports),'seconds':time.monotonic()-started,'reports':reports,
            'units':'dimensionless authored schematic units','cargo_ids':CARGO_IDS,
            'limits':'Geometric lumen connectivity and cargo/domain side checks. No proteins, lipid mass conservation, hemifusion leaflets, experimental kinetics or complete mechanism implemented.'}
    (output/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    scope=Path(__file__).resolve().parents[1]/'assets/vesicle-uptake-scope.json'
    (output/'source-scope.json').write_bytes(scope.read_bytes())
    manifest={'kind':'vesicle-topology-keys','version':1,'units':result['units'],
              'preparation_sha256':digest(__file__),'source_scope_sha256':digest(output/'source-scope.json'),
              'audit_sha256':digest(output/'audit.json'),'cargo_ids':CARGO_IDS,
              'cargo_radius':CARGO_RADIUS,'snapshots':snapshots,'limits':result['limits']}
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['passed','states','seconds']}))


if __name__=='__main__':main()
