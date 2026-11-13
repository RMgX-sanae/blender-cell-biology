"""Independent raw-source, CA contact and saved articulated-scene checks."""
import argparse
import hashlib
from itertools import product
import json
from pathlib import Path

import gemmi
import numpy as np
from molecular_asset_io import chain_file
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation


def require(condition,message):
    if not condition:raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_rows(path):
    rows={}
    for chain in gemmi.read_structure(str(path))[0]:
        rows[chain.name]=[(str(r.seqid),r.name,np.array(list(a.pos))*.1)
                         for r in chain for a in r if a.name=='CA' and a.occ>0 and a.altloc in ('\x00','A')]
    return rows


def audit(sources,config_path,motion,saved=None):
    record=json.loads((motion/'manifest.json').read_text());config=json.loads(config_path.read_text())
    require(record['config']==config and record['config_sha256']==digest(config_path),'Motion config differs')
    require(record['kind']=='clathrin-articulated-patch' and record['units']=='nm','Wrong motion kind/units')
    assets=Path(record['assembly_directory']);assembly=json.loads((assets/'assembly.json').read_text());source=assembly['source_config']
    require(source['source_id']=='3IYV' and digest(assets/'assembly.json')==record['assembly_sha256'],'Unqualified assembly')
    for name,expected in record['data_sha256'].items():require(digest(motion/name)==expected,'Changed motion data')
    for name,expected in assembly['data_sha256'].items():require(digest(assets/name)==expected,'Changed source assets')
    require(digest(sources/'3IYV.cif')==source['source_sha256'],'Raw asymmetric source changed')
    require(digest(sources/'3IYV-assembly1.cif')==source['reference_sha256'],'Raw complete assembly changed')
    rows=source_rows(sources/'3IYV.cif');reference=source_rows(sources/'3IYV-assembly1.cif')
    operators={o.name:o for g in gemmi.read_structure(str(sources/'3IYV.cif')).assemblies[0].generators for o in g.operators}
    units={u['id']:u for u in assembly['instances']};moving=config['moving_unit'];visible=[moving]+config['fixed_units']
    points={};identities={};numbers=[];owners=[]
    for key in visible:
        unit=units[key];operator=operators[unit['operator']];rotation=np.array(operator.transform.mat);shift=np.array(list(operator.transform.vec))*.1
        coordinates=[];ids=[]
        for chain in unit['chains']:
            selected=rows[chain];xyz=np.array([r[2] for r in selected])@rotation.T+shift
            other=reference[unit['reference_chains'][chain]]
            require([r[:2] for r in selected]==[r[:2] for r in other] and np.max(abs(xyz-np.array([r[2] for r in other])))<1e-4,'Reference identity/coordinates differ')
            expected=list(range(1,1631)) if len(selected)>100 else list(range(95,165))
            require([r[0] for r in selected]==list(map(str,expected)),'Source is not the reviewed whole modeled chain')
            require(np.max(np.linalg.norm(np.diff(xyz,axis=0),axis=1))<1,'Disconnected source chain')
            coordinates.append(xyz);ids.extend(f'3IYV:assembly1:{unit["reference_chains"][chain]}:{n}:{r}:CA' for n,r,_ in selected)
            if key==moving:numbers.extend(expected);owners.extend([chain]*len(selected))
        points[key]=np.concatenate(coordinates);identities[key]=ids
    require(record['atom_ids']==identities,'Source identities/order changed')
    source_ca=points[moving];obstacle=np.concatenate([points[k] for k in config['fixed_units']]);numbers=np.array(numbers);owners=np.array(owners)
    data=dict(np.load(motion/'motion.npz'));frames=record['frames']
    # The shared asset writer stores CA coordinates as float32 before unit conversion.
    source_error=max(float(np.max(abs(data['source_ca']-source_ca))),float(np.max(abs(data['obstacle_ca']-obstacle))))
    require(source_error<5e-6,'Prepared source coordinates differ beyond float32 conversion')
    expected_joints=[(c,h['pivot_residue'],h['maximum_degrees']) for c in units[moving]['chains'][:3] for h in config['hinges']]
    require([(j['chain'],j['pivot_residue']) for j in record['joints']]==[r[:2] for r in expected_joints] and
            np.allclose([np.rad2deg(j['maximum_radians']) for j in record['joints']],[r[2] for r in expected_joints],atol=1e-10,rtol=0),'Unreviewed hinge selection')
    pivots=[source_ca[np.flatnonzero((owners==c)&(numbers==n))[0]] for c,n,_ in expected_joints]
    require(np.allclose([j['pivot_nm'] for j in record['joints']],pivots,atol=5e-6,rtol=0),'Hinge pivot does not follow raw CA')
    if saved is None:
        poses=data['moving_ca'];eulers=Rotation.from_rotvec(data['angles'].reshape(-1,3)).as_euler('xyz').reshape(-1,len(pivots),3)
        matrices=np.repeat(np.eye(4)[None],frames+1,axis=0);matrices[:,:3,:3]=data['rotations']
        matrices[:,:3,3]=data['centers']-np.einsum('fij,j->fi',data['rotations'],source_ca.mean(0));step=1
        export_record=None;sample_data=None
    else:
        export_record=json.loads((saved/'manifest.json').read_text());sample_data=dict(np.load(saved/'samples.npz'))
        require(export_record['motion_sha256']==digest(motion/'manifest.json') and export_record['samples_sha256']==digest(saved/'samples.npz'),'Changed saved-scene sample binding')
        poses=sample_data['moving'];matrices=sample_data['matrices'];eulers=sample_data['eulers'];step=4
        require(np.array_equal(sample_data['times'],np.arange(frames*4+1)/4+1),'Sample times differ')
        require(np.max(abs(sample_data['fixed']-obstacle))<2e-5,'Fixed source units changed during motion')
        require(export_record['maximum_fixed_pose_change_nm']<2e-5,'Fixed scene geometry changed between samples')
        require(np.max(abs(poses[::4]-data['moving_ca']))<2e-5,'Saved integer poses differ from prepared motion')
    require(poses.shape==(frames*step+1,len(source_ca),3) and np.isfinite(poses).all(),'Missing/nonfinite pose inventory')
    require(np.max(abs(poses[0]-poses[-1]))<2e-5,'Motion does not close')
    require(np.max(abs(poses[frames*step//2]-source_ca))<2e-5,'Bound endpoint lost source coordinates')
    adjacent=np.flatnonzero((owners[:-1]==owners[1:])&(numbers[1:]==numbers[:-1]+1));lengths=np.linalg.norm(source_ca[adjacent+1]-source_ca[adjacent],axis=1)
    segments=np.zeros(len(source_ca),int)
    for leg,c in enumerate(units[moving]['chains'][:3]):
        for i,h in enumerate(config['hinges']):
            low=1 if i==0 else config['hinges'][i-1]['pivot_residue']
            segments[(owners==c)&(numbers>=low)&(numbers<h['pivot_residue'])]=leg*len(config['hinges'])+i+1
    def internal(p,radius):
        pairs=cKDTree(p).query_pairs(radius,output_type='ndarray');i,j=pairs.T
        return pairs[(segments[i]!=segments[j])&~((owners[i]==owners[j])&(abs(numbers[i]-numbers[j])<=3))]
    guard=config['guard_nm'];allowance=config['source_contact_tolerance_nm'];tree=cKDTree(obstacle)
    external_floor={(i,j):max(1e-6,np.linalg.norm(source_ca[i]-obstacle[j])-allowance)
                    for i,neighbors in enumerate(tree.query_ball_point(source_ca,guard)) for j in neighbors}
    internal_floor={(i,j):max(1e-6,np.linalg.norm(source_ca[i]-source_ca[j])-allowance) for i,j in internal(source_ca,guard)}
    maximum_shape_error=0.;maximum_length_error=0.;minimum_slack=float('inf');maximum_angles=np.zeros(len(pivots));margin=1.
    bounds=[]
    if saved is not None:
        for chain in units[moving]['chains']:
            vertices=sample_data['surface_'+chain];raw=dict(np.load(chain_file(assets/'prototypes', chain)))
            require(np.allclose(vertices,raw['vertices'],atol=2e-6,rtol=0),'Exported display surface differs')
            weights=sample_data.get('weights_'+chain)
            if weights is not None:
                nearest=raw['residue_numbers'][cKDTree(raw['atoms']).query(vertices)[1]]
                require(np.array_equal(weights,nearest),'Surface residue weights differ from source CA')
            for label in np.unique(np.digitize(weights,[h['pivot_residue'] for h in config['hinges']])) if weights is not None else [0]:
                selection=np.ones(len(vertices),bool) if weights is None else np.digitize(weights,[h['pivot_residue'] for h in config['hinges']])==label
                v=vertices[selection];corners=np.array(list(product(*zip(v.min(0),v.max(0)))))
                bounds.append((chain,int(weights[selection][0]) if weights is not None else 2000,corners))
        camera=np.array(export_record['camera_matrix']);aspect=export_record['resolution'][0]/export_record['resolution'][1];scale=export_record['ortho_scale']
        for key in config['fixed_units']:
            unit=units[key];matrix=np.array(unit['matrix_nm'])
            for chain in unit['chains']:
                v=sample_data['surface_'+chain];corners=np.array(list(product(*zip(v.min(0),v.max(0)))))
                view=(corners@matrix[:3,:3].T+matrix[:3,3]-camera[:3,3])@camera[:3,:3]
                coordinates=view[:,:2]/[scale,scale/aspect]+.5
                margin=min(margin,float(np.min(np.c_[coordinates,1-coordinates])))
                require(np.all(view[:,2]<0) and margin>0,'Whole fixed surface leaves camera')
    for pose,matrix,euler in zip(poses,matrices,eulers):
        body_rotation=matrix[:3,:3]
        require(np.allclose(body_rotation.T@body_rotation,np.eye(3),atol=2e-6,rtol=0) and abs(np.linalg.det(body_rotation)-1)<2e-6,'Global scale/reflection detected')
        expected=source_ca.copy();rotations=Rotation.from_euler('xyz',euler);joint_angles=rotations.magnitude()
        maximum_angles=np.maximum(maximum_angles,joint_angles)
        require(np.all(joint_angles<=np.deg2rad([r[2] for r in expected_joints])+2e-6),'Joint angle exceeds reviewed bound')
        for (chain,number,_),pivot,rotation in zip(expected_joints,pivots,rotations.as_matrix()):
            selected=(owners==chain)&(numbers<number);expected[selected]=(expected[selected]-pivot)@rotation.T+pivot
        expected=expected@body_rotation.T+matrix[:3,3];shape_error=float(np.max(abs(pose-expected)));maximum_shape_error=max(maximum_shape_error,shape_error)
        require(shape_error<2e-5,'Actual deformation differs from rigid-segment source motion')
        error=float(np.max(abs(np.linalg.norm(pose[adjacent+1]-pose[adjacent],axis=1)-lengths)));maximum_length_error=max(maximum_length_error,error)
        require(error<2e-5,'Consecutive source CA length changed')
        ext=np.array([(i,j) for i,neighbors in enumerate(tree.query_ball_point(pose,guard)) for j in neighbors],int).reshape(-1,2);inter=internal(pose,guard)
        slack=np.r_[[np.linalg.norm(pose[i]-obstacle[j])-external_floor.get((i,j),guard) for i,j in ext],
                     [np.linalg.norm(pose[i]-pose[j])-internal_floor.get((i,j),guard) for i,j in inter]]
        if len(slack):minimum_slack=min(minimum_slack,float(slack.min()))
        require(np.all(slack>=-2e-5),'New/worsened contact outside exact source-pair allowance')
        if saved is not None:
            for chain,number,corners in bounds:
                corners=corners.copy()
                for (owner,pivot_number,_),pivot,rotation in zip(expected_joints,pivots,rotations.as_matrix()):
                    if chain==owner and number<pivot_number:corners=(corners-pivot)@rotation.T+pivot
                view=(corners@body_rotation.T+matrix[:3,3]-camera[:3,3])@camera[:3,:3]
                coordinates=view[:,:2]/[scale,scale/aspect]+.5
                margin=min(margin,float(np.min(np.c_[coordinates,1-coordinates])))
                require(np.all(view[:,2]<0) and margin>0,'Whole moving surface leaves camera')
    origin=source_ca.mean(0);body_travel=float(np.linalg.norm(matrices[0,:3,:3]@origin+matrices[0,:3,3]-origin))
    require(body_travel>=config['travel_nm']-2e-5,'Insufficient body-reference release travel')
    require(tree.query(poses[0])[0].min()>=config['park_clearance_nm'],'Parked molecule remains in contact')
    return {'passed':True,'mode':record['mode'],'poses':len(poses),'source_CA_points':len(source_ca),'fixed_CA_points':len(obstacle),
        'raw_source_float32_error_nm':source_error,'maximum_actual_kinematic_error_nm':maximum_shape_error,'maximum_CA_step_change_nm':maximum_length_error,
        'maximum_hinge_degrees':np.rad2deg(maximum_angles).tolist(),'minimum_contact_slack_nm':minimum_slack,
        'body_reference_travel_nm':body_travel,'CA_centroid_travel_nm':float(np.linalg.norm(poses[0].mean(0)-origin)),
        'source_external_pairs':len(external_floor),'source_internal_pairs':len(internal_floor),
        'saved_scene':saved is not None,'moving_surface_camera_margin':margin if saved is not None else None,
        'source_hashes':{p.name:digest(p) for p in [sources/'3IYV.cif',sources/'3IYV-assembly1.cif']},
        'limits':'CA geometry with exact source-pair allowances, sampled at stated times; no all-atom energy, continuous-time or physiological pathway claim.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('sources','config','motion','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--saved-samples',type=Path);a=p.parse_args();result=audit(a.sources,a.config,a.motion,a.saved_samples)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
