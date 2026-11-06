"""Qualify a coarse articulated clathrin admission/release from reviewed source units."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from articulated_ca_motion import Contacts, plan, sample, shape
from interaction_motion import progress


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def model(assets, config):
    manifest=json.loads((assets/'assembly.json').read_text())
    if manifest['source_config']['source_id']!='3IYV' or not manifest['grouping_check']['trace_check']['passed']:
        raise ValueError('Use the reviewed full-triskelion 3IYV input')
    for file,expected in manifest['data_sha256'].items():
        if digest(assets/file)!=expected:raise ValueError('Prepared assembly data changed: '+file)
    units={r['id']:r for r in manifest['instances']};moving=units[config['moving_unit']]
    if not np.allclose(moving['matrix_nm'],np.eye(4)):raise ValueError('This scene adapter requires an identity-source moving unit')
    if len(set(config['fixed_units']+[config['moving_unit']]))!=len(config['fixed_units'])+1:
        raise ValueError('Moving and fixed source units must be distinct')
    identities={};points={};residues={};chains={}
    for key in [config['moving_unit']]+config['fixed_units']:
        unit=units[key];matrix=np.array(unit['matrix_nm']);coordinates=[];ids=[];numbers=[];owners=[]
        for index,c in enumerate(unit['chains']):
            data=np.load(assets/'prototypes'/('chain_'+c+'.npz'))
            coordinates.append(data['atoms']@matrix[:3,:3].T+matrix[:3,3])
            ids.extend(f'3IYV:assembly1:{unit["reference_chains"][c]}:{int(n)}:{r}:CA'
                       for n,r in zip(data['residue_numbers'],data['residue_names']))
            numbers.extend(data['residue_numbers']);owners.extend([index]*len(data['atoms']))
        points[key]=np.concatenate(coordinates);identities[key]=ids
        residues[key]=np.array(numbers);chains[key]=np.array(owners)
    key=config['moving_unit'];source=points[key];origin=source.mean(0);joints=[];segments=np.zeros(len(source),int)
    group=next(g for g in manifest['source_config']['groups'] if g['id']==moving['group'])
    for leg,c in enumerate(group['heavy']):
        previous=1
        for ordinal,spec in enumerate(config['hinges']):
            pivot=spec['pivot_residue'];mask=(chains[key]==leg)&(residues[key]<pivot)
            index=np.flatnonzero((chains[key]==leg)&(residues[key]==pivot))
            if len(index)!=1 or pivot<=previous:raise ValueError('Require ordered present hinge pivots')
            joints.append({'indices':np.flatnonzero(mask),'pivot':source[index[0]]-origin,
                           'maximum_radians':np.deg2rad(spec['maximum_degrees']),
                           'chain':c,'pivot_residue':pivot,'pivot_nm':source[index[0]].tolist()})
            segments[(chains[key]==leg)&(residues[key]>=previous)&(residues[key]<pivot)]=leg*len(config['hinges'])+ordinal+1
            previous=pivot
    obstacle=np.concatenate([points[k] for k in config['fixed_units']])
    contacts=Contacts(source,obstacle,segments,chains[key],residues[key],config['guard_nm'],config['source_contact_tolerance_nm'])
    return manifest,points,identities,joints,segments,chains[key],residues[key],contacts


def prepare(assets, config_path, output, frames, fps, mode):
    if output.exists():raise ValueError('Use a fresh output directory')
    if frames<24 or frames%4 or fps<=0:raise ValueError('Use positive fps and a frame count divisible by four')
    config=json.loads(config_path.read_text())
    manifest,points,identities,joints,segments,chains,residues,contacts=model(assets,config)
    source=points[config['moving_unit']];origin=source.mean(0)
    unit=next(u for u in manifest['instances'] if u['id']==config['moving_unit'])
    path=plan(source,joints,contacts,unit['hub_nm'],config['travel_nm'],config['park_clearance_nm'])
    fractions=1-progress(np.arange(frames*4+1)/(frames*4),mode)[0]
    centers,rotations,angles=sample(path,fractions);maximum_bond_error=0.;poses=[]
    adjacent=np.flatnonzero((chains[:-1]==chains[1:])&(residues[1:]==residues[:-1]+1))
    original_lengths=np.linalg.norm(source[adjacent+1]-source[adjacent],axis=1)
    for index,(center,rotation,angle) in enumerate(zip(centers,rotations,angles)):
        pose=shape(source-origin,joints,angle)@rotation.T+center
        if not contacts.safe(pose):raise ValueError('Final articulated schedule violates a source-pair contact bound')
        error=float(np.max(abs(np.linalg.norm(pose[adjacent+1]-pose[adjacent],axis=1)-original_lengths)))
        maximum_bond_error=max(maximum_bond_error,error)
        if error>1e-8:raise ValueError('Consecutive source CA distances changed')
        if index%4==0:poses.append(pose)
    poses=np.array(poses)
    if np.max(abs(poses[0]-poses[-1]))>1e-8 or np.max(abs(poses[frames//2]-source))>1e-8:
        raise ValueError('Loop closure or source-bound endpoint differs')
    output.mkdir(parents=True)
    np.savez_compressed(output/'path.npz',**path)
    np.savez_compressed(output/'motion.npz',centers=centers[::4],rotations=rotations[::4],angles=angles[::4],
                        moving_ca=poses,source_ca=source,obstacle_ca=contacts.obstacle,segments=segments,chains=chains,residues=residues)
    record={'kind':'clathrin-articulated-patch','passed':True,'units':'nm','frames':frames,'fps':fps,'mode':mode,
            'config':config,'config_sha256':digest(config_path),'assembly_directory':str(assets.resolve()),
            'assembly_sha256':digest(assets/'assembly.json'),'atom_ids':identities,'origin_nm':origin.tolist(),
            'joints':[{k:v for k,v in j.items() if k not in {'indices','pivot'}} for j in joints],
            'checks':{'passed':True,'sampled_poses':len(fractions),'maximum_consecutive_CA_distance_error_nm':maximum_bond_error,
                      'source_external_pairs':len(contacts.external_floor),'source_internal_pairs':len(contacts.internal_floor),
                      'policy':'All moving-fixed and changing internal CA pairs. Same rigid segments are invariant; same-chain neighbors within3 residues are excluded. Identical reviewed source pairs only may deepen by the declared tolerance.'},
            'path_waypoints':len(path['centers']),'scope':config['scope'],'limits':config['limits'],
            'data_sha256':{name:digest(output/name) for name in ('path.npz','motion.npz')}}
    (output/'manifest.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'passed':True,'poses':len(poses),'sampled_poses':len(fractions),'waypoints':len(path['centers'])}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('assets','config','output'):parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--frames',type=int,default=240);parser.add_argument('--fps',type=int,default=24)
    parser.add_argument('--mode',choices=['loop','teaching'],default='teaching');args=parser.parse_args()
    prepare(args.assets,args.config,args.output,args.frames,args.fps,args.mode)
