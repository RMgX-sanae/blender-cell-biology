"""Check every source member, shared transform, restraint and candidate pose.

Assembly selections are also compared with independent deposited coordinates.
This validates geometric handoff, not affinity or biological binding modes.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

from interaction_assets import verify_assets
from prepare_binding_candidates import load_partner, residue_id, rigid, sha


def require(value,message):
    if not value:raise ValueError(message)


def check(sources,prepared):
    report_path=prepared/'candidates.json';report=json.loads(report_path.read_text())
    folder=prepared/'interaction';m=json.loads((folder/'manifest.json').read_text())
    verify_assets(folder,m);config=report['config']
    require(m['source_config']==config and m['binding_candidate']['report_sha256']==sha(report_path),
            'Candidate report/config binding differs')
    candidates=[c for c in report['candidates'] if c['id']==m['binding_candidate']['id']]
    require(len(candidates)==1,'Candidate identity is ambiguous')
    selected=rigid(candidates[0]['matrix_nm'])
    coordinates={};polymers={};polymer_points={};errors={};proofs={};counts={}
    for role in ('anchor','moving'):
        members,proof=load_partner(sources,config[role],role)
        if proof:proofs[role]=proof
        declared=[a for a in m['actors'] if a['role']==role]
        require({a['id'] for a in declared}=={a['id'] for a in members},'Partner member inventory differs')
        base=rigid(report['monomer_to_base_matrices_angstrom'][role]);all_points=[];poly=[];records=[]
        for member in members:
            actor=next(a for a in declared if a['id']==member['id']);rows=member['rows']
            require(all(actor.get(k)==v for k,v in member.items() if k!='rows'),'Source member declaration differs')
            data=np.load(folder/actor['file'],allow_pickle=False)
            require(len(rows)==actor['atoms']==len(data['atoms']),'Source member atom count differs')
            for key,index in [('elements',1),('atom_names',2)]:
                require(np.array_equal(data[key],[r[index] for r in rows]),'Source atom identity differs: '+key)
            if member['kind']=='polymer':
                for key,index in [('residue_numbers',3),('residue_names',4)]:
                    require(np.array_equal(data[key],[r[index] for r in rows]),'Source residue identity differs: '+key)
                codes=data['insertion_codes'] if 'insertion_codes' in data else np.full(len(rows),'')
                require(np.array_equal(codes,[r[5] for r in rows]),'Source insertion codes differ')
            raw=np.array([r[0] for r in rows])
            expected=(raw@base[:3,:3].T+base[:3,3])*.1-report['shared_origin_nm']
            if role=='moving':expected=expected@selected[:3,:3].T+selected[:3,3]
            error=float(np.abs(data['atoms']-expected).max())
            require(np.isfinite(error) and error<=1e-10,'Source member no longer follows the common partner transform')
            errors[member['id']]=error;all_points.append(expected)
            if member['kind']=='polymer':
                poly.append(expected);records.extend((member.get('instance_id'),r) for r in rows)
        coordinates[role]=np.concatenate(all_points);polymer_points[role]=np.concatenate(poly);polymers[role]=records
        counts[role]={'polymer_atoms':len(records),'companion_atoms':len(coordinates[role])-len(records),'members':len(members)}
    require(m.get('partner_assembly_checks',{})==proofs and report.get('partner_assembly_checks',{})==proofs,
            'Independent partner assembly evidence differs')
    poses=np.load(folder/'motion.npz',allow_pickle=False)['matrices'];minimum=float('inf')
    require(len(poses)==m['frames']+1,'Candidate motion pose inventory differs')
    for pose in poses:
        rigid(pose);moving=coordinates['moving']@pose[:3,:3].T+pose[:3,3]
        # Dense all-pair distances in blocks bound memory, without KD trees or
        # exempting contacts inherited from a template.
        for start in range(0,len(moving),256):
            minimum=min(minimum,float(cdist(coordinates['anchor'],moving[start:start+256]).min()))
    require(minimum>=config['search']['atom_guard_nm'],'Candidate path has deep overlaps')
    distances=[]
    for item in config['restraints']:
        selected_points=[]
        for role in ('anchor','moving'):
            mask=[instance==item.get(role+'_instance') and residue_id(row)==item[role]
                  and (not item.get(role+'_atom') or row[2]==item[role+'_atom']) for instance,row in polymers[role]]
            require(any(mask),'Restraint has no matching source-instance atom')
            selected_points.append(polymer_points[role][mask])
        distance=float(cdist(*selected_points).min());lo,hi=item['range_nm']
        require(lo-1e-4<=distance<=hi+1e-4,'Candidate misses a reviewed instance restraint')
        distances.append(distance)
    return {'passed':True,'source_members':counts,'source_coordinate_error_max_nm':errors,
            'independent_dense_path_minimum_nm':minimum,'checked_poses':len(poses),
            'restraint_distances_nm':distances,'assembly_checks':proofs,
            'scope':'Source identities and complete-partner rigid placement/path; not binding affinity, kinetics or a complete mechanism.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources',type=Path,required=True);p.add_argument('--prepared',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=check(a.sources,a.prepared);a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':True,'checked_poses':result['checked_poses'],'source_members':result['source_members']}))
