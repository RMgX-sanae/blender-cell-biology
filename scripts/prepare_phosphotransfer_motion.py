"""Sample two PKA turnovers with stable atom IDs, source-guided opening and QC."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix
from scipy.spatial.transform import Rotation

from interaction_motion import matrices, progress, source_contacts, contacts


STATES=['closed','product','apo','bound_open','product_bound_open']
INTERVALS={'atp_in':(.04,.20),'m1':(.23,.32),'peptide_in':(.35,.48),
           'peptide_out':(.69,.80),'m1_out':(.82,.87),'adp_out':(.89,.98)}


def smooth(x):
    u=np.clip(x,0,1);return u*u*u*(10+u*(-15+6*u))


def pose_at(t,name,assets,manifest):
    start,end=INTERVALS[name];q=np.clip((t-start)/(end-start),0,1)
    departure=name.endswith('_out');fraction=.5+q/2 if departure else q/2
    spec=manifest['paths'][name];mobile=np.load(assets/(name+'-path.npz'))['moving']
    pose=matrices(mobile,spec['direction'],spec['distance_nm'],np.array([fraction]),
                  'loop',spec['turn_radians'],spec['bend_fraction'])[0]
    binding,_=progress(np.array([fraction]),'loop');free=1-float(binding[0])
    # Extend only the distant portion beyond both view windows. The screened
    # pocket approach is unchanged; the final full timeline is checked again.
    pose[:3,3]+=np.array(spec['direction'])*13*smooth((free-.6)/.4)
    opening=min(1,free/(1-spec['obstacle_opening']['closure_starts_at_progress']))
    return pose,opening


def apply(points,pose):return points@pose[:3,:3].T+pose[:3,3]


def rotation_between(first,second):
    x=first/np.linalg.norm(first);y=second/np.linalg.norm(second)
    axis=np.cross(x,y);length=np.linalg.norm(axis)
    if length<1e-9:return Rotation.identity()
    return Rotation.from_rotvec(axis/length*np.arctan2(length,np.dot(x,y)))


def group_interpolate(first,second,q):
    x=first-first[0];y=second-second[0]
    u,_,vt=np.linalg.svd(x.T@y);r=vt.T@np.diag([1,1,np.linalg.det(vt.T@u.T)])@u.T
    motion=Rotation.from_rotvec(Rotation.from_matrix(r).as_rotvec()*q)
    residual=y-x@r.T
    return first[0]*(1-q)+second[0]*q+motion.apply(x)+residual*q


def covalent_interpolate(first,second,bonds,q,initial=None,max_nfev=160,ftol=1e-8):
    """Geometric interpolation constrained by bonds and adjacent bond angles.

    This repairs collapsed Cartesian intermediates, not a molecular dynamics
    trajectory. The source endpoints and atom correspondence remain unchanged.
    """
    if q<=0:return first.copy()
    if q>=1:return second.copy()
    neighbors=[set() for _ in first]
    pairs={tuple(sorted(edge)) for edge in bonds}
    for i,j in bonds:neighbors[i].add(j);neighbors[j].add(i)
    for connected in neighbors:
        for i in connected:
            for j in connected:
                if i<j:pairs.add((i,j))
    edges=np.array(sorted(pairs));target=first*(1-q)+second*q
    lengths=(1-q)*np.linalg.norm(first[edges[:,0]]-first[edges[:,1]],axis=1)+q*np.linalg.norm(second[edges[:,0]]-second[edges[:,1]],axis=1)
    count=first.size;sparsity=lil_matrix((count+len(edges),count),dtype=int)
    sparsity.setdiag(1)
    for k,(i,j) in enumerate(edges):
        sparsity[count+k,3*i:3*i+3]=1;sparsity[count+k,3*j:3*j+3]=1
    def residual(flat):
        xyz=flat.reshape(-1,3)
        distance=np.linalg.norm(xyz[edges[:,0]]-xyz[edges[:,1]],axis=1)
        return np.r_[(xyz-target).ravel()*.15,(distance-lengths)*12]
    seed=target if initial is None else np.asarray(initial)
    fit=least_squares(residual,seed.ravel(),jac_sparsity=sparsity.tocsr(),max_nfev=max_nfev,ftol=ftol)
    if not fit.success:raise ValueError('Peptide geometric interpolation did not converge')
    return fit.x.reshape(-1,3)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cycles',type=int,default=2);p.add_argument('--frames-per-cycle',type=int,default=288)
    p.add_argument('--fps',type=int,default=24)
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh motion directory')
    if a.cycles<1 or a.frames_per_cycle<96 or a.fps<1:p.error('Use positive cycles/fps and at least 96 frames per cycle')
    m=json.loads((a.asset_dir/'manifest.json').read_text())
    if not m['path_screen_passed']:p.error('Fix the failed source-guided paths before animation')
    d={k:dict(np.load(a.asset_dir/(k+'.npz'))) for k in ['kinase','peptide','nucleotide_core','phosphoryl','metals']}
    kinase=d['kinase'];states=np.array([kinase[k] for k in ['atoms','product_atoms','open_atoms','bound_open_atoms','product_bound_open_atoms']])
    peptide=d['peptide'];core=d['nucleotide_core'];phosphoryl=d['phosphoryl'];metals=d['metals']
    atom_index={(int(res),str(name)):i for i,(res,name) in enumerate(zip(peptide['residue_numbers'],peptide['atom_names']))}
    cb=atom_index[621,'CB'];og=atom_index[621,'OG']
    rotation=rotation_between(peptide['atoms'][og]-peptide['atoms'][cb],peptide['product_atoms'][og]-peptide['product_atoms'][cb])
    center=peptide['product_atoms'][cb]
    # Non-bridging oxygen labels do not establish measured atom trajectories
    # between crystals. Choose the equivalent labeling with the smallest rigid
    # rotation before the substrate hydroxyl turns; record this inferred map.
    candidates=[];first=phosphoryl['atoms']-phosphoryl['atoms'][0]
    for permutation in itertools.permutations([1,2,3]):
        indices=[0,*permutation]
        endpoint=center+rotation.inv().apply(phosphoryl['product_atoms'][indices]-center)
        target=endpoint-endpoint[0]
        u,_,vt=np.linalg.svd(first.T@target)
        matrix=vt.T@np.diag([1,1,np.linalg.det(vt.T@u.T)])@u.T
        candidates.append((float(np.linalg.norm(Rotation.from_matrix(matrix).as_rotvec())),indices))
    _,oxygen_order=min(candidates)
    phosphoryl['product_atoms']=phosphoryl['product_atoms'][oxygen_order]
    prephosphate=center+rotation.inv().apply(phosphoryl['product_atoms']-center)
    atom_sets={k:[] for k in ['peptide','nucleotide_core','phosphoryl','M1','M2']};weights=[];transfers=[];turns=[]
    for frame in range(a.frames_per_cycle+1):
        t=frame/a.frames_per_cycle;poses={};opening={}
        for name in INTERVALS:poses[name],opening[name]=pose_at(t,name,a.asset_dir,m)
        chemical=smooth((t-.50)/.04);transfer=smooth((t-.54)/.06);turn=smooth((t-.60)/.06)
        # Common peptide coordinates interpolate between snapshots. The Ser21
        # hydroxyl follows a rotation, rather than collapsing through its C-beta.
        pep=covalent_interpolate(peptide['atoms'],peptide['product_atoms'],peptide['bonds'],chemical)
        v0=peptide['atoms'][og]-peptide['atoms'][cb];v1=peptide['product_atoms'][og]-peptide['product_atoms'][cb]
        rotated=Rotation.from_rotvec(rotation.as_rotvec()*turn).apply(v0/np.linalg.norm(v0))
        pep[og]=pep[cb]+rotated*((1-chemical)*np.linalg.norm(v0)+chemical*np.linalg.norm(v1))
        nucleotide=core['atoms']*(1-transfer)+core['product_atoms']*transfer
        group=group_interpolate(phosphoryl['atoms'],prephosphate,transfer)
        if turn>0:group=center+Rotation.from_rotvec(rotation.as_rotvec()*turn).apply(prephosphate-center)
        ions=metals['atoms']*(1-transfer)+metals['product_atoms']*transfer
        pep_pose=poses['peptide_in'] if t<.69 else poses['peptide_out']
        nt_pose=poses['atp_in'] if t<.89 else poses['adp_out']
        atom_sets['peptide'].append(apply(pep,pep_pose))
        atom_sets['nucleotide_core'].append(apply(nucleotide,nt_pose))
        group_pose=poses['atp_in'] if t<.54 else poses['peptide_out'] if t>=.69 else np.eye(4)
        atom_sets['phosphoryl'].append(apply(group,group_pose))
        atom_sets['M1'].append(apply(ions[:1],poses['m1'] if t<.82 else poses['m1_out']))
        atom_sets['M2'].append(apply(ions[1:],nt_pose))
        w=np.zeros(5)
        if t<.20:
            w[0]=1-opening['atp_in'];w[2]=opening['atp_in']
        elif t<.33:w[0]=1
        elif t<.35:
            q=smooth((t-.33)/.02);w[0]=1-q;w[3]=q
        elif t<.48:w[0]=1-opening['peptide_in'];w[3]=opening['peptide_in']
        elif t<.69:w[0]=1-transfer;w[1]=transfer
        elif t<.80:w[1]=1-opening['peptide_out'];w[4]=opening['peptide_out']
        elif t<.89:w[4]=1
        else:w[4]=1-opening['adp_out'];w[2]=opening['adp_out']
        weights.append(w);transfers.append(transfer);turns.append(turn)
    arrays={k:np.array(v) for k,v in atom_sets.items()};weights=np.array(weights)
    # Check moving components against the actual frame-dependent enzyme.
    # Exceptions are only atom pairs already present at the appropriate source
    # endpoint; accommodating an apo rotamer does not add new exceptions.
    protein=np.einsum('fs,snc->fnc',weights,states);poses=np.repeat(np.eye(4)[None],len(weights),axis=0)
    checks={}
    for key,points in arrays.items():
        initial=d[key]['atoms'] if key in d else metals['atoms'][:1] if key=='M1' else metals['atoms'][1:]
        final=d[key]['product_atoms'] if key in d else metals['product_atoms'][:1] if key=='M1' else metals['product_atoms'][1:]
        allowed={**source_contacts(initial,states[0],.2),**source_contacts(final,states[1],.2)}
        bad=[];minimum=1
        for f,pts in enumerate(points):
            count,dist=contacts(pts,protein[f],poses[f:f+1],.2,allowed)
            if count[0]:bad.append({'frame':f+1,'atoms':int(count[0])})
            minimum=min(minimum,float(dist[0]))
        checks[key]={'flagged_frames':bad,'minimum_distance_nm':minimum}
    # Covalent geometry of the rendered peptide and nucleotide cores must stay
    # finite and reasonably close to their two source endpoint bond lengths.
    bond_checks={}
    for key in ['peptide','nucleotide_core','phosphoryl']:
        edges=d[key]['bonds'];pts=arrays[key]
        lengths=np.linalg.norm(pts[:,edges[:,0]]-pts[:,edges[:,1]],axis=2)
        bond_checks[key]={'min_nm':float(lengths.min()),'max_nm':float(lengths.max()),
            'within_display_bounds':bool(lengths.min()>.08 and lengths.max()<.22)}
    # Reaction partner contacts include the source donor and acceptor covalent
    # bonds. Only those exact atom pairs, at their endpoint distances, are exempt.
    partner_checks={};keys=list(arrays)
    def endpoint(key,product):
        field='product_atoms' if product else 'atoms'
        if key in d:return d[key][field]
        return metals[field][:1] if key=='M1' else metals[field][1:]
    for i,key in enumerate(keys):
        for other in keys[i+1:]:
            allowed={**source_contacts(endpoint(key,False),endpoint(other,False),.2),
                     **source_contacts(endpoint(key,True),endpoint(other,True),.2)}
            bad=[];minimum=1
            for f in range(len(weights)):
                count,dist=contacts(arrays[key][f],arrays[other][f],poses[f:f+1],.2,allowed)
                if count[0]:bad.append({'frame':f+1,'atoms':int(count[0])})
                minimum=min(minimum,float(dist[0]))
            partner_checks[key+' / '+other]={'flagged_frames':bad,'minimum_distance_nm':minimum}
    identity=m['chemical_identity'];passed=not any(v['flagged_frames'] for v in [*checks.values(),*partner_checks.values()]) and all(v['within_display_bounds'] for v in bond_checks.values())
    a.output.mkdir(parents=True)
    np.savez_compressed(a.output/'cycle.npz',**arrays,enzyme_weights=weights,transfer=np.array(transfers),product_turn=np.array(turns))
    result={'asset_dir':str(a.asset_dir.resolve()),'frames':a.cycles*a.frames_per_cycle,'fps':a.fps,'cycles':a.cycles,
        'frames_per_cycle':a.frames_per_cycle,'mode':'teaching','intervals':INTERVALS,'enzyme_states':STATES,
        'identity':identity,'enzyme_contacts':checks,'partner_contacts':partner_checks,'bond_geometry':bond_checks,'passed':passed,
        'phosphoryl_atom_map':dict(zip(phosphoryl['atom_names'].tolist(),np.array(['P','O1P','O2P','O3P'])[oxygen_order].tolist())),
        'scope':'Two independent substrate/ATP packets; products leave the view. No on-screen product is reset to substrate.',
        'limits':['Atomic-center and display-bond checks are not a force field or dynamics validation.',
            'Peptide intermediates use geometric bond and adjacent-angle constraints, not an energy model. Final saved Blender views still require inspection.',
            'Equivalent non-bridging phosphate oxygen labels are matched for a smaller illustrative rotation; their atom trajectories are not measured.',
            'This is a catalytic-subunit example. It does not include cAMP production or PKA holoenzyme activation.']}
    (a.output/'manifest.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
    print(json.dumps({'passed':passed,'flagged_frames':{k:len(v['flagged_frames']) for k,v in checks.items()},
        'partner_flagged_frames':{k:len(v['flagged_frames']) for k,v in partner_checks.items()},'bonds':bond_checks}))


if __name__=='__main__':main()
