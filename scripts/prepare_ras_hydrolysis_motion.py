"""Prepare source-guided GAP association, GTP hydrolysis and GDP/Pi release states.

Heavy-atom bookkeeping and geometric screening; not a molecular dynamics model.
The selected Ras G-domain and GAP catalytic constructs omit membrane recruitment.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import gemmi
import numpy as np
from scipy.interpolate import BSpline, PchipInterpolator
from scipy.optimize import minimize, minimize_scalar
from scipy.spatial.transform import Rotation

from interaction_motion import matrices, source_contacts
from prepare_molecular_asset import atom_table, xyz
from prepare_ras_switch import fit_points
from prepare_ras_exchange_motion import pair_check, apply
from prepare_phosphotransfer_motion import covalent_interpolate, group_interpolate, rotation_between, smooth
from rigid_path_detour import refine


def merge_contacts(*records):
    result={}
    for record in records:
        for pair,distance in record.items():result[pair]=min(distance,result.get(pair,np.inf))
    return result


def gap_samples(data,route,q):
    p=data['atoms'];lookup={(int(r),str(n)):i for i,(r,n) in enumerate(zip(data['residue_numbers'],data['atom_names']))}
    specs=[(949,['CG','CD','CE','NZ']),(903,['CG','CD','NE','CZ','NH1','NH2'])]
    mobile={lookup[r,n] for r,names in specs for n in names}
    core=p[[i for i in range(len(p)) if i not in mobile]]
    poses=matrices(core,route['direction'],route['distance_nm'],q*.5,'loop',route['turn_radians'],route['bend_fraction'])
    local=np.repeat(p[None],len(q),axis=0)
    for (residue,names),record in zip(specs,route['torsions']):
        angles=np.unwrap(np.deg2rad(record['degrees']),axis=0)
        values=PchipInterpolator(np.linspace(0,1,len(angles)),angles)(q)
        selected=[lookup[residue,n] for n in names];distal=[i for i,n in enumerate(names) if n!='CG']
        ca,cb,cg=(p[lookup[residue,n]] for n in ['CA','CB','CG'])
        a1=cb-ca;a1/=np.linalg.norm(a1);a2=cg-cb;a2/=np.linalg.norm(a2)
        for f,(first,second) in enumerate(values):
            points=p[selected].copy()
            points[distal]=(points[distal]-cg)@Rotation.from_rotvec(a2*second).as_matrix().T+cg
            local[f,selected]=(points-cb)@Rotation.from_rotvec(a1*first).as_matrix().T+cb
    world=np.einsum('fij,fnj->fni',poses[:,:3,:3],local)+poses[:,None,:3,3]
    return local,world,poses,specs


def chemistry(nt,gd,source_path):
    names=list(nt['atom_names'])
    r,t=fit_points([nt['atoms'][names.index(n)] for n in gd['atom_names']],gd['atoms'])
    gtp=nt['atoms']@r.T+t
    origin=xyz(atom_table(source_path)['R']).mean(0)
    model=gemmi.read_structure(str(source_path))[0]
    water=next(res for res in model['R'] if res.name=='HOH' and res.seqid.num==230)
    w=(np.array(list(next(iter(water)).pos))-origin)*.1
    ids=[names.index(n) for n in ['PG','O1G','O2G','O3G']]
    p=gtp[ids[0]];axis=w-p;axis/=np.linalg.norm(axis);product_p=w-axis*.155
    tangent=gtp[ids[1]]-p;tangent-=axis*(tangent@axis);tangent/=np.linalg.norm(tangent)
    cross=np.cross(axis,tangent)
    def endpoint(phi,order):
        angles=phi+np.array(order)*2*np.pi/3
        directions=-axis/3+np.sqrt(8/9)*(np.cos(angles)[:,None]*tangent+np.sin(angles)[:,None]*cross)
        return np.vstack([product_p,product_p+.153*directions,w])
    best=None
    for order in itertools.permutations(range(3)):
        fit=minimize_scalar(lambda phi:float(np.sum((endpoint(phi,order)[1:4]-gtp[ids[1:]])**2)),bounds=(-np.pi,np.pi),method='bounded')
        if best is None or fit.fun<best[0]:best=(fit.fun,fit.x,order)
    return gtp,w,ids,endpoint(best[1],best[2])


def chemical_samples(gtp,water,ids,product,fractions):
    initial=np.vstack([gtp[ids],water]);output=[]
    for q in fractions:
        coordinates=initial.copy();coordinates[0]=initial[0]*(1-q)+product[0]*q
        for j in range(1,4):
            first=initial[j]-initial[0];last=product[j]-product[0]
            turn=rotation_between(first,last)
            direction=Rotation.from_rotvec(turn.as_rotvec()*q).apply(first/np.linalg.norm(first))
            coordinates[j]=coordinates[0]+direction*((1-q)*np.linalg.norm(first)+q*np.linalg.norm(last))
        output.append(coordinates)
    return np.array(output)


def refine_magnesium(actors,allowed):
    """Smooth local ion accommodation while retaining observed coordination pairs."""
    groups=['ras','gap','core','gamma','water'];points=actors['mg']
    counts=np.zeros(len(points),dtype=int);limits=[]
    for group in groups:
        allow=allowed.get('mg_'+group,{})
        bad,_=pair_check(points,actors[group],allow);counts+=bad
        limit=np.full(actors[group].shape[1],.20)
        for (i,j),distance in allow.items():
            assert i==0
            limit[j]=distance-.01
        limits.extend(limit)
    bad=np.flatnonzero(counts)
    if not len(bad):return points,None
    margin=max(10,round(len(points)*.035));first=max(1,int(bad[0])-margin);last=min(len(points)-2,int(bad[-1])+margin)
    frames=np.arange(first,last+1);u=(frames-first)/(last-first)
    knots=np.r_[np.zeros(4),np.linspace(0,1,7)[1:-1],np.ones(4)];n=len(knots)-4
    basis=BSpline(knots,np.eye(n),3)(u)[:,2:-2]
    fixed=np.concatenate([actors[group][frames] for group in groups],axis=1)
    target=np.asarray(limits)+.003
    def objective(flat):
        controls=flat.reshape(-1,3);offset=basis@controls
        delta=points[frames,0,None,:]+offset[:,None,:]-fixed
        distance=np.linalg.norm(delta,axis=2);penetration=np.maximum(target[None]-distance,0)
        loss=float(np.sum(controls**2))*.05+float(np.sum(penetration**2))*10000
        gradient=controls*.1
        differences=np.diff(controls,n=2,axis=0);loss+=float(np.sum(differences**2))*.2
        gradient[:-2]+=.4*differences;gradient[1:-1]-=.8*differences;gradient[2:]+=.4*differences
        force=-20000*np.sum(penetration[:,:,None]*delta/np.maximum(distance[:,:,None],1e-8),axis=1)
        gradient+=basis.T@force
        return loss,gradient.ravel()
    fit=minimize(objective,np.zeros((n-4)*3),jac=True,method='L-BFGS-B',bounds=[(-.1,.1)]*((n-4)*3),
                 options={'maxiter':500,'ftol':1e-12,'gtol':1e-7,'maxls':40})
    controls=fit.x.reshape(-1,3);offset=basis@controls;result=points.copy();result[frames,0]+=offset
    distances=np.linalg.norm(result[frames,0,None,:]-fixed,axis=2)
    record={'first_index':first,'last_index':last,'knots':knots.tolist(),
        'controls_nm':np.vstack([np.zeros((2,3)),controls,np.zeros((2,3))]).tolist(),
        'maximum_translation_nm':float(np.linalg.norm(offset,axis=1).max()),
        'minimum_clearance_above_pair_limit_nm':float(np.min(distances-np.asarray(limits)[None])),
        'optimizer_success':bool(fit.success),'optimizer_message':str(fit.message),
        'scope':'Local spline accommodation of one ion, with zero displacement and velocity at endpoints; coordination is not molecular dynamics.'}
    return result,record


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path,required=True);p.add_argument('--sources',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--frames',type=int,default=576)
    p.add_argument('--fps',type=int,default=24);p.add_argument('--mode',choices=['teaching','continuous'],default='continuous')
    a=p.parse_args()
    if a.output.exists():p.error('Use a fresh output directory.')
    if a.frames<240 or a.fps<1:p.error('Use at least 240 frames and positive fps.')
    source=json.loads((a.asset_dir/'manifest.json').read_text())
    route=json.loads((Path(__file__).resolve().parents[1]/'assets/ras-gap-route.json').read_text())
    for code,digest in route['source_sha256'].items():
        if source['sources'][code]['sha256']!=digest:raise ValueError('Route source changed; re-evaluate the pathway: '+code)
    source_path=a.sources/'1WQ1.cif'
    if hashlib.sha256(source_path.read_bytes()).hexdigest()!=route['source_sha256']['1WQ1']:
        raise ValueError('Water-coordinate source hash mismatch.')
    d={k:dict(np.load(a.asset_dir/(k+'.npz'))) for k in ['ras','gap','new_gtp','old_gdp','gap_gdp']}
    nt=d['new_gtp'];gd=d['gap_gdp'];protein=d['ras'];gap=d['gap']
    fitted,water,gamma_ids,product=chemistry(nt,gd,source_path)
    core_ids=[i for i in range(len(nt['atoms'])) if i not in gamma_ids]
    names=list(nt['atom_names']);core_names=[str(names[i]) for i in core_ids]
    target_indices=[list(d['old_gdp']['atom_names']).index(n) for n in core_names]
    source_indices=[list(gd['atom_names']).index(n) for n in core_names]
    target=d['old_gdp']['atoms'][target_indices]
    index_map={old:new for new,old in enumerate(core_ids)}
    bonds=np.array([[index_map[int(i)],index_map[int(j)]] for i,j in nt['bonds'] if i in index_map and j in index_map])
    intervals={'binding':[0,.36],'hydrolysis':[.36,.52],'release':[.58,1]} if a.mode=='continuous' else {'binding':[.06,.36],'hydrolysis':[.43,.57],'release':[.65,.95]}
    t=np.linspace(0,1,a.frames+1)
    phases={k:np.clip((t-start)/(end-start),0,1) for k,(start,end) in intervals.items()}
    bind=phases['binding'];release=phases['release'];transfer=smooth(phases['hydrolysis'])
    opening=smooth((bind-.08)/.62);relax=smooth((release-.08)/.72)
    gap_q=bind.copy();leaving=t>=intervals['release'][0];gap_q[leaving]=1-smooth(release[leaving]/.9)
    gap_local,gap_world,gap_pose,sidechains=gap_samples(gap,route,gap_q)
    ras=[];core=[];gamma=[];waters=[];metals=[]
    chemical=chemical_samples(fitted,water,gamma_ids,product,transfer)
    for f in range(len(t)):
        if leaving[f]:
            value=float(relax[f])
            r=covalent_interpolate(protein['gap'],protein['atoms'],protein['bonds'],value,initial=ras[-1] if ras else None,max_nfev=600,ftol=1e-5)
            c=covalent_interpolate(fitted[core_ids],target,bonds,value,initial=core[-1] if core else None,max_nfev=300,ftol=1e-6)
            mg=gd['magnesium']*(1-value)+d['old_gdp']['magnesium']*value
        else:
            value=float(opening[f])
            r=covalent_interpolate(protein['gtp'],protein['gap'],protein['bonds'],value,initial=ras[-1] if ras else None,max_nfev=600,ftol=1e-5)
            nucleotide=group_interpolate(nt['atoms'],fitted,value);c=nucleotide[core_ids]
            mg=nt['magnesium']*(1-value)+gd['magnesium']*value
        ras.append(r);core.append(c);metals.append(mg)
        gamma.append(nucleotide[gamma_ids] if t[f]<intervals['hydrolysis'][0] else chemical[f,:4])
        waters.append(chemical[f,4:])
        if f%48==0:print('constrained hydrolysis frame',f+1,flush=True)
    ras=np.array(ras);core=np.array(core);gamma=np.array(gamma);waters=np.array(waters);metals=np.array(metals)
    # Pi moves as one complete five-heavy-atom product after the bound intermediate.
    pi_phase=smooth((release-.22)/.56)
    pi_direction=np.array([-.951720469911,.149261341797,.268229166667])
    pi_pose=matrices(product,pi_direction,18,.5+pi_phase*.5,'loop',0,0)
    pi=apply(product,pi_pose)
    gamma[leaving]=pi[leaving,:4];waters[leaving]=pi[leaving,4:]
    report={'kind':'ras-gap-hydrolysis','asset_dir':str(a.asset_dir.resolve()),'frames':a.frames,'fps':a.fps,
            'mode':a.mode,'intervals':intervals,'passed':False,'status':'checking_all_partners',
            'threshold_nm':.20,'core_atom_names':core_names,'gamma_atom_names':['PG','O1G','O2G','O3G'],
            'water_source':{'id':'1WQ1','chain':'R','residue':230,'atom':'O'},
            'gtp_extra_translation_nm':[0,0,0],'pi_release':{'direction':pi_direction.tolist(),'phase':[.22,.78]},
            'source_contact_policy':'Only matching observed atom pairs, allowing at most 0.01 nm further approach; covalent reaction pairs listed separately.'}
    a.output.mkdir(parents=True)
    fixed=np.concatenate([ras,core,gamma,waters,metals],axis=1)
    counts,_=pair_check(gap_world,fixed)
    bad=np.flatnonzero(counts)
    if len(bad):
        if np.any(~leaving[bad]):
            report['early_gap_clashes']=(bad+1).tolist()
        else:
            margin=max(8,round((a.frames+1)*.035))
            gap_world,detour=refine(gap_world,[fixed],max(1,int(bad[0])-margin),min(a.frames-1,int(bad[-1])+margin),maximum=.4)
            report['gap_release_detour']=detour
    mg_native=merge_contacts(source_contacts(nt['magnesium'],nt['atoms'][core_ids],.2),
        source_contacts(gd['magnesium'],gd['atoms'][source_indices],.2),source_contacts(d['old_gdp']['magnesium'],target,.2))
    allowed={
        'gap_ras':source_contacts(gap['atoms'],protein['gap'],.2),
        'core_ras':merge_contacts(source_contacts(nt['atoms'][core_ids],protein['gtp'],.2),source_contacts(gd['atoms'][source_indices],protein['gap'],.2),source_contacts(target,protein['atoms'],.2)),
        'gamma_ras':source_contacts(nt['atoms'][gamma_ids],protein['gtp'],.2),
        'water_ras':source_contacts(water[None],protein['gap'],.2),
        'mg_ras':merge_contacts(source_contacts(nt['magnesium'],protein['gtp'],.2),source_contacts(gd['magnesium'],protein['gap'],.2),source_contacts(d['old_gdp']['magnesium'],protein['atoms'],.2)),
        'mg_core':mg_native,'mg_gamma':source_contacts(nt['magnesium'],nt['atoms'][gamma_ids],.2),
        'gamma_core':source_contacts(nt['atoms'][gamma_ids],nt['atoms'][core_ids],.2),
        'gamma_water':{(0,0):.155}}
    actors={'gap':gap_world,'ras':ras,'core':core,'gamma':gamma,'water':waters,'mg':metals}
    actors['mg'],metal_detour=refine_magnesium(actors,allowed);metals=actors['mg']
    if metal_detour:report['magnesium_return_detour']=metal_detour
    checks={}
    for first,second in [('gap','ras'),('core','ras'),('gamma','ras'),('water','ras'),('mg','ras'),
            ('core','gap'),('gamma','gap'),('water','gap'),('mg','gap'),('mg','core'),('mg','gamma'),
            ('mg','water'),('gamma','core'),('water','core'),('gamma','water')]:
        label=first+'_'+second
        counts,distance=pair_check(actors[first],actors[second],allowed.get(label))
        checks[label]={'maximum_overlap_atoms':int(counts.max()),'bad_frames':(np.flatnonzero(counts)+1).tolist(),'minimum_nm':float(distance.min())}
    lookup={(int(r),str(n)):i for i,(r,n) in enumerate(zip(gap['residue_numbers'],gap['atom_names']))}
    for residue,names in sidechains:
        selected=[lookup[residue,n] for n in names];other=[i for i in range(len(gap['atoms'])) if i not in selected]
        counts,dist=pair_check(gap_local[:,selected],gap_local[:,other],source_contacts(gap['atoms'][selected],gap['atoms'][other],.2))
        checks['gap_self_'+str(residue)]={'maximum_overlap_atoms':int(counts.max()),'bad_frames':(np.flatnonzero(counts)+1).tolist(),'minimum_nm':float(dist.min())}
    lengths={}
    for label,points,edges in [('ras',ras,protein['bonds']),('gap',gap_world,gap['bonds']),('core',core,bonds)]:
        values=np.linalg.norm(points[:,edges[:,0]]-points[:,edges[:,1]],axis=2)
        lengths[label]=[float(values.min()),float(values.max())]
    report.update(full_timeline_checks=checks,covalent_length_ranges_nm=lengths,
        heavy_atoms=34,atom_balance='GTP(32) + water oxygen(1) + Mg(1) -> GDP(28) + Pi(5) + Mg(1)',
        forming_covalent_pair=['GTP:PG','water:O'],leaving_covalent_pair=['GTP:PG','GTP:O3B'],
        chemistry_limits=['Hydrogens, proton transfers, protonation states, solvent and electronic energies are omitted.',
            'True GTP from 1QRA is fitted to the GDP core of the 1WQ1 transition-state analogue. No Al or F atom becomes phosphate.',
            'Pi is an idealized tetrahedral product, not an observed product crystal.',
            'Protein and ligand motion is source-guided geometry, not kinetics or molecular dynamics.'],
        passed=all(c['maximum_overlap_atoms']==0 for c in checks.values()) and all(.10<v[0] and v[1]<.22 for v in lengths.values()))
    report['status']='validated_geometry' if report['passed'] else 'failed_full_timeline'
    np.savez_compressed(a.output/'hydrolysis.npz',ras=ras,gap=gap_world,gap_local=gap_local,gap_pose=gap_pose,
        nucleotide_core=core,phosphate=gamma,water_oxygen=waters,magnesium=metals,core_bonds=bonds,
        transfer=transfer,pi_release=pi_phase)
    (a.output/'manifest.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({'passed':report['passed'],'checks':checks,'bond_ranges':lengths}),flush=True)
    if not report['passed']:raise SystemExit('Hydrolysis motion failed; keep diagnostic output for correction.')


if __name__=='__main__':main()
