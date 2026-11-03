"""Prepare a finite PTP1B two-step net reaction with persistent heavy-atom identities."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np

from ptp_geometry import PtpGeometry
from interaction_motion import matrices
from prepare_phosphotransfer_motion import smooth
from event_timing import phase_progress
from reaction_ledger import validate_sequence
from prepare_molecular_asset import atom_table,surface,xyz

ORDER=['transfer','peptide_release','q_turn','water_entry','hydrolysis']


def captions(config):
    timing=config['timing'];texts=config['caption_texts']
    boundaries=[0,timing['transfer'][0],timing['transfer'][1],timing['peptide_release'][0],
                timing['q_turn'][0],timing['water_entry'][0],timing['hydrolysis'][0],timing['hydrolysis'][1],1]
    if len(texts)!=8:raise ValueError('Supply eight PTP1B stage captions')
    return [[start,end,body] for start,end,body in zip(boundaries,boundaries[1:],texts) if end>start]


def sample(model,config,frames):
    time=np.arange(frames+1)/frames
    phases=phase_progress(time,config['timing'],before=zip(ORDER,ORDER[1:]))
    progress={key:smooth(value) for key,value in phases.items()}
    chemical=model.chemistry(progress['transfer'],progress['hydrolysis'])
    enzyme=np.array([model.q_pose(model.q_angles*q) for q in progress['q_turn']])
    def move(key,points,entry=False):
        route=config['routes'][key]
        fractions=progress[key]*.5+(0 if entry else .5)
        pose=matrices(points*.1,route['direction'],route['distance_nm'],fractions,'loop',route['turn_radians'],route['bend_fraction'])
        return np.einsum('fij,nj->fni',pose[:,:3,:3],points)+pose[:,None,:3,3]*10
    peptide=move('peptide_release',model.source['peptide'])
    water=move('water_entry',model.water[None],True)
    phosphate=chemical[:,:4]
    return {'enzyme':enzyme,'peptide':peptide,'phosphoryl':phosphate,'water':water},progress


def audit(model,data,progress):
    frames=len(data['enzyme']);errors=[];minima={};bond_error=0.;angle_error=0.
    for key in ['enzyme','peptide']:
        pairs=set(map(tuple,model.bonds[key]))
        angle_pairs={tuple(sorted((a,b))) for neighbors in model.neighbors[key] for a in neighbors for b in neighbors if a!=b}
        for edges,is_angle in [(pairs,False),(angle_pairs,True)]:
            edges=np.array(sorted(edges));source=model.source[key]
            expected=np.linalg.norm(source[edges[:,0]]-source[edges[:,1]],axis=1)
            actual=np.linalg.norm(data[key][:,edges[:,0]]-data[key][:,edges[:,1]],axis=2)
            error=float(np.max(abs(actual-expected)))
            if is_angle:angle_error=max(angle_error,error)
            else:bond_error=max(bond_error,error)
    sg=model.indices['enzyme'][215,'SG'];oh=model.indices['peptide'][992,'OH']
    # Static source atoms remain intact; only the four qualified Q262 side-chain atoms move internally.
    fixed=np.array([i for i in range(len(model.rows['enzyme'])) if i not in model.q_atoms])
    fixed_error=float(np.max(abs(data['enzyme'][:,fixed]-model.source['enzyme'][fixed])))
    self_excluded={}
    for i in model.q_atoms:
        neighbors=model.neighbors['enzyme'][i]
        self_excluded[int(i)]={int(i),*neighbors,*[j for neighbor in neighbors for j in model.neighbors['enzyme'][neighbor]]}
    pair_specs=[('enzyme','peptide'),('enzyme','phosphoryl'),('enzyme','water'),
                ('peptide','phosphoryl'),('peptide','water'),('phosphoryl','water')]
    flagged=0
    for frame in range(frames):
        for first,second in pair_specs:
            distances=np.linalg.norm(data[first][frame,:,None]-data[second][frame,None],axis=2)
            # Exempt only the explicit forming/breaking direct bonds, not arbitrary source contacts.
            if (first,second)==('enzyme','phosphoryl'):distances[sg,0]=np.inf
            if (first,second)==('peptide','phosphoryl'):distances[oh,0]=np.inf
            if (first,second)==('phosphoryl','water') and progress['hydrolysis'][frame]>0:distances[0,0]=np.inf
            key=first+'/'+second;value=float(distances.min());minima[key]=min(minima.get(key,float('inf')),value)
            bad=np.argwhere(distances<2.38);flagged+=len(bad)
            if len(bad) and len(errors)<30:
                i,j=bad[np.argmin([distances[i,j] for i,j in bad])]
                errors.append([frame+1,key,int(i),int(j),float(distances[i,j])])
        positions=data['enzyme'][frame]
        self_distances=np.linalg.norm(positions[model.q_atoms,None]-positions[None],axis=2)
        for row,i in enumerate(model.q_atoms):self_distances[row,list(self_excluded[int(i)])]=np.inf
        value=float(self_distances.min());minima['mobile_Q262_self']=min(minima.get('mobile_Q262_self',float('inf')),value)
        if value<2.38:
            flagged+=int(np.count_nonzero(self_distances<2.38))
            if len(errors)<30:errors.append([frame+1,'mobile_Q262_self',value])
    po=np.linalg.norm(data['phosphoryl'][:,1:]-data['phosphoryl'][:,:1],axis=2)
    oo=np.array([np.linalg.norm(data['phosphoryl'][:,i]-data['phosphoryl'][:,j],axis=1) for i,j in [(1,2),(1,3),(2,3)]])
    ranges={'P_O_A':[float(po.min()),float(po.max())],'nonbridging_O_O_A':[float(oo.min()),float(oo.max())]}
    geometric=po.min()>1.4 and po.max()<1.7 and oo.min()>2.2
    finite=all(np.isfinite(value).all() for value in data.values())
    return {'passed':finite and not flagged and geometric and max(bond_error,angle_error,fixed_error)<1e-7,
            'frames_checked':frames,'minimum_nonbonded_A':minima,'flagged_pairs':flagged,'first_failures':errors,
            'maximum_polymer_bond_length_change_A':bond_error,'maximum_polymer_adjacent_angle_pair_change_A':angle_error,
            'fixed_enzyme_atoms':len(fixed),'maximum_fixed_enzyme_change_A':fixed_error,'phosphoryl_geometry':ranges,
            'contact_policy':'2.38 A heavy-atom-center screen. Only explicit P-SG, P-TyrOH and forming P-water bonds are excluded; Q262 self 1-2/1-3 covalent neighbors excluded. No grandfathered novel collisions.',
            'limits':'Geometric teaching constraints, not energy, proton relay, measured dynamics or biological rate.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources',type=Path,required=True);parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--frames',type=int,default=720);parser.add_argument('--fps',type=int,default=24)
    args=parser.parse_args()
    if args.output.exists():parser.error('Use a fresh preparation directory')
    if args.frames<96 or args.fps<1:parser.error('Use at least 96 frames and positive fps')
    model=PtpGeometry(args.sources);config=json.loads(args.config.read_text())
    if config['version']!=1 or config['kind']!='ptp1b-two-step-net-reaction':parser.error('Use qualified PTP1B controls')
    data,progress=sample(model,config,args.frames);checks=audit(model,data,progress)
    if not checks['passed']:
        print(json.dumps(checks,indent=2));raise SystemExit('Combined PTP1B geometry failed; no assets written')
    sequence=model.ledger();book=validate_sequence(sequence)
    background={'ptp1b':model.source['enzyme']};background_records={}
    for name,code,expected in [('ubiquitin','1UBQ','056f98710cb2b36f633c45e41902a02eb446e82871da21ff2dd44f74a56ca0f6'),
                               ('mdm2','1YCR','03f62994ea1ff55d98c0d9d835b706319e523367fce56d9deb51b850c83e47a8')]:
        path=args.sources/(code+'.cif');record=json.loads(path.with_suffix('.cif.source.json').read_text())
        if record['id'].upper()!=code or record['sha256']!=expected or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError('Background source changed; requalify '+code)
        rows=atom_table(path)['A'];background[name]=xyz(rows)
        background_records[name]={'source_id':code,'author_chain':'A','sha256':expected,'selected_heavy_atoms':len(rows),
                                  'entry_url':record['entry_url'],'scope':'Human source construct used as an illustrative crowding shape, not measured composition or a PTP1B binding partner.'}
    args.output.mkdir(parents=True)
    np.savez_compressed(args.output/'motion.npz',**{k:v*.1 for k,v in data.items()},**{k+'_fraction':v for k,v in progress.items()})
    atom_ids={key:[f'{key}:{r[3]}:{r[2]}' for r in rows] for key,rows in model.rows.items()}
    elements={key:[r[1] for r in rows] for key,rows in model.rows.items()}
    atom_ids.update(phosphoryl=['phosphoryl:'+n for n in ['P','O1','O2','O3']],water=['water:O'])
    elements.update(phosphoryl=['P','O','O','O'],water=['O'])
    for key,rows in model.rows.items():
        np.savez_compressed(args.output/(key+'.npz'),atoms=model.source[key]*.1,bonds=model.bonds[key],
                            atom_ids=np.array(atom_ids[key]),elements=np.array(elements[key]))
    # Display guides only: preserve every source point, draw the whole C-alpha
    # trace and the complete selected active-site residues with distinct semantics.
    ca={r[3]:i for i,r in enumerate(model.rows['enzyme']) if r[2]=='CA'}
    context={(ca[n],ca[n+1]) for n in ca if n+1 in ca}
    active={181,215,221,262}
    context.update(tuple(sorted((int(a),int(b)))) for a,b in model.bonds['enzyme']
                   if model.rows['enzyme'][a][3] in active and model.rows['enzyme'][b][3] in active)
    np.savez_compressed(args.output/'enzyme_context.npz',bonds=np.array(sorted(context),dtype=np.int32))
    (args.output/'reaction-sequence.json').write_text(json.dumps(sequence,indent=2)+'\n')
    for name,points in background.items():
        vertices,faces=surface(points-points.mean(0),.9,1.15,.55)
        np.savez_compressed(args.output/('background_'+name+'.npz'),vertices=vertices*.1,faces=faces)
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.output.iterdir() if p.is_file()}
    manifest={'kind':'ptp1b-two-step-net-reaction','passed':True,'frames':args.frames,'fps':args.fps,'units':'nm',
              'atom_ids':atom_ids,'elements':elements,'data_sha256':hashes,'reaction_sequence':book,'captions':captions(config),
              'reaction_atom_indices':{'SG':model.indices['enzyme'][215,'SG'],'TyrOH':model.indices['peptide'][992,'OH'],
                                       'Q262':model.q_atoms.tolist()},'source_record':model.record,'checks':checks,'config':config,
              'scope':config['scope'],'background_sources':background_records,
              'presentation':'Finite continuous motion with observation pauses; not a seamless catalytic loop.'}
    (args.output/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'passed':True,'heavy_atoms':book['mapped_atoms'],'frames_checked':checks['frames_checked'],'output':str(args.output)}))


if __name__=='__main__':main()
