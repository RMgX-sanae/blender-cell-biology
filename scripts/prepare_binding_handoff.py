"""Compose two reviewed binding poses into a source-bounded occupancy handoff.

Reuse source qualification, rigid registration, obstacle-aware release and the
molecular-sequence display contract. This is a geometrical teaching model, not
an inferred rate law or a complete initiation pathway.
"""
import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist

from prepare_binding_candidates import fit_template, residue_id, sha, transform
from prepare_molecular_asset import surface, xyz
from prepare_source_states import qualify
from rigid_escape import escape, sample


def load_endpoints(sources, config):
    qualified = {s['pdb']:qualify(sources,s) for s in config['sources']}
    reference = config['reference']; chain = config['host_chain']
    host_rows = qualified[reference][2][chain]; origin = xyz(host_rows).mean(0)
    points = {'host':(xyz(host_rows)-origin)*.1}
    identities = {'host':[f'{reference}:{chain}:{residue_id(r)}:{r[4]}:{r[2]}' for r in host_rows]}
    fits = {}; provenance = {}
    for code, data in qualified.items():
        provenance[code] = {'source':data[1], 'chains':data[3], 'assembly':data[5]}
        if code == reference: matrix = np.eye(4)
        else:
            matrix, fits[code] = fit_template(data[2][chain], host_rows, config['alignments'][code])
        for guest in config['guests']:
            if guest['source'] != code: continue
            rows = data[2][guest['chain']]; key = guest['id']
            if key in points: raise ValueError('Duplicate component name')
            points[key] = (transform(xyz(rows),matrix)-origin)*.1
            identities[key] = [f'{code}:{guest["chain"]}:{residue_id(r)}:{r[4]}:{r[2]}' for r in rows]
        if code != reference and any(data[4].values()):
            raise ValueError('Qualify moving-companion ownership before using this handoff input')
    ligand_records = []
    for chain_name, residue, atoms in qualified[reference][4]['anchor']:
        key = f'ligand_{chain_name}_{residue.name}_{residue.seqid}'
        points[key] = (np.array([list(a.pos) for a in atoms])-origin)*.1
        identities[key] = [f'{reference}:{chain_name}:{residue.seqid}:{residue.name}:{a.name}' for a in atoms]
        ligand_records.append({'component':key,'source':reference,'author_chain':chain_name,
            'author_residue':str(residue.seqid),'ccd':residue.name,'atom_names':[a.name for a in atoms],
            'elements':[a.element.name for a in atoms], 'owner':'host'})
    if qualified[reference][4]['moving']:
        raise ValueError('This input supports fixed host companions; qualify moving companions separately')
    return points, identities, {'origin_angstrom':origin.tolist(), 'source_states':provenance,
                                'fits':fits,'host_companions':ligand_records}


def fractions(time, mode):
    def ease(q):
        q = np.clip(q,0,1); return q*q*q*(10+q*(-15+6*q))
    if mode == 'loop':
        first = ease(time/.25)-ease((time-.75)/.25)
        second = 1-ease((time-.25)/.25)+ease((time-.5)/.25)
        windows = [[0,.25],[.25,.5],[.5,.75],[.75,1]]
    elif mode == 'teaching':
        first = ease((time-.1)/.2)-ease((time-.8)/.2)
        second = 1-ease((time-.3)/.2)+ease((time-.6)/.2)
        windows = [[0,.1],[.1,.3],[.3,.5],[.5,.6],[.6,.8],[.8,1]]
    else: raise ValueError('Use loop or teaching timing')
    return [np.clip(first,0,1),np.clip(second,0,1)], windows


def screen(tracks, guard):
    checks = {}
    for left,right in combinations(tracks,2):
        distances = [float(cKDTree(a).query(b)[0].min()) for a,b in zip(tracks[left],tracks[right])]
        checks[left+'--'+right] = {'minimum_nm':min(distances),'closest_frame':int(np.argmin(distances))+1,
            'violating_frames':[i+1 for i,d in enumerate(distances) if d < guard]}
    return {'passed':not any(r['violating_frames'] for r in checks.values()), 'guard_nm':guard,
            'pairs':checks, 'scope':'All retained heavy atoms of all components at every output pose; no source-contact exemptions. Not a van der Waals or continuous-time proof.'}


def prepare(sources, config_path, output, frames, fps, mode):
    if output.exists(): raise ValueError('Use a fresh output directory')
    config = json.loads(config_path.read_text()); guests = config['guests']
    if len(guests) != 2 or not config.get('competition_evidence') or not config.get('scope'):
        raise ValueError('Supply two reviewed guest poses and evidence for the occupancy example')
    if frames < 24 or frames%4 or fps <= 0: raise ValueError('Use a frame count divisible by four and positive fps')
    guard = config['guard_nm']; points, identities, provenance = load_endpoints(sources,config)
    names = [g['id'] for g in guests]; fixed_names = [k for k in points if k not in names]
    fixed = np.concatenate([points[k] for k in fixed_names])
    overlap = cdist(points[names[0]],points[names[1]])
    if overlap.min() >= guard:
        raise ValueError('These specific binding poses do not establish the declared geometric exclusion')
    paths = {}; reports = {}; obstacles = fixed
    for guest in guests:
        key = guest['id']
        paths[key], reports[key] = escape(points[key],obstacles,guest['departure_direction'],guard,
            config['minimum_travel_nm'],config['park_clearance_nm'])
        # The second moves while the first is parked; final screening also checks
        # the first moving around the second's actual final parking pose.
        obstacles = np.concatenate([fixed,sample(points[key],paths[key],[1])[0]])
    time = np.arange(frames+1)/frames; amounts, windows = fractions(time,mode)
    tracks = {key:np.repeat(points[key][None],frames+1,axis=0) for key in fixed_names}
    for key, amount in zip(names,amounts): tracks[key] = sample(points[key],paths[key],amount)
    checks = screen(tracks,guard)
    if not checks['passed']: raise ValueError('Simultaneous-partner screening failed: '+json.dumps(checks))
    if max(float(np.max(abs(data[0]-data[-1]))) for data in tracks.values()) > 1e-9:
        raise ValueError('The declared cyclic handoff does not return the same atoms to the same poses')
    if len(config['captions'][mode]) != len(windows): raise ValueError('Caption count must match the selected timing')
    output.mkdir(parents=True)
    np.savez_compressed(output/'motion.npz',**tracks)
    np.savez_compressed(output/'paths.npz',**{key+'_'+field:array for key,path in paths.items() for field,array in path.items()})
    vertices, faces = surface(points['host']*10,.75,1.15,.55)
    np.savez_compressed(output/'background_host.npz',vertices=vertices*.1,faces=faces)
    captions = [[begin,end,text] for (begin,end),text in zip(windows,config['captions'][mode])]
    record = {'kind':'binding-occupancy-handoff','passed':True,'checks':checks,'units':'nm','frames':frames,'fps':fps,
        'mode':mode,'scope':config['scope'],'config':config,'config_sha256':sha(config_path),
        'atom_ids':identities,'provenance':provenance,'paths':reports,'captions':captions,
        'simultaneous_bound_poses':{'minimum_nm':float(overlap.min()),'atom_pairs_below_guard':int(np.count_nonzero(overlap<guard)),
            'interpretation':'The two aligned complete bound poses overlap; this does not exclude other partial-binding conformations or establish a kinetic mechanism.'},
        'data_sha256':{name:sha(output/name) for name in ['motion.npz','paths.npz','background_host.npz']},
        'limits':'Source-bound fragments remain rigid during illustrative departure/arrival. Registered alternative pose is inferred. No unfolding, phosphorylation, RNA recruitment, kinetic competition, affinity or full translation initiation is modeled.'}
    (output/'manifest.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'output':str(output),'frames':frames,'atoms':{k:len(v) for k,v in identities.items()},'passed':True}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('sources','config','output'): parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--frames',type=int,default=288); parser.add_argument('--fps',type=int,default=24)
    parser.add_argument('--mode',choices=['loop','teaching'],default='loop')
    args = parser.parse_args()
    prepare(args.sources,args.config,args.output,args.frames,args.fps,args.mode)
