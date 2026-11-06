"""Independently compare a prepared occupancy handoff with its raw mmCIF atoms.

Does not call the preparation, registration or motion planner. Checks retained
identity, rigid source geometry, fixed companions and every simultaneous pair.
"""
import argparse
import hashlib
from itertools import combinations
import json
from pathlib import Path

import gemmi
import numpy as np
from scipy.spatial import cKDTree


def require(condition, message):
    if not condition: raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def register(first, second):
    a = first-first.mean(0); b = second-second.mean(0)
    u, _, vt = np.linalg.svd(a.T@b)
    rotation = u@np.diag([1,1,np.linalg.det(u@vt)])@vt
    return rotation, second.mean(0)-first.mean(0)@rotation


def audit(sources, config_path, assets):
    config = json.loads(config_path.read_text()); record = json.loads((assets/'manifest.json').read_text())
    require(record['config'] == config and record['config_sha256'] == digest(config_path), 'Config changed')
    require(record['kind'] == 'binding-occupancy-handoff' and record['units'] == 'nm', 'Wrong asset kind/units')
    for name, expected in record['data_sha256'].items():
        require(digest(assets/name) == expected, 'Prepared data changed: '+name)
    models = {}; atoms = {}; chains = {}
    for source in config['sources']:
        code = source['pdb']; path = sources/(code+'.cif')
        require(digest(path) == source['sha256'], 'Raw source changed: '+code)
        models[code] = gemmi.read_structure(str(path))[0]
        for chain in models[code]:
            polymer = []
            for residue in chain:
                for atom in residue:
                    if atom.element.is_hydrogen or atom.occ <= 0 or atom.altloc not in ('\x00','A'): continue
                    identity = f'{code}:{chain.name}:{residue.seqid}:{residue.name}:{atom.name}'
                    require(identity not in atoms, 'Duplicate retained source identity')
                    atoms[identity] = np.array(list(atom.pos))
                    if residue.entity_type == gemmi.EntityType.Polymer: polymer.append(identity)
            chains[(code,chain.name)] = polymer
    reference = config['reference']; host_chain = config['host_chain']
    host_ids = chains[(reference,host_chain)]
    raw_host = np.array([atoms[i] for i in host_ids]); origin = raw_host.mean(0)
    expected = {'host':(raw_host-origin)*.1}; identities = {'host':host_ids}; fits = {}
    for guest in config['guests']:
        code = guest['source']; key = guest['id']; identities[key] = chains[(code,guest['chain'])]
        guest_points = np.array([atoms[i] for i in identities[key]])
        if code != reference:
            mapping = config['alignments'][code]
            def ca(entry, number):
                found = [i for i in chains[(entry,host_chain)] if i.split(':')[2] == str(number) and i.endswith(':CA')]
                require(len(found) == 1, 'Ambiguous source CA map'); return found[0]
            pairs = [(ca(code,a),ca(reference,b)) for a,b in mapping['residue_pairs']]
            differences = [{'source_residue':a,'template_residue':b,'source_name':x.split(':')[3],
                            'template_name':y.split(':')[3]}
                           for (a,b),(x,y) in zip(mapping['residue_pairs'],pairs) if x.split(':')[3] != y.split(':')[3]]
            require(differences == mapping['expected_substitutions'], 'Unreviewed sequence substitution')
            first = np.array([atoms[a] for a,b in pairs]); second = np.array([atoms[b] for a,b in pairs])
            rotation, shift = register(first,second)
            rmsd = float(np.sqrt(np.mean(np.sum((first@rotation+shift-second)**2,axis=1))))
            require(rmsd <= mapping['max_CA_rmsd_angstrom'], 'Independent registration failed')
            guest_points = guest_points@rotation+shift; fits[code] = {'CA_pairs':len(pairs),'rmsd_angstrom':rmsd}
        expected[key] = (guest_points-origin)*.1
    source = next(s for s in config['sources'] if s['pdb'] == reference)
    for ligand in source['ligand_instances']:
        number = str(ligand['residue_number'])+ligand.get('insertion_code','')
        key = f'ligand_{ligand["chain"]}_{ligand["component"]}_{number}'
        prefix = f'{reference}:{ligand["chain"]}:{number}:{ligand["component"]}:'
        identities[key] = [i for i in atoms if i.startswith(prefix)]
        require(identities[key] and ligand['owner_role'] == 'anchor', 'Missing or unqualified companion')
        expected[key] = (np.array([atoms[i] for i in identities[key]])-origin)*.1
    require(record['atom_ids'] == identities, 'Retained source atom identity/order differs')
    data = dict(np.load(assets/'motion.npz')); frames = record['frames']
    require(set(data) == set(expected), 'Missing or extra component')
    guests = [g['id'] for g in config['guests']]; rigidity = {}; subframe_rigidity = {}
    for key, coordinates in data.items():
        require(coordinates.shape == (frames+1,len(identities[key]),3) and np.isfinite(coordinates).all(), 'Invalid coordinate inventory')
        require(np.max(abs(coordinates[0]-coordinates[-1])) < 1e-8, 'Loop does not close')
        if key not in guests:
            require(np.max(abs(coordinates-expected[key])) < 1e-8, 'Fixed source/companion moved: '+key)
        else:
            bound_frame = 0 if key == guests[0] else frames//2
            require(np.max(abs(coordinates[bound_frame]-expected[key])) < 1e-8, 'Bound source pose differs: '+key)
        errors = []
        for pose in coordinates:
            rotation, shift = register(expected[key],pose)
            errors.append(float(np.max(np.linalg.norm(expected[key]@rotation+shift-pose,axis=1))))
        rigidity[key] = max(errors); require(rigidity[key] < 1e-7, 'Retained body deformed: '+key)
        errors = []
        for pose in (coordinates[:-1]+coordinates[1:])/2:
            rotation, shift = register(expected[key],pose)
            errors.append(float(np.max(np.linalg.norm(expected[key]@rotation+shift-pose,axis=1))))
        subframe_rigidity[key] = max(errors)
        require(subframe_rigidity[key] < .005, 'Linear display interpolation deforms a body excessively')
    minima = {}; guard = config['guard_nm']
    for left,right in combinations(data,2):
        minimum = float('inf')
        for fraction in (0,.25,.5,.75,1):
            a = data[left][:-1]*(1-fraction)+data[left][1:]*fraction
            b = data[right][:-1]*(1-fraction)+data[right][1:]*fraction
            minimum = min(minimum,min(float(cKDTree(x).query(y)[0].min()) for x,y in zip(a,b)))
        require(minimum >= guard, 'Simultaneous component overlap: '+left+'--'+right)
        minima[left+'--'+right] = minimum
    return {'passed':True,'mode':record['mode'],'poses':frames+1,'raw_source_fits':fits,
            'atom_counts':{k:len(v) for k,v in identities.items()},'maximum_rigid_fit_error_nm':rigidity,
            'linear_midframe_rigid_fit_error_nm':subframe_rigidity,'guard_nm':guard,'pair_minimum_nm':minima,
            'scope':'Independent raw source identities/registration, whole-body geometry, bound endpoints, fixed host/cap, closed loop and all simultaneous pairs at integer plus quarter-frame samples.',
            'limits':'Display interpolation has the measured small chord error; no exact continuous-time collision proof, energetic docking score, kinetics or biochemical displacement claim.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('sources','config','assets','output'): parser.add_argument('--'+name,type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists(): parser.error('Use a fresh audit report')
    result = audit(args.sources,args.config,args.assets)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'passed':True,'poses':result['poses'],'minimum_nm':min(result['pair_minimum_nm'].values())}))
