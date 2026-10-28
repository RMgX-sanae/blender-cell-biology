"""Screen intact regulatory-fragment teaching paths from their source geometry."""
import hashlib
import itertools
import json

import numpy as np
from scipy.interpolate import BSpline
from scipy.spatial import cKDTree

from prepare_contact_repairs import repair_path
from project_context import ROOT


def displacement(fraction, record):
    offset = np.array(record['offset']);side = np.cross(offset,[0,0,1])
    side *= 1.2/np.linalg.norm(side)
    far = max(0,(fraction-.5)*2)
    shift = offset*fraction+side*np.sin(np.pi*far)**2
    detour = record['detour']
    if detour and detour['first_fraction']<fraction<detour['last_fraction']:
        u = (fraction-detour['first_fraction'])/(detour['last_fraction']-detour['first_fraction'])
        shift += BSpline(detour['knots'],detour['controls_nm'],3)(u)
    return shift


def main():
    assets = ROOT/'assets';output = assets/'actin-path-repairs.json'
    record = {'passed':False,'motions':{},
        'source':{'entry':'2BAP','url':'https://www.rcsb.org/structure/2BAP',
            'sha256':hashlib.sha256((ROOT/'sources/2BAP.cif').read_bytes()).hexdigest()},
        'limits':'Complete source regulatory fragments and local translation corrections. The DAD release is an authored teaching path, not a measured trajectory or an atomic RhoA displacement mechanism.'}
    output.write_text(json.dumps(record,indent=2)+'\n')
    atoms = {c:np.load(assets/'2BAP'/(c+'.npz'))['atoms'].astype(float) for c in 'ABCD'}
    plans = json.loads((assets/'approach-plans.json').read_text())
    for c in 'CD':
        obstacles = np.concatenate([p for k,p in atoms.items() if k!=c])
        record['motions']['dad_'+c] = repair_path(atoms[c],obstacles,plans['dad_'+c]['offset'])
    minima = {a+'/'+b:float('inf') for a,b in itertools.combinations(atoms,2)}
    for fraction in np.linspace(0,1,2001):
        state = {c:p+(displacement(fraction,record['motions']['dad_'+c]) if c in 'CD' else 0) for c,p in atoms.items()}
        trees = {c:cKDTree(p) for c,p in state.items()}
        for a,b in itertools.combinations(atoms,2):
            value = float(trees[b].query(state[a])[0].min())
            minima[a+'/'+b] = min(minima[a+'/'+b],value)
    record['screen'] = {'samples':2001,'threshold_nm':.20,'pair_minima_nm':minima,
        'scope':'Both DAD fragments move together; the DID/dimerization source fragments stay fixed. Heavy atom centers only, not molecular dynamics.'}
    record['passed'] = min(minima.values())>=.20
    output.write_text(json.dumps(record,indent=2)+'\n')
    if not record['passed']:raise ValueError('Regulatory-fragment paths still overlap; inspect actin-path-repairs.json')
    print('REGULATORY_PATHS_READY',min(minima.values()),flush=True)


if __name__=='__main__':main()
