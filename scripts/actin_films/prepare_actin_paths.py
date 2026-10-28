"""Screen intact regulatory-fragment teaching paths from their source geometry."""
import hashlib
import itertools
import json

import numpy as np
from scipy.interpolate import BSpline
from scipy.spatial import cKDTree

from prepare_contact_repairs import repair_path
from project_context import ROOT


def smooth(value):
    value = np.clip(value,0,1)
    return value**3*(10+value*(-15+6*value))


def displacement(fraction, record):
    offset = np.array(record['offset']);side = np.cross(offset,[0,0,1])
    side *= 1.2/np.linalg.norm(side)
    far = max(0,(fraction-.5)*2)
    shift = offset*fraction+side*np.sin(np.pi*far)**2+np.array(record['far_shift'])*smooth(far)
    detour = record['detour']
    if detour and detour['first_fraction']<fraction<detour['last_fraction']:
        u = (fraction-detour['first_fraction'])/(detour['last_fraction']-detour['first_fraction'])
        shift += BSpline(detour['knots'],detour['controls_nm'],3)(u)
    return shift


def screen_translations(atoms, trajectory, samples=2001):
    """Check complete chains together; trajectory omits explicitly hidden chains."""
    trees = {c:cKDTree(p) for c,p in atoms.items()}
    bounds = {c:(p.min(0),p.max(0)) for c,p in atoms.items()}
    minima = {};worst = {};previous = {};radius = .5
    for index,fraction in enumerate(np.linspace(0,1,samples)):
        shifts = trajectory(fraction)
        for a,b in itertools.combinations(shifts,2):
            pair = a+'/'+b;relative = tuple(np.array(shifts[a])-shifts[b])
            if pair in previous and relative==previous[pair][0]:continue
            low,high = bounds[b];points = atoms[a]+relative
            selected = ((points>=low-radius)&(points<=high+radius)).all(1)
            value = min(radius,float(trees[b].query(points[selected],distance_upper_bound=radius)[0].min())) if selected.any() else radius
            previous[pair] = (relative,value)
            if value<minima.get(pair,float('inf')):minima[pair] = value;worst[pair] = float(fraction)
        if index%500==0:print('TRANSLATION_SCREEN',index,'of',samples-1,flush=True)
    return {'samples':samples,'threshold_nm':.20,'search_radius_nm':radius,
        'pair_minima_nm_capped_at_radius':minima,'worst_fractions':worst,
        'passed':min(minima.values())>=.20,
        'limits':'Source heavy atom centers and authored whole-chain translations. Hidden chains omitted until their recorded reveal. No force-field or kinetics claim.'}


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
