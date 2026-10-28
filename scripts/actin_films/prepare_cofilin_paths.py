"""Prepare a whole-subunit cofilactin boundary and illustrative disassembly.

6UC4 is a 9.2 A barbed-side boundary, not a severing trajectory. E/H and F/G
are complementary actin domains in the deposited model, not four monomers.
"""
import hashlib
import itertools
import json
import sys

import numpy as np
from scipy.interpolate import BSpline
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from prepare_actin_paths import displacement,smooth
from prepare_contact_repairs import fit_placement,repair_path,repair_rigid_release,transform
from plan_clear_approaches import choose
from project_context import ROOT,SKILL
sys.path.insert(0,str(SKILL/'scripts'))
from prepare_molecular_asset import atom_table,surface,xyz


def rigid_pose(fraction,record):
    shift = displacement(fraction,{**record,'detour':None})
    detour = record['detour'];q = np.zeros(6)
    if detour['first_fraction']<fraction<detour['last_fraction']:
        u = (fraction-detour['first_fraction'])/(detour['last_fraction']-detour['first_fraction'])
        q = BSpline(detour['knots'],detour['controls_nm_radians'],3)(u)
    pivot = np.array(detour['pivot_nm'])
    matrix = np.eye(4);matrix[:3,:3] = Rotation.from_rotvec(q[3:]).as_matrix()
    matrix[:3,3] = shift+q[:3]+pivot-matrix[:3,:3]@pivot
    return matrix


def main():
    assets = ROOT/'assets';output = assets/'cofilin-path-repairs.json'
    source = ROOT/'sources/6UC4.cif'
    provenance = json.loads(source.with_suffix('.cif.source.json').read_text())
    assert hashlib.sha256(source.read_bytes()).hexdigest()==provenance['sha256']
    units = {c:[c] for c in ['D','J','K','L','A','C','B','I','M','N','O','P']}
    units.update(EH=['E','H'],FG=['F','G'])
    cofilins = list('IMNOP')
    groups = {'left':['D','J','K','L','EH','FG']+cofilins,
              'right':['A','C','B'],'departing':['FG','P']}
    record = {'passed':False,'motions':{},'placements':{},'units':units,'groups':groups,
        'intervals':{'sever':[.36,.51],'monomer_release':[.66,.84]},
        'source':{'entry':'6UC4','url':'https://www.rcsb.org/structure/6UC4',
            'sha256':provenance['sha256'],'resolution_angstrom':9.2,
            'paper_doi':'10.1073/pnas.1915987117'},
        'limits':'The source is the barbed-end side of a cluster. It illustrates one possible boundary, not the preferentially disrupted pointed-side geometry or event frequencies. Early binding, severing and terminal release are authored paths. Actin conformations remain those of this partially decorated source; no naked-to-decorated transition or F-to-G conversion is measured. Whole cofilins receive recorded bounded rigid fits to remove low-resolution atomic-center overlaps; no atom is deleted or molecule sliced.'}
    output.write_text(json.dumps(record,indent=2)+'\n')
    table = atom_table(source)
    assert sorted(c for names in units.values() for c in names)==sorted(table)
    rows = {u:[r for c in chains for r in table[c]] for u,chains in units.items()}
    for u in ['EH','FG']:
        identities = [(r[3],r[2]) for r in rows[u]]
        assert len(identities)==len(set(identities)), 'Overlapping source domains: '+u
    origin = xyz(rows['FG']).mean(0)
    x = xyz(rows['B']).mean(0)-xyz(rows['D']).mean(0);x /= np.linalg.norm(x)
    z = xyz(rows['P']).mean(0)-origin;z -= x*np.dot(z,x);z /= np.linalg.norm(z)
    basis = np.array([x,np.cross(z,x),z]).T
    folder = assets/'cofilin_boundary';folder.mkdir(exist_ok=True)
    atoms = {};ownership = []
    for u,selected in rows.items():
        coords = (xyz(selected)-origin)@basis
        vertices,faces = surface(coords,1.15,1.25,.55)
        path = folder/(u+'.npz')
        np.savez_compressed(path,vertices=vertices*.1,faces=faces,atoms=coords.astype(np.float32)*.1,
            atom_names=np.array([r[2] for r in selected]),residue_numbers=np.array([r[3] for r in selected]),
            source_chains=np.array([c for c in units[u] for _ in table[c]]))
        atoms[u] = np.load(path)['atoms'].astype(float)
        ownership.append({'unit':u,'auth_chains':units[u],'atoms':len(selected),
            'vertices':len(vertices),'faces':len(faces),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    record['geometry'] = {'units':ownership,'polymer_atom_count':sum(len(v) for v in atoms.values()),
        'source_polymer_atom_count':sum(len(v) for v in table.values()),
        'origin_angstrom':origin.tolist(),'basis_columns':basis.tolist(),
        'representation':'Gaussian envelopes of all modeled polymer heavy atoms. ADP/Mg and waters are not displayed in this overview; nucleotide chemistry is a separate chapter. Source gaps are not reconstructed.'}
    original = dict(atoms)
    for c in cofilins:
        fixed = np.concatenate([p for k,p in atoms.items() if k!=c])
        matrix,fit = fit_placement(atoms[c],fixed)
        record['placements'][c] = fit
        atoms[c] = transform(atoms[c],matrix)
        print('BOUNDARY_COFILIN_FIT',c,fit['maximum_atom_displacement_nm'],flush=True)
    combine = lambda names:np.concatenate([atoms[c] for c in names])
    for i,c in enumerate(cofilins):
        key = 'cofilin_'+c
        fixed = combine([k for k in atoms if k!=c])
        center = atoms[c].mean(0)
        plan = choose(atoms[c],fixed,[0,center[1],center[2]],length=4,clearance=.21)
        record['motions'][key] = repair_path(atoms[c],fixed,plan['offset'])
        record['intervals'][key] = [.08+i*.022,.19+i*.022]
        print('COFILIN_PATH_READY',key,record['motions'][key]['after_minimum_nm'],flush=True)
    output.write_text(json.dumps(record,indent=2)+'\n')
    relative = np.array([18.,0,0])
    record['motions']['sever'] = repair_rigid_release(combine(groups['right']),combine(groups['left']),
        atoms['A'].mean(0),relative,samples=81,iterations=800)
    print('BOUNDARY_RELEASE_READY',record['motions']['sever']['after_minimum_nm'],flush=True)
    remaining = [k for k in groups['left'] if k not in groups['departing']]
    fixed = np.concatenate([combine(remaining),combine(groups['right'])+relative])
    plan = choose(combine(groups['departing']),fixed,[1,-.25,.1],length=9,clearance=.21)
    record['motions']['monomer_release'] = repair_path(combine(groups['departing']),fixed,plan['offset'])

    def progress(fraction,key):
        first,last = record['intervals'][key]
        return smooth((fraction-first)/(last-first))

    def trajectory(fraction):
        relative_pose = rigid_pose(progress(fraction,'sever'),record['motions']['sever'])
        left = np.eye(4);left[:3,3] = -.5*relative_pose[:3,3]
        right = left@relative_pose
        depart = displacement(progress(fraction,'monomer_release'),record['motions']['monomer_release'])
        poses = {}
        for c in atoms:
            local = np.eye(4)
            if c in groups['departing']:local[:3,3] += depart
            if c in cofilins:
                local[:3,3] += displacement(1-progress(fraction,'cofilin_'+c),record['motions']['cofilin_'+c])
                local = local@np.array(record['placements'][c]['correction'])
            poses[c] = (left if c in groups['left'] else right)@local
        return poses

    trees = {c:cKDTree(p) for c,p in original.items()}
    bounds = {c:(p.min(0),p.max(0)) for c,p in original.items()}
    minima = {};worst = {};previous = {};samples = [];radius = .5
    for index,fraction in enumerate(np.linspace(0,1,2001)):
        matrices = trajectory(fraction);samples.append([matrices[c] for c in units])
        inverses = {c:np.linalg.inv(m) for c,m in matrices.items()}
        for a,b in itertools.combinations(units,2):
            pair = a+'/'+b;relative_pose = inverses[b]@matrices[a]
            key = np.round(relative_pose,12).tobytes()
            if previous.get(pair)==key:continue
            previous[pair] = key
            points = transform(original[a],relative_pose);low,high = bounds[b]
            selected = ((points>=low-radius)&(points<=high+radius)).all(1)
            value = min(radius,float(trees[b].query(points[selected],distance_upper_bound=radius)[0].min())) if selected.any() else radius
            if value<minima.get(pair,float('inf')):minima[pair] = value;worst[pair] = float(fraction)
        if index%500==0:print('BOUNDARY_SCREEN',index,'of',2000,flush=True)
    record['screen'] = {'samples':2001,'threshold_nm':.20,'search_radius_nm':radius,
        'pair_minima_nm_capped_at_radius':minima,'worst_fractions':worst,
        'scope':'All complete source-derived polymer units, with rigid fits, binding, fragment separation and terminal release evaluated together. No force-field or kinetics claim.'}
    record['passed'] = min(minima.values())>=.20
    output.write_text(json.dumps(record,indent=2)+'\n')
    if not record['passed']:raise ValueError('Cofilactin boundary motion overlaps; inspect cofilin-path-repairs.json')
    np.savez_compressed(folder/'motion.npz',matrices=np.array(samples),units=np.array(list(units)))
    record['motion_sha256'] = hashlib.sha256((folder/'motion.npz').read_bytes()).hexdigest()
    output.write_text(json.dumps(record,indent=2)+'\n')
    print('COFILIN_BOUNDARY_READY',flush=True)


if __name__=='__main__':main()
