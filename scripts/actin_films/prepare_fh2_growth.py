"""Prepare source-guided FH2 elongation geometry and authored, screened motion.

The two 9B27 functional lobes exchange endpoint conformations. Their trajectories,
the G/F interpolation and their timing are teaching choices, not measured paths.
No modeled atoms, surface vertices or faces are removed to create loading space.
"""
import hashlib
import itertools
import json

import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation, Slerp

from prepare_contact_repairs import fit_placement, transform
from project_context import ROOT


def arrays(path):
    with np.load(path) as data:
        return {key:data[key] for key in data.files}


def smooth(value):
    value = np.clip(value, 0, 1)
    return value**3 * (10 + value * (-15 + 6 * value))


def rigid_fit(first, second):
    u, _, vt = np.linalg.svd((first-first.mean(0)).T @ (second-second.mean(0)))
    matrix = np.eye(4)
    matrix[:3, :3] = vt.T @ np.diag([1, 1, np.linalg.det(vt.T @ u.T)]) @ u.T
    matrix[:3, 3] = second.mean(0)-matrix[:3, :3] @ first.mean(0)
    return matrix


def prepare_conformers(assets, folder, seed):
    data = {c:arrays(assets/'9B27'/(c+'.npz')) for c in 'GH'}
    ca = {c:dict(zip(data[c]['ca_residues'],data[c]['ca'])) for c in 'GH'}
    records = {}
    for chain, other in [('G', 'H'), ('H', 'G')]:
        source = arrays(assets/'fh2_seed'/(chain+'.npz'))
        ownership = seed['lobes'][chain]['source_atom_ids']
        first_ca = [];second_ca = [];ca_parts = [];target_ids = []
        for name, index in ownership:
            target = other if name == chain else chain
            residue = int(data[name]['residue_numbers'][index]);atom = str(data[name]['atom_names'][index])
            candidates = np.flatnonzero((data[target]['residue_numbers'] == residue) & (data[target]['atom_names'] == atom))
            target_ids.append((target, int(candidates[0])) if len(candidates) else None)
            if atom == 'CA' and len(candidates):
                first_ca.append(ca[name][residue]);second_ca.append(ca[target][residue])
                ca_parts.append(name)
        first_ca = np.array(first_ca);second_ca = np.array(second_ca)
        fit = rigid_fit(first_ca, second_ca);inverse = np.linalg.inv(fit)
        delta = transform(second_ca, inverse)-first_ca

        def warp(points, part=None):
            selected = np.ones(len(first_ca),dtype=bool) if part is None else np.array(ca_parts)==part
            distance, index = cKDTree(first_ca[selected]).query(points, k=12)
            weights = 1/(distance**2+.02);weights /= weights.sum(1)[:, None]
            return points+np.sum(delta[selected][index]*weights[:, :, None], axis=1)

        alternate_atoms = source['atoms'].copy()
        for part in [chain,other]:
            selected = np.array([name==part for name,_ in ownership])
            alternate_atoms[selected] = warp(source['atoms'][selected],part)
        for i, target in enumerate(target_ids):
            if target is not None:
                name, index = target
                alternate_atoms[i] = transform(data[name]['atoms'][index:index+1], inverse)[0]
        alternate_vertices = warp(source['vertices'])
        np.savez_compressed(folder/(chain+'.npz'), f_vertices=source['vertices'],
            g_vertices=alternate_vertices, faces=source['faces'], f_atoms=source['atoms'],
            g_atoms=alternate_atoms, material_indices=source['material_indices'])
        records[chain] = {**seed['lobes'][chain], 'alternate_fit':fit.tolist(),
            'alternate_source_atom_ids':target_ids, 'unmatched_atoms':target_ids.count(None),
            'matched_CA':len(first_ca),
            'CA_RMSD_nm':float(np.sqrt(np.mean(np.sum(delta**2, axis=1)))),
            'maximum_surface_displacement_nm':float(np.linalg.norm(alternate_vertices-source['vertices'], axis=1).max()),
            'atom_displacement_percentiles_nm':np.percentile(np.linalg.norm(alternate_atoms-source['atoms'], axis=1), [0,50,95,100]).tolist()}
    actin = arrays(assets/'actin_morph.npz')
    actin['g_atoms'] = actin['g_atoms'].astype(float)
    source = arrays(assets/'2BTF/A.npz');target = arrays(assets/'9B27/F.npz')
    aligned = source['atoms'] @ actin['g_to_f_rotation'].T+actin['g_to_f_shift']
    by_name = {(int(n), str(k)):(i, point) for i,(n,k,point) in enumerate(zip(source['residue_numbers'], source['atom_names'], aligned))}
    matched = []
    for i, (residue, atom) in enumerate(zip(target['residue_numbers'], target['atom_names'])):
        value = by_name.get((int(residue), str(atom)))
        if value is not None:
            index, point = value;actin['g_atoms'][i] = point;matched.append([i,index])
    np.savez_compressed(folder/'actin.npz', **actin)
    records['actin'] = {'source_F':'9B27 F', 'source_G':'2BTF A', 'matched_atom_indices_F_G':matched,
        'unmatched_F_atoms_retained':len(target['atoms'])-len(matched),
        'limits':'Cross-isoform correspondence, not a same-molecule experimental trajectory. G surface and unmatched atoms use the CA-guided approximation; matched G atoms use the aligned source coordinates.'}
    return records


def fit_actin_ends(points, obstacles):
    centers = [p.mean(0) for p in points];trees = [cKDTree(p) for p in obstacles]

    def positions(values):
        return [(p-c) @ Rotation.from_rotvec(q[3:]).as_matrix().T+c+q[:3]
                for p,c,q in zip(points, centers, values.reshape(-1,6))]

    def objective(values):
        posed = positions(values);q = values.reshape(-1,6)
        loss = .03*np.sum(q[:,:3]**2)+np.sum(q[:,3:]**2)
        for i, p in enumerate(posed):
            for tree in trees+[cKDTree(v) for v in posed[:i]]:
                distance, _ = tree.query(p)
                loss += 10000*np.sum(np.maximum(.225-distance, 0)**2)
        return loss

    fit = minimize(objective, np.zeros(6*len(points)), method='L-BFGS-B',
        bounds=([(-.4,.4)]*3+[(-.09,.09)]*3)*len(points),
        options={'maxiter':350, 'ftol':1e-13, 'gtol':1e-7, 'maxls':40})
    corrections = []
    for center,q in zip(centers,fit.x.reshape(-1,6)):
        matrix = np.eye(4);matrix[:3,:3] = Rotation.from_rotvec(q[3:]).as_matrix()
        matrix[:3,3] = center+q[:3]-matrix[:3,:3] @ center;corrections.append(matrix)
    posed = positions(fit.x)
    minimum = min(float(cKDTree(b).query(p)[0].min()) for i,p in enumerate(posed) for b in obstacles+posed[:i])
    displacement = max(float(np.linalg.norm(p-q,axis=1).max()) for p,q in zip(points,posed))
    if minimum < .21 or displacement > .5:raise ValueError('FH2 growth endpoint fit requires source review')
    return corrections, {'minimum_nm':minimum, 'maximum_atom_displacement_nm':displacement,
        'optimizer_success':bool(fit.success), 'optimizer_message':str(fit.message)}


def prepare_endpoints(assets, folder, conformers):
    data = {c:np.load(assets/'9B27'/(c+'.npz'))['atoms'] for c in 'ABCDEFGH'}
    step = np.array(json.loads((assets/'9B27/manifest.json').read_text())['step'])
    targets = [step, step @ step]
    corrections, report = fit_actin_ends([transform(data['F'],m) for m in targets],
        [data[c] for c in 'ABCDEF']+[transform(data[c],targets[1]) for c in 'GH'])
    targets = [a @ b for a,b in zip(corrections,targets)]
    result = {'actin_targets':[m.tolist() for m in targets], 'actin_end_fit':report, 'head_targets':{}}
    lobes = {c:arrays(folder/(c+'.npz')) for c in 'GH'}
    core = np.concatenate([data[c] for c in 'ABCDEF'])
    for cycle in [1,2]:
        poses = {c:np.linalg.matrix_power(step,cycle) @ (np.array(conformers[c]['alternate_fit']) if cycle==1 else np.eye(4)) for c in 'GH'}
        parts = {c:transform(lobes[c]['g_atoms' if cycle==1 else 'f_atoms'],poses[c]) for c in 'GH'}
        actins = [transform(data['F'],m) for m in targets[:cycle]]
        for c in 'GH':
            correction, record = fit_placement(parts[c],np.concatenate([core]+actins+[parts['H' if c=='G' else 'G']]))
            record['matrix'] = (correction @ poses[c]).tolist();result['head_targets'][c+str(cycle)] = record
    return result


def refine_motion(fractions, mobile, pivots, fixed, label, *, target=.23, iterations=750, maximum_shift=.6):
    centered = mobile-pivots[:,None];trees = [cKDTree(p) for p in fixed]
    knots = np.r_[np.zeros(4), np.linspace(0,1,11)[1:-1], np.ones(4)]
    basis = BSpline(knots,np.eye(len(knots)-4),3)((fractions-.10)/.82)[:,2:-2]
    count = basis.shape[1];iteration = 0

    def progress(_):
        nonlocal iteration
        iteration += 1
        if iteration%100==0:print('FH2_GROWTH_OPTIMIZING',label,iteration,flush=True)

    def objective(flat):
        controls = flat.reshape(count,6);poses = basis @ controls;gradient = np.zeros_like(poses)
        loss = .05*np.sum(controls[:,:3]**2)+2*np.sum(controls[:,3:]**2)
        for i,q in enumerate(poses):
            moved = centered[i] @ Rotation.from_rotvec(q[3:]).as_matrix().T+pivots[i]+q[:3]
            distance, neighbor = trees[i].query(moved, distance_upper_bound=target)
            selected = distance < target;penetration = target-distance[selected]
            loss += 10000*(penetration @ penetration)
            force = -20000*penetration[:,None]*(moved[selected]-fixed[i][neighbor[selected]])/np.maximum(distance[selected,None],1e-9)
            gradient[i,:3] = force.sum(0)
            if selected.any():
                for axis in range(3):
                    h = np.zeros(3);h[axis] = 1e-6
                    derivative = (Rotation.from_rotvec(q[3:]+h).as_matrix()-Rotation.from_rotvec(q[3:]-h).as_matrix())/2e-6
                    gradient[i,3+axis] = np.sum(force*(centered[i,selected] @ derivative.T))
        result = basis.T @ gradient;result[:,:3] += .1*controls[:,:3];result[:,3:] += 4*controls[:,3:]
        return loss, result.ravel()

    fit = minimize(objective, np.zeros(count*6), jac=True, method='L-BFGS-B',
        bounds=([(-maximum_shift,maximum_shift)]*3+[(-.15,.15)]*3)*count,
        callback=progress,
        options={'maxiter':iterations, 'ftol':1e-12, 'gtol':1e-7, 'maxls':40})
    controls = np.vstack([np.zeros((2,6)),fit.x.reshape(count,6),np.zeros((2,6))])
    record = {'first_fraction':.10, 'last_fraction':.92, 'knots':knots.tolist(),
        'controls_nm_radians':controls.tolist(), 'optimizer_success':bool(fit.success),
        'optimizer_message':str(fit.message), 'iterations':int(fit.nit),
        'scope':'Rigid correction of complete modeled geometry; not molecular dynamics'}
    print('FH2_GROWTH_PATH',label,fit.message,flush=True)
    return record


class Growth:
    def __init__(self, assets, folder, endpoints):
        self.core = {c:np.load(assets/'9B27'/(c+'.npz'))['atoms'].astype(float) for c in 'ABCDEF'}
        self.actin = arrays(folder/'actin.npz');self.lobes = {c:arrays(folder/(c+'.npz')) for c in 'GH'}
        self.profilin = np.load(assets/'profilin_on_morph.npz')['atoms']
        self.escape = self.profilin.mean(0)-self.actin['g_atoms'].mean(0)
        self.escape *= 8/np.linalg.norm(self.escape)
        self.step = np.array(json.loads((assets/'9B27/manifest.json').read_text())['step'])
        self.offset = np.array(json.loads((assets/'approach-plans.json').read_text())['formin_add_1']['offset'])
        self.offset *= 10/np.linalg.norm(self.offset)
        self.endpoints = endpoints;self.repairs = {};self.splines = {}
        # Authored broad approach arc in the prepared 9B27 coordinate frame.
        # These are choreography controls, not experimental coordinates. The
        # source-dependent rigid corrections and dense screen run below.
        self.approach_controls = np.array([[0,0,0],[0,0,0],
            [.7544785047587277,-1.899359217060999,-3.5],
            [-.8198989524625281,-.9373020014506968,-3.5],
            [-.22448495224442636,-.25212069718343255,-.1632472873572147],
            [.057870982868324566,.046767614847914446,.040549335474331254],
            [-.022392597769487255,-.027430618019937264,-.012678781724506816],
            [0,0,0],[0,0,0]])
        self.approach_knots = np.r_[np.zeros(4),np.linspace(0,1,7)[1:-1],np.ones(4)]
        self.approach = BSpline(self.approach_knots,self.approach_controls,3)

    def add_repair(self, name, record):
        self.repairs[name] = record
        self.splines[name] = BSpline(record['knots'],record['controls_nm_radians'],3)

    def correction(self, name, u):
        return self.splines[name]((u-.10)/.82) if name in self.splines and .10<u<.92 else np.zeros(6)

    @staticmethod
    def corrected(matrix, local_points, q, turn=None):
        turn = np.eye(3) if turn is None else turn
        pivot = transform(local_points,matrix).mean(0)
        correction = np.eye(4);correction[:3,:3] = Rotation.from_rotvec(turn @ q[3:]).as_matrix()
        correction[:3,3] = pivot-correction[:3,:3] @ pivot+turn @ q[:3]
        return correction @ matrix

    def incoming(self, u, cycle):
        turn = np.eye(3) if cycle==0 else self.step[:3,:3]
        travel = 1-.82*smooth((u-.08)/.22)-.18*smooth((u-.42)/.23)
        matrix = np.array(self.endpoints['actin_targets'][cycle])
        matrix[:3,3] += turn @ (self.offset*travel)
        if .125<u<.9:matrix[:3,3] += turn @ self.approach((u-.125)/.775)
        g = 1-smooth((u-.49)/.16)
        local = self.actin['f_atoms']*(1-g)+self.actin['g_atoms']*g
        matrix = self.corrected(matrix,local,self.correction('incoming_first',u),turn)
        if cycle==1:matrix = self.corrected(matrix,local,self.correction('incoming_second',u))
        return matrix, g, self.escape*smooth((u-.48)/.17)

    def head(self, chain, u, cycle):
        opening = chain==('G' if cycle==0 else 'H')
        q = smooth((u-.41)/.34) if opening else smooth((u-.18)/.55)
        start = np.eye(4) if cycle==0 else np.array(self.endpoints['head_targets'][chain+'1']['matrix'])
        end = np.array(self.endpoints['head_targets'][chain+str(cycle+1)]['matrix'])
        matrix = np.eye(4)
        matrix[:3,:3] = Slerp([0,1],Rotation.from_matrix([start[:3,:3],end[:3,:3]]))(q).as_matrix()
        matrix[:3,3] = (1-q)*start[:3,3]+q*end[:3,3]
        if opening:
            turn = np.eye(3) if cycle==0 else self.step[:3,:3]
            matrix[:3,3] += turn @ np.array([0,-1.,-3.])*smooth((u-.25)/.15)*(1-smooth((u-.65)/.18))
        g = q if cycle==0 else 1-q
        if cycle==1 and chain=='G':
            local = self.lobes[chain]['f_atoms']*(1-g)+self.lobes[chain]['g_atoms']*g
            matrix = self.corrected(matrix,local,self.correction('head_second',u))
        return matrix, g

    def poses(self, fraction):
        cycle = 0 if fraction<=.5 else 1
        u = fraction*2 if cycle==0 else (fraction-.5)*2
        incoming = [self.incoming(u if i==cycle else (1 if i<cycle else 0),i) for i in range(2)]
        heads = [self.head(c,u,cycle) for c in 'GH']
        return {'matrices':np.array([v[0] for v in incoming+heads]),
            'g_fractions':np.array([v[1] for v in incoming+heads]),
            'profilin_shifts':np.array([v[2] for v in incoming])}

    def state(self, fraction):
        pose = self.poses(fraction);out = dict(self.core)
        for i in range(2):
            matrix = pose['matrices'][i];g = pose['g_fractions'][i]
            out['actin'+str(i+1)] = transform(self.actin['f_atoms']*(1-g)+self.actin['g_atoms']*g,matrix)
            out['profilin'+str(i+1)] = transform(self.profilin+pose['profilin_shifts'][i],matrix)
        for i,chain in enumerate('GH',2):
            g = pose['g_fractions'][i]
            out['lobe'+chain] = transform(self.lobes[chain]['f_atoms']*(1-g)+self.lobes[chain]['g_atoms']*g,pose['matrices'][i])
        return out


def linker_records(assets, seed, conformers, model):
    owner = {tuple(source):(chain,i) for chain in 'GH'
             for i,source in enumerate(conformers[chain]['source_atom_ids'])}
    records = []
    for original in seed['linkers']:
        record = dict(original);data = np.load(assets/'9B27'/(record['chain']+'.npz'))
        for label, residue in zip(['start','end'],record['residue_interval']):
            index = int(np.flatnonzero((data['residue_numbers']==residue) & (data['atom_names']=='CA'))[0])
            chain, i = owner[(record['chain'],index)]
            record[label+'_anchor'] = {'lobe':chain, 'atom_index':i,
                'position':model.lobes[chain]['f_atoms'][i].tolist(),
                'alternate_position':model.lobes[chain]['g_atoms'][i].tolist()}
        records.append(record)
    return records


def screen_sequence(model, linkers, samples=2001):
    minima = {};lengths = {r['chain']:0. for r in linkers};poses = []
    maximum_nearest_lobe = 0.;detached = []
    for index,fraction in enumerate(np.linspace(0,1,samples)):
        state = model.state(fraction);poses.append(model.poses(fraction))
        bounds = {name:(points.min(0),points.max(0)) for name,points in state.items()};distances = {}
        trees = {name:cKDTree(points) for name,points in state.items()}
        for a,b in itertools.combinations(state,2):
            if index and a in 'ABCDEF' and b in 'ABCDEF':continue
            pair = a+'/'+b;radius = .45
            lo,hi = bounds[b];selected = ((state[a]>=lo-radius) & (state[a]<=hi+radius)).all(1)
            value = float(trees[b].query(state[a][selected],distance_upper_bound=radius)[0].min()) if selected.any() else radius
            distances[frozenset([a,b])] = min(radius,value)
            minima[pair] = min(minima.get(pair,radius),value)
        connected = list('ABCDEF')
        for name in ['actin1','actin2']:
            if min(distances[frozenset([name,c])] for c in connected)<.45:connected.append(name)
        nearest_lobe = min(distances[frozenset(['lobe'+c,actin])] for c in 'GH' for actin in connected)
        maximum_nearest_lobe = max(maximum_nearest_lobe,nearest_lobe)
        if nearest_lobe>=.45:detached.append(float(fraction))
        for record in linkers:
            a,b = [state['lobe'+record[label+'_anchor']['lobe']][record[label+'_anchor']['atom_index']] for label in ['start','end']]
            radial = np.array([0,(a[1]+b[1])/2,(a[2]+b[2])/2]);radial /= np.linalg.norm(radial)
            t = np.linspace(0,1,21)[:,None]
            points = (1-t)*a+t*b+.4*np.sin(np.pi*t)*radial
            length = float(np.linalg.norm(np.diff(points,axis=0),axis=1).sum())
            lengths[record['chain']] = max(lengths[record['chain']],length)
        if index%500==0:print('FH2_GROWTH_SCREEN',index,'of',samples-1,flush=True)
    report = {'samples':samples, 'threshold_nm':.20, 'distance_search_radius_nm':.45,
        'pair_minima_nm_capped_at_search_radius':minima,
        'maximum_linker_lengths_nm':lengths,
        'association_screen':{'geometric_contact_radius_nm':.45,
            'maximum_nearest_lobe_distance_nm':maximum_nearest_lobe,
            'both_lobes_detached_at_fractions':detached},
        'limits':'Heavy atom centers, source-guided interpolated coordinates and schematic linker contour bounds; not a steric energy or dynamics calculation.'}
    report['passed'] = not detached and min(minima.values())>=.20 and all(lengths[r['chain']]<=r['contour_upper_bound_nm'] for r in linkers)
    return report, {key:np.array([pose[key] for pose in poses]) for key in poses[0]}


def main():
    assets = ROOT/'assets';folder = assets/'fh2_growth';folder.mkdir(exist_ok=True)
    seed = json.loads((assets/'fh2-seed-repairs.json').read_text())
    record = {'schema_version':1, 'passed':False,
        'sources':{code:{'url':'https://www.rcsb.org/structure/'+code,
            'sha256':hashlib.sha256((ROOT/'sources'/(code+'.cif')).read_bytes()).hexdigest()} for code in ['9B27','2BTF']},
        'references':['https://doi.org/10.1038/s41586-024-07637-0','https://doi.org/10.1126/science.adn9560'],
        'limits':'9B27 is a mature mDia1 end and 2BTF is a different actin isoform. Loading order, lobe stepping, G/F interpolation, disordered linkers and durations are authored teaching choices, not observed intermediates or molecular dynamics. Different formins need separate source qualification.'}
    manifest = assets/'fh2-growth.json';manifest.write_text(json.dumps(record,indent=2)+'\n')
    conformers = prepare_conformers(assets,folder,seed)
    endpoints = prepare_endpoints(assets,folder,conformers)
    print('FH2_GROWTH_ENDPOINTS_READY',flush=True)
    model = Growth(assets,folder,endpoints);fractions = np.linspace(.10,.92,165)
    first = [model.state(u*.5) for u in fractions]
    moving = np.array([np.concatenate([s['actin1'],s['profilin1']]) for s in first])
    pivots = np.array([s['actin1'].mean(0) for s in first])
    fixed = [np.concatenate([s[c] for c in list('ABCDEF')+['lobeG','lobeH']]) for s in first]
    model.add_repair('incoming_first',refine_motion(fractions,moving,pivots,fixed,'incoming_first',target=.225,iterations=650,maximum_shift=.5))
    second = [model.state(.5+u*.5) for u in fractions]
    moving = np.array([s['lobeG'] for s in second]);pivots = moving.mean(1)
    fixed = [np.concatenate([s[c] for c in list('ABCDEF')+['actin1','profilin1','lobeH']]) for s in second]
    model.add_repair('head_second',refine_motion(fractions,moving,pivots,fixed,'head_second'))
    second = [model.state(.5+u*.5) for u in fractions]
    moving = np.array([np.concatenate([s['actin2'],s['profilin2']]) for s in second])
    pivots = np.array([s['actin2'].mean(0) for s in second])
    fixed = [np.concatenate([s[c] for c in list('ABCDEF')+['actin1','profilin1','lobeG','lobeH']]) for s in second]
    model.add_repair('incoming_second',refine_motion(fractions,moving,pivots,fixed,'incoming_second'))
    linkers = linker_records(assets,seed,conformers,model)
    report, poses = screen_sequence(model,linkers)
    record.update({'conformers':conformers, 'endpoints':endpoints, 'repairs':model.repairs,
        'authored_approach':{'first_fraction':.125,'last_fraction':.9,
            'controls_nm':model.approach_controls.tolist(),'knots':model.approach_knots.tolist(),
            'initial_offset_nm':model.offset.tolist()},
        'linkers':linkers, 'screen':report, 'pose_names':['actin1','actin2','lobeG','lobeH'], 'passed':report['passed']})
    manifest.write_text(json.dumps(record,indent=2)+'\n')
    if not report['passed']:raise ValueError('Coupled FH2 growth failed clearance or linker checks; inspect fh2-growth.json')
    np.savez_compressed(folder/'motion.npz',**poses)
    record['motion_sha256'] = hashlib.sha256((folder/'motion.npz').read_bytes()).hexdigest()
    manifest.write_text(json.dumps(record,indent=2)+'\n')
    print('FH2_GROWTH_READY',min(report['pair_minima_nm_capped_at_search_radius'].values()),flush=True)


if __name__=='__main__':main()
