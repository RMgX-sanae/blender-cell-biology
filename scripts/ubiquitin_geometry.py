"""Source-specific RNF4/Ubc13/Ube2V2 geometry, with explicit inferred torsions.

This is a display construction for 5AIT, not docking, MD, or a reaction solver.
Internal calculations use Angstrom; exported scene coordinates use nm.
"""
import hashlib
import itertools
from pathlib import Path

import gemmi
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from prepare_molecular_asset import atom_table, xyz
from prepare_phosphotransfer import peptide_bonds


SOURCE_HASHES = {
    '5AIT': 'b932d163add2f529115311f4b65b9cd5b094c0736d5841e08bdfc621d4cb61a9',
    '1J7D': '6c7a6877e39fb5a335fa5d7c7f099d6384500b0d71b906032a384ba031ab7e66',
    '1UBQ': '056f98710cb2b36f633c45e41902a02eb446e82871da21ff2dd44f74a56ca0f6',
    '2JF5': 'a54c5de3c6fcd87185c4b03ebca6c355467f2bf206c5bbcb9628b6c512e7e807',
}
PROTEINS = ('ring', 'e2', 'donor', 'uev', 'acceptor')


def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def index(rows, residue, name):
    return next(i for i, r in enumerate(rows) if r[3] == residue and r[2] == name)


def rotate(points, first, second, selection, angle):
    axis = points[second]-points[first]; axis /= np.linalg.norm(axis)
    points[selection] = (points[selection]-points[first])@Rotation.from_rotvec(axis*angle).as_matrix().T+points[first]


def load_sources(directory):
    directory = Path(directory)
    for name, expected in SOURCE_HASHES.items():
        if digest(directory/(name+'.cif')) != expected:
            raise ValueError('Source version changed; requalify before using this recipe: '+name)
    source = atom_table(directory/'5AIT.cif'); native = atom_table(directory/'1J7D.cif')
    e2 = [list(r) for r in source['B'] if r[3] not in (87, 92)]; changes = []
    for residue, name, added, removed in [(87, 'CYS', ['SG'], ['CG', 'CD', 'CE', 'NZ']), (92, 'LYS', ['CG', 'CD', 'CE', 'NZ'], [])]:
        a = {r[2]: r for r in native['B'] if r[3] == residue}
        b = {r[2]: r for r in source['B'] if r[3] == residue}; anchors = ['N', 'CA', 'C', 'CB']
        aa = np.array([a[n][0] for n in anchors]); bb = np.array([b[n][0] for n in anchors])
        u, _, vt = np.linalg.svd((aa-aa.mean(0)).T@(bb-bb.mean(0)))
        rotation = vt.T@np.diag([1, 1, np.linalg.det(vt.T@u.T)])@u.T
        shift = bb.mean(0)-rotation@aa.mean(0)
        e2.extend([[r[0], r[1], r[2], residue, name] for r in source['B'] if r[3] == residue and r[2] not in removed])
        e2.extend([[(rotation@a[n][0]+shift).tolist(), a[n][1], n, residue, name] for n in added])
        changes.append({'residue': residue, 'native_name': name, 'added': added, 'removed': removed,
                        'anchor_rmsd_A': float(np.sqrt(np.mean(np.sum((aa@rotation.T+shift-bb)**2, axis=1))))})
    e2.sort(key=lambda r: (r[3], r[2]))
    structure = gemmi.read_structure(str(directory/'5AIT.cif')); op = gemmi.Op('-x+y,-x,z+1/3')
    acceptor = []
    for row in source['F']:
        fractional = np.array(op.apply_to_xyz(list(structure.cell.fractionalize(gemmi.Position(*row[0])))))+[1, 1, 0]
        acceptor.append([list(structure.cell.orthogonalize(gemmi.Fractional(*fractional))), *row[1:]])
    # Free acceptor, not another E2-linked crystallographic copy. 1UBQ supports
    # terminal atom identity; this planar OXT placement is geometrically inferred.
    terminal = {r[2]: np.array(r[0]) for r in acceptor if r[3] == 76}
    ca = terminal['CA']-terminal['C']; ca /= np.linalg.norm(ca)
    oxygen = terminal['O']-terminal['C']; oxygen /= np.linalg.norm(oxygen)
    direction = -ca-oxygen; direction /= np.linalg.norm(direction)
    acceptor.append([(terminal['C']+1.25*direction).tolist(), 'O', 'OXT', 76, 'GLY'])
    assert any(r[3] == 76 and r[2] == 'OXT' for r in atom_table(directory/'1UBQ.cif')['A'])
    product = atom_table(directory/'2JF5.cif')
    product_length = np.linalg.norm(xyz(product['B'])[index(product['B'], 76, 'C')]-xyz(product['A'])[index(product['A'], 63, 'NZ')])
    assert 1.2 < product_length < 1.5
    zinc = [[list(atom.pos), atom.element.name, atom.name, residue.seqid.num, residue.name]
            for chain in structure[0] if chain.name == 'A' for residue in chain if residue.name == 'ZN' for atom in residue]
    assert len(zinc) == 4
    rows = {'ring': source['A'], 'e2': e2, 'donor': source['C'], 'uev': source['D'], 'acceptor': acceptor, 'zinc': zinc}
    observed = {k: sorted({r[3] for r in rows[k]}) for k in PROTEINS}
    return rows, {'native_restoration': changes, 'free_acceptor_OXT': 'inferred planar terminal oxygen, not a measured pose',
                  'product_reference_CN_A': float(product_length), 'acceptor_symmetry': '-x+y,-x,z+1/3; lattice [1,1,0]',
                  'observed_author_residue_numbers': observed,
                  'unmodelled_termini': 'Ub Met1, Ubc13 outside 4-150, and UEV outside 6-145 are not reconstructed. Whole-source views retain all selected resolved atoms, not complete native sequences.'}


class UbiquitinGeometry:
    def __init__(self, directory):
        self.rows, self.source_record = load_sources(directory)
        self.source = {k: xyz(v) for k, v in self.rows.items()}; r = self.rows
        self.ids = {'sg': index(r['e2'], 87, 'SG'), 'cb': index(r['e2'], 87, 'CB'),
                    'c': index(r['donor'], 76, 'C'), 'o': index(r['donor'], 76, 'O'), 'ca': index(r['donor'], 76, 'CA'),
                    'nz': index(r['acceptor'], 63, 'NZ'), 'ce': index(r['acceptor'], 63, 'CE')}
        self.donor_axes = {}
        for residue in (74, 75, 76):
            rows = r['donor']
            self.donor_axes[residue] = [
                (index(rows, residue, 'N'), index(rows, residue, 'CA'), np.array([i for i, row in enumerate(rows) if row[3] > residue or (row[3] == residue and row[2] != 'N')])),
                (index(rows, residue, 'CA'), index(rows, residue, 'C'), np.array([i for i, row in enumerate(rows) if row[3] > residue or (row[3] == residue and row[2] in ('O', 'OXT'))]))]
        self.lys_axes = [(index(r['acceptor'], 63, a), index(r['acceptor'], 63, b), np.array([index(r['acceptor'], 63, n) for n in mobile]))
                         for a, b, mobile in [('CA', 'CB', ['CG', 'CD', 'CE', 'NZ']), ('CB', 'CG', ['CD', 'CE', 'NZ']), ('CG', 'CD', ['CE', 'NZ']), ('CD', 'CE', ['NZ'])]]
        self.lys92 = [index(r['e2'], 92, n) for n in ('CD', 'CE', 'NZ')]
        self.pivot = self.source['acceptor'][cKDTree(self.source['uev']).query(self.source['acceptor'])[0] < 4.5].mean(0)
        self.uev_pivot = self.source['uev'][cKDTree(self.source['e2']).query(self.source['uev'])[0] < 4.5].mean(0)
        self.bonds, self.neighbors = {}, {}
        for key in PROTEINS:
            bonds = peptide_bonds(r[key], self.source[key]*.1); neighbors = [set() for _ in r[key]]
            for a, b in bonds: neighbors[a].add(int(b)); neighbors[b].add(int(a))
            self.bonds[key], self.neighbors[key] = bonds, neighbors
        self.exclusions = {}
        for a, b, i, j in [('e2', 'donor', self.ids['sg'], self.ids['c']), ('donor', 'acceptor', self.ids['c'], self.ids['nz'])]:
            self.exclusions[a, b] = {(i, j)}|{(k, j) for k in self.neighbors[a][i]}|{(i, k) for k in self.neighbors[b][j]}
        self.self_pairs = {}
        mobile = {'e2': [self.ids['sg'], *[index(r['e2'], 92, n) for n in ('CG', 'CD', 'CE', 'NZ')]], 'donor': np.concatenate([a[2] for a in self.donor_axes[74]]),
                  'acceptor': [index(r['acceptor'], 63, n) for n in ('CG', 'CD', 'CE', 'NZ')]}
        for key, indices in mobile.items():
            pairs = set()
            for i in set(indices):
                close = {i}|self.neighbors[key][i]
                for n in self.neighbors[key][i]: close |= self.neighbors[key][n]
                pairs.update(tuple(sorted([i, j])) for j in range(len(r[key])) if j not in close)
            self.self_pairs[key] = np.array(sorted(pairs))

    def pose(self, parameters):
        x = np.asarray(parameters, float)
        if x.shape != (24,) or not np.isfinite(x).all(): raise ValueError('Expected 24 finite source-specific torsion/pose parameters')
        p = {k: v.copy() for k, v in self.source.items()}
        rotate(p['e2'], index(self.rows['e2'], 87, 'CA'), self.ids['cb'], [self.ids['sg']], x[0])
        for residue, angles in [(75, x[1:3]), (76, x[3:5]), (74, x[15:17])]:
            for args, angle in zip(self.donor_axes[residue], angles): rotate(p['donor'], *args, angle)
        for args, angle in zip(self.lys_axes, x[5:9]): rotate(p['acceptor'], *args, angle)
        p['acceptor'] = (p['acceptor']-self.pivot)@Rotation.from_rotvec(x[9:12]).as_matrix().T+self.pivot+x[12:15]
        p['uev'] = (p['uev']-self.uev_pivot)@Rotation.from_rotvec(x[17:20]).as_matrix().T+self.uev_pivot+x[20:23]
        rotate(p['e2'], self.lys92[0], self.lys92[1], [self.lys92[2]], x[23])
        return p

    def check(self, points, state):
        if state not in ('reactant', 'transition', 'product'): raise ValueError('Unknown chemical state')
        minima = {}
        for a, b in itertools.combinations(PROTEINS, 2):
            distance, near = cKDTree(points[a]).query(points[b], k=4)
            excluded = self.exclusions.get((a, b), set())
            if (a, b) == ('e2', 'donor') and state == 'product': excluded = set()
            if (a, b) == ('donor', 'acceptor') and state == 'reactant': excluded = set()
            for i, j in excluded: distance[j, near[j] == i] = np.inf
            minima[a+'/'+b] = float(distance.min())
        for key in ('e2', 'donor', 'uev', 'acceptor'):
            minima[key+'/zinc'] = float(cKDTree(points['zinc']).query(points[key])[0].min())
        internal = {k: float(np.linalg.norm(points[k][e[:, 0]]-points[k][e[:, 1]], axis=1).min()) for k, e in self.self_pairs.items()}
        c = points['donor'][self.ids['c']]
        o = points['donor'][self.ids['o']]; ca = points['donor'][self.ids['ca']]
        n = points['acceptor'][self.ids['nz']]; ce = points['acceptor'][self.ids['ce']]
        s = points['e2'][self.ids['sg']]; cb = points['e2'][self.ids['cb']]
        def angle(a, b, d):
            u, v = a-b, d-b
            return float(np.degrees(np.arccos(np.clip(np.dot(u, v)/(np.linalg.norm(u)*np.linalg.norm(v)), -1, 1))))
        def torsion(a, b, d, e):
            axis = d-b; axis /= np.linalg.norm(axis)
            u, v = a-b, e-d; u -= np.dot(u, axis)*axis; v -= np.dot(v, axis)*axis
            return float(np.degrees(np.arctan2(np.dot(np.cross(axis, u), v), np.dot(u, v))))
        normal = np.cross(o-c, ca-c); normal /= np.linalg.norm(normal)
        chemistry = {'O_C_S_degrees': angle(o, c, s), 'CB_S_C_degrees': angle(cb, s, c),
            'thioester_carbonyl_plane_error_degrees': float(np.degrees(np.arcsin(np.clip(abs(np.dot(normal, s-c))/np.linalg.norm(s-c), 0, 1)))),
            'O_C_N_degrees': angle(o, c, n), 'C_N_CE_degrees': angle(c, n, ce),
            'O_C_N_CE_degrees': torsion(o, c, n, ce), 'CA_C_N_CE_degrees': torsion(ca, c, n, ce)}
        return {'C_S_A': float(np.linalg.norm(c-points['e2'][self.ids['sg']])),
                'C_N_A': float(np.linalg.norm(c-points['acceptor'][self.ids['nz']])),
                'nonbonded_minima_A': minima, 'mobile_self_nonbonded_minima_A': internal, 'reaction_angles': chemistry}
