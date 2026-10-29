"""Inferred source-specific K63 product release with intact E2/UEV partners.

Rigid whole-protein departures and declared side-chain torsions preserve the
formed ubiquitin bond. This is not an experimental trajectory or MD model.
"""
import numpy as np
from scipy.spatial import cKDTree

from interaction_motion import matrices
from ubiquitin_geometry import UbiquitinGeometry, index, rotate


def smooth(t):
    x = np.clip(t, 0, 1)
    return x*x*x*(10+x*(-15+6*x))


class UbiquitinReleaseGeometry(UbiquitinGeometry):
    def __init__(self, directory):
        super().__init__(directory)
        self.release_axes = []
        for key, residue, definitions in [
            ('e2', 119, [('CA', 'CB', ['CG', 'OD1', 'OD2']), ('CB', 'CG', ['OD1', 'OD2'])]),
            ('uev', 32, [('CA', 'CB', ['OG'])]),
            ('acceptor', 62, [('CA', 'CB', ['CG', 'CD', 'OE1', 'NE2']), ('CB', 'CG', ['CD', 'OE1', 'NE2']), ('CG', 'CD', ['OE1', 'NE2'])]),
        ]:
            for first, second, selection in definitions:
                r = self.rows[key]
                self.release_axes.append((key, index(r, residue, first), index(r, residue, second), [index(r, residue, a) for a in selection]))
        # Add every newly mobile side-chain atom to the existing self-contact
        # screen. Only chemical graph 1-2 and 1-3 neighbors are excluded.
        for key in ('e2', 'uev', 'acceptor'):
            mobile = {i for k, _, _, selected in self.release_axes if k == key for i in selected}
            pairs = set(map(tuple, self.self_pairs.get(key, [])))
            for i in mobile:
                close = {i}|self.neighbors[key][i]
                for j in self.neighbors[key][i]: close |= self.neighbors[key][j]
                pairs.update(tuple(sorted((i, j))) for j in range(len(self.rows[key])) if j not in close)
            self.self_pairs[key] = np.array(sorted(pairs))

    def release_pose(self, product, spec, time):
        p = {k: v.copy() for k, v in product.items()}
        begin, end = spec['partner_interval']; q = np.clip((time-begin)/(end-begin), 0, 1)
        partner = spec['partner']; angles = np.asarray(partner['sidechain_torsions_radians'], float)
        if angles.shape != (6,) or not np.isfinite(angles).all(): raise ValueError('Expected six finite interface torsions')
        fraction = smooth(q/.15)*(1-smooth((q-.60)/.35))
        for (key, first, second, selected), angle in zip(self.release_axes, angles): rotate(p[key], first, second, selected, angle*fraction)
        travel = np.clip((q-.15)/.75, 0, 1)
        mobile = np.concatenate([product[k] for k in ('e2', 'uev')])*.1
        transform = matrices(mobile, partner['direction'], partner['distance_nm'], np.array([.5+.5*travel]), 'loop', partner['turn_radians'], partner['bend_fraction'])[0]
        side = smooth((q-partner['separation_start'])/(1-partner['separation_start']))
        transform[:3, 3] += np.asarray(partner['separation_offset_nm'])*side
        for key in ('e2', 'uev'): p[key] = p[key]@transform[:3, :3].T+transform[:3, 3]*10
        begin, end = spec['product_interval']; q = np.clip((time-begin)/(end-begin), 0, 1)
        departure = spec['product']; mobile = np.concatenate([product[k] for k in ('donor', 'acceptor')])*.1
        travel = np.clip(q/departure['primary_end'], 0, 1)
        transform = matrices(mobile, departure['direction'], departure['distance_nm'], np.array([.5+.5*travel]), 'loop', departure['turn_radians'], departure['bend_fraction'])[0]
        side = smooth((q-departure['separation_start'])/(1-departure['separation_start']))
        transform[:3, 3] += np.asarray(departure['separation_offset_nm'])*side
        for key in ('donor', 'acceptor'): p[key] = p[key]@transform[:3, :3].T+transform[:3, 3]*10
        return p

    def release_checks(self, arrays, states):
        selected = np.flatnonzero(np.asarray(states) == 'product')
        ca = {k: np.array([i for i, row in enumerate(self.rows[k]) if row[2] == 'CA']) for k in ('e2', 'uev')}
        first = selected[0]
        neighbors = cKDTree(arrays['e2'][first, ca['e2']]).query_ball_point(arrays['uev'][first, ca['uev']], 12.)
        pairs = np.array([(ca['e2'][i], ca['uev'][j]) for j, group in enumerate(neighbors) for i in group])
        if len(pairs) < 3: raise ValueError('Missing source E2/UEV backbone interface')
        distances = np.linalg.norm(arrays['e2'][selected][:, pairs[:, 0]]-arrays['uev'][selected][:, pairs[:, 1]], axis=2)
        error = float(np.max(abs(distances-distances[0])))
        product = np.concatenate([arrays[k][-1] for k in ('donor', 'acceptor')])
        others = np.concatenate([arrays[k][-1] for k in ('ring', 'zinc', 'e2', 'uev')])
        gap = float(cKDTree(others).query(product)[0].min())
        return {'passed': error < 1e-6 and gap > 8., 'retained_partner_CA_pairs': len(pairs),
                'maximum_partner_CA_distance_change_A': error, 'final_product_to_other_heavy_atom_gap_A': gap,
                'limit': 'The final gap is a display-separation criterion, not an affinity or biological unbinding threshold.'}
