"""Reviewed angles and handedness across physical atom aliases and coordinate frames."""
import copy
import unittest

import numpy as np

from reaction_ledger import validate
from reaction_pose_geometry import check_pose_geometry, measure_nm


class PoseGeometry(unittest.TestCase):
    def setUp(self):
        names = ['A', 'B', 'D', 'X']
        new = ['p/' + n for n in names]
        self.ledger = dict(version=2, scope='technical', bond_scope='three center bonds',
                           reactant_atoms=[(n, 'C') for n in names], product_atoms=[(n, 'C') for n in new],
                           reactant_bonds=[['X', n, 'sing'] for n in names[:3]],
                           product_bonds=[['p/X', n, 'sing'] for n in new[:3]],
                           atom_map=list(zip(names, new)), form=[], change_order=[], **{'break': []})
        self.ids = {'left': names[:2], 'right': names[2:]}
        # Neighbors in z=0 and center at z=+0.5 Angstrom. The second pose is a proper rigid transform.
        xyz = np.array([[.1, 0, 0], [0, .1, 0], [0, 0, 0], [0, 0, .05]])
        turn = np.array([[0., -1, 0], [1, 0, 0], [0, 0, 1]])
        self.xyz = np.stack([xyz, xyz @ turn + [2., -1., 4.]])
        self.review = {'version': 2, 'poses': []}
        for i, (a, b, d, x) in enumerate((names, new)):
            self.review['poses'].append(dict(state=i, sample=i,
                angles=[dict(atoms=[a, x, b], angle_deg=[78., 79.], basis='analytic technical tripod')],
                signed_heights=[dict(atoms=[a, b, d, x], height_A=[.45, .55], basis='ordered technical reference')]))

    def check(self):
        self.assertTrue(validate(self.ledger)['passed'])
        motion = {'left': self.xyz[:, :2], 'right': self.xyz[:, 2:]}
        return check_pose_geometry(self.review, self.ledger, self.ids, motion, 1)

    def test_units_proper_rotation_and_product_aliases(self):
        first, last = self.check()
        self.assertAlmostEqual(first['angles'][0]['angle_deg'], np.degrees(np.arccos(.2)))
        self.assertEqual(last['signed_heights'][0]['height_A'], .5)
        self.assertEqual(first['signed_heights'][0]['height_A'], last['signed_heights'][0]['height_A'])
        self.assertEqual(last['signed_heights'][0]['physical_atoms'], ['A', 'B', 'D', 'X'])
        self.assertEqual(last['signed_heights'][0]['bond_orders'], ['sing'] * 3)

    def test_mirror_preserves_angles_but_reverses_signed_height(self):
        self.xyz[1, 3, 2] -= .1
        with self.assertRaisesRegex(ValueError, 'signed_heights.*outside'):
            self.check()
        del self.review['poses'][1]['signed_heights']
        self.assertAlmostEqual(self.check()[1]['angles'][0]['angle_deg'], np.degrees(np.arccos(.2)))

    def test_bending_a_single_arm_keeps_bond_length(self):
        # Rotate arm A around X until it points along arm B.
        self.xyz[0, 0] = self.xyz[0, 1]
        pose = self.review['poses'][0]
        pose['bonds'] = [dict(atoms=['A', 'X'], order='sing', length_A=[1.11, 1.12], basis='fixed length')]
        with self.assertRaisesRegex(ValueError, 'angles.*outside'):
            self.check()

    def test_neighbor_order_controls_sign_and_duplicate_detection(self):
        pose = self.review['poses'][0]
        row = pose['signed_heights'][0]
        row['atoms'] = ['B', 'D', 'A', 'X']
        self.assertEqual(self.check()[0]['signed_heights'][0]['height_A'], .5)
        row['atoms'] = ['A', 'D', 'B', 'X']
        with self.assertRaisesRegex(ValueError, 'outside'):
            self.check()
        row['height_A'] = [-.55, -.45]
        self.assertEqual(self.check()[0]['signed_heights'][0]['height_A'], -.5)
        pose['signed_heights'].append(dict(row, atoms=['A', 'B', 'D', 'X'], height_A=[.45, .55]))
        with self.assertRaisesRegex(ValueError, 'Repeated signed_heights'):
            self.check()
        pose['signed_heights'].pop()
        angle = pose['angles'][0]
        pose['angles'].append(dict(angle, atoms=angle['atoms'][::-1]))
        with self.assertRaisesRegex(ValueError, 'Repeated angles'):
            self.check()

    def test_state_graph_must_contain_each_center_bond(self):
        self.ledger['reactant_bonds'].pop()
        self.ledger['product_bonds'].pop()
        with self.assertRaisesRegex(ValueError, 'covalent neighbor.*absent'):
            self.check()

    def test_unknown_and_invalid_reviews_cannot_silently_pass(self):
        original = self.review
        for kind, field, value, message in [('angles', 'angle_deg', [90, 181], 'interval'),
                ('angles', 'angle_deg', [-1, 2], 'interval'), ('signed_heights', 'height_A', [0, float('nan')], 'interval'),
                ('angles', 'atoms', ['X','B','X'], 'distinct'), ('signed_heights','basis','','basis')]:
            with self.subTest(kind=kind, field=field, value=value):
                self.review = copy.deepcopy(original)
                self.review['poses'][0][kind][0][field] = value
                with self.assertRaisesRegex(ValueError, message):
                    self.check()
        self.review = copy.deepcopy(original)
        self.review['poses'][0]['angel'] = []
        with self.assertRaisesRegex(ValueError, 'Unknown pose'):
            self.check()
        self.review = dict(original, version=1)
        with self.assertRaisesRegex(ValueError, 'require a version-2'):
            self.check()

    def test_planar_center_and_degenerate_anchors(self):
        self.xyz[0, 3] = [0.02, 0.02, 0]
        pose = self.review['poses'][0]
        pose.pop('angles')
        pose['signed_heights'][0]['height_A'] = [-.01, .01]
        self.assertEqual(self.check()[0]['signed_heights'][0]['height_A'], 0)
        self.xyz[0, 2] = (self.xyz[0, 0] + self.xyz[0, 1]) / 2
        with self.assertRaisesRegex(ValueError, 'collinear'):
            self.check()
        with self.assertRaisesRegex(ValueError, 'coincident'):
            measure_nm('angles', [[0,0,0], [0,0,0], [1,0,0]])
        self.assertEqual(measure_nm('angles', [[-1,0,0], [0,0,0], [1,0,0]]), 180)


if __name__ == '__main__':
    unittest.main()
