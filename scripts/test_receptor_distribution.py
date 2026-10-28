"""Numerical regression checks for geometry sampling, not biological validation."""
import unittest
import numpy as np
from prepare_receptor_distribution import surface_table, draw_surface, sample, ellipsoid


class DistributionTests(unittest.TestCase):
    def setUp(self):
        vertices = np.array([[0., 0, 0], [2, 0, 0], [0, 1, 0], [10, 0, 0], [18, 0, 0], [10, 1, 0]])
        self.table = surface_table(vertices, np.array([[0, 1, 2], [3, 4, 5]]))

    def test_world_area_not_face_or_uv_count(self):
        points, normals, ids = draw_surface(self.table, 100000, np.random.default_rng(1))
        self.assertLess(abs(np.mean(ids == 1)-.8), .004)
        self.assertTrue(np.allclose(normals, [0, 0, 1]))
        self.assertTrue(np.allclose(points[:, 2], 0))

    def test_side_rigidity_seed_and_last_attempt(self):
        kwargs = dict(min_anchor_distance=0, normal_sign=-1, seed=7, max_attempts=1)
        result = sample(self.table, 1, [[0, 0, 0]], lambda _: 1, 1, **kwargs)
        repeat = sample(self.table, 1, [[0, 0, 0]], lambda _: 1, 1, **kwargs)
        self.assertTrue(np.array_equal(result['positions'], repeat['positions']))
        rotation = result['rotations'][0]
        self.assertTrue(np.allclose(rotation.T@rotation, np.eye(3)))
        self.assertAlmostEqual(np.linalg.det(rotation), 1)
        self.assertTrue(np.allclose(rotation[:, 2], [0, 0, -1]))

    def test_no_silent_count_or_bad_weight(self):
        with self.assertRaisesRegex(ValueError, 'Packing not reached'):
            sample(self.table, 2, [[0, 0, 0]], lambda _: 1, 1, min_anchor_distance=100, max_attempts=10)
        with self.assertRaisesRegex(ValueError, 'Intensity field'):
            sample(self.table, 1, [[0, 0, 0]], lambda _: 2, 1, min_anchor_distance=0)
        with self.assertRaisesRegex(ValueError, 'prototype'):
            sample(self.table, 1, [[np.nan, 0, 0]], lambda _: 1, 1, min_anchor_distance=0)

    def test_whole_body_exclusion_not_just_anchors(self):
        # Every copy includes an atom near the same world origin; anchors alone
        # could fit, but a complete-body screen must reject a second copy.
        tiny = surface_table(np.array([[0., 0, 0], [.02, 0, 0], [0, .02, 0]]), np.array([[0, 1, 2]]))
        with self.assertRaisesRegex(ValueError, '1/2'):
            sample(tiny, 2, [[0, 0, 0], [0, 0, 5]], lambda _: 1, 1,
                   min_anchor_distance=0, atom_clearance=.2, max_attempts=20)

    def test_ellipsoid_outward_after_nonuniform_scale(self):
        vertices, faces = ellipsoid([4, 2, 1], 8, 16)
        table = surface_table(vertices, faces)
        self.assertTrue(np.all(np.einsum('ij,ij->i', table[0].mean(1), table[1]) > 0))


if __name__ == '__main__': unittest.main()
