"""Numerical regression checks for geometry sampling, not biological validation."""
import unittest
import numpy as np
from prepare_receptor_distribution import surface_table, draw_surface, sample, ellipsoid, validate_membrane_surface


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


class MeshInputTests(unittest.TestCase):
    def setUp(self):
        self.vertices, self.faces = ellipsoid([4, 2, 1], 8, 16)

    def test_closed_and_open_surface_contract(self):
        audit = validate_membrane_surface(self.vertices, self.faces, True)
        self.assertEqual(audit['boundary_edges'], 0)
        self.assertGreater(audit['components'][0]['signed_volume_nm3'], 0)
        # Removing one triangle leaves a hole without orphaning any vertex.
        audit = validate_membrane_surface(self.vertices, self.faces[1:])
        self.assertEqual(audit['boundary_edges'], 3)
        self.assertIsNone(audit['components'][0]['signed_volume_nm3'])
        with self.assertRaisesRegex(ValueError, 'boundary edges'):
            validate_membrane_surface(self.vertices, self.faces[1:], True)

    def test_inconsistent_or_reversed_winding(self):
        faces = self.faces.copy(); faces[0] = faces[0, ::-1]
        with self.assertRaisesRegex(ValueError, 'Inconsistent'):
            validate_membrane_surface(self.vertices, faces)
        with self.assertRaisesRegex(ValueError, 'outward face winding'):
            validate_membrane_surface(self.vertices, self.faces[:, ::-1])
        reflected = self.vertices*np.array([-1, 1, 1])
        with self.assertRaisesRegex(ValueError, 'outward face winding'):
            validate_membrane_surface(reflected, self.faces)
        validate_membrane_surface(reflected, self.faces[:, ::-1], True)

    def test_nonmanifold_and_degenerate_faces(self):
        with self.assertRaisesRegex(ValueError, 'Nonmanifold'):
            validate_membrane_surface(self.vertices, np.vstack([self.faces, self.faces[0]]))
        with self.assertRaisesRegex(ValueError, 'degenerate'):
            validate_membrane_surface(self.vertices, np.vstack([self.faces, [0, 0, 1]]))

    def test_each_disconnected_component_requires_outward_winding(self):
        # The combined volume stays positive; a single global sum would miss this.
        vertices = np.vstack([self.vertices, self.vertices*.5+[20, 0, 0]])
        faces = np.vstack([self.faces, self.faces[:, ::-1]+len(self.vertices)])
        with self.assertRaisesRegex(ValueError, 'outward face winding'):
            validate_membrane_surface(vertices, faces, True)
        faces[len(self.faces):] = self.faces+len(self.vertices)
        self.assertEqual(len(validate_membrane_surface(vertices, faces, True)['components']), 2)

    def test_affine_mesh_uses_world_normals_without_stretching_receptors(self):
        angle = .6
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                             [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        transform = rotation@np.diag([1.4, .7, 1.2])
        vertices = self.vertices@transform.T+[11, -7, 3]
        validate_membrane_surface(vertices, self.faces, True)
        expected = surface_table(self.vertices, self.faces)[1]@np.linalg.inv(transform)
        expected /= np.linalg.norm(expected, axis=1)[:, None]
        placed = sample(surface_table(vertices, self.faces), 16, [[0, 0, 0]], lambda _: 1, 1,
                        min_anchor_distance=.3, normal_sign=-1)
        for matrix, face in zip(placed['rotations'], placed['face_ids']):
            self.assertTrue(np.allclose(matrix[:, 2], -expected[face]))
            self.assertTrue(np.allclose(matrix.T@matrix, np.eye(3)))
            self.assertAlmostEqual(np.linalg.det(matrix), 1)


if __name__ == '__main__': unittest.main()
