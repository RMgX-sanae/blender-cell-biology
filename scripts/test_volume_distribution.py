"""Behavior checks for complete-body placement, including an independent oracle."""
import unittest
import numpy as np
from volume_distribution import sample,validate,ball_inside,enclosing_sphere


class VolumePlacementTests(unittest.TestCase):
    def test_anisotropic_grid_seed_and_rigid_bodies(self):
        mask=np.zeros((21,29,35),bool);mask[1:-1,1:-1,1:-1]=True
        mask[6:15,8:21,11:24]=False
        spacing=np.array([.13,.19,.31]);origin=np.array([4.,-3.,8.])
        body=np.array([[-.09,-.06,-.08],[.09,.06,.08],[-.09,.06,-.08],[.09,-.06,.08]])+[10.,-20.,15.]
        a=sample(mask,spacing,origin,body,30,seed=86,boundary_gap=.015,pair_gap=.03)
        b=sample(mask,spacing,origin,body,30,seed=86,boundary_gap=.015,pair_gap=.03)
        self.assertTrue(np.array_equal(a['transforms'],b['transforms']))
        self.assertEqual(validate(mask,spacing,origin,body,a['transforms'],boundary_gap=.015,pair_gap=.03)['count'],30)
        self.assertTrue(all(ball_inside(mask,spacing,origin,c,a['body_radius']+.015) for c in a['centers_xyz']))

    def test_center_inside_does_not_establish_body_containment(self):
        mask=np.ones((9,9,9),bool)
        self.assertTrue(mask[4,4,0])
        self.assertFalse(ball_inside(mask,[1,1,1],[0,0,0],[0,4,4],.6))
        self.assertTrue(ball_inside(mask,[1,1,1],[0,0,0],[0,4,4],.4))

    def test_vertices_and_center_can_miss_internal_obstacle(self):
        mask=np.ones((9,9,9),bool);mask[4,4,5]=False
        body=np.array([[-3.,0,0],[3.,0,0]])
        matrix=np.eye(4);matrix[:3,3]=[4,4,4]
        self.assertTrue(mask[4,4,1] and mask[4,4,7] and mask[4,4,4])
        with self.assertRaisesRegex(ValueError,'forbidden'):validate(mask,[1,1,1],[0,0,0],body,[matrix])

    def test_voxel_half_width_and_jitter_matter(self):
        mask=np.ones((7,7,7),bool);mask[3,3,4]=False
        # Distance to the forbidden center is .6, but distance to its box is .1.
        self.assertFalse(ball_inside(mask,[1,1,1],[0,0,0],[3.4,3,3],.2))
        self.assertTrue(ball_inside(mask,[1,1,1],[0,0,0],[3.4,3,3],.05))

    def test_bad_transforms_and_pair_overlap_rejected(self):
        mask=np.ones((11,11,11),bool);body=np.array([[-.2,0,0],[.2,0,0]])
        matrix=np.eye(4);matrix[:3,3]=[5,5,5]
        for factor in [2.,-1.]:
            broken=matrix.copy();broken[0,0]=factor
            with self.assertRaisesRegex(ValueError,'stretches or reflects'):validate(mask,[1,1,1],[0,0,0],body,[broken])
        with self.assertRaisesRegex(ValueError,'pair separation'):validate(mask,[1,1,1],[0,0,0],body,[matrix,matrix])

    def test_unreachable_placement_and_input_contract(self):
        mask=np.ones((5,5,5),bool);body=np.array([[-10.,0,0],[10.,0,0]])
        with self.assertRaisesRegex(ValueError,'No conservative'):sample(mask,[1,1,1],[0,0,0],body,1)
        with self.assertRaisesRegex(ValueError,'boolean'):sample(mask.astype('uint8'),[1,1,1],[0,0,0],body,1)
        with self.assertRaisesRegex(ValueError,'spacing'):sample(mask,[1,0,1],[0,0,0],body,1)
        self.assertGreater(enclosing_sphere(body)[1],0)


if __name__=='__main__':unittest.main()
