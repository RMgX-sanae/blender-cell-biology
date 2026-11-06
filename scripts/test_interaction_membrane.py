"""Whole-body membrane checks that do not need a Blender render."""
import unittest
import numpy as np
from interaction_membrane import plane_clearance


class MembraneChecks(unittest.TestCase):
    def setUp(self):
        self.membrane={'normal':[0,0,1],'origin_nm':[0,0,0],'hydrophobic_half_width_nm':1.2}

    def test_center_can_be_clear_while_body_crosses(self):
        body=np.array([[0,0,-1.],[0,0,1.]])
        pose=np.eye(4)[None];pose[0,2,3]=2
        self.assertGreater(plane_clearance(body.mean(0)[None],pose,self.membrane)[0],0)
        self.assertLess(plane_clearance(body,pose,self.membrane)[0],0)

    def test_envelope_and_intermediate_pose_matter(self):
        poses=np.repeat(np.eye(4)[None],3,axis=0);poses[:,2,3]=[2,1.5,2]
        atoms=np.array([[0.,0,0]])
        self.assertTrue(np.all(plane_clearance(atoms,poses,self.membrane)>0))
        np.testing.assert_allclose(plane_clearance(atoms,poses,self.membrane,.5),[.3,-.2,.3])

    def test_common_rigid_frame_change_preserves_clearance(self):
        body=np.array([[.2,-.5,.4],[1.,.6,-.2],[-.4,.1,.3]])
        poses=np.repeat(np.eye(4)[None],2,axis=0);poses[:,:3,3]=[[0,0,2],[1,2,1.5]]
        turn=np.array([[0,-1,0],[0,0,-1],[1,0,0]])
        shift=np.array([13.,-7.,4.]);changed=poses.copy()
        changed[:,:3,:3]=turn@poses[:,:3,:3];changed[:,:3,3]=poses[:,:3,3]@turn.T+shift
        plane={**self.membrane,'normal':turn@np.array([0.,0,1]),'origin_nm':shift}
        np.testing.assert_allclose(plane_clearance(body,poses,self.membrane,.2),plane_clearance(body,changed,plane,.2))

    def test_bad_plane_is_not_a_pass(self):
        for field,value in [('normal',[0,0,0]),('normal',[0,0,2]),('normal',[0,float('nan'),1]),
                            ('origin_nm',[0,0,float('inf')]),('hydrophobic_half_width_nm',0)]:
            with self.subTest(field=field,value=value),self.assertRaises(ValueError):
                plane_clearance(np.zeros((1,3)),np.eye(4)[None],{**self.membrane,field:value})


if __name__=='__main__':unittest.main()
