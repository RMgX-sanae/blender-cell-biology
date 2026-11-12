"""Run with Blender --background --python-exit-code 1 --python this_file.py."""
import sys,unittest
from pathlib import Path
import bpy,numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
from coordinate_animation import keyed_points,verify_point_animation,_curves


class CoordinateHandoff(unittest.TestCase):
    def setUp(self):
        for obj in list(bpy.data.objects):bpy.data.objects.remove(obj,do_unlink=True)
        self.points=np.array([[[0,0,0],[1,2,3],[-2,1,.5]],
                              [[1,3,0],[4,2,-1],[-2,-2,.5]],
                              [[5,4,2],[4,7,-1],[3,-2,1]]],float)

    def make(self,frames=None):
        return keyed_points('Checked points',self.points,bpy.context.scene.collection,[(0,1)],frames=frames)

    def check_evaluated(self,obj,frames,points):
        for t in np.r_[frames[0]-1,np.linspace(frames[0],frames[-1],41),frames[-1]+1]:
            integer=int(np.floor(t));bpy.context.scene.frame_set(integer,subframe=float(t-integer))
            evaluated=obj.evaluated_get(bpy.context.evaluated_depsgraph_get());mesh=evaluated.to_mesh()
            actual=np.array([tuple(v.co) for v in mesh.vertices]);evaluated.to_mesh_clear()
            expected=np.array([np.interp(t,frames,column) for column in points.reshape(len(points),-1).T]).reshape(-1,3)
            np.testing.assert_allclose(actual,expected,atol=3e-6,rtol=0)

    def test_default_playback_is_linear_between_samples(self):
        obj=self.make();verify_point_animation(obj,self.points)
        self.check_evaluated(obj,[1,2,3],self.points)
        self.assertEqual([tuple(edge.vertices) for edge in obj.data.edges],[(0,1)])

    def test_irregular_sample_times_and_endpoint_holds(self):
        frames=[1.,2.25,8.5];obj=self.make(frames)
        verify_point_animation(obj,self.points,frames);self.check_evaluated(obj,frames,self.points)

    def test_single_sample_is_static(self):
        points=self.points[:1];obj=keyed_points('Static',points,bpy.context.scene.collection,frames=[2.5])
        verify_point_animation(obj,points,[2.5]);self.check_evaluated(obj,[2.5],points)

    def test_separate_object_visibility_slot_does_not_change_coordinate_contract(self):
        obj=self.make();obj.hide_render=False;obj.keyframe_insert('hide_render',frame=1)
        self.assertEqual(obj.animation_data.action,obj.data.shape_keys.animation_data.action)
        self.assertNotEqual(obj.animation_data.action_slot,obj.data.shape_keys.animation_data.action_slot)
        verify_point_animation(obj,self.points);self.check_evaluated(obj,[1,2,3],self.points)

    def test_invalid_coordinates_and_times(self):
        bad=[(self.points[:0],None),(self.points[0],None),(self.points*np.nan,None),
             (self.points*1e40,None),(self.points,[1,2]),(self.points,[1,1,2]),
             (self.points,[3,2,1]),(self.points,[1,float('nan'),3]),
             (self.points,[1,1+1e-9,2])]
        for points,frames in bad:
            with self.subTest(frames=frames),self.assertRaises(ValueError):
                keyed_points('Invalid',points,bpy.context.scene.collection,frames=frames)

    def test_blender_time_key_coalescing_cannot_pass(self):
        # Distinct float32 times can still fall inside Blender's insertion threshold.
        with self.assertRaisesRegex(ValueError,'merged coordinate time keys'):
            self.make([1,1.005,2])

    def test_invalid_edges(self):
        for edges in [[(0,3)],[(0,0)],[(0.,1.)],[(-1,1)]]:
            with self.subTest(edges=edges),self.assertRaises(ValueError):
                keyed_points('Invalid edge',self.points,bpy.context.scene.collection,edges)

    def test_changed_coordinates_cannot_pass(self):
        obj=self.make();obj.data.shape_keys.key_blocks[1].data[0].co.x+=.01
        with self.assertRaisesRegex(ValueError,'Coordinate sample changed'):verify_point_animation(obj,self.points)

    def test_reordered_identities_cannot_pass(self):
        obj=self.make();swapped=self.points.copy();swapped[1,[0,1]]=swapped[1,[1,0]]
        with self.assertRaisesRegex(ValueError,'Coordinate sample changed'):verify_point_animation(obj,swapped)

    def test_bezier_time_interpolation_cannot_pass(self):
        obj=self.make();_curves(obj.data.shape_keys)[0].keyframe_points[0].interpolation='BEZIER'
        with self.assertRaisesRegex(ValueError,'linear interpolation'):verify_point_animation(obj,self.points)

    def test_changed_time_key_cannot_pass(self):
        obj=self.make();_curves(obj.data.shape_keys)[0].keyframe_points[1].co.x+=.1
        with self.assertRaisesRegex(ValueError,'timing differs'):verify_point_animation(obj,self.points)

    def test_curve_modifier_cannot_pass(self):
        obj=self.make();_curves(obj.data.shape_keys)[0].modifiers.new('NOISE')
        with self.assertRaisesRegex(ValueError,'linear interpolation'):verify_point_animation(obj,self.points)

    def test_muted_time_curve_cannot_pass(self):
        obj=self.make();_curves(obj.data.shape_keys)[0].mute=True
        with self.assertRaisesRegex(ValueError,'linear interpolation'):verify_point_animation(obj,self.points)


if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CoordinateHandoff))
    if not result.wasSuccessful():raise RuntimeError('Coordinate handoff regression failed')
