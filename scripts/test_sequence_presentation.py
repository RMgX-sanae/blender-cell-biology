"""Presentation changes retain the qualified path, event order and observation states."""
import copy
import unittest
import numpy as np
from sequence_presentation import presentation,source_positions,coordinates_at,caption_windows


class PresentationTests(unittest.TestCase):
    def setUp(self):
        self.source={'frames':4,'fps':24}
        self.profile={'version':1,'mode':'teaching','source_sha256':'source-A','fps':24,'frames':8,
                      'anchors':[[0,1],[2,3],[2,7],[4,9]]}

    def resolve(self,profile=None):
        return presentation(self.source,self.profile if profile is None else profile,'source-A')

    def test_default_samples_and_caption_windows_unchanged(self):
        result=presentation(self.source)
        self.assertEqual(result['source_samples'],[0,1,2,3,4])
        self.assertEqual(result['sample_frames'],[1,2,3,4,5])
        self.assertEqual(caption_windows([[0,.5,'before'],[.5,1,'after']],result,4),[[1,2,'before'],[3,4,'after']])

    def test_teaching_hold_retains_every_sample_and_corner(self):
        result=self.resolve()
        self.assertEqual(result['source_samples'],[0,1,2,2,3,4])
        np.testing.assert_array_equal(source_positions(result,[1,2,3,4,5,6,7,8,9]),[0,1,2,2,2,2,2,3,4])
        points=np.array([[[0,0,0]],[[1,0,0]],[[1,1,0]],[[1,1,1]],[[2,1,1]]])
        sampled=coordinates_at(points,source_positions(result,[2.5,3,5,7,7.5]))
        np.testing.assert_allclose(sampled[:,0],[[1,.5,0],[1,1,0],[1,1,0],[1,1,0],[1,1,.5]])
        self.assertEqual(caption_windows([[0,.5,'before'],[.5,1,'after']],result,4),[[1,2,'before'],[3,8,'after']])

    def test_continuous_retains_nonuniform_keys_without_teaching_text(self):
        profile={**self.profile,'mode':'continuous','frames':6,'anchors':[[0,1],[1,2],[4,7]]}
        result=self.resolve(profile)
        self.assertEqual(result['source_samples'],[0,1,2,3,4])
        np.testing.assert_allclose(result['sample_frames'],[1,2,11/3,16/3,7])
        self.assertEqual(result['captions'],'none')
        self.assertEqual(caption_windows([],result,4),[])

    def test_source_mismatch_and_loop_claim_rejected(self):
        for change,diagnostic in [({'source_sha256':'other'},'qualified source'),({'mode':'loop'},'not a closed loop')]:
            with self.subTest(change=change),self.assertRaisesRegex(ValueError,diagnostic):self.resolve({**self.profile,**change})

    def test_backtracking_missing_endpoints_and_noninteger_anchors_rejected(self):
        for anchors in [[[1,1],[4,9]],[[0,1],[3,9]],[[0,1],[4,8]],[[0,1],[3,4],[2,5],[4,9]],
                        [[0,1],[2,4],[3,4],[4,9]],[[0,1],[2.2,4],[4,9]],[[0,1],[True,4],[4,9]],[[0,1],[4,float('nan')]]]:
            with self.subTest(anchors=anchors),self.assertRaises(ValueError):self.resolve({**self.profile,'anchors':anchors})

    def test_continuous_cannot_insert_teaching_hold(self):
        with self.assertRaisesRegex(ValueError,'holds require teaching'):self.resolve({**self.profile,'mode':'continuous'})

    def test_invalid_timing_controls_rejected(self):
        for change in [{'fps':0},{'fps':True},{'frames':8.5},{'version':True},{'captions':'automatic'},{'frame_rate':24}]:
            with self.subTest(change=change),self.assertRaises(ValueError):self.resolve({**self.profile,**change})

    def test_caption_gaps_overlaps_and_invisible_events_rejected(self):
        for captions in [[[0,.4,'a'],[.6,1,'b']],[[0,.8,'a'],[.5,1,'b']],[[0,.01,'a'],[.01,.02,'b'],[.02,1,'c']]]:
            with self.subTest(captions=captions),self.assertRaises(ValueError):caption_windows(captions,self.resolve(),4)

    def test_out_of_range_queries_are_not_silently_clamped(self):
        for values in [[0],[10],[float('nan')]]:
            with self.subTest(values=values),self.assertRaises(ValueError):source_positions(self.resolve(),values)
        with self.assertRaises(ValueError):coordinates_at(np.zeros((5,1,3)),[4.01])

    def test_input_profiles_are_not_modified(self):
        original=copy.deepcopy(self.profile);self.resolve();self.assertEqual(original,self.profile)


if __name__=='__main__':unittest.main()
