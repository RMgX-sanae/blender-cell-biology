"""Event order and caption synchronization checks; no kinetic validation."""
import copy
import json
from pathlib import Path
import unittest

import numpy as np
from event_timing import continuous_intervals, phase_progress, validate_intervals, frame_window
from prepare_dub_hydrolysis import presentation_config


class EventTimingTests(unittest.TestCase):
    def test_frame_window_does_not_announce_an_event_one_frame_early(self):
        self.assertEqual(frame_window(.52,.7,720), (376,504))
        self.assertEqual(frame_window(0,.52,720), (1,375))
        self.assertEqual(frame_window(0,1,720), (1,720))
        # A narrow interval may contain no sampled frame; never invent one.
        first,last=frame_window(.5001,.5002,720)
        self.assertGreater(first,last)

    def test_ras_existing_phases_remain_bit_identical(self):
        time = np.linspace(0, 1, 577)
        for intervals in ({'binding':[0,.36], 'hydrolysis':[.36,.52], 'release':[.58,1]},
                          {'binding':[.06,.36], 'hydrolysis':[.43,.57], 'release':[.65,.95]}):
            result = phase_progress(time, intervals, [('binding','hydrolysis'), ('hydrolysis','release')])
            for name, (start, end) in intervals.items():
                self.assertTrue(np.array_equal(result[name], np.clip((time-start)/(end-start), 0, 1)))

    def test_invalid_intervals_and_event_order_are_rejected(self):
        for pair in ([0,0], [.8,.2], [-.1,.2], [0,1.1], [0,np.nan], [0,np.inf], [0,.5,1]):
            with self.subTest(pair=pair), self.assertRaises(ValueError):
                validate_intervals({'reaction':pair})
        with self.assertRaisesRegex(ValueError, 'finish before'):
            validate_intervals({'reaction':[0,.6], 'release':[.5,1]}, [('reaction','release')])
        with self.assertRaisesRegex(ValueError, 'Unknown event'):
            validate_intervals({'reaction':[0,1]}, [('binding','reaction')])
        with self.assertRaisesRegex(ValueError, 'nondecreasing'):
            phase_progress([0,.8,.2], {'reaction':[0,1]})

    def test_continuous_removes_only_holds_and_preserves_duration_ratios(self):
        result = continuous_intervals({'a':[.1,.3], 'b':[.5,.9]}, ['a','b'])
        self.assertEqual(result['a'][0], 0)
        self.assertEqual(result['b'][1], 1)
        self.assertEqual(result['a'][1], result['b'][0])
        self.assertAlmostEqual(result['a'][1], 1/3)
        with self.assertRaisesRegex(ValueError, 'every event exactly once'):
            continuous_intervals({'a':[0,.3], 'b':[.5,1]}, ['a','a'])


class DubPresentationTests(unittest.TestCase):
    def setUp(self):
        self.config = json.loads((Path(__file__).resolve().parents[1]/'assets/dub-hydrolysis.json').read_text())

    def test_default_teaching_caption_boundaries_unchanged(self):
        result = presentation_config(self.config, 'teaching')
        expected = [[0,.08], [.08,.48], [.48,.52], [.52,.7], [.7,.98], [.98,1]]
        self.assertEqual([row[:2] for row in result['captions']], expected)
        self.assertNotIn('captions', self.config)

    def test_changed_timing_moves_captions_and_supports_old_project_texts(self):
        config = copy.deepcopy(self.config)
        config['timing'] = {'reaction':[.08,.60], 'proximal_release':[.64,.76], 'distal_release':[.78,.98]}
        expected_texts = config.pop('caption_texts')
        config['captions'] = [[0,1,text] for text in expected_texts]
        result = presentation_config(config, 'teaching')
        def label(time): return next(text for start,end,text in result['captions'] if start <= time < end)
        self.assertEqual(label(.56), expected_texts[1])
        self.assertEqual(label(.62), expected_texts[2])
        self.assertEqual(label(.65), expected_texts[3])
        self.assertEqual(result['captions'][3][:2], [.64,.78])

    def test_continuous_caps_have_no_overlap_gaps_or_zero_length(self):
        result = presentation_config(self.config, 'continuous')
        self.assertEqual(len(result['captions']), 3)
        for t in np.arange(720)/720:
            self.assertEqual(sum(start <= t < end for start,end,_ in result['captions']), 1)
        for i, name in enumerate(['reaction','proximal_release','distal_release']):
            self.assertEqual(result['captions'][i][:2], result['timing'][name])

    def test_dub_chemical_order_and_loop_boundary_are_explicit(self):
        self.config['timing']['proximal_release'] = [.4,.69]
        with self.assertRaisesRegex(ValueError, 'reaction -> proximal_release'):
            presentation_config(self.config, 'teaching')
        self.config['timing']['proximal_release'] = [.52,.69]
        with self.assertRaisesRegex(ValueError, 'not a closed catalytic loop'):
            presentation_config(self.config, 'loop')


if __name__ == '__main__': unittest.main()
