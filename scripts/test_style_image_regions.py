"""Image-region regressions for display-space crops and signed linear luminance."""
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from PIL import Image, ImageCms

from analyze_style_image import local_contrast, measure


class ColorRegions(unittest.TestCase):
    def setUp(self):
        temp = TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.file = Path(temp.name) / 'appearance.png'
        self.im = Image.new('RGBA', (4, 1))
        self.im.putdata([(255, 255, 255, 255), (128, 128, 128, 255), (0, 0, 0, 255), (255, 255, 255, 0)])
        self.im.save(self.file)
        self.regions = {'subject': [0, 0, 1, 1], 'background': [1, 0, 2, 1]}

    def test_local_linear_measurement_with_either_polarity(self):
        midgray = ((128 / 255 + .055) / 1.055) ** 2.4
        for src, bg, expected in ((0, 1, 1 - midgray), (1, 0, midgray - 1), (2, 1, -midgray)):
            areas = dict(subject=[src, 0, src + 1, 1], background=[bg, 0, bg + 1, 1])
            self.assertAlmostEqual(local_contrast(self.file, areas)['subject_minus_background_Y'], expected)
        self.assertEqual(measure(self.file, [0, 0, 1, 1])['strong_edge_fraction'], 0.)

    def test_alpha_and_out_of_bounds_do_not_supply_a_fictitious_background(self):
        for bad in ([-1, 0, 1, 1], [0, 0, 5, 1], [1, 0, 1, 1]):
            with self.subTest(bad=bad), self.assertRaisesRegex(ValueError, 'rectangle'):
                measure(self.file, bad)
        with self.assertRaisesRegex(ValueError, 'translucent'):
            local_contrast(self.file, dict(self.regions, background=[2, 0, 4, 1]))

    def test_embedded_srgb_matches_an_opaque_source(self):
        before = local_contrast(self.file, self.regions)
        self.im.save(self.file, icc_profile=ImageCms.ImageCmsProfile(ImageCms.createProfile('sRGB')).tobytes())
        after = local_contrast(self.file, self.regions)
        self.assertEqual(after['regions'], before['regions'])
        self.assertIn('ICC', after['color_input'])


if __name__ == '__main__':
    unittest.main()
