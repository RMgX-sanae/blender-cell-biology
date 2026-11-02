"""Explicit project paths shared by the portable actin-film chapter scripts."""
import json
import os
from pathlib import Path

SKILL = Path(__file__).resolve().parents[2]
ROOT = Path(os.environ['BIO_FILM_PROJECT']).resolve()
if ROOT.is_relative_to(SKILL):
    raise ValueError('Generated film projects must be outside the skill directory')
OPTIONS = json.loads((ROOT / 'project.json').read_text())
FFMPEG = OPTIONS['ffmpeg']
FONT = OPTIONS['font_zh']
LATIN_FONT = OPTIONS['font_en']


def validate_render_settings(settings, shot):
    expected = {'engine': OPTIONS['engine'], 'width': OPTIONS['width'],
                'height': 2*round(OPTIONS['width']*9/32), 'percentage': 100,
                'samples': OPTIONS['samples'], 'fps': OPTIONS['fps'], 'fps_base': 1,
                'compositor': OPTIONS['compositor'], 'frame_start': 1, 'frames': shot['frames']}
    differences = {key: {'expected': value, 'actual': settings.get(key)}
                   for key, value in expected.items() if settings.get(key) != value}
    if shot['fps'] != OPTIONS['fps']:
        differences['storyboard_fps'] = {'expected': OPTIONS['fps'], 'actual': shot['fps']}
    if differences:
        raise ValueError('Chapter '+shot['id']+' does not match project render settings; '
                         'rebuild the affected chapter before rendering: '+json.dumps(differences))
