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
