"""Check that both molecular panels and the title contribute to a saved sequence render."""
import argparse
import json
from pathlib import Path
import sys
import tempfile

import bpy
import numpy as np


def render_pixels(scene, path):
    scene.render.filepath = str(path)
    bpy.ops.render.render(write_still=True, scene=scene.name)
    image = bpy.data.images.load(str(path), check_existing=False)
    try:
        image.colorspace_settings.name = 'Non-Color'
        width, height = image.size
        pixels = np.empty(width * height * 4, dtype=np.float32)
        image.pixels.foreach_get(pixels)
        return pixels.reshape(height, width, 4)[:, :, :3]
    finally:
        bpy.data.images.remove(image)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--frame', type=int, help='Use the saved opening frame by default')
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    scene = bpy.context.scene
    if not bpy.data.filepath or 'sequence_display_sha256' not in scene:
        parser.error('Open a saved two-panel molecular sequence first')
    frame = scene.frame_current if args.frame is None else args.frame
    if not scene.frame_start <= frame <= scene.frame_end:
        parser.error('Review frame lies outside the saved presentation range')
    out = args.directory.resolve()
    out.mkdir(parents=True, exist_ok=False)
    collections = [scene.collection.children[name] for name in
                   ('Sequence overview', 'Sequence detail', scene['compositor_label_collection'])]
    original_flags = [c.hide_render for c in collections]
    settings = scene.render.image_settings
    original_output = scene.render.filepath, settings.file_format, settings.color_depth
    original_time = scene.frame_current, scene.frame_subframe
    try:
        scene.frame_set(frame)
        settings.file_format = 'PNG'
        settings.color_depth = '8'
        intact = render_pixels(scene, out / 'frame.png')
        for collection in collections:
            collection.hide_render = True
        with tempfile.TemporaryDirectory(dir=out) as cache:
            control = render_pixels(scene, Path(cache) / 'without_subjects.png')
    finally:
        for collection, hidden in zip(collections, original_flags):
            collection.hide_render = hidden
        scene.render.filepath, settings.file_format, settings.color_depth = original_output
        scene.frame_set(original_time[0], subframe=original_time[1])

    # Blender image rows run from bottom to top. Leave the divider, caption
    # and legend outside the panel regions. Markers stay in the control so a
    # dot or a bond annotation alone cannot stand in for the molecular body.
    height, width = intact.shape[:2]
    delta = np.max(np.abs(intact - control), axis=2)
    bounds = {'overview': (.045, .46, .23, .78), 'detail': (.54, .955, .23, .78),
              'title': (.10, .90, .88, .99)}
    regions = {}
    for label, (left, right, bottom, top) in bounds.items():
        area = delta[round(bottom*height):round(top*height), round(left*width):round(right*width)]
        changed = int(np.count_nonzero(area >= 8 / 255))
        minimum = max(8, int(np.ceil(area.size * .001)))
        regions[label] = {'rectangle_fraction': [left, right, bottom, top], 'changed_pixels': changed,
                          'minimum_changed_pixels': minimum, 'passed': changed >= minimum}
    report = {'passed': all(r['passed'] for r in regions.values()), 'frame': frame,
              'scene_file': bpy.data.filepath, 'blender': bpy.app.version_string,
              'size': [width, height], 'display_sha256': scene['sequence_display_sha256'],
              'channel_difference_threshold_8bit': 8, 'regions': regions,
              'scope': 'Pixel contribution in the standard two molecular panels and central title. '
                       'Inspect shape, identity, occlusion, individual glyphs and scientific meaning separately.'}
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    if not report['passed']:
        raise ValueError('Missing rendered contribution in ' + ', '.join(k for k,v in regions.items() if not v['passed']))


if __name__ == '__main__':
    main()
