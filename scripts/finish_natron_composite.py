"""Verify Natron's actual linear pixels and apply the saved Blender display once.

Run inside Blender with --factory-startup. A saved/reopened graph alone is not
proof of rendering: every requested frame must exist and satisfy the recipe.
"""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from export_composite_layers import color_state, digest


def read_linear(path, resolution):
    image = bpy.data.images.load(str(path), check_existing=False)
    if list(image.size) != resolution or not image.is_float:
        raise ValueError('Expected full-size floating-point image: '+str(path))
    image.colorspace_settings.name = 'Linear Rec.709'; image.alpha_mode = 'PREMUL'
    pixels = np.empty(resolution[0]*resolution[1]*4, dtype=np.float32)
    image.pixels.foreach_get(pixels)
    pixels = pixels.reshape(resolution[1], resolution[0], 4)
    if not np.isfinite(pixels).all(): raise ValueError('Nonfinite pixels: '+str(path))
    return image, pixels


def matches(actual, expected, context):
    error = float(np.max(np.abs(actual-expected)))
    if not np.allclose(actual, expected, atol=2e-6, rtol=2e-6):
        raise ValueError(context+' pixel mismatch: '+str(error))
    return error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--')+1:])
    job = json.loads(args.job.read_text()); output = Path(job['output'])
    handoff_path = Path(job['handoff']); inputs = handoff_path.parent
    handoff = json.loads(handoff_path.read_text()); colors = handoff['colors']
    prepared = json.loads((output/'natron-graph.json').read_text())
    audit = json.loads((output/'reopened-audit.json').read_text())
    if (not audit['passed'] or audit['project_sha256'] != digest(output/'scene.ntp')
            or audit['handoff_sha256'] != digest(handoff_path)
            or prepared['handoff_sha256'] != digest(handoff_path)
            or prepared['graph'] != audit['graph'] or job['grade'] != prepared['grade']):
        raise ValueError('Saved project, recipe or handoff changed after graph verification')
    if digest(job['ocio_config']) != prepared['ocio_sha256']:
        raise ValueError('Natron OCIO configuration changed')
    for name, expected in handoff['files_sha256'].items():
        if digest(inputs/name) != expected: raise ValueError('Input changed: '+name)
    scene = bpy.context.scene
    runtime = color_state(scene)
    if runtime['ocio_files_sha256'] != colors['ocio_files_sha256'] or runtime['blender'] != colors['blender']:
        raise ValueError('Use the source Blender version and color resources for display conversion')
    scene.display_settings.display_device = colors['display']
    for property_name, key in [('view_transform','view'), ('look','look'), ('exposure','exposure'), ('gamma','gamma')]:
        setattr(scene.view_settings, property_name, colors[key])
    scene.render.dither_intensity = colors['dither']
    scene.render.image_settings.file_format = 'PNG'; scene.render.image_settings.color_mode = 'RGBA'
    scene.render.image_settings.color_depth = '8'
    gain = np.array(job['grade']['multiply'], dtype=np.float32)
    gamma = np.array(job['grade']['gamma'], dtype=np.float32)
    rows = []
    for frame in handoff['frames']:
        images = {}; arrays = {}
        for role in ('Beauty','Labels','Image'):
            path = inputs/handoff['patterns'][role].replace('####', f'{frame:04d}')
            images[role], arrays[role] = read_linear(path, handoff['resolution'])
        image_path = output/f'linear_{frame:04d}.exr'
        actual_image, actual = read_linear(image_path, handoff['resolution'])
        beauty, labels = arrays['Beauty'], arrays['Labels']
        base = labels+beauty*(1-labels[..., 3:4])
        base_error = matches(base, arrays['Image'], 'Published layers')
        alpha = beauty[..., 3:4]
        straight = np.divide(beauty[..., :3], alpha, out=np.zeros_like(beauty[..., :3]), where=alpha > 0)
        if np.any(straight < 0) and np.any(gamma != 1):
            raise ValueError('Qualify negative RGB grading separately before using non-unit gamma')
        graded = beauty.copy(); graded[..., :3] = np.power(straight*gain, 1/gamma)*alpha
        expected = labels+graded*(1-labels[..., 3:4])
        error = matches(actual, expected, 'Natron grade and merge')
        label_error = matches(actual[..., :3]-graded[..., :3]*(1-labels[..., 3:4]), labels[..., :3], 'Protected labels')
        preview = output/f'display_{frame:04d}.png'
        actual_image.save_render(str(preview), scene=scene)
        rows.append({'frame': frame, 'base_merge_max_error': base_error,
                     'recipe_max_error': error, 'label_contribution_max_error': label_error,
                     'original_difference_max': float(np.max(np.abs(actual-arrays['Image']))),
                     'linear_sha256': digest(image_path), 'display_sha256': digest(preview)})
        for image in list(images.values())+[actual_image]: bpy.data.images.remove(image)
    record = {'passed': True, 'frames': rows, 'resolution': handoff['resolution'],
              'project_sha256': digest(output/'scene.ntp'), 'handoff_sha256': digest(handoff_path),
              'display': colors, 'grade': job['grade'],
              'scope': 'Verified sampled linear frames, premultiplied merge and protected label contribution; display transform applied once. No full-film or motion validation.'}
    (output/'composite-check.json').write_text(json.dumps(record, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'passed': True, 'frames': handoff['frames'], 'output': str(output)}))


if __name__ == '__main__': main()
