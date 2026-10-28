"""Apply an editable compositor stack to the scene currently open in Blender."""
import argparse
import json
from pathlib import Path
import sys
import bpy

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compositor_stack import build_stack


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--recipes', type=Path, default=Path(__file__).resolve().parents[1]/'assets/render-recipes.json')
    p.add_argument('--recipe', required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--engine', choices=('cycles', 'eevee'))
    p.add_argument('--samples', type=int)
    p.add_argument('--width', type=int)
    p.add_argument('--frame', type=int, default=1)
    p.add_argument('--subject', action='append', default=[], help='Exact object or collection name; assigns whole objects index 1')
    p.add_argument('--render', action='store_true')
    a = p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.output.exists(): p.error('Use a fresh output directory')
    recipe = next((r for r in json.loads(a.recipes.read_text())['recipes'] if r['id'] == a.recipe), None)
    if recipe is None: p.error('Recipe not found in the specified file')
    scene = bpy.context.scene
    if a.engine: scene.render.engine = 'CYCLES' if a.engine == 'cycles' else 'BLENDER_EEVEE'
    if a.width:
        ratio = scene.render.resolution_y/scene.render.resolution_x
        scene.render.resolution_x = a.width; scene.render.resolution_y = 2*round(a.width*ratio/2)
    if a.samples:
        scene.cycles.samples = a.samples
        if hasattr(scene, 'eevee'): scene.eevee.taa_render_samples = a.samples
    selected=set()
    for name in a.subject:
        if name in bpy.data.collections:
            objects = bpy.data.collections[name].all_objects
        elif name in scene.objects:
            objects = [scene.objects[name]]
        else: p.error('Unknown subject object/collection: ' + name)
        for obj in objects: obj.pass_index = 1; selected.add(obj.name)
    for obj in scene.objects:
        if obj.get('legend_for_objects') and selected.intersection(json.loads(obj['legend_for_objects'])):
            obj.pass_index=1
    scene.frame_set(a.frame); bpy.context.view_layer.update()
    source = bpy.data.filepath
    report = build_stack(scene, recipe)
    a.output.mkdir(parents=True)
    scene.render.image_settings.file_format = 'PNG'
    scene.render.filepath = str(a.output/f'frame_{a.frame:04d}.png')
    bpy.ops.wm.save_as_mainfile(filepath=str(a.output/'scene.blend'))
    if a.render: bpy.ops.render.render(write_still=True)
    report.update(source_blend=source, recipe=recipe, blender=bpy.app.version_string, frame=a.frame,
                  resolution=[scene.render.resolution_x,scene.render.resolution_y], rendered=a.render)
    (a.output/'render-recipe-report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__': main()
