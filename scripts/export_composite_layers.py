"""Export a published layout and protected labels as linear EXR for compositing.

Run inside Blender with an existing scene. No geometry or source file is changed.
The current interface supports Blender 5.2's File Output node and built-in OCIO.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

import bpy

sys.path.insert(0, str(Path(__file__).resolve().parent))
from render_saved import selected_frames


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def color_state(scene):
    root = Path(bpy.utils.resource_path('LOCAL'))/'datafiles/colormanagement'
    config = os.environ.get('OCIO')
    if (config and Path(config).resolve() != (root/'config.ocio').resolve()) or scene.view_settings.use_curve_mapping:
        raise ValueError('Qualify custom OCIO/curves before exporting this handoff')
    return {'view': scene.view_settings.view_transform, 'look': scene.view_settings.look,
            'exposure': scene.view_settings.exposure, 'gamma': scene.view_settings.gamma,
            'display': scene.display_settings.display_device, 'dither': scene.render.dither_intensity,
            'linear_space': 'Linear Rec.709', 'blender': bpy.app.version_string,
            'ocio_files_sha256': {str(p.relative_to(root)): digest(p) for p in sorted(root.rglob('*')) if p.is_file()}}


def published_sources(tree):
    if tree is None or tree.get('display_contract_version') != 1 or not tree.get('display_has_labels'):
        raise ValueError('Require a published Beauty layout and separate protected Labels')
    outputs = [n for n in tree.nodes if n.bl_idname == 'NodeGroupOutput' and n.is_active_output]
    if len(outputs) != 1: raise ValueError('Require one active compositor output')
    ports = {}
    for name in ('Beauty', 'Labels', 'Image'):
        socket = outputs[0].inputs.get(name)
        if socket is None or len(socket.links) != 1: raise ValueError('Missing published output: '+name)
        ports[name] = socket.links[0].from_socket
    merge = ports['Image'].node
    if (merge.bl_idname != 'CompositorNodeAlphaOver' or merge.inputs['Type'].default_value != 'Over'
            or merge.inputs['Straight Alpha'].default_value or merge.inputs['Factor'].is_linked
            or merge.inputs['Factor'].default_value != 1):
        raise ValueError('Require a premultiplied Labels-over-Beauty final composition')
    for name, port in [('Background', ports['Beauty']), ('Foreground', ports['Labels'])]:
        if len(merge.inputs[name].links) != 1 or merge.inputs[name].links[0].from_socket != port:
            raise ValueError('Published layers do not describe the final composition')
    return ports


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--frames', required=True, help='Explicit frames/ranges; use all only for a full deliverable')
    p.add_argument('--width', type=int)
    p.add_argument('--samples', type=int)
    a = p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.output.exists(): p.error('Choose a fresh export directory')
    if not bpy.data.filepath: p.error('Load the qualified source scene first')
    scene = bpy.context.scene; source = Path(bpy.data.filepath); source_hash = digest(source)
    frames = selected_frames(a.frames, scene.frame_start, scene.frame_end)
    colors = color_state(scene)
    ports = published_sources(scene.compositing_node_group)
    tree = scene.compositing_node_group.copy(); tree.name = 'External composite export'
    scene.compositing_node_group = tree; ports = published_sources(tree)
    for node in tree.nodes:
        if node.bl_idname == 'CompositorNodeOutputFile': node.mute = True
    if a.width:
        ratio = scene.render.resolution_y/scene.render.resolution_x
        scene.render.resolution_x = a.width; scene.render.resolution_y = 2*round(a.width*ratio/2)
    if a.samples:
        scene.cycles.samples = a.samples
        if hasattr(scene, 'eevee'): scene.eevee.taa_render_samples = a.samples
    scene.render.resolution_percentage = 100
    nodes = {}
    for name, socket in ports.items():
        node = tree.nodes.new('CompositorNodeOutputFile'); node.name = 'Export '+name
        node.directory = str(a.output/name); node.file_name = 'frame_####'
        node.format.media_type = 'IMAGE'; node.format.file_format = 'OPEN_EXR'; node.format.color_mode = 'RGBA'
        node.format.color_depth = '32'; node.format.exr_codec = 'ZIP'
        node.format.linear_colorspace_settings.name = colors['linear_space']
        node.save_as_render = False
        node.file_output_items.new('RGBA', 'rgba'); tree.links.new(socket, node.inputs['rgba'])
        nodes[name] = node
    a.output.mkdir(parents=True)
    scene.render.image_settings.file_format = 'PNG'; scene.render.image_settings.color_mode = 'RGBA'
    scene.render.image_settings.color_depth = '8'
    for frame in frames:
        scene.frame_set(frame); scene.render.filepath = str(a.output/f'reference_{frame:04d}.png')
        bpy.ops.render.render(write_still=True)
    files = {str(f.relative_to(a.output)): digest(f) for f in sorted(a.output.rglob('*')) if f.is_file()}
    expected = {f'{name}/frame_{frame:04d}rgba.exr' for name in ports for frame in frames}
    expected.update(f'reference_{frame:04d}.png' for frame in frames)
    if set(files) != expected: raise ValueError('Unexpected exported frame/layer filenames')
    if digest(source) != source_hash: raise ValueError('Source scene file changed during export')
    record = {'kind': 'linear-composite-handoff', 'version': 1, 'source_blend': str(source),
              'source_sha256': source_hash, 'source_unchanged': True, 'frames': frames,
              'fps': scene.render.fps, 'fps_base': scene.render.fps_base,
              'resolution': [scene.render.resolution_x, scene.render.resolution_y],
              'colors': colors, 'alpha': 'premultiplied', 'files_sha256': files,
              'patterns': {name: f'{name}/frame_####rgba.exr' for name in ports},
              'scope': 'Existing complete layout, protected labels and original composite; no geometry, motion or new mechanism.',
              'external_compositor_executed': False}
    (a.output/'handoff.json').write_text(json.dumps(record, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'output': str(a.output), 'frames': frames, 'files': list(files)}))


if __name__ == '__main__': main()
