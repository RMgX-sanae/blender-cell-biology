"""Render selected/all frames from the .blend loaded by Blender, without resaving it."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import bpy
sys.path.insert(0,str(Path(__file__).resolve().parent))
from render_cache import valid_png


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--frames',default='all',help='all or comma-separated integer frames')
    p.add_argument('--width',type=int)
    p.add_argument('--compositor',help='Exact saved COMP_, COMP_STYLE_ or COMP_STACK_ suffix; does not change the .blend')
    p.add_argument('--exr',action='store_true',help='Scene-linear multilayer EXR for external compositing')
    p.add_argument('--resume',action='store_true',help='Reuse only matching scene/settings and fully validated PNG frames')
    args=p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    scene=bpy.context.scene
    if args.compositor:
        group=bpy.data.node_groups.get('COMP_'+args.compositor)
        if group is None:group=bpy.data.node_groups.get('COMP_STYLE_'+args.compositor)
        if group is None:group=bpy.data.node_groups.get('COMP_STACK_'+args.compositor)
        if group is None:raise SystemExit('Requested compositor branch is missing from this saved scene.')
        scene.compositing_node_group=group
        scene['selected_compositor']=args.compositor
    frames=list(range(scene.frame_start,scene.frame_end+1)) if args.frames=='all' else [int(f) for f in args.frames.split(',')]
    if any(f<scene.frame_start or f>scene.frame_end for f in frames):
        raise SystemExit('Requested frame outside saved scene range.')
    if args.resume and args.exr:raise SystemExit('Resume validation currently supports PNG; export EXR to a fresh directory.')
    if args.directory.exists() and not args.resume:raise SystemExit('Choose a fresh render directory or explicitly request --resume.')
    args.directory.mkdir(parents=True,exist_ok=True)
    if args.exr:
        if hasattr(scene.render.image_settings,'media_type'):
            scene.render.image_settings.media_type='MULTI_LAYER_IMAGE'
        scene.render.image_settings.file_format='OPEN_EXR_MULTILAYER'
        scene.render.image_settings.color_depth='32'
        scene.render.image_settings.exr_codec='ZIP'
        if hasattr(scene.render.image_settings,'use_exr_interleave'):
            scene.render.image_settings.use_exr_interleave=True
        scene.render.film_transparent=True
        scene.render.use_compositing=False
        for layer in scene.view_layers:
            layer.use_pass_z=True
            layer.use_pass_diffuse_color=True
            layer.use_pass_diffuse_direct=True
            layer.use_pass_diffuse_indirect=True
            layer.use_pass_glossy_direct=True
            layer.use_pass_glossy_indirect=True
            layer.use_pass_emit=True
            layer.use_pass_cryptomatte_object=True
            layer.pass_cryptomatte_depth=6
    if args.width:
        ratio=scene.render.resolution_y/scene.render.resolution_x
        scene.render.resolution_x=args.width
        scene.render.resolution_y=round(args.width*ratio)
    if scene.render.engine=='CYCLES':
        prefs=bpy.context.preferences.addons['cycles'].preferences
        try:
            prefs.compute_device_type='METAL'
            prefs.get_devices()
            metal=[d for d in prefs.devices if d.type=='METAL']
            for device in prefs.devices:device.use=device in metal
            if metal:scene.cycles.device='GPU'
        except (TypeError,RuntimeError):pass
    signature={'scene_sha256':hashlib.sha256(Path(bpy.data.filepath).read_bytes()).hexdigest(),
               'render_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'blender':bpy.app.version_string,'engine':scene.render.engine,
               'resolution':[scene.render.resolution_x,scene.render.resolution_y,scene.render.resolution_percentage],
               'fps':scene.render.fps,'compositor':scene.get('selected_compositor'),'exr':args.exr}
    receipt=args.directory/'render-signature.json'
    if receipt.exists():
        if json.loads(receipt.read_text())!=signature:raise SystemExit('Cached scene/settings differ; use a fresh render directory.')
    elif any(args.directory.iterdir()):
        raise SystemExit('Existing files have no render signature; use a fresh render directory.')
    receipt.write_text(json.dumps(signature,indent=2))
    started=time.monotonic();skipped=[];rendered=[]
    for frame in frames:
        scene.frame_set(frame)
        extension='exr' if args.exr else 'png'
        target=args.directory/f'frame_{frame:04}.{extension}'
        width=scene.render.resolution_x*scene.render.resolution_percentage//100
        height=scene.render.resolution_y*scene.render.resolution_percentage//100
        if args.resume and valid_png(target,width,height):
            skipped.append(frame);continue
        scene.render.filepath=str(target)
        bpy.ops.render.render(write_still=True)
        rendered.append(frame)
    (args.directory/'render-report.json').write_text(json.dumps({
        'source_blend':bpy.data.filepath,'frames':frames,'fps':scene.render.fps,
        'resolution':[scene.render.resolution_x,scene.render.resolution_y],
        'seconds':time.monotonic()-started,'blender':bpy.app.version_string,
        'rendered_this_run':rendered,'validated_cached_frames':skipped,
        'linear_multilayer_exr':args.exr,'external_compositor_executed':False,
        'compositor':scene.get('selected_compositor'),
        'exr_storage':'legacy interleaved for external compatibility' if args.exr else None,
        'preview_display_transform':scene.view_settings.view_transform,
        'preview_look':scene.view_settings.look,'preview_exposure':scene.view_settings.exposure,
        'handoff_note':'EXR is scene-linear, with premultiplied alpha. Match the saved Blender OCIO configuration and apply the display transform once.' if args.exr else None},indent=2))


if __name__=='__main__':main()
