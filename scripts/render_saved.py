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


def selected_frames(spec, first, last):
    if spec=='all':return list(range(first,last+1))
    selected=set()
    for part in spec.split(','):
        values=[int(x) for x in part.split(':')]
        if len(values)==1:selected.add(values[0])
        elif len(values) in [2,3]:
            start,end=values[:2];step=values[2] if len(values)==3 else 1
            if step<=0 or start>end:raise ValueError('Require increasing ranges and a positive step')
            if start<first or end>last:raise ValueError('Range outside saved scene')
            selected.update(range(start,end+1,step))
        else:raise ValueError('Use frames or inclusive START:END[:STEP] ranges')
    if not selected or min(selected)<first or max(selected)>last:raise ValueError('Frame outside saved scene')
    return sorted(selected)


def brief_frames(path, scene):
    raw = path.read_bytes()
    brief = json.loads(raw)
    rows = brief.get('review_frames') if isinstance(brief, dict) else None
    if not isinstance(rows, list) or not rows:
        raise ValueError('Give the project brief a nonempty review_frames list')
    questions = {}
    for row in rows:
        if not isinstance(row, dict) or not isinstance(row.get('look_for'), str) or not row['look_for'].strip():
            raise ValueError('Each review frame needs a look_for question')
        frame = row.get('frame')
        if frame == 'saved':
            if scene.frame_subframe:
                raise ValueError('Saved pose is at a subframe; specify an integer review frame')
            frame = scene.frame_current
        if type(frame) is not int or not scene.frame_start <= frame <= scene.frame_end:
            raise ValueError('Review frame must be saved or an integer inside the scene range')
        questions.setdefault(frame, []).append(row['look_for'].strip())
    frames = sorted(questions)
    return frames, {'sha256': hashlib.sha256(raw).hexdigest(),
                    'frames': [{'frame': f, 'look_for': questions[f]} for f in frames]}


def saved_compositor(scene, suffix):
    group = scene.get('detail_compositor_branches', {}).get(suffix)
    if group is not None:
        return group
    for prefix in ('COMP_', 'COMP_STYLE_', 'COMP_STACK_'):
        group = bpy.data.node_groups.get(prefix + suffix)
        if group is not None:
            return group
    return None


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory',type=Path,required=True)
    choice=p.add_mutually_exclusive_group()
    choice.add_argument('--frames',default='all',help='all, integer frames or inclusive START:END[:STEP] ranges, comma-separated')
    choice.add_argument('--review-brief',type=Path,help='Render only the poses and visual questions in a project brief')
    p.add_argument('--width',type=int)
    p.add_argument('--samples',type=int,help='Temporary Eevee/Cycles sample count for review; does not resave the scene')
    p.add_argument('--compositor',help='Exact saved COMP_, COMP_STYLE_ or COMP_STACK_ suffix; does not change the .blend')
    p.add_argument('--exr',action='store_true',help='Scene-linear multilayer EXR for external compositing')
    p.add_argument('--resume',action='store_true',help='Reuse only matching scene/settings and fully validated PNG frames')
    args=p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    scene=bpy.context.scene
    locator_script = None
    if args.compositor:
        group=saved_compositor(scene,args.compositor)
        if group is None:raise SystemExit('Requested compositor branch is missing from this saved scene.')
        scene.compositing_node_group=group
        scene['selected_compositor']=args.compositor
        if scene.get('linked_view_locator'):
            import build_linked_views
            locator = json.loads(scene['linked_view_locator'])
            build_linked_views.add_locator(scene, locator['at'], locator['radius'])
            locator_script = Path(build_linked_views.__file__)
    review = None
    try:
        if args.review_brief:
            frames, review = brief_frames(args.review_brief, scene)
        else:
            frames = selected_frames(args.frames, scene.frame_start, scene.frame_end)
    except (ValueError, OSError) as error:
        p.error(str(error))
    if args.width is not None and args.width<=0:p.error('Width must be positive')
    if args.samples is not None:
        if args.samples<=0:p.error('Samples must be positive')
        if scene.render.engine=='CYCLES':scene.cycles.samples=args.samples
        elif scene.render.engine=='BLENDER_EEVEE':scene.eevee.taa_render_samples=args.samples
        else:p.error('Sample override supports Eevee and Cycles')
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
               'blender':bpy.app.version_string,'engine':scene.render.engine,'scene_name':scene.name,
               'resolution':[scene.render.resolution_x,scene.render.resolution_y,scene.render.resolution_percentage],
               'fps':scene.render.fps,'fps_base':scene.render.fps_base,
               'samples':scene.cycles.samples if scene.render.engine=='CYCLES' else scene.eevee.taa_render_samples if scene.render.engine=='BLENDER_EEVEE' else None,
               'compositor':scene.get('selected_compositor'),'exr':args.exr}
    if locator_script:
        signature['locator_script_sha256'] = hashlib.sha256(locator_script.read_bytes()).hexdigest()
    receipt=args.directory/'render-signature.json'
    if receipt.exists():
        if json.loads(receipt.read_text())!=signature:raise SystemExit('Cached scene/settings differ; use a fresh render directory.')
    elif any(args.directory.iterdir()):
        raise SystemExit('Existing files have no render signature; use a fresh render directory.')
    receipt.write_text(json.dumps(signature,indent=2))
    started=time.monotonic();skipped=[];rendered=[]
    for frame in frames:
        scene.frame_set(frame, subframe=scene.get('linked_view_subframe', 0.0))
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
        'fps_base':scene.render.fps_base,'samples':signature['samples'],
        'scene_frame_range':[scene.frame_start,scene.frame_end],
        'source_seconds':[(frame-scene.frame_start)*scene.render.fps_base/scene.render.fps for frame in frames],
        'resolution':[scene.render.resolution_x,scene.render.resolution_y],
        'seconds':time.monotonic()-started,'blender':bpy.app.version_string,
        'rendered_this_run':rendered,'validated_cached_frames':skipped,
        **({'review_brief':review} if review is not None else {}),
        'linear_multilayer_exr':args.exr,'external_compositor_executed':False,
        'compositor':scene.get('selected_compositor'),
        'exr_storage':'legacy interleaved for external compatibility' if args.exr else None,
        'preview_display_transform':scene.view_settings.view_transform,
        'preview_look':scene.view_settings.look,'preview_exposure':scene.view_settings.exposure,
        'handoff_note':'EXR is scene-linear, with premultiplied alpha. Match the saved Blender OCIO configuration and apply the display transform once.' if args.exr else None},indent=2))


if __name__=='__main__':main()
