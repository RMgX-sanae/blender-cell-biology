"""Assemble existing Blender camera scenes into an editable, source-timed teaching sequence."""
import json
import bpy

import detail_scene as ds
from build_translation_scene import emission, text_object


def compose(name, chapters, reference):
    """Each chapter is (title, source scene, inclusive first frame, inclusive last frame)."""
    if not chapters or name in bpy.data.scenes:
        raise ValueError('Choose chapters and a fresh master-scene name')
    fps = (reference.render.fps,reference.render.fps_base)
    for title,source,start,end in chapters:
        if not title or not source.camera or not source.frame_start<=start<=end<=source.frame_end:
            raise ValueError('Select an in-range, nonempty camera shot')
        if (source.render.fps,source.render.fps_base)!=fps:
            raise ValueError('Match scene playback rates before joining shots')
    master=bpy.data.scenes.new(name)
    master.render.engine=reference.render.engine
    for field in ('resolution_x','resolution_y','resolution_percentage','pixel_aspect_x','pixel_aspect_y','fps','fps_base'):
        setattr(master.render,field,getattr(reference.render,field))
    master.display_settings.display_device=reference.display_settings.display_device
    for field in ('view_transform','look','exposure','gamma'):
        setattr(master.view_settings,field,getattr(reference.view_settings,field))
    master.sequencer_colorspace_settings.name=reference.sequencer_colorspace_settings.name
    master.render.use_sequencer=True
    strips=master.sequence_editor_create().strips
    cursor=1
    record=[]
    for title,source,start,end in chapters:
        left=start-source.frame_start
        strip=strips.new_scene(title,source,1,cursor-left)
        strip.scene_input='CAMERA'
        strip.frame_offset_start=left
        strip.frame_offset_end=source.frame_end-end
        stop=cursor+end-start
        if (strip.frame_final_start,strip.frame_final_end)!=(cursor,stop+1):
            raise ValueError('The scene strip did not retain its requested time range')
        master.timeline_markers.new(title,frame=cursor)
        for marker in source.timeline_markers:
            if start<marker.frame<=end:
                master.timeline_markers.new(marker.name,frame=cursor+marker.frame-start)
        record.append({'title':title,'scene':source.name,'source':[start,end],'master':[cursor,stop]})
        cursor=stop+1
    master.frame_start,master.frame_end=1,cursor-1
    master['native_chapters']=json.dumps(record,ensure_ascii=False)
    bpy.context.window.scene=master
    for area in bpy.context.screen.areas:
        if area.type=='VIEW_3D':
            area.type='SEQUENCE_EDITOR'
            area.spaces.active.view_type='SEQUENCER_PREVIEW'
            break
    master.frame_set(1)
    return master


def visible(obj, start, end):
    for frame,flag in ((start-1,True),(start,False),(end+1,True)):
        obj.hide_viewport=obj.hide_render=flag
        obj.keyframe_insert(data_path='hide_viewport',frame=frame)
        obj.keyframe_insert(data_path='hide_render',frame=frame)


def empty(name,group):
    obj=bpy.data.objects.new(name,None)
    obj.empty_display_size=.3
    group.objects.link(obj)
    return obj


def balls(name,positions,radius,mat,group,parent=None):
    parts=ds.Parts()
    for point in positions:
        parts.ball(point,radius,seg=12,rings=8)
    obj=parts.object(name,[mat],group); obj.parent=parent
    return obj


def new_shot(name,frames,view,reference,profile):
    from types import SimpleNamespace
    shot=bpy.data.scenes.new(name)
    bpy.context.window.scene=shot
    args=SimpleNamespace(width=reference.render.resolution_x,frames=frames,fps=reference.render.fps,
                         samples=reference.cycles.samples if reference.render.engine=='CYCLES' else reference.eevee.taa_render_samples,
                         engine='cycles' if reference.render.engine=='CYCLES' else 'eevee',compositor='clean')
    studio=dict(profile,background=profile.get('background','091A24'),camera={},
                compositor=profile.get('compositor',{'preset':'clean','contrast':2}))
    cam=ds.studio(studio,args,view)
    shot.render.resolution_y=reference.render.resolution_y
    shot.render.fps_base=reference.render.fps_base
    return shot,cam


def caption_set(scene,camera,fonts,lang,stages,scales,*,outline=False):
    group=ds.collection('Process labels · '+scene.name)
    overlay=empty('Camera typography · '+scene.name,group)
    overlay.parent=camera; overlay.location=(0,0,-1)
    for frame,scale in scales:
        overlay.scale=(scale,)*3;overlay.keyframe_insert(data_path='scale',frame=frame)
    ink=emission('Caption ink · '+scene.name,'EDF3F0')
    edge=emission('Caption outline · '+scene.name,'102028') if outline else None
    zh,latin=fonts
    for index,(start,cn,en) in enumerate(stages):
        body=cn if lang=='zh' else en
        text=text_object(body,'Phase · '+body,.0275,ink,group,zh if lang=='zh' else latin)
        text.data.font_bold=latin
        for char,formatting in zip(body,text.data.body_format):
            formatting.use_bold=(lang=='zh' and char.isascii())
        text.parent=overlay; text.location=(-.452,.255,0)
        for name in ('visible_shadow','visible_diffuse','visible_glossy','visible_transmission','visible_volume_scatter'):
            setattr(text,name,False)
        end=stages[index+1][0]-1 if index+1<len(stages) else scene.frame_end
        visible(text,start,end)
        if edge:
            border=text.copy();border.data=text.data.copy()
            border.name='Contrast · '+body;border.data.offset=.0011
            border.data.materials.clear();border.data.materials.append(edge)
            group.objects.link(border)
            border.location.z-=.00025
        scene.timeline_markers.new(body,frame=start)


def font_pair(args,lang):
    def load_font(given,names):
        if given:
            font=bpy.data.fonts.load(str(given.resolve()),check_existing=True)
            if not font.packed_file:font.pack()
            return font
        return next((f for f in bpy.data.fonts if f.packed_file and any(n in f.name.lower() for n in names)),None)
    zh=load_font(args.font_zh,('songti','simsun'))
    en=load_font(args.font_latin,('times new roman','timesnewroman'))
    if en is None or (lang=='zh' and zh is None):
        raise ValueError('Supply the missing Songti and/or Times New Roman font path for this project')
    return zh,en


def animated_path(name,frames,coordinates,radius,mat,group):
    obj=ds.path(name,coordinates(frames[0]),radius,mat,group)
    obj.shape_key_add(name='First pose')
    keys=obj.data.shape_keys;keys.use_relative=False
    for index,frame in enumerate(frames):
        key=keys.key_blocks[0] if index==0 else obj.shape_key_add(name='Frame '+str(frame))
        for p,co in zip(key.data,coordinates(frame)):p.co=co
        key.interpolation='KEY_LINEAR'
        keys.eval_time=key.frame;keys.keyframe_insert(data_path='eval_time',frame=frame)
    action=keys.animation_data.action
    for layer in action.layers:
        for strip in layer.strips:
            for bag in strip.channelbags:
                for curve in bag.fcurves:
                    for point in curve.keyframe_points:point.interpolation='LINEAR'
    return obj
