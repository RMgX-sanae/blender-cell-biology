"""Assemble existing Blender camera scenes into an editable, source-timed teaching sequence."""
import json
import bpy


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
