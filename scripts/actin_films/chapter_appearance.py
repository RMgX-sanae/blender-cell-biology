"""Reference appearance on unchanged chapter geometry and authored events."""
import hashlib,json
from pathlib import Path
import bpy,numpy as np
from mathutils import Vector
from bpy_extras.object_utils import world_to_camera_view
import detail_scene as ds
import interaction_environment as environment
from interaction_appearance import build_compositor
from compositor_stack import publish_composite_outputs


def studio_profile(original,appearance):
    result=dict(original)
    for key in ('background','world','lighting','view_transform','look','exposure'):
        if key in appearance:result[key]=appearance[key]
    result['background']=(appearance.get('environment') or {}).get('compositor',{}).get('background',result['background'])
    if 'lighting' not in appearance:
        result['lighting']=[{**row,'position':[v*1.4 for v in row['position']],'power':row['power']*1.4**2,'size':row['size']*1.4} for row in result['lighting']]
    result['camera']={}
    return result


def finish(shot,profile,options,studio):
    scene=bpy.context.scene;points=[]
    bodies=[o for group in ('Molecules','Schematics') for o in shot.groups[group].objects if o.type in ('MESH','CURVE')]
    for frame in range(1,shot.n+1):
        scene.frame_set(frame);bpy.context.view_layer.update();graph=bpy.context.evaluated_depsgraph_get()
        for obj in bodies:
            if obj.hide_render:continue
            evaluated=obj.evaluated_get(graph)
            points.extend([evaluated.matrix_world@Vector(p) for p in evaluated.bound_box])
    scene.frame_set(1);bpy.context.view_layer.update()
    env=None;groups=[]
    if profile.get('environment'):
        cloud=np.array(points)
        if not len(cloud):raise ValueError('Source-shaped context requires a modeled foreground')
        groups,env=environment.build(options['environment_dir'],profile['environment'],profile,(cloud.min(0),cloud.max(0)))
    # Diagram text keeps its intentional dark-on-light ink. Convert any legacy
    # colored caption to unlit ink before moving it out of the beauty lighting.
    converted={};fonts=[]
    for obj in shot.labels:
        original=obj.data.materials[0]
        if original.name not in converted:
            emits=[n for n in original.node_tree.nodes if n.type=='EMISSION']
            if emits:converted[original.name]=original
            else:
                color=original.node_tree.nodes['Principled BSDF'].inputs['Base Color'].default_value[:]
                mat=bpy.data.materials.new('Protected caption '+original.name);mat.use_nodes=True;mat.node_tree.nodes.clear()
                emission=mat.node_tree.nodes.new('ShaderNodeEmission');emission.inputs['Color'].default_value=color
                output=mat.node_tree.nodes.new('ShaderNodeOutputMaterial');mat.node_tree.links.new(emission.outputs[0],output.inputs['Surface'])
                converted[original.name]=mat
        obj.data.materials[0]=converted[original.name];obj['compositor_preserve_text_color']=True
        emit=next(n for n in obj.data.materials[0].node_tree.nodes if n.type=='EMISSION')
        kind='foreground' if original==shot.ink else 'diagram' if original==shot.diagram_ink else 'palette'
        fonts.append({'object':obj.name,'body':obj.data.body,'ink_role':kind,'palette_key':original.get('palette_key'),
                      'color':list(emit.inputs['Color'].default_value),'strength':emit.inputs['Strength'].default_value})
    compositor=build_compositor(scene,profile,studio['background'],shot.width)
    if groups:environment.compose(scene,groups,profile['environment'],studio['background'],shot.width)
    fade=[]
    if shot.c['id']=='a03':
        tree=scene.compositing_node_group;out=next(n for n in tree.nodes if n.type=='GROUP_OUTPUT' and n.is_active_output)
        final=out.inputs['Image'].links[0].from_node;sources={}
        for name,port in [('beauty','Background'),('labels','Foreground')]:
            source=final.inputs[port].links[0].from_socket;over=tree.nodes.new('CompositorNodeAlphaOver');over.label='Authored chapter fade '+name
            over.inputs['Background'].default_value=(*ds.linear(studio['background']),1) if name=='beauty' else (0,0,0,0)
            tree.links.new(source,over.inputs['Foreground']);tree.links.new(over.outputs['Image'],final.inputs[port]);sources[name]=over.outputs['Image']
            for t,value in [(0,1),(.285,1),(.300,0),(.307,0),(.324,1),(1,1)]:
                over.inputs['Factor'].default_value=value;over.inputs['Factor'].keyframe_insert('default_value',frame=shot.frame(t))
            fade.append(over.label)
        publish_composite_outputs(tree,sources['beauty'],sources['labels'])
    scene.frame_set(1);bpy.context.view_layer.update()
    record={'profile':profile,'profile_path':options['style_profile'],'profile_sha256':options['style_profile_sha256'],
            'environment_sha256':options.get('environment_sha256'),'environment':env,'studio':studio,
            'palette':shot.palette,'fonts':fonts,'font_zh_sha256':hashlib.sha256(Path(options['font_zh']).read_bytes()).hexdigest(),
            'font_en_sha256':hashlib.sha256(Path(options['font_en']).read_bytes()).hexdigest(),
            'camera':{'matrix':np.array(scene.camera.matrix_world).tolist(),'scale':scene.camera.data.ortho_scale},
            'base_scene_scale_nm':shot.original_scene_scale,
            'compositor':compositor,'editorial_fade_nodes':fade,
            'scope':'Appearance of the existing authored frontal teaching layout; unchanged source geometry/events, not additional biological or optical evidence.'}
    scene['actin_appearance']=json.dumps(record)
    return record


def label_state(scene):
    """Evaluate animated labels in their own layer, including excluded parents."""
    original=bpy.context.window.view_layer;states={}
    try:
        bpy.context.window.view_layer=scene.view_layers[scene['compositor_label_layer']]
        bpy.context.view_layer.update();graph=bpy.context.evaluated_depsgraph_get()
        for obj in scene.objects:
            if obj.type!='FONT':continue
            evaluated=obj.evaluated_get(graph)
            states[obj.name]={'matrix':np.array(evaluated.matrix_world),'hidden':obj.hide_render,
                              'bounds':[world_to_camera_view(scene,scene.camera,evaluated.matrix_world@Vector(p))[:] for p in evaluated.bound_box]}
    finally:
        bpy.context.window.view_layer=original;bpy.context.view_layer.update()
    return states
