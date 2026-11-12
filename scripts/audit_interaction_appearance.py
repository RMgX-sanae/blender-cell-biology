"""Check effective saved appearance independently of prepared interaction poses."""
import hashlib,json,re
from pathlib import Path
import bpy,numpy as np
from mathutils import Vector
from bpy_extras.object_utils import world_to_camera_view
from audit_surface_materials import check_materials
from detail_scene import linear


def check_appearance(scene,manifest):
    if 'interaction_appearance' not in scene:return None
    from audit_interaction_scene import displayed_materials
    record=json.loads(scene['interaction_appearance']);profile=record['profile'];path=Path(record['path'])
    if hashlib.sha256(path.read_bytes()).hexdigest()!=record['sha256'] or json.loads(path.read_text())!=profile:raise ValueError('Appearance input changed')
    saved=json.loads((Path(bpy.data.filepath).parent/'manifest.json').read_text())
    if saved['appearance']!=record:raise ValueError('Appearance record differs from project manifest')
    if json.loads(scene['style_profile'])!=record['studio']:raise ValueError('Resolved studio differs')
    def equal(a,b,message):
        if np.shape(a)!=np.shape(b) or not np.allclose(a,b,atol=2e-6,rtol=1e-5):raise ValueError(message)
    required=set();colors=record['actor_colors']
    if colors!=json.loads(scene['interaction_resolved_colors']):raise ValueError('Actor identity color record differs')
    for name,color in profile.get('actor_colors',{}).items():
        if colors[name]!=color:raise ValueError('Requested actor color was ignored')
    for actor in manifest['actors']:
        role='cofactor' if actor['kind']!='polymer' else 'nucleic_acid' if actor.get('polymer_type','').startswith(('polyribonucleotide','polydeoxyribonucleotide')) else 'protein'
        if role not in profile.get('surface_roles',{}):continue
        assigned=displayed_materials(scene.objects[actor['id']])
        if len(assigned)!=1 or None in assigned or next(iter(assigned)).get('surface_role')!=role:raise ValueError('Actor uses the wrong surface role')
        mat=next(iter(assigned));socket=mat.node_tree.nodes['Principled BSDF'].inputs['Base Color']
        if socket.is_linked:raise ValueError('Actor identity color overridden')
        equal(socket.default_value,(*linear(colors[actor['id']]),1),'Actor identity color differs')
        required.add(role)
    if manifest.get('membrane'):
        defaults={'heads':'#638fa0','tails':'#365f72','tether':json.loads((Path(__file__).resolve().parents[1]/'assets/style-library.json').read_text())['profiles'][record['style']]['palette'][0]}
        defaults.update(profile.get('membrane_palette',{}))
        for name,key,role in [('Membrane heads','heads','membrane_heads'),('Membrane tails','tails','membrane_tails'),('Cargo tether','tether','schematic_tether')]:
            mats=displayed_materials(scene.objects[name])
            if len(mats)!=1 or None in mats or next(iter(mats)).get('surface_role')!=role:raise ValueError('Membrane surface role differs')
            mat=next(iter(mats));equal(mat.node_tree.nodes['Principled BSDF'].inputs['Base Color'].default_value,(*linear(defaults[key]),1),'Membrane identity color differs')
            if role in profile.get('surface_roles',{}):required.add(role)
    env=record['environment']
    if env and 'environment'in profile.get('surface_roles',{}):required.add('environment')
    materials=check_materials(profile,required) if required else None
    for obj in scene.objects:
        if obj.get('legend_for_objects'):
            identities=json.loads(obj['legend_for_objects']);expected={colors[i] for i in identities}
            if len(expected)!=1:raise ValueError('Legend merges different identity colors')
            mats=displayed_materials(obj)
            if len(mats)!=1 or None in mats:raise ValueError('Legend material disconnected')
            mat=next(iter(mats));rgb=np.array(linear(next(iter(expected))))
            ramps=[n for n in mat.node_tree.nodes if n.type=='VALTORGB'];emits=[n for n in mat.node_tree.nodes if n.type=='EMISSION']
            if len(ramps)!=1 or len(emits)!=1 or len(emits[0].inputs['Color'].links)!=1 or emits[0].inputs['Color'].links[0].from_node!=ramps[0]:raise ValueError('Legend color path changed')
            if len(ramps[0].color_ramp.elements)!=2:raise ValueError('Legend color ramp changed')
            for stop,factor in zip(ramps[0].color_ramp.elements,[.92,1.08]):equal(stop.color,(*(rgb*factor),1),'Legend color differs')
        if obj.type!='FONT':continue
        chinese=bool(re.search(r'[\u3400-\u9fff]',obj.data.body));expected=record['font_zh_sha256'] if chinese else record['font_en_sha256']
        if not obj.data.font.packed_file or hashlib.sha256(bytes(obj.data.font.packed_file.data)).hexdigest()!=expected:raise ValueError('Packed language font differs')
        if chinese:
            if not obj.data.font_bold.packed_file or hashlib.sha256(bytes(obj.data.font_bold.packed_file.data)).hexdigest()!=record['font_en_sha256']:raise ValueError('Mixed-language Times font differs')
            if any(form.use_bold!=(ord(char)<128 or char in 'αβγμσ') for char,form in zip(obj.data.body,obj.data.body_format)):raise ValueError('Mixed-language font routing differs')
        for corner in obj.bound_box:
            point=world_to_camera_view(scene,scene.camera,obj.matrix_world@Vector(corner))
            if not (0<=point.x<=1 and 0<=point.y<=1 and point.z>0):raise ValueError('Label leaves the image')
    for label in profile.get('actor_labels',{}).values():
        if not any(o.type=='FONT' and o.data.body==label for o in scene.objects):raise ValueError('Requested actor label missing')
    studio=record['studio'];lights={o.name:o for o in scene.objects if o.type=='LIGHT'}
    if set(lights)!={r['name'] for r in studio['lighting']}:raise ValueError('Light inventory differs')
    camera=scene.camera;expected=record['camera']
    if camera.animation_data or camera.constraints:raise ValueError('Appearance camera must remain static')
    equal(np.array(camera.matrix_world),expected['matrix'],'Saved camera transform differs')
    if camera.data.type!=expected['type'] or camera.data.dof.use_dof!=expected['dof']:raise ValueError('Camera projection or focus differs')
    for key in ('ortho_scale','lens','sensor_width'):equal(getattr(camera.data,key),expected[key],'Saved camera '+key+' differs')
    equal(camera.data.dof.focus_distance,expected['focus_distance'],'Camera focus distance differs')
    equal(camera.data.dof.aperture_fstop,expected['fstop'],'Camera aperture differs')
    requested=profile.get('camera',{})
    if 'target' in requested:
        direction=np.array(requested['target'])-requested['position'];distance=np.linalg.norm(direction)
        equal(-np.array(camera.matrix_world.to_3x3())[:,2],direction/distance,'Requested camera target ignored')
        equal(camera.data.dof.focus_distance,distance,'Requested camera target focus ignored')
    for key,actual in [('scale',camera.data.ortho_scale),('lens',camera.data.lens),('position',camera.location)]:
        if key in requested:equal(actual,requested[key],'Requested camera '+key+' ignored')
    if 'type'in requested and camera.data.type!=requested['type']:raise ValueError('Requested camera type ignored')
    if 'fstop'in requested and (camera.data.dof.use_dof!=bool(requested['fstop']) or requested['fstop'] and not np.isclose(camera.data.dof.aperture_fstop,requested['fstop'])):raise ValueError('Requested camera focus ignored')
    for row in studio['lighting']:
        obj=lights[row['name']]
        if obj.data.type!='AREA' or obj.hide_render or obj.animation_data or obj.constraints:raise ValueError('Light type or visibility differs')
        equal(obj.location,row['position'],'Light position differs');equal(obj.data.energy,row['power'],'Light power differs')
        equal(obj.data.size,row['size'],'Light size differs');equal(obj.data.color,linear(row['color']),'Light color differs')
        target=np.array(camera.location)-np.array(camera.matrix_world.to_3x3())[:,2]*camera.data.dof.focus_distance
        direction=target-np.array(obj.location);direction/=np.linalg.norm(direction)
        equal(-np.array(obj.matrix_world.to_3x3())[:,2],direction,'Light aim differs')
    for key in ('lighting','world','view_transform','look','exposure'):
        if key in profile and studio[key]!=profile[key]:raise ValueError('Requested studio setting ignored: '+key)
    world=scene.world.node_tree.nodes['Background'];equal(world.inputs['Color'].default_value,(*linear(studio['world']),1),'World color differs');equal(world.inputs['Strength'].default_value,.35,'World strength differs')
    if scene.view_settings.view_transform!=studio.get('view_transform','AgX') or scene.view_settings.look!=studio.get('look','AgX - Medium High Contrast' if scene.view_settings.view_transform=='AgX' else 'None'):raise ValueError('View transform or look differs')
    equal(scene.view_settings.exposure,studio.get('exposure',0),'Exposure differs')
    if scene['selected_compositor']!=profile.get('render_recipe','neutral'):raise ValueError('Render recipe differs')
    tree=scene.compositing_node_group
    def node(label):
        rows=[n for n in tree.nodes if n.label==label]
        if len(rows)!=1:raise ValueError('Missing compositor node: '+label)
        return rows[0]
    def linked(socket,expected,message):
        if len(socket.links)!=1 or socket.links[0].from_socket!=expected:raise ValueError(message)
    base=node('Straight background plus premultiplied rendered foreground')
    beauty=scene.view_layers[scene['compositor_beauty_layer']]
    if any(a['id'] not in beauty.objects for a in manifest['actors']):raise ValueError('Actor excluded from sharp beauty layer')
    links=base.inputs['Foreground'].links
    if len(links)!=1 or links[0].from_node.type!='R_LAYERS' or links[0].from_node.layer!=beauty.name:raise ValueError('Sharp subject render layer differs')
    if env:
        from interaction_environment import audit_environment
        audit_environment(scene,env,manifest,node,linked,equal,base,studio['background'])
    else:equal(base.inputs['Background'].default_value,(*linear(studio['background']),1),'Flat background color differs')
    output=next(n for n in tree.nodes if n.type=='GROUP_OUTPUT' and n.is_active_output)
    if len(output.inputs['Image'].links)!=1:raise ValueError('Final image disconnected')
    final=output.inputs['Image'].links[0].from_node
    if final.type!='ALPHAOVER' or len(final.inputs['Foreground'].links)!=1:raise ValueError('Protected labels no longer last')
    labels=final.inputs['Foreground'].links[0].from_node
    if labels.type!='R_LAYERS' or labels.layer!=scene['compositor_label_layer']:raise ValueError('Labels routed through image filters')
    if scene['selected_compositor']=='neutral':linked(final.inputs['Background'],base.outputs['Image'],'Neutral beauty output bypassed')
    return {'passed':True,'role_materials':len(materials['materials']) if materials else 0,'lights':len(lights),
            'background_instances':len(env['objects']) if env else 0,'depth_layers':3 if env else 0,
            'scope':'Actual role/color/light/camera/font handoff and protected depth routing. Source geometry and poses are checked separately; inspect final pixels and mechanism suitability.'}
