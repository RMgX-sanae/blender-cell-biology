"""Check effective V-ATPase appearance and protect the authored dual-view layout."""
import hashlib,json,math,sys
from pathlib import Path
import bpy,numpy as np
from mathutils import Vector
from bpy_extras.object_utils import world_to_camera_view
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_surface_materials import check_materials
from interaction_environment import audit_environment
from export_composite_layers import published_sources
from detail_scene import linear
from vatpase_appearance import digest,read_profile,MEMBRANE


def text_bounds(scene):
    original=bpy.context.window.view_layer
    try:
        bpy.context.window.view_layer=scene.view_layers[scene['compositor_label_layer']]
        bpy.context.view_layer.update();graph=bpy.context.evaluated_depsgraph_get()
        return {o.name:[world_to_camera_view(scene,scene.camera,o.evaluated_get(graph).matrix_world@Vector(p))[:] for p in o.evaluated_get(graph).bound_box]
                for o in scene.objects if o.type=='FONT'}
    finally:bpy.context.window.view_layer=original;bpy.context.view_layer.update()


def audit(scene):
    if 'vatpase_appearance' not in scene:return None
    record=json.loads(scene['vatpase_appearance']);folder=Path(bpy.data.filepath).parent
    manifest=json.loads((folder/'manifest.json').read_text())
    if digest(record['profile_path'])!=record['profile_sha256'] or read_profile(record['profile_path'])!=record['profile'] or manifest['appearance']!=record:
        raise ValueError('Saved V-ATPase appearance differs from its immutable request')
    p=record['profile'];palette=record['palette'];studio=record['studio']
    if manifest['labels']!=[row['body'] for row in record['fonts']]:raise ValueError('Manifest lost the protected caption inventory')
    def equal(a,b,message):
        if np.shape(a)!=np.shape(b) or not np.allclose(a,b,atol=2e-6,rtol=1e-5):raise ValueError(message)
    if any(palette[k]!=v for k,v in p.get('palette',{}).items()):raise ValueError('Requested semantic palette ignored')
    info=json.loads((Path(manifest['asset_dir'])/'manifest.json').read_text());roles=set();bodies=[]
    for group in ('Pump','C ring detail','Membrane','Flow'):
        for obj in bpy.data.collections[group].objects:
            if obj.type not in ('MESH','CURVE'):continue
            if obj.get('source')=='6VQ6':
                c=obj['source_chain'];key='A' if c in 'ABC' else 'B' if c in 'DEF' else 'a' if c=='a' else 'ring_tracer' if c=='b' else 'ring' if c in info['c_ring_chains'] else 'rotor' if c in info['rotor_chains'] else 'stalk'
                expected={('protein',key)}
            elif group=='Membrane':expected={('membrane_heads','lipid'),('membrane_tails','tail')}
            elif obj.name.startswith(('H transport token ','Left H ')):expected={('proton_marker','proton')}
            elif obj.name=='Adenosine and two phosphates':expected={('chemical_token','atp'),('chemical_token','h')}
            elif obj.name in ('Terminal phosphate becomes Pi','ATP terminal bond'):expected={('chemical_token','h')}
            else:
                key='in' if obj.name.startswith(('Cytoplasmic','in flow','Whole-pump cytoplasmic')) else 'out' if obj.name.startswith(('Luminal','out flow','Whole-pump luminal')) else 'background' if obj.name=='ATP approach guide' else 'stalk'
                expected={('flow_guide',key)}
            if obj.get('source')=='6VQ6':
                modifiers=list(obj.modifiers)
                if len(modifiers)!=1 or modifiers[0].type!='SMOOTH' or not modifiers[0].show_render or not modifiers[0].show_viewport or modifiers[0].iterations!=2:
                    raise ValueError('Source envelope modifier differs')
                equal(modifiers[0].factor,.25,'Source envelope smoothing differs')
            elif obj.modifiers:raise ValueError('Schematic geometry modifier changed')
            if obj.type=='MESH':
                evaluated=obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
                mats={obj.material_slots[p.material_index].material if p.material_index<len(obj.material_slots) else None for p in evaluated.data.polygons}
            else:mats=set(obj.data.materials)
            if not mats or None in mats or {(m.get('surface_role'),m.get('palette_key')) for m in mats}!=expected:raise ValueError('V-ATPase object material identity differs: '+obj.name)
            for mat in mats:
                shader=mat.node_tree.nodes['Principled BSDF'];socket=shader.inputs['Base Color']
                if socket.is_linked:raise ValueError('Semantic color overridden')
                equal(socket.default_value,(*linear(palette[mat['palette_key']]),1),'Semantic color differs')
                if mat['surface_role']=='proton_marker':
                    equal(shader.inputs['Emission Color'].default_value,socket.default_value,'Proton color differs')
                    equal(shader.inputs['Emission Strength'].default_value,.28,'Proton glow differs')
            roles.update(r for r,_ in expected);bodies.append(obj.name)
    if record['environment']:roles.add('environment')
    if roles&p.get('surface_roles',{}).keys():check_materials(p,roles&p['surface_roles'].keys())
    camera=scene.camera;controls=p.get('camera',{});pan=controls.get('center_offset_nm',[0,0])
    if camera.data.type!='ORTHO' or camera.data.dof.use_dof or camera.parent or camera.constraints or camera.animation_data or camera.data.animation_data:raise ValueError('Authored dual-view camera changed')
    equal(np.array(camera.matrix_world),record['camera']['matrix'],'Camera transform differs')
    equal(camera.location,[pan[0],-100,8.5+pan[1]],'Requested pan ignored')
    equal(camera.data.ortho_scale,58*controls.get('scale_multiplier',1),'Requested scale ignored')
    lights={o.name:o for o in scene.objects if o.type=='LIGHT'}
    if set(lights)!={r['name'] for r in studio['lighting']}:raise ValueError('Light inventory differs')
    for row in studio['lighting']:
        obj=lights[row['name']]
        if obj.data.type!='AREA' or obj.hide_render or obj.parent or obj.constraints or obj.animation_data or obj.data.animation_data:raise ValueError('Light display differs')
        equal(obj.location,row['position'],'Light position differs');equal(obj.data.energy,row['power'],'Light power differs')
        equal(obj.data.size,row['size'],'Light size differs');equal(obj.data.color,linear(row['color']),'Light color differs')
        direction=np.array([pan[0],0,8.5+pan[1]])-np.array(obj.location);direction/=np.linalg.norm(direction)
        equal(-np.array(obj.matrix_world.to_3x3())[:,2],direction,'Light aim differs')
    for key in ('world','lighting','view_transform','look','exposure'):
        if key in p and p[key]!=studio[key]:raise ValueError('Requested studio setting ignored')
    world=scene.world.node_tree.nodes['Background'];equal(world.inputs['Color'].default_value,(*linear(studio['world']),1),'World color differs');equal(world.inputs['Strength'].default_value,.35,'World strength differs')
    if scene.view_settings.view_transform!=studio.get('view_transform','AgX') or scene.view_settings.look!=studio.get('look','AgX - Medium High Contrast' if scene.view_settings.view_transform=='AgX' else 'None'):raise ValueError('Display transform differs')
    equal(scene.view_settings.exposure,studio.get('exposure',0),'Exposure differs')
    tree=scene.compositing_node_group;ports=published_sources(tree)
    def node(label):
        found=[n for n in tree.nodes if n.label==label]
        if len(found)!=1:raise ValueError('Missing compositor node: '+label)
        return found[0]
    def linked(socket,expected,message):
        if len(socket.links)!=1 or socket.links[0].from_socket!=expected:raise ValueError(message)
    base=node('Straight background plus premultiplied rendered foreground');beauty=scene.view_layers[scene['compositor_beauty_layer']];text=scene.view_layers[scene['compositor_label_layer']]
    if any(name not in beauty.objects for name in bodies):raise ValueError('Foreground excluded from beauty')
    if not beauty.layer_collection.children[scene['compositor_label_collection']].exclude:raise ValueError('Text leaks into beauty')
    if {c.name for c in text.layer_collection.children if not c.exclude}!={scene['compositor_label_collection']}:raise ValueError('Foreground leaks into text layer')
    if record['membrane']!=MEMBRANE:raise ValueError('Cytosolic boundary differs')
    if record['environment']:
        audit_environment(scene,record['environment'],record,node,linked,equal,base,studio['background'])
        for row in record['environment']['objects']:
            obj=scene.objects[row['object']]
            if max((obj.matrix_world@v.co).x for v in obj.data.vertices)>2:raise ValueError('Context crosses into magnified inset')
    else:equal(base.inputs['Background'].default_value,(*linear(studio['background']),1),'Background differs')
    if scene.get('selected_compositor')!=p.get('render_recipe','neutral'):raise ValueError('Render recipe differs')
    render=ports['Labels'].node
    if render.type!='R_LAYERS' or render.layer!=text.name or ports['Labels']!=render.outputs['Image']:raise ValueError('Text routed through filters')
    if scene['selected_compositor']=='neutral' and ports['Beauty']!=base.outputs['Image']:raise ValueError('Neutral beauty bypassed')
    pending=[ports['Beauty'].node];reached=set()
    while pending:
        current=pending.pop()
        if current in reached:continue
        reached.add(current);pending.extend(link.from_node for socket in current.inputs for link in socket.links)
    if base not in reached:raise ValueError('Foreground composite does not reach output')
    fonts=[o for o in scene.objects if o.type=='FONT']
    if {o.name:o.data.body for o in fonts}!={r['object']:r['body'] for r in record['fonts']}:raise ValueError('Caption identity differs')
    for obj in fonts:
        if not obj.get('compositor_preserve_text_color'):raise ValueError('Text protection removed')
        for font,sha in zip([obj.data.font,obj.data.font_bold],record['font_sha256']):
            if not font.packed_file or hashlib.sha256(bytes(font.packed_file.data)).hexdigest()!=sha:raise ValueError('Packed bilingual font differs')
        if any(fmt.use_bold!=(ord(c)<0x400 or 0x2010<=ord(c)<=0x22ff) for c,fmt in zip(obj.data.body,obj.data.body_format)):raise ValueError('Mixed font routing differs')
        mat=obj.data.materials[0];out=next(n for n in mat.node_tree.nodes if n.type=='OUTPUT_MATERIAL' and n.is_active_output)
        if len(out.inputs['Surface'].links)!=1 or out.inputs['Surface'].links[0].from_node.type!='EMISSION':raise ValueError('Protected text is not unlit')
        emit=out.inputs['Surface'].links[0].from_node;equal(emit.inputs['Color'].default_value,(*linear(palette['ink']),1),'Text ink differs');equal(emit.inputs['Strength'].default_value,1,'Text strength differs')
    try:
        for f in range(1,386):
            scene.frame_set(f)
            for i,c in enumerate('ABC'):
                shader=scene.objects['6VQ6 chain '+c].data.materials[0].node_tree.nodes['Principled BSDF'];u=((f-1)/384*3-i)%3
                equal(shader.inputs['Emission Strength'].default_value,.13*max(0,math.sin(math.pi*u))**2 if u<1 else 0,'ATP highlight timing differs')
                equal(shader.inputs['Emission Color'].default_value,(*linear(palette['A']),1),'ATP highlight color differs')
    finally:scene.frame_set(1)
    return {'passed':True,'roles':sorted(roles),'objects':len(bodies),'lights':len(lights),'labels':len(fonts),'highlight_frames':385,
            'context_instances':len(record['environment']['objects']) if record['environment'] else 0,
            'scope':'Actual appearance, source context and protected captions; geometry/flow synchronization and pixels are separate checks.'}
