"""Inspect effective receptor materials, source background and protected depth composition."""
import hashlib,json
from pathlib import Path
import bpy,numpy as np
from audit_surface_materials import check_materials
from audit_interaction_scene import displayed_materials
from build_receptor_distribution_scene import appearance_profile,digest
from detail_scene import linear


def check_appearance(scene,record):
    if record.get('appearance') is None:return None
    source_path=Path(record['build_parameters']['style_profile'])
    if digest(source_path)!=record['appearance_sha256'] or json.loads(source_path.read_text())!=record['appearance']:
        raise ValueError('Selected appearance source changed')
    profile,settings=appearance_profile(record['build_parameters']['style'],record['appearance'])
    if json.loads(scene['style_profile'])!=profile or record['style']!=profile:raise ValueError('Resolved appearance differs')
    prepared=record['preparation'];config=prepared['config'];source_dir=Path(prepared['asset_dir'])
    if digest(source_dir/'manifest.json')!=prepared['source_manifest_sha256']:raise ValueError('Source manifest changed')
    source=json.loads((source_dir/'manifest.json').read_text())
    roles={r['chain']:'nucleic_acid' if r['nucleic_acid'] else 'protein' for r in source['chains']}
    used_roles=set(roles[c] for c in config['chains'])|{'membrane_heads','membrane_tails','membrane_surface'}
    if record['environment']:used_roles.add('environment')
    shaders=check_materials(profile,used_roles&set(profile.get('surface_roles',{}))) if profile.get('surface_roles') else None
    colors={c:profile['palette'][key] for c,key in [('R','receptor'),('A','alpha'),('B','beta'),('G','gamma')]}
    colors.update(config.get('chain_colors',{}));colors.update(profile.get('chain_colors',{}))
    if record['resolved_chain_colors']!={c:colors[c] for c in config['chains']}:raise ValueError('Identity color record differs')
    def equal(a,b,message):
        if np.shape(a)!=np.shape(b) or not np.allclose(a,b,atol=1e-6,rtol=1e-5):raise ValueError(message)
    def material(mat,color,role):
        if mat is None or mat.get('surface_role')!=role:raise ValueError('Active surface role differs')
        socket=mat.node_tree.nodes['Principled BSDF'].inputs['Base Color']
        if socket.is_linked:raise ValueError('Identity color overridden')
        equal(socket.default_value,(*linear(color),1),'Active identity color differs')
    for row in record['objects']:
        obj=scene.objects[row['object']]
        if len(obj.material_slots)!=1 or any(p.material_index!=0 for p in obj.data.polygons):raise ValueError('Receptor material slots differ')
        material(obj.material_slots[0].material,colors[row['chain']],roles[row['chain']])
        if len(obj.modifiers)!=1 or obj.modifiers[0].type!='SMOOTH' or obj.modifiers[0].iterations!=2:
            raise ValueError('Source smoothing changed')
        equal(obj.modifiers[0].factor,.35,'Source smoothing changed')
    membrane={'heads':profile['palette']['heads'],'tails':profile['palette']['tails'],'core':profile['palette']['heads']}
    membrane.update(config['display'].get('membrane_palette',{}).get(record['build_parameters']['style'],{}))
    membrane.update(profile.get('membrane_palette',{}))
    if record['resolved_membrane_palette']!=membrane:raise ValueError('Membrane color record differs')
    lipid=scene.objects['Intact two-tail lipid instances'];heads=bpy.data.materials['Lipid heads'];tails=bpy.data.materials['Lipid tails']
    if displayed_materials(lipid)!={heads,tails}:raise ValueError('Lipid material path disconnected or replaced')
    material(heads,membrane['heads'],'membrane_heads');material(tails,membrane['tails'],'membrane_tails')
    core=scene.objects['Continuum bilayer core']
    if len(core.material_slots)!=1:raise ValueError('Membrane core material slots differ')
    material(core.material_slots[0].material,membrane['core'],'membrane_surface')
    data=np.load(Path(record['build_parameters']['distribution_dir'])/'distribution.npz')
    points=np.empty((len(core.data.vertices),3),np.float32);core.data.vertices.foreach_get('co',points.ravel())
    if not np.array_equal(points,data['shell_vertices'].astype(np.float32)) or not np.array_equal([list(p.vertices) for p in core.data.polygons],data['shell_faces']):
        raise ValueError('Membrane core geometry changed')
    equal(np.array(core.matrix_world),np.eye(4),'Membrane core moved')
    if core.modifiers or core.hide_render:raise ValueError('Membrane core modifier or visibility changed')
    # Legend swatches are instanced paint materials, with the same identity
    # colors before the selected scene view transform and artistic recipe.
    for item in config['display']['legend']:
        obj=scene.objects[item['label']+' legend'];assigned=displayed_materials(obj)
        if len(assigned)!=1 or None in assigned:raise ValueError('Legend material disconnected')
        mat=next(iter(assigned));rgb=np.array(linear(colors[item['chain']]))
        ramps=[n for n in mat.node_tree.nodes if n.type=='VALTORGB'];emits=[n for n in mat.node_tree.nodes if n.type=='EMISSION']
        if len(ramps)!=1 or len(emits)!=1 or not emits[0].inputs['Color'].links or emits[0].inputs['Color'].links[0].from_node!=ramps[0]:
            raise ValueError('Legend color shader path differs')
        for stop,factor in zip(ramps[0].color_ramp.elements,[.92,1.08]):equal(stop.color,(* (rgb*factor),1),'Legend identity color differs')
    for obj in scene.objects:
        if obj.type!='FONT':continue
        for font,key in [(obj.data.font,'font_zh_sha256'),(obj.data.font_bold,'font_en_sha256')]:
            if not font.packed_file or hashlib.sha256(bytes(font.packed_file.data)).hexdigest()!=record[key]:raise ValueError('Packed language font changed')
    scale=config['display']['camera_scale_nm']/29.375
    target=np.array(config['display']['camera_target_nm'])
    def aimed(obj):
        direction=target-np.array(obj.location);direction/=np.linalg.norm(direction)
        equal(-np.array(obj.matrix_world.to_3x3())[:,2],direction,'Requested aim differs')
    camera=scene.camera
    equal(camera.location,config['display']['camera_position_nm'],'Prepared camera position differs')
    equal(camera.data.ortho_scale,config['display']['camera_scale_nm'],'Prepared camera scale differs')
    if camera.data.type!='ORTHO' or camera.data.dof.use_dof:raise ValueError('Prepared sharp camera changed')
    aimed(camera)
    for row in profile['lighting']:
        obj=scene.objects.get(row['name'])
        if obj is None or obj.type!='LIGHT' or obj.data.type!='AREA':raise ValueError('Requested area light missing')
        equal(obj.location,np.array(row['position'])*scale,'Light position differs')
        equal(obj.data.energy,row['power']*scale**2,'Light power differs')
        equal(obj.data.size,row['size']*scale,'Light size differs')
        equal(obj.data.color,linear(row['color']),'Light color differs')
        aimed(obj)
    equal(scene.world.node_tree.nodes['Background'].inputs['Color'].default_value,(*linear(profile['world']),1),'World color differs')
    equal(scene.world.node_tree.nodes['Background'].inputs['Strength'].default_value,.35,'World strength differs')
    if scene.view_settings.view_transform!=profile.get('view_transform','AgX'):raise ValueError('View transform differs')
    if scene.view_settings.look!=profile.get('look','AgX - Medium High Contrast' if scene.view_settings.view_transform=='AgX' else 'None'):raise ValueError('View look differs')
    equal(scene.view_settings.exposure,profile.get('exposure',0),'Exposure differs')
    tree=scene.compositing_node_group
    def node(label):
        result=[n for n in tree.nodes if n.label==label]
        if len(result)!=1:raise ValueError('Missing or ambiguous compositor node: '+label)
        return result[0]
    def linked(socket,expected,message):
        if len(socket.links)!=1 or socket.links[0].from_socket!=expected:raise ValueError(message)
    env=record['environment'];mesh_cache={}
    if env:
        if {o.name for o in scene.objects if o.get('environment_only')}!={r['object'] for r in env['objects']}:raise ValueError('Background inventory changed')
        # Excluded objects have stale matrix_world after reopening. Inspect each
        # object in the evaluated depth layer that actually renders it.
        matrices={};original_layer=bpy.context.window.view_layer
        try:
            for index,depth in enumerate(['far','middle','near']):
                bpy.context.window.view_layer=scene.view_layers['Cytosol '+depth+' render']
                bpy.context.view_layer.update();graph=bpy.context.evaluated_depsgraph_get()
                for row in env['objects']:
                    if row['depth_layer']==index:matrices[row['object']]=np.array(scene.objects[row['object']].evaluated_get(graph).matrix_world)
        finally:bpy.context.window.view_layer=original_layer
        palette={'ribosome':'#638e8d','actin':'#8b829a','ubiquitin':'#728d9b'};palette.update(profile.get('environment_palette',{}))
        for row in env['objects']:
            obj=scene.objects[row['object']];name=row['asset'];path=Path(env['asset_directory'])/(name+'.npz')
            if name not in mesh_cache:
                if digest(path)!=env['asset_sha256'][name]:raise ValueError('Background source changed')
                mesh_cache[name]=dict(np.load(path))
            reference=mesh_cache[name];points=np.empty((len(obj.data.vertices),3),np.float32);obj.data.vertices.foreach_get('co',points.ravel())
            if not np.array_equal(points,reference['vertices']) or not np.array_equal([list(p.vertices) for p in obj.data.polygons],reference['faces']):raise ValueError('Background source geometry changed')
            equal(matrices[row['object']],row['matrix_nm'],'Background source placement changed')
            if obj.hide_render or obj.modifiers or len(obj.material_slots)!=1:raise ValueError('Background display changed')
            material(obj.material_slots[0].material,palette[name],'environment')
        base=node('Cellular environment base color');previous=base.outputs[0]
        equal(previous.default_value,(*linear(profile['background']),1),'Background base color differs')
        beauty=scene.view_layers[scene['compositor_beauty_layer']]
        for index,depth in enumerate(['far','middle','near']):
            group='Cytosol '+depth;layer=scene.view_layers[group+' render']
            if not beauty.layer_collection.children[group].exclude:raise ValueError('Background leaks into sharp subject layer')
            allowed={group,scene.camera.users_collection[0].name}
            if {c.name for c in layer.layer_collection.children if not c.exclude}!=allowed:raise ValueError('Subject or labels leak into blurred background')
            blur=node(layer.name+' defocus');over=node(layer.name+' background strength')
            equal(blur.inputs['Size'].default_value,[settings.get('blur_pixels_1280',[18,9,28])[index]*scene.render.resolution_x/1280]*2,'Depth blur differs')
            equal(over.inputs['Factor'].default_value,settings.get('opacity',[.26,.24,.14])[index],'Depth opacity differs')
            if len(blur.inputs['Image'].links)!=1:raise ValueError('Depth render input disconnected')
            render=blur.inputs['Image'].links[0].from_node
            if render.type!='R_LAYERS' or render.layer!=layer.name:raise ValueError('Wrong depth render layer')
            linked(over.inputs['Foreground'],blur.outputs['Image'],'Depth blur bypassed')
            linked(over.inputs['Background'],previous,'Depth order changed');previous=over.outputs['Image']
        linked(node('Straight background plus premultiplied rendered foreground').inputs['Background'],previous,'Depth backdrop not used')
    output=next(n for n in tree.nodes if n.type=='GROUP_OUTPUT' and n.is_active_output)
    texts=[n for n in tree.nodes if n.type=='R_LAYERS' and n.layer==scene['compositor_label_layer']]
    if len(texts)!=1:raise ValueError('Protected text layer changed')
    text=texts[0]
    if len(output.inputs['Image'].links)!=1:raise ValueError('Protected text output disconnected')
    final=output.inputs['Image'].links[0].from_node
    if final.type!='ALPHAOVER':raise ValueError('Protected text no longer last')
    linked(final.inputs['Foreground'],text.outputs['Image'],'Protected text routed through filters')
    return {'passed':True,'role_materials':len(shaders['materials']) if shaders else 0,
            'source_background_instances':len(env['objects']) if env else 0,'depth_layers':3 if env else 0,
            'lights':len(profile['lighting']),'scope':'Effective shaders/colors, prepared background/core geometry, requested lights, packed fonts and actual protected depth graph; not biological abundance, calibrated optics or final visual quality.'}
