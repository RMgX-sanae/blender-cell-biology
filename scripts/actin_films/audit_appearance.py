"""Inspect effective chapter appearance, timed text and preserved editorial dip."""
import argparse,hashlib,json,sys
from pathlib import Path
import bpy,numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from appearance_config import load_saved,chapter_profile
from chapter_appearance import label_state
from audit_interaction_scene import displayed_materials
from audit_surface_materials import check_materials
from detail_scene import linear
from interaction_environment import audit_environment
from export_composite_layers import published_sources


def audit(scene):
    if 'actin_appearance' not in scene:return None
    folder=Path(bpy.data.filepath).parent;chapter=json.loads((folder/'manifest.json').read_text())
    options=json.loads((folder.parents[1]/'project.json').read_text());record=json.loads(scene['actin_appearance'])
    profile=chapter_profile(load_saved(options),scene['chapter_id'])
    if profile!=record['profile'] or chapter['appearance']!=record:raise ValueError('Saved chapter appearance differs from its project request')
    if record['profile_sha256']!=options['style_profile_sha256'] or record['environment_sha256']!=options.get('environment_sha256'):raise ValueError('Appearance rebuild signature differs')
    def equal(a,b,message):
        if np.shape(a)!=np.shape(b) or not np.allclose(a,b,atol=2e-6,rtol=1e-5):raise ValueError(message)
    palette=record['palette']
    if any(palette[k]!=v for k,v in profile.get('palette',{}).items()):raise ValueError('Requested semantic palette ignored')
    roles=set();bodies=[]
    for group in ('Molecules','Schematics'):
        for obj in bpy.data.collections[group].objects:
            if obj.type not in ('MESH','CURVE'):continue
            role=('cofactor' if obj.get('ccd') else 'protein') if group=='Molecules' else 'schematic_linker' if obj.type=='CURVE' else 'membrane_symbol' if obj.name=='Membrane lipids schematic' else 'schematic_symbol'
            if group=='Schematics' and obj.type=='MESH':
                evaluated=obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
                mats={evaluated.material_slots[p.material_index].material if p.material_index<len(evaluated.material_slots) else None for p in evaluated.data.polygons}
            else:mats=displayed_materials(obj) if obj.type=='MESH' else set(obj.data.materials)
            if not mats or None in mats or any(m.get('surface_role')!=role for m in mats):raise ValueError('Wrong object surface role: '+obj.name)
            for mat in mats:
                key=mat.get('palette_key')
                if key not in palette:raise ValueError('Unknown object palette identity')
                socket=mat.node_tree.nodes['Principled BSDF'].inputs['Base Color']
                if socket.is_linked:raise ValueError('Object identity color overridden')
                equal(socket.default_value,(*linear(palette[key]),1),'Object identity color differs')
            roles.add(role);bodies.append(obj.name)
    if record['environment']:roles.add('environment')
    check_materials(profile,roles&profile.get('surface_roles',{}).keys())
    camera=scene.camera;settings=profile.get('camera',{});pan=settings.get('center_offset_nm',[0,0])
    if camera.data.type!='ORTHO' or camera.data.dof.use_dof or camera.parent or camera.animation_data or camera.data.animation_data or camera.constraints:raise ValueError('Authored frontal camera changed')
    equal(np.array(camera.matrix_world),record['camera']['matrix'],'Camera transform differs')
    equal(camera.location,(pan[0],-100,pan[1]),'Requested camera pan ignored')
    equal(camera.data.ortho_scale,record['base_scene_scale_nm']*settings.get('scale_multiplier',1),'Requested camera scale ignored')
    studio=record['studio'];lights={o.name:o for o in scene.objects if o.type=='LIGHT'}
    if set(lights)!={row['name'] for row in studio['lighting']}:raise ValueError('Light inventory differs')
    for row in studio['lighting']:
        obj=lights[row['name']]
        if obj.data.type!='AREA' or obj.hide_render or obj.parent or obj.animation_data or obj.data.animation_data or obj.constraints:raise ValueError('Light type or visibility differs')
        equal(obj.location,row['position'],'Light position differs');equal(obj.data.energy,row['power'],'Light power differs')
        equal(obj.data.size,row['size'],'Light size differs');equal(obj.data.color,linear(row['color']),'Light color differs')
        direction=np.array([pan[0],0,pan[1]])-np.array(obj.location);direction/=np.linalg.norm(direction)
        equal(-np.array(obj.matrix_world.to_3x3())[:,2],direction,'Light aim differs')
    for key in ('lighting','world','view_transform','look','exposure'):
        if key in profile and studio[key]!=profile[key]:raise ValueError('Requested studio field ignored: '+key)
    world=scene.world.node_tree.nodes['Background'];equal(world.inputs['Color'].default_value,(*linear(studio['world']),1),'World color differs')
    equal(world.inputs['Strength'].default_value,.35,'World strength differs')
    if scene.view_settings.view_transform!=studio.get('view_transform','AgX') or scene.view_settings.look!=studio.get('look','AgX - Medium High Contrast' if scene.view_settings.view_transform=='AgX' else 'None'):raise ValueError('Display transform differs')
    equal(scene.view_settings.exposure,studio.get('exposure',0),'Exposure differs')
    if scene.get('selected_compositor')!=profile.get('render_recipe','neutral'):raise ValueError('Render recipe differs')
    tree=scene.compositing_node_group;ports=published_sources(tree)
    def node(label):
        found=[n for n in tree.nodes if n.label==label]
        if len(found)!=1:raise ValueError('Missing compositor node: '+label)
        return found[0]
    def linked(socket,expected,message):
        if len(socket.links)!=1 or socket.links[0].from_socket!=expected:raise ValueError(message)
    base=node('Straight background plus premultiplied rendered foreground')
    beauty=scene.view_layers[scene['compositor_beauty_layer']];typography=scene.view_layers[scene['compositor_label_layer']]
    if any(name not in beauty.objects for name in bodies):raise ValueError('Foreground excluded from sharp beauty')
    if not beauty.layer_collection.children[scene['compositor_label_collection']].exclude:raise ValueError('Typography leaks into beauty filters')
    if {c.name for c in typography.layer_collection.children if not c.exclude}!={scene['compositor_label_collection']}:raise ValueError('Foreground leaks into typography')
    if record['environment']:audit_environment(scene,record['environment'],{},node,linked,equal,base,studio['background'])
    else:equal(base.inputs['Background'].default_value,(*linear(studio['background']),1),'Plain background differs')
    text_nodes=[n for n in tree.nodes if n.type=='R_LAYERS' and n.layer==typography.name]
    if len(text_nodes)!=1:raise ValueError('Text render layer differs')
    text=text_nodes[0]
    if record['editorial_fade_nodes']:
        if scene['chapter_id']!='a03' or len(record['editorial_fade_nodes'])!=2:raise ValueError('Unexpected editorial cut')
        for name,port in [('beauty','Beauty'),('labels','Labels')]:
            fade=node('Authored chapter fade '+name)
            if ports[port]!=fade.outputs['Image']:raise ValueError('Editorial fade bypassed')
            if name=='labels':linked(fade.inputs['Foreground'],text.outputs['Image'],'Labels routed through image filters')
            for t,value in [(0,1),(.285,1),(.300,0),(.307,0),(.324,1),(1,1)]:
                frame=min(scene.frame_end,max(1,1+round(t*(scene.frame_end-1))));scene.frame_set(frame)
                equal(fade.inputs['Factor'].default_value,value,'Editorial cut timing differs')
    else:
        if scene['chapter_id']=='a03':raise ValueError('Authored editorial cut missing')
        if ports['Labels']!=text.outputs['Image']:raise ValueError('Labels routed through image filters')
    beauty_source=node('Authored chapter fade beauty').inputs['Foreground'].links[0].from_socket if record['editorial_fade_nodes'] else ports['Beauty']
    if scene['selected_compositor']=='neutral' and beauty_source!=base.outputs['Image']:raise ValueError('Neutral beauty output bypassed')
    pending=[beauty_source.node];reached=set()
    while pending:
        current=pending.pop()
        if current in reached:continue
        reached.add(current);pending.extend(link.from_node for socket in current.inputs for link in socket.links)
    if base not in reached:raise ValueError('Foreground composite does not reach beauty output')
    if {o.name for o in scene.objects if o.type=='FONT'}!={r['object'] for r in record['fonts']}:raise ValueError('Text inventory differs')
    for row in record['fonts']:
        obj=scene.objects[row['object']]
        if obj.data.body!=row['body'] or not obj.get('compositor_preserve_text_color'):raise ValueError('Caption identity or protected ink changed')
        for font,sha in [(obj.data.font,record['font_zh_sha256']),(obj.data.font_bold,record['font_en_sha256'])]:
            if not font.packed_file or hashlib.sha256(bytes(font.packed_file.data)).hexdigest()!=sha:raise ValueError('Packed bilingual font differs')
        if any(fmt.use_bold!=(ord(char)<0x400 or 0x2010<=ord(char)<=0x22ff) for char,fmt in zip(obj.data.body,obj.data.body_format)):raise ValueError('Mixed language font routing differs')
        mat=obj.data.materials[0];emits=[n for n in mat.node_tree.nodes if n.type=='EMISSION'];outputs=[n for n in mat.node_tree.nodes if n.type=='OUTPUT_MATERIAL' and n.is_active_output]
        if len(emits)!=1 or len(outputs)!=1 or len(outputs[0].inputs['Surface'].links)!=1 or outputs[0].inputs['Surface'].links[0].from_node!=emits[0]:raise ValueError('Protected text is not unlit ink')
        expected=(.012,.025,.032,1) if row['ink_role']=='diagram' else (*linear(palette['ink' if row['ink_role']=='foreground' else row['palette_key']]),1)
        equal(emits[0].inputs['Color'].default_value,expected,'Protected text color differs')
        equal(emits[0].inputs['Strength'].default_value,row['strength'],'Protected text strength differs')
    boundaries=chapter.get('beat_boundaries',[0,1/3,2/3,1]);frames=scene.frame_end;checked=0
    try:
        for frame in range(1,frames+1):
            scene.frame_set(frame);states=label_state(scene)
            for name,state in states.items():
                obj=scene.objects[name]
                if obj.get('chapter_title'):
                    if state['hidden'] or obj.data.body!=chapter['title']:raise ValueError('Chapter title missing')
                elif 'chapter_beat_index'in obj:
                    i=obj['chapter_beat_index'];first=1+round(boundaries[i]*(frames-1));last=frames if i==2 else 1+round((boundaries[i+1]-1/frames)*(frames-1))
                    if state['hidden']==(first<=frame<=last) or obj.data.body!=chapter['beats'][i]:raise ValueError('Timed teaching beat differs')
                else:continue
                if not state['hidden'] and any(not (0<=p[0]<=1 and 0<=p[1]<=1 and p[2]>0) for p in state['bounds']):raise ValueError('Chapter heading leaves the image')
            checked+=1
    finally:scene.frame_set(1);bpy.context.view_layer.update()
    return {'passed':True,'roles':sorted(roles),'source_and_diagram_objects':len(bodies),'lights':len(lights),
            'context_instances':len(record['environment']['objects']) if record['environment'] else 0,
            'timed_label_frames':checked,'scope':'Effective appearance, protected ink and authored heading/cut timing. Source geometry, motion, contact and final pixels require their separate checks.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(sys.argv[sys.argv.index('--')+1:]);result=audit(bpy.context.scene)
    args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
