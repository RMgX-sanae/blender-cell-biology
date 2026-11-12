"""Check a saved generic molecular sequence against its qualified prepared data."""
import argparse,json,sys
from pathlib import Path
import bpy,numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_phosphotransfer_scene import vertices
from build_ubiquitin_transfer_scene import digest
from coordinate_animation import verify_point_animation
from molecular_sequence_assets import load_assets
from audit_surface_materials import check_materials
from audit_interaction_scene import displayed_materials
from detail_scene import linear
from sequence_presentation import presentation,source_positions,coordinates_at,caption_windows


def check_appearance(scene,record,root):
    appearance=record.get('appearance')
    if appearance is None:return None
    assert appearance==json.loads((root/'style-profile.json').read_text()),'Saved appearance record differs'
    assert digest(root/'style-profile.json')==record['appearance_sha256']==scene['sequence_style_profile_sha256'],'Appearance fingerprint differs'
    roles=set(appearance['component_roles'].values())
    if 'environment' in appearance['surface_roles']:roles.add('environment')
    shaders=check_materials(appearance,roles)
    colors=dict(record['display']['colors'][record['style']]);colors.update(appearance.get('semantic_colors',{}))
    assert record['resolved_colors']==colors,'Resolved semantic colors differ'
    def equal(actual,expected,label):
        assert np.allclose(actual,expected,rtol=1e-5,atol=1e-6),label
    def color_check(obj,material,color):
        assert displayed_materials(obj)=={material},'Actual material assignment differs: '+obj.name
        equal(material.diffuse_color,color,'Display color differs: '+obj.name)
        socket=material.node_tree.nodes.get('Principled BSDF').inputs['Base Color']
        assert not socket.is_linked,'Unreviewed base-color override: '+obj.name
        equal(socket.default_value,color,'Surface color differs: '+obj.name)
    for view in record['display']['views']:
        for key,settings in view['actors'].items():
            obj=scene.objects[key+' '+view['id']];material=bpy.data.materials[view['id']+' '+key]
            assert material['surface_role']==appearance['component_roles'][key],'Component material role differs'
            color=settings.get('color',colors[settings.get('color_key',key)])
            color_check(obj,material,(*linear(color),1))
    palette=appearance.get('environment_palette',{})
    fallback=dict(zip(record['display']['background_components'],['#778f9a','#789994','#8f839e']))
    for entry in record['environment']:
        obj=scene.objects[entry['object']];key=entry['component'];material=bpy.data.materials['Background '+key]
        rgb=linear(palette.get(key,fallback[key]))
        if record['style']=='dark' and not palette:rgb=[value*.55 for value in rgb]
        color_check(obj,material,(*rgb,1))
        layer='Molecular background '+('far','middle','near')[entry['depth_layer']]+' render'
        bpy.context.window.view_layer=scene.view_layers[layer];bpy.context.view_layer.update()
        evaluated=obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
        equal(np.asarray(evaluated.matrix_world),entry['matrix'],'Background source placement differs')
        data=np.load(Path(record['asset_directory'])/('background_'+key+'.npz'))
        points=np.empty(len(obj.data.vertices)*3);obj.data.vertices.foreach_get('co',points)
        equal(points.reshape(-1,3),data['vertices'],'Background source shape differs')
        assert np.array_equal([list(p.vertices) for p in obj.data.polygons],data['faces']),'Background source faces differ'
    if 'atom_subdivisions' in appearance:
        for obj in scene.objects:
            for modifier in obj.modifiers:
                if modifier.type!='NODES':continue
                for node in modifier.node_group.nodes:
                    if node.bl_idname=='GeometryNodeMeshIcoSphere':
                        assert node.inputs['Subdivisions'].default_value==appearance['atom_subdivisions'],'Atom display tessellation differs'
    studio=json.loads(scene['style_profile'])
    for key in ('lighting','world','view_transform','look','exposure','background'):
        if key in appearance:assert studio[key]==appearance[key],'Selected studio field differs: '+key
    equal(scene.world.node_tree.nodes['Background'].inputs['Color'].default_value,(*linear(studio['world']),1),'World color differs')
    assert scene.view_settings.view_transform==studio.get('view_transform','AgX'),'View transform differs'
    default_look='AgX - Medium High Contrast' if scene.view_settings.view_transform=='AgX' else 'None'
    assert scene.view_settings.look==studio.get('look',default_look),'View look differs'
    equal(scene.view_settings.exposure,studio.get('exposure',0),'Exposure differs')
    for item in studio['lighting']:
        light=scene.objects[item['name']]
        assert light.type=='LIGHT' and light.data.type=='AREA','Selected area light missing'
        equal(light.location,np.array(item['position'])*1.2,'Light position differs')
        equal(light.data.energy,item['power']*1.2**2,'Light power differs')
        equal(light.data.size,item['size']*1.2,'Light size differs')
        equal(light.data.color,linear(item['color']),'Light color differs')
    settings=appearance.get('environment_compositor',{}).get(record['style'],{})
    assert record['environment_compositor']==settings,'Environment settings record differs'
    def node(label):
        found=[n for n in scene.compositing_node_group.nodes if n.label==label]
        assert len(found)==1,'Missing or ambiguous appearance node: '+label
        return found[0]
    equal(node('Cellular environment base color').outputs[0].default_value,
          (*linear(settings.get('background',studio['background'])),1),'Environment base color differs')
    for index,depth in enumerate(('far','middle','near')):
        prefix='Molecular background '+depth+' render'
        blur=node(prefix+' defocus');strength=node(prefix+' background strength')
        pixels=settings.get('blur_pixels_1280',[18,9,28])[index]*scene.render.resolution_x/1280
        equal(blur.inputs['Size'].default_value,[pixels,pixels],'Environment blur differs')
        equal(strength.inputs['Factor'].default_value,settings.get('opacity',[.26,.24,.14])[index],'Environment opacity differs')
        assert strength.inputs['Foreground'].links[0].from_node==blur,'Environment blur disconnected'
    return {'passed':True,'role_materials':len(shaders['materials']),'source_background_instances':len(record['environment']),
            'lighting_checked':len(studio['lighting']),'depth_layers_checked':3,'profile_sha256':record['appearance_sha256']}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
    if args.output.exists():parser.error('Use a fresh audit output')
    root=Path(bpy.data.filepath).parent;record=json.loads((root/'manifest.json').read_text())
    scene=bpy.context.scene;assets=Path(record['asset_directory']);display=record['display']
    source,data,input_checks=load_assets(assets,display);source_count=source['frames'];frames=record['frames']
    timing=record.get('timing_profile')
    if timing is not None:
        assert timing==json.loads((root/'timing-profile.json').read_text()),'Saved presentation profile differs'
        assert digest(root/'timing-profile.json')==record['timing_profile_sha256']==scene['sequence_timing_profile_sha256'],'Presentation fingerprint differs'
    schedule=presentation(source,timing,digest(assets/'manifest.json'))
    assert record.get('presentation',schedule)==schedule,'Resolved presentation schedule differs'
    assert frames==schedule['frames'] and record['fps']==schedule['fps'],'Presentation duration differs'
    samples=np.array(schedule['source_samples']);sample_frames=np.array(schedule['sample_frames'])
    check_frames=np.unique(np.r_[sample_frames,np.arange(1,frames+2)])
    positions=source_positions(schedule,check_frames)
    render_positions=source_positions(schedule,np.arange(1,frames+2))
    windows=caption_windows(display['captions'],schedule,source_count)
    assert record.get('caption_windows',windows)==windows,'Resolved caption windows differ'
    assert source==record['preparation'] and display==json.loads((root/'display.json').read_text()),'Saved input records differ'
    assert record['kind']=='molecular-sequence-scene' and source['passed'] and source['checks']['passed']
    assert scene.frame_end==frames and scene.render.fps==schedule['fps']
    assert digest(assets/'manifest.json')==record['source_sha256']==scene['sequence_source_sha256']
    assert digest(root/'display.json')==record['display_sha256']==scene['sequence_display_sha256']
    appearance_check=check_appearance(scene,record,root)
    physical_ids=[identity for values in source['atom_ids'].values() for identity in values]
    expected_actors={key+' '+view['id'] for view in display['views'] for key in view['actors']}
    assert set(record['actors'])==expected_actors and len(record['actors'])==len(expected_actors)
    # Coordinate checks read the saved source meshes; generated display geometry
    # is inspected in rendered outputs, not rebuilt thousands of times here.
    for name in record['actors']+[marker['object'] for marker in record['markers']]:
        for modifier in scene.objects[name].modifiers:modifier.show_viewport=False
    error=0.;checked=0
    for view,layer_name in zip(display['views'],['Sequence overview','Sequence detail']):
        bpy.context.window.view_layer=scene.view_layers[layer_name]
        objects=[scene.objects[key+' '+view['id']] for key in view['actors']]
        for obj in objects:
            assert obj['view']==view['id'] and obj['component'] in view['actors']
            assert json.loads(obj['atom_ids'])==source['atom_ids'][obj['component']]
            verify_point_animation(obj,data[obj['component']][samples],sample_frames)
            assert not any(mod.type=='BOOLEAN' for mod in obj.modifiers)
            settings=view['actors'][obj['component']];edges=settings.get('bonds',[])
            if 'bond_file' in settings:edges=np.load(assets/settings['bond_file'])['bonds']
            expected_edges={tuple(sorted(map(int,pair))) for pair in edges}
            assert {tuple(sorted(edge.vertices)) for edge in obj.data.edges}==expected_edges,(obj.name,'display edges differ')
        matrix=np.asarray(record['transforms'][view['id']]);offset=-5.5 if view['id']==display['views'][0]['id'] else 5.5
        for frame,position in zip(check_frames,positions):
            scene.frame_set(int(frame),subframe=float(frame%1));graph=bpy.context.evaluated_depsgraph_get()
            for obj in objects:
                observed=vertices(obj);expected=coordinates_at(data[obj['component']],[position])[0];evaluated=obj.evaluated_get(graph)
                assert observed.shape==expected.shape and np.isfinite(observed).all()
                difference=float(np.max(abs(observed-expected)));error=max(error,difference);assert difference<2e-4,(frame,obj.name,difference)
                actual=np.asarray(evaluated.matrix_world)
                assert np.allclose(actual,matrix,atol=1e-5) and not evaluated.hide_render,(frame,obj.name,'transform/visibility')
                assert np.ptp(np.linalg.svd(actual[:3,:3],compute_uv=False))<1e-5,'Stretched molecular view'
                if obj['component'] in view.get('keep_visible',[]):
                    projected=observed@actual[:3,:3].T+actual[:3,3]
                    assert abs(projected[:,0]-offset).max()<4.3 and abs(projected[:,1]).max()<3.45,'Main source construct clipped'
                checked+=1
            if frame%200==0:print(view['id'],frame,'coordinate samples checked',flush=True)
    bpy.context.window.view_layer=scene.view_layers['Sequence annotations']
    marker_specs={settings['id']+' '+view['id']:(settings,view) for view in display['views']
                  for settings in display.get('highlights',[])+display.get('bond_annotations',[]) if view['id'] in settings['views']}
    assert set(marker_specs)=={marker['object'] for marker in record['markers']},'Display marker inventory differs'
    for marker in record['markers']:
        settings,view=marker_specs[marker['object']]
        for key in ('atoms','component','indices','interval'):
            if key in settings:assert marker[key]==settings[key],'Display marker identity or interval differs'
        expected=(np.stack([data[k][:,index] for k,index in marker['atoms']],axis=1) if 'atoms' in marker
                  else data[marker['component']][:,marker['indices']])
        verify_point_animation(scene.objects[marker['object']],expected[samples],sample_frames)
        matrix=np.asarray(record['transforms'][marker['view']]);offset=-5.5 if view['id']==display['views'][0]['id'] else 5.5
        projected=coordinates_at(expected,render_positions)@matrix[:3,:3].T+matrix[:3,3]
        radius=settings.get('radius_nm',.024)*8.6/view['span_nm']
        visible=(abs(projected[:,:,0]-offset).max(1)+radius<4.3)&(abs(projected[:,:,1]).max(1)+radius<3.45)
        if 'interval' in settings:
            time=render_positions/source_count;visible&=(time>=settings['interval'][0])&(time<settings['interval'][1])
        assert marker['visible']==visible.tolist(),'Annotation visibility differs from source event/window'
        for mod in scene.objects[marker['object']].modifiers:mod.show_viewport=False
    for frame in range(1,frames+2):
        scene.frame_set(frame);graph=bpy.context.evaluated_depsgraph_get()
        for marker in record['markers']:
            obj=scene.objects[marker['object']];matrix=np.asarray(record['transforms'][marker['view']]);evaluated=obj.evaluated_get(graph)
            source_points=(np.stack([data[k][:,index] for k,index in marker['atoms']],axis=1) if 'atoms' in marker
                           else data[marker['component']][:,marker['indices']])
            expected=coordinates_at(source_points,[render_positions[frame-1]])[0]
            assert vertices(obj).shape==expected.shape and np.max(abs(vertices(obj)-expected))<2e-4,(frame,obj.name,'annotation detached')
            assert np.allclose(np.asarray(evaluated.matrix_world),matrix,atol=1e-5)
            assert evaluated.hide_render==(not marker['visible'][frame-1]),(frame,obj.name,'annotation window')
    bpy.context.window.view_layer=scene.view_layers[scene['compositor_label_layer']]
    fonts=[obj for obj in scene.objects if obj.type=='FONT'];captions={row[2] for row in display['captions']}
    for obj in fonts:
        assert obj.data.font.packed_file and obj.data.font_bold.packed_file
        for char,form in zip(obj.data.body,obj.data.body_format):assert form.use_bold==(ord(char)<128 or char in '·—')
    for frame in range(1,frames+1):
        scene.frame_set(frame);graph=bpy.context.evaluated_depsgraph_get()
        active=[obj.data.body for obj in fonts if obj.data.body in captions and not obj.evaluated_get(graph).hide_render]
        expected=[body for start,end,body in windows if start<=frame<=end]
        assert active==expected,(frame,'caption mismatch')
    tree=scene.compositing_node_group
    mask=next(node for node in tree.nodes if node.label=='Select matching premultiplied panel alpha')
    assert mask.inputs['Type'].default_value=='Replace Alpha' and mask.inputs['Alpha'].is_linked
    result={'passed':True,'render_frames':frames,'coordinate_samples_checked_per_view':len(samples),'evaluated_times_per_view':len(check_frames),'actor_frame_checks':checked,
            'presentation_mode':schedule['mode'],'source_samples_retained':source_count+1,'observation_hold_keys':len(samples)-source_count-1,
            'physical_heavy_atom_identities':len(physical_ids),'maximum_coordinate_error_nm':error,'ledger_checks':input_checks['ledger_checks'],
            'required_input_hashes_checked':True,'stored_coordinate_keys_and_linear_timing_checked':True,
            'source_and_display_fingerprints_checked':True,'both_views_and_annotations_checked':True,
            'packed_fonts_and_caption_timing_checked':True,'scope':source['scope'],
            'appearance_check':appearance_check,
            'limits':'Saved coordinate/identity/display checks; reaction-source qualification and rendered clarity remain separate.'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='ledger_checks'},ensure_ascii=False))


if __name__=='__main__':main()
