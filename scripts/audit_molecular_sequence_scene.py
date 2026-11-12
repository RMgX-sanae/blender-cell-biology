"""Check a saved generic molecular sequence against its qualified prepared data."""
import argparse,json,sys
from pathlib import Path
import bpy,numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_phosphotransfer_scene import vertices
from build_ubiquitin_transfer_scene import digest
from coordinate_animation import verify_point_animation
from molecular_sequence_assets import load_assets


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
    if args.output.exists():parser.error('Use a fresh audit output')
    root=Path(bpy.data.filepath).parent;record=json.loads((root/'manifest.json').read_text())
    scene=bpy.context.scene;assets=Path(record['asset_directory']);display=record['display']
    source,data,input_checks=load_assets(assets,display);frames=source['frames']
    assert source==record['preparation'] and display==json.loads((root/'display.json').read_text()),'Saved input records differ'
    assert record['kind']=='molecular-sequence-scene' and source['passed'] and source['checks']['passed']
    assert scene.frame_end==frames and scene.render.fps==source['fps']
    assert digest(assets/'manifest.json')==record['source_sha256']==scene['sequence_source_sha256']
    assert digest(root/'display.json')==record['display_sha256']==scene['sequence_display_sha256']
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
            verify_point_animation(obj,data[obj['component']])
            assert not any(mod.type=='BOOLEAN' for mod in obj.modifiers)
            settings=view['actors'][obj['component']];edges=settings.get('bonds',[])
            if 'bond_file' in settings:edges=np.load(assets/settings['bond_file'])['bonds']
            expected_edges={tuple(sorted(map(int,pair))) for pair in edges}
            assert {tuple(sorted(edge.vertices)) for edge in obj.data.edges}==expected_edges,(obj.name,'display edges differ')
        matrix=np.asarray(record['transforms'][view['id']]);offset=-5.5 if view['id']==display['views'][0]['id'] else 5.5
        for frame in range(1,frames+2):
            scene.frame_set(frame);graph=bpy.context.evaluated_depsgraph_get()
            for obj in objects:
                observed=vertices(obj);expected=data[obj['component']][frame-1];evaluated=obj.evaluated_get(graph)
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
    for marker in record['markers']:
        expected=(np.stack([data[k][:,index] for k,index in marker['atoms']],axis=1) if 'atoms' in marker
                  else data[marker['component']][:,marker['indices']])
        verify_point_animation(scene.objects[marker['object']],expected)
        for mod in scene.objects[marker['object']].modifiers:mod.show_viewport=False
    for frame in range(1,frames+2):
        scene.frame_set(frame);graph=bpy.context.evaluated_depsgraph_get()
        for marker in record['markers']:
            obj=scene.objects[marker['object']];matrix=np.asarray(record['transforms'][marker['view']]);evaluated=obj.evaluated_get(graph)
            expected=(np.array([data[k][frame-1,index] for k,index in marker['atoms']]) if 'atoms' in marker
                      else data[marker['component']][frame-1,marker['indices']])
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
        expected=[body for start,end,body in display['captions'] if start<=(frame-1)/frames<end]
        assert len(active)==len(expected)==1 and active==expected,(frame,'caption mismatch')
    tree=scene.compositing_node_group
    mask=next(node for node in tree.nodes if node.label=='Select matching premultiplied panel alpha')
    assert mask.inputs['Type'].default_value=='Replace Alpha' and mask.inputs['Alpha'].is_linked
    result={'passed':True,'render_frames':frames,'coordinate_samples_checked_per_view':frames+1,'actor_frame_checks':checked,
            'physical_heavy_atom_identities':len(physical_ids),'maximum_coordinate_error_nm':error,'ledger_checks':input_checks['ledger_checks'],
            'required_input_hashes_checked':True,'stored_coordinate_keys_and_linear_timing_checked':True,
            'source_and_display_fingerprints_checked':True,'both_views_and_annotations_checked':True,
            'packed_fonts_and_caption_timing_checked':True,'scope':source['scope'],
            'limits':'Saved coordinate/identity/display checks; reaction-source qualification and rendered clarity remain separate.'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='ledger_checks'},ensure_ascii=False))


if __name__=='__main__':main()
