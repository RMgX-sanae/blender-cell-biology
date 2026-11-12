"""Exercise invalid requests and native scene changes without saving mutations."""
import argparse,hashlib,json,sys,tempfile
from pathlib import Path
import bpy
sys.path.insert(0,str(Path(__file__).resolve().parent))
from vatpase_appearance import read_profile
from audit_vatpase_appearance import audit


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args(sys.argv[sys.argv.index('--')+1:]);path=Path(bpy.data.filepath);original=hashlib.sha256(path.read_bytes()).hexdigest()
    baseline=audit(bpy.context.scene);record=json.loads(bpy.context.scene['vatpase_appearance'])
    if not record['environment']:p.error('Use a styled scene with source context')
    manifest_path=path.parent/'manifest.json';manifest_bytes=manifest_path.read_bytes()
    try:
        broken=json.loads(manifest_bytes);broken['labels']=[];manifest_path.write_text(json.dumps(broken))
        try:audit(bpy.context.scene)
        except ValueError as error:assert 'protected caption inventory' in str(error)
        else:raise AssertionError('Accepted missing caption inventory')
    finally:manifest_path.write_bytes(manifest_bytes)
    invalid=[{'camera':{'type':'PERSP'}},{'camera':{'scale_multiplier':0}},{'camera':{'center_offset_nm':[0,float('nan')]}},
             {'palette':{'unknown':'FFFFFF'}},{'palette':{'A':'invalid'}},{'surface_roles':{'unknown':{}}},{'lighting':[]},
             {'environment':{'compartment':'lumen'}},{'exposure':float('inf')}]
    with tempfile.TemporaryDirectory() as folder:
        profile=Path(folder)/'style.json'
        for case in invalid:
            profile.write_text(json.dumps(case))
            try:read_profile(profile)
            except ValueError:pass
            else:raise AssertionError('Accepted invalid appearance '+str(case))
    cases={'role':'material identity differs','color':'Semantic color differs','shader':'Saved shader differs',
           'smoothing':'Source envelope smoothing differs','camera':'Camera transform differs','scale':'Requested scale ignored',
           'light':'Light power differs','font':'Packed bilingual font differs','ink':'Text ink differs','caption':'Caption identity differs',
           'typography':'Foreground leaks into text layer','background':'Environment placement changed','blur':'Context blur width differs',
           'highlight':'ATP highlight timing differs'}
    if bpy.context.scene['selected_compositor']=='neutral':cases['beauty']='Neutral beauty bypassed'
    results=[]
    for case,diagnostic in cases.items():
        bpy.ops.wm.open_mainfile(filepath=str(path));scene=bpy.context.scene;tree=scene.compositing_node_group
        obj=scene.objects['6VQ6 chain A'];mat=obj.data.materials[0];title=scene.objects['V 型 ATPase']
        def node(label):return next(n for n in tree.nodes if n.label==label)
        if case=='role':mat['surface_role']='environment'
        if case=='color':mat.node_tree.nodes['Principled BSDF'].inputs['Base Color'].default_value[0]+=.2
        if case=='shader':mat.node_tree.nodes['Principled BSDF'].inputs['Roughness'].default_value=.9
        if case=='smoothing':obj.modifiers[0].factor=.6
        if case=='camera':scene.camera.location.x+=1
        if case=='scale':scene.camera.data.ortho_scale*=.7
        if case=='light':scene.objects[record['studio']['lighting'][0]['name']].data.energy+=1000
        if case=='font':title.data.font=title.data.font_bold
        if case=='ink':next(n for n in title.data.materials[0].node_tree.nodes if n.type=='EMISSION').inputs['Color'].default_value[0]+=.2
        if case=='caption':title.data.body+='?'
        if case=='typography':scene.view_layers[scene['compositor_label_layer']].layer_collection.children['Pump'].exclude=False
        if case=='background':scene.objects[record['environment']['objects'][0]['object']].location.x+=1
        if case=='blur':node('Interaction context far render defocus').inputs['Size'].default_value=(1,1)
        if case=='highlight':mat.node_tree.animation_data_clear();mat.node_tree.nodes['Principled BSDF'].inputs['Emission Strength'].default_value=0
        if case=='beauty':
            render=next(n for n in tree.nodes if n.type=='R_LAYERS' and n.layer==scene['compositor_beauty_layer'])
            out=next(n for n in tree.nodes if n.type=='GROUP_OUTPUT' and n.is_active_output)
            tree.links.new(render.outputs['Image'],out.inputs['Beauty']);tree.links.new(render.outputs['Image'],out.inputs['Image'].links[0].from_node.inputs['Background'])
        bpy.context.view_layer.update()
        try:audit(scene)
        except ValueError as error:
            assert diagnostic in str(error),(case,str(error));results.append({'case':case,'rejected':True,'diagnostic':str(error)})
        else:raise AssertionError('Accepted native mutation '+case)
    bpy.ops.wm.open_mainfile(filepath=str(path));assert audit(bpy.context.scene)==baseline
    assert hashlib.sha256(path.read_bytes()).hexdigest()==original
    assert manifest_path.read_bytes()==manifest_bytes
    report={'passed':True,'invalid_profile_cases':len(invalid),'native_cases':results,'missing_caption_inventory_rejected':True,'restored_positive':True,'saved_file_unchanged':True,'manifest_unchanged':True}
    args.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'passed':True,'invalid_profiles':len(invalid),'native_cases':len(results)}))


if __name__=='__main__':main()
