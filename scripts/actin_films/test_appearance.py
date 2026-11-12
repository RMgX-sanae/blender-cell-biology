"""Reject changed native appearance; restore the saved chapter without edits."""
import argparse,hashlib,json,sys
from pathlib import Path
import bpy
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_appearance import audit


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args(sys.argv[sys.argv.index('--')+1:]);path=Path(bpy.data.filepath)
    original=hashlib.sha256(path.read_bytes()).hexdigest();baseline=audit(bpy.context.scene)
    record=json.loads(bpy.context.scene['actin_appearance'])
    if not record['environment']:p.error('Use a styled chapter with source context')
    cases={'role':'Wrong object surface role','color':'Object identity color differs','shader':'Saved shader differs',
           'light':'Light power differs','camera':'Camera transform differs','scale':'Requested camera scale ignored','camera_animation':'Authored frontal camera changed',
           'font':'Packed bilingual font differs','ink':'Protected text color differs','caption':'Caption identity',
           'typography':'Foreground leaks into typography','background_shape':'Environment source geometry changed',
           'background_move':'Environment placement changed','blur':'Context blur width differs',
           'beat':'Timed teaching beat differs','labels':'Published layers do not describe'}
    if bpy.context.scene['chapter_id']=='a03':cases['fade']='Editorial cut timing differs'
    if bpy.context.scene['selected_compositor']=='neutral':cases['beauty']='Neutral beauty output bypassed'
    results=[]
    for case,diagnostic in cases.items():
        bpy.ops.wm.open_mainfile(filepath=str(path));scene=bpy.context.scene;tree=scene.compositing_node_group
        actor=next(o for o in scene.objects if o.get('atoms_file'));mat=actor.data.materials[0]
        title=next(o for o in scene.objects if o.get('chapter_title'))
        beat=next(o for o in scene.objects if o.get('chapter_beat_index',-1)==0)
        def node(label):return next(n for n in tree.nodes if n.label==label)
        if case=='role':mat['surface_role']='environment'
        if case=='color':mat.node_tree.nodes['Principled BSDF'].inputs['Base Color'].default_value[0]+=.2
        if case=='shader':mat.node_tree.nodes['Principled BSDF'].inputs['Roughness'].default_value=.9
        if case=='light':scene.objects[record['studio']['lighting'][0]['name']].data.energy+=1000
        if case=='camera':scene.camera.location.x+=1
        if case=='scale':scene.camera.data.ortho_scale*=.7
        if case=='camera_animation':scene.camera.data.keyframe_insert('ortho_scale',frame=1)
        if case=='font':title.data.font=title.data.font_bold
        if case=='ink':next(n for n in title.data.materials[0].node_tree.nodes if n.type=='EMISSION').inputs['Color'].default_value[0]+=.2
        if case=='caption':title.data.body+='?'
        if case=='typography':scene.view_layers[scene['compositor_label_layer']].layer_collection.children['Molecules'].exclude=False
        if case=='background_shape':scene.objects[record['environment']['objects'][0]['object']].data.vertices[0].co.x+=.2
        if case=='background_move':scene.objects[record['environment']['objects'][0]['object']].location.x+=1
        if case=='blur':node('Interaction context far render defocus').inputs['Size'].default_value=(1,1)
        if case=='beat':beat.animation_data_clear();beat.hide_render=True
        if case=='labels':
            output=next(n for n in tree.nodes if n.type=='GROUP_OUTPUT' and n.is_active_output)
            tree.links.new(node('Interaction context far render defocus').outputs['Image'],output.inputs['Image'].links[0].from_node.inputs['Foreground'])
        if case=='fade':tree.animation_data_clear()
        if case=='beauty':
            render=next(n for n in tree.nodes if n.type=='R_LAYERS' and n.layer==scene['compositor_beauty_layer'])
            if record['editorial_fade_nodes']:tree.links.new(render.outputs['Image'],node('Authored chapter fade beauty').inputs['Foreground'])
            else:
                output=next(n for n in tree.nodes if n.type=='GROUP_OUTPUT' and n.is_active_output)
                tree.links.new(render.outputs['Image'],output.inputs['Beauty']);tree.links.new(render.outputs['Image'],output.inputs['Image'].links[0].from_node.inputs['Background'])
        bpy.context.view_layer.update()
        try:audit(scene)
        except ValueError as error:
            assert diagnostic in str(error),(case,str(error));results.append({'case':case,'rejected':True,'diagnostic':str(error)})
        else:raise AssertionError('Accepted '+case)
    bpy.ops.wm.open_mainfile(filepath=str(path))
    assert audit(bpy.context.scene)==baseline and hashlib.sha256(path.read_bytes()).hexdigest()==original
    args.output.write_text(json.dumps({'passed':True,'cases':results,'restored_positive':True,'saved_file_unchanged':True},indent=2)+'\n')
    print(json.dumps({'passed':True,'native_negative_cases':len(results),'saved_file_unchanged':True}))


if __name__=='__main__':main()
