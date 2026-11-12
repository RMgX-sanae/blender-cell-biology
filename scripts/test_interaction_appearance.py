"""Mutate a saved styled membrane interaction and restore without writing changes."""
import argparse,hashlib,json,sys
from pathlib import Path
import bpy
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_interaction_scene import audit
from audit_interaction_appearance import check_appearance


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);path=Path(bpy.data.filepath);original=hashlib.sha256(path.read_bytes()).hexdigest()
    baseline=audit();assert baseline['passed'];record=json.loads(bpy.context.scene['interaction_appearance'])
    if not record['environment'] or not record['font_zh_sha256']:p.error('Use a styled scene with context and mixed Chinese labels')
    manifest=json.loads((Path(bpy.context.scene['interaction_asset_dir'])/'manifest.json').read_text());rows=[]
    cases={'role':'wrong surface role','color':'Actor identity color differs','shader':'Saved shader differs',
           'light':'Light power differs','light_aim':'Light aim differs','camera':'Saved camera transform differs',
           'scale':'Saved camera ortho_scale differs','legend':'Legend color differs','font':'Packed language font differs',
           'background_shape':'Environment source geometry changed','background_move':'Environment placement changed',
           'blur':'Context blur width differs','layer':'Subject or labels leak','labels':'Labels routed through image filters'}
    for case,diagnostic in cases.items():
        bpy.ops.wm.open_mainfile(filepath=str(path));scene=bpy.context.scene;tree=scene.compositing_node_group
        actor=scene.objects[manifest['actors'][0]['id']];mat=bpy.data.materials[actor.name]
        light=scene.objects[record['studio']['lighting'][0]['name']]
        def node(label):return next(n for n in tree.nodes if n.label==label)
        if case=='role':actor.data.materials[0]=bpy.data.materials['Membrane heads']
        if case=='color':mat.node_tree.nodes['Principled BSDF'].inputs['Base Color'].default_value[0]+=.2
        if case=='shader':mat.node_tree.nodes['Principled BSDF'].inputs['Roughness'].default_value=.9
        if case=='light':light.data.energy+=1000
        if case=='light_aim':light.rotation_euler.x+=.3
        if case=='camera':scene.camera.location.x+=1
        if case=='scale':scene.camera.data.ortho_scale*=.7
        if case=='legend':
            dot=next(o for o in scene.objects if o.get('legend_for_objects'))
            next(n for n in dot.data.materials[0].node_tree.nodes if n.type=='VALTORGB').color_ramp.elements[0].color[0]+=.2
        if case=='font':
            text=next(o for o in scene.objects if o.type=='FONT' and '亚基'in o.data.body);text.data.font=text.data.font_bold
        if case=='background_shape':scene.objects[record['environment']['objects'][0]['object']].data.vertices[0].co.x+=.2
        if case=='background_move':scene.objects[record['environment']['objects'][0]['object']].location.x+=1
        if case=='blur':node('Interaction context far render defocus').inputs['Size'].default_value=(1,1)
        if case=='layer':scene.view_layers['Interaction context far render'].layer_collection.children['Interaction partners'].exclude=False
        if case=='labels':
            output=next(n for n in tree.nodes if n.type=='GROUP_OUTPUT' and n.is_active_output)
            final=output.inputs['Image'].links[0].from_node
            tree.links.new(node('Interaction context far render defocus').outputs['Image'],final.inputs['Foreground'])
        bpy.context.view_layer.update()
        try:check_appearance(scene,manifest)
        except ValueError as error:
            assert diagnostic in str(error),(case,str(error));rows.append({'case':case,'rejected':True,'diagnostic':str(error)})
        else:raise AssertionError('Accepted '+case)
    bpy.ops.wm.open_mainfile(filepath=str(path));restored=audit()
    assert baseline==restored and hashlib.sha256(path.read_bytes()).hexdigest()==original
    a.output.write_text(json.dumps({'passed':True,'cases':rows,'restored_positive':True,'saved_file_unchanged':True},indent=2)+'\n')
    print(json.dumps({'passed':True,'native_negative_cases':len(rows),'saved_file_unchanged':True}))


if __name__=='__main__':main()
