"""Reject actual saved receptor appearance drift; restore without saving mutations."""
import argparse,hashlib,json,sys
from pathlib import Path
import bpy

sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_receptor_distribution_scene import audit
from receptor_appearance_audit import check_appearance


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);path=Path(bpy.data.filepath)
    original=hashlib.sha256(path.read_bytes()).hexdigest();record=json.loads((a.project/'manifest.json').read_text())
    if not record.get('appearance') or not record.get('environment'):p.error('Use a styled scene with all three environment layers')
    baseline=audit(a.project);assert baseline['passed'];rows=[]
    cases={'role':'Active surface role differs','shader':'Saved shader differs',
           'lipid':'Lipid material path disconnected or replaced','legend':'Legend identity color differs',
           'core':'Membrane core geometry changed','background_shape':'Background source geometry changed',
           'background_move':'Background source placement changed','blur':'Depth blur differs',
           'depth_order':'Depth order changed','layer_leak':'Subject or labels leak',
           'light':'Light power differs','camera':'Prepared camera scale differs',
           'font':'Packed language font changed','text':'Protected text routed through filters'}
    for case,diagnostic in cases.items():
        bpy.ops.wm.open_mainfile(filepath=str(path));scene=bpy.context.scene;tree=scene.compositing_node_group
        obj=scene.objects[record['objects'][0]['object']]
        def node(label):return next(n for n in tree.nodes if n.label==label)
        if case=='role':obj.data.materials[0]=bpy.data.materials['Lipid heads']
        if case=='shader':obj.data.materials[0].node_tree.nodes['Principled BSDF'].inputs['Roughness'].default_value=.9
        if case=='lipid':
            group=scene.objects['Intact two-tail lipid instances'].modifiers[0].node_group
            next(n for n in group.nodes if n.type=='SET_MATERIAL' and n.inputs['Material'].default_value==bpy.data.materials['Lipid heads']).inputs['Material'].default_value=bpy.data.materials['Lipid tails']
        if case=='legend':
            name=record['preparation']['config']['display']['legend'][0]['label']+' identity'
            next(n for n in bpy.data.materials[name].node_tree.nodes if n.type=='VALTORGB').color_ramp.elements[0].color[0]+=.2
        if case=='core':scene.objects['Continuum bilayer core'].data.vertices[0].co.x+=.2
        if case=='background_shape':scene.objects[record['environment']['objects'][0]['object']].data.vertices[0].co.x+=.2
        if case=='background_move':scene.objects[record['environment']['objects'][0]['object']].location.x+=1
        if case=='blur':node('Cytosol far render defocus').inputs['Size'].default_value=(1,1)
        if case=='depth_order':tree.links.new(node('Cellular environment base color').outputs[0],node('Cytosol middle render background strength').inputs['Background'])
        if case=='layer_leak':scene.view_layers['Cytosol far render'].layer_collection.children['Source receptor complexes'].exclude=False
        if case=='light':scene.objects[record['style']['lighting'][0]['name']].data.energy+=10000
        if case=='camera':scene.camera.data.ortho_scale+=10
        if case=='font':
            text=next(o for o in scene.objects if o.type=='FONT');text.data.font=text.data.font_bold
        if case=='text':
            output=next(n for n in tree.nodes if n.type=='GROUP_OUTPUT' and n.is_active_output)
            final=output.inputs['Image'].links[0].from_node
            tree.links.new(node('Cytosol far render defocus').outputs['Image'],final.inputs['Foreground'])
        bpy.context.view_layer.update()
        try:check_appearance(scene,record)
        except ValueError as error:
            assert diagnostic in str(error),(case,str(error));rows.append({'case':case,'rejected':True,'diagnostic':str(error)})
        else:raise AssertionError('Accepted '+case)
    bpy.ops.wm.open_mainfile(filepath=str(path));restored=audit(a.project)
    assert baseline==restored and hashlib.sha256(path.read_bytes()).hexdigest()==original
    a.output.write_text(json.dumps({'passed':True,'cases':rows,'restored_positive':True,'saved_file_unchanged':True},indent=2)+'\n')
    print(json.dumps({'passed':True,'negative_cases':len(rows),'saved_file_unchanged':True}))


if __name__=='__main__':main()
