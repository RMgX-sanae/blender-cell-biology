"""Mutate a real saved complex scene, reject drift, and restore it in memory."""
import argparse,hashlib,json,sys
from pathlib import Path
import bpy

sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_complex_scene import audit


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--asset-dir',type=Path,required=True)
    p.add_argument('--profile',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);path=Path(bpy.data.filepath)
    original=hashlib.sha256(path.read_bytes()).hexdigest();profile=json.loads(a.profile.read_text())
    baseline=audit(a.asset_dir,profile);rows=[]
    cases={'role':'Wrong active role','color':'Wrong active identity color','surface':'Source surface changed',
           'atom':'Source atom centers changed','camera':'Camera clips','light':'Requested light differs',
           'inventory':'Scene component inventory differs','modifier':'Source envelope modifiers changed'}
    if {'position','target','roll_degrees'}.issubset(profile.get('camera',{})):cases['roll']='Explicit camera roll ignored'
    for case,diagnostic in cases.items():
        bpy.ops.wm.open_mainfile(filepath=str(path))
        components=json.loads(bpy.context.scene['complex_display'])['components']
        surface=next(bpy.data.objects[r['object']] for r in components.values() if r['surface_role']=='protein')
        atom=next(bpy.data.objects[r['object']] for r in components.values() if r['representation']=='atoms')
        if case=='role':surface.data.materials[0]=atom.data.materials[0]
        if case=='color':surface.data.materials[0].node_tree.nodes['Principled BSDF'].inputs['Base Color'].default_value[0]+=.2
        if case=='surface':surface.data.vertices[0].co.x+=.2
        if case=='atom':atom.data.vertices[0].co.x+=.2
        if case=='camera':bpy.context.scene.camera.data.ortho_scale=.1
        if case=='light':bpy.data.objects[profile['lighting'][0]['name']].data.energy+=100
        if case=='inventory':bpy.data.objects.remove(surface,do_unlink=True)
        if case=='modifier':surface.modifiers[0].factor=.8
        if case=='roll':bpy.context.scene.camera.rotation_euler.rotate_axis('Z',.1)
        bpy.context.view_layer.update()
        try:audit(a.asset_dir,profile)
        except ValueError as error:
            assert diagnostic in str(error),(case,str(error));rows.append({'case':case,'rejected':True,'diagnostic':str(error)})
        else:raise AssertionError('Accepted '+case)
    bpy.ops.wm.open_mainfile(filepath=str(path));restored=audit(a.asset_dir,profile)
    assert baseline==restored and hashlib.sha256(path.read_bytes()).hexdigest()==original
    a.output.write_text(json.dumps({'passed':True,'cases':rows,'restored_positive':True,'saved_file_unchanged':True},indent=2)+'\n')
    print(json.dumps({'passed':True,'negative_cases':len(rows),'saved_file_unchanged':True}))


if __name__=='__main__':main()
