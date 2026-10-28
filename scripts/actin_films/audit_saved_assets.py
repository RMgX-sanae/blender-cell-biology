"""Check intact source envelopes, embedded fonts, and saved compositor choices."""
import json
from pathlib import Path
import sys
import bpy
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from project_context import ROOT, SKILL, FFMPEG, FONT, LATIN_FONT, OPTIONS
results=[]
for folder in sorted((ROOT/'shots').iterdir()):
    if not (folder/'scene.blend').exists():continue
    bpy.ops.wm.open_mainfile(filepath=str(folder/'scene.blend'))
    scene=bpy.context.scene;meshes=[];labels=[]
    for obj in scene.objects:
        if obj.type=='FONT':
            labels.append({'text':obj.data.body,'chinese_font':obj.data.font.filepath,
                           'latin_font':obj.data.font_bold.filepath,
                           'embedded':bool(obj.data.font.packed_file and obj.data.font_bold.packed_file),
                           'mixed_selection_ok':all(fmt.use_bold==(ord(char)<0x400 or 0x2010<=ord(char)<=0x22ff) for char,fmt in zip(obj.data.body,obj.data.body_format))})
        if not obj.get('atoms_file'):continue
        d=np.load(obj['atoms_file']);vertices=d['vertices'] if 'vertices' in d else d['f_vertices'];faces=d['faces']
        actual=np.empty(len(obj.data.vertices)*3,dtype=np.float32);obj.data.vertices.foreach_get('co',actual)
        actual=actual.reshape(-1,3)
        f=np.empty(len(obj.data.polygons)*3,dtype=np.int32);obj.data.polygons.foreach_get('vertices',f)
        whole=(actual.shape==vertices.shape and np.allclose(actual,vertices,atol=2e-5) and np.array_equal(f.reshape(-1,3),faces))
        meshes.append({'object':obj.name,'source':obj['atoms_file'],'intact_basis_mesh':bool(whole),'shape_keys':[] if not obj.data.shape_keys else [k.name for k in obj.data.shape_keys.key_blocks],
                       'no_boolean_cut':not any(m.type=='BOOLEAN' for m in obj.modifiers)})
    rec={'shot':folder.name,'meshes':meshes,'labels':labels,'compositors':[g.name for g in bpy.data.node_groups if g.name.startswith('COMP_')],
         'resolution':[scene.render.resolution_x,scene.render.resolution_y],'fps':scene.render.fps}
    results.append(rec)
    print('ASSETS_CHECKED',folder.name,len(meshes),len(labels),flush=True)
(ROOT/'validation/saved-assets.json').write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n')
assert all(m['intact_basis_mesh'] and m['no_boolean_cut'] for r in results for m in r['meshes'])
assert all(x['embedded'] and x['mixed_selection_ok'] for r in results for x in r['labels'])
print('INTACT_MESHES_AND_FONTS_PASS',len(results),flush=True)
