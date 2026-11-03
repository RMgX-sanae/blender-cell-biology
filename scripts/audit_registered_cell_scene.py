"""Check saved cell surface identities and their shared uniform placement in Blender."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Vector

sys.path.insert(0,str(Path(__file__).resolve().parent))
from prepare_allen_environment import load_cell_bundle


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    scene=bpy.context.scene;record=json.loads(scene['registered_source_cell'])
    if hashlib.sha256((a.asset_dir/'manifest.json').read_bytes()).hexdigest()!=record['asset_manifest_sha256']:
        raise ValueError('Source manifest differs from the saved scene')
    data,sources=load_cell_bundle(a.asset_dir,record['cell'])
    if sources!=record['source_assets']:raise ValueError('Saved source identity differs')
    expected=np.asarray(record['common_matrix']);basis=expected[:3,:3]
    scales=np.linalg.svd(basis,compute_uv=False)
    if np.linalg.det(basis)<=0 or not np.allclose(scales,scales[0],rtol=2e-6,atol=1e-8):
        raise ValueError('Expected a proper rotation with one common positive scale')
    reports={}
    for role,source in data.items():
        obj=scene.objects[record['objects'][role]]
        if obj.type!='MESH' or obj.modifiers or obj.hide_render or obj.data.shape_keys:
            raise ValueError('Source mesh hidden or deformed: '+role)
        if obj.get('source_sha256')!=sources[role]['sha256']:
            raise ValueError('Wrong object source identity: '+role)
        if not np.allclose(np.asarray(obj.matrix_world),expected,atol=1e-7,rtol=1e-7):
            raise ValueError('Registered components moved independently: '+role)
        coords=np.empty(len(obj.data.vertices)*3,np.float32);obj.data.vertices.foreach_get('co',coords)
        coords=coords.reshape(-1,3)
        if not np.array_equal(coords,source['vertices']):raise ValueError('Source vertex coordinates changed: '+role)
        faces=np.array([list(poly.vertices) for poly in obj.data.polygons])
        if not np.array_equal(faces,source['faces']):raise ValueError('Source faces changed: '+role)
        projected=np.array([world_to_camera_view(scene,scene.camera,obj.matrix_world@Vector(v)) for v in coords])
        if np.any(projected[:,:2]<.01) or np.any(projected[:,:2]>.99) or np.any(projected[:,2]<=0):
            raise ValueError('Source geometry leaves the camera frame: '+role)
        reports[role]={'vertices':len(coords),'triangles':len(faces),'projected_bounds':[
            projected[:,:2].min(0).tolist(),projected[:,:2].max(0).tolist()]}
    labels=[o for o in scene.objects if o.type=='FONT']
    if len(labels)!=4:raise ValueError('Expected one title and three labels')
    for obj in labels:
        for key,font in [('zh',obj.data.font),('en',obj.data.font_bold)]:
            if not font.packed_file or hashlib.sha256(font.packed_file.data).hexdigest()!=record['font_sha256'][key]:
                raise ValueError('Saved font differs: '+key)
        if any(form.use_bold!=(ord(char)<128) for char,form in zip(obj.data.body,obj.data.body_format)):
            raise ValueError('Chinese/English font routing changed')
    result={'passed':True,'cell':record['cell'],'source_assets':sources,'objects':reports,
        'common_uniform_scale':float(scales[0]),'fonts_checked':len(labels),
        'limits':'Source surface/face identity, common uniform placement, camera bounds and packed fonts only; no new segmentation, topology, contact, physiological distribution or visual legibility proof.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':True,'cell':record['cell'],'objects':len(reports)}))


if __name__=='__main__':main()
