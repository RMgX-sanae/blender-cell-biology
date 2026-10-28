"""Check actual saved interaction geometry, transforms, framing and loop seam."""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Vector


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    scene=bpy.context.scene;folder=Path(scene['interaction_asset_dir'])
    manifest=json.loads((folder/'manifest.json').read_text());expected=np.load(folder/'motion.npz')['matrices']
    errors=[];meshes=[]
    for actor in manifest['actors']:
        obj=scene.objects.get(actor['id'])
        if obj is None:errors.append('Missing actor '+actor['id']);continue
        data=np.load(folder/actor['file'])
        reference=data['atoms'] if obj.modifiers else data['vertices']
        actual=np.empty(len(obj.data.vertices)*3);obj.data.vertices.foreach_get('co',actual);actual=actual.reshape(-1,3)
        intact=actual.shape==reference.shape and np.allclose(actual,reference,atol=2e-6)
        if not obj.modifiers and 'faces' in data:
            faces=np.array([list(face.vertices) for face in obj.data.polygons])
            intact=bool(intact and np.array_equal(faces,data['faces']))
        meshes.append({'id':actor['id'],'intact_source_geometry':bool(intact),'source_chain':actor['chain']})
        if not intact:errors.append('Geometry changed '+actor['id'])
        if obj.get('biological_id')!=actor['id']:errors.append('Identity changed '+actor['id'])
        if any(m.type=='BOOLEAN' for m in obj.modifiers):errors.append('Cutting modifier on '+actor['id'])
    maximum_error=0.;minimum_margin=1.;clipped=[]
    for frame,pose in enumerate(expected,1):
        scene.frame_set(frame);bpy.context.view_layer.update()
        actual=np.array(scene.objects['moving'].matrix_world)
        maximum_error=max(maximum_error,float(abs(actual-pose).max()))
        for actor in manifest['actors']:
            obj=scene.objects[actor['id']]
            for point in obj.bound_box:
                v=world_to_camera_view(scene,scene.camera,obj.matrix_world@Vector(point))
                margin=min(v.x,v.y,1-v.x,1-v.y);minimum_margin=min(minimum_margin,margin)
                if margin<0 or v.z<=0:clipped.append({'frame':frame,'actor':actor['id']});break
    if maximum_error>2e-5:errors.append('Saved animation differs from prepared poses')
    if clipped:errors.append('Main actor bounding boxes leave the camera frame')
    if not manifest['motion']['collision_screen_passed']:errors.append('Prepared path contact screen failed')
    if scene.frame_end!=manifest['frames'] or scene.render.fps!=manifest['fps']:errors.append('Timing differs from manifest')
    result={'scene':bpy.data.filepath,'actors':meshes,'pose_error_max':maximum_error,
            'camera_minimum_fractional_margin':minimum_margin,'clipped_examples':clipped[:15],
            'checked_frames_including_virtual_seam':len(expected),'prepared_contact_screen':manifest['motion'],
            'errors':errors,'passed':not errors,
            'limits':'Technical and geometric checks only. Inspect encoded video, readability, source suitability and biological scope separately.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2))
    print(json.dumps({'passed':not errors,'errors':errors,'minimum_margin':minimum_margin}))
    if errors:raise SystemExit(1)


if __name__=='__main__':main()
