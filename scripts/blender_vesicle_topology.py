"""Build or audit one editable topology snapshot in Blender; no continuous film."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import bpy
import numpy as np
from mathutils import Vector

sys.path.insert(0,str(Path(__file__).resolve().parent))
import detail_scene as ds
from apply_scale_environment import surface
from blender_volume import evaluated_mesh


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_key(directory,t):
    record=json.loads((directory/'manifest.json').read_text())
    if record.get('kind')!='vesicle-topology-keys' or record.get('version')!=1 or record.get('units')!='dimensionless authored schematic units':
        raise ValueError('Expected prepared vesicle topology keys')
    key=next((k for k in record['snapshots'] if abs(k['time']-t)<1e-9),None)
    if key is None:raise ValueError('Time is not a prepared snapshot')
    for name,expected in [(key['file'],key['sha256']),('audit.json',record['audit_sha256']),('source-scope.json',record['source_scope_sha256'])]:
        if digest(directory/name)!=expected:raise ValueError('Prepared topology files differ from the manifest')
    if not json.loads((directory/'audit.json').read_text()).get('passed'):raise ValueError('Topology checks did not pass')
    with np.load(directory/key['file'],allow_pickle=False) as data:data=dict(data)
    return record,key,data


def build(directory,t,output,width,samples):
    if output.exists():raise ValueError('Use a fresh scene output file')
    record,key,data=load_key(directory,t)
    bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
    group=ds.collection('Topology snapshot');scene=bpy.context.scene
    for role,color in [('lumen_surface','#63d7eb'),('cytosolic_surface','#6ca8d7')]:
        obj=ds.mesh(role,data[role+'_vertices'].tolist(),data[role+'_faces'].tolist(),
                    [surface(role,color,.18,bump=0)],group)
        obj['topology_role']=role
    colors=['#ffe16c','#f5a24b','#f27463']
    for identifier,position,color in zip(record['cargo_ids'],data['cargo'],colors):
        parts=ds.Parts();parts.ball(position,record['cargo_radius'],seg=24,rings=16)
        obj=parts.object(identifier,[surface(identifier,color,bump=0)],group);obj['topology_role']=identifier
    parts=ds.Parts();parts.ball(data['ligand'],.08,seg=24,rings=16)
    obj=parts.object('ligand',[surface('ligand','#ef8ff2',bump=0)],group);obj['topology_role']='ligand'
    # An explicitly idealized single-pass receptor, not a molecular structure.
    normal=Vector(data['normal']);anchor=Vector(data['receptor'])
    bpy.ops.mesh.primitive_cylinder_add(vertices=24,radius=.025,depth=.45,
                                        location=anchor+normal*.055)
    obj=bpy.context.object;obj.name='receptor';obj.rotation_mode='QUATERNION'
    obj.rotation_quaternion=normal.to_track_quat('Z','Y');obj['topology_role']='receptor'
    obj.data.materials.append(surface('receptor','#b5ef8a',bump=0))
    profile=json.loads((Path(__file__).resolve().parents[1]/'assets/style-organelle.json').read_text())
    profile.update(background='#071825',camera={},exposure=0)
    args=SimpleNamespace(engine='cycles',samples=samples,width=width,fps=24,frames=1,compositor='clean')
    view=({'position':(2,-20,-4),'target':(1.8,0,-6),'scale':8.6} if t>=.62 else
          {'position':(7,-22,4),'target':(.3,0,-3.5),'scale':15})
    camera=ds.studio(profile,args,{**view,'light_scale':.45})
    if t<.62:ds.fit_orthographic(camera,[o for o in scene.objects if o.type=='MESH'])
    scene.cycles.transparent_max_bounces=24
    scene['vesicle_topology']=json.dumps({'manifest_sha256':digest(directory/'manifest.json'),
         'time':key['time'],'units':record['units'],'scope':'Static authored topology; molecular machinery and continuous mechanism pending.'})
    output.parent.mkdir(parents=True,exist_ok=True)
    bpy.ops.wm.save_as_mainfile(filepath=str(output))


def audit(directory):
    scene=bpy.context.scene;binding=json.loads(scene['vesicle_topology'])
    if digest(directory/'manifest.json')!=binding['manifest_sha256']:raise ValueError('Saved scene bound to a different prepared case')
    record,key,data=load_key(directory,binding['time'])
    expected={'lumen_surface','cytosolic_surface',*record['cargo_ids'],'ligand','receptor'}
    objects=[obj for obj in scene.objects if obj.type=='MESH']
    roles=[obj.get('topology_role') for obj in objects]
    if len(roles)!=len(expected) or set(roles)!=expected:raise ValueError('Saved topology object inventory differs')
    for obj in objects:
        role=obj['topology_role']
        if obj.hide_render or obj.hide_get() or not obj.visible_get():raise ValueError('Tracked topology object is hidden')
        vertices,faces,matrix=evaluated_mesh(obj)
        if role.endswith('_surface'):
            if not np.array_equal(vertices,data[role+'_vertices']) or not np.array_equal(faces,data[role+'_faces']) or not np.allclose(matrix,np.eye(4),atol=1e-8,rtol=0):
                raise ValueError('Saved membrane surface differs from qualified geometry')
        else:
            points=vertices@matrix[:3,:3].T+matrix[:3,3]
            edges=np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1)
            if np.any(np.unique(edges,axis=0,return_counts=True)[1]!=2) or np.linalg.det(matrix[:3,:3])<=0:
                raise ValueError('Saved probe is open or reflected')
            if role=='receptor':
                offset=points-data['receptor'];axial=offset@data['normal']
                radial=np.linalg.norm(offset-axial[:,None]*data['normal'],axis=1)
                direction=matrix[:3,2]/np.linalg.norm(matrix[:3,2])
                if not np.allclose(direction,data['normal'],atol=2e-6,rtol=0) or not np.allclose([axial.min(),axial.max()],[-.17,.28],atol=2e-6,rtol=0) or not np.allclose(radial,.025,atol=2e-6,rtol=0):
                    raise ValueError('Saved receptor domains changed')
            else:
                position=data['ligand'] if role=='ligand' else data['cargo'][record['cargo_ids'].index(role)]
                radius=.08 if role=='ligand' else record['cargo_radius']
                expected_probe=ds.Parts();expected_probe.ball(position,radius,seg=24,rings=16)
                if points.shape!=np.asarray(expected_probe.v).shape or not np.allclose(points,expected_probe.v,atol=2e-6,rtol=0) or [tuple(p.vertices) for p in obj.data.polygons]!=expected_probe.f:
                    raise ValueError('Saved complete cargo or ligand geometry changed')
    return {'passed':True,'time':key['time'],'objects':len(objects),'saved_geometry_matches_prepared':True,
            'limits':'Static snapshot geometry only; camera/materials and biological completeness are separate checks.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['build','audit'])
    p.add_argument('--input',type=Path,required=True);p.add_argument('--time',type=float,default=.8)
    p.add_argument('--output',type=Path);p.add_argument('--width',type=int,default=1280);p.add_argument('--samples',type=int,default=24)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.command=='build':
        if a.output is None:p.error('build requires --output')
        build(a.input,a.time,a.output,a.width,a.samples)
    result=audit(a.input)
    if a.command=='audit' and a.output:
        if a.output.exists():raise ValueError('Use a fresh audit output file')
        a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
