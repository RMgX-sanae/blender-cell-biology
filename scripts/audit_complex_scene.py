"""Check a saved static complex's source geometry and active artistic appearance."""
import argparse,hashlib,json,math,sys
from pathlib import Path
import bpy
import numpy as np
from mathutils import Quaternion,Vector

sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_surface_materials import check_materials
from detail_scene import linear


def audit(assets,profile):
    scene=bpy.context.scene;source=json.loads((assets/'manifest.json').read_text())
    record=json.loads(scene['complex_display']);components=record['components']
    if record['asset_manifest_sha256']!=hashlib.sha256((assets/'manifest.json').read_bytes()).hexdigest():
        raise ValueError('Source manifest changed')
    expected={row['chain']:row for row in source['chains']}
    if set(components)!=set(expected):raise ValueError('Missing source component')
    names={v['object'] for v in components.values()}
    if names!={o.name for o in scene.objects if o.type=='MESH'}:raise ValueError('Scene component inventory differs')
    check_materials(profile,{'protein','nucleic_acid'})
    bounds=[];rows=[]
    for chain,row in components.items():
        obj=scene.objects[row['object']];data=np.load(assets/f'chain_{chain}.npz')
        if not np.allclose(np.array(obj.matrix_world),np.eye(4),atol=1e-7,rtol=0):
            raise ValueError('Source component moved or stretched: '+chain)
        role='nucleic_acid' if expected[chain]['nucleic_acid'] else 'protein'
        if len(obj.data.materials)!=1 or not obj.data.materials[0] or obj.data.materials[0].get('surface_role')!=role:
            raise ValueError('Wrong active role material: '+chain)
        mat=obj.data.materials[0];shader=mat.node_tree.nodes['Principled BSDF']
        color=(*linear(profile['palette'][row['palette_key']]),1)
        if (shader.inputs['Base Color'].is_linked or not np.allclose(shader.inputs['Base Color'].default_value,color,atol=1e-7,rtol=1e-6)
                or any(p.material_index!=0 for p in obj.data.polygons)):
            raise ValueError('Wrong active identity color: '+chain)
        xyz=np.empty((len(obj.data.vertices),3),np.float32);obj.data.vertices.foreach_get('co',xyz.ravel())
        if row['representation']=='surface':
            loops=np.empty(len(obj.data.loops),np.int32);obj.data.loops.foreach_get('vertex_index',loops)
            if not np.array_equal(xyz,data['vertices']) or not np.array_equal(loops,data['faces'].ravel()):
                raise ValueError('Source surface changed: '+chain)
            if (len(obj.modifiers)!=1 or obj.modifiers[0].type!='SMOOTH' or abs(obj.modifiers[0].factor-.35)>1e-6
                    or obj.modifiers[0].iterations!=2 or not obj.modifiers[0].show_render or not obj.modifiers[0].show_viewport):
                raise ValueError('Source envelope modifiers changed: '+chain)
        elif row['representation']=='atoms':
            # The existing atom display stores 62 vertices per sphere before
            # the inferred bond tubes. Symmetric sphere means recover centers.
            n=len(data['atoms']);atom_vertices=xyz[:n*62].reshape(n,62,3)
            if obj.modifiers or not np.allclose(atom_vertices.mean(1),data['atoms'],atol=1e-5,rtol=0):
                raise ValueError('Source atom centers changed: '+chain)
            radii={'C':.105,'N':.103,'O':.10,'P':.14,'S':.13}
            radius=np.array([radii.get(e,.105) for e in data['elements']])
            if not np.allclose(np.linalg.norm(atom_vertices-data['atoms'][:,None,:],axis=2),radius[:,None],atol=2e-6,rtol=0):
                raise ValueError('Source atom display shape changed: '+chain)
        else:raise ValueError('Unsupported complex representation')
        bounds.extend(camera_point for camera_point in (scene.camera.matrix_world.inverted()@(obj.matrix_world@Vector(p)) for p in obj.bound_box))
        rows.append({'chain':chain,'role':role,'representation':row['representation'],'hidden_complete_subunit':bool(obj.hide_render)})
    camera=scene.camera;points=np.array(bounds);aspect=scene.render.resolution_x/scene.render.resolution_y
    if camera.data.type!='ORTHO' or np.any(points[:,2]>=-camera.data.clip_start):raise ValueError('Invalid complex camera')
    settings=profile.get('camera',{})
    if {'position','target','roll_degrees'}.issubset(settings):
        rotation=(Vector(settings['target'])-Vector(settings['position'])).to_track_quat('-Z','Y')
        rotation=rotation@Quaternion((0,0,1),math.radians(settings['roll_degrees']))
        if abs(rotation.dot(camera.rotation_euler.to_quaternion()))<1-1e-7:raise ValueError('Explicit camera roll ignored')
    if max(np.abs(points[:,0]).max(),np.abs(points[:,1]).max()*aspect)>camera.data.ortho_scale*.42+1e-4:
        raise ValueError('Camera clips the source complex or required margin')
    if 'scale' in settings and not np.isclose(camera.data.ortho_scale,settings['scale'],rtol=1e-6,atol=1e-6):
        raise ValueError('Explicit camera scale ignored')
    if 'position' in settings and not np.allclose(camera.location,settings['position'],atol=1e-6,rtol=0):
        raise ValueError('Explicit camera position ignored')
    if 'target' in settings:
        direction=(Vector(settings['target'])-camera.location).normalized()
        if direction.dot(camera.rotation_euler.to_quaternion()@Vector((0,0,-1)))<1-1e-6:raise ValueError('Explicit camera aim ignored')
    for light in profile['lighting']:
        obj=scene.objects.get(light['name'])
        if (obj is None or obj.type!='LIGHT' or obj.data.type!='AREA'
                or not np.allclose(obj.location,np.array(light['position'])*1.2,atol=1e-6,rtol=0)
                or not np.isclose(obj.data.energy,light['power']*1.2**2)
                or not np.isclose(obj.data.size,light['size']*1.2)
                or not np.allclose(obj.data.color,linear(light['color']),atol=1e-7,rtol=1e-6)):
            raise ValueError('Requested light differs: '+light['name'])
    if (scene.view_settings.view_transform!=profile.get('view_transform','AgX')
            or not np.isclose(scene.view_settings.exposure,profile.get('exposure',0))):raise ValueError('Display transform changed')
    return {'passed':True,'components':rows,'camera_scale':float(camera.data.ortho_scale),
            'limits':'Retained surface meshes and atom centers/radii, active role shaders, identity colors, requested lights and framing. Inferred bond-tube mesh, full chemical assembly, optical calibration and visual quality require separate review.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--asset-dir',type=Path,required=True)
    p.add_argument('--profile',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);result=audit(a.asset_dir,json.loads(a.profile.read_text()))
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':True,'components':len(result['components']),'camera_scale':result['camera_scale']}))


if __name__=='__main__':main()
