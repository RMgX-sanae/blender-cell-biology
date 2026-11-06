"""Export complete mesh geometry, place a voxel layout, or audit its saved scene.

Run in Blender. Single mesh snapshots only; realize compound/instanced geometry
and match viewport/render modifier settings before export. Lengths use micrometers.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import subprocess
import tempfile

import bpy
import numpy as np
from mathutils import Matrix

sys.path.insert(0,str(Path(__file__).resolve().parent))
from volume_distribution import ball_inside,enclosing_sphere


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def uniform_matrix(matrix):
    matrix=np.asarray(matrix,float)
    if not np.isfinite(matrix).all():raise ValueError('Expected a finite display matrix')
    scales=np.linalg.svd(matrix[:3,:3],compute_uv=False)
    if (np.linalg.det(matrix[:3,:3])<=0 or
        not np.allclose(scales,scales[0],atol=1e-7,rtol=2e-6)):
        raise ValueError('Expected one positive uniform scale and a proper rotation')
    return matrix


def evaluated_mesh(obj):
    if obj.type!='MESH' or obj.instance_type!='NONE' or obj.children or obj.particle_systems:
        raise ValueError('Realize the complete body as a single mesh before export or checking')
    for modifier in obj.modifiers:
        if modifier.type=='NODES' or modifier.show_viewport!=modifier.show_render:
            raise ValueError('Realize geometry nodes and match viewport/render modifiers first')
        if modifier.type in {'SUBSURF','MULTIRES'} and modifier.levels!=modifier.render_levels:
            raise ValueError('Viewport and render subdivision levels differ')
        if modifier.type=='SCREW' and modifier.steps!=modifier.render_steps:
            raise ValueError('Viewport and render screw steps differ')
    if getattr(obj.cycles,'use_adaptive_subdivision',False):raise ValueError('Bake adaptive subdivision first')
    for slot in obj.material_slots:
        mat=slot.material
        if mat and mat.use_nodes and any(n.type=='OUTPUT_MATERIAL' and n.inputs['Displacement'].is_linked for n in mat.node_tree.nodes):
            raise ValueError('Bake material displacement into the mesh first')
    graph=bpy.context.evaluated_depsgraph_get();evaluated=obj.evaluated_get(graph);mesh=evaluated.to_mesh()
    try:
        vertices=np.empty(len(mesh.vertices)*3,np.float32);mesh.vertices.foreach_get('co',vertices)
        mesh.calc_loop_triangles();faces=np.asarray([tuple(t.vertices) for t in mesh.loop_triangles],np.int32)
        if not len(vertices) or not len(faces):raise ValueError('Expected complete mesh vertices and faces')
        return vertices.reshape(-1,3),faces,np.asarray(evaluated.matrix_world,dtype=float)
    finally:evaluated.to_mesh_clear()


def export_body(obj, output, micrometers_per_unit):
    output=Path(output);source=Path(bpy.data.filepath)
    if output.exists():raise ValueError('Use a fresh body output directory')
    if not source.is_file():raise ValueError('Save the source model before exporting its snapshot')
    if not np.isfinite(micrometers_per_unit) or micrometers_per_unit<=0:raise ValueError('Declare positive micrometers per source unit')
    vertices,faces,matrix=evaluated_mesh(obj);uniform_matrix(matrix)
    source_hash=digest(source)
    # Blender's dirty flag misses some Python edits. Reopen the file in another
    # process and compare the actual snapshot instead of trusting that flag.
    with tempfile.TemporaryDirectory(prefix='volume-source-check-') as directory:
        directory=Path(directory);probe=directory/'probe.py';result=directory/'snapshot.npz'
        probe.write_text("import sys, bpy, numpy as np\nfrom pathlib import Path\n"
                         "a=sys.argv[sys.argv.index('--')+1:]\nsys.path.insert(0,a[0])\n"
                         "from blender_volume import evaluated_mesh\nbpy.context.scene.frame_set(int(a[2]))\n"
                         "v,f,m=evaluated_mesh(bpy.context.scene.objects[a[1]])\n"
                         "np.savez(a[3],vertices=v,faces=f,matrix=m)\n")
        process=subprocess.run([bpy.app.binary_path,'--disable-autoexec','-b',str(source),'--python-exit-code','1',
                                '--python',str(probe),'--',str(Path(__file__).resolve().parent),obj.name,
                                str(bpy.context.scene.frame_current),str(result)],capture_output=True,text=True)
        if process.returncode or not result.is_file():raise ValueError('Could not evaluate the saved source model: '+process.stderr[-1000:])
        with np.load(result,allow_pickle=False) as saved:
            if (digest(source)!=source_hash or not np.array_equal(vertices,saved['vertices']) or
                not np.array_equal(faces,saved['faces']) or not np.allclose(matrix,saved['matrix'],atol=1e-8,rtol=0)):
                raise ValueError('Current geometry differs from the saved source; save the intended model before export')
    vertices=((vertices@matrix[:3,:3].T+matrix[:3,3])*micrometers_per_unit).astype(np.float32)
    enclosing_sphere(vertices)
    output.mkdir(parents=True);np.savez_compressed(output/'body.npz',vertices=vertices,faces=faces)
    record={'kind':'evaluated-volume-body','version':1,'units':'micrometer',
            'source_blend_sha256':source_hash,'source_geometry_reopened_and_compared':True,
            'source_object':obj.name,'source_frame':bpy.context.scene.frame_current,
            'micrometers_per_source_unit':micrometers_per_unit,'blender':bpy.app.version_string,
            'vertices':len(vertices),'triangles':len(faces),'data_sha256':digest(output/'body.npz'),
            'limits':'Static evaluated triangle mesh snapshot at float32 precision; materials, animation and other objects are not exported.'}
    (output/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    return record


def load_layout(path):
    path=Path(path);record=json.loads((path/'manifest.json').read_text())
    if record.get('kind')!='complete-body-volume-layout' or record.get('version')!=1 or record.get('units')!='micrometer':
        raise ValueError('Expected a complete-body micrometer layout')
    names=['source-masks.npz','body.npz','layout.npz']
    if set(record['files'])!=set(names) or any(digest(path/name)!=record['files'][name] for name in names):
        raise ValueError('Prepared volume files differ from the layout manifest')
    if (record['files']['body.npz']!=record['body']['data_sha256'] or
        record['files']['source-masks.npz']!=record['region']['data_sha256']['source-masks.npz'] or
        record['body'].get('source_geometry_reopened_and_compared') is not True):
        raise ValueError('Prepared volume files do not retain their qualified input binding')
    data={name:dict(np.load(path/name,allow_pickle=False)) for name in names}
    return record,data


def place(path, anchor, collection_name='Volume layout'):
    scene=bpy.context.scene;record,data=load_layout(path)
    if 'volume_layout' in scene or collection_name in bpy.data.collections:raise ValueError('Use a scene without an existing volume layout')
    if anchor.type!='EMPTY' or anchor.instance_type!='NONE':raise ValueError('Use an empty representing the shared source coordinate frame')
    bpy.context.view_layer.update();common=uniform_matrix(anchor.matrix_world)
    group=bpy.data.collections.new(collection_name);scene.collection.children.link(group)
    body=data['body.npz'];mesh=bpy.data.meshes.new('Complete volume body')
    mesh.from_pydata(body['vertices'].tolist(),[],body['faces'].tolist());mesh.update()
    objects=[]
    for index,matrix in enumerate(data['layout.npz']['transforms']):
        obj=bpy.data.objects.new(f'Placed body {index+1:03d}',mesh);group.objects.link(obj)
        obj.parent=anchor;obj.matrix_parent_inverse=Matrix.Identity(4);obj.matrix_basis=Matrix(matrix)
        obj['volume_layout_index']=index;objects.append(obj.name)
    scene['volume_layout']=json.dumps({'manifest_sha256':digest(Path(path)/'manifest.json'),
        'anchor':anchor.name,'common_matrix':common.tolist(),'objects':objects,'collection':group.name,
        'frame':scene.frame_current,'limits':record['limits']})
    bpy.context.view_layer.update()
    return [scene.objects[name] for name in objects]


def audit(path):
    scene=bpy.context.scene;record,data=load_layout(path);binding=json.loads(scene['volume_layout'])
    if digest(Path(path)/'manifest.json')!=binding['manifest_sha256']:raise ValueError('Scene and prepared manifest differ')
    if scene.frame_current!=binding['frame']:raise ValueError('Audit the declared static snapshot frame')
    anchor=scene.objects[binding['anchor']];common=uniform_matrix(anchor.matrix_world)
    if not np.allclose(common,binding['common_matrix'],rtol=0,atol=2e-6):raise ValueError('Shared source frame moved after binding')
    inverse=np.linalg.inv(common);group=scene.collection.children[binding['collection']]
    if group.hide_render or set(o.name for o in group.objects)!=set(binding['objects']):raise ValueError('Layout collection hidden or membership changed')
    body=data['body.npz'];region=data['source-masks.npz'];matrices=data['layout.npz']['transforms'];spheres=[]
    for index,name in enumerate(binding['objects']):
        obj=scene.objects[name]
        if obj.hide_render or not obj.visible_get() or obj.parent!=anchor or obj.get('volume_layout_index')!=index:
            raise ValueError('Placed body hidden, detached or reidentified')
        vertices,faces,world=evaluated_mesh(obj)
        if not np.array_equal(vertices,body['vertices']) or not np.array_equal(faces,body['faces']):
            raise ValueError('Complete evaluated body geometry changed')
        actual=inverse@world
        if not np.allclose(actual,matrices[index],atol=3e-6,rtol=2e-6):raise ValueError('Placed body transform changed')
        # Bound the ACTUAL float-precision saved geometry, not the intended poses.
        local_center,_=enclosing_sphere(body['vertices'])
        center=actual[:3,:3]@local_center+actual[:3,3]
        radius=float(np.linalg.norm((vertices-local_center)@actual[:3,:3].T,axis=1).max())
        if not ball_inside(region['available'],region['spacing_xyz_um'],region['origin_xyz_um'],
                           center,radius+record['boundary_gap_um']):raise ValueError('Saved body intersects forbidden region')
        if any(np.linalg.norm(center-c)<=radius+r+record['pair_gap_um'] for c,r in spheres):
            raise ValueError('Saved whole-body bounds overlap')
        spheres.append((center,radius))
    if len(spheres)!=record['count'] or len(spheres)!=len(matrices):raise ValueError('Placed body count differs')
    return {'passed':True,'count':len(spheres),'vertices_per_body':len(body['vertices']),
            'triangles_per_body':len(body['faces']),'actual_saved_geometry_checked':True,
            'limits':'Declared static mesh snapshots, shared source frame and voxel containment only; context geometry provenance, camera, materials, biological distribution and motion require separate checks.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);commands=p.add_subparsers(dest='action',required=True)
    export=commands.add_parser('export');export.add_argument('--object',required=True)
    export.add_argument('--micrometers-per-unit',type=float,required=True);export.add_argument('--output',type=Path,required=True)
    add=commands.add_parser('place');add.add_argument('--layout',type=Path,required=True)
    add.add_argument('--anchor',required=True);add.add_argument('--output',type=Path,required=True)
    check=commands.add_parser('audit');check.add_argument('--layout',type=Path,required=True);check.add_argument('--output',type=Path,required=True)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.output.exists():p.error('Use a fresh output path')
    if a.action=='export':result=export_body(bpy.context.scene.objects[a.object],a.output,a.micrometers_per_unit)
    elif a.action=='place':
        objects=place(a.layout,bpy.context.scene.objects[a.anchor]);result=audit(a.layout)
        a.output.parent.mkdir(parents=True,exist_ok=True);bpy.ops.wm.save_as_mainfile(filepath=str(a.output))
    else:
        result=audit(a.layout);a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
