"""Build editable, variable-topology membrane motion from the prepared teaching interval."""
import json
from hashlib import sha256
from pathlib import Path
from types import SimpleNamespace

import bpy
import numpy as np
from mathutils import Vector

import detail_scene as ds
from apply_scale_environment import surface
from membrane_field_nodes import membrane_field


def read_interval(directory):
    record = json.loads((directory/'motion.json').read_text())
    if (record.get('kind') != 'vesicle-topology-motion' or record.get('version') != 1
            or record.get('units') != 'dimensionless authored schematic units'
            or sha256((directory/'source-scope.json').read_bytes()).hexdigest() != record['source_scope_sha256']):
        raise ValueError('Use a prepared motion interval with its matching source scope')
    ids, rows = record['cargo_ids'], record['samples']
    if (len(ids) != 3 or len(set(ids)) != len(ids) or set(ids) & {'ligand','receptor','Membrane'}
            or any(row['frame'] != frame or len(row['cargo']) != len(ids) for frame,row in enumerate(rows,1))):
        raise ValueError('Every playback frame needs the same three identified cargo probes')
    if record['event'] not in {'scission','fusion','uptake'} or not 1 < record['boundary_frame'] <= len(rows):
        raise ValueError('Choose a prepared event with a connection change inside the playback interval')
    if record['event']=='uptake':
        for name,key,rises in [('scission','neck_connected',False),('fusion','pore_connected',True)]:
            changes=[r['frame'] for i,r in enumerate(rows[1:],1) if r[key]!=rows[i-1][key] and r[key]==rises]
            if changes!=[record['boundaries'][name]]:raise ValueError('Unexpected '+name+' boundary')
    return record


def shell_material(name, color):
    mat = surface(name,color,bump=0,settings={'roughness':.39},emission=.13)
    tree = mat.node_tree
    shader = tree.nodes.get('Principled BSDF')
    fresnel = tree.nodes.new('ShaderNodeFresnel')
    fresnel.inputs['IOR'].default_value = 1.38
    opacity = tree.nodes.new('ShaderNodeValToRGB')
    opacity.label = 'Clear face, visible rim'
    opacity.color_ramp.elements[0].color = (.078,.078,.078,1)
    opacity.color_ramp.elements[1].color = (.82,.82,.82,1)
    tree.links.new(fresnel.outputs['Fac'],opacity.inputs['Fac'])
    clear = tree.nodes.new('ShaderNodeBsdfTransparent')
    mix = tree.nodes.new('ShaderNodeMixShader')
    tree.links.new(opacity.outputs['Color'],mix.inputs[0])
    tree.links.new(clear.outputs[0],mix.inputs[1])
    tree.links.new(shader.outputs[0],mix.inputs[2])
    tree.links.new(mix.outputs[0],tree.nodes.get('Material Output').inputs['Surface'])
    return mat


def build(directory,output,width,samples,options=None):
    if output.exists():
        raise ValueError('Use a fresh scene output file')
    record = read_interval(directory)
    rows = record['samples']
    scene = bpy.context.scene
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)
    group = ds.collection('Membrane and tracked cargo')
    scission = record['event']=='scission'
    bounds = ((-4.5,-4.5,-3.75),(4.5,4.5,1.5)) if scission else ((-2.25,-2.5,-8.5),(5.875,2.5,-3.5))
    full = record['event']=='uptake'
    if full:
        if options is None or any(getattr(options,k,None) is None for k in ('coat','ap2','machinery')):
            raise ValueError('Full uptake needs qualified --coat, --ap2 and --machinery bundles')
        bounds=((-4.5,-3.25,-8.75),(5.875,3.25,1.5))
    # Skim the donor sheet at fission. Orbit about the fusion axis to separate
    # the two depth-offset tracers without foreshortening the pore.
    view = ({'position':(9.4,-10.1,-1.28),'target':(0,0,-1.23),'scale':6.6} if scission else
            {'position':(1.8,-10.,5.4),'target':(1.8,0,-6.),'scale':8.2})
    profile = json.loads((Path(__file__).resolve().parents[1]/'assets/style-organelle.json').read_text())
    profile.update(background='#071B27',camera={},exposure=0)
    args = SimpleNamespace(engine='cycles',samples=samples,width=width,fps=record['fps'],frames=len(rows),compositor='clean')
    ds.studio(profile,args,{**view,'light_scale':.46})
    scene.cycles.transparent_max_bounces=24

    membrane = bpy.data.objects.new('Membrane',bpy.data.meshes.new('Native membrane field'))
    group.objects.link(membrane)
    mats = [shell_material('Cytosolic membrane','#428E99'),shell_material('Lumen-facing membrane','#87DBCA')]
    for mat in mats:
        membrane.data.materials.append(mat)
    modifier, controls = membrane_field(membrane,record,bounds,mats)
    modifier.show_viewport = False
    ids = record['cargo_ids']
    actors = {}
    for key,color in zip([*ids,'ligand'],['#FFE36A','#FFAE64','#FF806C','#ED9BDD']):
        parts = ds.Parts()
        parts.ball((0,0,0),.08 if key=='ligand' else record['cargo_radius'],seg=24,rings=16)
        actors[key] = parts.object(key,[surface(key,color,bump=0,emission=.50)],group)
        actors[key]['topology_role'] = key
    parts = ds.Parts()
    parts.tube([(0,0,-.17),(0,0,.28)],.025,sides=24)
    receptor = parts.object('receptor',[surface('receptor','#B7EDA7',bump=0,emission=.3)],group)
    receptor['topology_role'] = 'receptor'
    receptor.rotation_mode = 'QUATERNION'
    membrane['topology_role'] = 'membrane'
    channels = {'Center':'center','Radius':'radius','Neck radius':'neck_radius','Neck open':'neck_connected',
                'Pore radius':'pore_radius','Pore open':'pore_connected'}
    interface = getattr(modifier,'properties',None)
    paths = {name:modifier.path_from_id() + ('.properties.inputs.'+key+'.value' if interface else '["'+key+'"]')
             for name,key in controls.items()}
    for row in rows:
        f = row['frame']
        for name,key in channels.items():
            prop = controls[name]
            value = row[key]
            if interface:
                getattr(interface.inputs,prop).value = value
            else:
                modifier[prop] = tuple(value) if isinstance(value,list) else value
            membrane.keyframe_insert(data_path=paths[name],frame=f,group=name)
        for key,point in zip([*ids,'ligand'],[*row['cargo'],row['ligand']]):
            obj = actors[key]
            obj.location = point
            obj.keyframe_insert(data_path='location',frame=f)
        receptor.location = row['receptor']
        receptor.rotation_quaternion = Vector(row['normal']).to_track_quat('Z','Y')
        receptor.keyframe_insert(data_path='location',frame=f)
        receptor.keyframe_insert(data_path='rotation_quaternion',frame=f)
    switches = {paths[name] for name in ('Neck open','Pore open')}
    for action in [obj.animation_data.action for obj in [membrane,receptor,*actors.values()]]:
        for layer in action.layers:
            for strip in layer.strips:
                for bag in strip.channelbags:
                    for curve in bag.fcurves:
                        mode = 'CONSTANT' if curve.data_path in switches else 'LINEAR'
                        for point in curve.keyframe_points:
                            point.interpolation = mode
    modifier.show_viewport = True
    labels = ('膜内陷','切断','囊泡分离') if scission else ('两腔分离','开孔','货物进入受体腔')
    if full:
        from vesicle_uptake_actors import animate
        animate(scene,record,options)
    else:
        for text,frame in zip(labels,[1,record['boundary_frame'],len(rows)]):
            scene.timeline_markers.new(text,frame=frame)
    scene['vesicle_motion'] = json.dumps({'motion_sha256':sha256((directory/'motion.json').read_bytes()).hexdigest(),'event':record['event'],
        'units':record['units'],'frame_of_connectivity_change':record['boundary_frame'],'limits':record['limits']})
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type=='VIEW_3D':
                space = area.spaces.active
                space.region_3d.view_perspective='CAMERA'
                if screen == bpy.context.screen:
                    space.shading.type='MATERIAL'
    scene.frame_set(1)
    output.parent.mkdir(parents=True,exist_ok=True)
    bpy.ops.wm.save_as_mainfile(filepath=str(output))


def shell_components(mesh):
    import bmesh
    graph = bmesh.new()
    try:
        graph.from_mesh(mesh)
        remaining = set(graph.verts)
        counts = []
        while remaining:
            stack = [remaining.pop()]
            size = 0
            while stack:
                vertex = stack.pop()
                size += 1
                for edge in vertex.link_edges:
                    adjacent = edge.other_vert(vertex)
                    if adjacent in remaining:
                        remaining.remove(adjacent)
                        stack.append(adjacent)
            counts.append(size)
        return sorted(counts)
    finally:
        graph.free()


def audit(directory):
    from mathutils.bvhtree import BVHTree
    record = read_interval(directory)
    scene = bpy.context.scene
    binding = json.loads(scene['vesicle_motion'])
    if binding['motion_sha256'] != sha256((directory/'motion.json').read_bytes()).hexdigest():
        raise ValueError('The saved animation and prepared interval differ')
    boundary = record['boundary_frame']
    if (scene.frame_start,scene.frame_end,scene.render.fps) != (1,len(record['samples']),record['fps']):
        raise ValueError('The stored playback range or frame rate changed')
    membrane = scene.objects['Membrane']
    modifier = next((m for m in membrane.modifiers if m.type=='NODES' and m.show_render and m.show_viewport),None)
    if modifier is None or not np.allclose(membrane.matrix_world,np.eye(4)):
        raise ValueError('The native field must remain enabled in its prepared coordinate frame')
    ids = record['cargo_ids']
    actors = {key:scene.objects[key] for key in [*ids,'ligand','receptor']}
    for key,obj in actors.items():
        if obj.get('topology_role')!=key or obj.hide_render or not obj.visible_get():
            raise ValueError('Tracked actor is missing or hidden: '+key)
        if key!='receptor':
            radius = .08 if key=='ligand' else record['cargo_radius']
            verts = np.array([list(v.co) for v in obj.data.vertices])
            if not np.allclose(np.linalg.norm(verts,axis=1),radius,rtol=0,atol=2e-6):
                raise ValueError('Tracked probe dimensions changed: '+key)
    checks = []
    scission = record['event']=='scission'
    full=record['event']=='uptake'
    frames={1,max(1,boundary-1),boundary,len(record['samples'])}
    if full:
        pore=record['boundaries']['fusion']
        frames.update((pore-1,pore))
    for frame in sorted(frames):
        row = record['samples'][frame-1]
        scene.frame_set(frame)
        graph = bpy.context.evaluated_depsgraph_get()
        mesh = membrane.evaluated_get(graph).data
        components = shell_components(mesh)
        # The cropped donor sheet has one closed rim. Each isolated, complete
        # vesicle has two boundary surfaces, one facing each aqueous side.
        at_neck=scission or (full and frame<=boundary)
        connected = row['neck_connected'] if at_neck else row['pore_connected']
        wanted = (1 if connected else 3) if scission else (2 if connected else 4)
        if full:wanted=3 if row['neck_connected'] or row['pore_connected'] else 5
        if len(components)!=wanted:
            raise ValueError(f'At frame {frame} the native membrane has {len(components)} surface components, expected {wanted}')
        tree = BVHTree.FromObject(membrane,graph)
        center = Vector(row['center'])
        other = Vector(record['neck_end'] if at_neck else record['target_center'])
        line = other-center
        hit = tree.ray_cast(center,line.normalized(),line.length)[0]
        if (hit is None)!=connected:
            raise ValueError(f'The rendered membrane blocks the wrong aqueous passage at frame {frame}')
        margin = float('inf')
        for key,point in zip([*ids,'ligand','receptor'],[*row['cargo'],row['ligand'],row['receptor']]):
            matrix = actors[key].matrix_world
            if not np.allclose(list(matrix.translation),point,rtol=0,atol=2e-6) or not np.allclose(matrix.to_scale(),[1,1,1]):
                raise ValueError(f'{key} moved or stretched at frame {frame}')
            if key=='receptor':
                if not np.allclose(matrix.to_3x3().col[2],row['normal'],rtol=0,atol=2e-6):
                    raise ValueError('The receptor changed membrane side')
                continue
            gap = tree.find_nearest(Vector(point))[3] - (.08 if key=='ligand' else record['cargo_radius'])
            ray = center-Vector(point)
            if gap<=0 or (ray.length>1e-6 and tree.ray_cast(Vector(point),ray.normalized(),ray.length)[0] is not None):
                raise ValueError(f'{key} touches membrane or is separated from the carrier lumen at frame {frame}')
            margin = min(margin,gap)
        checks.append({'frame':frame,'time':row['time'],'surface_components':len(components),
                       'aqueous_passage_open':connected,'minimum_probe_clearance':margin})
    molecular={}
    if full:
        from vesicle_uptake_actors import audit_machinery
        molecular={'molecular_actors':audit_machinery(scene,record),'boundaries':record['boundaries']}
    return {**molecular,'passed':True,'event':record['event'],'playback_frames':scene.frame_end,
            'frame_of_connectivity_change':boundary,'cargo_ids':ids,'native_mesh_checks':checks,
            'limits':record['limits']}
