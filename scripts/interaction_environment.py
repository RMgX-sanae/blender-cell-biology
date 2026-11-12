"""Native-size reviewed cytosolic shapes, kept outside the complete actor path."""
import hashlib,json
from pathlib import Path
import bpy,numpy as np
from mathutils import Euler,Matrix
import detail_scene as ds
from cellular_environment import composite
from interaction_membrane import plane_clearance


def build(asset_dir,settings,appearance,subject_bounds,membrane=None):
    allowed={'compartment','context_note','context_sources','objects','palette','compositor'}
    if set(settings)-allowed:raise ValueError('Unsupported interaction environment field')
    if settings.get('compartment')!='cytosol' or not settings.get('context_note') or not settings.get('context_sources'):
        raise ValueError('These source assets require an explicitly reviewed cytosolic context')
    root=Path(asset_dir);manifest=root/'manifest.json';source=json.loads(manifest.read_text())
    if source.get('kind')!='schematic-cytosolic-environment' or source.get('units')!='nm':raise ValueError('Require prepared native-nm cytosolic assets')
    palette={'ribosome':'#294F5A','actin':'#44506A','ubiquitin':'#4F6C76'};palette.update(settings.get('palette',{}))
    if set(palette)-{'ribosome','actin','ubiquitin'}:raise ValueError('Unknown environment source palette')
    groups=[ds.collection('Interaction context '+name) for name in ('far','middle','near')]
    meshes={};data={};hashes={};records=[]
    for row in settings['objects']:
        if set(row)-{'asset','depth','position_nm','rotation_degrees','repeat_indices'}:raise ValueError('Unsupported background placement field')
        name=row['asset'];depth=row['depth']
        if name not in palette or depth not in ('far','middle','near'):raise ValueError('Unknown background source or depth')
        for key in ('position_nm','rotation_degrees'):
            if np.shape(row.get(key))!=(3,) or not np.isfinite(row[key]).all():raise ValueError('Background transforms require finite triples')
        if name not in data:
            path=root/(name+'.npz');data[name]=dict(np.load(path));hashes[name]=hashlib.sha256(path.read_bytes()).hexdigest()
            mesh=bpy.data.meshes.new('Source context '+name);mesh.from_pydata(data[name]['vertices'],[],data[name]['faces']);mesh.update()
            mat=ds.material('Interaction background '+name,palette[name],.78,0,settings=appearance.get('surface_roles',{}).get('environment'));mat['surface_role']='environment';mesh.materials.append(mat)
            for polygon in mesh.polygons:polygon.use_smooth=True
            meshes[name]=mesh
        repeats=row.get('repeat_indices',[0])
        if not repeats or any(type(i)is not int for i in repeats) or len(set(repeats))!=len(repeats) or name!='actin' and repeats!=[0]:raise ValueError('Only actin supports distinct integer source repeat indices')
        base=np.eye(4);base[:3,:3]=np.array(Euler(np.radians(row['rotation_degrees'])).to_matrix());base[:3,3]=row['position_nm']
        for index in repeats:
            matrix=base@np.linalg.matrix_power(data[name]['repeat_nm'],index) if name=='actin' else base
            points=data[name]['vertices']@matrix[:3,:3].T+matrix[:3,3]
            lo,hi=points.min(0),points.max(0)
            if np.all(hi>=subject_bounds[0]) and np.all(lo<=subject_bounds[1]):
                raise ValueError(f'Background bounds overlap the complete actor path: {name}, placement {row["position_nm"]}, bounds {lo.tolist()} to {hi.tolist()}, subject {np.asarray(subject_bounds).tolist()}')
            clearance=None
            if membrane:
                guard=membrane['head_offset_nm']+membrane['head_radius_nm']-membrane['hydrophobic_half_width_nm']
                clearance=float(plane_clearance(data[name]['vertices'],matrix[None],membrane,guard)[0])
                if clearance<0:raise ValueError(f'Background source crosses the cytosolic membrane boundary: {name}, clearance {clearance:.6g} nm')
            obj=bpy.data.objects.new('Interaction context '+name+' '+str(len(records)),meshes[name]);groups[('far','middle','near').index(depth)].objects.link(obj)
            obj.matrix_world=Matrix(matrix);obj['environment_only']=True;obj['source_structure']=source['sources'][name]['id'];obj['depth_layer']=depth
            records.append({'object':obj.name,'asset':name,'depth':depth,'matrix_nm':matrix.tolist(),'cytosolic_clearance_nm':clearance})
    if not records:raise ValueError('Provide at least one reviewed environment placement')
    return groups,{'asset_dir':str(root.resolve()),'manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest(),
                   'asset_sha256':hashes,'source_manifest':source,'objects':records,'palette':palette,'settings':settings,
                   'scope':'Intact native-size source shapes with authored static positions and separate depth defocus; not measured abundance or diffusion.'}


def compose(scene,groups,settings,background,width):
    beauty=scene.view_layers[scene['compositor_beauty_layer']];layers=[]
    lights={c.name for o in scene.objects if o.type in ('CAMERA','LIGHT') for c in o.users_collection if c.name!=scene['compositor_label_collection']}
    for group in groups:
        beauty.layer_collection.children[group.name].exclude=True
        layer=scene.view_layers.new(group.name+' render')
        for child in layer.layer_collection.children:child.exclude=child.name not in lights|{group.name}
        layers.append(layer)
    tree=scene.compositing_node_group;image=composite(tree,layers,background,width,settings=settings.get('compositor',{}))
    base=next(n for n in tree.nodes if n.label=='Straight background plus premultiplied rendered foreground')
    tree.links.new(image,base.inputs['Background'])


def audit_environment(scene,record,manifest,node,linked,equal,base,background):
    root=Path(record['asset_dir'])
    if hashlib.sha256((root/'manifest.json').read_bytes()).hexdigest()!=record['manifest_sha256']:raise ValueError('Environment source manifest changed')
    if {o.name for o in scene.objects if o.get('environment_only')}!={r['object'] for r in record['objects']}:raise ValueError('Environment object inventory differs')
    data={};original=bpy.context.window.view_layer;previous=node('Cellular environment base color').outputs[0]
    equal(previous.default_value,(*ds.linear(background),1),'Environment background color differs')
    settings=record['settings'].get('compositor',{});beauty=scene.view_layers[scene['compositor_beauty_layer']]
    try:
        for i,depth in enumerate(('far','middle','near')):
            group='Interaction context '+depth;layer=scene.view_layers[group+' render']
            if not beauty.layer_collection.children[group].exclude:raise ValueError('Environment leaks into sharp beauty')
            allowed={group,scene.camera.users_collection[0].name}
            if {c.name for c in layer.layer_collection.children if not c.exclude}!=allowed:raise ValueError('Subject or labels leak into context blur')
            bpy.context.window.view_layer=layer;bpy.context.view_layer.update();graph=bpy.context.evaluated_depsgraph_get()
            for row in record['objects']:
                if row['depth']!=depth:continue
                name=row['asset'];path=root/(name+'.npz')
                if name not in data:
                    if hashlib.sha256(path.read_bytes()).hexdigest()!=record['asset_sha256'][name]:raise ValueError('Environment source mesh changed')
                    data[name]=dict(np.load(path))
                obj=scene.objects[row['object']];reference=data[name]
                if obj.parent or obj.modifiers or obj.constraints or obj.animation_data or obj.hide_render:raise ValueError('Environment whole-object display changed')
                points=np.empty((len(obj.data.vertices),3),np.float32);obj.data.vertices.foreach_get('co',points.ravel())
                if not np.array_equal(points,reference['vertices'].astype(np.float32)) or not np.array_equal([list(p.vertices) for p in obj.data.polygons],reference['faces']):raise ValueError('Environment source geometry changed')
                matrix=np.array(obj.evaluated_get(graph).matrix_world);equal(matrix,row['matrix_nm'],'Environment placement changed')
                if len(obj.material_slots)!=1 or obj.material_slots[0].material.get('surface_role')!='environment':raise ValueError('Environment material assignment differs')
                mat=obj.material_slots[0].material;equal(mat.node_tree.nodes['Principled BSDF'].inputs['Base Color'].default_value,(*ds.linear(record['palette'][name]),1),'Environment color differs')
                if manifest.get('membrane'):
                    membrane=manifest['membrane'];guard=membrane['head_offset_nm']+membrane['head_radius_nm']-membrane['hydrophobic_half_width_nm']
                    clearance=float(plane_clearance(points,matrix[None],membrane,guard)[0])
                    if clearance<0:raise ValueError('Saved environment crosses membrane boundary')
                    equal(clearance,row['cytosolic_clearance_nm'],'Saved environment clearance differs')
            blur=node(layer.name+' defocus');over=node(layer.name+' background strength')
            equal(blur.inputs['Size'].default_value,[settings.get('blur_pixels_1280',[18,9,28])[i]*scene.render.resolution_x/1280]*2,'Context blur width differs')
            equal(over.inputs['Factor'].default_value,settings.get('opacity',[.26,.24,.14])[i],'Context opacity differs')
            if len(blur.inputs['Image'].links)!=1:raise ValueError('Context image disconnected')
            render=blur.inputs['Image'].links[0].from_node
            if render.type!='R_LAYERS' or render.layer!=layer.name:raise ValueError('Context render layer differs')
            linked(over.inputs['Foreground'],blur.outputs['Image'],'Context blur bypassed')
            linked(over.inputs['Background'],previous,'Context depth order changed');previous=over.outputs['Image']
        linked(base.inputs['Background'],previous,'Context layers not used')
    finally:bpy.context.window.view_layer=original
