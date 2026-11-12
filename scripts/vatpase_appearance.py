"""Reference appearance for the existing frontal whole-pump/c-ring layout."""
import hashlib,json,math,re
from pathlib import Path
import bpy,numpy as np
from mathutils import Vector
import detail_scene as ds
from interaction_appearance import build_compositor,studio_profile
import interaction_environment as environment

PALETTE={'A','B','rotor','ring','ring_tracer','a','stalk','lipid','tail','h','atp','in','out','background','proton','ink'}
ROLES={'protein','membrane_heads','membrane_tails','chemical_token','proton_marker','flow_guide','environment'}
MEMBRANE={'origin_nm':[0,0,0],'normal':[0,0,1],'hydrophobic_half_width_nm':1.9,'head_offset_nm':2.05,'head_radius_nm':.25}


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_profile(path):
    p=json.loads(Path(path).read_text())
    allowed={'name','reference','reference_page','usage','observations','implementation_assumptions','palette','surface_roles',
             'background','world','lighting','view_transform','look','exposure','camera','render_recipe','environment'}
    if set(p)-allowed:raise ValueError('Unsupported V-ATPase appearance fields: '+str(set(p)-allowed))
    if set(p.get('palette',{}))-PALETTE:raise ValueError('Unknown V-ATPase semantic palette key')
    if set(p.get('surface_roles',{}))-ROLES:raise ValueError('Unknown V-ATPase surface role')
    colors=list(p.get('palette',{}).values())+[p[k] for k in ('world','background') if k in p]
    if any(not isinstance(c,str) or not re.fullmatch(r'#?[0-9a-fA-F]{6}',c) for c in colors):raise ValueError('Appearance colors must be sRGB hex')
    camera=p.get('camera',{})
    if set(camera)-{'scale_multiplier','center_offset_nm'}:raise ValueError('Authored dual view supports orthographic scale_multiplier/center_offset_nm; other cameras require layout reauthoring')
    scale=camera.get('scale_multiplier',1);pan=camera.get('center_offset_nm',[0,0])
    if not isinstance(scale,(int,float)) or not math.isfinite(scale) or scale<=0:raise ValueError('Camera scale must be finite and positive')
    if np.shape(pan)!=(2,) or not np.isfinite(pan).all():raise ValueError('Camera pan must be a finite x/z pair')
    if 'lighting'in p:
        if not p['lighting'] or len({v['name'] for v in p['lighting']})!=len(p['lighting']):raise ValueError('Provide uniquely named area lights')
        for light in p['lighting']:
            if set(light)!={'name','position','power','color','size'}:raise ValueError('Area lights require name/position/power/color/size')
            if np.shape(light['position'])!=(3,) or not np.isfinite(light['position']).all():raise ValueError('Light position must be finite')
            if any(not math.isfinite(light[k]) or light[k]<=0 for k in ('power','size')):raise ValueError('Light power and size must be positive')
            if not re.fullmatch(r'#?[0-9a-fA-F]{6}',light['color']):raise ValueError('Invalid light color')
    if 'exposure'in p and not math.isfinite(p['exposure']):raise ValueError('Exposure must be finite')
    env=p.get('environment')
    if env is not None and (env.get('compartment')!='cytosol' or not env.get('context_note') or not env.get('context_sources')):
        raise ValueError('Source context requires a reviewed cytosolic environment')
    return p


def shade(profile,palette,key,role,name=None):
    mat=ds.material(name or role+' '+key,palette[key],.5,0,settings=profile.get('surface_roles',{}).get(role))
    mat['surface_role']=role;mat['palette_key']=key
    return mat


def finish(scene,args,appearance,palette,studio):
    # Keep the original geometry and baked flow, including off-screen pool paths.
    points=[]
    bodies=[o for o in scene.objects if o.type in ('MESH','CURVE')]
    for frame in range(1,386):
        scene.frame_set(frame);bpy.context.view_layer.update();graph=bpy.context.evaluated_depsgraph_get()
        for obj in bodies:
            if obj.hide_render:continue
            eo=obj.evaluated_get(graph);points.extend(eo.matrix_world@Vector(p) for p in eo.bound_box)
    scene.frame_set(1);bpy.context.view_layer.update()
    groups=[];env=None
    if appearance.get('environment'):
        cloud=np.array(points)
        groups,env=environment.build(args.environment_dir,appearance['environment'],appearance,(cloud.min(0),cloud.max(0)),MEMBRANE)
        # The right panel is a magnified isolated ring, not another physical pump.
        for row in env['objects']:
            obj=scene.objects[row['object']]
            if max((obj.matrix_world@v.co).x for v in obj.data.vertices)>2:
                raise ValueError('Native-size background must remain in the whole-pump panel (x <= 2 nm)')
    ink=bpy.data.materials.new('Protected V-ATPase ink');ink.use_nodes=True;ink.node_tree.nodes.clear()
    emit=ink.node_tree.nodes.new('ShaderNodeEmission');emit.inputs['Color'].default_value=(*ds.linear(palette['ink']),1)
    out=ink.node_tree.nodes.new('ShaderNodeOutputMaterial');ink.node_tree.links.new(emit.outputs[0],out.inputs['Surface'])
    fonts=[]
    for obj in scene.objects:
        if obj.type!='FONT':continue
        obj.data.materials.clear();obj.data.materials.append(ink);obj['compositor_preserve_text_color']=True
        fonts.append({'object':obj.name,'body':obj.data.body})
    compositor=build_compositor(scene,appearance,studio['background'],args.width)
    if groups:environment.compose(scene,groups,appearance['environment'],studio['background'],args.width)
    snapshot=args.output.resolve()/'appearance.json';snapshot.write_bytes(args.style_profile.read_bytes())
    record={'profile':appearance,'profile_path':str(snapshot),'profile_sha256':digest(snapshot),'palette':palette,'studio':studio,
            'environment':env,'membrane':MEMBRANE,'compositor':compositor,'fonts':fonts,
            'font_sha256':[digest(args.font_zh),digest(args.font_latin)],
            'camera':{'matrix':np.array(scene.camera.matrix_world).tolist(),'scale':scene.camera.data.ortho_scale},
            'scope':'Appearance of the existing source-bounded dual-view loop; no new geometry, catalytic conformation, proton trajectory or measured context.'}
    scene['vatpase_appearance']=json.dumps(record)
    return record
