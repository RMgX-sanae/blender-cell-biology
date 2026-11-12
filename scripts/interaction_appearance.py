"""Explicit reference appearance for the existing interaction geometry and motion."""
import hashlib,json,re
from pathlib import Path
import bpy,numpy as np
from mathutils import Vector
import detail_scene as ds
from compositor_stack import build_stack


def read_profile(path,actors):
    profile=json.loads(Path(path).read_text())
    allowed={'name','reference','reference_page','usage','observations','implementation_assumptions',
             'surface_roles','actor_colors','actor_labels','background','world','lighting','view_transform','look','exposure',
             'camera','render_recipe','environment','representation','membrane_palette'}
    if set(profile)-allowed:raise ValueError('Unsupported interaction appearance fields: '+str(set(profile)-allowed))
    if set(profile.get('actor_colors',{}))-{a['id'] for a in actors}:raise ValueError('Appearance color refers to an unknown actor')
    if set(profile.get('actor_labels',{}))-{a['id'] for a in actors} or any(not isinstance(v,str) or not v.strip() for v in profile.get('actor_labels',{}).values()):raise ValueError('Actor labels require known actors and nonempty text')
    if set(profile.get('membrane_palette',{}))-{'heads','tails','tether'}:raise ValueError('Membrane palette accepts heads/tails/tether')
    colors=list(profile.get('actor_colors',{}).values())+[profile[k] for k in ('background','world') if k in profile]
    if any(not isinstance(c,str) or not re.fullmatch(r'#?[0-9a-fA-F]{6}',c) for c in colors):raise ValueError('Appearance colors must be sRGB hex')
    camera=profile.get('camera',{})
    if set(camera)-{'position','target','scale','type','lens','fstop','roll_degrees'}:raise ValueError('Unsupported interaction camera field')
    if ('position'in camera)!=('target'in camera):raise ValueError('Provide camera position and target together')
    for key in ('position','target'):
        if key in camera and (np.shape(camera[key])!=(3,) or not np.isfinite(camera[key]).all()):raise ValueError('Camera vectors must be finite triples')
    for key in ('scale','lens','fstop','roll_degrees'):
        if key in camera and (np.shape(camera[key])!=() or not np.isfinite(camera[key]) or key!='roll_degrees' and camera[key]<(0 if key=='fstop' else 1e-8)):raise ValueError('Invalid camera '+key)
    if camera.get('type','ORTHO') not in ('ORTHO','PERSP'):raise ValueError('Camera type must be ORTHO or PERSP')
    if camera.get('type')=='PERSP' and 'scale'in camera:raise ValueError('Perspective camera uses position/lens; scale is orthographic only')
    if 'representation'in profile and profile['representation'] not in ('surface','atoms'):raise ValueError('Interaction representation must be surface or atoms')
    if 'exposure'in profile and (np.shape(profile['exposure'])!=() or not np.isfinite(profile['exposure'])):raise ValueError('Exposure must be finite')
    for light in profile.get('lighting',[]):
        if set(light)!={'name','position','power','color','size'}:raise ValueError('Area lights require name/position/power/color/size')
        if np.shape(light['position'])!=(3,) or not np.isfinite(light['position']).all():raise ValueError('Light position must be a finite triple')
        if any(np.shape(light[k])!=() or not np.isfinite(light[k]) or light[k]<=0 for k in ('power','size')):raise ValueError('Light power and size must be positive')
        if not re.fullmatch(r'#?[0-9a-fA-F]{6}',light['color']):raise ValueError('Invalid light color')
    if 'lighting'in profile and (not profile['lighting'] or len({r['name'] for r in profile['lighting']})!=len(profile['lighting'])):raise ValueError('Provide uniquely named area lights')
    return profile


def camera_basis(camera,basis,orbit):
    if 'position'in camera:
        if orbit:raise ValueError('Explicit camera position cannot be combined with --orbit-degrees')
        direction=np.array(camera['target'])-camera['position']
        if np.linalg.norm(direction)<1e-8:raise ValueError('Camera position and target must differ')
        basis=np.array(Vector(direction).to_track_quat('-Z','Y').to_matrix())
    if camera.get('roll_degrees'):
        angle=np.radians(camera['roll_degrees']);c,s=np.cos(angle),np.sin(angle)
        basis=basis@np.array([[c,-s,0],[s,c,0],[0,0,1]])
    return basis


def studio_profile(default,appearance):
    result=dict(default)
    for key in ('background','world','lighting','view_transform','look','exposure'):
        if key in appearance:result[key]=appearance[key]
    result['camera']={}
    return result


def configure_camera(scene,settings,basis,center,width,style):
    camera=scene.camera;data=camera.data
    data.type=settings.get('type','PERSP' if style=='cinematic' else 'ORTHO')
    if data.type=='PERSP' and 'scale'in settings:raise ValueError('Perspective camera uses position/lens; scale is orthographic only')
    if data.type=='ORTHO' and 'lens'in settings:raise ValueError('Use scale for an orthographic camera; lens controls perspective')
    data.lens=settings.get('lens',72);data.sensor_fit='HORIZONTAL';data.sensor_width=36
    data.ortho_scale=settings.get('scale',width)
    if data.type=='ORTHO' and data.ortho_scale+1e-6<width:raise ValueError('Requested camera scale clips the complete path or label clearance')
    position=np.array(settings.get('position',center+basis[:,2]*width*2))
    matrix=np.eye(4);matrix[:3,:3]=basis;matrix[:3,3]=position
    from mathutils import Matrix
    camera.matrix_world=Matrix(matrix)
    distance=float(np.linalg.norm(position-center));data.dof.focus_distance=distance
    fstop=settings.get('fstop',0);data.dof.use_dof=bool(fstop)
    if fstop:data.dof.aperture_fstop=fstop
    span=data.ortho_scale if data.type=='ORTHO' else distance*data.sensor_width/data.lens
    return float(span),distance


def build_compositor(scene,appearance,background,width):
    name=appearance.get('render_recipe','neutral')
    recipes=json.loads((Path(__file__).resolve().parents[1]/'assets/render-recipes.json').read_text())['recipes']
    recipe=next((dict(r) for r in recipes if r['id']==name),None)
    if recipe is None:raise ValueError('Unknown interaction render_recipe')
    recipe.update(background=background,label_min_pixels=26*width/1280)
    return build_stack(scene,recipe)


def save_request(scene,path,appearance,studio,camera,environment,compositor):
    record={'path':str(Path(path).resolve()),'sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            'profile':appearance,'studio':studio,'camera':camera,'environment':environment,'compositor':compositor}
    scene['interaction_appearance']=json.dumps(record)
    return record
