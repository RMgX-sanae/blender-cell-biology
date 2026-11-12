"""Explicit appearance input for the authored frontal actin chapter layout."""
import hashlib,json,math,re
from pathlib import Path

PALETTE={'teal','blue','gold','pink','green','red','cream','dim','ink'}
ROLES={'protein','cofactor','schematic_symbol','schematic_linker','membrane_symbol','environment'}
CHAPTERS={prefix+str(i).zfill(2) for prefix,count in [('a',8),('f',7),('c',6)] for i in range(1,count+1)}


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_profile(path):
    profile=json.loads(Path(path).read_text())
    allowed={'name','reference','reference_page','usage','observations','implementation_assumptions',
             'surface_roles','palette','background','world','lighting','view_transform','look','exposure',
             'camera','render_recipe','environment','chapters'}
    if set(profile)-allowed:raise ValueError('Unsupported actin appearance fields: '+str(set(profile)-allowed))
    if set(profile.get('palette',{}))-PALETTE:raise ValueError('Unknown actin semantic palette key')
    if set(profile.get('surface_roles',{}))-ROLES:raise ValueError('Unknown actin material role')
    colors=list(profile.get('palette',{}).values())+[profile[k] for k in ('background','world') if k in profile]
    if any(not isinstance(c,str) or not re.fullmatch(r'#?[0-9a-fA-F]{6}',c) for c in colors):raise ValueError('Actin appearance colors must be sRGB hex')
    chapters=profile.get('chapters',{})
    if set(chapters)-CHAPTERS:raise ValueError('Unknown actin chapter override')
    if any(set(row)-{'camera','environment'} for row in chapters.values()):raise ValueError('Chapter overrides accept camera/environment only')
    for row in [profile,*chapters.values()]:
        camera=row.get('camera',{})
        if set(camera)-{'scale_multiplier','center_offset_nm'}:raise ValueError('Authored frontal layout supports camera scale_multiplier/center_offset_nm; perspective/orbit require reauthoring the diagram layout')
        scale=camera.get('scale_multiplier',1);offset=camera.get('center_offset_nm',[0,0])
        if not isinstance(scale,(int,float)) or not math.isfinite(scale) or scale<=0:raise ValueError('Camera scale multiplier must be finite and positive')
        if len(offset)!=2 or any(not isinstance(v,(int,float)) or not math.isfinite(v) for v in offset):raise ValueError('Camera center offset is a finite x/z pair in nm')
        env=row.get('environment')
        if env is not None and (not isinstance(env,dict) or env.get('compartment')!='cytosol' or not env.get('context_note') or not env.get('context_sources')):raise ValueError('Actin background requires an explicitly reviewed cytosolic context')
    lights=profile.get('lighting')
    if lights is not None:
        if not lights or len({l['name'] for l in lights})!=len(lights):raise ValueError('Provide uniquely named area lights')
        for light in lights:
            if set(light)!={'name','position','power','color','size'}:raise ValueError('Area lights require name/position/power/color/size')
            if len(light['position'])!=3 or any(not math.isfinite(v) for v in light['position']):raise ValueError('Light position must be finite')
            if any(not math.isfinite(light[k]) or light[k]<=0 for k in ('power','size')):raise ValueError('Light power and size must be positive')
            if not re.fullmatch(r'#?[0-9a-fA-F]{6}',light['color']):raise ValueError('Invalid light color')
    if 'exposure'in profile and not math.isfinite(profile['exposure']):raise ValueError('Exposure must be finite')
    return profile


def chapter_profile(profile,chapter):
    if profile is None:return None
    result={k:v for k,v in profile.items() if k!='chapters'}
    override=profile.get('chapters',{}).get(chapter,{})
    result['camera']={**result.get('camera',{}),**override.get('camera',{})}
    if 'environment'in override:result['environment']=override['environment']
    return result


def load_saved(options):
    path=options.get('style_profile')
    if not path:return None
    if digest(path)!=options['style_profile_sha256']:raise ValueError('Project appearance snapshot changed; provide a new explicit profile and rebuild')
    return read_profile(path)


def configure(options,old,path,environment_dir,root,stages):
    if path is None and old.get('style_profile'):load_saved(old)
    source=path or old.get('style_profile');profile=read_profile(source) if source else None
    options.update(style_profile=None,style_profile_sha256=None,appearance_recipe=None,environment_dir=None,environment_sha256=None)
    if profile is not None:
        sha=digest(source);options.update(style_profile=str(root/'appearance'/(sha+'.json')),style_profile_sha256=sha,appearance_recipe=profile.get('render_recipe','neutral'))
    required=profile is not None and any(chapter_profile(profile,c).get('environment') for c in CHAPTERS)
    env=environment_dir or (old.get('environment_dir') if required else None)
    if bool(env)!=required:raise ValueError('Provide --environment-dir exactly when the profile declares source-shaped backgrounds')
    if env:
        env=Path(env).resolve();files=[env/'manifest.json',*sorted(env.glob('*.npz'))]
        manifest=json.loads(files[0].read_text())
        if manifest.get('kind')!='schematic-cytosolic-environment' or manifest.get('units')!='nm':raise ValueError('Use prepared native-nm cytosolic environment assets')
        hashes={p.name:digest(p) for p in files}
        options.update(environment_dir=str(env),environment_sha256=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest())
    changed=any(old.get(k)!=options[k] for k in ('style_profile_sha256','environment_sha256','appearance_recipe'))
    if old and changed and 'build' not in stages and set(stages)&{'render','encode','audit'}:raise ValueError('Changed appearance requires rebuilding affected chapters before rendering or encoding')
    if profile is not None:
        target=Path(options['style_profile']);target.parent.mkdir(exist_ok=True)
        if target.exists():
            if digest(target)!=options['style_profile_sha256']:raise ValueError('Existing immutable appearance snapshot is corrupted')
        else:target.write_bytes(Path(source).read_bytes())
    return profile
