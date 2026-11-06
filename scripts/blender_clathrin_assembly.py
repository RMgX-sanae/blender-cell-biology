"""Build or audit intact modeled clathrin units using shared source meshes."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Matrix,Vector
from bpy_extras.object_utils import world_to_camera_view

sys.path.insert(0,str(Path(__file__).resolve().parent))
import detail_scene as ds


def load_assets(folder):
    path=folder/'assembly.json';m=json.loads(path.read_text())
    if m['units']!='nm' or not m['grouping_check']['passed']:raise ValueError('Require a qualified nm assembly')
    expected={'source-selection.json','prototypes/manifest.json'}|{'prototypes/chain_'+c+'.npz' for u in m['instances'] for c in u['chains']}
    if set(m['data_sha256'])!=expected:raise ValueError('Incomplete prepared file inventory')
    for name,sha in m['data_sha256'].items():
        if hashlib.sha256((folder/name).read_bytes()).hexdigest()!=sha:raise ValueError('Prepared file changed: '+name)
    return m,hashlib.sha256(path.read_bytes()).hexdigest()


def build(args):
    if args.output.exists():raise ValueError('Use a fresh output path')
    m,fingerprint=load_assets(args.assets)
    ids={u['id'] for u in m['instances']};visible=set(args.show or ids)
    if not visible or not visible<=ids or args.highlight not in ids:raise ValueError('Unknown or empty whole-unit selection')
    profile_path=args.style_profile or Path(__file__).resolve().parents[1]/'assets/style-molecular-reference.json'
    profile=json.loads(profile_path.read_text())
    if args.style_profile is None:
        profile.update(background='0b1420',bump=0,exposure=0,roughness=.6,roughness_background=.65)
        profile['palette'].update(alpha='#447786',protein='#ffc34d',dna1='#f2a5bf')
    bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
    group=ds.collection('Complete modeled clathrin units');prototypes={};shown=[]
    mats={role:ds.material('Clathrin '+role,profile['palette'][key],
                          profile.get('roughness_background',profile['roughness']) if role=='background' else profile['roughness'],profile['bump'])
          for role,key in [('background','alpha'),('heavy','protein'),('light','dna1')]}
    heavy={c for g in m['source_config']['groups'] for c in g['heavy']}
    for unit in m['instances']:
        rig=bpy.data.objects.new(unit['id'],None);group.objects.link(rig);rig.matrix_world=Matrix(unit['matrix_nm'])
        rig['assembly_unit_id']=unit['id']
        for c in unit['chains']:
            name=unit['id']+'/'+c
            if c not in prototypes:
                data=np.load(args.assets/'prototypes'/('chain_'+c+'.npz'))
                obj=ds.mesh(name,data['vertices'],data['faces'],[mats['background']],group);prototypes[c]=obj
            else:obj=ds.duplicate(prototypes[c],name,group)
            obj.parent=rig;obj.matrix_basis=Matrix.Identity(4)
            obj['source_chain']=c;obj['assembly_unit_id']=unit['id'];obj['source_reference_chain']=unit['reference_chains'][c]
            obj.hide_render=unit['id'] not in visible
            obj.material_slots[0].link='OBJECT'
            obj.material_slots[0].material=mats[('heavy' if c in heavy else 'light') if unit['id']==args.highlight else 'background']
            if unit['id'] in visible:shown.append(obj)
    bpy.context.view_layer.update()
    target=np.mean([np.array(obj.matrix_world)[:3,:3]@np.mean([list(v) for v in obj.bound_box],axis=0)+np.array(obj.matrix_world)[:3,3] for obj in shown],axis=0)
    hub=np.array(next(u for u in m['instances'] if u['id']==args.highlight)['hub_nm']);front=hub/np.linalg.norm(hub)
    args.frames=1;args.fps=24
    camera=ds.studio(profile,args,{'position':(target+front*180).tolist(),'target':target.tolist(),'scale':110,'light_scale':4})
    ds.fit_orthographic(camera,shown,margin=.09)
    scene=bpy.context.scene;scene['assembly_asset_dir']=str(args.assets.resolve());scene['assembly_sha256']=fingerprint
    scene['assembly_visible_units']=json.dumps(sorted(visible));scene['assembly_highlight']=args.highlight
    scene['assembly_scope']=m['state'];scene['source_limits']=m['source_config']['limits']
    scene['assembly_compositor']=args.compositor
    scene['assembly_style_source']=str(profile_path.resolve())
    args.output.parent.mkdir(parents=True,exist_ok=True);bpy.ops.wm.save_as_mainfile(filepath=str(args.output))


def audit():
    scene=bpy.context.scene;folder=Path(scene['assembly_asset_dir']);m,fingerprint=load_assets(folder)
    if fingerprint!=scene['assembly_sha256']:raise ValueError('Prepared assembly manifest changed after build')
    visible=set(json.loads(scene['assembly_visible_units']));ids={u['id'] for u in m['instances']}
    if not visible or not visible<=ids:raise ValueError('Invalid whole-unit visibility declaration')
    expected={u['id']+'/'+c for u in m['instances'] for c in u['chains']}
    actual={o.name for o in scene.objects if o.type=='MESH'}
    if actual!=expected:raise ValueError('Saved chain inventory differs from complete source units')
    meshes={};cache={};minimum=1.;maximum=0.
    for unit in m['instances']:
        rig=scene.objects.get(unit['id'])
        if rig is None or rig.get('assembly_unit_id')!=unit['id'] or not np.allclose(np.asarray(rig.matrix_world),unit['matrix_nm'],atol=2e-6,rtol=0):raise ValueError('Source unit transform changed')
        for c in unit['chains']:
            obj=scene.objects[unit['id']+'/'+c]
            if c not in cache:cache[c]=dict(np.load(folder/'prototypes'/('chain_'+c+'.npz')))
            data=cache[c]
            if obj.parent!=rig or obj.get('source_chain')!=c or obj.get('assembly_unit_id')!=unit['id'] or obj.get('source_reference_chain')!=unit['reference_chains'][c]:raise ValueError('Source chain ownership changed')
            if obj.hide_render!=(unit['id'] not in visible) or obj.modifiers:raise ValueError('Whole chain hidden or modified')
            error=float(abs(np.asarray(obj.matrix_world)-unit['matrix_nm']).max());maximum=max(maximum,error)
            if error>2e-6:raise ValueError('Chain no longer follows its complete unit')
            if obj.data.as_pointer() not in meshes:
                coords=np.empty(len(obj.data.vertices)*3);obj.data.vertices.foreach_get('co',coords)
                faces=np.array([list(f.vertices) for f in obj.data.polygons])
                if coords.size!=data['vertices'].size or not np.allclose(coords.reshape(-1,3),data['vertices'],atol=2e-6,rtol=0) or not np.array_equal(faces,data['faces']):raise ValueError('Source mesh changed or cut')
                meshes[obj.data.as_pointer()]=c
            elif meshes[obj.data.as_pointer()]!=c:raise ValueError('Different chains share the wrong mesh')
            if unit['id'] in visible:
                for point in obj.bound_box:
                    v=world_to_camera_view(scene,scene.camera,obj.matrix_world@Vector(point));minimum=min(minimum,v.x,v.y,1-v.x,1-v.y)
                    if min(v.x,v.y,1-v.x,1-v.y)<0 or v.z<=0:raise ValueError('Whole visible chain leaves camera')
    if scene.frame_start!=1 or scene.frame_end!=1:raise ValueError('This source assembly is a static model')
    return {'passed':True,'scene':bpy.data.filepath,'units':len(ids),'chains':len(expected),'shared_source_meshes':len(meshes),
       'visible_units':sorted(visible),'hidden_whole_units':sorted(ids-visible),'transform_maximum_error':maximum,'camera_minimum_margin':minimum,
       'prepared_reference_check':m['independent_reference_check'],'limits':m['source_config']['limits']}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['build','audit'])
    p.add_argument('--assets',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--show',nargs='+');p.add_argument('--highlight',default='1-ABC')
    p.add_argument('--width',type=int,default=1280);p.add_argument('--samples',type=int,default=24)
    p.add_argument('--engine',choices=['cycles','eevee'],default='cycles')
    p.add_argument('--style-profile',type=Path,help='Project profile using the detail_scene lighting/palette schema')
    p.add_argument('--compositor',choices=['clean','reference','monochrome'],default='clean')
    args=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if args.action=='build':
        if args.assets is None:p.error('Build requires --assets')
        build(args)
    else:
        result=audit();args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
