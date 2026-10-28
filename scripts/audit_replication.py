"""Check the saved fork: intact units, complementary daughters, growth and polarity."""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Vector
from bpy_extras.object_utils import world_to_camera_view


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path,required=True)
    p.add_argument('--replisome-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    s=bpy.context.scene;m=json.loads(Path(bpy.data.filepath).with_name('manifest.json').read_text())
    errors=[];cropped=[];count=m['nucleotide_count_per_parent']
    objects={(ob['strand_role'],int(ob['index'])):ob for ob in s.objects if 'strand_role' in ob}
    if len(objects)!=4*count:errors.append('Missing strand units')
    radii={'C':.105,'N':.103,'O':.10,'P':.14,'S':.13}
    checked_meshes=set();intact_meshes=0
    for ob in objects.values():
        if ob['representation']=='RNA primer symbol' or ob.data.name in checked_meshes:continue
        checked_meshes.add(ob.data.name);d=np.load(a.asset_dir/(ob['prototype']+'.npz'))
        v=np.empty(len(ob.data.vertices)*3,dtype=np.float32);ob.data.vertices.foreach_get('co',v)
        centers=v.reshape(-1,3)[:len(d['atoms'])*62:62].copy()
        centers[:,2]-=[radii.get(e,.105) for e in d['elements']]
        if not (np.allclose(centers,d['atoms'],atol=2e-6) and len(ob.modifiers)==0 and ob.data.shape_keys is None):
            errors.append('Prototype geometry altered '+ob.name)
        intact_meshes+=1
    cmg=[ob for ob in s.objects if 'source_chain' in ob]
    if {ob['source_chain'] for ob in cmg}!=set('234567CDEFG'):errors.append('Incorrect CMG chain set')
    front=s.camera.rotation_euler.to_matrix()@Vector((0,0,1))
    cmg_max={}
    for ob in cmg:
        d=np.load(a.replisome_dir/f'chain_{ob["source_chain"]}.npz')
        v=np.empty(len(ob.data.vertices)*3,dtype=np.float32);ob.data.vertices.foreach_get('co',v)
        f=np.empty(len(ob.data.loops),dtype=np.int32);ob.data.loops.foreach_get('vertex_index',f)
        if not (np.array_equal(v.reshape(-1,3),d['vertices']) and np.array_equal(f,d['faces'].ravel())
                and all(mod.type=='SMOOTH' for mod in ob.modifiers)):
            errors.append('CMG chain changed '+ob['source_chain'])
        cmg_max[ob.name]=float((d['vertices']@np.array(front)).max())
    previous_lead=set();previous_lag=[set() for _ in range(count//12)]
    final_counts={};minimum_depth_gap=float('inf');rigid_error=0.;pair_error=0.
    for frame in range(s.frame_start,s.frame_end+1):
        s.frame_set(frame);state=m['frame_poses'][frame-1]
        roles={role:{i for (r,i),ob in objects.items() if r==role and not ob.hide_render}
               for role in ['parent_A','parent_B','leading','lagging']}
        if roles['parent_A']!=set(range(count)) or roles['parent_B']!=set(range(count)):
            errors.append(f'Parent identity lost at {frame}')
        lead=roles['leading']
        if not lead or lead!=set(range(max(lead)+1)) or not previous_lead<=lead:
            errors.append(f'Leading synthesis not a growing prefix at {frame}')
        previous_lead=lead
        if lead and max(lead)>state['fork_index']-5+1e-6:errors.append(f'Leading DNA ahead of fork at {frame}')
        for block in range(count//12):
            present={i for i in roles['lagging'] if i//12==block}
            primer={i for i in present if i%12>=9};dna=present-primer;j=block*12+11
            if dna and (primer!={j-2,j-1,j} or dna!=set(range(min(dna),j-2))):
                errors.append(f'Lagging DNA lacks its primer or grows in the wrong order at {frame}')
            if not previous_lag[block]<=dna:errors.append(f'Lagging identity lost at {frame}')
            previous_lag[block]=dna
            if present and j+5>state['fork_index']+1e-6:errors.append(f'Primer ahead of exposed template at {frame}')
        min_actor=float('inf')
        for (role,i),ob in objects.items():
            if ob.hide_render:continue
            r=np.array(ob.matrix_world.to_3x3());rigid_error=max(rigid_error,float(np.max(np.abs(r.T@r-np.eye(3)))))
            if ob.hide_viewport:errors.append(f'Visible object hidden in viewport at {frame}')
            corners=[ob.matrix_world@Vector(v) for v in ob.bound_box]
            min_actor=min(min_actor,min(v.dot(front) for v in corners))
            xy=[world_to_camera_view(s,s.camera,v) for v in corners]
            if min(v.x for v in xy)<0 or max(v.x for v in xy)>1 or min(v.y for v in xy)<0 or max(v.y for v in xy)>1:
                cropped.append([frame,ob.name])
            if role in ['leading','lagging']:
                parent=objects['parent_B' if role=='leading' else 'parent_A',i]
                pair_error=max(pair_error,float(np.max(np.abs(np.array(ob.matrix_world)-np.array(parent.matrix_world)))))
        positions=[np.array(ob.location) for ob in cmg]
        if not all(np.allclose(x,positions[0],atol=2e-6) for x in positions):errors.append(f'CMG subunits move apart at {frame}')
        max_cmg=max(cmg_max[ob.name]+ob.location.dot(front) for ob in cmg)
        minimum_depth_gap=min(minimum_depth_gap,min_actor-max_cmg)
        for ob in [o for o in s.objects if o.type=='FONT' and not o.hide_render and not o.name.startswith('Outline ')]:
            xy=[world_to_camera_view(s,s.camera,ob.matrix_world@Vector(v)) for v in ob.bound_box]
            if min(v.x for v in xy)<0 or max(v.x for v in xy)>1 or min(v.y for v in xy)<0 or max(v.y for v in xy)>1:
                cropped.append([frame,ob.name])
        final_counts={k:len(v) for k,v in roles.items()}
    if rigid_error>3e-6:errors.append('Nucleotide unit scaled or deformed')
    if pair_error>2e-6:errors.append('Daughter complement does not share its parent track')
    if minimum_depth_gap<=0:errors.append('CMG reference layer overlaps foreground depth')
    if cropped:errors.append('Molecules or labels outside frame')
    report={'blend':bpy.data.filepath,'frames':s.frame_end,'intact_unique_nucleotide_meshes':intact_meshes,
        'intact_CMG_chains':len(cmg),'final_visible_units':final_counts,
        'maximum_rigid_error':rigid_error,'maximum_pair_transform_error':pair_error,
        'minimum_CMG_depth_gap_nm':minimum_depth_gap,'cropped':cropped,'errors':errors,
        'scope':'Actual saved identities, rigid units, complement pairing, synthesis order, symbolic primer prerequisites, intact offset CMG and camera. Not a biological time calibration or atomic simulation.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['frames','final_visible_units','minimum_CMG_depth_gap_nm','errors']}))
    if errors:raise RuntimeError('Replication audit failed; see report.')


if __name__=='__main__':main()
