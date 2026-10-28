"""Audit the saved one-C teaching scene; does not certify catalytic intermediates."""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Vector
from bpy_extras.object_utils import world_to_camera_view


def centers(obj, count, elements, depsgraph=None):
    mesh = obj.data if depsgraph is None else obj.evaluated_get(depsgraph).data
    coords = np.empty(len(mesh.vertices)*3,dtype=np.float32)
    mesh.vertices.foreach_get('co',coords)
    result = coords.reshape(-1,3)[:count*62:62].copy()
    radii={'C':.105,'N':.103,'O':.10,'P':.14,'S':.13}
    result[:,2] -= [radii.get(e,.105) for e in elements]
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path,required=True)
    p.add_argument('--substrate-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    s=bpy.context.scene
    manifest=json.loads(Path(bpy.data.filepath).with_name('manifest.json').read_text())
    errors=[];cropped=[];occluded=[]
    intact=[]
    mapping={'T':'Template DNA — all modeled residues',
             'N':'Non-template DNA — hidden whole chain','P':'RNA primer P residues 7–19'}
    data={ch:np.load(a.asset_dir/f'chain_{ch}.npz') for ch in mapping}
    for ch,name in mapping.items():
        d=data[ch];mask=d['residue_numbers']<20 if ch=='P' else np.ones(len(d['atoms']),bool)
        ob=bpy.data.objects[name]
        ok=np.allclose(centers(ob,int(mask.sum()),d['elements'][mask]),d['atoms'][mask],atol=2e-6)
        ok=ok and len(ob.modifiers)==0 and ob.data.shape_keys is None
        intact.append({'chain':ch,'modeled_atoms':int(mask.sum()),'source_retained':bool(ok)})
        if not ok:errors.append('Scaffold geometry mismatch '+ch)
    offset=manifest['presentation']['Pol_II_offset_nm']
    context_count=0
    for ob in s.objects:
        if 'source_chain' not in ob:continue
        d=np.load(a.asset_dir/f'chain_{ob["source_chain"]}.npz')
        v=np.empty(len(ob.data.vertices)*3,dtype=np.float32);ob.data.vertices.foreach_get('co',v)
        f=np.empty(len(ob.data.loops),dtype=np.int32);ob.data.loops.foreach_get('vertex_index',f)
        if not (np.array_equal(v.reshape(-1,3),d['vertices']) and np.array_equal(f,d['faces'].ravel())
                and np.allclose(ob.location,offset) and
                all(m.type=='SMOOTH' and abs(m.factor-.35)<1e-6 and m.iterations==2 for m in ob.modifiers)):
            errors.append('Pol II integrity/offset mismatch '+ob.name)
        context_count+=1
    if context_count!=12:errors.append('Expected twelve complete Pol II protein chains')
    ctp=np.load(a.substrate_dir/'ctp.npz');count=int(ctp['product_atom_count'])
    core=bpy.data.objects['Incoming CTP core → RNA C20']
    ppi=bpy.data.objects['CTP beta/gamma phosphates → PPi']
    if not np.allclose(centers(ppi,len(ctp['atoms'])-count,ctp['elements'][count:]),ctp['atoms'][count:],atol=2e-6):
        errors.append('PPi atoms changed')
    atom_count=len(ctp['atoms'])
    if atom_count!=29 or count!=20:errors.append('CTP atom bookkeeping failed')
    stages=[ob for ob in s.objects if ob.name.startswith('Stage ')]
    base_indices=[list(ctp['atom_names']).index(n) for n in ['C4','C5','C6']]
    substrate_error=product_error=None
    for frame in range(s.frame_start,s.frame_end+1):
        s.frame_set(frame);dg=bpy.context.evaluated_depsgraph_get()
        current=centers(core,count,ctp['elements'][:count],dg)
        if frame==s.frame_start:
            substrate_error=float(np.max(np.abs(current-ctp['atoms'][:count])))
        if frame==s.frame_end:
            product_error=float(np.max(np.abs(current-ctp['product_atoms'])))
        if sum(not ob.hide_render for ob in stages)!=1:errors.append(f'Caption count at {frame}')
        if not bpy.data.objects[mapping['N']].hide_render:errors.append('Non-template chain unexpectedly visible')
        for ob in [core,ppi]:
            corners=[world_to_camera_view(s,s.camera,ob.matrix_world@Vector(v)) for v in ob.bound_box]
            if min(v.x for v in corners)<0 or max(v.x for v in corners)>1 or min(v.y for v in corners)<0 or max(v.y for v in corners)>1:
                cropped.append([frame,ob.name])
        if frame==s.frame_end or .40<(frame-1)/(s.frame_end-1)<.60:
            point=core.matrix_world@Vector(current[base_indices].mean(0))
            front=s.camera.rotation_euler.to_matrix()@Vector((0,0,1))
            hit=s.ray_cast(dg,point+front*100,-front)
            if hit[0] and hit[4]!=core and (hit[1]-point).dot(front)>.15:
                occluded.append({'frame':frame,'object':hit[4].name})
        for ob in [o for o in s.objects if o.type=='FONT' and not o.hide_render and not o.name.startswith('Outline ')]:
            corners=[world_to_camera_view(s,s.camera,ob.matrix_world@Vector(v)) for v in ob.bound_box]
            if min(v.x for v in corners)<0 or max(v.x for v in corners)>1 or min(v.y for v in corners)<0 or max(v.y for v in corners)>1:
                cropped.append([frame,ob.name])
    if max(substrate_error,product_error)>2e-6:errors.append('Substrate/product endpoint coordinates differ')
    if cropped:errors.append('Objects or labels outside camera frame')
    if occluded:errors.append('Added C base hidden during explanation')
    source=data['P']
    source_point=lambda res,name:source['atoms'][(source['residue_numbers']==res)&(source['atom_names']==name)][0]
    length=float(np.linalg.norm(source_point(19,"O3'")-source_point(20,'P')))
    if not .14<length<.18:errors.append('Unexpected product linkage length')
    result={'blend':bpy.data.filepath,'frames':s.frame_end,'intact_scaffold':intact,
        'intact_Pol_II_chains':context_count,'CTP_to_residue_plus_PPi_atoms':[atom_count,count,atom_count-count],
        'substrate_endpoint_error_nm':substrate_error,'product_endpoint_error_nm':product_error,
        'product_bond_length_nm':length,'cropped':cropped,'occluded':occluded,'errors':errors,
        'scope':'Saved geometry, atom bookkeeping, observed endpoint, complete subunits, camera and captions. Authored interpolation is not a validated catalytic trajectory.'}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['frames','CTP_to_residue_plus_PPi_atoms','product_endpoint_error_nm','errors']}))
    if errors:raise RuntimeError('Transcription scene audit failed; see report.')


if __name__=='__main__':main()
