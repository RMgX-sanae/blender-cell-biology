"""Inspect saved two-cycle geometry, timing, connectivity, visibility and bookkeeping."""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Vector
from bpy_extras.object_utils import world_to_camera_view

sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_transcription import centers
from build_translation_scene import atom


def world(ob,points):
    matrix=np.array(ob.matrix_world)
    return points @ matrix[:3,:3].T+matrix[:3,3]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path,required=True)
    p.add_argument('--substrate-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    s=bpy.context.scene
    info=json.loads(Path(bpy.data.filepath).with_name('manifest.json').read_text())
    dna=np.load(a.asset_dir/'chain_T.npz');rna=np.load(a.asset_dir/'chain_P.npz')
    ctp=np.load(a.substrate_dir/'ctp.npz');utp=np.load(a.substrate_dir/'utp.npz')
    core1=bpy.data.objects['Incoming CTP core → RNA C20']
    core2=bpy.data.objects['Incoming UTP core → modeled RNA U21']
    ppi2=bpy.data.objects['UTP beta/gamma phosphates → second PPi']
    template=bpy.data.objects['Template DNA except complete A16']
    a16=bpy.data.objects['Complete template A16 — modeled loading pose']
    primer=bpy.data.objects['RNA primer P residues 7–19']
    parent=bpy.data.objects['Authored one-position scaffold translocation']
    errors=[];cropped=[];occluded=[]
    intact=[]
    for ob,d,mask in [(template,dna,dna['residue_numbers']!=16),
                      (a16,dna,dna['residue_numbers']==16),
                      (primer,rna,rna['residue_numbers']<20)]:
        actual=centers(ob,int(mask.sum()),d['elements'][mask])
        error=float(np.max(np.abs(actual-d['atoms'][mask])))
        intact.append({'object':ob.name,'atoms':int(mask.sum()),'local_coordinate_error_nm':error})
        if error>2e-6 or ob.modifiers or ob.data.shape_keys:errors.append('Native local geometry altered: '+ob.name)
    protein_count=0
    for ob in s.objects:
        if 'source_chain' not in ob:continue
        d=np.load(a.asset_dir/f'chain_{ob["source_chain"]}.npz')
        vertices=np.empty(len(ob.data.vertices)*3,dtype=np.float32)
        ob.data.vertices.foreach_get('co',vertices)
        faces=np.empty(len(ob.data.loops),dtype=np.int32);ob.data.loops.foreach_get('vertex_index',faces)
        if not (np.array_equal(vertices.reshape(-1,3),d['vertices']) and
                np.array_equal(faces,d['faces'].ravel()) and
                np.allclose(ob.location,info['presentation']['Pol_II_offset_nm'])):
            errors.append('Pol II chain mismatch: '+ob.name)
        protein_count+=1
    if protein_count!=12:errors.append('Expected 12 complete Pol II chains')
    names=list(utp['atom_names']);pnames=list(utp['product_names'])
    count=int(utp['product_atom_count'])
    expected=utp['atoms'] @ utp['translocation_rotation'].T+utp['translocation_translation']
    core_error=float(np.max(np.abs(centers(core2,count,utp['elements'][:count])-expected[:count])))
    ppi_error=float(np.max(np.abs(centers(ppi2,len(utp['atoms'])-count,utp['elements'][count:])-expected[count:])))
    if max(core_error,ppi_error)>2e-6:errors.append('UTP core/PPi atom positions changed')
    if (len(ctp['atoms']),int(ctp['product_atom_count']),len(utp['atoms']),count)!=(29,20,29,20):
        errors.append('Two-cycle heavy-atom accounting failed')
    stages=[ob for ob in s.objects if ob.name.startswith('Cycle stage ')]
    dna_links=[];first_endpoint=None;last_bond=None;pair_distances=None;scaffold_cropped=[]
    min_entry=(float('inf'),None);min_nonbond=(float('inf'),None)
    base_ids=[names.index(n) for n in ['C4','C5','C6']]
    mask=np.ones((count,len(dna['atoms'])+len(rna['atoms'])),bool)
    for i,n in enumerate(names[:count]):
        if n in ['PA','O1A','O2A',"O5'"]:
            mask[i,len(dna['atoms']):] &= ~((rna['residue_numbers']==20)&np.isin(rna['atom_names'],["O3'","C3'"]))
    for frame in range(s.frame_start,s.frame_end+1):
        s.frame_set(frame);dg=bpy.context.evaluated_depsgraph_get();time=(frame-1)/s.render.fps
        if sum(not ob.hide_render for ob in stages)!=1:errors.append(f'Caption count at {frame}')
        if not bpy.data.objects['Non-template DNA — hidden whole chain'].hide_render:
            errors.append('Non-template strand became visible')
        t16=world(a16,np.array([atom(dna,16,'P'),atom(dna,16,"O3'")]))
        neighbors=world(parent,np.array([atom(dna,15,"O3'"),atom(dna,17,'P')]))
        dna_links.extend(np.linalg.norm(t16-neighbors,axis=1).tolist())
        points=world(parent,np.concatenate([dna['atoms'],rna['atoms']]))
        points[np.flatnonzero(dna['residue_numbers']==16)]=world(a16,dna['atoms'][dna['residue_numbers']==16])
        camera=np.array(s.camera.matrix_world.inverted())
        cp=points @ camera[:3,:3].T+camera[:3,3]
        half_width=s.camera.data.ortho_scale*.5
        half_height=half_width*s.render.resolution_y/s.render.resolution_x
        if np.abs(cp[:,0]).max()+.14>half_width or np.abs(cp[:,1]).max()+.14>half_height:
            scaffold_cropped.append(frame)
        if frame==265:
            first_endpoint=float(np.max(np.abs(world(core1,centers(core1,20,ctp['elements'][:20],dg))-ctp['product_atoms'])))
        if time>=13.5:
            cur=world(core2,expected[:count])
            env=world(parent,np.concatenate([dna['atoms'],rna['atoms']]))
            env[np.flatnonzero(dna['residue_numbers']==16)]=world(a16,dna['atoms'][dna['residue_numbers']==16])
            distance=np.linalg.norm(cur[:,None,:]-env[None,:,:],axis=2)
            if time<=17.5 and distance.min()<min_entry[0]:min_entry=(float(distance.min()),frame)
            distance[~mask]=np.inf
            if distance.min()<min_nonbond[0]:min_nonbond=(float(distance.min()),frame)
            if time>=16.5:
                point=Vector(cur[base_ids].mean(0));front=s.camera.rotation_euler.to_matrix() @ Vector((0,0,1))
                hit=s.ray_cast(dg,point+front*100,-front)
                if hit[0] and hit[4]!=core2 and (hit[1]-point).dot(front)>.15:
                    occluded.append({'frame':frame,'object':hit[4].name})
        # The first PPi deliberately crosses the picture boundary while leaving.
        visible=[ob for ob in [core1,core2,ppi2] if not ob.hide_render and
                 (ob==core1 or time>=16.5)]
        visible += [ob for ob in s.objects if ob.type=='FONT' and not ob.hide_render and not ob.name.startswith('Outline ')]
        for ob in visible:
            corners=[world_to_camera_view(s,s.camera,ob.matrix_world @ Vector(v)) for v in ob.bound_box]
            if min(v.x for v in corners)<0 or max(v.x for v in corners)>1 or min(v.y for v in corners)<0 or max(v.y for v in corners)>1:
                cropped.append([frame,ob.name])
        if frame==s.frame_end:
            last_bond=float(np.linalg.norm(world(parent,np.array([atom(rna,20,"O3'")]))[0]-cur[names.index('PA')]))
            a_atoms={n:world(a16,np.array([atom(dna,16,n)]))[0] for n in ['N1','N6']}
            pair_distances=[float(np.linalg.norm(cur[names.index('N3')]-a_atoms['N1'])),
                            float(np.linalg.norm(cur[names.index('O4')]-a_atoms['N6']))]
    if first_endpoint is None or first_endpoint>2e-6:errors.append('First experimental C20 endpoint differs')
    if not (.14<min(dna_links) and max(dna_links)<.18):errors.append('DNA junction length during modeled loading')
    if not .15<last_bond<.17:errors.append('Second RNA connection length')
    if min_entry[0]<.18:errors.append('UTP entry intersects an existing nucleotide')
    if min_nonbond[0]<.18:errors.append('UTP core nonbonded clash during joining')
    if cropped:errors.append('Caption or retained reaction actor cropped')
    if scaffold_cropped:errors.append('Native nucleic-acid scaffold cropped')
    if occluded:errors.append('U base obscured during its explanation')
    result={'blend':bpy.data.filepath,'frames':s.frame_end,'fps':s.render.fps,'seconds':s.frame_end/s.render.fps,
        'native_geometry':intact,'intact_Pol_II_chains':protein_count,
        'CTP_and_UTP_heavy_atoms_each':[29,20,9],'UTP_core_error_nm':core_error,'PPi_error_nm':ppi_error,
        'first_experimental_endpoint_error_nm':first_endpoint,'second_RNA_bond_length_nm':last_bond,
        'DNA_link_length_range_nm':[min(dna_links),max(dna_links)],
        'modeled_A_U_pair_contact_distances_nm':pair_distances,
        'minimum_entry_heavy_atom_distance_nm_and_frame':min_entry,
        'minimum_nonbonded_core_distance_nm_and_frame':min_nonbond,
        'cropped':cropped,'scaffold_cropped':scaffold_cropped,'occluded':occluded,'errors':errors,
        'scope':'Saved geometry and teaching-animation checks only; U21/loading/translocation are modeled, not validated catalytic dynamics.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False))
    if errors:raise RuntimeError('Two-cycle audit failed; inspect report.')


if __name__=='__main__':main()
