"""Audit the reopened translation scene: identities, intact atoms, codon and visibility."""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Vector


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    scene = bpy.context.scene
    manifest = json.loads(Path(bpy.data.filepath).with_name('manifest.json').read_text())
    errors = []
    intact = []
    mesh_vertices = {}
    data = {c: np.load(args.asset_dir/f'chain_{c}.npz') for c in ['At','Pt','Et','Mr']}
    # nucleic_atoms creates 62 vertices per atom before its bond geometry.
    # Compare the actual sphere centers to every modeled heavy-atom position.
    radii = {'C':.105,'N':.103,'O':.10,'P':.14,'S':.13}
    for ch in ['At','Pt','Et']:
        obj = bpy.data.objects['tRNA identity '+ch]
        count = len(data[ch]['atoms'])
        coords = np.empty(len(obj.data.vertices)*3, dtype=np.float32)
        obj.data.vertices.foreach_get('co', coords)
        mesh_vertices[ch] = coords.reshape(-1,3).copy()
        centers = coords.reshape(-1,3)[:count*62:62].copy()
        centers[:,2] -= [radii.get(e,.105) for e in data[ch]['elements']]
        good = (len(obj.data.vertices) == count*62+len(data[ch]['bonds'])*12 and
                np.allclose(centers,data[ch]['atoms'],atol=2e-6) and
                obj.data.shape_keys is None and len(obj.modifiers)==0)
        intact.append({'chain':ch,'modeled_atoms':count,'intact':bool(good)})
        if not good: errors.append('Source-derived tRNA geometry changed: '+ch)
    offset = np.array(manifest['presentation']['large_subunit_offset_nm'])
    context_count = 0
    for obj in scene.objects:
        if obj.type!='MESH' or 'source_chain' not in obj: continue
        chain = obj['source_chain']
        d = np.load(args.asset_dir/f'chain_{chain}.npz')
        coords = np.empty(len(obj.data.vertices)*3, dtype=np.float32)
        obj.data.vertices.foreach_get('co',coords)
        loops = np.empty(len(obj.data.loops),dtype=np.int32)
        obj.data.loops.foreach_get('vertex_index',loops)
        expected_offset = offset if chain.startswith('L') else np.zeros(3)
        if not (np.array_equal(coords.reshape(-1,3),d['vertices']) and
                np.array_equal(loops,d['faces'].ravel()) and
                np.allclose(obj.location,expected_offset)):
            errors.append('Context geometry/offset mismatch: '+chain)
        context_count += 1
    def point(chain,res,name):
        d=data[chain]
        return Vector(d['atoms'][(d['residue_numbers']==res)&(d['atom_names']==name)][0])
    cropped=[]
    occluded=[]
    rigid_errors=[]
    attachment_errors=[]
    entry_mrna_nearest={'distance_nm':float('inf'),'frame':None}
    mr_objects=[o for o in scene.objects if o.get('source_nucleic_chain')=='Mr']
    assert scene.camera.data.type=='ORTHO', 'This scene audit expects its orthographic camera.'
    camera_frame=scene.camera.data.view_frame(scene=scene)
    left,right=min(v.x for v in camera_frame),max(v.x for v in camera_frame)
    bottom,top=min(v.y for v in camera_frame),max(v.y for v in camera_frame)
    stages=[o for o in scene.objects if o.name.startswith('Stage ')]
    for frame in range(scene.frame_start,scene.frame_end+1):
        scene.frame_set(frame)
        state=manifest['frame_poses'][frame-1]
        for ch in ['At','Pt','Et']:
            ob=bpy.data.objects['tRNA identity '+ch]
            r=np.array(ob.matrix_world.to_3x3())
            if not np.allclose(r.T@r,np.eye(3),atol=2e-6) or np.linalg.det(r)<.99999:
                rigid_errors.append([frame,ch])
        for ch in ['At','Pt']:
            ob=bpy.data.objects['tRNA identity '+ch]
            matrix=np.array(scene.camera.matrix_world.inverted()@ob.matrix_world)
            xyz=mesh_vertices[ch]@matrix[:3,:3].T+matrix[:3,3]
            # Rotated bounding-box corners can leave the frame even while
            # every visible vertex is inside. Check the rendered geometry.
            if (xyz[:,0].min()<left or xyz[:,0].max()>right or
                    xyz[:,1].min()<bottom or xyz[:,1].max()>top):
                cropped.append([frame,ch])
        if state['dock']<1:
            matrix=np.array(bpy.data.objects['tRNA identity At'].matrix_world)
            a=data['At']['atoms']@matrix[:3,:3].T+matrix[:3,3]
            mr_points=[]
            for ob in mr_objects:
                matrix=np.array(ob.matrix_world)
                xyz=data['Mr']['atoms'][data['Mr']['residue_numbers']==ob['residue_number']]
                mr_points.append(xyz@matrix[:3,:3].T+matrix[:3,3])
            b=np.concatenate(mr_points)
            distance=float(np.sqrt(np.sum((a[:,None,:]-b[None,:,:])**2,axis=2).min()))
            if distance<entry_mrna_nearest['distance_nm']:
                entry_mrna_nearest={'distance_nm':distance,'frame':frame}
        if sum(not o.hide_render for o in stages)!=1: errors.append(f'Caption count at {frame}')
        a=bpy.data.objects['tRNA identity At'].matrix_world@point('At',76,"O3'")
        aa=bpy.data.objects['Peptide residue Ile'].matrix_world.translation
        if abs((aa-a).length-.42)>2e-5: attachment_errors.append([frame,'Ile-to-A-tRNA'])
        if state['peptide_transfer']==0:
            pp=bpy.data.objects['tRNA identity Pt'].matrix_world@point('Pt',76,"O3'")
            pc=bpy.data.objects['Peptide residue Arg'].matrix_world.translation
            if abs((pc-pp).length-.42)>2e-5: attachment_errors.append([frame,'initial-peptide-to-P-tRNA'])
        # Inspect actual rendered geometry during the transfer and final hold.
        t=(frame-1)/(scene.frame_end-1)
        if .34<=t<=.49 or t>=.90:
            dg=bpy.context.evaluated_depsgraph_get()
            front=scene.camera.rotation_euler.to_matrix()@Vector((0,0,1))
            for name in ['Gly','Ala','Arg','Ile']:
                ob=bpy.data.objects['Peptide residue '+name]
                position=ob.matrix_world.translation
                hit=scene.ray_cast(dg,position+front*100,-front)
                if hit[0] and hit[4]!=ob and (hit[1]-position).dot(front)>.20:
                    occluded.append({'frame':frame,'residue':name,'occluder':hit[4].name})
    scene.frame_set(scene.frame_end)
    codon_errors=[]
    for res in [5,6,7]:
        obj=next(o for o in scene.objects if o.get('source_nucleic_chain')=='Mr' and o.get('residue_number')==res)
        final=obj.matrix_world@point('Mr',res,"C4'")
        expected=point('Mr',res-3,"C4'")
        codon_errors.append(float((final-expected).length))
    if max(codon_errors)>2e-5: errors.append('Final mRNA codon displacement is not three residues')
    if rigid_errors: errors.append('Non-rigid tRNA transform')
    if attachment_errors: errors.append('Peptide detaches unexpectedly')
    if cropped: errors.append('Primary tRNA leaves the camera frame')
    if occluded: errors.append('Peptide residues obscured during explanation')
    if entry_mrna_nearest['distance_nm']<.20:
        errors.append('Incoming tRNA crosses mRNA during approach')
    result={'blend':bpy.data.filepath,'frames':scene.frame_end,'intact_tRNAs':intact,
            'intact_context_surface_chains':context_count,'final_codon_position_errors_nm':codon_errors,
            'rigid_errors':rigid_errors,'attachment_errors':attachment_errors,
            'primary_tRNA_cropped_frames':cropped,'peptide_occlusion':occluded,'errors':errors,
            'entry_tRNA_mRNA_nearest_atom_centers':entry_mrna_nearest,
            'scope':'Source geometry, recorded entities, attachment, three-residue transport and direct visibility. Not biochemical kinetics or complete molecular collision certification.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'frames':scene.frame_end,'errors':errors,'occlusion_count':len(occluded),'cropped_count':len(cropped)}))
    if errors: raise RuntimeError('Translation audit failed; see report.')


if __name__=='__main__': main()
