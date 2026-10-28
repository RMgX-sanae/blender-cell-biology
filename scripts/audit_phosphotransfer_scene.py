"""Verify the reopened PKA scene against its independently prepared atom motion."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import bpy
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))


def vertices(obj):
    mesh=obj.evaluated_get(bpy.context.evaluated_depsgraph_get()).to_mesh()
    values=np.empty(len(mesh.vertices)*3,dtype=np.float32);mesh.vertices.foreach_get('co',values)
    obj.evaluated_get(bpy.context.evaluated_depsgraph_get()).to_mesh_clear()
    return values.reshape(-1,3)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.output.exists():p.error('Use a new audit file')
    scene=bpy.context.scene;root=Path(bpy.data.filepath).parent
    m=json.loads((root/'manifest.json').read_text());motion=m['motion'];count=motion['frames_per_cycle']
    samples=np.load(Path(scene['motion_directory'])/'cycle.npz');assets=Path(motion['asset_dir'])
    assert scene.frame_end==motion['frames'] and scene.render.fps==motion['fps']
    assert motion['passed'],'Motion preflight failed'
    groups={i:{} for i in range(1,motion['cycles']+1)};identities=[]
    transforms={}
    for name,t in m['transforms'].items():
        matrix=np.eye(4);matrix[:3,:3]=np.array(t['basis']).T*t['scale']
        matrix[:3,3]=np.array(t['offset'])-np.array(t['center'])*t['scale'];transforms[name]=matrix
    for name in m['actors']:
        obj=scene.objects[name]
        bpy.context.window.view_layer=scene.view_layers['Overview render' if obj['view']=='overview' else 'Detail render']
        scene.frame_set((obj['packet']-1)*count+count//2+1);bpy.context.view_layer.update()
        evaluated=obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
        assert np.allclose(np.array(evaluated.matrix_world),transforms[obj['view']],atol=1e-5),name
        assert not any(mod.type=='BOOLEAN' for mod in obj.modifiers)
        identities.extend(json.loads(obj['atom_ids']))
        if obj['view']=='overview':groups[obj['packet']][obj['component']]=obj
        else:
            original=scene.objects[name.removesuffix('_detail')]
            assert obj.data==original.data,'Views do not share atom coordinates'
    expected_atoms=m['sources']['chemical_identity']['heavy_atoms_initial']
    assert expected_atoms==m['sources']['chemical_identity']['heavy_atoms_final']
    assert set(Counter(identities).values())=={2},'Each stable atom must appear once in each view'
    assert len(set(identities))==expected_atoms*motion['cycles']
    bpy.context.window.view_layer=scene.view_layers['Overview render']
    for packet in groups.values():
        for obj in packet.values():
            for modifier in obj.modifiers:modifier.show_viewport=False
    maximum=0.;checked=0;active_counts=[]
    for frame in range(1,motion['frames']+1):
        scene.frame_set(frame);index=(frame-1)%count;packet_index=(frame-1)//count+1
        active=0
        for number,packet in groups.items():
            for key,obj in packet.items():
                evaluated=obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
                assert evaluated.hide_render==(number!=packet_index),'Packet reuse/reset visibility error'
                if number!=packet_index:continue
                observed=vertices(obj);expected=samples[key][index]
                assert observed.shape==expected.shape
                error=float(np.max(np.abs(observed-expected)));maximum=max(maximum,error)
                # Float32 absolute-key time conversion adds up to about
                # 0.0001 nm here. This tolerance is separate from contact QC.
                assert error<2e-4,(frame,obj.name,error)
                assert np.all(np.isfinite(observed));active+=len(observed);checked+=1
        assert active==expected_atoms;active_counts.append(active)
    bpy.context.window.view_layer=scene.view_layers['Enzyme envelope render']
    enzyme=scene.objects['PKAc intact overview'];d=np.load(assets/'kinase.npz')
    assert len(enzyme.data.vertices)==len(d['atoms'])
    assert any(n.bl_idname=='GeometryNodePointsToVolume' for mod in enzyme.modifiers if mod.type=='NODES' for n in mod.node_group.nodes)
    for modifier in enzyme.modifiers:modifier.show_viewport=False
    states=np.array([d[k] for k in ['atoms','product_atoms','open_atoms','bound_open_atoms','product_bound_open_atoms']])
    protein_error=0.;longest_connection=0.;bonds=d['bonds']
    for frame in range(1,motion['frames']+1):
        scene.frame_set(frame);index=(frame-1)%count
        expected=np.einsum('s,snc->nc',samples['enzyme_weights'][index],states)
        coordinates=vertices(enzyme)
        error=float(np.max(np.abs(coordinates-expected)));protein_error=max(protein_error,error)
        assert error<2e-5,(frame,'enzyme geometry',error)
        longest_connection=max(longest_connection,float(np.linalg.norm(coordinates[bonds[:,0]]-coordinates[bonds[:,1]],axis=1).max()))
    assert longest_connection<.22,('Excessively stretched kinase connection',longest_connection)
    bpy.context.window.view_layer=scene.view_layers['Typography render']
    captions={label for _,_,label in __import__('build_phosphotransfer_scene').CAPTIONS}
    for frame in range(1,motion['frames']+1):
        scene.frame_set(frame)
        depsgraph=bpy.context.evaluated_depsgraph_get()
        visible=[obj for obj in scene.objects if obj.type=='FONT' and obj.data.body in captions and not obj.evaluated_get(depsgraph).hide_render]
        assert len(visible)==1,(frame,'overlapping or absent stage captions')
    result={'passed':True,'source_blend':bpy.data.filepath,'frames_checked':motion['frames'],
        'actor_frame_checks':checked,'max_atom_coordinate_error_nm':maximum,'max_enzyme_atom_error_nm':protein_error,
        'atom_coordinate_tolerance_nm':2e-4,
        'stable_atoms_per_packet':expected_atoms,'independent_packets':motion['cycles'],
        'views_share_coordinates':True,'all_source_enzyme_atoms_preserved':True,'envelope_rebuilt_from_atoms':True,'stage_captions_checked':True,
        'maximum_enzyme_connection_nm':longest_connection,
        'limitations':'Coordinate/identity checks do not validate molecular kinetics, biological generality or visual readability.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2))
    print(json.dumps(result,ensure_ascii=False))


if __name__=='__main__':main()
