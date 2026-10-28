"""Audit reopened Ras exchange geometry, atom identity, timing and two views."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import bpy
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_phosphotransfer_scene import vertices


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.output.exists():p.error('Use a fresh audit file')
    scene=bpy.context.scene;root=Path(bpy.data.filepath).parent
    manifest=json.loads((root/'manifest.json').read_text());motion=manifest['motion']
    samples=np.load(Path(scene['motion_directory'])/'exchange.npz');asset=Path(motion['asset_dir'])
    assert motion['passed'] and 'full_timeline_checks' in motion
    assert motion['maximum_ras_connection_nm']<.22 and motion['minimum_ras_connection_nm']>.10
    assert scene.frame_end==motion['frames'] and scene.render.fps==motion['fps']
    transforms={}
    for name,t in manifest['transforms'].items():
        matrix=np.eye(4);matrix[:3,:3]=np.array(t['basis']).T*t['scale']
        matrix[:3,3]=np.array(t['offset'])-np.array(t['center'])*t['scale'];transforms[name]=matrix
    identities=[];actors={}
    for name in manifest['actors']:
        obj=scene.objects[name];identities.extend(json.loads(obj['atom_ids']))
        assert not any(m.type=='BOOLEAN' for m in obj.modifiers)
        if obj['view']=='overview':actors[obj['component']]=obj
        else:assert obj.data==scene.objects[obj['component']].data
    assert set(Counter(identities).values())=={2}
    assert len(set(identities))==motion['heavy_atoms_in_nucleotides_and_metals']==62
    bpy.context.window.view_layer=scene.view_layers['Overview render']
    for obj in actors.values():
        for modifier in obj.modifiers:modifier.show_viewport=False
    chemical_error=0.;checked=0;bond_error=0.;nucleotides={}
    for key in ['old_gdp','new_gtp']:
        d=np.load(asset/(key+'.npz'));nucleotides[key]=d
        observed={tuple(sorted(e.vertices)) for e in actors[key].data.edges}
        expected={tuple(sorted(map(int,e))) for e in d['bonds']}
        assert observed==expected,(key,'Saved chemical bonds differ from CCD preparation')
        reached={0}
        for _ in d['atom_names']:
            reached.update(j for i,j in expected if i in reached)
            reached.update(i for i,j in expected if j in reached)
        assert len(reached)==len(d['atoms']),(key,'Disconnected heavy-atom graph')
    for frame in range(1,motion['frames']+1):
        scene.frame_set(frame)
        for key,obj in actors.items():
            observed=vertices(obj);expected=samples[key][frame-1]
            assert observed.shape==expected.shape
            error=float(np.max(np.abs(observed-expected)));chemical_error=max(chemical_error,error)
            assert error<2e-4,(frame,key,error)
            evaluated=obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
            assert np.allclose(np.array(evaluated.matrix_world),transforms['overview'],atol=1e-5)
            assert not evaluated.hide_render;checked+=1
            if key in nucleotides:
                d=nucleotides[key];bonds=d['bonds']
                baseline=np.linalg.norm(d['atoms'][bonds[:,0]]-d['atoms'][bonds[:,1]],axis=1)
                lengths=np.linalg.norm(observed[bonds[:,0]]-observed[bonds[:,1]],axis=1)
                bond_error=max(bond_error,float(np.max(np.abs(lengths-baseline))))
    assert bond_error<2e-4,('Nucleotide bond lengths changed',bond_error)
    bpy.context.window.view_layer=scene.view_layers['Protein envelopes render']
    ras=scene.objects['Ras G-domain'];sos=scene.objects['SOS catalytic region'];source=np.load(asset/'sos.npz')
    assert len(ras.data.vertices)==1322
    assert any(n.bl_idname=='GeometryNodePointsToVolume' for mod in ras.modifiers if mod.type=='NODES' for n in mod.node_group.nodes)
    original=np.array([list(v.co) for v in sos.data.vertices])
    assert np.allclose(original,source['vertices'],atol=1e-6),'SOS source mesh was distorted'
    for modifier in ras.modifiers:modifier.show_viewport=False
    ras_error=0.;sos_error=0.
    for frame in range(1,motion['frames']+1):
        scene.frame_set(frame)
        error=float(np.max(np.abs(vertices(ras)-samples['ras'][frame-1])));ras_error=max(ras_error,error)
        assert error<2e-4,(frame,'Ras',error)
        expected=transforms['overview']@samples['sos_pose'][frame-1]
        actual=np.array(sos.evaluated_get(bpy.context.evaluated_depsgraph_get()).matrix_world)
        error=float(np.max(np.abs(actual-expected)));sos_error=max(sos_error,error)
        assert error<1e-4,(frame,'SOS',error)
    bpy.context.window.view_layer=scene.view_layers['Typography render']
    captions={row[2] for row in manifest['captions']}
    for frame in range(1,motion['frames']+1):
        scene.frame_set(frame);deps=bpy.context.evaluated_depsgraph_get()
        visible=[o for o in scene.objects if o.type=='FONT' and o.data.body in captions and not o.evaluated_get(deps).hide_render]
        assert len(visible)==1,(frame,'caption overlap or absence')
    empty=(motion['intervals'][0][1]+motion['intervals'][1][0])/2
    empty_frame=round(empty*motion['frames'])
    for view,matrix in transforms.items():
        x0,x1=(-11,0) if view=='overview' else (0,11)
        scale=manifest['transforms'][view]['scale']
        for key in actors:
            points=samples[key][empty_frame]@matrix[:3,:3].T+matrix[:3,3];padding=.09*scale
            lo=points.min(0)-padding;hi=points.max(0)+padding
            assert hi[0]<x0 or lo[0]>x1 or hi[1]<-3.52 or lo[1]>3.52,(view,key,'Nucleotide-free state still shows ligand')
    result={'passed':True,'frames':motion['frames'],'chemical_actor_frame_checks':checked,
        'stable_chemical_heavy_atoms':62,'max_chemical_coordinate_error_nm':chemical_error,
        'max_ras_coordinate_error_nm':ras_error,'max_sos_transform_error':sos_error,
        'sos_mesh_preserved':True,'ras_density_from_all_atoms':True,'nucleotide_free_views_checked':True,
        'captions_checked':True,'chemical_bonds_unchanged':True,'nucleotide_graphs_connected':True,'maximum_bond_length_error_nm':bond_error,
        'limits':'Source-guided teaching geometry; does not validate solvent, kinetics, biological generality or visual readability.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2))
    print(json.dumps(result))


if __name__=='__main__':main()
