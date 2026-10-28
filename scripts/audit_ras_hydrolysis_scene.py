"""Audit saved Ras hydrolysis/state-cycle identities, graphs and two views."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import bpy
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_phosphotransfer_scene import vertices


def component_sizes(atom_count,edges):
    remaining=set(range(atom_count));sizes=[]
    while remaining:
        reached={min(remaining)}
        while True:
            expanded=reached|{j for i,j in edges if i in reached}|{i for i,j in edges if j in reached}
            if expanded==reached:break
            reached=expanded
        remaining-=reached;sizes.append(len(reached))
    return sorted(sizes)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if a.output.exists():p.error('Use a fresh audit file')
    scene=bpy.context.scene;root=Path(bpy.data.filepath).parent
    m=json.loads((root/'manifest.json').read_text());motion=m['motion'];count=motion['frames']
    cycle=motion['kind']=='ras-state-cycle'
    samples=dict(np.load(Path(scene['motion_directory'])/('state-cycle.npz' if cycle else 'hydrolysis.npz')));assets=Path(motion['asset_dir'])
    assert motion['passed'] and motion['kind'] in ['ras-gap-hydrolysis','ras-state-cycle']
    assert scene.frame_end==count and scene.render.fps==motion['fps']
    assert all(c['maximum_overlap_atoms']==0 for c in motion['full_timeline_checks'].values())
    # Read the saved atom coordinates and transforms without regenerating every
    # display surface. Rendering/visual review checks those surfaces separately.
    for obj in scene.objects:
        for modifier in obj.modifiers:modifier.show_viewport=False
    transforms={}
    for name,t in m['transforms'].items():
        basis=np.array(t['basis']);assert np.allclose(basis.T@basis,np.eye(3),atol=1e-7)
        matrix=np.eye(4);matrix[:3,:3]=basis.T*t['scale']
        matrix[:3,3]=np.array(t['offset'])-np.array(t['center'])*t['scale'];transforms[name]=matrix
    identities=[];actors={}
    for name in m['actors']:
        obj=scene.objects[name];identities.extend(json.loads(obj['atom_ids']))
        assert not any(mod.type=='BOOLEAN' for mod in obj.modifiers)
        bpy.context.window.view_layer=scene.view_layers['Overview render' if obj['view']=='overview' else 'Detail render']
        scene.frame_set(count//2);bpy.context.view_layer.update();evaluated=obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
        assert np.allclose(np.array(evaluated.matrix_world),transforms[obj['view']],atol=1e-5),name
        if obj['view']=='overview':actors[obj['component']]=obj
        else:assert obj.data==scene.objects[obj['component']].data
    expected_ids={('GTP_B:' if cycle else 'GTP:')+n for n in np.load(assets/'new_gtp.npz')['atom_names']}|{'water:O','Mg_B' if cycle else 'Mg:1'}
    if cycle:expected_ids|={'GDP_A:'+n for n in np.load(assets/'old_gdp.npz')['atom_names']}|{'Mg_A'}
    atom_count=63 if cycle else 34
    assert set(identities)==expected_ids and len(expected_ids)==atom_count
    assert set(Counter(identities).values())=={2}
    core_names=motion['core_atom_names'];gamma_names=motion['gamma_atom_names']
    names=core_names+gamma_names+['water:O','Mg:1'];indices={n:i for i,n in enumerate(names)}
    edges=[]
    for key,offset in [('nucleotide_core',0),('phosphate',28)]:
        edges.extend(tuple(sorted((e.vertices[0]+offset,e.vertices[1]+offset))) for e in actors[key].data.edges)
    leaving=tuple(sorted((indices['PG'],indices['O3B'])));forming=(indices['PG'],indices['water:O'])
    reactant_edges=set(edges)|{leaving};product_edges=set(edges)|{forming}
    gtp=np.load(assets/'new_gtp.npz');gdp=np.load(assets/'old_gdp.npz')
    expected_gtp={tuple(sorted((indices[str(gtp['atom_names'][i])],indices[str(gtp['atom_names'][j])]))) for i,j in gtp['bonds']}
    assert reactant_edges==expected_gtp,'Reactant graph differs from the CCD GTP graph'
    expected_gdp={tuple(sorted((indices[str(gdp['atom_names'][i])],indices[str(gdp['atom_names'][j])]))) for i,j in gdp['bonds']}
    assert {edge for edge in product_edges if max(edge)<28}==expected_gdp
    assert component_sizes(34,reactant_edges)==[1,1,32]
    assert component_sizes(34,product_edges)==[1,5,28]
    assert len([e for e in product_edges if indices['PG'] in e])==4
    if cycle:
        old_edges={tuple(sorted(map(int,e.vertices))) for e in actors['old_gdp'].data.edges}
        assert old_edges=={tuple(sorted(map(int,e))) for e in gdp['bonds']},'Old GDP source bonds changed'
        old_edges={tuple(v+34 for v in edge) for edge in old_edges}
        reactant_edges|=old_edges;product_edges|=old_edges
        assert component_sizes(63,reactant_edges)==[1,1,1,28,32]
        assert component_sizes(63,product_edges)==[1,1,5,28,28]
    assert all('Al' not in atom and ':F' not in atom for atom in identities)
    bpy.context.window.view_layer=scene.view_layers['Overview render']
    bonds={name:scene.objects[v['object']] for name,v in m['reaction_bonds'].items()}
    for obj in [*actors.values(),*bonds.values()]:
        for modifier in obj.modifiers:modifier.show_viewport=False
    max_error=0.;max_bond_display_error=0.;checked=0;bond_display_checks=0
    for frame in range(1,count+1):
        scene.frame_set(frame);deps=bpy.context.evaluated_depsgraph_get()
        for key,obj in actors.items():
            actual=vertices(obj);expected=samples[key][frame-1]
            assert actual.shape==expected.shape
            error=float(np.max(np.abs(actual-expected)));max_error=max(error,max_error)
            assert error<2e-4,(frame,key,error)
            evaluated=obj.evaluated_get(deps)
            assert not evaluated.hide_render
            assert np.allclose(np.array(evaluated.matrix_world),transforms['overview'],atol=1e-5)
            checked+=1
        for name,obj in bonds.items():
            kind=m['reaction_bonds'][name]['kind'];value=samples['transfer'][frame-1]
            radius=.026*(1-value if kind=='leaving' else value)
            group=obj.modifiers[0].node_group
            circle=next(n for n in group.nodes if n.bl_idname=='GeometryNodeCurvePrimitiveCircle')
            error=abs(float(circle.inputs['Radius'].default_value)-radius)
            max_bond_display_error=max(error,max_bond_display_error);assert error<1e-7,(frame,name,error)
            second=samples['nucleotide_core'][frame-1,indices['O3B']] if kind=='leaving' else samples['water_oxygen'][frame-1,0]
            expected=np.array([samples['phosphate'][frame-1,0],second])
            assert np.allclose(vertices(obj),expected,atol=2e-4),(frame,name,'Bond endpoint drift')
            copy=scene.objects[m['reaction_bonds'][name]['detail']]
            assert copy.data==obj.data and copy.modifiers[0].node_group==group
            assert np.allclose(np.array(copy.matrix_world),transforms['detail'],atol=1e-5)
            bond_display_checks+=1
    bpy.context.window.view_layer=scene.view_layers['Protein envelopes render']
    protein_errors={}
    for key,name in m['proteins'].items():
        obj=scene.objects[name];source=np.load(assets/(key+'.npz'))
        assert len(obj.data.vertices)==len(source['atoms'])
        assert any(n.bl_idname=='GeometryNodePointsToVolume' for mod in obj.modifiers if mod.type=='NODES' for n in mod.node_group.nodes)
        for modifier in obj.modifiers:modifier.show_viewport=False
        maximum=0.
        for frame in range(1,count+1):
            scene.frame_set(frame);error=float(np.max(np.abs(vertices(obj)-samples[key][frame-1])))
            maximum=max(error,maximum);assert error<2e-4,(frame,key,error)
        protein_errors[key]=maximum
    sos_error=0.
    if cycle:
        sos=scene.objects[m['sos_object']];source=np.load(assets/'sos.npz')
        actual=np.array([list(v.co) for v in sos.data.vertices])
        assert np.allclose(actual,source['vertices'],atol=1e-6),'SOS source mesh changed'
        for frame in range(1,count+1):
            scene.frame_set(frame);observed=np.array(sos.evaluated_get(bpy.context.evaluated_depsgraph_get()).matrix_world)
            expected=transforms['overview']@samples['sos_pose'][frame-1]
            error=float(np.max(np.abs(observed-expected)));sos_error=max(error,sos_error)
            assert error<1e-4,(frame,'SOS transform',error)
        start=np.load(assets/'ras.npz')['atoms']
        assert np.allclose(samples['ras'][0],start,atol=1e-6) and np.allclose(samples['ras'][-1],start,atol=1e-6)
        assert max(motion['handoff_coordinate_errors_nm'].values())<1e-6
        matrix=transforms['overview'];scale=m['transforms']['overview']['scale']
        for frame in [0,count]:
            for key,points in [('gap',samples['gap'][frame]),('sos',source['vertices']@samples['sos_pose'][frame,:3,:3].T+samples['sos_pose'][frame,:3,3])]:
                projected=points@matrix[:3,:3].T+matrix[:3,3];lo=projected.min(0)-.30*scale;hi=projected.max(0)+.30*scale
                assert hi[0]<-11 or lo[0]>0 or hi[1]<-3.52 or lo[1]>3.52,(frame,key,'Waiting/departed protein leaks into overview')
        # An empty nucleotide-binding site must remain visually empty in both
        # panels around the middle of the exchange stage, despite the new view.
        empty=motion['exchange_frames']//2
        for view,matrix in transforms.items():
            x0,x1=(-11,0) if view=='overview' else (0,11);scale=m['transforms'][view]['scale']
            for key in actors:
                points=samples[key][empty]@matrix[:3,:3].T+matrix[:3,3];lo=points.min(0)-.09*scale;hi=points.max(0)+.09*scale
                assert hi[0]<x0 or lo[0]>x1 or hi[1]<-3.52 or lo[1]>3.52,(view,key,'Empty state shows chemical actor')
        matrix=transforms['detail'];scale=m['transforms']['detail']['scale']
        water=samples['water_oxygen'][:motion['exchange_frames']+1,0]@matrix[:3,:3].T+matrix[:3,3]
        radius=.055*scale;inside=(water[:,0]+radius>=0)&(water[:,0]-radius<=11)&(water[:,1]+radius>=-3.52)&(water[:,1]-radius<=3.52)
        assert not np.any(inside&(water[:,2]+radius>=scene.camera.location.z-scene.camera.data.clip_start)), 'Water appears through the camera near plane'
    # Check final identity-preserving product coordinates against the independent
    # prepared GDP source, not just agreement with the animated cache.
    lookup={str(n):i for i,n in enumerate(gdp['atom_names'])}
    assert np.allclose(samples['nucleotide_core'][-1],gdp['atoms'][[lookup[n] for n in core_names]],atol=2e-4)
    assert np.allclose(samples['magnesium'][-1],gdp['magnesium'],atol=1e-7)
    product=np.concatenate([samples['phosphate'],samples['water_oxygen']],axis=1)
    after=np.flatnonzero(samples['transfer']>=1-1e-8)
    bonds_length=np.linalg.norm(product[after,1:]-product[after,:1],axis=2)
    assert np.max(np.ptp(bonds_length,axis=0))<1e-6,'Pi product distorts during release'
    assert np.allclose(bonds_length[0],[.153,.153,.153,.155],atol=1e-6)
    bpy.context.window.view_layer=scene.view_layers['Typography render']
    caption_bodies={row[2] for row in m['captions']}
    for frame in range(1,count+1):
        scene.frame_set(frame);deps=bpy.context.evaluated_depsgraph_get()
        visible=[o for o in scene.objects if o.type=='FONT' and o.data.body in caption_bodies and not o.evaluated_get(deps).hide_render]
        assert len(visible)==1,(frame,'caption overlap or absence')
    result={'passed':True,'frames':count,'chemical_actor_frame_checks':checked,
        'stable_chemical_heavy_atoms':atom_count,'max_chemical_coordinate_error_nm':max_error,
        'max_protein_coordinate_errors_nm':protein_errors,'all_source_protein_atoms_retained':True,
        'reactant_components_heavy_atoms':component_sizes(atom_count,reactant_edges),'product_components_heavy_atoms':component_sizes(atom_count,product_edges),
        'reaction_bond_display_frame_checks':bond_display_checks,'max_bond_radius_error_nm':max_bond_display_error,
        'final_gdp_coordinates_source_matched':True,'pi_rigid_after_formation':True,'captions_checked':True,
        'integrated_state_cycle':cycle,'max_sos_transform_error':sos_error if cycle else None,
        'empty_exchange_state_checked':cycle,
        'limits':'Checks saved geometry and bookkeeping, not proton transfer, energetics, kinetics or visual readability.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2));print(json.dumps(result))


if __name__=='__main__':main()
