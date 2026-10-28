"""Audit saved hydrolysis atom identities, chemical graphs and synchronized views."""
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
    samples=dict(np.load(Path(scene['motion_directory'])/'hydrolysis.npz'));assets=Path(motion['asset_dir'])
    assert motion['passed'] and motion['kind']=='ras-gap-hydrolysis'
    assert scene.frame_end==count and scene.render.fps==motion['fps']
    assert all(c['maximum_overlap_atoms']==0 for c in motion['full_timeline_checks'].values())
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
    expected_ids={'GTP:'+n for n in np.load(assets/'new_gtp.npz')['atom_names']}|{'water:O','Mg:1'}
    assert set(identities)==expected_ids and len(expected_ids)==34
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
        'stable_chemical_heavy_atoms':34,'max_chemical_coordinate_error_nm':max_error,
        'max_protein_coordinate_errors_nm':protein_errors,'all_source_protein_atoms_retained':True,
        'reactant_components_heavy_atoms':[1,1,32],'product_components_heavy_atoms':[1,5,28],
        'reaction_bond_display_frame_checks':bond_display_checks,'max_bond_radius_error_nm':max_bond_display_error,
        'final_gdp_coordinates_source_matched':True,'pi_rigid_after_formation':True,'captions_checked':True,
        'limits':'Checks saved geometry and bookkeeping, not proton transfer, energetics, kinetics or visual readability.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2));print(json.dumps(result))


if __name__=='__main__':main()
