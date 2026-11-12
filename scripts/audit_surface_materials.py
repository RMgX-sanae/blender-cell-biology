"""Check actual saved material nodes against the selected artistic role profile.

Run inside Blender on a saved scene. This checks shader handoff, not source
geometry, calibrated optical properties, or the beauty/readability of a render.
"""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np


def check_materials(profile, required=()):
    roles=profile.get('surface_roles',{})
    checked={}
    def equal(actual,expected,message):
        a,b=np.asarray(actual,float),np.asarray(expected,float)
        if a.shape!=b.shape or not np.isfinite(a).all() or not np.isfinite(b).all() or not np.allclose(a,b,atol=1e-7,rtol=1e-6):
            raise ValueError(message)
    for mat in bpy.data.materials:
        role=mat.get('surface_role')
        if role not in roles:continue
        if not mat.use_nodes:raise ValueError('Role material has no node tree: '+mat.name)
        values=json.loads(mat.get('appearance_parameters','{}'))
        for key,value in roles[role].items():
            if key not in values:raise ValueError('Missing resolved setting: '+mat.name+'/'+key)
            equal(values[key],value,'Profile setting was not applied: '+mat.name+'/'+key)
        shader=mat.node_tree.nodes.get('Principled BSDF')
        if shader is None:raise ValueError('Missing role shader: '+mat.name)
        for key,socket in [('roughness','Roughness'),('specular','Specular IOR Level'),
                           ('subsurface_weight','Subsurface Weight'),('subsurface_radius','Subsurface Radius'),
                           ('subsurface_scale','Subsurface Scale')]:
            if shader.inputs[socket].is_linked:raise ValueError('Unexpected driver overrides role setting: '+mat.name+'/'+key)
            equal(shader.inputs[socket].default_value,values[key],'Saved shader differs: '+mat.name+'/'+key)
        outputs=[n for n in mat.node_tree.nodes if n.type=='OUTPUT_MATERIAL' and n.is_active_output]
        if not outputs or any(n.inputs['Displacement'].is_linked for n in outputs):
            raise ValueError('Missing output or geometry displacement in appearance material: '+mat.name)
        reachable=set();pending=[link.from_node for n in outputs for link in n.inputs['Surface'].links]
        while pending:
            node=pending.pop()
            if node in reachable:continue
            reachable.add(node)
            pending.extend(link.from_node for socket in node.inputs for link in socket.links)
        if shader not in reachable:raise ValueError('Role shader does not reach the active surface output: '+mat.name)
        normal=shader.inputs['Normal']
        if values['bump']:
            if len(normal.links)!=1 or normal.links[0].from_node.type!='BUMP' or normal.links[0].from_socket.name!='Normal':
                raise ValueError('Bump does not reach shader normal: '+mat.name)
            bump=normal.links[0].from_node
            if bump.inputs['Distance'].is_linked or bump.inputs['Strength'].is_linked:
                raise ValueError('Unexpected bump setting override: '+mat.name)
            equal(bump.inputs['Distance'].default_value,values['bump'],'Wrong bump distance: '+mat.name)
            equal(bump.inputs['Strength'].default_value,values['bump_strength'],'Wrong bump strength: '+mat.name)
            if len(bump.inputs['Height'].links)!=1 or bump.inputs['Height'].links[0].from_node.type!='TEX_NOISE' or bump.inputs['Height'].links[0].from_socket.identifier!='Fac':
                raise ValueError('Missing surface texture: '+mat.name)
            noise=bump.inputs['Height'].links[0].from_node
            if noise.inputs['Scale'].is_linked or noise.inputs['Detail'].is_linked:
                raise ValueError('Unexpected texture setting override: '+mat.name)
            equal(noise.inputs['Scale'].default_value,values['texture_scale'],'Wrong texture frequency: '+mat.name)
            equal(noise.inputs['Detail'].default_value,values['texture_detail'],'Wrong texture detail: '+mat.name)
            links=noise.inputs['Vector'].links
            if len(links)!=1 or links[0].from_node.type!='TEX_COORD' or links[0].from_socket.name!='Object':
                raise ValueError('Expected local Object texture coordinates: '+mat.name)
        elif normal.is_linked:raise ValueError('Unexpected surface perturbation: '+mat.name)
        checked[mat.name]={'role':role,'parameters':values}
    seen={entry['role'] for entry in checked.values()}
    if not checked or set(required)-seen:
        raise ValueError('Required role materials missing: '+str(set(required)-seen))
    return {'passed':True,'materials':checked,'roles':sorted(seen),
            'limits':'Saved shader values and texture connections only; independently inspect source geometry, colors, lighting, composition and actual rendered pixels.'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--profile',type=Path,required=True)
    p.add_argument('--require-roles',nargs='+',default=[])
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    report=check_materials(json.loads(a.profile.read_text()),a.require_roles)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'passed':True,'materials':len(report['materials']),'roles':report['roles']}))


if __name__=='__main__':main()
