"""Reopen a generated style scene and check source coordinates and live loop state.

Run inside Blender with a saved .blend loaded. This is a technical audit; visual
reference comparison and biological interpretation remain separate.
"""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    scene=bpy.context.scene
    errors=[];checked=0
    for obj in scene.objects:
        if obj.type=='MESH':
            points=np.array([v.co[:] for v in obj.data.vertices])
            if not np.isfinite(points).all():errors.append('Nonfinite geometry: '+obj.name)
            if args.asset_dir and obj.get('element'):
                source=np.load(args.asset_dir/f"chain_{obj['source_chain']}.npz",allow_pickle=False)
                expected=source['atoms'][source['elements']==obj['element']]
                if points.shape!=expected.shape or not np.allclose(points,expected,atol=1e-6,rtol=0):
                    errors.append('Source atom coordinates changed: '+obj.name)
                checked+=len(points)
    missing=[]
    grain=[]
    for im in bpy.data.images:
        if im.name.startswith('Generated fixed compositor grain'):
            pixels=np.array(im.pixels[:]).reshape(-1,4)[:,:3]
            grain.append({'name':im.name,'mean':float(pixels.mean()),'std':float(pixels.std()),'packed':bool(im.packed_file)})
            if pixels.mean()<.8 or not im.packed_file:errors.append('Missing or black grain: '+im.name)
        if im.source=='FILE' and not im.packed_file and im.filepath and not Path(bpy.path.abspath(im.filepath)).exists():
            missing.append(im.name)
    errors.extend('Missing external image: '+x for x in missing)
    trees=[t for t in bpy.data.node_groups if t.bl_idname in {'GeometryNodeTree','CompositorNodeTree'}]
    invalid=[t.name for t in trees if any(not link.is_valid for link in t.links)]
    errors.extend('Invalid node link: '+x for x in invalid)
    length=scene.frame_end-scene.frame_start+1
    frames=[1,1+length//4,1+length//2,1+3*length//4,scene.frame_end+1]
    reference_points=None
    if args.asset_dir:
        chains=sorted({obj.get('source_chain') for obj in scene.objects if obj.get('source_chain')})
        if chains:
            reference_points=np.concatenate([np.load(args.asset_dir/f'chain_{c}.npz',allow_pickle=False)['vertices'] for c in chains])
    states=[]
    margins=[]
    for frame in frames:
        scene.frame_set(frame);bpy.context.view_layer.update()
        root=bpy.data.objects.get('Rigid style-preview motion')
        state=[]
        if root is not None:
            state.extend(np.array(root.matrix_world).ravel())
            if reference_points is not None and scene.camera.data.type=='ORTHO':
                homogeneous=np.column_stack((reference_points,np.ones(len(reference_points))))
                projected=homogeneous@np.array(root.matrix_world).T@np.array(scene.camera.matrix_world.inverted()).T
                width=scene.camera.data.ortho_scale
                height=width*scene.render.resolution_y/scene.render.resolution_x
                x=.5+projected[:,0]/width;y=.5+projected[:,1]/height
                margin=float(min(x.min(),y.min(),1-x.max(),1-y.max()))
                margins.append({'frame':frame,'minimum_fractional_margin':margin})
                if margin<.01:errors.append(f'Subject reaches frame edge at frame {frame}.')
        deps=bpy.context.evaluated_depsgraph_get()
        particles=[]
        for instance in deps.object_instances:
            if instance.is_instance and instance.parent and instance.parent.original.name=='Intracellular puncta':
                particles.append(instance.matrix_world.translation[:])
        if particles:
            points=np.array(particles)
            state.extend(points.mean(0))
            state.append(len(points))
        states.append(np.array(state))
    has_motion=bool(len(states[0]))
    expected_motion=bpy.data.objects.get('Rigid style-preview motion') is not None or scene.get('style_id')=='fluorescence'
    seam=float(np.max(abs(states[0]-states[-1]))) if has_motion else 0.
    movement=float(np.max(abs(states[0]-states[1]))) if has_motion else 0.
    if seam>1e-5:errors.append('Periodic state does not close.')
    if expected_motion and movement<1e-5:errors.append('Style preview expected motion but remains static.')
    scene.frame_set(1)
    report={'style':scene.get('style_id'),'blender':bpy.app.version_string,'saved_scene':bpy.data.filepath,
            'source_atoms_checked':checked,'node_groups':len(trees),'grain':grain,
            'loop_state_seam_max':seam,'quarter_cycle_state_change':movement,'motion_state_present':has_motion,'motion_expected':expected_motion,
            'frames_checked':frames,'camera_margins':margins,'errors':errors,'passed':not errors,
            'limits':'Technical invariants only; inspect actual rendered images and decoded video separately.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    if errors:raise SystemExit(1)


if __name__=='__main__':main()
