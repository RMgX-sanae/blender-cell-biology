"""Add synchronized, readable proton traffic to the whole-pump view."""
import json
import math
from pathlib import Path

import bpy
import numpy as np
from mathutils import Vector
from detail_scene import Parts, material, path


def add_main_protons(root):
    root=Path(root)
    scene=bpy.context.scene
    group=bpy.data.collections['Flow']
    info=json.loads((root/'manifest.json').read_text())
    mat=material('H+ motion in both views','F27B36',.36,0)
    shader=mat.node_tree.nodes.get('Principled BSDF')
    shader.inputs['Emission Color'].default_value=mat.diffuse_color
    shader.inputs['Emission Strength'].default_value=.28
    for obj in bpy.data.objects:
        if obj.name.startswith('H transport token '):
            obj.data.materials[0]=mat

    center=np.array([-11.,-2.,0.])
    radius=4.8
    inlet_angle=math.radians(-12)
    carry=math.tau-math.radians(24)
    entry_duration=math.tau*2/16
    exit_duration=math.tau*3.2/16
    def position(a):return center+np.array([radius*math.cos(a),radius*math.sin(a),0])
    inlet=position(inlet_angle)
    outlet=position(math.radians(12))
    incoming=np.array([[-33,-9,5.5],[-23,-9,5.3],[-9,-9,4.3],inlet])
    outgoing=np.array([outlet,[-6.6,-9,-2.7],[-9.5,-10,-5.9],[-17,-10,-10.8]])
    def bezier(points,u):
        v=1-u
        return v*v*v*points[0]+3*v*v*u*points[1]+3*v*u*u*points[2]+u*u*u*points[3]

    for name,points,color,span in [
        ('Whole-pump cytoplasmic entry',incoming,'6CD5D0',(.2,1)),
        ('Whole-pump luminal exit',outgoing,'E9A56D',(0,.78))]:
        guide=material(name+' color',color,.6,0)
        path(name,[bezier(points,u) for u in np.linspace(*span,80)],.045,guide,group)
        u=.52 if 'entry' in name else .61
        p=Vector(bezier(points,u));v=Vector(bezier(points,u+.015)-bezier(points,u-.015)).normalized()
        side=v.cross(Vector((0,-1,0))).normalized()
        path(name+' direction',[p-v*.55+side*.25,p,p-v*.55-side*.25],.055,guide,group)

    particles=[]
    for i,c in enumerate(info['carriers']):
        base=math.atan2(c['site_nm'][1],c['site_nm'][0])
        triple=[]
        for role in ['bound','inlet','outlet']:
            parts=Parts();parts.ball((0,0,0),.34 if i==0 else .29,seg=20,rings=12)
            obj=parts.object('Left H '+role+' '+c['chain'],[mat],group)
            obj['carrier_chain']=c['chain'];obj['role']=role
            obj['evidence']='Proton occupancy/flow graphic; foreground displacement for visibility'
            triple.append(obj)
        particles.append((base,*triple))
    for frame in range(1,386):
        angle=-math.tau*(frame-1)/384
        for base,bound,entry,exit in particles:
            a=base+angle;phase=(inlet_angle-a)%math.tau
            entry_progress=(phase-(math.tau-entry_duration))/entry_duration
            exit_progress=((phase-carry)%math.tau)/exit_duration
            states=[(bound,position(a),phase>carry),
                    (entry,bezier(incoming,np.clip(entry_progress,0,1)),not 0<=entry_progress<=1),
                    (exit,bezier(outgoing,np.clip(exit_progress,0,1)),not 0<=exit_progress<=1)]
            for obj,pos,hide in states:
                obj.location=pos;obj.hide_render=hide
                obj.keyframe_insert('location',frame=frame);obj.keyframe_insert('hide_render',frame=frame)
    for _,*objects in particles:
        for obj in objects:
            action=obj.animation_data.action
            for layer in action.layers:
                for strip in layer.strips:
                    bag=strip.channelbag(obj.animation_data.action_slot)
                    for fc in bag.fcurves:
                        for key in fc.keyframe_points:
                            key.interpolation='CONSTANT' if fc.data_path=='hide_render' else 'LINEAR'
    scene['left_proton_motion']=True
    scene['left_proton_display']='Same c-ring carrier phases as the right view; bound markers offset 2 nm toward viewer. Entry and exit streams have longer visible approach/departure paths; loading/unloading phases match both views.'
    scene['left_proton_pool_ids']='Bound markers represent carrier occupancy; separate inlet/outlet particles represent successive protons in open pools, not a claim of individually measured atom trajectories.'
    scene.frame_set(1)

