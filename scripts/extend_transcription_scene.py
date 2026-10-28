"""Slow a saved one-C scene and add the subsequent U cycle and authored translocation.

Run with Blender -b first-cycle/scene.blend --python this_file -- ... .
Preserves the original lighting/materials/compositor. Output is a fresh scene.
"""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Matrix, Quaternion, Vector

sys.path.insert(0,str(Path(__file__).resolve().parent))
from build_transcription_scene import subset
from build_molecular_scene import nucleic_atoms
from build_translation_scene import atom, smooth, text_object, update_line
from detail_scene import path, save_and_render


def curves(action):
    if hasattr(action,'fcurves'):yield from action.fcurves
    for layer in getattr(action,'layers',[]):
        for strip in layer.strips:
            for bag in getattr(strip,'channelbags',[]):yield from bag.fcurves


def visibility(ob,visible,frame):
    ob.hide_render=ob.hide_viewport=not visible
    ob.keyframe_insert(data_path='hide_render',frame=frame)
    ob.keyframe_insert(data_path='hide_viewport',frame=frame)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path,required=True)
    p.add_argument('--substrate-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--frame',type=int,default=457)
    p.add_argument('--render',action='store_true')
    p.add_argument('--animate',action='store_true')
    args=p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    if args.output.exists():p.error('Choose a fresh output directory.')
    source_blend=bpy.data.filepath
    old=json.loads(Path(source_blend).with_name('manifest.json').read_text())
    info=json.loads((args.substrate_dir/'manifest.json').read_text())
    if 'next_cycle' not in info:p.error('Run prepare_transcription_next.py first.')
    if old['source_asset']['source_sha256']!=info['source_asset']['source_sha256']:
        p.error('The saved first cycle and new substrate must use the same source.')
    s=bpy.context.scene
    args.frames,args.fps,args.width,args.samples=552,24,old['width'],old['samples']
    s.frame_start,s.frame_end,s.render.fps=1,args.frames,args.fps
    if not 1<=args.frame<=args.frames:p.error('Preview frame outside animation.')
    # Retiming all original actions keeps the first cycle's authored geometry.
    old_times=[0,.12,.40,.43,.58,.84,1]
    seconds=[0,1,4,5,6.5,9.5,11]
    for action in list(bpy.data.actions):
        for fc in curves(action):
            for key in fc.keyframe_points:
                for value in [key.co,key.handle_left,key.handle_right]:
                    value.x=1+args.fps*np.interp((value.x-1)/(old['frames']-1),old_times,seconds)
            fc.update()
    scaffold=bpy.data.collections['02 Native modeled DNA and RNA']
    reaction=bpy.data.collections['03 One-C addition — explicit teaching reconstruction']
    annotation=bpy.data.collections['04 Captions and polarity']
    reaction.name='03 Two successive additions — teaching reconstruction'
    dna=np.load(args.asset_dir/'chain_T.npz');rna=np.load(args.asset_dir/'chain_P.npz')
    utp=np.load(args.substrate_dir/'utp.npz');count=int(utp['product_atom_count'])
    template=bpy.data.objects['Template DNA — all modeled residues']
    dna_mat=template.data.materials[0]
    bpy.data.objects.remove(template,do_unlink=True)
    template=nucleic_atoms(subset(dna,dna['residue_numbers']!=16),'Template DNA except complete A16',dna_mat,scaffold)
    a16=nucleic_atoms(subset(dna,dna['residue_numbers']==16),'Complete template A16 — modeled loading pose',dna_mat,scaffold)
    a16['evidence']='Native nucleotide rigidly repositioned; no atoms removed or sliced.'
    parent=bpy.data.objects.new('Authored one-position scaffold translocation',None)
    scaffold.objects.link(parent)
    core1=bpy.data.objects['Incoming CTP core → RNA C20']
    ppi1=bpy.data.objects['CTP beta/gamma phosphates → PPi']
    for ob in [template,a16,bpy.data.objects['Non-template DNA — hidden whole chain'],
               bpy.data.objects['RNA primer P residues 7–19'],core1,
               bpy.data.objects['New RNA O3′–P bond']]:ob.parent=parent
    # Recalculate anchors in world space; rotating a label's foreground offset
    # with the molecule would tilt its text and displace it far from the anchor.
    polarity=[('DNA 模板链 3′',atom(dna,39,"C4'"),(.9,.8)),
              ('5′',atom(dna,1,"C4'"),(.9,.4)),
              ('RNA 片段 5′',atom(rna,7,"C4'"),(.7,.7))]
    link1=path('Template T15 O3′ to T16 P',[(0,0,0)]*2,.047,dna_mat,scaffold)
    link2=path('Template T16 O3′ to T17 P',[(0,0,0)]*2,.047,dna_mat,scaffold)
    link1.parent=link2.parent=parent
    front=np.array(s.camera.rotation_euler.to_matrix() @ Vector((0,0,1)))
    right=np.array(s.camera.rotation_euler.to_matrix() @ Vector((1,0,0)))
    up=np.array(s.camera.rotation_euler.to_matrix() @ Vector((0,1,0)))
    trans_r,trans_s=utp['translocation_rotation'],utp['translocation_translation']
    local_r,local_s=utp['template_rotation'],utp['template_translation']
    trans_q=Matrix(trans_r).to_quaternion();local_q=Matrix(local_r).to_quaternion()
    origin=atom(dna,16,"C1'");delta=local_r @ origin+local_s-origin
    identity=Quaternion((1,0,0,0))
    camera_start=np.array(s.camera.location)
    camera_scale=s.camera.data.ortho_scale
    parent.rotation_mode=a16.rotation_mode='QUATERNION'
    world_atoms=utp['atoms'] @ trans_r.T+trans_s
    udata={**dict(utp),'atoms':world_atoms}
    core2=nucleic_atoms(subset(udata,np.arange(len(world_atoms))<count),
        'Incoming UTP core → modeled RNA U21',core1.data.materials[0],reaction)
    ppi2=nucleic_atoms(subset(udata,np.arange(len(world_atoms))>=count),
        'UTP beta/gamma phosphates → second PPi',ppi1.data.materials[0],reaction)
    core2['evidence']='Modeled U21; not deposited in 5FLM.'
    ppi2['evidence']='CCD UTP beta/gamma phosphate atoms; charge/protonation omitted.'
    # Once the first addition is complete, C20 takes the existing RNA color.
    # Gold then continues to identify the nucleotide currently being added.
    incorporated=core1.data.materials[0].copy();incorporated.name='Incorporated C20 becomes RNA orange'
    core1.data.materials[0]=incorporated
    base_color=incorporated.node_tree.nodes.get('Principled BSDF').inputs['Base Color']
    gold=np.array(base_color.default_value)
    rna_color=np.array(bpy.data.objects['RNA primer P residues 7–19'].data.materials[0].node_tree.nodes.get('Principled BSDF').inputs['Base Color'].default_value)
    names=list(utp['atom_names']);base=['N1','C2','N3','C4','C5','C6','O2','O4']
    a16_data=subset(dna,dna['residue_numbers']==16)
    a16_final=(a16_data['atoms'] @ local_r.T+local_s) @ trans_r.T+trans_s
    outward=world_atoms[[names.index(n) for n in base]].mean(0)-a16_final.mean(0)
    outward/=np.linalg.norm(outward)
    pa=world_atoms[names.index('PA')];o3a=world_atoms[names.index('O3A')]
    c20_o3=atom(rna,20,"O3'");u21_o3=world_atoms[names.index("O3'")]
    c20_final=trans_r @ c20_o3+trans_s
    addition2=path('Second RNA O3′–P bond',[(0,0,0)]*2,.035,
        bpy.data.objects['New RNA O3′–P bond'].data.materials[0],reaction)
    cleavage2=path('Second alpha–beta phosphate connection',[(0,0,0)]*2,.035,ppi2.data.materials[0],reaction)
    terminal=bpy.data.objects['RNA growing 3′ end']
    # Retain the first-cycle 3′ tracking until it reaches C20.
    original_terminal=[]
    for frame in range(1,266):
        s.frame_set(frame);original_terminal.append(tuple(terminal.location))
    terminal.animation_data_clear()
    ppi1_extra=np.array(right*3.2-up*.9+front*1.8)
    poses=[]
    for frame in range(1,args.frames+1):
        time=(frame-1)/args.fps
        advance=smooth(time,11,13.5)
        base_color.default_value=(1-advance)*gold+advance*rna_color
        base_color.keyframe_insert(data_path='default_value',frame=frame)
        q=identity.slerp(trans_q,advance);r=np.array(q.to_matrix());shift=trans_s*advance
        parent.rotation_quaternion=q;parent.location=shift
        parent.keyframe_insert(data_path='rotation_quaternion',frame=frame)
        parent.keyframe_insert(data_path='location',frame=frame)
        s.camera.location=camera_start+up*.30*advance
        s.camera.data.ortho_scale=camera_scale+1.5*advance
        s.camera.keyframe_insert(data_path='location',frame=frame)
        s.camera.data.keyframe_insert(data_path='ortho_scale',frame=frame)
        for body,point,(dx,dy) in polarity:
            anchor=r @ point+shift
            label=bpy.data.objects['Polarity '+body]
            label.location=anchor+right*dx+up*dy+front*20
            label.keyframe_insert(data_path='location',frame=frame)
            update_line(bpy.data.objects['Leader '+body],
                [anchor+front*20,np.array(label.location)-right*.1],frame)
        aq=identity.slerp(local_q,advance);ar=np.array(aq.to_matrix());ashift=origin-ar @ origin+delta*advance
        a16.rotation_quaternion=aq;a16.location=ashift
        a16.keyframe_insert(data_path='rotation_quaternion',frame=frame)
        a16.keyframe_insert(data_path='location',frame=frame)
        update_line(link1,[atom(dna,15,"O3'"),ar @ atom(dna,16,'P')+ashift],frame)
        update_line(link2,[ar @ atom(dna,16,"O3'")+ashift,atom(dna,17,'P')],frame)
        docking=smooth(time,13.5,16.5);joining=smooth(time,17.5,19);release=smooth(time,19,21.5)
        approach=(1-docking)*(right*9+up*1.6+front*4)
        hold=(front*.30+outward*.30)*(1-joining)
        core2.location=approach+hold
        ppi2.location=approach+hold+release*(right*3.2-up*.9+front*1.8)
        for ob in [core2,ppi2]:
            ob.keyframe_insert(data_path='location',frame=frame)
            visibility(ob,time>=13.5,frame)
        update_line(cleavage2,[pa+np.array(core2.location),o3a+np.array(ppi2.location)],frame)
        update_line(addition2,[c20_final,pa+np.array(core2.location)],frame)
        for ob,value in [(cleavage2,1-joining),(addition2,joining)]:
            ob.data.bevel_depth=.035*value;ob.data.keyframe_insert(data_path='bevel_depth',frame=frame)
            visibility(ob,time>=13.5,frame)
        if time<11:terminal.location=original_terminal[frame-1]
        elif time<17.5:terminal.location=r @ c20_o3+shift+right*.7-up*.75+front*20
        else:terminal.location=(1-joining)*c20_final+joining*u21_o3+right*.7-up*.75+front*20
        terminal.keyframe_insert(data_path='location',frame=frame)
        # The first PPi leaves the view before the next substrate arrives.
        if time>=9.5:
            ppi1.location=ppi1_extra+right*10*smooth(time,9.5,11)
            ppi1.keyframe_insert(data_path='location',frame=frame)
        visibility(ppi1,time<11,frame)
        poses.append({'frame':frame,'seconds':time,'advance':advance,'second_docking':docking,
            'second_joining':joining,'second_release':release,
            'scaffold_rotation_wxyz':list(q),'scaffold_translation_nm':list(shift),
            'A16_rotation_wxyz':list(aq),'A16_translation_nm':list(ashift),
            'U_core_translation_nm':list(core2.location),'second_PPi_translation_nm':list(ppi2.location)})
    # Replace captions rather than leave the first scene's final caption on screen.
    old_stage=[ob for ob in list(annotation.objects) if ob.name.startswith('Stage ') or
               ob.name.startswith('Outline Stage ') or ob.name in ['Substrate accounting','Outline Substrate accounting']]
    for ob in old_stage:bpy.data.objects.remove(ob,do_unlink=True)
    footer=bpy.data.objects['Presentation scope'];font=footer.data.font
    ink=footer.data.materials[0];halo=bpy.data.materials['Annotation halo']
    footer.data.body='完整分子分层展示｜第一轮 C20 来自 5FLM；转位与第二轮 U21 为教学重构'
    bpy.data.objects['Outline Presentation scope'].data.body=footer.data.body
    stages=[(0,1,'第一轮：RNA 的 3′端等待 CTP'),(1,5,'第一轮：CTP 对准模板 G'),
        (5,6.5,'第一轮：连接变化示意，加入 C'),(6.5,9.5,'第一轮：焦磷酸 PPi 离开'),
        (9.5,11,'第一轮完成：C20 产物坐标来自 5FLM'),
        (11,13.5,'转位示意：模板与 RNA 前进一个位置'),
        (13.5,17.5,'第二轮：UTP 对准下一个模板 A'),
        (17.5,19,'第二轮：加入 U，RNA 继续在 3′端延长'),
        (19,21.5,'第二轮：焦磷酸 PPi 离开'),
        (21.5,23.1,'已连续加入 C 和 U；第二轮构象为教学重构')]
    new_labels=[]
    for i,(start,end,body) in enumerate(stages):
        ob=text_object(body,f'Cycle stage {i+1}',.25,ink,annotation,font)
        ob.parent=s.camera;ob.location=(-6.65,3.99,-15)
        for state in poses:visibility(ob,start<=state['seconds']<end,state['frame'])
        new_labels.append(ob)
    formula=text_object('RNA + CTP + UTP → RNA 增加 C、U + 2 PPi','Two-cycle accounting',.23,ink,annotation,font)
    formula.parent=s.camera;formula.location=(-3,-3.3,-15);new_labels.append(formula)
    for name,obj,points,kind in [('UTP',core2,world_atoms,'substrate'),('PPi',ppi2,world_atoms[count:],'product')]:
        label=text_object(name,'Second '+name+' label',.20,ink,annotation,font)
        label.rotation_euler=s.camera.rotation_euler;label.parent=obj
        label.location=points.mean(0)+right*.5+up*.4+front*20
        for state in poses:
            show=state['seconds']>=15.25 and (state['second_joining']<=.05 if kind=='substrate' else state['second_joining']>=.95)
            visibility(label,show,state['frame'])
        new_labels.append(label)
    for name in ['PPi label','Outline PPi label']:
        ob=bpy.data.objects[name]
        for state in poses:
            if state['seconds']>=9.5:visibility(ob,False,state['frame'])
    for original in new_labels:
        d=original.data.copy();d.offset=.010;d.materials.clear();d.materials.append(halo)
        ob=bpy.data.objects.new('Outline '+original.name,d);annotation.objects.link(ob)
        ob.parent=original;ob.location=(0,0,-.02);ob.visible_shadow=False
        for state in poses:
            s.frame_set(state['frame']);visibility(ob,not original.hide_render,state['frame'])
    # Keep camera-attached captions at the same screen margin during the small pullback.
    for ob in list(annotation.objects):
        if ob.parent!=s.camera:continue
        position=np.array(ob.location);scale=np.array(ob.scale)
        for state in poses:
            factor=(camera_scale+1.5*state['advance'])/camera_scale
            ob.location=(position[0]*factor,position[1]*factor,position[2])
            ob.scale=scale*factor
            ob.keyframe_insert(data_path='location',frame=state['frame'])
            ob.keyframe_insert(data_path='scale',frame=state['frame'])
    manifest={**old,'scene':'Two successive RNA additions C then U with slower pacing',
        'source_first_cycle_blend':source_blend,'substrate_reconstruction':info,
        'frame_poses':poses,'first_cycle_retime':{'original_normalized_times':old_times,'new_seconds':seconds},
        'stages':[{'start_seconds':a,'end_seconds':b,'text':body} for a,b,body in stages],
        'next_cycle_geometry':info['next_cycle'],
        'color_semantics':'Gold marks the current added nucleotide; incorporated C20 becomes RNA orange during the inter-cycle pause/advance.',
        'second_atom_bookkeeping':{'UTP_heavy_atoms':len(utp['atoms']),'added_U_residue_heavy_atoms':count,'PPi_heavy_atoms':len(utp['atoms'])-count},
        'sources':old['sources']+['https://www.rcsb.org/ligand/UTP'],
        'limitations':old['limitations']+info['next_cycle']['limitations']}
    # The first-cycle substrate limitation describes that asset, not this new authored motion.
    manifest['limitations']=[x for x in manifest['limitations'] if 'No hydrogen, protonation' not in x]
    manifest['limitations'].append('No hydrogen, protonation, metal chemistry or trigger-loop closure is simulated.')
    bpy.ops.file.pack_all()
    save_and_render(args,manifest)


if __name__=='__main__':main()
