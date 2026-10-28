"""One-nucleotide RNA extension: complete molecules and explicit teaching reconstruction.

5FLM supplies the product scaffold. CCD ideal CTP supplies a chemically connected
substrate. Morphing between them explains bookkeeping, not a catalytic pathway.
"""
import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Vector

sys.path.insert(0, str(Path(__file__).resolve().parent))
from detail_scene import collection, load_chain, material, path, save_and_render, studio
from build_molecular_scene import nucleic_atoms
from build_translation_scene import atom, emission, smooth, text_object, update_line


def subset(data, mask):
    ids = np.flatnonzero(mask)
    mapping = {int(old): new for new, old in enumerate(ids)}
    bonds = [(mapping[int(a)],mapping[int(b)]) for a,b in data['bonds']
             if int(a) in mapping and int(b) in mapping]
    return {'atoms':data['atoms'][ids], 'elements':data['elements'][ids], 'bonds':bonds}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path,required=True)
    p.add_argument('--substrate-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--font',type=Path)
    p.add_argument('--style-profile',type=Path,
                   default=Path(__file__).resolve().parent.parent/'assets/style-molecular-reference.json')
    p.add_argument('--engine',choices=['cycles','eevee'],default='cycles')
    p.add_argument('--compositor',choices=['clean','reference','monochrome'])
    p.add_argument('--width',type=int,default=1920)
    p.add_argument('--samples',type=int,default=64)
    p.add_argument('--frames',type=int,default=168)
    p.add_argument('--fps',type=int,default=24)
    p.add_argument('--frame',type=int,default=140)
    p.add_argument('--render',action='store_true')
    p.add_argument('--animate',action='store_true')
    args = p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    if args.output.exists(): p.error('Choose a fresh output directory.')
    if args.frames<24 or not 1<=args.frame<=args.frames: p.error('Check frame count and preview frame.')
    source = json.loads((args.asset_dir/'manifest.json').read_text())
    substrate_info = json.loads((args.substrate_dir/'manifest.json').read_text())
    if (not source['source_file'].upper().startswith('5FLM') or
        source['source_sha256'] != substrate_info['source_asset']['source_sha256']):
        p.error('Use the matching 5FLM asset and its prepared CTP reconstruction.')
    profile = json.loads(args.style_profile.read_text())
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)
    mats = {k:material(k,v,profile['roughness'],profile['bump']) for k,v in profile['palette'].items()}
    target = np.array((-.1,-2.3,1.2))
    front = np.array((.9396926,.3420201,0)); front /= np.linalg.norm(front)
    right = np.cross((0,0,1),front); right /= np.linalg.norm(right)
    up = np.cross(front,right)
    context = collection('01 Intact Pol II — displaced only for a layered view')
    scaffold = collection('02 Native modeled DNA and RNA')
    reaction = collection('03 One-C addition — explicit teaching reconstruction')
    annotation = collection('04 Captions and polarity')
    pol_mat = material('Pol II contextual surface','276E78',.64,.006)
    shader = pol_mat.node_tree.nodes.get('Principled BSDF')
    color = shader.inputs['Base Color'].default_value
    shader.inputs['Base Color'].default_value = tuple(v*.38 for v in color[:3])+(1,)
    pol_objects, maximum_depth = [], -float('inf')
    for ch in 'ABCDEFGHIJKL':
        ob, d = load_chain(args.asset_dir,ch,pol_mat,context,'Intact Pol II chain '+ch)
        pol_objects.append(ob)
        maximum_depth = max(maximum_depth,float((d['vertices']@front).max()))
    data = {ch:np.load(args.asset_dir/f'chain_{ch}.npz') for ch in 'NPT'}
    template = nucleic_atoms(data['T'],'Template DNA — all modeled residues',mats['dna2'],scaffold)
    template['source_nucleic_chain'] = 'T'
    non_template = nucleic_atoms(data['N'],'Non-template DNA — hidden whole chain',mats['dna1'],scaffold)
    non_template['source_nucleic_chain'] = 'N'
    non_template.hide_render = non_template.hide_viewport = True
    primer_data = subset(data['P'],data['P']['residue_numbers']<20)
    primer = nucleic_atoms(primer_data,'RNA primer P residues 7–19',mats['rna'],scaffold)
    primer['source_nucleic_chain'] = 'P'; primer['residue_range'] = '7–19'
    ctp = np.load(args.substrate_dir/'ctp.npz')
    count = int(ctp['product_atom_count'])
    core_data = subset(ctp,np.arange(len(ctp['atoms']))<count)
    extra_data = subset(ctp,np.arange(len(ctp['atoms']))>=count)
    core = nucleic_atoms(core_data,'Incoming CTP core → RNA C20',mats['protein'],reaction)
    extra = nucleic_atoms(extra_data,'CTP beta/gamma phosphates → PPi',mats['accent'],reaction)
    # Each endpoint has intact local geometry. Only the authored reaction segment
    # interpolates the substrate core to its observed product coordinates.
    product_data = {**core_data,'atoms':ctp['product_atoms']}
    product = nucleic_atoms(product_data,'Temporary product geometry',mats['protein'],reaction)
    core.shape_key_add(name='CCD substrate')
    final_key = core.shape_key_add(name='Observed 5FLM C20 product')
    values = np.empty(len(product.data.vertices)*3,dtype=np.float32)
    product.data.vertices.foreach_get('co',values)
    final_key.data.foreach_set('co',values)
    product_mesh = product.data
    bpy.data.objects.remove(product,do_unlink=True)
    bpy.data.meshes.remove(product_mesh)
    names = list(ctp['atom_names'])
    pa, o3a = ctp['atoms'][names.index('PA')], ctp['atoms'][names.index('O3A')]
    product_p = ctp['product_atoms'][list(ctp['product_names']).index('P')]
    primer_o3 = atom(data['P'],19,"O3'")
    product_o3 = ctp['product_atoms'][list(ctp['product_names']).index("O3'")]
    base_ids = [names.index(n) for n in ['N1','C2','N3','C4','C5','C6','O2','N4']]
    outward = ctp['product_atoms'][base_ids].mean(0)-data['T']['atoms'][data['T']['residue_numbers']==17].mean(0)
    outward /= np.linalg.norm(outward)
    cleavage = path('Breaking alpha–beta phosphate connection',[(0,0,0)]*2,.035,mats['accent'],reaction)
    addition = path('New RNA O3′–P bond',[(0,0,0)]*2,.035,mats['rna'],reaction)
    poses = []
    for frame in range(1,args.frames+1):
        t = (frame-1)/(args.frames-1)
        docking = smooth(t,.12,.40)
        joining = smooth(t,.43,.58)
        release = smooth(t,.58,.84)
        approach = (1-docking)*(right*3+up*.8+front*2)
        hold_offset = (front*.30+outward*.30)*(1-joining)
        core.location = approach+hold_offset
        extra.location = approach+hold_offset+release*(right*3.2-up*.9+front*1.8)
        core.keyframe_insert(data_path='location',frame=frame)
        extra.keyframe_insert(data_path='location',frame=frame)
        final_key.value = joining
        final_key.keyframe_insert(data_path='value',frame=frame)
        current_p = (1-joining)*pa+joining*product_p+np.array(core.location)
        update_line(cleavage,[current_p,o3a+np.array(extra.location)],frame)
        update_line(addition,[primer_o3,current_p],frame)
        cleavage.data.bevel_depth = .035*(1-joining)
        cleavage.data.keyframe_insert(data_path='bevel_depth',frame=frame)
        addition.data.bevel_depth = .035*joining
        addition.data.keyframe_insert(data_path='bevel_depth',frame=frame)
        poses.append({'frame':frame,'docking':docking,'joining':joining,'release':release,
                      'core_translation_nm':list(core.location),'ppi_translation_nm':list(extra.location)})
    actor_points = np.concatenate([data['T']['atoms'],data['P']['atoms'],ctp['atoms']])
    offset = -front*max(0.,maximum_depth-float((actor_points@front).min())+1.8)
    for ob in pol_objects: ob.location = offset
    camera = studio(profile,args,{'position':target+front*50,'target':target,'scale':14.5,'light_scale':.55})
    light = bpy.data.lights.new('Readable active-site fill','AREA')
    light.energy,light.size = 4200,10
    lob = bpy.data.objects.new(light.name,light)
    bpy.data.collections['90 Camera and lighting'].objects.link(lob)
    lob.location = target+front*14+up*6
    lob.rotation_euler = (Vector(target)-lob.location).to_track_quat('-Z','Y').to_euler()
    font = bpy.data.fonts.load(str(args.font)) if args.font else None
    zh = bool(font)
    bg = profile['background'].lstrip('#')
    light_background = sum(int(bg[i:i+2],16) for i in (0,2,4))>480
    ink = emission('Annotation ink','15252C' if light_background else 'E4EAEC')
    stages = [(0,.12,'RNA 在 3′端延长' if zh else 'RNA extends at its 3′ end'),
              (.12,.43,'CTP 接近 RNA 末端' if zh else 'CTP approaches the RNA end'),
              (.43,.58,'连接变化示意：增加一个 C 残基' if zh else 'Illustrative bond change: one C residue added'),
              (.58,.84,'焦磷酸 PPi 离开' if zh else 'Pyrophosphate, PPi, departs'),
              (.84,1.01,'RNA 增加一个残基；末态来自 5FLM' if zh else 'RNA gains one residue; product from 5FLM')]
    for i,(a,b,body) in enumerate(stages):
        ob = text_object(body,f'Stage {i+1}',.25,ink,annotation,font)
        ob.parent = camera; ob.location = (-6.65,3.99,-15)
        for frame in range(1,args.frames+1):
            t=(frame-1)/(args.frames-1)
            ob.hide_render = not a<=t<b
            ob.keyframe_insert(data_path='hide_render',frame=frame)
    footer = ('分层观察：完整 Pol II 后移，非模板链隐藏｜底物接近与成键为教学重构' if zh else
              'Layered view: Pol II offset, non-template strand hidden | Authored approach and bond change')
    ob = text_object(footer,'Presentation scope',.145,ink,annotation,font)
    ob.parent = camera; ob.location = (-6.65,-4.16,-15)
    labels = [
        ('DNA 模板链 3′' if zh else 'Template DNA 3′',atom(data['T'],39,"C4'"),(.9,.8)),
        ('5′',atom(data['T'],1,"C4'"),(.9,.4)),
        ('RNA 片段 5′' if zh else 'RNA fragment 5′',atom(data['P'],7,"C4'"),(.7,.7)),
    ]
    for body,point,(dx,dy) in labels:
        ob = text_object(body,'Polarity '+body,.19,ink,annotation,font)
        ob.rotation_euler = camera.rotation_euler
        ob.data.align_x = 'LEFT' if dx>0 else 'RIGHT'
        ob.location = point+right*dx+up*dy+front*20
        path('Leader '+body,[point+front*20,np.array(ob.location)-right*np.sign(dx)*.1],.007,ink,annotation)
    terminal = text_object('3′','RNA growing 3′ end',.22,ink,annotation)
    terminal.rotation_euler = camera.rotation_euler
    substrate_label = text_object('RNA + CTP → 延长的 RNA + PPi' if zh else 'RNA + CTP → extended RNA + PPi',
                                  'Substrate accounting',.23,ink,annotation,font)
    substrate_label.rotation_euler = camera.rotation_euler
    substrate_label.location = target-right*3-up*3.3+front*20
    for state in poses:
        t=state['joining']
        location=(1-t)*primer_o3+t*product_o3
        terminal.location=location+right*.7-up*.75+front*20
        terminal.keyframe_insert(data_path='location',frame=state['frame'])
    for name,obj,points,phase in [('CTP',core,ctp['atoms'],'substrate'),
                                 ('PPi',extra,ctp['atoms'][count:],'product')]:
        label=text_object(name,name+' label',.20,ink,annotation)
        label.rotation_euler=camera.rotation_euler
        label.parent=obj
        label.location=points.mean(0)+right*.5+up*.4+front*20
        for state in poses:
            label.hide_render=state['joining']>.05 if phase=='substrate' else state['joining']<.95
            label.keyframe_insert(data_path='hide_render',frame=state['frame'])
    halo = emission('Annotation halo','F7F8F9' if light_background else '10212B')
    for original in list(annotation.objects):
        if original.type!='FONT':continue
        d=original.data.copy();d.offset=.010;d.materials.clear();d.materials.append(halo)
        ob=bpy.data.objects.new('Outline '+original.name,d);annotation.objects.link(ob)
        ob.parent=original;ob.location=(0,0,-.02);ob.visible_shadow=False
        if original.name.startswith('Stage '):
            a,b,_=stages[int(original.name.split()[-1])-1]
            for frame in range(1,args.frames+1):
                t=(frame-1)/(args.frames-1);ob.hide_render=not a<=t<b
                ob.keyframe_insert(data_path='hide_render',frame=frame)
        elif original.name in ['CTP label','PPi label']:
            for state in poses:
                ob.hide_render=state['joining']>.05 if original.name=='CTP label' else state['joining']<.95
                ob.keyframe_insert(data_path='hide_render',frame=state['frame'])
    bpy.ops.file.pack_all()
    manifest={'scene':'One-C RNA extension — product-based teaching reconstruction',
        'source_asset':source,'substrate_reconstruction':substrate_info,
        'coordinate_units':'1 Blender unit = 1 nm',
        'presentation':{'molecules_sliced':False,'hidden_complete_chains':['N'],
                        'Pol_II_offset_nm':offset.tolist(),'meaning':'Layered display, not Pol II dissociation.'},
        'frame_poses':poses,'stages':[{'start':a,'end':b,'text':body} for a,b,body in stages],
        'atom_bookkeeping':{'CTP_heavy_atoms':len(ctp['atoms']),'added_C_residue_heavy_atoms':count,
                            'PPi_heavy_atoms':len(ctp['atoms'])-count},
        'new_product_bond_length_nm':float(np.linalg.norm(primer_o3-product_p)),
        'sources':['https://www.rcsb.org/structure/5FLM','https://www.rcsb.org/ligand/CTP',
                   'https://doi.org/10.1038/nature16482'],
        'limitations':substrate_info['limitations']+[
            'The scene explains nucleotide addition, not an entire initiation/elongation/termination cycle.',
            'Product-based interpolation is not an atomistic transition-state or reaction-energy calculation.',
            'The non-template source gap is neither filled nor joined; that entire chain is hidden.',
            'Representations use artistic atom radii and single cylinders for bonds; bond orders remain in CCD data.']}
    save_and_render(args,manifest)


if __name__=='__main__': main()
