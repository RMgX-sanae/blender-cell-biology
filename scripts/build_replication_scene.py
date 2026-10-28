"""A detailed teaching fork: intact nucleotide units, explicit strand polarity.

1BNA supplies local atomic shapes; 7PFO supplies an intact CMG reference layer.
The fork, sequence, primer symbols and movement are authored teaching geometry.
"""
import argparse
import json
import math
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Matrix,Vector

sys.path.insert(0,str(Path(__file__).resolve().parent))
from detail_scene import Parts,collection,duplicate,load_chain,material,path,save_and_render,studio
from build_molecular_scene import nucleic_atoms
from build_translation_scene import emission,smooth,text_object,update_line


def pose(index,side,fork,count):
    u=min(1.,max(0.,(fork-index)/14))
    y=side*3.8*u*u*(3-2*u)
    slope=-side*3.8*6*u*(1-u)/(14*.34)
    a=math.atan(slope);h=index*math.tau/10.5
    rz=np.array([[math.cos(a),-math.sin(a),0],[math.sin(a),math.cos(a),0],[0,0,1]])
    rx=np.array([[1,0,0],[0,math.cos(h),-math.sin(h)],[0,math.sin(h),math.cos(h)]])
    return rz@rx,np.array(((index-(count-1)/2)*.34,y,0))


def set_pose(ob,rot,loc,frame,visible=True):
    ob.rotation_mode='QUATERNION';ob.rotation_quaternion=Matrix(rot).to_quaternion();ob.location=loc
    ob.keyframe_insert(data_path='rotation_quaternion',frame=frame)
    ob.keyframe_insert(data_path='location',frame=frame)
    ob.hide_render=ob.hide_viewport=not visible
    ob.keyframe_insert(data_path='hide_render',frame=frame)
    ob.keyframe_insert(data_path='hide_viewport',frame=frame)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path,required=True)
    p.add_argument('--replisome-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--font',type=Path)
    p.add_argument('--style-profile',type=Path,default=Path(__file__).resolve().parent.parent/'assets/style-molecular-reference.json')
    p.add_argument('--engine',choices=['cycles','eevee'],default='cycles')
    p.add_argument('--compositor',choices=['clean','reference','monochrome'])
    p.add_argument('--width',type=int,default=1920)
    p.add_argument('--samples',type=int,default=64)
    p.add_argument('--frames',type=int,default=144)
    p.add_argument('--fps',type=int,default=24)
    p.add_argument('--frame',type=int,default=100)
    p.add_argument('--render',action='store_true')
    p.add_argument('--animate',action='store_true')
    a=p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    if a.output.exists():p.error('Choose a fresh output directory.')
    if a.frames<24 or not 1<=a.frame<=a.frames:p.error('Check frame range.')
    source=json.loads((a.asset_dir/'manifest.json').read_text())
    cmg_source=json.loads((a.replisome_dir/'manifest.json').read_text())
    if not cmg_source['source_file'].upper().startswith('7PFO'):p.error('CMG chain selection requires 7PFO.')
    profile=json.loads(a.style_profile.read_text())
    bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
    parents=collection('01 Parent strands — identities retained')
    new=collection('02 Newly synthesized DNA and primer symbols')
    links=collection('03 Backbone guides and pairing symbols')
    context=collection('04 Intact CMG — offset reference layer')
    proto_group=collection('80 Reusable complete nucleotide prototypes')
    annotation=collection('05 Polarity and explanation')
    palette=profile['palette']
    mats={role:material(role,palette[key],profile['roughness'],min(profile['bump'],.003))
          for role,key in {'parent':'dna2','leading':'protein','lagging':'dna1','primer':'accent'}.items()}
    count=72;sequence=('CGAT'*18)
    data={(b,ch):np.load(a.asset_dir/f'{b}_{ch}.npz') for b in 'CGAT' for ch in 'AB'}
    prototypes={}
    for key,d in data.items():
        ob=nucleic_atoms(d,'Prototype '+''.join(key),mats['parent'],proto_group)
        ob.hide_render=ob.hide_viewport=True;prototypes[key]=ob
    prototype_primer={}
    for b in 'CGAT':
        d=data[b,'B'];pt=d['atoms'][list(d['atom_names']).index("C4'")]
        part=Parts();part.ball(pt,.18,seg=14,rings=9)
        ob=part.object('Primer symbol prototype '+b,[mats['primer']],proto_group)
        ob.hide_render=ob.hide_viewport=True;prototype_primer[b]=ob
    objects={};backbones={}
    for role,ch,mat in [('parent_A','A','parent'),('parent_B','B','parent'),
                        ('leading','A','leading'),('lagging','B','lagging')]:
        for i,b in enumerate(sequence):
            primer=role=='lagging' and i%12>=9
            ob=duplicate(prototype_primer[b] if primer else prototypes[b,ch],f'{role} nucleotide {i:02}',parents if role.startswith('parent') else new)
            ob.hide_render=ob.hide_viewport=False
            ob.material_slots[0].link='OBJECT';ob.material_slots[0].material=mats['primer' if primer else mat]
            ob['strand_role']=role;ob['index']=i;ob['prototype']=b+'_'+ch
            ob['representation']='RNA primer symbol' if primer else 'complete 1BNA nucleotide prototype'
            objects[role,i]=ob
        for i in range(count-1):
            if role=='lagging' and i//12!=(i+1)//12:continue
            backbones[role,i]=path(f'{role} backbone guide {i:02}',[(0,0,0)]*2,.038,mats[mat],links)
    cmg=[];cmg_points=[]
    cmg_mat=material('Intact CMG depth layer',palette['heads'],.62,profile['bump'])
    shader=cmg_mat.node_tree.nodes.get('Principled BSDF');rgba=shader.inputs['Base Color'].default_value
    shader.inputs['Base Color'].default_value=tuple(v*.65 for v in rgba[:3])+(1,)
    for ch in '234567CDEFG':
        ob,d=load_chain(a.replisome_dir,ch,cmg_mat,context,'CMG complete chain '+ch)
        cmg.append(ob);cmg_points.extend(d['vertices'])
    cmg_points=np.array(cmg_points);cmg_center=cmg_points.mean(0)
    front=np.array((0,-.35,1.));front/=np.linalg.norm(front)
    right=np.array((1.,0,0));up=np.cross(front,right)
    depth=float(((cmg_points-cmg_center)@front).max())+5
    camera=studio(profile,a,{'position':front*65,'target':(0,0,0),'scale':33,'light_scale':1.1})
    font=bpy.data.fonts.load(str(a.font)) if a.font else None;zh=bool(font)
    bg=profile['background'].lstrip('#');light=sum(int(bg[i:i+2],16) for i in (0,2,4))>480
    ink=emission('Annotation ink','15252C' if light else 'E4EAEC')
    def label(body,name,size,x,y,mat=ink):
        ob=text_object(body,name,size,mat,annotation,font)
        ob.rotation_euler=camera.rotation_euler
        ob.location=right*x+up*y+front*20
        return ob
    label('后随链：分段延伸' if zh else 'Lagging strand: discontinuous','Lagging explanation',.44,-13.7,7.3)
    label('前导链：连续延伸' if zh else 'Leading strand: continuous','Leading explanation',.44,-13.7,-6.7)
    label('RNA 引物符号' if zh else 'RNA primer symbols','Primer legend',.32,2.5,7.3)
    label('亲代 DNA' if zh else 'Parental DNA','Parent legend',.34,5,9)
    path('Parent color key',[right*3.8+up*9.2+front*20,right*4.7+up*9.2+front*20],.055,mats['parent'],annotation)
    path('Primer color key',[right*1.3+up*7.5+front*20,right*2.2+up*7.5+front*20],.055,mats['primer'],annotation)
    label('分层示意：CMG 后移；复制叉与运动为构造｜未显示引物去除和连接' if zh else
          'Layered CMG reference; authored fork and motion | Primer removal and ligation omitted',
          'Presentation scope',.24,-14.8,-9.3)
    for body,x,y in [('5′',-13.4,4.6),('3′',-13.4,-4.9),('3′',12.6,.7),('5′',12.6,-1.4)]:
        label(body,f'Parent polarity {body} {x} {y}',.32,x,y)
    label('5′','Leading new 5′',.31,-12.1,-5.5)
    lead_end=label('3′','Leading new 3′',.31,0,-5.5)
    lead_arrow=path('Leading synthesis direction',[(0,0,0)]*2,.024,mats['leading'],annotation)
    lead_tip=path('Leading direction arrowhead',[(0,0,0)]*3,.024,mats['leading'],annotation)
    lag_arrows={};lag_labels={}
    for j in [11,23,35,47,59,71]:
        lag_arrows[j]=(path(f'Lagging fragment {j//12+1} direction',[(0,0,0)]*2,.024,mats['lagging'],annotation),
                       path(f'Lagging fragment {j//12+1} arrowhead',[(0,0,0)]*3,.024,mats['lagging'],annotation))
        lag_labels[j]=(label('5′',f'Lagging fragment {j//12+1} 5′',.26,0,5.4),
                       label('3′',f'Lagging fragment {j//12+1} 3′',.26,0,5.4))
    # Each daughter duplex uses the same matched local nucleotide-pair geometry.
    frames=[];guide_lengths=[]
    def point(role,i,name,poses):
        ch='A' if role in ('parent_A','leading') else 'B'
        d=data[sequence[i],ch]
        if role=='lagging' and i%12>=9:name="C4'"
        xyz=d['atoms'][list(d['atom_names']).index(name)]
        rot,loc=poses[role,i]
        return rot@xyz+loc
    for frame in range(1,a.frames+1):
        t=(frame-1)/(a.frames-1);fork=16+35*smooth(t,.06,.88)
        poses={};visible={};leading_end=math.floor(fork-5)
        fragments=[]
        for j in [11,23,35,47,59,71]:
            active=fork>=j+5
            dna=min(9,max(0,math.floor((fork-j-5)*1.4))) if active else 0
            fragments.append({'right_index':j,'active':active,'DNA_indices':list(range(j-2-dna,j-2)) if dna else [],
                              'primer_indices':[j-2,j-1,j] if active else []})
        lag_present={i for f in fragments for i in f['DNA_indices']+f['primer_indices']}
        for role in ['parent_A','parent_B','leading','lagging']:
            for i in range(count):
                side=1 if role in ['parent_A','lagging'] else -1
                rot,loc=pose(i,side,fork,count);poses[role,i]=(rot,loc)
                show=role.startswith('parent') or (i<=leading_end if role=='leading' else i in lag_present)
                visible[role,i]=show;set_pose(objects[role,i],rot,loc,frame,show)
            for i in range(count-1):
                if (role,i) not in backbones:continue
                ob=backbones[role,i]
                endpoints=([point(role,i,"O3'",poses),point(role,i+1,'P',poses)] if role in ['parent_A','leading'] else
                           [point(role,i+1,"O3'",poses),point(role,i,'P',poses)])
                update_line(ob,endpoints,frame)
                show=visible[role,i] and visible[role,i+1]
                ob.hide_render=ob.hide_viewport=not show
                ob.keyframe_insert(data_path='hide_render',frame=frame)
                ob.keyframe_insert(data_path='hide_viewport',frame=frame)
                if show:guide_lengths.append(float(np.linalg.norm(endpoints[1]-endpoints[0])))
        fork_x=(fork-(count-1)/2)*.34
        cmg_offset=np.array((fork_x,0,0))-cmg_center-front*depth
        for ob in cmg:
            ob.location=cmg_offset;ob.keyframe_insert(data_path='location',frame=frame)
        end_x=(leading_end-(count-1)/2)*.34
        aa=right*(-12.1)+up*(-5.9)+front*20
        bb=right*end_x+up*(-5.9)+front*20
        update_line(lead_arrow,[aa,bb],frame)
        update_line(lead_tip,[bb-right*.3+up*.18,bb,bb-right*.3-up*.18],frame)
        lead_end.location=right*(end_x-.1)+up*(-5.5)+front*20
        lead_end.keyframe_insert(data_path='location',frame=frame)
        for f in fragments:
            j=f['right_index'];start_x=(j-(count-1)/2)*.34
            end_x=(j-2-len(f['DNA_indices'])-(count-1)/2)*.34
            arrow_y=5.1+(j//12%2)*.75
            start=right*start_x+up*arrow_y+front*20;end=right*end_x+up*arrow_y+front*20
            arrow,tip=lag_arrows[j]
            update_line(arrow,[start,end],frame)
            update_line(tip,[end+right*.25+up*.16,end,end+right*.25-up*.16],frame)
            for ob in [arrow,tip,*lag_labels[j]]:
                ob.hide_render=ob.hide_viewport=not f['active']
                ob.keyframe_insert(data_path='hide_render',frame=frame)
                ob.keyframe_insert(data_path='hide_viewport',frame=frame)
            for ob,x in zip(lag_labels[j],[start_x-.1,end_x-.1]):
                ob.location=right*x+up*(arrow_y+.3)+front*20;ob.keyframe_insert(data_path='location',frame=frame)
        frames.append({'frame':frame,'fork_index':fork,'fork_x_nm':fork_x,'leading_DNA_indices':list(range(leading_end+1)),
                       'lagging_fragments':fragments,'CMG_offset_nm':cmg_offset.tolist()})
    label('两条新链均以 5′ → 3′ 合成' if zh else 'Both new strands grow 5′ → 3′','Synthesis rule',.48,-14.8,9)
    halo=emission('Annotation halo','F8F9FA' if light else '10212B')
    for original in list(annotation.objects):
        if original.type!='FONT':continue
        d=original.data.copy();d.offset=.016;d.materials.clear();d.materials.append(halo)
        ob=bpy.data.objects.new('Outline '+original.name,d);annotation.objects.link(ob)
        ob.parent=original;ob.location=(0,0,-.02);ob.visible_shadow=False
        if original.name.startswith('Lagging fragment '):
            number=int(original.name.split()[2])-1
            for state in frames:
                ob.hide_render=ob.hide_viewport=not state['lagging_fragments'][number]['active']
                ob.keyframe_insert(data_path='hide_render',frame=state['frame'])
                ob.keyframe_insert(data_path='hide_viewport',frame=state['frame'])
    bpy.ops.file.pack_all()
    save_and_render(a,{'scene':'Semiconservative replication fork — explicit teaching construction',
        'source_prototypes':source,'CMG_source':cmg_source,'sequence_parent_A_5_to_3':sequence,
        'complement_parent_B_3_to_5':sequence.translate(str.maketrans('ATCG','TAGC')),
        'units':'1 Blender unit = 1 nm','nucleotide_count_per_parent':count,
        'helix':{'rise_nm':.34,'bp_per_turn':10.5,'interpretation':'canonical display parameters'},
        'presentation':{'molecules_sliced':False,'CMG_chains':'234567CDEFG',
                        'meaning':'Complete CMG is offset behind the authored fork. Its placement is not an experimentally reconstructed DNA–CMG interface.'},
        'frame_poses':frames,'backbone_guide_length_range_nm':[min(guide_lengths),max(guide_lengths)],
        'sources':['https://www.rcsb.org/structure/1BNA','https://www.rcsb.org/structure/7PFO',
                   'https://doi.org/10.15252/embj.2021108819'],
        'limitations':['A local elongation window, not a complete replication cycle or experimental DNA trajectory.',
            'RNA primers are three-unit green symbols; their atomic ribose/base chemistry is not represented.',
            'Leading-strand initiation is outside this window. Polymerases/primase, nucleotide binding, ATP chemistry, proofreading, primer removal and ligation are omitted.',
            'Fragment lengths, fork speed and helix curvature are authored for readability, not measured biological parameters.',
            'Local nucleotide atoms remain intact; connecting tubes bridge the authored path and are not validated covalent bond geometry.',
            'All parent identities persist; new DNA appears as complete units. This does not show diffusion or binding of every nucleotide.']})


if __name__=='__main__':main()
