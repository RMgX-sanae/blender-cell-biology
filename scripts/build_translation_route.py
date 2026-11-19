"""Join registered mammalian initiation, the native elongation shot, termination and recycling."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Matrix, Vector

sys.path.insert(0,str(Path(__file__).resolve().parent))
import detail_scene as ds
from build_translation_scene import smooth
from molecular_asset_io import chain_file
from scene_chapters import animated_path, balls, caption_set, compose, empty, font_pair, new_shot, visible


INIT,STOP,RECYCLE=192,216,168
ZERO=np.zeros(3)
RECIPE=Path(__file__).resolve().parents[1]/'assets/translation-route-sources.json'


def unit(a):return a/np.linalg.norm(a)


def keys(obj,events):
    for frame,position in events:
        obj.location=position
        obj.keyframe_insert(data_path='location',frame=frame)


def show_family(obj,first,last):
    for node in [obj,*obj.children_recursive]:
        if node.type in ('MESH','CURVE'):visible(node,first,last)


def materials(profile):
    pal=profile['palette']
    colors={'small':pal['beta'],'large':pal['alpha'],'protein':pal['heads'],
            'mRNA':pal['dna2'],'codon':pal['dna1'],'P':pal['rna'],'E':pal['accent'],
            'peptide':pal['peptide'],'Cterm':pal['amino_acid'],'factor1':'DBAF62',
            'factor3':'69AB9D','recycler':'A990CB','cap':'B76CA6','ion':'9EE2DA','water':'7AD9EF'}
    mats={}
    for role,color in colors.items():
        mat=ds.material('Translation route · '+role,color,roughness=.53,bump=.005)
        shader=mat.node_tree.nodes['Principled BSDF']
        shader.inputs['Emission Color'].default_value=mat.diffuse_color
        shader.inputs['Emission Strength'].default_value=.09 if role in ('small','large','protein') else .28
        mats[role]=mat
    return mats


def camera_pose(cam,frame,center,front,right,scale):
    right=unit(right-front*np.dot(right,front));up=np.cross(front,right)
    cam.rotation_euler=Matrix(np.column_stack((right,up,front)).tolist()).to_euler()
    cam.location=center+front*66
    cam.data.ortho_scale=scale
    for prop in ('location','rotation_euler'):cam.keyframe_insert(data_path=prop,frame=frame)
    cam.data.keyframe_insert(data_path='ortho_scale',frame=frame)


def local_atom(data,residue,atom):
    indices=np.flatnonzero((data['residue_numbers']==residue)&(data['atom_names']==atom))
    if len(indices)!=1:raise ValueError(f'Ambiguous or absent landmark {residue} {atom}')
    return data['atoms'][indices[0]].copy()


def molecular_actor(folder,entry,chain,role,group,parent,record,mats):
    name=entry+' · '+chain+' · '+role
    if role=='mRNA':
        with np.load(chain_file(folder,chain),allow_pickle=False) as data:
            start=record['codon'][1];parts=ds.Parts()
            for pos,element,number in zip(data['atoms'],data['elements'],data['residue_numbers']):
                parts.ball(pos,ds.SPACEFILL_RADII_NM.get(str(element),.16),
                           mat=int(start<=number<start+3),seg=7,rings=5)
            obj=parts.object(name,[mats['mRNA'],mats['codon']],group)
            obj['codon_author_start']=start
    else:
        obj,archive=ds.load_chain(folder,chain,mats[role],group,name,
                                  representation='spacefill' if role in ('P','E') else 'surface')
        archive.close()
    obj.parent=parent;obj['pdb']=entry;obj['source_chain']=chain
    for lig in record['companions']:
        if lig['chain']!=chain:continue
        points=[atom['xyz_nm'] for atom in lig['atoms']]
        marker=balls(entry+' · '+lig['component']+' '+str(lig['number']),points,.13,mats['ion'],group,obj)
        marker['source_component']=lig['component']
    return obj


def intact_subunit(reference,prefix,group,parent,mats):
    result=[]
    for native in reference.objects:
        if not native.name.startswith(prefix+' complete chain '):continue
        chain=native['source_chain']
        role=('small' if prefix=='40S' else 'large') if chain in ('S2','L5','L7','L8') else 'protein'
        obj=bpy.data.objects.new(prefix+' · '+chain,native.data)
        group.objects.link(obj);obj.parent=parent
        obj.material_slots[0].link='OBJECT';obj.material_slots[0].material=mats[role]
        obj['pdb']='9SPI';obj['source_chain']=chain
        result.append(obj)
    if not result:raise ValueError('Missing intact donor subunit '+prefix)
    return result


def overview_axes(donor):
    centers={row['chain']:np.array(row['centroid_nm']) for row in donor['source_asset']['chains']}
    separation=unit(centers['L5']-centers['S2'])
    travel=centers['At']-centers['Pt']
    right=unit(travel-separation*np.dot(travel,separation))
    front=unit(np.cross(right,separation)+.30*separation)
    return separation,right,front


def shot(name,n,center,front,right,scale,reference,profile):
    scene,cam=new_shot(name,n,{'position':center+front*66,'target':center,'scale':scale,'light_scale':1.3},reference,profile)
    camera_pose(cam,1,center,front,right,scale)
    return scene,cam,ds.collection('Mechanism · '+name)


def initiation(prepared,reference,record,axes,profile,fonts,lang,mats):
    sep,right,front=axes
    scene,cam,g=shot('Initiation · 48S to 80S',INIT,np.array([-1.,-1.,0.]),front,right,46,reference,profile)
    folder=prepared/'48s-initiation'
    small=empty('Scanning 40S',g);small['entity_id']='40s'
    rna=empty('Coding mRNA at AUG',g);rna['entity_id']='coding-mRNA'
    early=empty('eIF2 and initiation-state ABCE1',g);early.parent=small
    eif1a=empty('Joining factor eIF1A',g);eif1a.parent=small
    role=record['roles']
    for item in json.loads((folder/'manifest.json').read_text())['chains']:
        ch=item['chain'];part=small;material='small' if item['nucleic_acid'] else 'protein'
        if ch==role['mRNA']:part,material=rna,'mRNA'
        elif ch==role['initiator_tRNA']:material='P'
        elif ch==role['eIF1A']:part,material=eif1a,'factor1'
        elif ch in role['eIF2']:part,material=early,'factor1'
        elif ch==role['ABCE1']:part,material=early,'recycler'
        actor=molecular_actor(folder,'6YAL',ch,material,g,part,record,mats)
        if ch=='2':
            # From this camera only the intact 18S chain masks the codon/anticodon.
            for f,hidden in [(1,False),(55,True),(101,False)]:
                actor.hide_render=actor.hide_viewport=hidden
                actor.keyframe_insert(data_path='hide_render',frame=f)
                actor.keyframe_insert(data_path='hide_viewport',frame=f)
    with np.load(chain_file(folder,role['mRNA']),allow_pickle=False) as data:
        start=local_atom(data,36,"C4'");next_site=local_atom(data,37,"C4'")
        aug=np.mean([local_atom(data,i,"C4'") for i in range(52,55)],axis=0)
    direction=unit(start-next_site)
    leader=np.array([start,start+direction*.9,start+direction*1.9-.7*right,start+direction*2.5-2.2*right])
    stub=ds.path('Schematic upstream leader',leader,.11,mats['mRNA'],g);stub.parent=rna
    cap=balls('5′ cap · teaching mark',[leader[-1]],.32,mats['cap'],g,rna)
    cap['scope']='5′ extension and cap are schematic; the modeled mRNA begins at residue 36'
    keys(small,[(1,-5.4*right+.45*front),(24,-4.9*right+.18*front),(54,ZERO)])
    keys(early,[(1,ZERO),(93,ZERO),(116,1.9*front+.9*right),(132,5.3*front+2*right)])
    show_family(early,1,132)
    keys(eif1a,[(1,ZERO),(150,ZERO),(177,3*front-1.5*right)])
    show_family(eif1a,1,177)
    large=empty('Joining 60S',g);large['entity_id']='60s'
    components=intact_subunit(reference,'60S',g,large,mats)
    keys(large,[(102,15*sep+4*right),(116,13*sep+3.5*right),(145,4*sep+.55*right),(164,ZERO)])
    for obj in components:visible(obj,102,INIT)
    close=unit(np.array([.572,.777,.262]))
    views=[(1,np.array([-1.,-1.,0.]),46,front),(54,np.array([-2.,-2.,0.]),39,front),
           (55,aug+.15*sep,15,close),(100,aug+.15*sep,15,close),
           (101,np.array([-1.,1.,0.]),55,front),(INIT,ZERO,49,front)]
    for f,target,size,direction in views:camera_pose(cam,f,target,direction,right,size)
    caption_set(scene,cam,fonts,lang,[
        (1,'40S 携起始 tRNA 从 5′ 端扫描','The 40S scans from the 5′ end with initiator tRNA'),
        (55,'起始 tRNA 识别 AUG','Initiator tRNA recognizes AUG'),
        (101,'eIF2 释放，准备接纳大亚基','eIF2 leaves after start-codon recognition'),
        (133,'60S 加入，eIF5B 水解 GTP 后离开','60S joins, then eIF5B hydrolyzes GTP and leaves')],
        [(f,size) for f,_,size,_ in views],outline=True)
    scene.frame_set(1)
    return scene


def peptide_motion(prepared,group,mats):
    guide=json.loads((prepared/'exit-guide.json').read_text())
    path=np.asarray(guide['world_path_nm']);trace=np.asarray(guide['C_to_N_trace_nm'])
    arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(path,axis=0),axis=1))]
    base=arc[:len(trace)];travel=arc[-2]
    def coords(f):
        position=base+travel*smooth(f,112,199)
        return np.column_stack([np.interp(position,arc,path[:,i]) for i in range(3)])
    frames=sorted({1,107,112,199,STOP,*range(114,200,3)})
    ribbon=animated_path('Nascent peptide · C-terminal trace',frames,coords,.125,mats['peptide'],group)
    ribbon['entity_id']='translated-peptide';ribbon['residues']='5LZT Val68 carbonyl and C-alpha 68–54'
    tip=balls('Released C terminus · Val68',[ZERO],.22,mats['Cterm'],group)
    tip['entity_id']='peptide-C68'
    keys(tip,[(f,coords(f)[0]) for f in frames])
    return coords,trace[0]


def termination(prepared,reference,donor,records,profile,fonts,lang,mats):
    direction=reference.camera.rotation_euler.to_matrix()@Vector((0,0,1))
    front=np.array(direction);right=np.array(reference.camera.rotation_euler.to_matrix()@Vector((1,0,0)))
    scene,cam,g=shot('Termination · stop to released peptide',STOP,np.array([-3.,0.,.4]),front,right,26,reference,profile)
    back=empty('Whole 60S behind exposed reaction site',g)
    back.location=donor['presentation']['large_subunit_offset_nm']
    intact_subunit(reference,'60S',g,back,mats)
    substrate=empty('Stop substrate · maintained identities',g)
    pre=records['5LZT'];post=records['5LZU']
    for ch,role in [('hh','mRNA'),('2','P'),('3','E')]:
        obj=molecular_actor(prepared/'5LZT','5LZT',ch,role,g,substrate,pre,mats)
        obj['entity_id']={'hh':'coding-mRNA','2':'terminal-P-tRNA','3':'terminal-E-tRNA'}[ch]
    recognize=molecular_actor(prepared/'5LZT','5LZT','ii','factor1',g,substrate,pre,mats)
    show_family(recognize,1,63)
    escort=molecular_actor(prepared/'5LZT','5LZT','jj','factor3',g,substrate,pre,mats)
    keys(escort,[(1,ZERO),(39,ZERO),(57,1.6*right+front),(73,5*right+3*front)])
    show_family(escort,1,73)
    active=molecular_actor(prepared/'5LZU','5LZU','ii','factor1',g,substrate,post,mats)
    show_family(active,64,STOP)
    coords,start=peptide_motion(prepared,g,mats)
    with np.load(chain_file(prepared/'5LZT','2'),allow_pickle=False) as data:
        oxygen=local_atom(data,76,"O3'")
    ester=ds.path('Peptidyl–tRNA ester · schematic', [oxygen,start],.065,mats['Cterm'],g)
    visible(ester,1,108)
    water=balls('Water oxygen · hydrolysis teaching mark',[ZERO],.115,mats['water'],g)
    frames=[82,94,107,109,112,*range(115,201,3),STOP]
    keys(water,[(f,(start+front*(.24+1.1*(1-smooth(f,82,107)))+.6*right*(1-smooth(f,82,107)))
                  if f<=107 else coords(f)[0]+.16*front) for f in frames])
    visible(water,82,STOP)
    views=[(1,np.array([-3.,0.,.4]),26),(48,np.array([-3.,0.,.4]),26),
           (87,start+np.array([-1.0,0,0]),15),(116,start+np.array([-1.0,0,0]),15),
           (146,coords(146).mean(axis=0),18),(188,coords(188).mean(axis=0),24),(STOP,coords(STOP).mean(axis=0),24)]
    for f,center,size in views:camera_pose(cam,f,center,front,right,size)
    caption_set(scene,cam,fonts,lang,[
        (1,'后续延伸后，UGA 进入 A 位','After further elongation, UGA reaches the A site'),
        (30,'eRF1 识别 UGA，eRF3 水解 GTP 后离开','eRF1 reads UGA and eRF3 leaves after GTP hydrolysis'),
        (77,'eRF1 促进肽酰 tRNA 水解','eRF1 promotes peptidyl-tRNA hydrolysis'),
        (112,'肽链与 tRNA 分开，C 端沿出口离开','The freed C terminus follows the peptide exit tunnel')],
        [(f,size) for f,_,size in views],outline=True)
    scene.frame_set(1)
    return scene,coords(STOP)


def recycling(prepared,reference,records,axes,product,profile,fonts,lang,mats):
    sep,right,front=axes
    center=np.array([1.,.5,0.]);size=55
    scene,cam,g=shot('Recycling · separate and reuse',RECYCLE,center,front,right,size,reference,profile)
    small=empty('Recycled 40S',g);large=empty('Released 60S',g)
    small['entity_id']='40s';large['entity_id']='60s'
    intact_subunit(reference,'40S',g,small,mats);intact_subunit(reference,'60S',g,large,mats)
    keys(large,[(1,ZERO),(41,ZERO),(69,3.8*sep+1.5*right),(97,7*sep+3*right),(RECYCLE,7.6*sep+3.3*right)])
    keys(small,[(1,ZERO),(41,ZERO),(93,-2.5*sep-1.6*right),(RECYCLE,-2.8*sep-1.6*right)])
    for chain,role,ident,out in [('hh','mRNA','coding-mRNA',-11*right+9*front+1.5*sep),
                                 ('2','P','terminal-P-tRNA',7.5*right+10*front+1.3*sep)]:
        actor=molecular_actor(prepared/'5LZT','5LZT',chain,role,g,small,records['5LZT'],mats)
        actor['entity_id']=ident
        keys(actor,[(1,ZERO),(94,ZERO),(116,.4*out+front),(145,out),(RECYCLE,out+front)])
    factor=molecular_actor(prepared/'5LZU','5LZU','ii','factor1',g,small,records['5LZU'],mats)
    keys(factor,[(1,ZERO),(54,ZERO),(80,3*front+3*right),(95,5*front+5*right)])
    show_family(factor,1,95)
    ab=molecular_actor(prepared/'5LZV','5LZV','jj','recycler',g,small,records['5LZV'],mats)
    keys(ab,[(1,5*right+2*front),(22,1.4*right+front),(37,ZERO),(107,ZERO),(142,5*right+3*front)])
    show_family(ab,1,142)
    chain=ds.path('Released peptide · same C-terminal fragment',product,.125,mats['peptide'],g)
    chain['entity_id']='translated-peptide'
    end=balls('Released peptide C68',[product[0]],.22,mats['Cterm'],g,chain)
    keys(chain,[(1,ZERO),(RECYCLE,2*right+1.3*front)])
    views=[(1,center,size),(42,center,size),(94,center+1.2*sep,68),(RECYCLE,center+1.2*sep,68)]
    for f,target,z in views:camera_pose(cam,f,target,front,right,z)
    caption_set(scene,cam,fonts,lang,[
        (1,'ABCE1 结合终止后的核糖体','ABCE1 engages the post-termination ribosome'),
        (42,'ABCE1 利用 ATP 分开大小亚基','ABCE1 uses ATP to split the subunits'),
        (95,'去酰化 tRNA 与 mRNA 从 40S 释放','Deacylated tRNA and mRNA leave the 40S'),
        (144,'大小亚基可进入下一轮翻译','Both subunits are available for another round')],
        [(f,z) for f,_,z in views],outline=True)
    scene.frame_set(1)
    return scene


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--font-zh',type=Path);parser.add_argument('--font-latin',type=Path)
    args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
    if args.output.exists():parser.error('Choose a fresh output .blend')
    native=Path(bpy.data.filepath);donor=json.loads(native.with_name('manifest.json').read_text())
    reference=bpy.context.scene
    if (donor.get('scene')!='Human translation elongation substrate cycle' or len(donor['source_asset']['chains'])!=83
            or (reference.frame_start,reference.frame_end,reference.render.fps,reference.render.fps_base)!=(1,192,24,1)):
        raise ValueError('Open the source-qualified 9SPI 192-frame 24 fps elongation scene')
    records=json.loads((args.prepared/'registered/registration.json').read_text())
    recipe=json.loads(RECIPE.read_text())
    if donor['source_asset']['source_sha256']!=recipe['reference']['sha256']:
        raise ValueError('9SPI source changed')
    for spec in recipe['states']:
        pdb=spec['pdb']; folder=args.prepared/('48s-initiation' if pdb=='6YAL' else pdb)
        source=json.loads((folder/'manifest.json').read_text())
        if source['source_sha256']!=spec['sha256'] or set(records[pdb]['polymer_chains'])!={c['chain'] for c in source['chains']}:
            raise ValueError('Source preparation differs: '+pdb)
    appearance=donor.get('appearance',{});profile=appearance.get('profile') or json.loads((RECIPE.parent/'style-translation-reference.json').read_text())
    lang='zh' if any('\u4e00'<=c<='\u9fff' for c in reference.objects['Stage 1'].data.body) else 'en'
    fonts=font_pair(args,lang);mats=materials(profile);axes=overview_axes(donor)
    start=initiation(args.prepared,reference,records['6YAL'],axes,profile,fonts,lang,mats)
    stop,product=termination(args.prepared,reference,donor,records,profile,fonts,lang,mats)
    recycle=recycling(args.prepared,reference,records,axes,product,profile,fonts,lang,mats)
    for obj in reference.objects:
        if obj.name.startswith('Peptide residue '):obj['entity_id']='translated-peptide'
    for stage,cn,en in [(1,'eEF1A–GTP 递送氨酰 tRNA','eEF1A–GTP delivers aminoacyl-tRNA'),
                        (4,'eEF2–GTP 驱动一个密码子的转位','eEF2–GTP drives one-codon translocation')]:
        obj=reference.objects['Stage '+str(stage)];body=cn if lang=='zh' else en
        obj.data.body=body
        for c,fmt in zip(body,obj.data.body_format):fmt.use_bold=(lang=='zh' and c.isascii())
        halo=reference.objects.get('Annotation outline '+obj.name)
        if halo:halo.data.body=body
    titles=('扫描与起始','肽转移和转位','终止和肽链释放','亚基分离与回收') if lang=='zh' else ('Start recognition and joining','Elongation','Termination and release','Recycling')
    master=compose('Translation · full teaching route',[(titles[0],start,1,INIT),(titles[1],reference,1,192),
                   (titles[2],stop,1,STOP),(titles[3],recycle,1,RECYCLE)],reference)
    master['translation_route']=json.dumps({'source_blend_sha256':hashlib.sha256(native.read_bytes()).hexdigest(),
        'sources':[s['pdb'] for s in recipe['states']]+['9SPI'],'language':lang,
        'scope':'A mammalian teaching synthesis across registered initiation, elongation and stop specimens. Functional handoffs and peptide exit are schematic.'})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    bpy.ops.wm.save_as_mainfile(filepath=str(args.output.resolve()),compress=True)
    print('TRANSLATION_ROUTE '+master['native_chapters'])


if __name__=='__main__':main()
