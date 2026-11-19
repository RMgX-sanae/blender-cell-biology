"""Create the human promoter and mammalian RNA-fate shots around an existing 5FLM two-addition scene."""
import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Matrix, Vector

sys.path.insert(0,str(Path(__file__).resolve().parent))
import detail_scene as ds
from build_translation_scene import emission, smooth, text_object
from scene_chapters import compose

CONFIG=Path(__file__).resolve().parents[1]/'assets/pol2-transcription-cycle.json'
PIC_END=144
ESCAPE_END=72
FATE_END=336


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def visible(obj, start, end):
    for frame,flag in ((start-1,True),(start,False),(end+1,True)):
        obj.hide_viewport=obj.hide_render=flag
        obj.keyframe_insert(data_path='hide_viewport',frame=frame)
        obj.keyframe_insert(data_path='hide_render',frame=frame)


def empty(name,group):
    obj=bpy.data.objects.new(name,None)
    obj.empty_display_size=.3
    group.objects.link(obj)
    return obj


def balls(name,positions,radius,mat,group,parent=None):
    parts=ds.Parts()
    for point in positions:
        parts.ball(point,radius,seg=12,rings=8)
    obj=parts.object(name,[mat],group); obj.parent=parent
    return obj


def times(*events):
    return sorted(set(events))


def new_shot(name,frames,view,reference,profile):
    from types import SimpleNamespace
    shot=bpy.data.scenes.new(name)
    bpy.context.window.scene=shot
    args=SimpleNamespace(width=reference.render.resolution_x,frames=frames,fps=reference.render.fps,
                         samples=reference.cycles.samples if reference.render.engine=='CYCLES' else reference.eevee.taa_render_samples,
                         engine='cycles' if reference.render.engine=='CYCLES' else 'eevee',compositor='clean')
    studio=dict(profile,background=profile.get('background','091A24'),camera={},
                compositor=profile.get('compositor',{'preset':'clean','contrast':2}))
    cam=ds.studio(studio,args,view)
    shot.render.resolution_y=reference.render.resolution_y
    shot.render.fps_base=reference.render.fps_base
    return shot,cam


def caption_set(scene,camera,fonts,lang,stages,scales):
    group=ds.collection('Process labels · '+scene.name)
    overlay=empty('Camera typography · '+scene.name,group)
    overlay.parent=camera; overlay.location=(0,0,-1)
    for frame,scale in scales:
        overlay.scale=(scale,)*3;overlay.keyframe_insert(data_path='scale',frame=frame)
    ink=emission('Caption ink · '+scene.name,'EDF3F0')
    zh,latin=fonts
    for index,(start,cn,en) in enumerate(stages):
        body=cn if lang=='zh' else en
        text=text_object(body,'Phase · '+body,.0275,ink,group,zh if lang=='zh' else latin)
        text.data.font_bold=latin
        for char,formatting in zip(body,text.data.body_format):
            formatting.use_bold=(lang=='zh' and char.isascii())
        text.parent=overlay; text.location=(-.452,.255,0)
        for name in ('visible_shadow','visible_diffuse','visible_glossy','visible_transmission','visible_volume_scatter'):
            setattr(text,name,False)
        end=stages[index+1][0]-1 if index+1<len(stages) else scene.frame_end
        visible(text,start,end)
        scene.timeline_markers.new(body,frame=start)


def render_mats():
    colors={'polA':'42838B','polB':'357484','polSmall':'517D82','tbp':'E8B561','gtf':'819C8C',
            'xpb':'E4A75B','tfiih':'548B85','other':'647A85','X':'AB5893','Y':'E6D3B6',
            'rna':'F58D67','polyA':'FFBB84','ion':'ACCCA1','cpsf':'69BEAD','activity':'A3E0D1',
            'phosphate':'F9D773','ctd':'92ACBB','downstream':'F3C55F'}
    strengths={'rna':.44,'polyA':.42,'X':.23,'Y':.23,'activity':.35,'phosphate':.32,
               'ctd':.20,'cpsf':.15,'xpb':.13,'tbp':.13}
    mats={}
    for role,color in colors.items():
        mat=ds.material('Cycle · '+role,color,roughness=.47,bump=.006)
        shader=mat.node_tree.nodes['Principled BSDF']
        shader.inputs['Emission Color'].default_value=mat.diffuse_color
        shader.inputs['Emission Strength'].default_value=strengths.get(role,.075)
        mats[role]=mat
    return mats


def protein_role(ch):
    if ch=='A':return 'polA'
    if ch=='B':return 'polB'
    if ch in 'CDEFGHIJKL':return 'polSmall'
    if ch=='P':return 'tbp'
    if ch=='V':return 'xpb'
    if ch in 'W0123':return 'tfiih'
    return 'other' if ch=='U' else 'gtf'


def promoter(geo,meta,reference,profile,fonts,lang,mats):
    shot,cam=new_shot('Promoter · human PIC',PIC_END,
        {'position':(2.2,-52,3.7),'target':(2.2,0,3.7),'scale':37,'light_scale':1.35},reference,profile)
    for state,source in [('closed','5IY6'),('initial','5IY8')]:
        g=ds.collection('Promoter · '+source)
        for ch in meta['states'][state]['chains']:
            key=state+'_'+ch
            obj=ds.mesh(source+' · '+ch,geo[key+'_v'].tolist(),geo[key+'_f'].tolist(),[mats[protein_role(ch)]],g)
            obj['source_chain']=ch;obj['pdb']=source
            if ch in 'ABCDEFGHIJKL':obj['entity_id']='pol2'
            if state=='closed':
                begin,end=(1 if ch in 'MNOP' else 29),68
            else:
                begin,end=69,PIC_END
                if ch in meta['states'][state]['occluders_at_rna']:
                    obj.hide_render=obj.hide_viewport=True
                    obj['display_scope']='Whole chain hidden for the transcription-bubble view'
                    begin=end=0
            if begin:visible(obj,begin,end)
            if key+'_ions' in geo:
                ion=balls(source+' · ions '+ch,geo[key+'_ions'],.13,mats['ion'],g)
                if begin: visible(ion,begin,end)
                else: ion.hide_render=ion.hide_viewport=True
        for ch in 'XY':
            key=state+'_'+ch
            centers,bases,ids=(geo[key+'_'+suffix] for suffix in ('backbone','base','ids'))
            model=ds.Parts()
            runs=np.split(np.arange(len(ids)),np.flatnonzero((np.diff(ids)!=1)|(np.linalg.norm(np.diff(centers,axis=0),axis=1)>1.15))+1)
            for run in runs:
                if len(run)>1:model.tube(centers[run],.145,sides=10)
            for center,base in zip(centers,bases):
                model.tube([center,base],.075,sides=8);model.ball(base,.16,seg=10,rings=6)
            strand=model.object(source+' · promoter DNA '+ch,[mats[ch]],g)
            strand['source_chain']=ch;strand['pdb']=source
            visible(strand,1 if state=='closed' else 69,68 if state=='closed' else PIC_END)
    guide=geo['short_rna_guide']; tangent=np.gradient(guide,axis=0)/3
    u=np.linspace(0,1,8,endpoint=False)
    pieces=[cubic(a,a+tangent[i],b-tangent[i+1],b,u) for i,(a,b) in enumerate(zip(guide,guide[1:]))]
    rna=ds.path('Initiating RNA · backbone guide',np.vstack([*pieces,guide[-1]]),.125,mats['rna'],g)
    rna['entity_id']='upstream-rna'
    rna['coordinate_evidence']='Six 5FLM RNA C4′ positions fitted on RPB1; 5IY8 has DNA but no deposited RNA polymer chain'
    visible(rna,95,PIC_END)
    for f,v in [(95,.04),(120,1.)]:
        rna.data.bevel_factor_end=v;rna.data.keyframe_insert(data_path='bevel_factor_end',frame=f)
    for i,(point,base) in enumerate(zip(guide,geo['short_rna_bases'])):
        part=ds.Parts()
        part.tube([point,base],.067,sides=10);part.ball(base,.135,seg=12,rings=8)
        landmark=part.object('Initiating RNA · base '+str(i+1),[mats['rna']],g)
        visible(landmark,95+i*5,PIC_END)
    views=[(1,37,2.2,3.7),(48,37,2.2,3.7),(73,22,.7,1.4),(99,14.8,0,.35),(PIC_END,14.8,0,.35)]
    for f,scale,x,z in views:
        cam.location=(x,-52,z)
        cam.keyframe_insert(data_path='location',frame=f)
        cam.data.ortho_scale=scale;cam.data.keyframe_insert(data_path='ortho_scale',frame=f)
    caption_set(shot,cam,fonts,lang,[(1,'起始因子识别启动子','Initiation factors recognize the promoter'),
        (29,'Pol II 加入起始复合物','Pol II joins the initiation complex'),
        (69,'TFIIH 利用 ATP 打开启动子','TFIIH uses ATP to open the promoter'),
        (95,'短 RNA 在转录泡内形成','A short RNA forms in the transcription bubble')],[(f,v) for f,v,_,_ in views])
    shot.frame_set(1)
    return shot


def cubic(a,b,c,d,u):
    u=np.asarray(u).reshape(-1,1)
    return (1-u)**3*a+3*(1-u)**2*u*b+3*(1-u)*u**2*c+u**3*d


def strands(frame):
    f=frame
    x=-5.45+2.65*smooth(f,1,72)+5.65*smooth(f,73,170)+1.35*smooth(f,170,232)
    site=np.array([x,0.,.26*math.sin(x*.38)-.35])
    departure=smooth(f,232,267)
    pol=site+np.array([.8*departure,.22*departure,2.1*departure])
    feeding=smooth(f,1,133)
    release=smooth(f,141,180)
    fed=smooth(f,86,134)
    # After cleavage the upstream transcript separates from the downstream
    # substrate. Only the latter remains tethered to the active polymerase.
    cut=np.array([.06,-1.85,-2.08])
    join=(1-fed)*(site+np.array([-.25,-.7,-.3]))+fed*cut
    end=join+np.array([-1.85*release,-.17*release,-.78*release])
    start=end+np.array([-.65-4.5*feeding,-.25,.18+.65*feeding])
    phases=np.linspace(0,1,36)
    body=cubic(start,start+[-.2,-.2,1.1],end+[-1.35,-.3,-1.10],end,phases)
    end5=join+release*np.array([.22,-.05,.14])
    eat=smooth(f,183,234)
    follow=lambda u:cubic(end5,end5+[1.,0,.65],site+[-1.45,-1.3,-1.65],site,u)
    fragment=follow(eat+(1-eat)*phases)
    tail_len=smooth(f,154,209)
    tail=cubic(end,end+[.3,-.04,-.28],end+[1.24,-.18,-.75],end+[2.,-.14,-.43],phases*tail_len)
    return site,pol,body,fragment,tail,follow([eat])[0]


def animated_path(name,frames,coordinates,radius,mat,group):
    obj=ds.path(name,coordinates(frames[0]),radius,mat,group)
    obj.shape_key_add(name='First pose')
    keys=obj.data.shape_keys;keys.use_relative=False
    for index,frame in enumerate(frames):
        key=keys.key_blocks[0] if index==0 else obj.shape_key_add(name='Frame '+str(frame))
        for p,co in zip(key.data,coordinates(frame)):p.co=co
        key.interpolation='KEY_LINEAR'
        keys.eval_time=key.frame;keys.keyframe_insert(data_path='eval_time',frame=frame)
    action=keys.animation_data.action
    for layer in action.layers:
        for strip in layer.strips:
            for bag in strip.channelbags:
                for curve in bag.fcurves:
                    for point in curve.keyframe_points:point.interpolation='LINEAR'
    return obj


def fates(geo,meta,reference,src,profile,fonts,lang,mats):
    shot,cam=new_shot('Transcript and polymerase · overview',FATE_END,
        {'position':(0,-36,.5),'target':(0,0,.5),'scale':21.5,'light_scale':.75},reference,profile)
    g=ds.collection('One polymerase · transcript fates')
    pol=empty('Pol II · persistent protein complex',g);pol['entity_id']='pol2'
    rec=meta['elongation_frame'];basis=np.asarray(rec['basis_columns'])
    origin=np.asarray(src['source_asset']['source_origin_angstrom'])
    rotation=basis.T @ np.asarray(src['source_asset']['source_basis_columns'])
    transform=Matrix(rotation.tolist()).to_4x4()
    transform.translation=Vector(basis.T @ (origin-rec['origin_angstrom'])/10)
    for ch in 'ABCDEFGHIJKL':
        native=reference.objects['Intact Pol II chain '+ch]
        obj=bpy.data.objects.new('Overview · 5FLM '+ch,native.data)
        g.objects.link(obj);obj.parent=pol;obj.matrix_basis=transform
        obj.material_slots[0].link='OBJECT';obj.material_slots[0].material=mats[protein_role(ch)]
        obj['pdb']='5FLM';obj['source_chain']=ch
    pol.scale=(.40,)*3
    frames=times(*range(1,FATE_END,7),1,ESCAPE_END,73,133,140,141,154,180,183,209,232,234,235,267,FATE_END)
    states={f:strands(f) for f in frames}
    for f in frames:
        pol.location=states[f][1];pol.keyframe_insert(data_path='location',frame=f)
    xs=np.linspace(-10.6,10.6,300)
    for sign,role in ((1,'X'),(-1,'Y')):
        def locus(f):
            site=states[f][0][0]
            open_fraction=(1-smooth(f,232,267))*np.exp(-((xs-site)/1.12)**2)
            theta=xs*4.8
            y=.18*np.sin(.27*xs)+sign*.21*np.cos(theta)*(1-open_fraction)
            z=.26*np.sin(.38*xs)-.35+sign*(.22*np.sin(theta)+.39*open_fraction)
            return np.column_stack((xs,y,z))
        dna=animated_path('Gene DNA · '+role,frames,locus,.078,mats[role],g)
        dna['coordinate_evidence']='Schematic gene, not the experimental 5FLM template length'
    body=animated_path('Upstream RNA · coding transcript',frames,lambda f:states[f][2],.135,mats['rna'],g)
    body['entity_id']='upstream-rna'
    down=mats['rna'].copy();down.name='Downstream substrate · color follows cleavage'
    for channel in ('Base Color','Emission Color'):
        color=down.node_tree.nodes['Principled BSDF'].inputs[channel]
        for f,k in [(140,'F58D67'),(146,'F3C55F')]:
            color.default_value=(*ds.linear(k),1);color.keyframe_insert(data_path='default_value',frame=f)
    fragment=animated_path('Downstream RNA · polymerase-associated',frames,lambda f:states[f][3],.135,down,g)
    fragment['entity_id']='downstream-rna';visible(fragment,1,234)
    tail=animated_path('Poly(A) · grows from the upstream 3′ end',frames,lambda f:states[f][4],.105,mats['polyA'],g)
    tail['covalent_parent']='upstream-rna';visible(tail,155,FATE_END)
    cleavage=empty('CPSF73 · processing landmark',g)
    obj=ds.mesh('2I7V · CPSF73 catalytic segment',geo['cpsf_v'].tolist(),geo['cpsf_f'].tolist(),[mats['cpsf']],g)
    obj.parent=cleavage;obj['pdb']='2I7V'
    zinc=balls('2I7V · resolved zinc',geo['cpsf_zn'],.09,mats['ion'],g,cleavage)
    sulfate=balls('2I7V · crystallization sulfate',geo['cpsf_so4'],.085,mats['ion'],g,cleavage)
    cleavage.scale=(.41,)*3
    for f,loc in [(105,(-.5,-2.5,-4.0)),(136,(.05,-2.15,-2.41)),(146,(.05,-2.15,-2.41)),(169,(.75,-2.2,-3.65))]:
        cleavage.location=loc;cleavage.keyframe_insert(data_path='location',frame=f)
    for ob in [obj,zinc,sulfate]:visible(ob,105,170)
    a=np.linspace(0,2*math.pi,50,endpoint=False)
    loop=np.column_stack((.31*np.cos(a),np.zeros(len(a)),.31*np.sin(a)))
    xr=ds.path('XRN2 · 5′ degradation marker',loop,.072,mats['activity'],g,cyclic=True)
    xr['representation']='Activity marker, no molecular docking pose';visible(xr,177,236)
    for f in frames:
        xr.location=states[f][5]+[0,-.13,0];xr.keyframe_insert(data_path='location',frame=f)
    attach=np.asarray(rec['ctd_attachment_nm'])
    controls=attach+np.array([[0,0,0],[.5,-1.,1.],[2.2,-2.6,2.0],[3.1,-4.1,3.6],[4.85,-5.9,3.0],[6.25,-6.5,3.95]])
    ctd=ds.path('RPB1 CTD · schematic extension',controls,.15,mats['ctd'],g);ctd.parent=pol
    ctd['attachment']='5FLM RPB1 last modeled CA '+str(rec['rpb1_last_ca'])
    for index,offset in enumerate(controls[-2:]):
        phosphate=balls('CTD phosphate symbol '+str(index+1),[(0,0,0)],.48,mats['phosphate'],g,pol)
        start=274+index*7
        for f,location in [(1,offset),(start,offset),(start+20,offset+[.7,-.25,1.35])]:
            phosphate.location=location;phosphate.keyframe_insert(data_path='location',frame=f)
        visible(phosphate,30+index*8,start+22)
    # The reference mark stays on the promoter while the whole enzyme and RNA leave.
    promoter_mark=ds.path('Promoter position',[(x,-.42,.26*math.sin(.38*x)-.89) for x in np.linspace(-6.5,-4.9,16)],.10,mats['tbp'],g)
    visible(promoter_mark,1,ESCAPE_END)
    caption_set(shot,cam,fonts,lang,[(1,'Pol II 离开启动子，RNA 继续延伸','Pol II leaves the promoter as RNA grows'),
        (73,'延伸至转录本的 3′ 端加工区','Transcription reaches the 3′ processing region'),
        (108,'CPSF73 切开新生 RNA','CPSF73 cleaves the nascent RNA'),
        (154,'上游 RNA 加上 poly(A) 尾','Poly(A) polymerase extends the upstream RNA'),
        (182,'XRN2 从新露出的 5′ 端降解下游 RNA','XRN2 degrades the downstream RNA from its new 5′ end'),
        (235,'转录终止，Pol II 离开 DNA','Transcription terminates and Pol II leaves DNA'),
        (273,'FCP1 去磷酸化，Pol II 可再次起始','FCP1 removes CTD phosphate for Pol II reuse')],[(1,21.5)])
    shot.frame_set(1)
    return shot


def font_pair(args,lang):
    def load_font(given,names):
        if given:
            font=bpy.data.fonts.load(str(given.resolve()),check_existing=True)
            if not font.packed_file:font.pack()
            return font
        return next((f for f in bpy.data.fonts if f.packed_file and any(n in f.name.lower() for n in names)),None)
    zh=load_font(args.font_zh,('songti','simsun'))
    en=load_font(args.font_latin,('times new roman','timesnewroman'))
    if en is None or (lang=='zh' and zh is None):
        raise ValueError('Supply the missing Songti and/or Times New Roman font path for this project')
    return zh,en


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--font-zh',type=Path)
    parser.add_argument('--font-latin',type=Path)
    args=parser.parse_args(sys.argv[sys.argv.index('--')+1:])
    if args.output.exists():parser.error('Choose a fresh output .blend file')
    reference=bpy.context.scene
    blend=Path(bpy.data.filepath)
    src=json.loads(blend.with_name('manifest.json').read_text())
    if (src.get('scene')!='Two successive RNA additions C then U with slower pacing'
            or src['source_asset']['entry_url']!='https://www.rcsb.org/structure/5FLM'
            or (reference.frame_start,reference.frame_end,reference.render.fps,reference.render.fps_base)!=(1,552,24,1)
            or any('Intact Pol II chain '+c not in reference.objects for c in 'ABCDEFGHIJKL')):
        raise ValueError('Open the saved 24 fps two-addition 5FLM scene with its matching manifest')
    meta=json.loads((args.prepared/'manifest.json').read_text())
    archive=args.prepared/'geometry.npz'
    if (meta['kind']!='pol2-transcription-route-v1' or meta['config_sha256']!=digest(CONFIG)
            or meta['geometry_sha256']!=digest(archive)):
        raise ValueError('Prepared geometry and source scope differ')
    ref_appearance=src.get('appearance',{})
    profile=ref_appearance.get('profile') or json.loads((CONFIG.parent/'style-transcription-reference.json').read_text())
    lang=ref_appearance.get('language','zh' if any('songti' in f.name.lower() for f in bpy.data.fonts) else 'en')
    lang='zh' if lang.startswith('zh') else 'en'
    fonts=font_pair(args,lang)
    mats=render_mats()
    with np.load(archive,allow_pickle=False) as geo:
        start=promoter(geo,meta,reference,profile,fonts,lang,mats)
        overview=fates(geo,meta,reference,src,profile,fonts,lang,mats)
    reference.objects['RNA primer P residues 7–19']['entity_id']='upstream-rna'
    reference.objects['Incoming CTP core → RNA C20']['joins_entity']='upstream-rna'
    titles=('启动子识别与起始','启动子逃逸','两轮核苷酸加入','RNA 释放与聚合酶回收') if lang=='zh' else ('Promoter recognition','Promoter escape','Two nucleotide additions','RNA and polymerase fates')
    master=compose('Pol II · full teaching route',[(titles[0],start,1,PIC_END),(titles[1],overview,1,ESCAPE_END),
                   (titles[2],reference,1,552),(titles[3],overview,ESCAPE_END+1,FATE_END)],reference)
    master['pol2_cycle']=json.dumps({'source_blend_sha256':digest(blend),'config_sha256':digest(CONFIG),
        'geometry_sha256':meta['geometry_sha256'],'language':lang,'entities':['pol2','upstream-rna','downstream-rna'],
        'scope':json.loads(CONFIG.read_text())['scope']},ensure_ascii=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    bpy.ops.wm.save_as_mainfile(filepath=str(args.output.resolve()),compress=True)
    print(json.dumps({'master':master.name,'frames':master.frame_end,'chapters':json.loads(master['native_chapters'])},ensure_ascii=False))


if __name__=='__main__':main()
