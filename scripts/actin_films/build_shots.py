"""Build editable Blender chapters for two complete, source-backed actin films.

Molecular surfaces retain whole modeled chains. Motion is authored teaching
motion; source structures, constructs and schematic steps remain in metadata.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Matrix, Quaternion, Vector

sys.path.insert(0,str(Path(__file__).resolve().parent))
from project_context import ROOT, SKILL, FFMPEG, FONT, LATIN_FONT, OPTIONS
sys.path.insert(0,str(SKILL/'scripts'))
from detail_scene import Parts, collection, material, mesh, path, studio



OLD=ROOT/'assets/filaments'
STATES=ROOT/'assets/states'

def T(x=0,y=0,z=0):return Matrix.Translation((x,y,z))
def R(angle,axis='Y'):return Matrix.Rotation(math.radians(angle),4,axis)
def power(m,n):return Matrix(np.linalg.matrix_power(np.array(m),n))
def smooth(u):
    u=min(1,max(0,u));return u*u*u*(10+u*(-15+6*u))

def sample_detour(fraction,repair):
    if not repair or not repair['detour']:return np.zeros(6)
    data=repair['detour'];first=data['first_fraction'];last=data['last_fraction']
    if fraction<=first or fraction>=last:return np.zeros(6)
    u=(fraction-first)/(last-first);knots=data['knots']
    basis=[float(a<=u<b) for a,b in zip(knots,knots[1:])]
    for degree in range(1,4):
        updated=[]
        for i in range(len(basis)-1):
            left=knots[i+degree]-knots[i];right=knots[i+degree+1]-knots[i+1]
            updated.append((0 if left==0 else (u-knots[i])/left*basis[i])
                           +(0 if right==0 else (knots[i+degree+1]-u)/right*basis[i+1]))
        basis=updated
    controls=data.get('controls_nm_radians',data.get('controls_nm'))
    value=np.array(basis)@np.array(controls)
    return np.r_[value,np.zeros(6-len(value))]

class Shot:
    def __init__(self,config,width):
        self.c=config;self.n=config['frames'];self.events=[];self.objects=[];self.counter=0
        self.width=width;self.cached={};self.groups={}
        self.plans=json.loads((ROOT/'assets/approach-plans.json').read_text()) if (ROOT/'assets/approach-plans.json').exists() else {}
        repairs=ROOT/'assets/contact-repairs.json'
        self.repairs=json.loads(repairs.read_text()) if repairs.exists() else {}
        self.primary=[];self.labels=[];self.scene_scale=52
        bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
        bpy.data.orphans_purge(do_recursive=True)
        for group in list(bpy.data.collections):
            if not group.objects and not group.children:
                bpy.data.collections.remove(group)
        for n in ['Molecules','Schematics','Labels','Environment']:
            self.groups[n]=collection(n)
        palette={'teal':'4CACA9','blue':'4D83B5','gold':'E8B852','pink':'BF6A9C',
                 'green':'8DB788','red':'D98873','cream':'D5DCCB','dim':'45616B','ink':'ECF2F0'}
        self.mats={k:material(k,v,.53,0) for k,v in palette.items()}
        self.ink=bpy.data.materials.new('Emission typography');self.ink.use_nodes=True
        nodes=self.ink.node_tree.nodes;nodes.clear()
        emission=nodes.new('ShaderNodeEmission');emission.inputs['Color'].default_value=(.76,.88,.87,1)
        emission.inputs['Strength'].default_value=.85
        out=nodes.new('ShaderNodeOutputMaterial');self.ink.node_tree.links.new(emission.outputs[0],out.inputs['Surface'])
        self.diagram_ink=self.ink.copy();self.diagram_ink.name='Dark diagram typography'
        self.diagram_ink.node_tree.nodes.get('Emission').inputs['Color'].default_value=(.012,.025,.032,1)
        self.font=bpy.data.fonts.load(FONT);self.font.pack()
        self.latin_font=bpy.data.fonts.load(LATIN_FONT);self.latin_font.pack()

    def rig(self,name,parent=None,matrix=None):
        obj=bpy.data.objects.new(name,None);self.groups['Molecules'].objects.link(obj)
        obj.parent=parent
        if matrix is not None:obj.matrix_basis=matrix
        obj['motion_evidence']='schematic';return obj

    def asset(self,code,chain,name=None,color='teal',parent=None,mat=None,main=True):
        key=(code,chain)
        if key not in self.cached:
            p=ROOT/'assets'/code/(chain+'.npz')
            self.cached[key]=np.load(p,allow_pickle=False)
        d=self.cached[key]
        obj=mesh(name or code+' '+chain,d['vertices'].tolist(),d['faces'].tolist(),[self.mats[color]],self.groups['Molecules'])
        obj.parent=parent
        if mat is not None:obj.matrix_basis=mat
        obj['source_structure']=code;obj['source_chain']=chain
        obj['atoms_file']=str(ROOT/'assets'/code/(chain+'.npz'))
        obj['geometry_evidence']='experimental-coordinate-derived complete modeled chain envelope'
        if main:self.primary.append(obj)
        self.objects.append(obj);return obj

    def compound(self,code,chains,colors=None,name=None,matrix=None,main=True):
        rig=self.rig(name or code,matrix=matrix)
        obs=[self.asset(code,c,color=(colors or ['teal'])[i%len(colors or ['teal'])],parent=rig,main=main) for i,c in enumerate(chains)]
        return rig,obs

    def oldcap(self,site,parent=None):
        # Residue-number correspondence is recovered from cached source mmCIF
        # in the scientific preparation helper, not guessed inside Blender.
        tr=self.repairs['placements'][site]['matrix']
        rig=self.rig('CapZ complete heterodimer',parent=parent,matrix=Matrix(tr))
        rig['placement_evidence']='Source-guided bounded rigid adjustment to a different actin structure; not an observed complex'
        out=[]
        for c in ['G','H']:
            a=np.load(OLD/'cap'/('chain_'+c+'.npz'))
            obj=mesh('CapZ '+c,a['vertices'].tolist(),a['faces'].tolist(),[self.mats['red']],self.groups['Molecules'])
            obj.parent=rig;obj['source_structure']='8F8Q';obj['source_chain']=c
            obj['atoms_file']=str(OLD/'cap'/('chain_'+c+'.npz'));out.append(obj);self.primary.append(obj)
        return rig,out

    def actin(self,name,target=None,parent=None,tracked=False,g_fraction=0,asset='actin_morph.npz'):
        d=np.load(ROOT/'assets'/asset)
        obj=mesh(name,d['f_vertices'].tolist(),d['faces'].tolist(),[self.mats['gold' if tracked else 'teal']],self.groups['Molecules'])
        obj.shape_key_add(name='F reference')
        key=obj.shape_key_add(name='G reference guided approximation')
        key.data.foreach_set('co',d['g_vertices'].reshape(-1));key.value=g_fraction
        obj.parent=parent;obj.matrix_basis=target if target is not None else Matrix.Identity(4)
        obj['biological_id']=name;obj['geometry_evidence']='F source envelope / CA-guided G approximation; cross-isoform teaching morph'
        obj['atoms_file']=str(ROOT/'assets'/asset)
        self.primary.append(obj);self.objects.append(obj);return obj

    def profilin(self,name,parent=None,matrix=None):
        d=np.load(ROOT/'assets/profilin_on_morph.npz')
        obj=mesh(name,d['vertices'].tolist(),d['faces'].tolist(),[self.mats['pink']],self.groups['Molecules'])
        obj.parent=parent
        if matrix is not None:obj.matrix_basis=matrix
        obj['source_structure']='2BTF';obj['source_chain']='P'
        obj['atoms_file']=str(ROOT/'assets/profilin_on_morph.npz')
        self.primary.append(obj);self.objects.append(obj);return obj

    def frame(self,t):return min(self.n,max(1,1+round(t*(self.n-1))))
    def docking(self,obj,target,interval,offset,release=False,turn=False,far_shift=(0,0,0),repair=None):
        first,last=map(self.frame,interval);offset=Vector(offset)
        obj.rotation_mode='QUATERNION'
        pivot=Vector(np.mean([tuple(v) for v in obj.bound_box],axis=0)) if obj.type=='MESH' else Vector((0,0,0))
        side=offset.cross(Vector((0,0,1))).normalized()*1.2
        for f in [1]+list(range(first,last+1))+[self.n]:
            progress=smooth((f-first)/max(1,last-first))
            d=progress if release else 1-progress
            far=max(0,(d-.5)*2)
            travel=offset*d+side*math.sin(math.pi*far)**2+Vector(far_shift)*smooth(far)
            correction=sample_detour(d,repair)
            travel+=Vector(correction[:3])
            mat=T(*travel)@target
            angle=float(np.linalg.norm(correction[3:]))
            if angle>1e-12:
                anchor=Vector(repair['detour']['pivot_nm'])
                rotation=Quaternion(Vector(correction[3:]/angle),angle).to_matrix().to_4x4()
                mat=T(*travel)@T(*anchor)@rotation@T(*(-anchor))@target
            if turn:mat=mat@T(*pivot)@R(18*smooth(far),'Z')@T(*(-pivot))
            loc,quat,_=mat.decompose();obj.location=loc;obj.rotation_quaternion=quat
            obj.keyframe_insert('location',frame=f);obj.keyframe_insert('rotation_quaternion',frame=f)
        self.events.append({'object':obj.name,'first':first,'last':last,'action':'source-atom screened illustrative release' if release else 'source-atom screened illustrative approach'})
    def move(self,obj,start,end,interval,arc=(0,0,0),action='authored movement'):
        first,last=map(self.frame,interval)
        a,qa,_=start.decompose();b,qb,_=end.decompose()
        obj.rotation_mode='QUATERNION'
        for f in [1]+list(range(first,last+1))+[self.n]:
            u=smooth((f-first)/max(1,last-first))
            obj.location=a.lerp(b,u)+Vector(arc)*math.sin(math.pi*u)**2
            obj.rotation_quaternion=qa.slerp(qb,u);obj.scale=(1,1,1)
            obj.keyframe_insert('location',frame=f);obj.keyframe_insert('rotation_quaternion',frame=f)
        self.events.append({'object':obj.name,'first':first,'last':last,'action':action})

    def morph(self,obj,start,end,interval):
        key=obj.data.shape_keys.key_blocks[1]
        first,last=map(self.frame,interval)
        for f in [1]+list(range(first,last+1))+[self.n]:
            key.value=start+(end-start)*smooth((f-first)/max(1,last-first));key.keyframe_insert('value',frame=f)
        self.events.append({'object':obj.name,'first':first,'last':last,'action':'schematic G/F state morph'})

    def visibility(self,objects,first=0,last=1):
        if not isinstance(objects,(list,tuple)):objects=[objects]
        lo,hi=self.frame(first),self.frame(last)
        for obj in objects:
            for f,hide in [(1,lo>1),(max(1,lo-1),lo>1),(lo,False),(hi,False),(hi+1,hi<self.n)]:
                if f>self.n:continue
                obj.hide_render=hide;obj.keyframe_insert('hide_render',frame=f)

    def text(self,body,position,size=.85,color=None,align='CENTER'):
        data=bpy.data.curves.new(body,'FONT');data.body=body;data.font=self.font;data.size=size
        # The regular TNR face occupies Blender's alternate-face slot so mixed
        # Chinese/Latin strings retain exact layout without synthetic bolding.
        data.font_bold=self.latin_font
        for char,fmt in zip(body,data.body_format):fmt.use_bold=ord(char)<0x400 or 0x2010<=ord(char)<=0x22ff
        data.align_x=align;data.align_y='CENTER';data.space_line=1.15
        obj=bpy.data.objects.new(body,data);self.groups['Labels'].objects.link(obj)
        obj.location=position;obj.rotation_euler=(math.pi/2,0,0)
        data.materials.append(color or self.ink);self.labels.append(obj);return obj

    def line(self,points,color='cream',width=.08,name='Schematic connection'):
        return path(name,points,width,self.mats[color],self.groups['Schematics'])

    def arrow(self,a,b,color='cream'):
        a,b=Vector(a),Vector(b);v=(b-a).normalized()
        self.line([a,b],color)
        wing=v.cross(Vector((0,1,0)))
        self.line([b-v*.65+wing*.4,b,b-v*.65-wing*.4],color)

    def node(self,label,x,z,color='blue',width=7.3):
        bpy.ops.mesh.primitive_cube_add(size=1,location=(x,.5,z))
        obj=bpy.context.object;obj.name='Diagram node '+label
        obj.dimensions=(width,.35,2.2);bpy.ops.object.transform_apply(location=False,rotation=False,scale=True)
        for c in list(obj.users_collection):c.objects.unlink(obj)
        self.groups['Schematics'].objects.link(obj);obj.data.materials.append(self.mats[color])
        mod=obj.modifiers.new('Diagram corner radius','BEVEL');mod.width=.2;mod.segments=3
        self.text(label,(x,-40,z),1.2,color=self.diagram_ink);return obj

    def token(self,label,pos,color='gold',radius=.75):
        p=Parts();p.ball((0,0,0),radius,seg=20,rings=12)
        obj=p.object('Schematic '+label,[self.mats[color]],self.groups['Schematics']);obj.location=pos
        lab=self.text(label,(0,-40,0),radius*.8,color=self.diagram_ink)
        lab.parent=obj;lab.location=(0,-40,0)
        obj['representation']='chemical/signaling symbol, not atomic scale'
        return obj,lab

    def step(self,code):return Matrix(json.loads((ROOT/'assets'/code/'manifest.json').read_text())['step'])

    def center(self,objects,offset=(0,0,0),rotation=None):
        # Rigid composition transform; protein sizes remain unchanged.
        bpy.context.scene.frame_set(self.n)
        bpy.context.view_layer.update()
        points=[obj.matrix_world@Vector(v) for obj in objects for v in obj.bound_box]
        rotation=rotation or Matrix.Identity(4)
        rotated=np.array([rotation@Vector(p) for p in points])
        center=(np.min(rotated,axis=0)+np.max(rotated,axis=0))/2
        rig=self.rig('Camera composition rig')
        roots=set(objects)
        for obj in objects:
            while obj.parent:
                obj=obj.parent
            roots.add(obj)
        roots={obj for obj in roots if obj.parent is None}
        for obj in roots:obj.parent=rig
        rig.matrix_basis=T(*(np.array(offset)-center))@rotation
        bpy.context.view_layer.update()
        return rig

    def finish(self):
        # Linear interpolation of already smooth, frame-sampled authored paths.
        for action in bpy.data.actions:
            for layer in action.layers:
                for strip in layer.strips:
                    for bag in strip.channelbags:
                        for curve in bag.fcurves:
                            for k in curve.keyframe_points:
                                k.interpolation='CONSTANT' if 'hide_render' in curve.data_path else 'LINEAR'
        profile=json.loads((SKILL/'assets/style-molecular-reference.json').read_text())
        profile['background']='172B35';profile['bump']=0;profile['roughness']=.53
        profile['exposure']=.3;profile['compositor']['contrast']=1
        args=argparse.Namespace(engine=OPTIONS['engine'],width=self.width,samples=OPTIONS['samples'],fps=self.c['fps'],frames=self.n,compositor=OPTIONS['compositor'])
        studio(profile,args,{'position':(0,-100,0),'target':(0,0,0),'scale':self.scene_scale,'light_scale':1.4})
        scene=bpy.context.scene;scene.render.resolution_y=2*round(self.width*9/32)
        scene.render.image_settings.color_mode='RGB';scene.render.image_settings.compression=25
        scene.render.use_file_extension=True
        scene.render.film_transparent=True
        scene.render.fps=self.c['fps']
        scene.render.fps_base=1
        if self.c['id']=='a03':
            # An explicit editorial dip separates Arp recruitment from the
            # selected monomer-pool view. Whole units enter view at the cut.
            for tree in [t for t in bpy.data.node_groups if t.name.startswith('COMP_')]:
                over=next(n for n in tree.nodes if n.bl_idname=='CompositorNodeAlphaOver')
                for t,value in [(0,1),(.285,1),(.300,0),(.307,0),(.324,1),(1,1)]:
                    over.inputs['Factor'].default_value=value;over.inputs['Factor'].keyframe_insert('default_value',frame=self.frame(t))
        if hasattr(scene,'eevee') and hasattr(scene.eevee,'taa_render_samples'):
            scene.eevee.taa_render_samples=OPTIONS['samples']
        W=self.scene_scale
        readable=self.c['id'] in {'a03','a04','a08','c04','f02','f04','f07'}
        self.text(self.c['title'],(0,-75,W*.244),W*(.029 if readable else .025))
        for i,beat in enumerate(self.c['beats']):
            obj=self.text(beat,(0,-75,W*.192),W*(.026 if readable else .024))
            boundaries=self.c.get('beat_boundaries',[0,1/3,2/3,1])
            self.visibility(obj,boundaries[i],1 if i==2 else boundaries[i+1]-1/self.n)
        # A few whole source-derived background proteins give depth; their
        # illustrative placement/density is explicitly not a cell measurement.
        data=np.load(ROOT/'assets/ubiquitin/chain_A.npz')
        rng=np.random.default_rng(901)
        for i in range(12):
            x=rng.uniform(-W*.52,W*.52);z=rng.uniform(-W*.16,W*.16)
            if abs(x)<W*.32:continue
            obj=mesh('Context ubiquitin '+str(i),data['vertices'].tolist(),data['faces'].tolist(),[self.mats['dim']],self.groups['Environment'])
            obj.location=(x,24+rng.uniform(0,8),z);obj.rotation_euler=rng.uniform(0,math.tau,3)
            obj['source_structure']='1UBQ';obj['placement']='schematic environmental depth, not measured density'
        scene.frame_set(1);bpy.context.view_layer.update()
        folder=ROOT/'shots'/self.c['id'];folder.mkdir(exist_ok=True)
        scene.render.filepath=str(folder/'preview.png')
        scene['chapter_id']=self.c['id'];scene['biological_scope']=self.c['subtitle']
        scene['animation_evidence']='Authored teaching motion; not MD or observed reaction kinetics'
        scene['events']=json.dumps(self.events)
        bpy.ops.wm.save_as_mainfile(filepath=str(folder/'scene.blend'))
        record={**self.c,'events':self.events,'scene_scale_nm':self.scene_scale,
                'source_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'objects':len(scene.objects),'representation_limits':json.loads((ROOT/'assets/manifest.json').read_text())}
        (folder/'manifest.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
        return folder


def input_scene(s,formin=False):
    s.scene_scale=58
    labs=['GPCR','Gα13','p115RhoGEF','RhoA'] if formin else ['PDGF / 受体','PI3K / PIP3','Rac-GEF','Rac1']
    xs=[-21,-7,7,21]
    for x,lab in zip(xs,labs):s.node(lab,x,5,width=11)
    for a,b in zip(xs,xs[1:]):s.arrow((a+5.5,0,5),(b-5.7,0,5))
    signal,_=s.token('信号',(-21,-2,8),'green',.65)
    s.move(signal,T(-21,-2,8),T(21,-2,8),(.07,.32),(0,0,1.3),'upstream signal diagram')
    s.node(labs[-1]+'–GDP',-11,-3,'blue',10)
    s.node(labs[-1]+'–GTP',11,-3,'gold',10)
    s.arrow((-5,0,-3),(5,0,-3))
    s.text('GEF：核苷酸交换', (0,-1,-.4),.86)
    gdp,_=s.token('GDP',(-11,-1,-6),'blue')
    atp,_=s.token('GTP',(21,-1,-7),'gold')
    s.move(gdp,T(-11,-1,-6),T(-22,-1,-8),(.36,.49),(0,0,-1),'GDP release diagram')
    s.move(atp,T(21,-1,-7),T(11,-1,-6),(.51,.66),(0,0,-1),'new GTP binding diagram')

def wrc_scene(s):
    s.scene_scale=78
    rig,obs=s.compound('7USE',list('ABCDEFG'),['teal','blue','gold','green','pink','red','cream'])
    for obj,side in zip(obs[-2:],[-1,1]):
        s.docking(obj,Matrix.Identity(4),(.14,.40),Vector(s.plans['rac_'+('F' if side==-1 else 'G')]['offset'])*.4)
    s.center(obs,offset=(-2,0,-1),rotation=R(145,'Y')@R(15,'Z'))
    s.line([(3,-4,-3),(8,-4,-3),(10,-4,1)],'gold',.12,'WCA unresolved flexible segment schematic')
    wc=s.node('WCA 可用',11,2,'gold',7.7)
    s.visibility([wc]+[s.labels[-1]],.46,1)
    s.text('五亚基 WRC',(-17,-18,4),1.1)

def supply_scene(s):
    s.scene_scale=40
    r,obs=s.compound('2BTF',['A','P'],['gold','pink'])
    centers=[np.load(ROOT/'assets/2BTF'/f'{c}.npz')['atoms'].mean(0) for c in ['A','P']]
    rot=Vector(centers[1]-centers[0]).rotation_difference(Vector((-.7,-1,.1))).to_matrix().to_4x4()
    s.center(obs,offset=(-3,0,0),rotation=rot)
    obs[0]['biological_id']='actin_001'
    s.text('actin_001',(-4,-8,4.2),.8)
    s.text('profilin',(-8,-8,-3.6),.8)
    for lab,start,end,color,interval in [
        ('ADP',(-1,-8,-4),(-13,-8,-6),'blue',(.12,.30)),
        ('ATP',(13,-8,-6),(-1,-8,-4),'gold',(.40,.60))]:
        obj,_=s.token(lab,start,color,.85);s.move(obj,T(*start),T(*end),interval,(0,-1,1.5),'nucleotide exchange symbol')
    s.node('Tβ4 封存池',11,3,'blue',8)
    s.arrow((9,0,1.6),(5,0,-.4))

def branch_setup(s,mother=True,daughter=True,arp=True,parent=None):
    groups={}
    for label,chains,colors in [('mother',list('NOPQRS'),['blue','teal']),('arp',list('ABCDEFG'),['pink','red','green','cream']),('daughter',list('HIJK'),['teal','blue'])]:
        if not {'mother':mother,'arp':arp,'daughter':daughter}[label]:continue
        rig=s.rig(label,parent=parent);objs=[]
        for i,c in enumerate(chains):objs.append(s.asset('7TPT',c,color=colors[i%len(colors)],parent=rig))
        groups[label]=(rig,objs)
    return groups

def nucleate_branch(s):
    s.scene_scale=90
    record=json.loads((ROOT/'assets/branch-path-repairs.json').read_text())
    if not record['passed']:raise ValueError('Branch path preparation failed; prepare this project again')
    groups=branch_setup(s)
    arp,arpobs=groups['arp']
    path_record=record['motions']['arp']
    s.docking(arp,Matrix.Identity(4),record['intervals']['arp'],path_record['offset'],repair=path_record)
    for i,o in enumerate(groups['daughter'][1]):
        key='branch_'+list('HIJK')[i];path_record=record['motions'][key]
        s.docking(o,Matrix.Identity(4),record['intervals'][key],path_record['offset'],far_shift=path_record['far_shift'],repair=path_record)
        s.visibility(o,record['daughter_reveal'],1)
    allobs=sum([v[1] for v in groups.values()],[])
    comp=s.center(allobs,offset=(0,0,-1),rotation=R(20))
    s.text('已有母丝',(-20,-40,-8),1.8)
    s.text('Arp2/3 七亚基',(-19,-40,2),1.8)
    wca,wca_label=s.token('WCA',(-24,-40,12.2),'gold',2.2)
    wca_label.data.materials[0]=s.ink
    s.move(wca,T(-24,-40,12.2),T(-33,-40,12.2),(.67,.80),(0,0,1),'NPF release symbol')

def branch_growth(s,reentry=False):
    s.scene_scale=60
    groups=branch_setup(s,mother=False,arp=False)
    record=json.loads((ROOT/'assets/branch-growth.json').read_text())
    motion_file=ROOT/'assets/branch_growth/motion.npz'
    if not record['passed'] or hashlib.sha256(motion_file.read_bytes()).hexdigest()!=record['motion_sha256']:
        raise ValueError('Branch growth preparation failed or its motion file changed; prepare this project again')
    if hashlib.sha256((ROOT/'assets/actin_morph.npz').read_bytes()).hexdigest()!=record['actin_asset_sha256']:
        raise ValueError('Branch growth actin asset changed; prepare this project again')
    with np.load(motion_file) as data:motion={key:data[key] for key in data.files}
    added=[]
    for i in range(2):
        obj=s.actin('actin_001' if i==0 else 'actin_002',tracked=i==0,g_fraction=1)
        obj.rotation_mode='QUATERNION';added.append(obj)
        for f in range(1,s.n+1):
            index=(f-1)/(s.n-1)*(len(motion['matrices'])-1);lo=int(index);hi=min(lo+1,len(motion['matrices'])-1);q=index-lo
            a,qa,_=Matrix(motion['matrices'][lo,i]).decompose();b,qb,_=Matrix(motion['matrices'][hi,i]).decompose()
            obj.location=a.lerp(b,q);obj.rotation_quaternion=qa.slerp(qb,q)
            obj.keyframe_insert('location',frame=f);obj.keyframe_insert('rotation_quaternion',frame=f)
            key=obj.data.shape_keys.key_blocks[1];key.value=(1-q)*motion['g_fractions'][lo,i]+q*motion['g_fractions'][hi,i]
            key.keyframe_insert('value',frame=f)
    s.events.append({'object':'Two daughter-actin additions','first':1,'last':s.n,
        'action':'Source-fitted endpoints and coupled approach/G-F interpolation',
        'limits':record['limits'],'motion_preparation':str(ROOT/'assets/branch-growth.json')})
    s.center(sum([v[1] for v in groups.values()],[])+added,offset=(-3,0,2.5),rotation=R(75))
    s.text('开放的倒钩端 +',(14,-40,-11.8),1.4)
    s.text('金色：跟踪单体',(-16,-40,-11.8),1.4)

def autoinhib(s):
    s.scene_scale=38
    record=json.loads((ROOT/'assets/actin-path-repairs.json').read_text())
    if not record['passed']:raise ValueError('Regulatory release preparation failed; prepare this project again')
    rig,obs=s.compound('2BAP',list('ABCD'),['teal','blue','gold','gold'])
    for i,o in enumerate(obs[2:]):
        path_record=record['motions']['dad_'+list('CD')[i]]
        s.docking(o,Matrix.Identity(4),(.23,.43),path_record['offset'],release=True,repair=path_record)
    s.center(obs,offset=(-2,0,0),rotation=R(60,'Y')@R(-90,'X')@R(-20,'Z'))
    s.text('RhoA–GTP',(11,-14,3),1.2)
    s.arrow((7,0,1.5),(2.5,0,.7),'green')
    s.text('DID / 二聚化调控片段',(-7,-40,5.6),.95)
    s.text('金色：DAD 抑制片段',(-8,-40,-8.3),1.05)

def formin_setup(s,actin_chains='ABCDEF'):
    fr,heads=s.compound('9B27',['G','H'],['pink','green'],name='FH2 dimer')
    core=[s.asset('9B27',c,color='blue' if i%2 else 'teal') for i,c in enumerate(actin_chains)]
    return core,fr,heads

def formin_seed(s):
    s.scene_scale=64
    repairs=json.loads((ROOT/'assets/fh2-seed-repairs.json').read_text())
    core=[s.asset('9B27',c,color='blue' if i%2 else 'teal') for i,c in enumerate('DEF')]
    heads={}
    for chain,partner in [('G','H'),('H','G')]:
        obj=s.asset('fh2_seed',chain,name='FH2 body '+chain+' / partner lasso '+partner,
                    color='pink' if chain=='G' else 'green')
        obj.data.materials.append(s.mats['green' if chain=='G' else 'pink'])
        obj.data.polygons.foreach_set('material_index',s.cached[('fh2_seed',chain)]['material_indices'].astype(np.int32))
        obj['source_structure']='9B27';obj['source_chain']=chain+' body; '+partner+' lasso'
        obj['geometry_evidence']='All original source vertices and faces retained; only existing unresolved-linker gaps separate rigid groups'
        heads[chain]=obj
    for chain,obj in zip('DEF',core):
        record=repairs['motions'][chain]
        s.docking(obj,Matrix.Identity(4),repairs['intervals'][chain],record['offset'],repair=record)
    head=heads['H'];offset=Vector(repairs['motions']['H']['offset'])
    side=offset.cross(Vector((0,0,1))).normalized()*1.2
    for f in range(1,s.n+1):
        u=(f-1)/(s.n-1);first,last=repairs['opening'];a,b=repairs['closing']
        d=smooth((u-first)/(last-first))*(1-smooth((u-a)/(b-a)));far=max(0,(d-.5)*2)
        correction=sample_detour(d,repairs['motions']['H'])
        head.location=offset*d+side*math.sin(math.pi*far)**2+Vector(correction[:3])
        head.keyframe_insert('location',frame=f)
    linkers=[];linker_checks=[]
    for record in repairs['linkers']:
        obj=path('Unresolved FH2 linker '+record['chain'],[(0,0,0)]*21,.055,
                 s.mats['pink' if record['chain']=='G' else 'green'],s.groups['Schematics'])
        obj['geometry_evidence']=record['representation'];linkers.append(obj);lengths=[]
        obj['source_linker_contour_nm']=record['contour_upper_bound_nm']
        obj['linker_anchor_start']={'object':heads[record['start_lobe']].name,'position':record['start']}
        obj['linker_anchor_end']={'object':heads[record['end_lobe']].name,'position':record['end']}
        for f in range(1,s.n+1):
            bpy.context.scene.frame_set(f);bpy.context.view_layer.update()
            a=heads[record['start_lobe']].matrix_basis@Vector(record['start'])
            b=heads[record['end_lobe']].matrix_basis@Vector(record['end'])
            radial=Vector((0,(a.y+b.y)/2,(a.z+b.z)/2)).normalized()
            points=[a.lerp(b,t/20)+radial*(.4*math.sin(math.pi*t/20)) for t in range(21)]
            length=sum((y-x).length for x,y in zip(points,points[1:]));lengths.append(length)
            if length>record['contour_upper_bound_nm']:
                raise ValueError('Schematic FH2 linker exceeds its residue contour bound')
            for point,position in zip(obj.data.splines[0].points,points):
                point.co=(*position,1);point.keyframe_insert('co',frame=f)
        linker_checks.append({'chain':record['chain'],'maximum_curve_length_nm':max(lengths),
                              'contour_upper_bound_nm':record['contour_upper_bound_nm']})
    s.events.append({'object':'FH2 seed assembly','action':'Sequential intact monomers with one rigid functional lobe opening and returning',
                     'limits':repairs['limits'],'unresolved_linker_checks':linker_checks})
    s.center(core+list(heads.values())+linkers,offset=(-1,0,0),rotation=R(-12))
    s.text('FH2 二聚体',(17,-40,5.5),1.4)
    free=s.text('游离单体',(-12,-40,-11),1.4);s.visibility(free,0,.48)
    seed=s.text('极性短丝',(-9,-40,-8),1.4);s.visibility(seed,.48,1)

def formin_growth(s,reentry=False):
    s.scene_scale=72
    folder=ROOT/'assets/fh2_growth'
    record=json.loads((ROOT/'assets/fh2-growth.json').read_text())
    if not record['passed'] or hashlib.sha256((folder/'motion.npz').read_bytes()).hexdigest()!=record['motion_sha256']:
        raise ValueError('FH2 growth preparation failed or its motion file changed; prepare this project again')
    with np.load(folder/'motion.npz') as data:motion={key:data[key] for key in data.files}
    core=[s.asset('9B27',c,color='blue' if i%2 else 'teal') for i,c in enumerate('ABCDEF')]
    monomers=[];profilins=[];heads={}
    for i in range(2):
        obj=s.actin('actin_001' if i==0 else 'actin_002',tracked=i==0,g_fraction=1,asset='fh2_growth/actin.npz')
        obj['geometry_evidence']=record['conformers']['actin']['limits']
        monomers.append(obj);prof=s.profilin('profilin cycle '+str(i+1),parent=obj)
        prof.data.materials[0]=s.mats['cream'];profilins.append(prof)
    for chain,partner in [('G','H'),('H','G')]:
        data=np.load(folder/(chain+'.npz'))
        mats=[s.mats['pink' if chain=='G' else 'green'],s.mats['green' if chain=='G' else 'pink']]
        obj=mesh('FH2 body '+chain+' / partner lasso '+partner,data['f_vertices'].tolist(),data['faces'].tolist(),mats,s.groups['Molecules'])
        obj.data.polygons.foreach_set('material_index',data['material_indices'].astype(np.int32))
        obj.shape_key_add(name='Source functional lobe')
        alternate=obj.shape_key_add(name='Source-guided alternate lobe')
        alternate.data.foreach_set('co',data['g_vertices'].reshape(-1))
        obj['atoms_file']=str(folder/(chain+'.npz'));obj['source_structure']='9B27'
        obj['geometry_evidence']='Complete modeled body and partner lasso, source-guided alternate conformation; unresolved linkers are separate schematic curves'
        heads[chain]=obj;s.primary.append(obj);s.objects.append(obj)
    objects=monomers+list(heads.values())
    for obj in objects:obj.rotation_mode='QUATERNION'
    for f in range(1,s.n+1):
        index=(f-1)/(s.n-1)*(len(motion['matrices'])-1);lo=int(index);hi=min(lo+1,len(motion['matrices'])-1);q=index-lo
        for i,obj in enumerate(objects):
            a,qa,_=Matrix(motion['matrices'][lo,i]).decompose();b,qb,_=Matrix(motion['matrices'][hi,i]).decompose()
            obj.location=a.lerp(b,q);obj.rotation_quaternion=qa.slerp(qb,q)
            obj.keyframe_insert('location',frame=f);obj.keyframe_insert('rotation_quaternion',frame=f)
            key=obj.data.shape_keys.key_blocks[1];key.value=(1-q)*motion['g_fractions'][lo,i]+q*motion['g_fractions'][hi,i]
            key.keyframe_insert('value',frame=f)
        for i,obj in enumerate(profilins):
            obj.location=(1-q)*motion['profilin_shifts'][lo,i]+q*motion['profilin_shifts'][hi,i]
            obj.keyframe_insert('location',frame=f)
    linkers=[]
    for item in record['linkers']:
        obj=path('Unresolved FH2 linker '+item['chain'],[(0,0,0)]*21,.055,
            s.mats['pink' if item['chain']=='G' else 'green'],s.groups['Schematics'])
        obj['geometry_evidence']=item['representation'];obj['source_linker_contour_nm']=item['contour_upper_bound_nm']
        for label in ['start','end']:
            anchor=item[label+'_anchor']
            obj['linker_anchor_'+label]={'object':heads[anchor['lobe']].name,
                'position':anchor['position'],'alternate_position':anchor['alternate_position']}
        for f in range(1,s.n+1):
            bpy.context.scene.frame_set(f);bpy.context.view_layer.update();anchors=[]
            for label in ['start','end']:
                anchor=item[label+'_anchor'];head=heads[anchor['lobe']];g=head.data.shape_keys.key_blocks[1].value
                local=Vector(anchor['position']).lerp(Vector(anchor['alternate_position']),g)
                anchors.append(head.matrix_basis@local)
            a,b=anchors;radial=Vector((0,(a.y+b.y)/2,(a.z+b.z)/2)).normalized()
            points=[a.lerp(b,t/20)+radial*(.4*math.sin(math.pi*t/20)) for t in range(21)]
            if sum((y-x).length for x,y in zip(points,points[1:]))>item['contour_upper_bound_nm']:
                raise ValueError('FH2 growth linker exceeds its residue contour bound')
            for point,position in zip(obj.data.splines[0].points,points):
                point.co=(*position,1);point.keyframe_insert('co',frame=f)
        linkers.append(obj)
    s.events.append({'object':'Coupled FH2 elongation','first':1,'last':s.n,
        'action':'Two intact monomers, profilin release and alternating complete functional lobes',
        'limits':record['limits'],'motion_preparation':str(ROOT/'assets/fh2-growth.json')})
    s.center(core+objects+linkers,offset=(-3,0,0),rotation=R(-50,'X')@R(-15))
    s.text('FH2 / 倒钩端',(14,-40,-14),1.65)
    s.text('金色：跟踪单体',(-15,-40,-14),1.65)
    # FH1 has no resolved coordinates in this source. A timed graphic cue
    # indicates supply via profilin without inventing an atomic tether pose.
    center=Vector(np.load(ROOT/'assets/profilin_on_morph.npz')['vertices'].mean(0))
    for cycle,prof in enumerate(profilins):
        label=s.text('FH1 → profilin 供料',(26,-40,10),1.65)
        shaft=path('FH1 supply cue '+str(cycle),[(0,0,0)]*21,.055,s.ink,s.groups['Schematics'])
        tip=path('FH1 supply arrow '+str(cycle),[(0,0,0)]*3,.055,s.ink,s.groups['Schematics'])
        for obj in [shaft,tip]:obj['geometry_evidence']='Graphic functional supply cue; not a physical bond or resolved FH1 structure'
        for f in range(1,s.n+1):
            bpy.context.scene.frame_set(f);bpy.context.view_layer.update()
            point=prof.matrix_world@center;end=Vector((point.x,-40,point.z+1.6))
            start=Vector((26,-40,8.7));middle=(start+end)*.5+Vector((3,0,0))
            points=[start*(1-t/20)**2+middle*(2*(1-t/20)*(t/20))+end*(t/20)**2 for t in range(21)]
            for control,position in zip(shaft.data.splines[0].points,points):
                control.co=(*position,1);control.keyframe_insert('co',frame=f)
            direction=(end-points[-2]).normalized();wing=direction.cross(Vector((0,1,0)))
            for control,position in zip(tip.data.splines[0].points,[end-direction*.55+wing*.3,end,end-direction*.55-wing*.3]):
                control.co=(*position,1);control.keyframe_insert('co',frame=f)
        s.visibility([label,shaft,tip],cycle*.5,cycle*.5+.12)

def chemistry(s):
    s.scene_scale=44
    source=json.loads((STATES/'manifest.json').read_text())
    for i,record in enumerate(source['states']):
        x=(i-1)*12
        d=np.load(STATES/(record['id']+'.npz'))
        o=mesh(record['id']+' whole F subunit',d['vertices'].tolist(),d['faces'].tolist(),[s.mats['teal']],s.groups['Molecules'])
        o.location=(x,0,1.3);o['source_structure']=record['id'];s.primary.append(o)
        o['atoms_file']=str(STATES/(record['id']+'.npz'))
        s.text(['ATP 状态类似物','ADP · Pi','ADP'][i],(x,-8,5.7),.88)
        adp=next(l for l in record['ligands'] if l['ccd']=='ADP');center=np.mean(adp['positions_nm'],axis=0)
        for ligand in record['ligands']:
            parts=Parts()
            for pos,el in zip(ligand['positions_nm'],ligand['elements']):
                parts.ball((np.array(pos)-center)*2.1, .18 if el=='Mg' else .14,seg=10,rings=7)
            col='green' if ligand['ccd']=='MG' else ('pink' if ligand['ccd']=='BEF' else 'gold' if ligand['ccd']=='PO4' else 'cream')
            obj=parts.object(record['id']+' '+ligand['ccd'],[s.mats[col]],s.groups['Molecules'])
            obj.location=(x,-1,-4.1)
            obj['source_structure']=record['id'];obj['ccd']=ligand['ccd'];obj['magnification']=2.1
            if ligand['ccd']=='PO4':
                s.move(obj,T(x,-1,-4.1),T(x+4,-1,-6),(.56,.72),(0,-2,1),'Pi release in an enlarged ligand-only inset')
        s.text(['ADP + BeF3 + Mg','ADP + Pi + Mg','ADP + Mg'][i],(x,-9,-6.8),.82)
    s.arrow((-7,-3,1.5),(-5,-3,1.5),'gold');s.arrow((5,-3,1.5),(7,-3,1.5),'gold')
    s.text('配体放大 ×2.1',(0,-12,-9.1),.9)

def network(s,formin=False):
    s.scene_scale=80 if formin else 110
    if formin:
        for j in range(3):
            rig=s.rig('Complete filament '+str(j),matrix=T(0,5*j,(j-1)*6))
            for i,c in enumerate('ABCDEF'):
                obj=s.asset('9B27',c,color='teal' if j!=1 else 'blue',parent=rig)
            s.move(rig,T(0,5*j,(j-1)*10),T(0,5*j,(j-1)*6),(.10,.33),(0,0,0),'filament organization schematic')
        for x in [-4,2,8]:s.line([(x,-1,-6),(x+1,-1,6)],'cream',.12,'Schematic crosslinker connection')
        s.node('RhoA',-24,7,'green',8);s.node('mDia1',-24,0,'pink',8);s.node('ROCK',-24,-7,'red',8)
        s.arrow((-24,0,5.5),(-24,0,1.4));s.line([(-28,0,7),(-31,0,7),(-31,0,-7),(-28,0,-7)],'green')
        s.text('交联 / 极性 / 肌球蛋白 / 黏附条件',(4,-15,10),1)
    else:
        for j,(dx,dz,theta) in enumerate([(-14,-6,-25),(9,-5,25)]):
            root=s.rig('Network branch '+str(j),matrix=T(dx,6*j,dz)@R(theta))
            groups=branch_setup(s,parent=root)
            if j==1:
                fit=Matrix(json.loads((ROOT/'assets/actin-target-fits.json').read_text())['branch_K'])
                cap,_=s.oldcap('branch',parent=root)
                endpoint=cap.matrix_basis.copy()
                offset=fit.to_3x3()@Vector(s.plans['cap']['offset'])*.6
                s.docking(cap,endpoint,(.14,.37),offset,repair=s.repairs['motions']['cap_branch'])
        # Membrane shown as a patch of complete lipid symbols; no protein cuts.
        p=Parts()
        for x in np.linspace(-26,26,49):
            for y in [-2,0,2]:
                for sign in [-1,1]:
                    p.ball((x,y,16.5+sign*.65),.36,seg=8,rings=5)
                    p.tube([(x,y,16.5+sign*.42),(x+.15,y,16.5)],.08)
        mem=p.object('Membrane lipids schematic',[s.mats['dim']],s.groups['Schematics'])
        s.move(mem,Matrix.Identity(4),T(0,0,1.5),(.4,.66),(0,0,0),'membrane shape output, not simulated force')
        s.text('细胞膜',(0,-15,14),1.9)
        s.text('CapZ 封端',(23,-18,-15),2.0)
        s.text('分支接头',(-24,-18,-15),2.0)

def stop_signal(s):
    s.scene_scale=45
    for x,l,c in [(-13,'Rac1–GTP','gold'),(13,'Rac1–GDP','blue')]:s.node(l,x,2,c,10)
    s.arrow((-7,0,2),(7,0,2));s.text('GAP 促进水解',(0,-2,5.5),1)
    pi,_=s.token('Pi',(11,-2,1),'green',.8)
    s.move(pi,T(11,-2,1),T(16,-2,-4),(.24,.43),(1,0,0),'signal Pi release symbol')
    s.node('新成核输入减少',0,-2,'pink',15)
    s.text('已有细丝仍在，等待后续周转',(0,-10,-6),.92)

def cofilin_control(s):
    s.scene_scale=46
    rig,obs=s.compound('6VAO',['I'],['pink']);s.center(obs,offset=(0,0,-1),rotation=R(20,'Z'))
    s.node('PAK / ROCK',-15,5,'blue',10);s.node('LIMK',-3,5,'red',8);s.node('SSH',15,5,'green',8)
    s.arrow((-9,0,5),(-7.3,0,5));s.arrow((-2,0,3.6),(0,0,1.5),'red');s.arrow((13,0,3.6),(4,0,1),'green')
    p,_=s.token('P',(-5,-6,0),'gold',.7)
    s.move(p,T(-5,-6,0),T(1,-6,0),(.12,.28),(0,0,1),'Ser3 phosphorylation symbol')
    # Distinct ligand object represents the released phosphate after the switch.
    q,_=s.token('Pi',(1,-6,0),'green',.7)
    s.visibility([p,s.labels[-2]],0,.5)
    s.visibility([q,s.labels[-1]],.5,1)
    s.move(q,T(1,-6,0),T(11,-6,-4),(.51,.68),(0,0,1),'cofilin dephosphorylation symbol')
    s.text('cofilin：Ser3 调控状态示意',(0,-12,-5.8),.9)

def sever(s):
    s.scene_scale=78
    record=json.loads((ROOT/'assets/cofilin-path-repairs.json').read_text())
    motion_file=ROOT/'assets/cofilin_boundary/motion.npz'
    if not record['passed'] or record['source']['entry']!='6UC4':
        raise ValueError('Cofilactin boundary preparation missing or failed; prepare this project again')
    if hashlib.sha256(motion_file.read_bytes()).hexdigest()!=record['motion_sha256']:
        raise ValueError('Cofilactin motion file changed; prepare this project again')
    for unit in record['geometry']['units']:
        if hashlib.sha256((ROOT/'assets/cofilin_boundary'/(unit['unit']+'.npz')).read_bytes()).hexdigest()!=unit['sha256']:
            raise ValueError('Cofilactin source asset changed: '+unit['unit'])
    with np.load(motion_file) as data:
        poses=data['matrices'];units=data['units'].tolist()
    obs=[]
    for i,unit in enumerate(units):
        chains=record['units'][unit];cofilin=unit in record['placements']
        obj=s.asset('cofilin_boundary',unit,name='actin_001' if unit=='FG' else ('cofilin ' if cofilin else 'actin ')+unit,
            color='gold' if unit=='FG' else 'pink' if cofilin else 'teal')
        obj['source_structure']='6UC4';obj['source_chain']='+'.join(chains)
        obj['geometry_evidence']='Complete modeled biological subunit; complementary source domains joined before meshing'
        if unit=='FG':obj['biological_id']='actin_001'
        if cofilin:obj['placement_evidence']='Bounded whole-cofilin rigid fit in a 9.2 A boundary model; not an atomic interface measurement'
        obj.rotation_mode='QUATERNION'
        for f in range(1,s.n+1):
            index=(f-1)/(s.n-1)*(len(poses)-1);lo=int(index);hi=min(lo+1,len(poses)-1);q=index-lo
            a,qa,_=Matrix(poses[lo,i]).decompose();b,qb,_=Matrix(poses[hi,i]).decompose()
            obj.location=a.lerp(b,q);obj.rotation_quaternion=qa.slerp(qb,q)
            obj.keyframe_insert('location',frame=f);obj.keyframe_insert('rotation_quaternion',frame=f)
        obs.append(obj)
    s.events.append({'object':'Complete cofilactin boundary and tracked monomer','first':1,'last':s.n,
        'action':'Local cluster binding, rigid boundary separation, then terminal monomer release',
        'limits':record['limits'],'motion_preparation':str(ROOT/'assets/cofilin-path-repairs.json')})
    s.center(obs,offset=(0,0,0),rotation=R(-12,'Y'))
    s.text('局部 cofilin 簇',(-20,-40,-13.8),1.6)
    s.text('未结合区',(23,-40,-13.8),1.6)
    label=s.text('完整单体返回单体池',(2,-40,-9.5),1.6)
    s.visibility([label],.84,1)

def debranch(s):
    s.scene_scale=96
    groups=branch_setup(s)
    child,childobs=groups['daughter']
    s.docking(child,Matrix.Identity(4),(.27,.49),s.plans['debranch']['offset'],release=True,repair=s.repairs['motions']['debranch'])
    gm=s.asset('4JD2','H',color='gold',main=False)
    center=np.load(ROOT/'assets/4JD2/H.npz')['atoms'].mean(0)
    s.move(gm,T(-23,-16,4)@T(*(-center)),T(-18,-16,4)@T(*(-center)),(.08,.23),(0,-3,1),'GMF stand-off symbol; not a docking trajectory')
    allobs=sum([v[1] for v in groups.values()],[])
    s.center(allobs,offset=(-2,0,0),rotation=R(-20))
    s.text('GMF',(-17,-16,7),1.8)
    s.text('母丝',(-13,-16,-8),1.8)
    s.text('接头保留 / 再成核',(16,-16,-12),1.6)

def cap_control(s):
    s.scene_scale=84
    core,fr,heads=formin_setup(s)
    cap,capobs=s.oldcap('formin')
    s.docking(fr,Matrix.Identity(4),(.12,.32),s.plans['formin_release']['offset'],release=True,far_shift=s.repairs['motions']['formin_release']['far_shift'],repair=s.repairs['motions']['formin_release'])
    s.docking(cap,cap.matrix_basis.copy(),(.42,.65),s.plans['cap']['offset'],repair=s.repairs['motions']['cap_formin'])
    s.text('formin 离开 → CapZ 封端',(0,-16,8),1.8)

def recycle(s):
    s.scene_scale=60
    rig,obs=s.compound('6FM2-assembly1',['A','B','A-2','B-2'],['gold','pink','teal','pink'])
    obs[0]['biological_id']='actin_001'
    s.docking(obs[0],Matrix.Identity(4),(.08,.28),Vector(s.plans['cap_recycle']['offset'])*.6,repair=s.repairs['motions']['cap_recycle'])
    s.center(obs,offset=(-2,0,0),rotation=R(35,'Y')@R(15,'Z'))
    old,_=s.token('ADP',(-2,-10,-3),'blue',1.15)
    new,new_label=s.token('ATP',(15,-10,-5),'gold',1.15)
    new_label.data.materials[0]=s.mats['dim']
    s.move(old,T(-2,-10,-3),T(-15,-10,-6),(.36,.49),(0,0,-1),'old ADP leaves pool symbol')
    s.move(new,T(15,-10,-5),T(-2,-10,-3),(.57,.72),(0,0,1),'different ATP molecule enters symbol')
    s.text('CAP–actin',(0,-16,5),1.3)
    co=s.asset('6VAO','I',color='pink')
    center=np.load(ROOT/'assets/6VAO/I.npz')['atoms'].mean(0)
    s.move(co,T(-12,-12,-1)@T(*(-center)),T(-18,-12,2)@T(*(-center)),(.06,.24),(0,-3,1),'cofilin partner release shown at schematic stand-off positions')
    s.text('cofilin',(-13,-12,7),1.4)

def closing(s):
    s.scene_scale=52
    labels=['输入 / 许可','原料 / 成核','延伸 / 功能','终止 / 拆解','回收 / 再进入']
    loc=[(-16,5),(0,5),(16,5),(9,-4),(-9,-4)]
    for label,(x,z) in zip(labels,loc):s.node(label,x,z,'teal',13)
    for (x,z),(a,b) in zip(loc,loc[1:]+loc[:1]):
        direction=Vector((a-x,0,b-z)).normalized()
        edge=min(6.7/max(abs(direction.x),1e-8),1.3/max(abs(direction.z),1e-8))
        s.arrow(Vector((x,0,z))+direction*edge,Vector((a,0,b))-direction*edge,'gold')
    s.text('同一亚基的身份持续；核苷酸可以交换',(0,-10,-8),.85)

def main():
    p=argparse.ArgumentParser();p.add_argument('--shots',default='all');p.add_argument('--width',type=int,default=OPTIONS['width'])
    p.add_argument('--preview',action='store_true');a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    config=json.loads((ROOT/'storyboard.json').read_text());selected=a.shots.split(',')
    if (selected==['all'] or set(selected)&{'a05','a07','c05','f06'}) and not (ROOT/'assets/contact-repairs.json').exists():
        raise ValueError('Contact repair assets missing; run build_actin_films.py --steps prepare for this project first')
    if (selected==['all'] or 'f03' in selected) and not (ROOT/'assets/fh2-seed-repairs.json').exists():
        raise ValueError('FH2 seed assets missing; run build_actin_films.py --steps prepare for this project first')
    if (selected==['all'] or set(selected)&{'f04','f07'}) and not (ROOT/'assets/fh2-growth.json').exists():
        raise ValueError('FH2 growth assets missing; run build_actin_films.py --steps prepare for this project first')
    if (selected==['all'] or 'f02' in selected) and not (ROOT/'assets/actin-path-repairs.json').exists():
        raise ValueError('Regulatory path preparation missing; run build_actin_films.py --steps prepare for this project first')
    if (selected==['all'] or 'a03' in selected) and not (ROOT/'assets/branch-path-repairs.json').exists():
        raise ValueError('Branch path preparation missing; run build_actin_films.py --steps prepare for this project first')
    if (selected==['all'] or set(selected)&{'a04','a08'}) and not (ROOT/'assets/branch-growth.json').exists():
        raise ValueError('Branch growth preparation missing; run build_actin_films.py --steps prepare for this project first')
    if (selected==['all'] or 'c04' in selected) and not (ROOT/'assets/cofilin-path-repairs.json').exists():
        raise ValueError('Cofilin path preparation missing; run build_actin_films.py --steps prepare for this project first')
    methods={'a01':lambda s:input_scene(s),'f01':lambda s:input_scene(s,True),'a02':wrc_scene,
             'c01':supply_scene,'a03':nucleate_branch,'a04':branch_growth,'a08':lambda s:branch_growth(s,True),
             'f02':autoinhib,'f03':formin_seed,'f04':formin_growth,'f07':lambda s:formin_growth(s,True),
             'c02':chemistry,'a05':network,'f05':lambda s:network(s,True),'a06':stop_signal,
             'c03':cofilin_control,'c04':sever,'a07':debranch,'f06':cap_control,'c05':recycle,'c06':closing}
    for c in config['shots']:
        if selected!=['all'] and c['id'] not in selected:continue
        s=Shot(c,a.width);methods[c['id']](s);folder=s.finish()
        if a.preview:
            for fraction in [.05,.48,.85]:
                bpy.context.scene.frame_set(s.frame(fraction))
                bpy.context.scene.render.filepath=str(folder/('key_'+str(round(fraction*100))+'.png'))
                bpy.ops.render.render(write_still=True)
        print('CHAPTER_BUILT',c['id'],c['frames'],flush=True)

if __name__=='__main__':main()
