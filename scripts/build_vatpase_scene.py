"""Build a 16-second seamless V-ATPase loop, with a whole-pump and c-ring view.

Pump envelopes come from 6VQ6. Rotation, solvent-pool tokens, membrane and
half-channel guides are authored mechanism graphics, not atomistic trajectories.
"""
import argparse
import json
import math
from pathlib import Path
import sys
from types import SimpleNamespace

import bpy
import numpy as np
from mathutils import Vector

SKILL=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(SKILL/'scripts'))
from detail_scene import Parts, collection, material, mesh, path, studio

parser=argparse.ArgumentParser()
parser.add_argument('--asset-dir',type=Path,required=True)
parser.add_argument('--crowder-dir',type=Path,required=True)
parser.add_argument('--font-zh',type=Path,required=True)
parser.add_argument('--font-latin',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--compositor',choices=['clean','reference','monochrome'],default='reference')
parser.add_argument('--width',type=int,default=1600)
parser.add_argument('--samples',type=int,default=64)
parser.add_argument('--engine',choices=['eevee','cycles'],default='eevee')
parser.add_argument('--preview',action='store_true')
args=parser.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
ROOT=args.output.resolve();ASSET=args.asset_dir.resolve()
if ROOT.exists():parser.error('Choose a fresh project output directory.')
if ROOT.is_relative_to(SKILL):parser.error('Project output must be outside the skill directory.')
ROOT.mkdir(parents=True)
if args.width<480 or args.width%32:parser.error('Use width >=480 and divisible by 32 for even 16:9 dimensions.')
bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
bpy.data.orphans_purge(do_recursive=True)
groups={n:collection(n) for n in ['Pump','Membrane','C ring detail','Flow','Labels','Environment']}
palette={'A':'53B6AF','B':'78A9BC','rotor':'F1C168','ring':'D1A35D','ring_tracer':'B87962','a':'B97799',
         'stalk':'B3C9CD','lipid':'345362','tail':'466169','h':'F9B557','atp':'BCA7DA',
         'in':'6CD5D0','out':'E9A56D','background':'102932'}
mats={k:material(k,c,.5,0) for k,c in palette.items()}
def emission(name,color,strength=1):
    m=material(name,color,.5,0);p=m.node_tree.nodes.get('Principled BSDF')
    p.inputs['Emission Color'].default_value=m.diffuse_color
    p.inputs['Emission Strength'].default_value=strength
    return m
ink=emission('Readable large labels','DBE9E9',.8)
mats['h']=emission('Proton markers','FFC66B',.4)
font=bpy.data.fonts.load(str(args.font_zh.resolve()));font.pack()
latin=bpy.data.fonts.load(str(args.font_latin.resolve()));latin.pack()

def empty(name,group='Pump',parent=None,location=(0,0,0)):
    o=bpy.data.objects.new(name,None);groups[group].objects.link(o)
    o.parent=parent;o.location=location;return o
def label(body,x,z,size=.92,align='CENTER'):
    data=bpy.data.curves.new(body,'FONT');data.body=body;data.font=font;data.font_bold=latin
    for c,f in zip(body,data.body_format):f.use_bold=ord(c)<0x400 or 0x2010<=ord(c)<=0x22ff
    data.size=size;data.align_x=align;data.align_y='CENTER'
    o=bpy.data.objects.new(body,data);groups['Labels'].objects.link(o)
    o.location=(x,-17,z);o.rotation_euler=(math.pi/2,0,0);data.materials.append(ink)
    return o
def ball(name,r,mat,group='Flow',parent=None,location=(0,0,0)):
    p=Parts();p.ball((0,0,0),r,seg=20,rings=12)
    o=p.object(name,[mat],groups[group]);o.parent=parent;o.location=location;return o
def curve(name,pts,r,mat,group='Flow'):
    return path(name,pts,r,mat,groups[group])
def bezier(points,u):
    p=np.array(points);v=1-u
    return v*v*v*p[0]+3*v*v*u*p[1]+3*v*u*u*p[2]+u*u*u*p[3]
def smooth(u):
    u=np.clip(u,0,1);return u*u*(3-2*u)
def key(o,f,location=None,scale=None,rotation=None):
    if location is not None:o.location=location;o.keyframe_insert('location',frame=f)
    if scale is not None:o.scale=(scale,)*3;o.keyframe_insert('scale',frame=f)
    if rotation is not None:o.rotation_euler=rotation;o.keyframe_insert('rotation_euler',frame=f)

info=json.loads((ASSET/'manifest.json').read_text())
main=empty('Whole pump',location=(-11,0,0))
rotor=empty('Central rotor clockwise from cytosol',parent=main)
detail=empty('Isolated complete c ring','C ring detail',location=(12.8,0,7.6))
detail.rotation_euler=(math.pi/2,0,0);detail.scale=(1.6,)*3
detail_rotor=empty('c ring detail rotor','C ring detail',parent=detail)
primary=[];asset_lookup={};ring_assets=[]
for item in info['chains']:
    c=item['chain'];d=np.load(ASSET/item['file'])
    color='A' if c in 'ABC' else 'B' if c in 'DEF' else 'a' if c=='a' else 'ring_tracer' if c=='b' else 'ring' if c in info['c_ring_chains'] else 'rotor' if c in info['rotor_chains'] else 'stalk'
    mat=mats[color]
    if c in 'ABC':mat=mat.copy();mat.name='Catalytic A '+c
    o=mesh('6VQ6 chain '+c,d['vertices'].tolist(),d['faces'].tolist(),[mat],groups['Pump'])
    o.parent=rotor if c in info['rotor_chains'] else main
    mod=o.modifiers.new('Sub-angstrom envelope smoothing','SMOOTH');mod.factor=.25;mod.iterations=2
    o['source']='6VQ6';o['source_chain']=c;o['asset_file']=str(ASSET/item['file'])
    o['role']=item['role'];primary.append(o);asset_lookup[c]=o
    if c in info['c_ring_chains']:
        copy=o.copy();copy.data=o.data;groups['C ring detail'].objects.link(copy)
        copy.name='Complete ring chain '+c;copy.parent=detail_rotor;ring_assets.append(copy)

# All lipids are complete head/tail units. Front units are omitted as a viewing
# window; there is no Boolean cut through a protein or lipid.
rng=np.random.default_rng(26)
lip=Parts()
for x in np.arange(-24.5,1.5,.82):
    for y in np.arange(-5.0,7.5,.82):
        local=np.array([x+11,y])
        if (local[0]/8.2)**2+(local[1]/5.8)**2<1:continue
        if y<-1.5 and -19<x<-2:continue
        for s in [-1,1]:
            xx=x+rng.uniform(-.15,.15);yy=y+rng.uniform(-.15,.15)
            z=s*2.05
            lip.ball((xx,yy,z),(.25,.25,.22),0,seg=10,rings=6)
            for dx in [-.11,.11]:
                lip.tube([(xx+dx,yy,z-s*.16),(xx+dx*.8,yy+.06,z-s*.75),
                          (xx+dx*1.5,yy-.05,z-s*1.95)],.065,1,sides=5)
lip.object('Whole lipid molecules with omitted foreground units',[mats['lipid'],mats['tail']],groups['Membrane'])

# A few source-derived background macromolecules provide depth, without implying
# a measured cytosolic composition or concentration.
bg=np.load(args.crowder_dir/'chain_A.npz')
for i,pos in enumerate([(-25,14,13),(-1,15,20),(-3,15,4)]):
    o=mesh('Background ubiquitin '+str(i),bg['vertices'].tolist(),bg['faces'].tolist(),[mats['background']],groups['Environment'])
    o.location=pos;o.rotation_euler=tuple(rng.uniform(0,math.tau,3))

# Right-hand view: the ten full c-ring chains, viewed from the cytoplasmic side.
# Subunit a is omitted wholesale in this inset, while separately drawn curves
# explain the two half-channels. Dots are displaced above/outside the surface.
center=np.array([12.8,-8.5,7.6]);radius=7.25
in_angle=math.radians(-12);out_angle=math.radians(12)
carry=math.tau-math.radians(24)
def point(angle):return center+np.array([radius*math.cos(angle),0,radius*math.sin(angle)])
inlet=point(in_angle);outlet=point(out_angle)
incoming=[np.array([32,-8.5,-3]),np.array([24,-8.5,-2]),np.array([23,-8.5,4]),inlet]
outgoing=[outlet,np.array([23,-8.5,11]),np.array([24,-8.5,16]),np.array([32,-8.5,18])]
curve('Cytoplasmic aqueous half-channel guide',[bezier(incoming,u) for u in np.linspace(.10,1,70)],.10,mats['in'])
curve('Luminal aqueous half-channel guide',[bezier(outgoing,u) for u in np.linspace(0,.90,70)],.10,mats['out'])
for name,pts,mat in [('in',incoming,mats['in']),('out',outgoing,mats['out'])]:
    u=.63;v=Vector(bezier(pts,u)-bezier(pts,u-.02)).normalized();p=Vector(bezier(pts,u))
    wing=v.cross(Vector((0,-1,0))).normalized()
    curve(name+' flow arrow',[p-v*.7+wing*.33,p,p-v*.7-wing*.33],.09,mat)
curve('Ring direction arc',[center+np.array([radius*math.cos(t),.65,radius*math.sin(t)]) for t in np.linspace(math.radians(160),math.radians(100),55)],.045,mats['stalk'])
p=Vector(center+np.array([radius*math.cos(math.radians(100)),.65,radius*math.sin(math.radians(100))]))
v=Vector((math.sin(math.radians(100)),0,-math.cos(math.radians(100))))
wing=v.cross(Vector((0,-1,0))).normalized()
curve('Clockwise arrow',[p-v*.65+wing*.3,p,p-v*.65-wing*.3],.075,mats['stalk'])
protons=[]
for i,c in enumerate(info['carriers']):
    o=ball('H transport token '+c['chain'],.25 if i else .34,mats['h'])
    o['carrier_chain']=c['chain'];o['carrier_residue']=c['residue']
    protons.append((o,math.atan2(c['site_nm'][1],c['site_nm'][0])))

# One repeated ATP graphic splits into an ADP graphic and an independently
# tracked terminal phosphate. Three turnovers accompany each complete rotation.
adp=empty('Adenosine diphosphate token','Flow')
token=Parts();token.ball((-.48,0,0),(.48,.23,.37),0,seg=16,rings=10)
token.ball((.2,0,0),.24,1,seg=16,rings=10);token.ball((.78,0,0),.24,1,seg=16,rings=10)
token.tube([(-.1,0,0),(.2,0,0),(.78,0,0)],.07,1)
token.object('Adenosine and two phosphates',[mats['atp'],mats['h']],groups['Flow']).parent=adp
gamma=ball('Terminal phosphate becomes Pi',.24,mats['h'])
link=curve('ATP terminal bond',[(.78,0,0),(1.36,0,0)],.07,mats['h']);link.parent=adp
dock=np.array([-15.0,-10,14.5])
atp_in=[[-32,-10,19],[-23,-10,19],[-18,-10,18],dock]
adp_out=[dock,[-18,-10,12],[-22,-10,10],[-29,-10,9]]
pi_out=[dock+np.array([1.36,0,0]),[-18,-10,10],[-24,-10,7],[-30,-10,6]]
curve('ATP approach guide',[bezier(atp_in,u) for u in np.linspace(.27,.94,40)],.028,mats['background'])

label('V 型 ATPase',-11,23.2,1.18)
label('c 环俯视',12.8,19.6,1.06)
label('ATP',-22.2,20,.98)
label('ADP',-22.3,12.3,.98)
label('Pi',-22.3,6.8,.98)
label('胞质',-22.5,3.7,1.0)
label('囊泡腔',-21.8,-4.6,1.0)
label('V1',-1.1,15.1,1.0)
label('V0',1.0,-.2,1.0)
label('胞质入口',21.2,-3.1,1.05)
label('囊泡出口',21.2,17.6,1.05)
label('H+',6.3,14.6,1.0)

profile=json.loads((SKILL/'assets/style-molecular-reference.json').read_text())
profile.update(name='Continuous V-ATPase loop; large labels only',background='09242E',world='8DA7B2',exposure=0)
profile['compositor']={'preset':'reference','contrast':1}
profile['lighting']=[
    {'name':'Large soft key','position':[-20,-23,30],'power':21000,'color':'F6E6D1','size':16},
    {'name':'Soft frontal fill','position':[21,-23,15],'power':18000,'color':'C3E5EF','size':18},
    {'name':'Cool rim','position':[1,10,18],'power':16000,'color':'B3D7EF','size':12},
    {'name':'Membrane underside','position':[-13,-18,-12],'power':9000,'color':'D7E5E8','size':14}]
settings=SimpleNamespace(engine=args.engine,width=args.width,fps=24,frames=384,samples=args.samples,compositor=args.compositor)
cam=studio(profile,settings,{'position':(0,-100,8.5),'target':(0,0,8.5),'scale':58,'type':'ORTHO'})
scene=bpy.context.scene;scene.render.resolution_y=round(args.width*9/16)
if args.engine=='eevee':scene.eevee.taa_render_samples=args.samples
scene.render.image_settings.color_mode='RGB';scene.render.film_transparent=True
scene.render.fps=24
events=[]
for f in range(1,386):
    t=(f-1)/384;angle=-math.tau*t
    key(rotor,f,rotation=(0,0,angle));key(detail_rotor,f,rotation=(0,0,angle))
    for o,base in protons:
        a=base+angle;phase=(in_angle-a)%math.tau
        if phase<=carry:
            pos=point(a);scale=1
        else:
            gap=(phase-carry)/(math.tau-carry)
            if gap<.45:
                pos=bezier(outgoing,gap/.45);scale=1
            elif gap>.55:
                pos=bezier(incoming,(gap-.55)/.45);scale=1
            else:
                pos=np.array([34,-8.5,20]);scale=0
        key(o,f,location=pos,scale=scale)
    phase=(3*t)%1
    if phase<.42:
        p=bezier(atp_in,smooth(phase/.42));g=p+np.array([1.36,0,0]);s=1;gs=1;ls=1
    elif phase<.52:
        p=dock;g=p+np.array([1.36,0,0]);s=1;gs=1;ls=1
    elif phase<.90:
        u=smooth((phase-.52)/.38);p=bezier(adp_out,u);g=bezier(pi_out,u)
        s=1;gs=1;ls=0
    else:
        p=np.array([-30,-10,9]);g=np.array([-30,-10,6]);s=0;gs=0;ls=0
    key(adp,f,location=p,scale=s);key(gamma,f,location=g,scale=gs);key(link,f,scale=ls)
    for i,c in enumerate('ABC'):
        shader=asset_lookup[c].data.materials[0].node_tree.nodes.get('Principled BSDF')
        u=(t*3-i)%3
        strength=.13*max(0,math.sin(math.pi*u))**2 if u<1 else 0
        shader.inputs['Emission Color'].default_value=asset_lookup[c].data.materials[0].diffuse_color
        shader.inputs['Emission Strength'].default_value=strength
        shader.inputs['Emission Strength'].keyframe_insert('default_value',frame=f)

# Constant-rate rendering between the baked states, including the loop seam.
for action in bpy.data.actions:
    for layer in action.layers:
        for strip in layer.strips:
            for slot in action.slots:
                bag=strip.channelbag(slot,ensure=False)
                if bag:
                    for fc in bag.fcurves:
                        for k in fc.keyframe_points:k.interpolation='LINEAR'
scene.frame_end=384;scene.frame_set(1)
scene['animation_mode']='continuous seamless loop'
scene['loop_frames']=384;scene['full_turn_seconds']=16
scene['chemistry']='3 ATP turnovers and 10 cyclic carrier trajectories per rotor revolution; graphical tokens, no atomistic chemical simulation'
scene['inset']='Whole c-ring chains only, seen from V1/cytosolic side. Complete a and other subunits omitted in inset; guides represent half-channels, not reconstructed pore surfaces.'
scene['display_policy']='Large essential labels only. Sources and modeling limits are in the delivery document.'
scene['provenance_manifest']=str(ASSET/'manifest.json')
from vatpase_motion import add_main_protons
add_main_protons(ASSET)
scene.render.filepath=str(ROOT/'preview.png')
bpy.ops.wm.save_as_mainfile(filepath=str(ROOT/'scene.blend'))
manifest={'source':'6VQ6','asset_dir':str(ASSET),'crowder_dir':str(args.crowder_dir.resolve()),'compositor':args.compositor,'frames':384,'fps':24,'seconds':16,'width':args.width,'height':scene.render.resolution_y,
          'engine':scene.render.engine,'samples':args.samples,'mode':'continuous-loop','font_zh':font.name,'font_latin':latin.name,'font_files':[str(args.font_zh.resolve()),str(args.font_latin.resolve())],
          'labels':[o.data.body for o in groups['Labels'].objects if o.type=='FONT'],
          'atomistic_trajectory':False,'left_proton_motion':True,'inset':scene['inset'],'display_policy':scene['display_policy'],
          'stoichiometry':{'ATP_per_revolution':3,'H_carriers':10},'blender':bpy.app.version_string}
(ROOT/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
if args.preview:
    for f in [53,145,278]:
        scene.frame_set(f);scene.render.filepath=str(ROOT/f'preview_{f:03d}.png')
        bpy.ops.render.render(write_still=True)
