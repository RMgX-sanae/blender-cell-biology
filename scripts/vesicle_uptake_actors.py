"""Choreograph whole molecular actors on the prepared membrane and its persistent cargo."""
import hashlib
import json
from pathlib import Path

import bpy
import numpy as np
from mathutils import Matrix, Vector

import detail_scene as ds
from apply_scale_environment import surface
from blender_clathrin_assembly import load_assets
from molecular_asset_io import chain_file
from scene_chapters import animated_path, balls, caption_set, empty, font_pair, visible


def frame(t,n):return 1+round((n-1)*t)


def ramp(t,a,b):
    x=np.clip((t-a)/(b-a),0,1)
    return x*x*(3-2*x)


def norm(v):return v/np.linalg.norm(v)


def keys(obj,poses):
    for f,point in poses:
        obj.location=point
        obj.keyframe_insert(data_path='location',frame=f)


def show(actor,first,last):
    for obj in [actor,*actor.children_recursive]:
        if obj.type in ('MESH','CURVE'):visible(obj,first,last)


def source_mesh(name,arrays,mat,group,parent):
    obj=ds.mesh(name,arrays[0],arrays[1],[mat],group)
    obj.parent=parent
    return obj


def add_coat(folder,rows,g,mats):
    spec,_=load_assets(folder);units=spec['instances']
    if spec['source_config']['source_id']!='3IYV':raise ValueError('Expected a qualified 3IYV cage')
    center=np.mean([u['hub_nm'] for u in units],axis=0);data={}
    for c in {c for u in units for c in u['chains']}:
        with np.load(chain_file(folder/'prototypes',c)) as z:data[c]=(z['vertices'],z['faces'])
    # The cage is centered at the origin of its source symmetry operators.
    # Choose one isotropic display scale so its innermost modeled envelope clears the bilayer.
    radii=np.concatenate([np.linalg.norm(v-center,axis=1) for v,_ in data.values()])
    scale=1.13/np.min(radii)
    if np.linalg.norm(center)>1e-5:raise ValueError('The source cage symmetry center changed')
    root=empty('Clathrin cap on cytosolic side',g)
    keys(root,[(r['frame'],np.array([*r['center'][:2],r['receptor'][2]+1])) for r in rows if r['time']<=.46])
    root.scale=(scale,)*3
    heavy={c for group in spec['source_config']['groups'] for c in group['heavy']}
    protos={};used=[]
    for u in sorted(units,key=lambda u:u['hub_nm'][2]):
        mat=np.asarray(u['matrix_nm'])
        max_z=max(np.max(data[c][0]@mat[2,:3]+mat[2,3]) for c in u['chains'])*scale
        # Reveal only entire triskelia that fit below the donor sheet.
        first=next((r['frame'] for r in rows if .14<r['time']<.336 and r['receptor'][2]+1+max_z<-.14),None)
        if first is None:continue
        outward=norm(np.array(u['hub_nm'])-center)
        rig=empty('Triskelion '+u['id'],g);rig.parent=root;rig.matrix_basis=Matrix(u['matrix_nm'])
        rig['source_unit']=u['id'];rig['source_chains']=','.join(u['chains'])
        for c in u['chains']:
            name='3IYV '+u['id']+' '+c
            if c in protos:obj=ds.duplicate(protos[c],name,g)
            else:obj=source_mesh(name,data[c],mats['coat'],g,rig);protos[c]=obj
            obj.parent=rig;obj.matrix_basis=Matrix.Identity(4);obj['source_chain']=c
            obj.material_slots[0].link='OBJECT'
            role='light' if c not in heavy else ('coat_focus' if outward[1]<-.61 else
                                                'coat_dim' if outward[1]>.20 else 'coat')
            obj.material_slots[0].material=mats[role]
        base=mat[:3,3].copy()
        release=.357+.011*(outward[1]+1)/2
        stop=frame(release+.050,len(rows))
        keys(rig,[(first,base),(frame(release,len(rows)),base),
                  (stop,base+outward*.57/scale+np.array([.07,-.09,-.05])/scale)])
        show(rig,first,stop);used.append(u['id'])
    if len(used)<3:raise ValueError('Too few intact coat units clear the donor sheet')
    root['source']='3IYV';root['scale_from_nm']=scale;root['shown_units']=json.dumps(used)
    return used


def add_ap2(folder,rows,g,mats):
    manifest=json.loads((folder/'manifest.json').read_text())
    if (not manifest.get('membrane_screen_passed') or manifest['source']['id']!='6YAF'
            or manifest['membrane']['normal']!=[0,0,1]):raise ValueError('Use the registered 6YAF AP2/cargo membrane assets')
    names=['anchor_P.npz']+['moving_'+c+'.npz' for c in 'ABMS']
    for name in names:
        if hashlib.sha256((folder/name).read_bytes()).hexdigest()!=manifest['data_sha256'][name]:
            raise ValueError('Changed AP2 source mesh '+name)
    with np.load(folder/names[0]) as z:anchor=z['atoms'].mean(0)
    group=empty('AP2 capture on receptor tail',g)
    group.rotation_euler=(np.pi,0,0);size=.059;group.scale=(size,)*3
    for name in names:
        with np.load(folder/name) as z:
            obj=source_mesh('6YAF '+name[:-4],(z['vertices']-anchor,z['faces']),mats['motif' if name==names[0] else 'adaptor'],g,group)
            obj['source']='6YAF'
    n=len(rows)
    def point(row):
        t=row['time'];approach=(1-ramp(t,.015,.085))*np.array([-.46,-.26,-.39])
        release=ramp(t,.365,.417)*np.array([-.63,-.15,-.45])
        return np.asarray(row['receptor'])+np.asarray(row['normal'])*.27+approach+release
    keys(group,[(r['frame'],point(r)) for r in rows if r['time']<.45])
    show(group,1,frame(.423,n))


def add_dynamin(rows,g,mat):
    n=len(rows);frames=list(range(frame(.23,n),frame(.35,n),3))+[frame(.35,n)]
    def helix(f):
        r=rows[f-1];theta=np.linspace(.22,4.8*np.pi,128)
        radius=r['neck_radius']+.154
        return np.column_stack((radius*np.cos(theta),radius*np.sin(theta),np.linspace(-.34,-.55,len(theta))))
    obj=animated_path('Dynamin neck collar',frames,helix,.051,mat,g)
    obj['scope']='Schematic GTP-linked helix on the cytosolic neck, no copy-number or structural trajectory claim'
    visible(obj,frames[0],rows[next(i for i,r in enumerate(rows) if r['time']>=.35)-1]['frame'])


def coiled_tether(a,b,length):
    # Sample a bow of fixed discrete contour length at changing endpoint separation.
    delta=b-a;direction=norm(delta)
    reference=np.array([0.,0.,1.]) if abs(direction[2])<.9 else np.array([0.,1.,0.])
    side=norm(reference-direction*np.dot(reference,direction))
    u=np.linspace(0,1,23)[:,None]
    def points(h):return a+u*delta+h*np.sin(np.pi*u)*side
    if np.linalg.norm(delta)>length+1e-6:raise ValueError('Tether endpoints exceed the retained contour')
    low,high=0.,length
    for _ in range(24):
        mid=(low+high)/2
        if np.linalg.norm(np.diff(points(mid),axis=0),axis=1).sum()<length:low=mid
        else:high=mid
    return points((low+high)/2)


def target_machinery(folder,rows,target,g,mats):
    record=json.loads((folder/'structures.json').read_text())
    if hashlib.sha256((folder/'structures.npz').read_bytes()).hexdigest()!=record['geometry_sha256']:
        raise ValueError('Changed Rab5 or SNARE geometry')
    n=len(rows);F=lambda t:frame(t,n)
    z=np.load(folder/'structures.npz',allow_pickle=False)
    side=norm(np.array([.50,.07,.86]))
    lateral=norm(np.cross([0,0,1],side));up=np.cross(side,lateral)
    axes=np.column_stack((side,lateral,up));scale=.078
    def rab_origin(r):
        return np.asarray(r['center'])+1.15*side+(1-ramp(r['time'],.426,.462))*np.array([-.26,-.53,.30])
    rab=empty('Carrier Rab5 coordinate frame',g)
    rab.rotation_euler=Matrix(axes.tolist()).to_euler();rab.scale=(scale,)*3
    keys(rab,[(r['frame'],rab_origin(r)) for r in rows if r['time']>.419])
    core=source_mesh('3MJH Rab5A', (z['A_verts'],z['A_faces']),mats['rab'],g,rab)
    finger=source_mesh('3MJH EEA1 zinc finger',(z['B_verts'],z['B_faces']),mats['tether'],g,rab)
    zinc=balls('3MJH zinc',z['ZN'],.16,mats['ion'],g,rab)
    show(core,F(.426),n);show(finger,F(.62),F(.862));show(zinc,F(.62),F(.862))
    lig=balls('3MJH GTP from cytosol',z['GTP'],.155,mats['gtp'],g,rab)
    mg=balls('3MJH Mg',z['MG'],.15,mats['ion'],g,rab)
    keys(lig,[(F(.477),(1.1,-3.7,2.8)),(F(.508),(0,0,0))])
    show(lig,F(.477),n);show(mg,F(.508),n)
    site=z['GTP'].mean(0)
    gdp=balls('Outgoing GDP identity probe',[site],.53,mats['gdp'],g,rab)
    keys(gdp,[(F(.425),(0,0,0)),(F(.452),(0,0,0)),(F(.479),(-1.,-3.5,1.8))]);show(gdp,F(.426),F(.479))
    # EEA1 C-terminal membrane anchor and the incoming Rab5-bound N-terminal finger.
    anchor=target+2.13*norm(np.array([-.72,.10,.69]))
    end_local=np.array(z['finger_tail'])
    def far(f):
        t=rows[f-1]['time']
        bound=rab_origin(rows[f-1])+axes@(end_local*scale)
        captured=rab_origin(rows[F(.62)-1])+axes@(end_local*scale)
        if t<.62:return captured+(1-ramp(t,.575,.62))*np.array([.13,.34,.16])
        return bound+ramp(t,.862,.936)*(anchor+np.array([-.30,.03,.89])-bound)
    frames=sorted({F(t) for t in (.575,.62,.69,.77,.80,.862,.936,1)}|set(range(F(.575),n+1,4)))
    length=max(np.linalg.norm(far(f)-anchor) for f in frames)+.15
    tether=animated_path('EEA1 target-anchored contour',frames,lambda f:coiled_tether(anchor,far(f),length),.024,mats['tether'],g)
    tether['contour']=length;tether['anchored_end']=anchor.tolist()
    visible(tether,frames[0],n)
    anchor_obj=balls('Target PI3P site',[anchor],.055,mats['tether'],g);visible(anchor_obj,frames[0],n)
    # Source coordinates of the parallel four-helix core. Missing membrane attachments are drawn separately.
    core_axes=np.array([[0,1,0],[0,0,-1],[-1,0,0]],float)
    helix_scale=.072
    def sites(row,code):
        t=row['time'];cx=row['center'][0]
        radius=.62
        left=cx+np.sqrt(1-radius**2);right=target[0]-np.sqrt(4-radius**2)
        lane='ABCD'.index(code)*.035
        # After a pore forms, the R anchor can move on its continuous rim
        # onto the same membrane patch as the other three anchors.
        x=(left+(right+.035-left)*ramp(t,.812,.908)) if code=='A' else right+.035*lane
        side=np.array([0,-np.cos(.74+lane),np.sin(.74+lane)])
        rad=max(np.sqrt(max(0,1-(x-cx)**2)),np.sqrt(max(0,4-(x-target[0])**2)),
                row['pore_radius'] if row['pore_connected'] else 0.)
        return np.array([x,0.,target[2]])+rad*side
    for c in 'ABCD':
        trace=z['helix_'+c];native=(trace-trace[0])*helix_scale
        segment=ds.path('2NPS '+c+' '+record['sources']['2NPS']['chains'][c],native,.027,mats[c],g)
        segment['source']='2NPS';segment['snare_chain']=c
        base_offset=trace[0]*helix_scale
        def pose(f,code=c,offset=base_offset):
            r=rows[f-1];t=r['time']
            rim_c=sites(r,'A');rim_t=sites(r,'B')
            weight=.5;spread=1-ramp(t,.717,.796)+ramp(t,.938,.994)
            sign=-1 if code=='A' else 1
            pick=sites(r,code)
            base=(rim_c+rim_t)*weight+np.array([0,-.045,.75])+core_axes@offset
            base+=spread*.72*(pick-(rim_c+rim_t)*weight)
            # Swing the C-terminal ends apart about the source N termini.
            rotation=Matrix.Rotation(float(sign*spread*.43),3,'Y')
            orient=rotation@Matrix(core_axes.tolist())
            return base,orient,pick
        fs=sorted({F(t) for t in (.69,.717,.77,.796,.815,.889,.938,.994,1)}|set(range(F(.69),n+1,5)))
        endpoint=native[-1]
        for f in fs:
            base,orient,_=pose(f);segment.location=base;segment.rotation_euler=orient.to_euler()
            segment.keyframe_insert(data_path='location',frame=f);segment.keyframe_insert(data_path='rotation_euler',frame=f)
        visible(segment,fs[0],n)
        def stub(f):
            base,orient,pick=pose(f)
            a=base+np.asarray(orient)@endpoint
            return np.array([a,.5*(a+pick)+np.array([0,-.04,.02]),pick])
        link=animated_path('Membrane attachment '+c,fs,stub,.014,mats[c],g)
        link['scope']='schematic continuation of the source motif to its assigned membrane'
        visible(link,fs[0],n)
    z.close()
    return length


def animate(scene,record,opt):
    rows=record['samples'];n=len(rows);F=lambda t:frame(t,n)
    g=ds.collection('Uptake and early-endosomal machinery')
    colors={'coat':'A57C5A','coat_focus':'F0B978','coat_dim':'756757','light':'CEB89E','adaptor':'D69599','motif':'E7F3A5','collar':'AB8EE1',
            'rab':'58BFA6','tether':'F4C779','gtp':'F6DF78','gdp':'A6ABCC','ion':'E8F7EE',
            'A':'ED9A80','B':'83BADF','C':'95D19E','D':'D6A3D9'}
    mats={k:surface('Uptake '+k,'#'+v,bump=0,emission=.33) for k,v in colors.items()}
    shown=add_coat(opt.coat,rows,g,mats)
    add_ap2(opt.ap2,rows,g,mats);add_dynamin(rows,g,mats['collar'])
    length=target_machinery(opt.machinery,rows,np.array(record['target_center']),g,mats)
    stages=[(0,'配体结合受体，AP2 招募被覆','Ligand capture recruits AP2 and clathrin'),
            (.15,'胞质侧的网格蛋白随膜内陷组装','Clathrin assembles on the cytosolic face'),
            (.25,'Dynamin 利用 GTP 收紧膜颈','Dynamin uses GTP during neck scission'),
            (.36,'辅伴侣招募 Hsc70，ATP 驱动去被覆','Auxilin or GAK recruits Hsc70 for ATP-driven uncoating'),
            (.435,'去被覆后，GEF 促使 Rab5 释放 GDP','After uncoating, a GEF promotes GDP release from Rab5'),
            (.479,'Rab5 结合来自胞质的 GTP','Rab5 binds GTP from the cytosol'),
            (.575,'内体表面的 EEA1 捕获囊泡','Endosomal EEA1 captures the carrier'),
            (.675,'系留蛋白弯曲，两层膜靠近','The tether bends as the membranes approach'),
            (.731,'内体 SNARE 从 N 端向 C 端配对','Endosomal SNAREs zipper from the N to the C ends'),
            (.80,'融合孔打开，两个腔相通','The fusion pore opens between the two lumens'),
            (.892,'同一批货物通过融合孔进入内体','The same cargo moves through the pore into the endosome'),
            (.947,'NSF 与 α-SNAP 利用 ATP 拆开 SNARE','NSF and α-SNAP use ATP to disassemble the cis complex')]
    cam=scene.camera
    views=[(0,[0,0,-.05],[.37,-1,.12],4.9),(.12,[0,0,-.35],[.35,-1,.11],5.5),
           (.265,[0,0,-1.43],[.28,-1,.012],7.5),(.43,[0,0,-1.60],[.23,-1,.08],6.3),
           (.51,[-.65,0,-3.25],[.18,-1,.28],6.8),(.605,[.55,0,-5.75],[.09,-1,.48],8.2),
           (.70,[.85,0,-5.50],[.06,-1,.49],6.6),(.759,[1.21,-.15,-5.24],[.08,-1,.40],3.4),
           (.795,[1.21,-.15,-5.24],[.08,-1,.40],3.4),(.836,[1.89,0,-5.95],[.01,-1,.76],8.6),
           (.942,[1.89,0,-5.95],[.01,-1,.76],8.6),(1,[1.35,-.15,-5.31],[.08,-1,.45],4.9)]
    for t,focus,front,size in views:
        cam.location=np.array(focus)+norm(np.array(front))*16;ds.aim(cam,focus)
        cam.data.ortho_scale=size
        for prop in ('location','rotation_euler'):cam.keyframe_insert(data_path=prop,frame=F(t))
        cam.data.keyframe_insert(data_path='ortho_scale',frame=F(t))
    fonts=font_pair(opt,opt.language)
    caption_set(scene,cam,fonts,opt.language,[(F(t),cn,en) for t,cn,en in stages],[(F(t),s) for t,_,_,s in views],outline=True)
    scene['uptake_machinery']=json.dumps({'coat_units':shown,'target_tether_contour':length,'source_roles':['6YAF','3IYV','3MJH','2NPS'],
        'fates':{'coat':'cytosolic pool','Rab5-GTP':'common endosomal membrane','EEA1':'target-side anchor','SNAREs':'cis complex then separated motifs on common membrane','cargo':'target lumen'}})


def audit_machinery(scene,record):
    info=json.loads(scene['uptake_machinery']);n=len(record['samples']);F=lambda t:frame(t,n)
    scene.frame_set(F(.341))
    units=[obj for obj in scene.objects if obj.get('source_unit')]
    if {o['source_unit'] for o in units}!=set(info['coat_units']) or not all(len(o.children)==6 and all(not ch.hide_render for ch in o.children) for o in units):
        raise ValueError('Whole six-chain clathrin display changed')
    scene.frame_set(F(.51))
    if any(not c.hide_render for u in units for c in u.children) or scene.objects['3MJH GTP from cytosol'].hide_render:
        raise ValueError('Coat did not return to the pool before the Rab5-GTP state')
    tether=scene.objects['EEA1 target-anchored contour']
    lengths=[];spans=[]
    for key in list(tether.data.shape_keys.key_blocks)[::4]:
        p=np.array([list(k.co)[:3] for k in key.data]);lengths.append(float(np.linalg.norm(np.diff(p,axis=0),axis=1).sum()));spans.append(float(np.linalg.norm(p[-1]-p[0])))
        if not np.allclose(p[0],tether['anchored_end'],atol=2e-6):raise ValueError('EEA1 target-side anchor moved')
    if max(abs(l-info['target_tether_contour']) for l in lengths)>.00001 or max(spans)-min(spans)<.3:
        raise ValueError('Tether contour or changing end-to-end separation differs')
    attachment_error=0.
    target=np.array(record['target_center'])
    for t in (.77,.889,1.):
        f=F(t);r=record['samples'][f-1];c=np.array(r['center']);axis=target-c
        for code in 'ABCD':
            block=scene.objects['Membrane attachment '+code].data.shape_keys.key_blocks['Frame '+str(f)]
            point=np.array(block.data[-1].co[:3])
            offset=np.clip(np.dot(point-c,axis)/np.dot(axis,axis),0,1)*axis
            phi=min(-point[2],np.linalg.norm(point-c)-r['radius'],np.linalg.norm(point-target)-record['target_radius'])
            if r['pore_connected']:phi=min(phi,np.linalg.norm(point-c-offset)-r['pore_radius'])
            attachment_error=max(attachment_error,abs(float(phi)))
    if attachment_error>2e-5:raise ValueError('A SNARE anchor no longer meets the actual membrane mid-surface')
    helices=[o for o in scene.objects if o.get('snare_chain')]
    if sorted(o['snare_chain'] for o in helices)!=list('ABCD'):raise ValueError('Missing complete four-chain endosomal SNARE motif set')
    scene.frame_set(F(.819));bound={o.name:o.matrix_world.copy() for o in helices}
    scene.frame_set(n)
    if any(o.hide_render or abs(o.matrix_world.to_quaternion().rotation_difference(bound[o.name].to_quaternion()).angle)<.1 for o in helices):
        raise ValueError('Native cis-SNARE disassembly motion is absent')
    return {**info,'validated_whole_triskelia':len(units),'tether_sample_contours':[min(lengths),max(lengths)],
            'tether_end_to_end_range':[min(spans),max(spans)],'source_SNARE_chains':4,'max_SNARE_anchor_midsurface_error':attachment_error}
