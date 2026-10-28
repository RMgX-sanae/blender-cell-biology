"""Build editable reference-style comparisons from a prepared molecular asset.

Molecular coordinates retain their nm scale. Cell appearance examples are
schematic adherent-cell arrangements in arbitrary display units, not reconstructions.
Run with Blender 5.2: --python build_style_scene.py -- --styles all --asset-dir ...
"""
import argparse
import json
import math
from pathlib import Path
import sys
import time

import bpy
import numpy as np
from mathutils import Vector

sys.path.insert(0, str(Path(__file__).resolve().parent))
import detail_scene as ds
import style_nodes as sn


def smooth(points, samples=8):
    points = np.asarray(points)
    if len(points) < 2:
        return points
    padded = np.concatenate((points[:1], points, points[-1:]))
    result = []
    for i in range(1, len(padded)-2):
        a, b, c, d = padded[i-1:i+3]
        for t in np.linspace(0, 1, samples, endpoint=False):
            result.append(.5*((2*b)+(-a+c)*t+(2*a-5*b+4*c-d)*t*t+(-a+3*b-3*c+d)*t**3))
    return np.array(result+[points[-1]])


def curves(name, lines, mat, radius, group):
    data = bpy.data.curves.new(name, 'CURVE')
    data.dimensions = '3D'
    for points in lines:
        if len(points) < 2:
            continue
        s = data.splines.new('POLY')
        s.points.add(len(points)-1)
        s.points.foreach_set('co', np.column_stack((points, np.ones(len(points)))).ravel())
    obj = bpy.data.objects.new(name, data)
    group.objects.link(obj)
    sn.curve_tubes(obj, mat, radius)
    return obj


def ribbon(data, secondary, chain, profile, group, white=False):
    atoms, names, numbers = data['atoms'], data['atom_names'], data['residue_numbers']
    ca_ids = np.flatnonzero(names == 'CA')
    if len(ca_ids) < 3:
        raise ValueError('Ribbon example needs a protein chain with at least three CA atoms.')
    ca, residues = atoms[ca_ids], numbers[ca_ids]
    if np.any(np.diff(residues)!=1) or np.any(np.linalg.norm(np.diff(ca,axis=0),axis=1)>.5):
        raise ValueError('CA trace has a chain break or nonconsecutive author numbering; prepare explicit fragments before drawing ribbons.')
    points = smooth(ca)
    # The author annotations decide helix/sheet/loop, not a geometric guess.
    mats = [sn.material('Ribbon '+k, profile, '#e9f0f5' if white else c, 'pbr')
            for k, c in zip(('sheet', 'helix', 'loop'), (profile['palette'][1], profile['palette'][2], profile['palette'][0]))]
    loop = curves('Continuous CA trace', [points], mats[2], .045, group)
    loop['evidence'] = 'smoothed experimental CA trace; not a backbone atom bond representation'
    # Peptide carbonyl directions establish ribbon width, with sign continuity.
    directions = []
    for index, res in enumerate(residues):
        oi = np.flatnonzero((numbers == res) & (names == 'O'))
        ci = np.flatnonzero((numbers == res) & (names == 'C'))
        v = atoms[oi[0]]-atoms[ci[0]] if len(oi) and len(ci) else np.array([0., 0., 1.])
        if directions and np.dot(v, directions[-1]) < 0:
            v = -v
        directions.append(v)
    directions = smooth(directions)
    for seg in secondary:
        if seg['chain'] != chain:
            continue
        ids = np.flatnonzero((residues >= seg['start']) & (residues <= seg['end']))
        if len(ids) < 2:
            continue
        start, end = ids[0]*8, ids[-1]*8
        verts, faces = [], []
        for i in range(start, end+1):
            tangent = points[min(i+1, len(points)-1)]-points[max(0, i-1)]
            tangent /= np.linalg.norm(tangent)
            width = directions[i]-tangent*np.dot(tangent, directions[i])
            width /= np.linalg.norm(width)
            normal = np.cross(tangent, width)
            fraction = (i-start)/(end-start)
            w = .18
            if seg['kind'] == 'sheet':
                w = .16 if fraction < .75 else max(.012, .29*(1-fraction)/.25)
            for sign_w, sign_n in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
                verts.append(points[i]+sign_w*w*width+sign_n*.018*normal)
        for i in range(end-start):
            for j in range(4):
                faces.append((4*i+j, 4*i+(j+1)%4, 4*(i+1)+(j+1)%4, 4*(i+1)+j))
        faces += [(3, 2, 1, 0), tuple(4*(end-start)+j for j in range(4))]
        obj = ds.mesh(f"{seg['kind']} {seg['start']}-{seg['end']}", verts, faces,
                      [mats[seg['kind'] == 'helix']], group)
        obj['evidence'] = 'deposited secondary-structure annotation, smoothed CA trace and peptide carbonyl orientation'


def molecule(asset_dir, chain, profile, group):
    data = np.load(asset_dir / f'chain_{chain}.npz', allow_pickle=False)
    rep = profile['representation']
    mat = sn.material('Main molecular material', profile)
    if rep in {'surface', 'crowded_surface', 'shell_ribbon'}:
        obj = ds.mesh('Coordinate-derived envelope', data['vertices'], data['faces'], [mat], group)
        obj['evidence'] = 'Gaussian envelope from deposited heavy atoms; not experimental density'
        mod = obj.modifiers.new('Subtle envelope smoothing', 'SMOOTH')
        mod.factor, mod.iterations = .25, 2
    if rep in {'atoms', 'bonds', 'environment_atoms'}:
        elements = data['elements']
        colors = {'C': profile['palette'][0], 'N': '#497acb', 'O': '#c55a51', 'S': '#d5b557', 'P': '#b694cc'}
        vdw = {'C': .17, 'N': .155, 'O': .152, 'S': .18, 'P': .18}
        for element in np.unique(elements):
            select = elements == element
            c = colors.get(element, profile['palette'][0])
            if profile['material'] == 'toon' or rep == 'environment_atoms':
                c = profile['palette'][0] if element == 'C' else profile['palette'][1]
            m = sn.material('Atom '+element, profile, c)
            obj = sn.point_object('Atoms '+element, data['atoms'][select], group)
            obj['element'] = str(element)
            sn.instance_spheres(obj, m, .055 if rep == 'bonds' else vdw.get(element, .17), 3)
            obj['evidence'] = 'deposited heavy atoms; conventional radii; hydrogen omitted'
        if rep == 'bonds':
            xyz = data['atoms']
            # Candidate protein heavy-atom bonds, restricted to same/adjacent residues.
            cov = np.array([{'C': .076, 'N': .071, 'O': .066, 'S': .105}.get(e, .076) for e in elements])
            dist = np.linalg.norm(xyz[:, None]-xyz[None, :], axis=-1)
            mask = (dist > .06) & (dist < cov[:, None]+cov[None, :]+.035)
            mask &= abs(data['residue_numbers'][:, None]-data['residue_numbers'][None, :]) <= 1
            pairs = np.argwhere(np.triu(mask, 1))
            obj = curves('Distance-inferred heavy-atom bonds', [xyz[pair] for pair in pairs],
                         sn.material('Bonds', profile, '#9ba6af'), .018, group)
            obj['evidence'] = 'distance candidate bonds, not CCD chemical validation; inspect before chemistry teaching'
            obj['bond_count'] = len(pairs)
    if rep in {'ribbon', 'shell_ribbon', 'wire'}:
        secondary = json.loads((asset_dir/'secondary-structure.json').read_text())['segments']
        ribbon(data, secondary, chain, profile, group, white=profile['material'] == 'xray')
    if rep in {'wire', 'volume'}:
        obj = sn.point_object('Atom-derived density', data['atoms'], group)
        sn.density(obj, mat, profile, isosurface=rep == 'wire')
        if rep == 'wire':
            dec = obj.modifiers.new('Readable mesh density', 'DECIMATE')
            dec.ratio = .2
            wire = obj.modifiers.new('Thin density grid', 'WIREFRAME')
            wire.thickness = .009
            wire.use_replace = True
    if rep in {'crowded_surface', 'environment_atoms'}:
        env = ds.collection('Complete background molecules; schematic locations')
        rng = np.random.default_rng(34)
        camera=bpy.context.scene.camera
        rotation=camera.rotation_euler.to_matrix()
        right=np.array(rotation.col[0]);up=np.array(rotation.col[1]);back=np.array(rotation.col[2])
        for i in range(42 if rep == 'crowded_surface' else 24):
            row, col = divmod(i, 5 if rep == 'crowded_surface' else 4)
            if rep=='crowded_surface':
                layer=i//25;row=(i%25)//5
                u=(col-2)*3.45+layer*1.6;v=(row-2)*3.35+layer*1.4
                depth=4+layer*4.8
            else:
                u=(col-1.5)*3.0;v=(row-2.5)*2.8;depth=5+rng.random()*3
            u,v=np.array([u,v])+rng.uniform(-.45,.45,2)
            location=right*u+up*v-back*depth
            m = sn.material(f'Environment tone {i%2}', profile, profile['palette'][1+i%2])
            copy = ds.mesh(f'Complete ubiquitin context {i+1}', data['vertices'], data['faces'], [m], env)
            copy.location = location
            copy.rotation_euler = rng.uniform(-math.pi, math.pi, 3)
            copy['evidence'] = 'same source structure, arbitrary placement for style demonstration; not measured cell composition'
    return data


def cell_surface(name, center, scale, angle, mat, group, phase=0):
    verts, faces = [], []
    rings, sides = 24, 64
    for j in range(rings+1):
        phi = math.pi*j/rings
        for i in range(sides):
            theta = math.tau*i/sides
            r = 1+.1*math.cos(3*theta+phase)+.06*math.sin(5*theta+phase)
            x, y = r*math.sin(phi)*scale[0]*math.cos(theta), r*math.sin(phi)*scale[1]*math.sin(theta)
            verts.append((center[0]+x*math.cos(angle)-y*math.sin(angle),
                          center[1]+x*math.sin(angle)+y*math.cos(angle), center[2]+scale[2]*math.cos(phi)))
    for j in range(rings):
        for i in range(sides):
            a=j*sides+i; b=j*sides+(i+1)%sides
            faces.append((a,b,b+sides,a+sides))
    obj=ds.mesh(name, verts, faces, [mat], group)
    obj['evidence']='schematic cell form; not segmented microscopy data'
    return obj


def surface_granules(obj, mat):
    t, source, target = sn.geometry_tree(obj, 'GN_Surface_distribution')
    scatter = sn.node(t, 'GeometryNodeDistributePointsOnFaces')
    scatter.distribute_method = 'RANDOM'
    scatter.inputs['Density'].default_value = 75
    scatter.inputs['Seed'].default_value = 26
    t.links.new(source, scatter.inputs['Mesh'])
    ico = sn.node(t, 'GeometryNodeMeshIcoSphere')
    ico.inputs['Radius'].default_value = .027
    ico.inputs['Subdivisions'].default_value = 1
    setter = sn.node(t, 'GeometryNodeSetMaterial')
    setter.inputs['Material'].default_value = mat
    t.links.new(ico.outputs['Mesh'], setter.inputs['Geometry'])
    ins = sn.node(t, 'GeometryNodeInstanceOnPoints')
    t.links.new(scatter.outputs['Points'], ins.inputs['Points'])
    t.links.new(scatter.outputs['Rotation'], ins.inputs['Rotation'])
    t.links.new(setter.outputs['Geometry'], ins.inputs['Instance'])
    ins.inputs['Scale'].default_value = (1, 1, 1.5)
    join = sn.node(t, 'GeometryNodeJoinGeometry')
    t.links.new(source, join.inputs['Geometry'])
    t.links.new(ins.outputs['Instances'], join.inputs['Geometry'])
    t.links.new(join.outputs['Geometry'], target)


def cells(style, profile, group, frames):
    rng=np.random.default_rng(28)
    mats=[sn.material('Cell channel '+str(i), profile, c) for i,c in enumerate(profile['palette'])]
    if style=='histology':
        centers=[((j-2)*2.45+(i%2)*.6,(i-1.5)*2.1) for i in range(4) for j in range(5)]
    elif style=='sem':centers=[(0,0)]
    elif style=='tem':centers=[(-2.5,0),(2.6,.4)]
    else:
        centers=[(-2.4,1.7),(2.25,1.45),(-1.9,-1.85),(2.75,-2)]
    puncta=[]
    for idx,(cx,cy) in enumerate(centers):
        cx,cy=np.array([cx,cy])+rng.uniform(-.22,.22,2)
        rx,ry=(1.48,1.18) if style=='histology' else (2.8,2.0) if style=='sem' else (2.5,1.8) if style=='tem' else (1.8,1.18)
        rx,ry=np.array([rx,ry])*rng.uniform(.9,1.1,2)
        angle=rng.uniform(-.6,.6)
        def xy(u,v,z):
            return (cx+u*math.cos(angle)-v*math.sin(angle),cy+u*math.sin(angle)+v*math.cos(angle),z)
        if style!='fluorescence':
            body=cell_surface(f'Whole cell {idx}',(cx,cy,.22),(rx,ry,.2),angle,mats[0],group,idx)
            if style=='sem':surface_granules(body,mats[0])
        cell_surface(f'Nucleus {idx}',xy(.1,-.03,.48),(.51,.39,.25),angle,mats[1],group,idx)
        if style in {'fluorescence','histology','tem'}:
            points=[]
            for k in range(380):
                a=rng.uniform(0,math.tau);r=math.sqrt(rng.random())
                points.append(xy(.1+.48*r*math.cos(a),-.03+.36*r*math.sin(a),.48+.25*math.sqrt(1-r*r)+.008))
            chromatin=sn.point_object(f'Nuclear texture points {idx}',points,group)
            color='#87a0df' if style=='fluorescence' else '#40244e' if style=='histology' else '#7d817e'
            sn.instance_spheres(chromatin,sn.material(f'Nuclear texture {idx}',profile,color),.012,1)
            chromatin['evidence']='artistic nuclear texture, not resolved DNA localization'
        if style=='sem':
            projections=[]
            for k in range(260):
                a=rng.uniform(0,math.tau);r=rng.uniform(.32,.94)
                x,y=rx*r*math.cos(a),ry*r*math.sin(a)
                z=.22+.2*math.sqrt(1-r*r)
                length=rng.uniform(.08,.3);direction=rng.uniform(0,math.tau)
                projections.append(np.array([xy(x+length*.6*u*u*math.cos(direction),y+length*.6*u*u*math.sin(direction),z+length*u) for u in np.linspace(0,1,12)]))
            curves('Schematic cell surface projections',projections,mats[2],.016,group)
        if style=='fluorescence':
            # Flat adhesion scaffold with bent filaments outside the nuclear zone.
            lines=[]
            for i in range(58):
                a=rng.uniform(0,math.tau);b=a+rng.uniform(1.1,4.9)
                start=np.array([math.cos(a),math.sin(a)])*rng.uniform(.78,.98)
                end=np.array([math.cos(b),math.sin(b)])*rng.uniform(.78,.98)
                side=np.array([-(end-start)[1],(end-start)[0]])
                bend=rng.uniform(-.2,.2)
                pts=[]
                for u in np.linspace(0,1,28):
                    q=(1-u)*start+u*end+side*bend*math.sin(math.pi*u)
                    if np.linalg.norm(q)<.36:q=q/max(np.linalg.norm(q),.001)*.36
                    pts.append(xy(rx*q[0],ry*q[1],.27+.08*math.sin(u*math.pi)))
                lines.append(np.array(pts))
            curves(f'Cytoskeletal channel {idx}',lines,mats[0],.012,group)
        organelles=[]
        for k in range(32 if style=='histology' else 68):
            th=rng.uniform(0,math.tau); r=rng.uniform(.43,.9)
            x,y=rx*r*math.cos(th),ry*r*math.sin(th)
            z=.22+.2*math.sqrt(max(0,1-r*r))
            pos=xy(x,y,z+.035)
            puncta.append(pos)
            if style=='tem' and k%3==0:
                theta=np.linspace(0,math.tau,60)
                loop=np.array([xy(x+.22*math.cos(t),y+.075*math.sin(t),.46) for t in theta])
                organelles.append(loop)
        if style=='tem':
            curves(f'Closed membrane outlines {idx}',organelles,mats[2],.018,group)
    granule_mat=mats[2] if style in {'tem','histology','fluorescence'} else mats[0]
    obj=sn.point_object('Intracellular puncta',puncta,group)
    sn.instance_spheres(obj,granule_mat,.035 if style!='histology' else .024,2,
                        loop_frames=frames,motion=.05 if style=='fluorescence' else 0)
    obj['evidence']='schematic particles and periodic style-test motion; not measured localization or transport'


def setup(profile, width, samples, frames, is_cells):
    scene=bpy.context.scene
    scene.render.engine=profile['engine']
    scene.render.resolution_x=width
    scene.render.resolution_y=2*round(width/3)
    scene.render.resolution_percentage=100
    scene.render.image_settings.file_format='PNG'
    scene.render.film_transparent=True
    scene.render.fps=24
    scene.frame_start,scene.frame_end=1,frames
    scene.view_settings.view_transform=profile.get('view_transform','Standard' if profile['material'] in {'toon','paint','height','emission','stain'} else 'AgX')
    scene.view_settings.look='None'
    scene.view_layers[0].use_pass_normal=True
    scene.view_layers[0].use_pass_z=True
    scene.view_layers[0].use_pass_cryptomatte_object=True
    if scene.render.engine=='CYCLES':
        scene.cycles.samples=samples
        scene.cycles.use_denoising=True
        scene.cycles.seed=19
        scene.cycles.use_animated_seed=False
        scene.cycles.transparent_max_bounces=16
        try:
            prefs=bpy.context.preferences.addons['cycles'].preferences
            prefs.compute_device_type='METAL';prefs.get_devices()
            metal=[d for d in prefs.devices if d.type=='METAL']
            for device in prefs.devices:device.use=device in metal
            if metal:scene.cycles.device='GPU'
        except (TypeError,RuntimeError):pass
    world=bpy.data.worlds.new('Controlled fill')
    scene.world=world;world.use_nodes=True
    world.node_tree.nodes['Background'].inputs['Color'].default_value=(.35,.4,.5,1)
    world.node_tree.nodes['Background'].inputs['Strength'].default_value=.24
    bpy.ops.object.camera_add(location=(5,-7,15) if profile['material']=='sem' else (0,0,24) if is_cells or profile['representation']=='crowded_surface' else (6,-9,5.2))
    camera=bpy.context.object;camera.name='Reference comparison camera'
    ds.aim(camera,(0,0,0))
    camera.data.type='ORTHO';camera.data.ortho_scale=8 if profile['material']=='sem' else 13 if is_cells else 17 if profile['representation']=='crowded_surface' else 6.8
    scene.camera=camera
    for name,loc,power,size in [('Key',(-4,-5,8),650,5),('Fill',(5,-1,4),330,5),('Rim',(1,5,6),850,3)]:
        bpy.ops.object.light_add(type='AREA',location=loc)
        light=bpy.context.object;light.name=name;light.data.energy=power;light.data.shape='DISK';light.data.size=size
        ds.aim(light,(0,0,0))
    if profile['representation']=='environment_atoms':
        camera.data.dof.use_dof=True
        camera.data.dof.focus_distance=camera.location.length
        camera.data.dof.aperture_fstop=4
    # Both branches remain editable in the saved project.
    sn.compositor(scene,profile,'clean')
    if profile['compositor']!='clean':sn.compositor(scene,profile)
    return scene


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--asset-dir',type=Path)
    p.add_argument('--chain',default='A',help='author chain IDs separated by commas, or all')
    p.add_argument('--profiles',type=Path,default=Path(__file__).resolve().parents[1]/'assets/style-library.json')
    p.add_argument('--styles',default='clay',help='comma-separated IDs or all')
    p.add_argument('--engine',choices=['profile','eevee','cycles'],default='profile')
    p.add_argument('--compositor',default='profile',choices=['profile','clean','ink','filmic','fluorescence','grain','dic','phase','soft'])
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--width',type=int,default=960)
    p.add_argument('--samples',type=int,default=48)
    p.add_argument('--frames',type=int,default=96)
    p.add_argument('--render-frames',default='1',help='comma-separated frames; frame N+1 can check the seam')
    p.add_argument('--animate',action='store_true')
    args=p.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else [])
    if args.output.exists():p.error('Choose a fresh output directory.')
    if args.output.resolve().is_relative_to(Path(__file__).resolve().parents[1]):
        p.error('Write generated projects outside the skill directory.')
    profiles=json.loads(args.profiles.read_text())['profiles']
    styles=list(profiles) if args.styles=='all' else args.styles.split(',')
    if set(styles)-set(profiles):p.error('Unknown style IDs: '+str(set(styles)-set(profiles)))
    if any(profiles[s]['subject']=='molecule' for s in styles) and not args.asset_dir:
        p.error('Molecular styles require --asset-dir; see prepare_style_asset.py.')
    if args.width<128 or args.width%2 or args.samples<1 or args.frames<2:
        p.error('Use an even width >=128, positive samples and at least two frames.')
    args.output.mkdir(parents=True)
    for style in styles:
        start=time.monotonic()
        bpy.ops.wm.read_factory_settings(use_empty=True)
        profile=dict(profiles[style])
        if args.engine!='profile':profile['engine']='CYCLES' if args.engine=='cycles' else 'BLENDER_EEVEE'
        if args.compositor!='profile':profile['compositor']=args.compositor
        if profile['material']=='toon' and profile['engine']!='BLENDER_EEVEE':
            p.error('The illustrative Shader-to-RGB material requires Eevee; choose paint or PBR for Cycles.')
        directory=args.output/style;directory.mkdir()
        scene=setup(profile,args.width,args.samples,args.frames,profile['subject']=='cells')
        group=ds.collection('Biological subject')
        if profile['subject']=='molecule':
            source_manifest=json.loads((args.asset_dir/'manifest.json').read_text())
            chains=[c['chain'] for c in source_manifest['chains']] if args.chain=='all' else args.chain.split(',')
            source_points=[]
            for index,chain in enumerate(chains):
                chain_profile=dict(profile)
                palette=profile['palette']
                chain_profile['palette']=palette[index%len(palette):]+palette[:index%len(palette)]
                before=set(group.objects)
                source_data=molecule(args.asset_dir,chain,chain_profile,group)
                source_points.append(source_data['vertices'])
                for obj in set(group.objects)-before:obj['source_chain']=chain
            if profile['representation']!='crowded_surface':
                # Fit actual coordinate-derived points instead of inflated
                # rotated bounding boxes, which over-pad multichain complexes.
                camera=scene.camera
                rotation=np.array(camera.rotation_euler.to_matrix())
                points=np.concatenate(source_points)
                poses=[]
                for angle in np.linspace(-.18,.18,9):
                    turn=np.array([[math.cos(angle),-math.sin(angle),0],
                                   [math.sin(angle),math.cos(angle),0],[0,0,1]])
                    poses.append((points@turn.T-np.array(camera.location))@rotation)
                projected=np.concatenate(poses)
                low,high=projected.min(0),projected.max(0)
                center=(low+high)/2
                camera.location+=Vector(rotation@np.array([center[0],center[1],0]))
                aspect=scene.render.resolution_x/scene.render.resolution_y
                camera.data.ortho_scale=max(high[0]-low[0],(high[1]-low[1])*aspect)/.8
            root=bpy.data.objects.new('Rigid style-preview motion',None)
            scene.collection.objects.link(root)
            for obj in list(group.objects):obj.parent=root
            driver=root.driver_add('rotation_euler',2).driver
            driver.expression=f'0.18*sin(2*pi*(frame-1)/{args.frames})'
        else:cells(style,profile,group,args.frames)
        scene['style_id']=style
        scene['source_units']='nm' if profile['subject']=='molecule' else 'arbitrary display units; schematic cell morphology'
        scene['evidence']='structure-derived representations or explicitly schematic cell forms; artistic reference matching'
        scene['animation_scope']='periodic style stability preview, not a molecular mechanism or biological simulation'
        scene.frame_set(1)
        bpy.ops.wm.save_as_mainfile(filepath=str(directory/'scene.blend'))
        frames=list(range(1,args.frames+1)) if args.animate else [int(x) for x in args.render_frames.split(',')]
        render_dir=directory/'frames';render_dir.mkdir()
        for frame in frames:
            scene.frame_set(frame)
            scene.render.filepath=str(render_dir/f'frame_{frame:04d}.png')
            bpy.ops.render.render(write_still=True)
        report={'style':style,'blender':bpy.app.version_string,'engine':scene.render.engine,'profile':profile,
                'frames':frames,'seconds':time.monotonic()-start,'resolution':[scene.render.resolution_x,scene.render.resolution_y],
                'asset_source':json.loads((args.asset_dir/'manifest.json').read_text()) if profile['subject']=='molecule' else None,
                'nodes':[{'name':t.name,'type':t.bl_idname,'nodes':[n.bl_idname for n in t.nodes]} for t in bpy.data.node_groups],
                'visual_review':'required separately; a successful render is not reference-match validation'}
        (directory/'build-report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
        print('STYLE_RENDER_COMPLETE',style,round(report['seconds'],2),flush=True)


if __name__=='__main__':main()
