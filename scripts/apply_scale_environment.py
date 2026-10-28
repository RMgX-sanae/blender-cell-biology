"""Add scale-specific depth backgrounds to the skill's existing scene builders.

Background cells use Allen meshes; ER uses reference-guided 3D curves or an image backplate.
Membrane/Ras shapes are retained.
Run in Blender 5.x on a saved source scene; always write a fresh project folder.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Euler, Matrix, Vector

sys.path.insert(0, str(Path(__file__).resolve().parent))
import detail_scene as ds
import style_nodes as sn


def geometry_signature(objects):
    digest = hashlib.sha256()
    for obj in sorted(objects, key=lambda o: o.name):
        if obj.type not in {'MESH', 'CURVE', 'FONT', 'EMPTY'}:
            continue
        digest.update(obj.name.encode())
        digest.update(np.asarray(obj.matrix_world, dtype=np.float64).tobytes())
        if obj.type == 'MESH':
            coords = np.empty(len(obj.data.vertices) * 3, dtype=np.float32)
            obj.data.vertices.foreach_get('co', coords)
            digest.update(coords.tobytes())
    return digest.hexdigest()


def surface(name, color, opacity=1):
    mat = ds.material(name, color, .65, .008)
    shader = mat.node_tree.nodes.get('Principled BSDF')
    shader.inputs['Emission Color'].default_value = (*ds.linear(color), 1)
    shader.inputs['Emission Strength'].default_value = .18
    if opacity < 1:
        tree = mat.node_tree
        transparent = tree.nodes.new('ShaderNodeBsdfTransparent')
        mix = tree.nodes.new('ShaderNodeMixShader')
        mix.inputs[0].default_value = opacity
        tree.links.new(transparent.outputs[0], mix.inputs[1])
        tree.links.new(shader.outputs[0], mix.inputs[2])
        tree.links.new(mix.outputs[0], tree.nodes.get('Material Output').inputs['Surface'])
    return mat


def spatial_background(scene, scale, groups, asset_dir):
    asset_dir=Path(asset_dir)
    manifest=json.loads((asset_dir/'manifest.json').read_text())
    if manifest['kind']=='reference-guided-er-network':
        if scale!='organelle':raise ValueError('The traced ER environment is for the organelle scene.')
        return reference_network(scene,groups,asset_dir,manifest)
    if manifest['kind']!='allen-visual-guide-environment':
        raise ValueError('Run prepare_allen_environment.py first.')
    if scale!='cell':raise ValueError('Use prepare_reference_er_network.py for the fine organelle environment.')
    sources={entry['name']:entry for entry in manifest['assets']}
    data={name:dict(np.load(asset_dir/(name+'.npz'))) for name in sources}
    camera=scene.camera
    camera.data.ortho_scale*=1.20 if scale=='cell' else 1.12
    width=camera.data.ortho_scale;depth=camera.location.length
    records=[]
    def place(name, group, transform, material):
        d=data[name]
        obj=ds.mesh('Allen background '+name,d['vertices'].tolist(),d['faces'].tolist(),[material],group)
        obj.matrix_world=transform
        obj['source_sha256']=sources[name]['sha256']
        obj['evidence']='Public display-model geometry preserved; source methodology and resolution limits retained.'
        records.append({'object':obj.name,'asset':name,'matrix':list(map(list,transform))})
    def fit(names, angles, fraction, xy, behind):
        pose=Euler(angles).to_matrix().to_4x4()
        vertices=np.concatenate([data[name]['vertices'] for name in names])
        rotated=vertices @ np.asarray(pose)[:3,:3].T
        lo,hi=rotated.min(0),rotated.max(0);center=(lo+hi)/2
        factor=width*fraction/(hi[0]-lo[0])
        local=Matrix.Translation((xy[0]*width,xy[1]*width,-depth-behind*width))
        return camera.matrix_world @ local @ Matrix.Scale(factor,4) @ Matrix.Translation(-Vector(center)) @ pose
    if scale=='cell':
        for prefix,angles,fraction,xy,behind,layer in [
                ('er',(1.08,.12,.47),.64,(-.47,.19),.75,1),
                ('mito',(1.21,-.20,-.39),.62,(.23,.38),.98,0),
                ('golgi',(.96,.14,.65),.60,(.46,-.41),-.20,2)]:
            names=[prefix+'-'+suffix for suffix in ['cell','dna','structure']]
            transform=fit(names,angles,fraction,xy,behind)
            membrane=surface('Source cell membrane '+prefix,'88A5A0',.2)
            # A Fresnel-weighted shell reveals the paired source nucleus/marker.
            tree=membrane.node_tree;fresnel=tree.nodes.new('ShaderNodeFresnel')
            fresnel.inputs['IOR'].default_value=1.38
            opacity=tree.nodes.new('ShaderNodeMath');opacity.operation='MULTIPLY_ADD'
            opacity.inputs[1].default_value=.65;opacity.inputs[2].default_value=.10
            tree.links.new(fresnel.outputs[0],opacity.inputs[0])
            tree.links.new(opacity.outputs[0],next(n for n in tree.nodes if n.bl_idname=='ShaderNodeMixShader').inputs[0])
            mats=[membrane,surface('Source nuclear material '+prefix,'83859A',.7),
                  surface('Source intracellular marker '+prefix,'829E96',.65)]
            for name,mat in zip(names,mats):place(name,groups[layer],transform,mat)
        contents='Three distinct interphase cell display meshes from Allen Cell, retaining each membrane, DNA and organelle-marker shape and their shared source coordinates.'
    return contents,{'source_manifest':manifest,'placements':records,
                     'layout':'Source shapes retained under rigid transforms and uniform scale; camera placement is authored, not a measured neighborhood around this foreground model.'}


def reference_network(scene,groups,asset_dir,manifest):
    path=asset_dir/'network.npz'
    if hashlib.sha256(path.read_bytes()).hexdigest()!=manifest['geometry_sha256']:
        raise ValueError('Prepared ER network does not match its manifest.')
    source=np.load(path);points=source['points'];radii=source['radii'];offsets=source['offsets']
    camera=scene.camera;camera.data.ortho_scale*=1.12;width=camera.data.ortho_scale
    curve=bpy.data.curves.new('Reference-guided ER branches','CURVE')
    curve.dimensions='3D';curve.resolution_u=4;curve.bevel_depth=1;curve.bevel_resolution=3;curve.use_fill_caps=True
    for start,end in zip(offsets[:-1],offsets[1:]):
        spline=curve.splines.new('BEZIER');spline.bezier_points.add(int(end-start)-1)
        for i,(point,co) in enumerate(zip(spline.bezier_points,points[start:end])):
            point.co=co;point.radius=float(radii[start+i]);point.handle_left_type='AUTO';point.handle_right_type='AUTO'
    obj=bpy.data.objects.new('ER network with authored depth',curve);groups[1].objects.link(obj)
    transform=camera.matrix_world @ Matrix.Translation((.03*width,.01*width,-camera.location.length-width*.43)) @ Matrix.Scale(width*1.35,4) @ Euler((.25,-.15,-.10)).to_matrix().to_4x4()
    obj.matrix_world=transform;obj['evidence']=manifest['evidence']
    mat=surface('Reference-guided ER membrane','829784')
    mat.node_tree.nodes.get('Principled BSDF').inputs['Emission Strength'].default_value=.02
    curve.materials.append(mat)
    light_data=bpy.data.lights.new('ER volume side light','AREA');light_data.energy=480;light_data.shape='DISK';light_data.size=width*.8
    light=bpy.data.objects.new('ER volume side light',light_data);groups[1].objects.link(light)
    light.matrix_world=camera.matrix_world @ Matrix.Translation((-width*.4,width*.6,-camera.location.length-width*.16))
    return 'Reference-guided three-dimensional ER curves with local branch-density variation and authored depth.', {
        'source_manifest':manifest,'matrix':list(map(list,transform)),'curve_paths':len(curve.splines),
        'layout':'Illustrative 3D geometry derived from image XY morphology; not a measured mitochondrial neighborhood.'}


def membrane_background(source, groups):
    # Re-layer existing whole molecules without altering their coordinates.
    for obj in list(bpy.context.scene.objects):
        group = None
        if 'background complex' in obj.name:
            group = groups[0]
        elif obj.name.startswith('Foreground ribosome'):
            group = groups[2]
        elif obj.name.startswith(('Background ribosome', 'Ubiquitin')):
            group = groups[1]
        if group:
            for old in list(obj.users_collection):
                old.objects.unlink(obj)
            group.objects.link(obj)
    manifest = json.loads((source.parent/'manifest.json').read_text())
    heads = np.asarray(manifest['lipid_positions_nm_and_leaflet'])
    heads = heads[heads[:,3] == 1,:3]
    glycan = surface('Schematic extracellular glycan', 'A7BC9F')
    for i,(x,y) in enumerate([(-15,8),(-10,20),(-3,20),(5,25),(13,12),(17,24)]):
        center = heads[np.argmin(np.linalg.norm(heads[:,:2]-(x,y),axis=1))]
        sugars=ds.Parts()
        stem=[center+(.18*math.sin(j*2),0,j*.60) for j in range(5)]
        ds.path('Extracellular glycan stem '+str(i),stem,.07,glycan,groups[0])
        for point in stem[1:]:sugars.ball(point,.20,seg=10,rings=6)
        for branch in (-1,1):
            pts = [stem[-1]+(branch*.4*j,.1*math.sin(j*3+i),.57*j) for j in range(4)]
            ds.path('Extracellular glycan branch '+str(i),pts,.07,glycan,groups[0])
            for point in pts[1:]:sugars.ball(point,.20,seg=10,rings=6)
        sugars.object('Glycan residue motifs '+str(i),[glycan],groups[0])
    return 'Existing 3SN6 complexes and cytosolic 1UBQ/9SPI shapes; six schematic glycan scaffolds on the extracellular leaflet, not measured composition.'


def layer_compositor(scene, groups, style, width, scale):
    base = scene.view_layers[0]
    for child in base.layer_collection.children:
        if child.collection in groups:
            child.exclude = True
    layers = []
    for group in groups:
        layer = scene.view_layers.new(group.name+' render')
        for child in layer.layer_collection.children:
            is_studio = child.collection not in groups and any(o.type in {'LIGHT','CAMERA'} for o in child.collection.objects)
            child.exclude = child.collection != group and not is_studio
        layers.append(layer)
    tree = bpy.data.node_groups.new('COMP_scale_environment_'+style, 'CompositorNodeTree')
    tree.interface.new_socket(name='Image',in_out='OUTPUT',socket_type='NodeSocketColor')
    color = sn.node(tree,'CompositorNodeRGB','Scale environment base')
    color.outputs[0].default_value = (*ds.linear('DCE5E4' if style=='light' else '102632'),1)
    image = color.outputs[0]
    def subject_over(background):
        render=sn.node(tree,'CompositorNodeRLayers','Sharp subject');render.layer=base.name
        over=sn.node(tree,'CompositorNodeAlphaOver','Sharp subject between environmental planes')
        tree.links.new(background,over.inputs['Background']);tree.links.new(render.outputs['Image'],over.inputs['Foreground'])
        return over.outputs['Image']
    for layer,pixels,strength in zip(layers,[13,6,20], [.38,.42,.22] if style=='dark' else [.32,.35,.18]):
        if scale=='organelle':
            strength=.27 if style=='dark' else .23
            pixels=10
        if scale=='cell':
            pixels=[22,13,30][layers.index(layer)]
            strength*=.78
            if layer==layers[2]:image=subject_over(image)
        render = sn.node(tree,'CompositorNodeRLayers');render.layer=layer.name
        blur = sn.node(tree,'CompositorNodeBlur',layer.name+' defocus')
        blur.inputs['Size'].default_value = (pixels*width/1280, pixels*width/1280)
        tree.links.new(render.outputs['Image'],blur.inputs['Image'])
        over = sn.node(tree,'CompositorNodeAlphaOver',layer.name+' strength')
        over.inputs['Factor'].default_value = strength
        tree.links.new(image,over.inputs['Background'])
        tree.links.new(blur.outputs['Image'],over.inputs['Foreground']);image=over.outputs['Image']
    if scale!='cell':image=subject_over(image)
    output=sn.node(tree,'NodeGroupOutput');tree.links.new(image,output.inputs['Image'])
    scene.compositing_node_group=tree


def reference_background(scene, image_path, profile_path, style, width, scale):
    profile=json.loads(profile_path.read_text())
    digest=hashlib.sha256(image_path.read_bytes()).hexdigest()
    if digest!=profile['sha256']:raise ValueError('Reference image does not match its source profile.')
    image=bpy.data.images.load(str(image_path.resolve()),check_existing=True)
    x,y,w,h=profile['crop_bottom_left_xywh']
    if min(x,y)<0 or min(w,h)<=0 or x+w>image.size[0] or y+h>image.size[1]:
        raise ValueError('Reference crop is outside the image.')
    if abs(w/h-scene.render.resolution_x/scene.render.resolution_y)>.001:
        raise ValueError('Choose a crop with the output aspect ratio; do not stretch microscopy.')
    image.pack()
    scene.camera.data.ortho_scale*=1.20 if scale=='cell' else 1.12
    tree=bpy.data.node_groups.new('COMP_reference_environment_'+style,'CompositorNodeTree')
    tree.interface.new_socket(name='Image',in_out='OUTPUT',socket_type='NodeSocketColor')
    source=sn.node(tree,'CompositorNodeImage','Public microscopy reference');source.image=image
    crop=sn.node(tree,'CompositorNodeCrop','Context crop; original scale bar is not a scene scale')
    for key,value in zip(['X','Y','Width','Height'],[x,y,w,h]):crop.inputs[key].default_value=value
    tree.links.new(source.outputs['Image'],crop.inputs['Image'])
    resize=sn.node(tree,'CompositorNodeScale','Uniform image scale')
    resize.inputs['X'].default_value=width/w;resize.inputs['Y'].default_value=width/w
    tree.links.new(crop.outputs['Image'],resize.inputs['Image'])
    gray=sn.node(tree,'CompositorNodeRGBToBW');tree.links.new(resize.outputs['Image'],gray.inputs['Image'])
    gamma=sn.node(tree,'ShaderNodeMath','Artistic grade, not quantitative signal');gamma.operation='POWER'
    gamma.inputs[1].default_value=profile['gamma'];tree.links.new(gray.outputs[0],gamma.inputs[0])
    ramp=sn.node(tree,'ShaderNodeValToRGB','Muted environmental colors')
    for element,color in zip(ramp.color_ramp.elements,profile['palette'][style]):
        element.color=(*ds.linear(color),1)
    ramp.color_ramp.elements[1].position=profile['white_point']
    tree.links.new(gamma.outputs[0],ramp.inputs[0])
    blur=sn.node(tree,'CompositorNodeBlur','Fixed reference-plane defocus')
    blur.inputs['Size'].default_value=(profile['blur_at_1280']*width/1280,)*2
    tree.links.new(ramp.outputs[0],blur.inputs['Image'])
    foreground=sn.node(tree,'CompositorNodeRLayers','Sharp subject');foreground.layer=scene.view_layers[0].name
    over=sn.node(tree,'CompositorNodeAlphaOver','Subject over reference image')
    tree.links.new(blur.outputs['Image'],over.inputs['Background'])
    tree.links.new(foreground.outputs['Image'],over.inputs['Foreground'])
    output=sn.node(tree,'NodeGroupOutput');tree.links.new(over.outputs['Image'],output.inputs['Image'])
    scene.compositing_node_group=tree
    return 'Cropped, recolored and defocused microscopy image used as a fixed illustrative background.', {
        'kind':'reference-image-backplate','source':profile,'source_sha256':digest,
        'limitation':'Two-dimensional background; no recovered depth, measured foreground relationship, quantitative intensity or camera-orbit parallax.'}


def molecular_style(scene, style, width):
    tree=scene.compositing_node_group
    base=next(n for n in tree.nodes if n.label=='Cellular environment base color')
    base.outputs[0].default_value=(*ds.linear('DCE5E4' if style=='light' else '071F2B'),1)
    if style=='light':
        ink=sn.material('Dark labels for light background',{'material':'paint'},'#243943','paint')
        for obj in bpy.data.collections['Typography'].objects:
            if obj.type=='FONT' or obj.name=='Panel separator':
                obj.data.materials.clear();obj.data.materials.append(ink)
        for key,label,color in [('ras','Ras','#367B96'),('sos','SOS','#AF5973'),
                                ('old_gdp','GDP','#3A86AD'),('new_gtp','GTP','#238F6D'),
                                ('gamma','γ-P','#BA811F'),('old_mg','Mg2+','#846399'),
                                ('new_mg',None,'#846399')]:
            for name in [key]+([label+' legend'] if label else []):
                old=bpy.data.materials.get(name)
                if old:
                    replacement=sn.material('Light theme '+name,{'material':'pbr','roughness':.44,'bump':0},color)
                    old.user_remap(replacement)
        for n in tree.nodes:
            if n.label=='Intact envelopes with visibility overlay':n.inputs[0].default_value=.85
    for n in tree.nodes:
        if n.bl_idname=='CompositorNodeBlur' and 'defocus' in n.label:
            pixels=18 if 'far' in n.label else 9 if 'middle' in n.label else 28
            n.inputs['Size'].default_value=(pixels*width/1280,pixels*width/1280)
    return 'Source-derived 1UBQ, 9SPI and 8F8R shapes in three existing cytosolic background layers; independently magnified mechanism views.'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scale',choices=['cell','organelle','membrane','molecule'],required=True)
    p.add_argument('--style',choices=['light','dark'],required=True)
    p.add_argument('--engine',choices=['cycles','eevee'],default='cycles')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--context-assets',type=Path,help='Prepared Allen cell meshes or reference-guided ER network.')
    p.add_argument('--reference-image',type=Path,help='Optional fixed microscopy backplate for cell/organelle scenes.')
    p.add_argument('--reference-profile',type=Path,help='Source, hash, crop and grading parameters for --reference-image.')
    p.add_argument('--width',type=int,default=1280)
    p.add_argument('--samples',type=int,default=32)
    p.add_argument('--frame',type=int,default=1)
    p.add_argument('--render',action='store_true')
    args=p.parse_args(sys.argv[sys.argv.index('--')+1:])
    if args.output.exists():p.error('Choose a fresh output directory.')
    if bool(args.reference_image)!=bool(args.reference_profile):p.error('Supply both --reference-image and --reference-profile.')
    if args.reference_image and args.scale not in {'cell','organelle'}:p.error('The reference-image adapter supports the cell/organelle source builders.')
    if args.reference_image and args.context_assets:p.error('Choose a reference image or prepared 3D context assets.')
    if args.scale in {'cell','organelle'} and not (args.context_assets or args.reference_image):p.error('Supply --context-assets or a reference image/profile pair.')
    source=Path(bpy.data.filepath)
    if not source.is_file():p.error('Open a saved source scene before applying this adapter.')
    scene=bpy.context.scene;scene.frame_set(args.frame);bpy.context.view_layer.update()
    originals=list(scene.objects);before=geometry_signature(originals)
    scene.render.resolution_x=args.width
    scene.render.resolution_y=round(args.width*(9/16 if args.scale=='molecule' else 5/8))
    scene.render.resolution_percentage=100;scene.render.film_transparent=True
    context_record=None
    if args.reference_image:
        contents,context_record=reference_background(scene,args.reference_image,args.reference_profile,args.style,args.width,args.scale)
    elif args.scale=='molecule':
        contents=molecular_style(scene,args.style,args.width)
    else:
        groups=[ds.collection('Scale environment '+name) for name in ['far','middle','near']]
        if args.scale=='membrane':contents=membrane_background(source,groups)
        else:contents,context_record=spatial_background(scene,args.scale,groups,args.context_assets)
        layer_compositor(scene,groups,args.style,args.width,args.scale)
    scene.render.engine='CYCLES' if args.engine=='cycles' else 'BLENDER_EEVEE'
    if args.engine=='cycles':
        scene.cycles.samples=args.samples;scene.cycles.use_denoising=True
        prefs=bpy.context.preferences.addons['cycles'].preferences
        try:
            prefs.compute_device_type='METAL';prefs.get_devices()
            gpu=[d for d in prefs.devices if d.type=='METAL']
            for d in prefs.devices:d.use=d in gpu
            scene.cycles.device='GPU' if gpu else 'CPU'
        except TypeError:
            scene.cycles.device='CPU'
    after=geometry_signature(originals)
    if before!=after:raise RuntimeError('Existing subject geometry or placement changed.')
    args.output.mkdir(parents=True)
    record={'source_scene':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'scale':args.scale,'style':args.style,'contents':contents,'context':context_record,'original_geometry_unchanged':before==after,
            'geometry_signature':after,'frame':args.frame,'resolution':[args.width,scene.render.resolution_y],
            'engine':scene.render.engine,'samples':args.samples,'blender':bpy.app.version_string,
            'defocus':'Editable artistic depth-plane blur, not a calibrated microscope PSF.',
            'scope':'Scale-dependent illustration; counts and placement are authored, not measured.',
            'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    scene['scale_environment']=json.dumps(record)
    scene.render.filepath=str((args.output/'preview.png').resolve())
    bpy.ops.wm.save_as_mainfile(filepath=str((args.output/'scene.blend').resolve()))
    if args.render:bpy.ops.render.render(write_still=True)
    (args.output/'environment-style.json').write_text(json.dumps(record,ensure_ascii=False,indent=2))
    print(json.dumps({k:record[k] for k in ['scale','style','engine','resolution','original_geometry_unchanged']},ensure_ascii=False))


if __name__=='__main__':main()
