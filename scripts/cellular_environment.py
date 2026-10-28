"""Source-shaped depth layers; artistic defocus stays separate from foreground."""
import json
from pathlib import Path
import bpy
import numpy as np
from mathutils import Matrix,Euler
import detail_scene as ds
import style_nodes as sn


def build(asset_dir,scale_nm,frames):
    asset_dir=Path(asset_dir);source=json.loads((asset_dir/'manifest.json').read_text())
    if source['kind']!='schematic-cytosolic-environment':raise ValueError('Prepare the cellular environment assets first')
    groups=[ds.collection('Cellular environment '+name) for name in ['far','middle','near']]
    rng=np.random.default_rng(741)
    colors={'ribosome':'#638e8d','ubiquitin':'#728d9b','actin':'#8b829a'}
    materials={name:sn.material('Background '+name,{'material':'pbr','roughness':.62,'bump':0},color) for name,color in colors.items()}
    data={name:dict(np.load(asset_dir/(name+'.npz'))) for name in colors}
    meshes={}
    for name,d in data.items():
        mesh=bpy.data.meshes.new('Background intact '+name);mesh.from_pydata(d['vertices'].tolist(),[],d['faces'].tolist());mesh.update()
        mesh.materials.append(materials[name]);meshes[name]=mesh
        for polygon in mesh.polygons:polygon.use_smooth=True
    records=[]
    def place(name,layer,center,projection,angles,index,local=None):
        obj=bpy.data.objects.new(f'Environment {name} {index}',meshes[name]);groups[layer].objects.link(obj)
        matrix=Matrix.Translation(center)@Euler(angles).to_matrix().to_4x4()@Matrix.Scale(scale_nm*projection,4)
        if local is not None:matrix=matrix@Matrix(local)
        obj.matrix_world=matrix;obj['environment_only']=True;obj['source_structure']=source['sources'][name]['id']
        obj['depth_layer']=layer;obj['evidence']='schematic placement of source-derived shape'
        records.append({'object':obj.name,'asset':name,'depth_layer':layer,'projection_factor':projection,'matrix':np.array(matrix).tolist()})
        return obj
    for i,(center,projection,angles) in enumerate([
        ((-10,3,-18),.55,(.5,.7,.6)),((10.6,-3.1,-18),.55,(1.8,.1,-.3)),
        ((3.5,7.4,-25),.29,(.8,1.5,.2))]):place('ribosome',0,center,projection,angles,i)
    # One shared matrix transforms every actin subunit, preserving its helical assembly.
    for fiber,(center,angles,projection) in enumerate([
        ((-10,-4,-13),(.2,-.4,-.45),.42),((8.3,3.2,-20),(1.3,.8,.35),.32)]):
        steps=[np.linalg.matrix_power(data['actin']['repeat_nm'],i) for i in range(-6,7)]
        for i,step in enumerate(steps):place('actin',1,center,projection,angles,fiber*13+i,step)
    for i in range(22):
        layer=0 if i<10 else 1 if i<19 else 2
        center=(float(rng.uniform(-13,13)),float(rng.uniform(-7.5,7.5)),[-24,-14,-6][layer])
        projection=[.42,.70,1.2][layer]
        place('ubiquitin',layer,center,projection,rng.uniform(0,6.28,3),i)
    # Slow, closed authored drift. It is background motion, not Brownian dynamics.
    for index,group in enumerate(groups):
        parent=bpy.data.objects.new(group.name+' drift',None);group.objects.link(parent)
        for obj in list(group.objects):
            if obj!=parent:obj.parent=parent
        for frame in np.linspace(1,frames+1,17):
            phase=(frame-1)/frames*2*np.pi
            parent.location=(.13*(np.cos(phase)-1),.10*np.sin(phase),0)
            parent.keyframe_insert('location',frame=float(frame))
    return groups,{'asset_directory':str(asset_dir.resolve()),'source_manifest':source,'objects':records,
        'layout':'Illustrative cytosolic depth backdrop, separate from the enlarged mechanism panels; not measured concentration or a common physical scale bar.',
        'defocus':'Three depth planes composited with separately editable blur widths; an artistic depth-of-field approximation, not a microscope PSF.',
        'motion':'Slow closed authored drift, not a diffusion simulation.'}


def composite(tree,layers,background,width):
    color=sn.node(tree,'CompositorNodeRGB','Cellular environment base color')
    color.outputs[0].default_value=(*ds.linear(background),1);image=color.outputs[0]
    for layer,pixels,strength in zip(layers,[18,9,28],[.26,.24,.14]):
        render=sn.node(tree,'CompositorNodeRLayers');render.layer=layer.name
        blur=sn.node(tree,'CompositorNodeBlur',layer.name+' defocus')
        blur.inputs['Size'].default_value=(pixels*width/1280,pixels*width/1280)
        tree.links.new(render.outputs['Image'],blur.inputs['Image'])
        over=sn.node(tree,'CompositorNodeAlphaOver',layer.name+' background strength')
        over.inputs['Factor'].default_value=strength
        tree.links.new(image,over.inputs['Background']);tree.links.new(blur.outputs['Image'],over.inputs['Foreground']);image=over.outputs['Image']
    return image
