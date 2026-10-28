"""Editable Geometry Nodes, materials and compositors for reference-led styling.

Validated with Blender 5.2.3. These are appearance tools, not optical forward
models or biological dynamics. Import the small functions in a project builder;
build_style_scene.py is a reproducible demonstration, not a universal cell model.
"""
import math

import bpy
import numpy as np

from detail_scene import linear


def node(tree, kind, label=None):
    n = tree.nodes.new(kind)
    n.label = label or ''
    # Readable left-to-right layout; reusable groups retain their own labels.
    count = len(tree.nodes) - 1
    n.location = ((count % 6) * 210, -(count // 6) * 260)
    return n


def math_node(tree, operation, a=None, b=None):
    n = node(tree, 'ShaderNodeMath')
    n.operation = operation
    for i, v in enumerate((a, b)):
        if v is None:
            continue
        if isinstance(v, (int, float)):
            n.inputs[i].default_value = v
        else:
            tree.links.new(v, n.inputs[i])
    return n.outputs[0]


def ramp(tree, source, stops, label=''):
    n = node(tree, 'ShaderNodeValToRGB', label)
    r = n.color_ramp
    r.interpolation = 'LINEAR'
    for i, (position, color) in enumerate(stops):
        el = r.elements[i] if i < 2 else r.elements.new(position)
        el.position = position
        el.color = (*linear(color), 1) if isinstance(color, str) else color
    tree.links.new(source, n.inputs[0])
    return n.outputs['Color']


def material(name, profile, color=None, model=None):
    color = color or profile['palette'][0]
    model = model or profile['material']
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    m.diffuse_color = (*linear(color), 1)
    t = m.node_tree
    t.nodes.clear()
    out = node(t, 'ShaderNodeOutputMaterial')
    noise = node(t, 'ShaderNodeTexNoise', 'Object-space texture; not measured molecular detail')
    tex = node(t, 'ShaderNodeTexCoord')
    t.links.new(tex.outputs['Object'], noise.inputs['Vector'])
    noise.inputs['Scale'].default_value = 24 if model in {'sem', 'stain', 'emission'} else 8
    noise.inputs['Detail'].default_value = 3
    base = tuple(linear(color))
    if model == 'volume':
        info = node(t, 'ShaderNodeVolumeInfo', 'Illustrative density generated from atom coordinates')
        vol = node(t, 'ShaderNodeVolumePrincipled')
        t.links.new(math_node(t, 'MULTIPLY', info.outputs['Density'], 1.8), vol.inputs['Density'])
        vol.inputs['Color'].default_value = (*base, 1)
        vol.inputs['Anisotropy'].default_value = .12
        vol.inputs['Emission Color'].default_value = (*base, 1)
        t.links.new(math_node(t, 'MULTIPLY', info.outputs['Density'], .3), vol.inputs['Emission Strength'])
        t.links.new(vol.outputs[0], out.inputs['Volume'])
        return m
    if model in {'paint', 'stain', 'emission', 'height', 'toon'}:
        emit = node(t, 'ShaderNodeEmission')
        if model == 'height':
            geo = node(t, 'ShaderNodeNewGeometry')
            xyz = node(t, 'ShaderNodeSeparateXYZ')
            t.links.new(geo.outputs['Position'], xyz.inputs[0])
            height = math_node(t, 'MULTIPLY', xyz.outputs['Z'], .65)
            height = math_node(t, 'ADD', height, math_node(t, 'MULTIPLY', noise.outputs['Fac'], .055))
            t.links.new(height, emit.inputs['Color'])
        elif model == 'toon':
            diffuse = node(t, 'ShaderNodeBsdfDiffuse')
            diffuse.inputs['Color'].default_value = (1, 1, 1, 1)
            rgb = node(t, 'ShaderNodeShaderToRGB', 'Eevee-only controlled cel shading')
            t.links.new(diffuse.outputs[0], rgb.inputs[0])
            tone = ramp(t, rgb.outputs[0], [(0, tuple(v*.45 for v in base)+(1,)),
                                           (.45, tuple(v*.8 for v in base)+(1,)),
                                           (.8, (*base, 1))])
            t.links.new(tone, emit.inputs['Color'])
        elif model == 'stain':
            # Two scales make a visible stain texture without pretending it is
            # a segmented nucleus or an electron-scattering calculation.
            fine = node(t, 'ShaderNodeTexNoise', 'Fine stain granularity, artistic')
            fine.inputs['Scale'].default_value = 95
            fine.inputs['Detail'].default_value = 2
            t.links.new(tex.outputs['Object'], fine.inputs['Vector'])
            mixed = math_node(t, 'ADD', math_node(t, 'MULTIPLY', noise.outputs['Fac'], .72),
                              math_node(t, 'MULTIPLY', fine.outputs['Fac'], .28))
            tone = ramp(t, mixed, [(.25, tuple(v*.25 for v in base)+(1,)),
                                  (.52, (*base, 1)),
                                  (.74, tuple(min(1., v*1.6+.06) for v in base)+(1,))])
            t.links.new(tone, emit.inputs['Color'])
        else:
            strength = .08 if model == 'paint' else .55
            tone = ramp(t, noise.outputs['Fac'], [(0, tuple(v*(1-strength) for v in base)+(1,)),
                                                 (1, tuple(v*(1+strength) for v in base)+(1,))])
            t.links.new(tone, emit.inputs['Color'])
        emit.inputs['Strength'].default_value = profile.get('emission_strength', .85) if model == 'emission' else 1
        if model == 'emission' and profile.get('emission_depth',0):
            facing=node(t,'ShaderNodeLayerWeight','Artistic shape cue, not a fluorescence optical model')
            depth=profile['emission_depth']
            shade=math_node(t,'ADD',math_node(t,'MULTIPLY',math_node(t,'SUBTRACT',1,facing.outputs['Facing']),depth),1-depth)
            t.links.new(math_node(t,'MULTIPLY',shade,profile.get('emission_strength',.85)),emit.inputs['Strength'])
        t.links.new(emit.outputs[0], out.inputs['Surface'])
        return m
    bsdf = node(t, 'ShaderNodeBsdfPrincipled')
    bsdf.inputs['Base Color'].default_value = (*base, 1)
    bsdf.inputs['Roughness'].default_value = profile.get('roughness', .5)
    bsdf.inputs['Specular IOR Level'].default_value = .2
    if model == 'pbr':
        bsdf.inputs['Subsurface Weight'].default_value = .025
        bsdf.inputs['Subsurface Radius'].default_value = (.05, .03, .02)
    bump_size = .018 if model == 'sem' else profile.get('bump', 0)
    if bump_size:
        bump = node(t, 'ShaderNodeBump')
        bump.inputs['Strength'].default_value = .3
        bump.inputs['Distance'].default_value = bump_size
        t.links.new(noise.outputs['Fac'], bump.inputs['Height'])
        t.links.new(bump.outputs[0], bsdf.inputs['Normal'])
    surface = bsdf.outputs[0]
    if model in {'shell', 'xray', 'sem'}:
        facing = node(t, 'ShaderNodeLayerWeight')
        edge = math_node(t, 'POWER', facing.outputs['Facing'], 3)
        if model == 'sem':
            tone = ramp(t, edge, [(0, '#626669'), (1, '#e5e6e4')])
            t.links.new(tone, bsdf.inputs['Base Color'])
            bsdf.inputs['Specular IOR Level'].default_value = .03
        else:
            trans = node(t, 'ShaderNodeBsdfTransparent')
            mix = node(t, 'ShaderNodeMixShader')
            fac = math_node(t, 'ADD', math_node(t, 'MULTIPLY', edge, .32), .015 if model == 'xray' else .025)
            t.links.new(fac, mix.inputs[0])
            t.links.new(trans.outputs[0], mix.inputs[1])
            if model == 'xray':
                emit = node(t, 'ShaderNodeEmission')
                emit.inputs['Color'].default_value = (*base, 1)
                surface = emit.outputs[0]
            t.links.new(surface, mix.inputs[2])
            surface = mix.outputs[0]
    t.links.new(surface, out.inputs['Surface'])
    return m


def geometry_tree(obj, name):
    tree = bpy.data.node_groups.new(name, 'GeometryNodeTree')
    tree.interface.new_socket(name='Geometry', in_out='INPUT', socket_type='NodeSocketGeometry')
    tree.interface.new_socket(name='Geometry', in_out='OUTPUT', socket_type='NodeSocketGeometry')
    inp = node(tree, 'NodeGroupInput')
    out = node(tree, 'NodeGroupOutput')
    modifier = obj.modifiers.new(name, 'NODES')
    modifier.node_group = tree
    return tree, inp.outputs['Geometry'], out.inputs['Geometry']


def point_object(name, coordinates, group, radii=None):
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(np.asarray(coordinates).tolist(), [], [])
    mesh.update()
    if radii is not None:
        attr = mesh.attributes.new('radius_scene', 'FLOAT', 'POINT')
        attr.data.foreach_set('value', np.asarray(radii, dtype=float))
    obj = bpy.data.objects.new(name, mesh)
    group.objects.link(obj)
    obj['point_count'] = len(coordinates)
    return obj


def instance_spheres(obj, mat, radius=.15, subdivisions=2, loop_frames=96, motion=0):
    t, source, target = geometry_tree(obj, 'GN_Atom_or_punctum_instances')
    convert = node(t, 'GeometryNodeMeshToPoints')
    convert.mode = 'VERTICES'
    t.links.new(source, convert.inputs['Mesh'])
    points = convert.outputs['Points']
    if motion:
        time = node(t, 'GeometryNodeInputSceneTime', 'Periodic demonstration motion; not a diffusion simulation')
        phase = math_node(t, 'MULTIPLY', math_node(t, 'SUBTRACT', time.outputs['Frame'], 1), math.tau/loop_frames)
        xyz = node(t, 'ShaderNodeCombineXYZ')
        t.links.new(math_node(t, 'MULTIPLY', math_node(t, 'SINE', phase), motion), xyz.inputs['X'])
        t.links.new(math_node(t, 'MULTIPLY', math_node(t, 'SUBTRACT', math_node(t, 'COSINE', phase), 1), motion), xyz.inputs['Y'])
        position = node(t, 'GeometryNodeSetPosition')
        t.links.new(points, position.inputs['Geometry'])
        t.links.new(xyz.outputs[0], position.inputs['Offset'])
        points = position.outputs['Geometry']
    sphere = node(t, 'GeometryNodeMeshIcoSphere')
    sphere.inputs['Radius'].default_value = radius
    sphere.inputs['Subdivisions'].default_value = subdivisions
    smooth = node(t, 'GeometryNodeSetShadeSmooth')
    t.links.new(sphere.outputs['Mesh'], smooth.inputs['Geometry'])
    setter = node(t, 'GeometryNodeSetMaterial')
    setter.inputs['Material'].default_value = mat
    t.links.new(smooth.outputs['Geometry'], setter.inputs['Geometry'])
    instance = node(t, 'GeometryNodeInstanceOnPoints')
    t.links.new(points, instance.inputs['Points'])
    t.links.new(setter.outputs['Geometry'], instance.inputs['Instance'])
    if obj.data.attributes.get('radius_scene'):
        attr = node(t, 'GeometryNodeInputNamedAttribute')
        attr.data_type = 'FLOAT'
        attr.inputs['Name'].default_value = 'radius_scene'
        t.links.new(math_node(t, 'DIVIDE', attr.outputs['Attribute'], radius), instance.inputs['Scale'])
    t.links.new(instance.outputs['Instances'], target)
    return t


def curve_tubes(obj, mat, radius, spacing=.045):
    t, source, target = geometry_tree(obj, 'GN_Resampled_fiber_tubes')
    resample = node(t, 'GeometryNodeResampleCurve')
    resample.inputs['Mode'].default_value = 'Length'
    resample.inputs['Length'].default_value = spacing
    t.links.new(source, resample.inputs['Curve'])
    circle = node(t, 'GeometryNodeCurvePrimitiveCircle')
    circle.inputs['Resolution'].default_value = 8
    circle.inputs['Radius'].default_value = radius
    tube = node(t, 'GeometryNodeCurveToMesh')
    t.links.new(resample.outputs['Curve'], tube.inputs['Curve'])
    t.links.new(circle.outputs['Curve'], tube.inputs['Profile Curve'])
    tube.inputs['Fill Caps'].default_value = True
    setter = node(t, 'GeometryNodeSetMaterial')
    setter.inputs['Material'].default_value = mat
    t.links.new(tube.outputs['Mesh'], setter.inputs['Geometry'])
    smooth = node(t, 'GeometryNodeSetShadeSmooth')
    t.links.new(setter.outputs['Geometry'], smooth.inputs['Geometry'])
    t.links.new(smooth.outputs['Geometry'], target)


def density(obj, mat, profile, isosurface=False):
    t, source, target = geometry_tree(obj, 'GN_Illustrative_density_from_atoms')
    convert = node(t, 'GeometryNodeMeshToPoints')
    convert.mode = 'VERTICES'
    t.links.new(source, convert.inputs['Mesh'])
    vol = node(t, 'GeometryNodePointsToVolume', 'Synthetic density; never label this as an experimental map')
    vol.inputs['Resolution Mode'].default_value = 'Size'
    vol.inputs['Voxel Size'].default_value = profile.get('volume_voxel_nm', .06)
    vol.inputs['Radius'].default_value = profile.get('volume_radius_nm', .32)
    t.links.new(convert.outputs['Points'], vol.inputs['Points'])
    result = vol.outputs['Volume']
    if isosurface:
        iso = node(t, 'GeometryNodeVolumeToMesh')
        iso.inputs['Threshold'].default_value = .35
        iso.inputs['Adaptivity'].default_value = .08
        t.links.new(result, iso.inputs['Volume'])
        result = iso.outputs['Mesh']
    setter = node(t, 'GeometryNodeSetMaterial')
    setter.inputs['Material'].default_value = mat
    t.links.new(result, setter.inputs['Geometry'])
    t.links.new(setter.outputs['Geometry'], target)
    return t


def grain_image(width, height, amount):
    # A generated, packed image is self-contained in the project .blend, never
    # copied into the skill. Fixed seed prevents random texture flicker.
    rng = np.random.default_rng(217)
    im = bpy.data.images.new('Generated fixed compositor grain', width, height, alpha=True, float_buffer=True)
    a = np.ones((height, width, 4), dtype=np.float32)
    a[:, :, :3] = (1 + rng.normal(0, amount, (height, width, 1))).clip(.5, 1.5)
    im.colorspace_settings.name = 'Non-Color'
    im.pixels.foreach_set(a.ravel())
    im.update()
    im.pack()
    return im


def compositor(scene, profile, mode=None):
    if not hasattr(scene, 'compositing_node_group'):
        raise RuntimeError('This editable style library was validated on Blender 5.2; adapt the compositor API before using Blender 4.x.')
    mode = mode or profile['compositor']
    t = bpy.data.node_groups.new('COMP_STYLE_'+mode, 'CompositorNodeTree')
    t.use_fake_user = True
    t.interface.new_socket(name='Image', in_out='OUTPUT', socket_type='NodeSocketColor')
    layers = node(t, 'CompositorNodeRLayers')
    over = node(t, 'CompositorNodeAlphaOver', 'Composite premultiplied foreground before grading')
    over.inputs['Background'].default_value = (*linear(profile['background']), 1)
    t.links.new(layers.outputs['Image'], over.inputs['Foreground'])
    image = over.outputs['Image']
    if mode == 'ink':
        # Screen-space normal edges are an appearance effect, not mesh edits.
        edges = node(t, 'CompositorNodeFilter', 'Fine normal discontinuities; inspect at delivery resolution')
        edges.inputs['Type'].default_value = 'Sobel'
        t.links.new(layers.outputs['Normal'], edges.inputs['Image'])
        bw = node(t, 'CompositorNodeRGBToBW')
        t.links.new(edges.outputs[0], bw.inputs[0])
        inner = math_node(t, 'MULTIPLY', bw.outputs[0], profile.get('outline_strength', .2))
        silhouette = node(t, 'CompositorNodeFilter', 'Silhouette kept separate from surface creases')
        silhouette.inputs['Type'].default_value = 'Sobel'
        t.links.new(layers.outputs['Alpha'], silhouette.inputs['Image'])
        edge = math_node(t, 'ADD', inner, math_node(t, 'MULTIPLY', silhouette.outputs[0], .25))
        tone = math_node(t, 'MAXIMUM', math_node(t, 'SUBTRACT', 1, edge), .15)
        multiply = node(t, 'ShaderNodeMixRGB')
        multiply.blend_type = 'MULTIPLY'
        multiply.inputs[0].default_value = 1
        t.links.new(image, multiply.inputs[1])
        t.links.new(tone, multiply.inputs[2])
        image = multiply.outputs[0]
    if mode in {'dic', 'phase'}:
        bw = node(t, 'CompositorNodeRGBToBW', 'Geometric height proxy, not measured optical path')
        t.links.new(image, bw.inputs[0])
        if mode == 'dic':
            shifted = []
            for sign in (-1, 1):
                n = node(t, 'CompositorNodeTransform', 'Opposed shear samples')
                n.inputs['X'].default_value = sign*1.25
                n.inputs['Y'].default_value = sign*.8
                t.links.new(bw.outputs[0], n.inputs['Image'])
                shifted.append(n.outputs[0])
            value = math_node(t, 'ADD', math_node(t, 'MULTIPLY', math_node(t, 'SUBTRACT', *shifted), 5), .38)
        else:
            blur = node(t, 'CompositorNodeBlur', 'Uncalibrated halo kernel')
            blur.inputs['Size'].default_value = (11, 11)
            t.links.new(bw.outputs[0], blur.inputs['Image'])
            halo = math_node(t, 'SUBTRACT', blur.outputs[0], bw.outputs[0])
            value = math_node(t, 'ADD', math_node(t, 'MULTIPLY', halo, 2.8), .32)
        image = value
    if mode in {'soft', 'fluorescence'}:
        blur = node(t, 'CompositorNodeBlur', 'Artistic blur; no calibrated PSF implied')
        size = profile.get('blur_px', .6)
        blur.inputs['Size'].default_value = (size, size)
        t.links.new(image, blur.inputs['Image'])
        image = blur.outputs[0]
    if mode in {'fluorescence', 'filmic'}:
        glow = node(t, 'CompositorNodeGlare', 'Restrained bloom; inspect fine structures')
        glow.inputs['Type'].default_value = 'Fog Glow'
        glow.inputs['Threshold'].default_value = .65 if mode == 'fluorescence' else 1.4
        glow.inputs['Strength'].default_value = .12
        glow.inputs['Size'].default_value = .25
        t.links.new(image, glow.inputs['Image'])
        image = glow.outputs['Image']
    if mode in {'grain', 'dic', 'phase', 'soft'}:
        texture = node(t, 'CompositorNodeImage', 'Fixed grain; regenerated from seed on rebuild')
        texture.image = grain_image(scene.render.resolution_x, scene.render.resolution_y,
                                    profile.get('grain_amount', .035))
        mix = node(t, 'ShaderNodeMixRGB')
        mix.blend_type = 'MULTIPLY'
        mix.inputs[0].default_value = 1
        t.links.new(image, mix.inputs[1])
        t.links.new(texture.outputs['Image'], mix.inputs[2])
        image = mix.outputs[0]
    out = node(t, 'NodeGroupOutput')
    t.links.new(image, out.inputs['Image'])
    scene.compositing_node_group = t
    scene['selected_compositor'] = mode
    return t
