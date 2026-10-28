"""Editable teaching prototypes; run with Blender --background --python this.py -- ... ."""
import argparse
import json
import math
import random
import sys
from pathlib import Path

import bpy
from mathutils import Vector

sys.path.insert(0, str(Path(__file__).resolve().parent))
from distributions import patch_points

PALETTE = {
    "membrane": (0.17, 0.54, 0.60, 1), "dna": (0.12, 0.31, 0.61, 1),
    "new_dna": (0.91, 0.46, 0.13, 1), "rna": (0.78, 0.22, 0.38, 1),
    "protein": (0.55, 0.35, 0.65, 1), "motor": (0.92, 0.56, 0.13, 1),
    "track": (0.19, 0.54, 0.39, 1), "nucleus": (0.50, 0.43, 0.65, 1),
    "neutral": (0.60, 0.66, 0.68, 1), "text": (0.055, 0.075, 0.09, 1),
    "mitochondria": (0.78, 0.22, 0.38, 1), "golgi": (0.91, 0.46, 0.13, 1),
}
ARGS = None
MATERIALS = {}
NOTES = []


def tag(obj, compartment="schematic_space", evidence="schematic", role=None):
    obj["evidence"] = evidence
    obj["compartment"] = compartment
    if role:
        obj["biological_role"] = role
    return obj


def material(name):
    if name in MATERIALS:
        return MATERIALS[name]
    mat = bpy.data.materials.new(name)
    mat.diffuse_color = PALETTE[name]
    mat.use_nodes = True
    shader = mat.node_tree.nodes.get("Principled BSDF")
    shader.inputs["Base Color"].default_value = PALETTE[name]
    shader.inputs["Roughness"].default_value = ARGS.roughness
    if name == "text":
        emission = mat.node_tree.nodes.new("ShaderNodeEmission")
        emission.inputs["Color"].default_value = PALETTE[name]
        emission.inputs["Strength"].default_value = 1.0
        mat.node_tree.links.new(emission.outputs[0], mat.node_tree.nodes.get("Material Output").inputs["Surface"])
    MATERIALS[name] = mat
    return mat


def finish(obj, name, mat, compartment="schematic_space"):
    obj.name = name
    if hasattr(obj.data, "materials"):
        obj.data.materials.append(material(mat))
    if obj.type == "MESH":
        for poly in obj.data.polygons:
            poly.use_smooth = True
    return tag(obj, compartment)


def ball(name, location, scale, mat, parent=None):
    bpy.ops.mesh.primitive_uv_sphere_add(segments=20, ring_count=12, location=location)
    obj = finish(bpy.context.object, name, mat)
    obj.scale = (scale,)*3 if isinstance(scale, (int, float)) else scale
    obj.parent = parent
    return obj


def line(name, points, radius, mat, parent=None, cyclic=False):
    curve = bpy.data.curves.new(name, "CURVE")
    curve.dimensions = "3D"
    curve.resolution_u = 12
    curve.bevel_depth = radius
    curve.bevel_resolution = 3
    spline = curve.splines.new("POLY")
    spline.points.add(len(points)-1)
    for point, co in zip(spline.points, points):
        point.co = (*co, 1)
    spline.use_cyclic_u = cyclic
    obj = bpy.data.objects.new(name, curve)
    bpy.context.collection.objects.link(obj)
    obj.parent = parent
    return finish(obj, name, mat)


def empty(name):
    obj = bpy.data.objects.new(name, None)
    bpy.context.collection.objects.link(obj)
    return tag(obj)


def key(obj, path, frame, value):
    setattr(obj, path, value)
    obj.keyframe_insert(data_path=path, frame=max(1, round(frame)))


def visibility(obj, first, last):
    for frame, visible in [(1, first == 1), (max(1, first-1), False), (first, True), (last, True), (last+1, False)]:
        # The duplicate frame at 1 intentionally ends with the correct first state.
        key(obj, "hide_render", frame, not visible)
        key(obj, "hide_viewport", frame, not visible)


def text_label(body, location, size=0.28, name=None):
    bpy.ops.object.text_add(location=location)
    obj = finish(bpy.context.object, name or body, "text")
    obj.data.body = body
    obj.data.size = size
    obj.data.align_x = "CENTER"
    obj.rotation_euler = bpy.context.scene.camera.rotation_euler
    return obj


def camera(location=(0, -20, 11), target=(0, 0, 0), scale=14):
    bpy.ops.object.camera_add(location=location)
    obj = bpy.context.object
    obj.rotation_euler = (Vector(target)-obj.location).to_track_quat("-Z", "Y").to_euler()
    obj.data.type = "ORTHO"
    obj.data.ortho_scale = scale
    obj.data.lens = 50
    bpy.context.scene.camera = obj
    return obj


def shell(name, center, axes, mat="membrane", pinch=0.0):
    # Back half of an ellipsoid: a deliberate cutaway, not an intact transparent membrane.
    vertices, faces = [], []
    nu, nv = 40, 28
    for i in range(nu+1):
        theta = 0.015 + (math.pi-0.03)*i/nu
        x = axes[0]*math.cos(theta)
        contraction = 1-pinch*math.exp(-(x/(axes[0]*0.27))**2)
        for j in range(nv+1):
            phi = math.pi*j/nv
            vertices.append((center[0]+x,
                             center[1]+axes[1]*math.sin(theta)*math.sin(phi)*contraction,
                             center[2]+axes[2]*math.sin(theta)*math.cos(phi)*contraction))
    for i in range(nu):
        for j in range(nv):
            a = i*(nv+1)+j
            faces.append((a, a+1, a+nv+2, a+nv+1))
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)
    return finish(obj, name, mat)


def compositor(style):
    scene = bpy.context.scene
    scene.render.use_compositing = style != "none"
    if style == "none":
        return
    # Blender 5 uses a compositing node group; Blender 4 uses scene.node_tree.
    if hasattr(scene, "compositing_node_group"):
        tree = bpy.data.node_groups.new("Biology_Compositor", "CompositorNodeTree")
        tree.interface.new_socket(name="Image", in_out="OUTPUT", socket_type="NodeSocketColor")
        scene.compositing_node_group = tree
        output = tree.nodes.new("NodeGroupOutput")
    else:
        scene.use_nodes = True
        tree = scene.node_tree
        tree.nodes.clear()
        output = tree.nodes.new("CompositorNodeComposite")
    render = tree.nodes.new("CompositorNodeRLayers")
    render.location = (-400, 0)
    source = render.outputs["Image"]
    if style == "contrast":
        grade = tree.nodes.new("CompositorNodeBrightContrast")
        grade.inputs["Contrast"].default_value = 8.0
        tree.links.new(source, grade.inputs["Image"])
        source = grade.outputs["Image"]
    elif style == "glow":
        glow = tree.nodes.new("CompositorNodeGlare")
        if hasattr(glow, "glare_type"):
            glow.glare_type = "FOG_GLOW"
        else:
            glow.inputs["Type"].default_value = "Fog Glow"
        if hasattr(glow, "threshold"):
            glow.threshold = 1.0
        tree.links.new(source, glow.inputs["Image"])
        source = glow.outputs["Image"]
        NOTES.append("Glow is a compositing effect, not measured fluorescence or protein abundance.")
    output.location = (220, 0)
    if ARGS.style == "teaching" and not ARGS.transparent:
        over = tree.nodes.new("CompositorNodeAlphaOver")
        if "Background" in over.inputs:
            over.inputs["Factor"].default_value = 1.0
            over.inputs["Background"].default_value = (*ARGS.background, 1)
            tree.links.new(source, over.inputs["Foreground"])
        else:
            over.inputs[0].default_value = 1.0
            over.inputs[1].default_value = (*ARGS.background, 1)
            tree.links.new(source, over.inputs[2])
        source = over.outputs["Image"]
    tree.links.new(source, output.inputs["Image"])


def setup():
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete(use_global=False)
    scene = bpy.context.scene
    if ARGS.engine == "cycles":
        scene.render.engine = "CYCLES"
        scene.cycles.samples = 24
        scene.cycles.use_denoising = True
    else:
        for engine in ("BLENDER_EEVEE_NEXT", "BLENDER_EEVEE"):
            try:
                scene.render.engine = engine
                break
            except TypeError:
                pass
        else:
            raise RuntimeError("No supported Eevee engine; choose --engine cycles explicitly")
    scene.render.resolution_x = ARGS.width
    scene.render.resolution_y = ARGS.height
    scene.render.resolution_percentage = 100
    scene.render.image_settings.file_format = "PNG"
    scene.render.image_settings.color_mode = "RGBA"
    scene.render.fps = ARGS.fps
    scene.frame_start, scene.frame_end = 1, ARGS.frames
    scene.world.use_nodes = True
    bg = scene.world.node_tree.nodes.get("Background")
    bg.inputs["Color"].default_value = (*ARGS.background, 1)
    bg.inputs["Strength"].default_value = 0.35
    scene.view_settings.view_transform = "Standard" if ARGS.style == "teaching" else "AgX"
    scene.render.film_transparent = ARGS.transparent or (ARGS.style == "teaching" and ARGS.compositor != "none")
    camera()
    for name, pos, power, size in [("Key", (-5,-7,10),500,7),("Fill",(6,-2,4),220,6),("Rim",(0,5,7),350,5)]:
        bpy.ops.object.light_add(type="AREA", location=pos)
        light = bpy.context.object
        light.name = name
        light.data.energy, light.data.shape, light.data.size = power, "DISK", size
        light.rotation_euler = (-light.location).to_track_quat("-Z", "Y").to_euler()
    compositor(ARGS.compositor)
    scene["evidence"] = "schematic except imported coordinate objects explicitly tagged measured"
    scene["scale_mode"] = "schematic; no physical scale bar"
    scene["time_mode"] = "schematic; playback seconds are not biological seconds"


def make_cell():
    shell("Plasma_membrane_cutaway", (0,0,0), (5,3.4,3.4))
    shell("Nuclear_envelope_cutaway", (-1.5,0,0.4), (1.7,1.5,1.6), "nucleus")
    ball("Nucleolus_membraneless", (-1.6,0,0.6), 0.52, "protein")
    for i in range(5):
        z = -0.7 + i*0.35
        line(f"ER_cisterna_{i}", [(-3,0.8,z),(-2.4,1.3,z),(-1,1.8,z),(0.3,1.5,z)],0.12,"membrane")
    for i in range(4):
        line(f"Golgi_cisterna_{i}", [(1.1,0.2,0.2+i*.3),(2,0.45,0.5+i*.3),(2.8,0.25,0.2+i*.3)],.15,"golgi")
    for i, pos in enumerate([(2.8,0,-1.3),(-2.4,-.3,-1.6),(2.5,.5,2)]):
        mito = ball(f"Mitochondrion_{i}",pos,(.82,.4,.4),"mitochondria")
        mito["detail"] = "overview proxy; cristae and membrane topology omitted"
    for i, pos in enumerate([(3.8,.3,.2),(.6,-.2,-1.8),(-.2,.3,2.5)]):
        ball(f"Endolysosome_proxy_{i}", pos, .28, "motor")
    text_label("Animal cell | cutaway schematic", (0,0,4.8), .34)
    text_label("Nucleus", (-2.8,-1.6,1.9))
    text_label("Golgi", (2.2,-1.4,1.9))
    text_label("Organelle abundance and sizes are illustrative",(0,-1.2,-3.8),.24)
    NOTES.extend(["Generic animal interphase overview; not a particular cell reconstruction.",
                  "ER/Golgi are overview proxies. No claim of continuous lumen in this low-detail starter."])


def make_receptors():
    camera((0,-13,12), scale=14)
    # XY membrane midplane, extracellular +Z, cytosolic -Z.
    for z in (-.18,.18):
        bpy.ops.mesh.primitive_cube_add(size=1,location=(0,0,z))
        obj = finish(bpy.context.object, f"Membrane_leaflet_{z}", "membrane")
        obj.scale = (11,6.8,.09)
    points = patch_points(ARGS.count,10.2,6.0,ARGS.distribution,.62,ARGS.seed)
    for i,(x,y) in enumerate(points):
        root = empty(f"Receptor_{i:03d}")
        root.location=(x,y,0)
        root["position_model"]=ARGS.distribution
        root["outward_normal"]=[0.0,0.0,1.0]
        line(f"TM_{i}",[(0,0,-.35),(0,0,.45)],.10,"protein",root)
        ball(f"Ectodomain_{i}",(0,0,.65),(.23,.19,.31),"protein",root)
        line(f"Cytosolic_tail_{i}",[(0,0,-.36),(.1,0,-.65)],.055,"protein",root)
    text_label(f"{ARGS.distribution.capitalize()} receptor positions | hypothesis",(0,0,4.2),.32)
    text_label("Extracellular +Z",(-3.8,1.2,1.8),.25)
    text_label("Receptor shape, count and density are schematic",(0,-4,-.6),.24)
    NOTES.append("Single-pass generic receptors, not GPCRs. Exact count and footprint exclusion are preserved.")


def make_vesicle():
    camera((0,-18,8),target=(0,0,.5),scale=13)
    for z in (-.24,.24):
        line("Microtubule_outline",[(-5.5,0,z),(5.5,0,z)],.11,"track")
    for i in range(37):
        x=-5.4+i*.3
        line(f"Tubulin_repeat_proxy_{i}",[(x,0,-.24),(x,0,.24)],.07,"track")
    root=empty("Motor_vesicle_transport")
    v=shell("Vesicle_cutaway",(0,0,2.1),(1.1,.9,.9))
    v.parent=root
    for i,p in enumerate([(-.3,-.1,2),(.3,.1,2.3),(0,-.1,2.5)]):
        b=ball(f"Luminal_cargo_{i}",p,.17,"rna",root)
        b["compartment"]="vesicle_lumen"
    line("Motor_stalk",[(0,0,.34),(0,0,1.3)],.09,"motor",root)
    for i,x in enumerate((-.18,.18)):
        ball(f"Motor_head_{i}",(x,0,.37),(.18,.18,.14),"motor",root)
    key(root,"location",1,(-3.6,0,0))
    key(root,"location",ARGS.frames,(3.6,0,0))
    text_label("-",(-5.8,0,.6),.5)
    text_label("+",(5.8,0,.6),.5)
    text_label("Plus-end-directed motor transport",(0,0,4.0),.32)
    text_label("Cargo remains inside the vesicle | time not to scale",(0,-1,-1.2),.25)
    NOTES.extend(["Generic plus-end-directed kinesin-like transport; kinesin-family direction is not universal.",
                  "Microtubule is a schematic track, not its 13-protofilament structure. No motor stepping, budding or fusion."])


def make_mitosis():
    camera((0,-23,6),scale=15)
    total=ARGS.frames
    envelope=shell("Mother_cell_cutaway",(0,0,0),(5.6,2.2,3.1))
    envelope.shape_key_add(name="Basis")
    pinched=shell("Pinch_template",(0,0,0),(5.6,2.2,3.1),pinch=.98)
    shape=envelope.shape_key_add(name="Furrow")
    for dst,src in zip(shape.data,pinched.data.vertices):
        dst.co=src.co
    bpy.data.objects.remove(pinched,do_unlink=True)
    key(shape,"value",total*.5,0)
    key(shape,"value",total*.9,1)
    cut=round(total*.94)
    visibility(envelope,1,cut-1)
    for side in (-1,1):
        daughter=shell(f"Daughter_{side}",(side*2.9,0,0),(2.65,2.15,2.8))
        visibility(daughter,cut,total)
        nucleus=shell(f"Reformed_nucleus_{side}",(side*2.8,0,0),(1.5,1.2,1.8),"nucleus")
        visibility(nucleus,cut,total)
        for i,z in enumerate((-1.2,-.4,.4,1.2)):
            chrom=empty(f"Chromatid_{side}_{i}")
            chrom["homolog_pair"]=i//2
            chrom["chromosome_id"]=i
            # One sister chromatid each; paired sisters begin adjacent at metaphase.
            line(f"Chromatid_mesh_{side}_{i}",[(side*.3,0,z-.3),(0,0,z),(side*.3,0,z+.3)],.10,"dna",chrom)
            key(chrom,"location",1,(side*.18,0,0))
            key(chrom,"location",total*.35,(side*.18,0,0))
            key(chrom,"location",total*.72,(side*2.8,0,0))
            key(chrom,"location",total,(side*2.8,0,0))
            # Lines ending at kinetochores: animate endpoints via shape keys.
            mt=line(f"Kinetochore_MT_{side}_{i}",[(side*4.4,.12,0),(side*.18,.12,z)],.025,"track")
            # Curve point co keyframes preserve an attached endpoint throughout separation.
            point=mt.data.splines[0].points[-1]
            for f,x in [(1,side*.18),(round(total*.35),side*.18),(round(total*.72),side*2.8)]:
                point.co=(x,.12,z,1)
                point.keyframe_insert(data_path="co",frame=f)
            visibility(mt,1,cut-1)
        ball(f"Spindle_pole_{side}",(side*4.4,.12,0),.20,"motor")
    text_label("Metaphase > anaphase > cytokinesis",(0,-.8,4.1),.33)
    text_label("2n = 4 teaching model | sister chromatids separate",(0,-.8,-3.8),.26)
    NOTES.extend(["Starts after DNA replication and nuclear-envelope breakdown; not the full cell cycle.",
                  "Four replicated chromosomes at metaphase; four single-chromatid chromosomes per daughter.",
                  "Furrow uses a cutaway deformation; final abscission is a labeled topology switch, not a membrane simulation."])


def make_replication():
    camera((0,-20,6),scale=14)
    # Fork moves +X. Upper parental template is 3'->5' along +X.
    for z in (-.18,.18):
        line("Parental_duplex_ahead",[(1,0,z),(5.3,0,z)],.06,"dna")
    for side in (-1,1):
        line(f"Parental_template_{side}",[(-5,0,side*1.45),(0,0,side*1.45),(1,0,side*.18)],.07,"dna")
    lead=line("Leading_nascent_DNA_5_to_3",[(-4.9,-.1,1.15),(.05,-.1,1.15)],.07,"new_dna")
    key(lead.data,"bevel_factor_end",1,.06)
    key(lead.data,"bevel_factor_end",ARGS.frames,1)
    for i in range(4):
        right=-3.3+i*1.05
        # Each fragment grows toward -X: its 3' end is on the left.
        fragment=line(f"Okazaki_5_to_3_{i}",[(right,-.1,-1.15),(right-.8,-.1,-1.15)],.07,"new_dna")
        start=1+round(i*ARGS.frames*.18)
        key(fragment.data,"bevel_factor_end",1,0)
        key(fragment.data,"bevel_factor_end",start,0)
        key(fragment.data,"bevel_factor_end",min(ARGS.frames,start+ARGS.frames*.2),1)
        ball(f"Primer_5prime_{i}",(right,-.1,-1.15),.09,"rna")
    ball("Helicase_proxy",(.9,0,0),(.45,.32,.54),"motor")
    text_label("Leading: 5' > 3'",(-2.4,-.2,2.1),.30)
    text_label("Lagging fragments: 3' < 5'",(-2.4,-.2,-2.1),.30)
    text_label("Parental DNA",(3.2,-.2,1),.28)
    text_label("Fork-local schematic | both new strands grow 5' to 3'",(0,-1,3.6),.29)
    text_label("Blue: parental   Orange: new DNA   Rose: primers",(0,-.8,-3),.24)
    NOTES.extend(["Fork-local stationary schematic, not moving replication-fork kinematics.",
                  "Leading strand grows +X, lagging fragments -X. No polymerase stoichiometry, primer removal or ligation animation."])


def make_transcription():
    camera((0,-20,7),scale=14)
    for z in (-.22,.22):
        line("DNA_duplex_overview",[(-5,0,z),(5,0,z)],.08,"dna")
    polymerase=ball("RNA_polymerase_proxy",(-3,0,0),(.65,.48,.65),"protein")
    key(polymerase,"location",1,(-3,0,0))
    key(polymerase,"location",ARGS.frames,(3,0,0))
    # RNA 3' stays at the catalytic proxy; the emerged 5' end trails away.
    points=[(-3-.1*i,-.18,-.3-.07*i) for i in range(30)]
    rna=line("Nascent_RNA_5prime_to_3prime",list(reversed(points)),.07,"rna")
    for frame,px in [(1,-3),(ARGS.frames,3)]:
        length=.4+3.2*(frame-1)/(ARGS.frames-1)
        for i,p in enumerate(rna.data.splines[0].points):
            t=i/29
            p.co=(px-length*(1-t),-.18,-.3-1.5*(1-t),1)
            p.keyframe_insert(data_path="co",frame=frame)
    text_label("DNA template: 3' > 5'",(0,-.2,1.8),.31)
    text_label("RNA synthesis: 5' > 3'",(0,-.2,-2.7),.31)
    text_label("3' end remains at polymerase; 5' end emerges",(0,-.3,3),.28)
    NOTES.extend(["Linear template overview; transcription bubble and RNA-DNA hybrid are not geometrically resolved.",
                  "RNA is synthesized at its 3' end; RNA motion is a schematic kinematic construction."])


def make_translation():
    camera((0,-20,7),scale=14)
    line("mRNA_5_to_3",[(-5,0,0),(5,0,0)],.08,"rna")
    root=empty("Ribosome_translation")
    ball("Large_subunit",(0,0,.55),(.72,.53,.62),"protein",root)
    ball("Small_subunit",(0,0,-.32),(.65,.43,.3),"neutral",root)
    peptide=line("Nascent_peptide_N_to_C",[(0,0,2.8),(.2,0,2.4),(-.1,0,1.15),(0,0,.55)],.08,"new_dna",root)
    # C-terminal attachment is inside the ribosome; the visible N-terminal segment emerges.
    for frame,tip in [(1,1.3),(ARGS.frames,3.0)]:
        for i,point in enumerate(peptide.data.splines[0].points):
            t=i/3
            point.co=(.15*math.sin(i*2),-.1,tip*(1-t)+.55*t,1)
            point.keyframe_insert(data_path="co",frame=frame)
    key(root,"location",1,(-3,0,0))
    key(root,"location",ARGS.frames,(3,0,0))
    text_label("5'",(-5.5,0,.6),.36)
    text_label("3'",(5.5,0,.6),.36)
    text_label("Ribosome moves 5' > 3' on mRNA",(0,-.3,-1.5),.31)
    text_label("Peptide: N terminus emerges first; C terminus grows",(0,-.3,4.1),.27)
    NOTES.extend(["Elongation overview only; tRNA A/P/E cycle, codon decoding and termination require the detailed recipe.",
                  "C-terminal attachment is inside the ribosome proxy, not at the external tunnel exit; exact tunnel geometry is omitted."])


def make_structure():
    path=Path(__file__).resolve().parents[1]/"assets/1BNA.pdb"
    atoms=[]
    for record in path.read_text().splitlines():
        if record.startswith("ATOM  ") and record[16:17] in (" ","A"):
            atoms.append((record[21],tuple(float(record[k:k+8]) for k in (30,38,46)),record[76:78].strip()))
    if not atoms:
        raise ValueError("Bundled 1BNA coordinates contain no atoms")
    center=Vector([sum(a[1][i] for a in atoms)/len(atoms) for i in range(3)])
    for i,(chain,co,element) in enumerate(atoms):
        pos=(Vector(co)-center)*.13
        atom=ball(f"1BNA_atom_{i}",pos,.13*(1.8 if element=="P" else 1.5),"dna" if chain=="A" else "new_dna")
        atom["evidence"]="measured"
        atom["source"]="PDB 1BNA experimental coordinate model; atoms only; no waters"
        atom["source_chain"]=chain
    camera((0,-15,8),scale=10)
    bpy.context.scene["scale_mode"]="uniform atom-center scale; illustrative atom radii"
    bpy.context.scene["meters_per_blender_unit"]=1e-10/.13
    text_label("B-DNA | PDB 1BNA",(0,-.5,3.8),.34)
    text_label("Experimental coordinates; illustrative atom radii",(0,-1,-3.6),.24)
    NOTES.append("1 Blender unit = 1/0.13 angstrom for atom centers; atom radii are illustrative, no bonds inferred.")


def linearize():
    # Blender 4 legacy and Blender 4.4+/5 layered actions.
    for action in bpy.data.actions:
        curves=[]
        if hasattr(action,"fcurves"):
            curves.extend(action.fcurves)
        for layer in getattr(action,"layers",[]):
            for strip in layer.strips:
                for bag in getattr(strip,"channelbags",[]):
                    curves.extend(bag.fcurves)
        for curve in curves:
            for k in curve.keyframe_points:
                k.interpolation="CONSTANT" if "hide_" in curve.data_path else "LINEAR"


def main():
    global ARGS
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scene",choices=("cell","receptors","vesicle","mitosis","replication","transcription","translation","structure"),required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--style",choices=("teaching","realistic"),default="teaching")
    parser.add_argument("--style-profile",type=Path)
    parser.add_argument("--engine",choices=("auto","eevee","cycles"),default="auto")
    parser.add_argument("--compositor",choices=("none","clean","contrast","glow"),default="clean")
    parser.add_argument("--distribution",choices=("uniform","clustered","polarized"),default="clustered")
    parser.add_argument("--count",type=int,default=32)
    parser.add_argument("--seed",type=int,default=23)
    parser.add_argument("--frames",type=int,default=96)
    parser.add_argument("--fps",type=int,default=24)
    parser.add_argument("--width",type=int,default=960)
    parser.add_argument("--height",type=int,default=720)
    parser.add_argument("--transparent",action="store_true")
    parser.add_argument("--preview",action="store_true")
    parser.add_argument("--animate",action="store_true")
    ARGS=parser.parse_args(sys.argv[sys.argv.index("--")+1:] if "--" in sys.argv else [])
    if ARGS.frames<12 or ARGS.fps<1 or min(ARGS.width,ARGS.height)<64 or ARGS.width%2 or ARGS.height%2:
        parser.error("Use >=12 frames, positive fps and even resolution >=64")
    if ARGS.output.exists() and any(ARGS.output.iterdir()):
        parser.error("Output directory must be new or empty; do not overwrite existing projects")
    ARGS.background=(1,1,1) if ARGS.style=="teaching" else (.035,.05,.065)
    ARGS.roughness=.55 if ARGS.style=="teaching" else .30
    if ARGS.style=="realistic":
        PALETTE["text"]=(.88,.91,.94,1)
    profile={}
    if ARGS.style_profile:
        profile=json.loads(ARGS.style_profile.read_text())
        unknown=set(profile)-{"reference","observations","inferences","palette","background","roughness","compositor"}
        if unknown:
            parser.error(f"Unsupported profile keys: {sorted(unknown)}")
        for name,color in profile.get("palette",{}).items():
            if name not in PALETTE or len(color)!=4 or not all(0<=v<=1 for v in color):
                parser.error("Palette values must be known semantic keys and four linear RGBA values in [0,1]")
            PALETTE[name]=tuple(color)
        ARGS.background=tuple(profile.get("background",ARGS.background))
        ARGS.roughness=profile.get("roughness",ARGS.roughness)
        if len(ARGS.background)!=3 or not all(0<=x<=1 for x in ARGS.background) or not 0<=ARGS.roughness<=1:
            parser.error("Profile background/roughness outside [0,1]")
        if "compositor" in profile:
            if profile["compositor"] not in ("none","clean","contrast","glow"):
                parser.error("Unknown compositor in profile")
            ARGS.compositor=profile["compositor"]
    ARGS.output.mkdir(parents=True,exist_ok=True)
    random.seed(ARGS.seed)
    setup()
    globals()["make_"+ARGS.scene]()
    linearize()
    scene=bpy.context.scene
    scene.frame_set(1)
    scene.render.filepath=str(ARGS.output/"frames/frame_")
    bpy.ops.wm.save_as_mainfile(filepath=str(ARGS.output/"scene.blend"))
    manifest={"scene":ARGS.scene,"blender_version":bpy.app.version_string,"parameters":{k:str(v) if isinstance(v,Path) else v for k,v in vars(ARGS).items()},
              "style_profile":profile,"object_count":len(scene.objects),"simplifications":NOTES,
              "scale_mode":scene["scale_mode"],"time_mode":scene["time_mode"],
              "meters_per_blender_unit":scene.get("meters_per_blender_unit"),
              "evidence":"teaching prototypes; see per-object tags", "generator":str(Path(__file__).resolve()),
              "validation":"Generated, not a biological simulation; inspect actual renders separately."}
    (ARGS.output/"manifest.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+"\n")
    if ARGS.preview:
        for frame in sorted({1,round(ARGS.frames/2),ARGS.frames}):
            scene.frame_set(frame)
            scene.render.filepath=str(ARGS.output/f"preview_{frame:04d}.png")
            bpy.ops.render.render(write_still=True)
    if ARGS.animate:
        (ARGS.output/"frames").mkdir(exist_ok=True)
        scene.render.filepath=str(ARGS.output/"frames/frame_")
        bpy.ops.render.render(animation=True)
    print("CELL_BIOLOGY_BUILD_OK",str(ARGS.output))


if __name__=="__main__":
    main()
