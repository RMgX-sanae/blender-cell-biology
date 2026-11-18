"""Small Blender studio used by the three independent projects."""

import argparse
import os
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Vector

SKILL = Path(os.environ.get("CELL_BIOLOGY_SKILL", Path.home() / ".codex/skills/blender-cell-biology"))
sys.path.insert(0, str(SKILL / "scripts"))
from detail_scene import linear, material


def options():
    p = argparse.ArgumentParser()
    p.add_argument("--draft", action="store_true")
    p.add_argument("--width", type=int, default=1920)
    p.add_argument("--samples", type=int, default=56)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--save-only", action="store_true")
    return p.parse_args(sys.argv[sys.argv.index("--") + 1 :] if "--" in sys.argv else [])


def init():
    bpy.ops.wm.read_factory_settings(use_empty=True)
    bpy.context.preferences.filepaths.save_version = 0


def collection(name, parent):
    c = bpy.data.collections.new(name)
    parent.children.link(c)
    return c


def group(name, col):
    o = bpy.data.objects.new(name, None)
    col.objects.link(o)
    o.empty_display_size = 0.45
    return o


def tri_mesh(name, v, f, col, mat=None):
    v = np.ascontiguousarray(v, dtype=np.float32)
    f = np.ascontiguousarray(f, dtype=np.int32)
    data = bpy.data.meshes.new(name)
    data.vertices.add(len(v))
    data.vertices.foreach_set("co", v.ravel())
    data.loops.add(f.size)
    data.loops.foreach_set("vertex_index", f.ravel())
    data.polygons.add(len(f))
    data.polygons.foreach_set("loop_start", np.arange(0, f.size, 3, dtype="i"))
    data.polygons.foreach_set("use_smooth", np.ones(len(f), dtype=bool))
    data.update(calc_edges=True)
    obj = bpy.data.objects.new(name, data)
    col.objects.link(obj)
    if mat:
        data.materials.append(mat)
    return obj


def mesh(name, v, f, col, mat=None, smooth=True):
    d = bpy.data.meshes.new(name)
    d.from_pydata(v, [], f)
    d.update()
    if smooth:
        d.polygons.foreach_set("use_smooth", np.ones(len(d.polygons), dtype=bool))
    ob = bpy.data.objects.new(name, d)
    col.objects.link(ob)
    if mat:
        d.materials.append(mat)
    return ob


def shader(name, color, rough=0.40, metal=0, coat=0, emission=0, **opts):
    mat = material(name, color, settings=dict(roughness=rough, bump=0, subsurface_weight=0, **opts))
    p = mat.node_tree.nodes.get("Principled BSDF")
    p.inputs["Metallic"].default_value = metal
    p.inputs["Coat Weight"].default_value = coat
    if emission:
        p.inputs["Emission Color"].default_value = (*linear(color), 1)
        p.inputs["Emission Strength"].default_value = emission
    return mat


def scene(name, a):
    s = (
        bpy.data.scenes[0]
        if len(bpy.data.scenes) == 1 and bpy.data.scenes[0].name == "Scene"
        else bpy.data.scenes.new(name)
    )
    s.name = name
    bpy.context.window.scene = s
    s.render.engine = "CYCLES"
    s.cycles.samples = 8 if a.draft else a.samples
    s.cycles.use_denoising = True
    s.cycles.use_adaptive_sampling = True
    s.cycles.adaptive_threshold = 0.10 if a.draft else 0.032
    s.cycles.max_bounces = 5
    s.cycles.transparent_max_bounces = 6
    try:
        pref = bpy.context.preferences.addons["cycles"].preferences
        pref.compute_device_type = "METAL"
        pref.refresh_devices()
        devices = [d for d in pref.devices if d.type == "METAL"]
        if devices and os.environ.get("SHOWCASE_CPU") != "1":
            for dev in pref.devices:
                dev.use = dev.type == "METAL"
            s.cycles.device = "GPU"
    except (RuntimeError, TypeError):
        pass
    w = 800 if a.draft else a.width
    s.render.resolution_x = w
    s.render.resolution_y = round(w * 2 / 3)
    s.render.resolution_percentage = 100
    s.render.image_settings.file_format = "PNG"
    s.render.image_settings.color_mode = "RGB"
    s.view_settings.view_transform = "AgX"
    s.render.film_transparent = False
    return s


def camera(s, location, target, span, perspective=False, fstop=None):
    d = bpy.data.cameras.new(s.name + " lens")
    ob = bpy.data.objects.new("Camera " + s.name, d)
    s.collection.objects.link(ob)
    ob.location = Vector(location)
    v = Vector(target) - ob.location
    ob.rotation_euler = v.to_track_quat("-Z", "Y").to_euler()
    if perspective:
        d.type = "PERSP"
        d.lens = d.sensor_width * v.length / span
    else:
        d.type = "ORTHO"
        d.ortho_scale = span
    d.clip_start = 0.01
    d.clip_end = 3000
    if fstop:
        d.dof.use_dof = True
        d.dof.focus_distance = v.length
        d.dof.aperture_fstop = fstop
        d.dof.aperture_blades = 8
    s.camera = ob
    return ob


def studio(s, cam, target, span, palette, world_color="698697", intensity=1.0):
    span = float(span)
    w = bpy.data.worlds.new(s.name + " world")
    s.world = w
    w.use_nodes = True
    bg = w.node_tree.nodes.get("Background")
    bg.inputs["Color"].default_value = (*linear(world_color), 1)
    bg.inputs["Strength"].default_value = 0.35
    q = cam.rotation_euler.to_quaternion()
    right = q @ Vector((1, 0, 0))
    up = q @ Vector((0, 1, 0))
    toward = q @ Vector((0, 0, 1))
    target = Vector(target)
    lights = collection(s.name + " | Softboxes", s.collection)
    for name, xyz, coef, size, color in [
        ("Key", (-0.50, 0.78, 0.70), 55, 0.64, "FFF1DC"),
        ("Soft fill", (0.75, 0.18, 0.38), 21, 0.65, "C5E7F5"),
        ("Rim", (0.27, 0.55, -0.57), 82, 0.48, "C5E5F7"),
        ("Warm edge", (-0.68, -0.05, -0.29), 27, 0.36, "FFD5AF"),
    ]:
        ld = bpy.data.lights.new(name, "AREA")
        ld.energy = coef * span * span * intensity
        ld.shape = "DISK"
        ld.size = size * span
        ld.color = linear(color)
        ob = bpy.data.objects.new(name, ld)
        lights.objects.link(ob)
        ob.location = target + span * (xyz[0] * right + xyz[1] * up + xyz[2] * toward)
        ob.rotation_euler = (target - ob.location).to_track_quat("-Z", "Y").to_euler()
    # A camera-facing color field gives the background a quiet center-to-edge gradient.
    c = collection(s.name + " | Backdrop", s.collection)
    back = target - toward * (span * 2.8)
    sx, sy = span * 3.25, span * 2.3
    pos = [back + dx * sx * right + dy * sy * up for dx, dy in [(-1, -1), (1, -1), (1, 1), (-1, 1)]]
    plane = mesh("Studio gradient", pos, [(0, 1, 2, 3)], c, smooth=False)
    uv = plane.data.uv_layers.new()
    for i, pt in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)]):
        uv.data[i].uv = pt
    m = bpy.data.materials.new("Background " + s.name)
    m.use_nodes = True
    n = m.node_tree.nodes
    l = m.node_tree.links
    n.clear()
    tc = n.new("ShaderNodeTexCoord")
    dist = n.new("ShaderNodeVectorMath")
    dist.operation = "DISTANCE"
    dist.inputs[1].default_value = (0.466, 0.535, 0)
    l.new(tc.outputs["UV"], dist.inputs[0])
    ramp = n.new("ShaderNodeValToRGB")
    ramp.color_ramp.elements[0].color = (*linear(palette[0]), 1)
    ramp.color_ramp.elements[1].color = (*linear(palette[1]), 1)
    ramp.color_ramp.elements[1].position = 0.24
    l.new(dist.outputs["Value"], ramp.inputs[0])
    em = n.new("ShaderNodeEmission")
    l.new(ramp.outputs[0], em.inputs[0])
    out = n.new("ShaderNodeOutputMaterial")
    l.new(em.outputs[0], out.inputs["Surface"])
    plane.data.materials.append(m)
    for attr in ("visible_diffuse", "visible_glossy", "visible_shadow", "visible_transmission"):
        setattr(plane, attr, False)


def finish(a, scenes):
    """Save with both cameras, then render each still. Both use ordinary Blender scenes."""
    a.out.mkdir(parents=True, exist_ok=True)
    for s, stem in scenes:
        s.render.filepath = "//" + stem + ".png"
        s["delivery_stem"] = stem
    bpy.context.window.scene = scenes[0][0]
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type == "VIEW_3D":
                area.spaces.active.region_3d.view_perspective = "CAMERA"
                area.spaces.active.overlay.show_overlays = False
    if not a.draft:
        bpy.ops.wm.save_as_mainfile(filepath=str(a.out / "model.blend"), compress=True)
    if a.save_only:
        return
    for s, stem in scenes:
        bpy.context.window.scene = s
        s.render.filepath = str(a.out / (stem + ".png"))
        bpy.ops.render.render(write_still=True, scene=s.name)
    print("Rendered", [(s.name, stem) for s, stem in scenes])
