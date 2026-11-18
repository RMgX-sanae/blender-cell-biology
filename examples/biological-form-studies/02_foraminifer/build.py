"""Share one uncut CT-derived shell between cameras and preserve a calibrated virtual section."""

import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Vector

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import studio as st
from views import view
from detail_scene import compositors


def registered_plane(scene, info):
    data = bpy.data.images.load(str(HERE / "registered-slice.png"))
    data.name = "Original 16-bit CT | native Z 192"
    data.colorspace_settings.name = "Non-Color"
    data.pack()
    ny, nx = info["source_shape_zyx"][1:]
    step = info["source_voxel_um"] / info["unit_um"]
    x = nx * step * 0.5
    y = ny * step * 0.5
    z = (info["origin_um_xyz"][2] + info["raw_z_index"] * info["source_voxel_um"]) / 100
    ob = st.mesh(
        "Calibrated CT pixel plane",
        [(-x, -y, z), (x, -y, z), (x, y, z), (-x, y, z)],
        [(0, 1, 2, 3)],
        scene.collection,
        smooth=False,
    )
    uv = ob.data.uv_layers.new(name="Physical pixel orientation")
    for ref, xy in zip(uv.uv, [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ref.vector = xy
    m = bpy.data.materials.new("Window 27300 to 43000")
    m.use_nodes = True
    n = m.node_tree.nodes
    l = m.node_tree.links
    n.clear()
    tex = n.new("ShaderNodeTexImage")
    tex.image = data
    tex.interpolation = "Closest"
    remap = n.new("ShaderNodeMapRange")
    remap.inputs["From Min"].default_value = 27300 / 65535
    remap.inputs["From Max"].default_value = 43000 / 65535
    remap.clamp = True
    emit = n.new("ShaderNodeEmission")
    out = n.new("ShaderNodeOutputMaterial")
    l.new(tex.outputs["Color"], remap.inputs["Value"])
    l.new(remap.outputs[0], emit.inputs["Color"])
    l.new(emit.outputs[0], out.inputs["Surface"])
    ob.data.materials.append(m)
    cam = st.camera(scene, (0, 0, 26), (0, 0, z), y * 3.12)
    # Looking down native +Z with source +Y up; scene and PNG share the voxel-center origin.
    cam.rotation_euler = (0.0, 0.0, 0.0)
    scene.view_settings.view_transform = "Standard"
    scene.world = bpy.data.worlds.new("Black outside acquired field")
    scene.world.color = (0, 0, 0)
    scene.world.use_nodes = True
    scene.world.node_tree.nodes.get("Background").inputs["Strength"].default_value = 0
    compositors(scene, {"background": "000000", "compositor": {"preset": "clean", "contrast": 0}})
    ob["source_z_index"] = info["raw_z_index"]
    ob["pixel_size_um"] = info["source_voxel_um"]


def main():
    a = st.options()
    st.init()
    whole = st.scene("01  Complete calcite test", a)
    pores = st.scene("02  Shell wall", a)
    raw = st.scene("03  Registered virtual section", a)
    info = json.loads((HERE / "inventory.json").read_text())
    data = np.load(HERE / "geometry.npz")
    col = st.collection("Segmented specimen GE_GAD_2", whole.collection)
    pores.collection.children.link(col)
    # The scan supplies pores, surface roughness and chamber silhouettes; add no procedural displacement.
    m = st.material(
        "Calcite | illustrative reflected light",
        "D4C1A0",
        settings=dict(
            roughness=0.40,
            bump=0,
            subsurface_weight=0.070,
            subsurface_radius=[0.67, 0.39, 0.23],
            subsurface_scale=0.065,
            specular=0.33,
        ),
    )
    p = m.node_tree.nodes.get("Principled BSDF")
    p.inputs["Coat Weight"].default_value = 0.10
    shell = st.tri_mesh("Globigerinella adamsi | unsliced test", data["shell_v"], data["shell_f"], col, m)
    shell["specimen_id"] = info["specimen"]
    shell["threshold_raw"] = info["threshold"]
    shell["unit_um"] = 100.0
    d = Vector((0.08, -0.35, 1)).normalized()
    rotation = (-d).to_track_quat("-Z", "Y").to_euler()
    rotation.rotate_axis("Z", -0.68)
    q = rotation.to_quaternion()
    right = np.array(q @ Vector((1, 0, 0)))
    up = np.array(q @ Vector((0, 1, 0)))
    bounds = data["shell_v"] @ np.stack((right, up), axis=1)
    lo, hi = bounds.min(0), bounds.max(0)
    center = (lo + hi) / 2
    span = max(hi[0] - lo[0], 1.50 * (hi[1] - lo[1])) * 1.22
    view(whole, right * center[0] + up * center[1], d, span, ("33454B", "11212A"), roll=-0.68, light=0.48)
    view(pores, (0.55, 0.5, 0), (0.63, 0.56, 0.89), 11.3, ("35424A", "0E1C25"), roll=0.32, light=0.53)
    for sc in (whole, pores):
        sc.view_settings.look = "AgX - Medium High Contrast"
        for ob in sc.objects:
            if ob.type == "LIGHT" and ob.name.startswith("Soft fill"):
                ob.data.energy *= 0.51
    registered_plane(raw, info)
    for sc in (whole, pores, raw):
        sc["source"] = "Siccha et al. (2023), PANGAEA.949585, GE_GAD_2. One scene unit = 100 micrometres."
    st.finish(a, [(whole, "test"), (pores, "wall"), (raw, "ct-section")])


if __name__ == "__main__":
    main()
