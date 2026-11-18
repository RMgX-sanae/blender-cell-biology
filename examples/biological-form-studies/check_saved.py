"""Reopen with automatic execution disabled, check source registration and render the saved cameras."""

import argparse
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Vector

p = argparse.ArgumentParser()
p.add_argument("--out", type=Path, required=True)
p.add_argument("--render", action="store_true")
a = p.parse_args(sys.argv[sys.argv.index("--") + 1 :])
out = a.out.resolve()
scenes = {s["delivery_stem"]: s for s in bpy.data.scenes if "delivery_stem" in s}
assert scenes and all(s.camera and s.compositing_node_group for s in scenes.values())
assert all(
    any(n.type == "R_LAYERS" and n.scene in (None, s) for n in s.compositing_node_group.nodes)
    for s in scenes.values()
)
result = {"scenes": list(scenes), "version": bpy.app.version_string}


def coords(mesh):
    data = np.empty(len(mesh.vertices) * 3, dtype="f4")
    mesh.vertices.foreach_get("co", data)
    return data.reshape(-1, 3)


if "branch" in scenes:
    assert set(scenes) == {"assembly", "rna-core", "branch"}
    info = json.loads((out / "inventory.json").read_text())
    arr = np.load(out / "geometry.npz")
    raw = next(o for o in scenes["assembly"].objects if "All deposited heavy atoms" in o.name)
    ids = np.empty(len(raw.data.vertices), "i4")
    raw.data.attributes["atom_site_id"].data.foreach_get("value", ids)
    assert np.array_equal(ids, arr["atom_site_ids"]) and np.array_equal(coords(raw.data), arr["xyz_nm"])
    groups = [g for g in info["groups"] if "v_" + g["label"] in arr]
    for g in groups:
        obj = next(
            o
            for o in scenes["assembly"].objects
            if o.get("label_asym_id") == g["label"] and o.get("entity_id")
        )
        assert np.array_equal(coords(obj.data), arr["v_" + g["label"]])
    pair = info["substrand_link"]
    measured = float(np.linalg.norm(arr["xyz_nm"][pair[0]] - arr["xyz_nm"][pair[1]]) * 10)
    for s in scenes.values():
        bond = next(o for o in s.objects if "endpoint_indices" in o)
        assert list(bond["endpoint_indices"]) == pair and abs(bond["distance_angstrom"] - measured) < 1e-5
    bpy.context.window.scene = scenes["branch"]
    bpy.context.view_layer.update()
    annotation = next(o for o in scenes["branch"].objects if o.type == "FONT")
    assert annotation.data.body == "2′–5′" and annotation.data.font.packed_file
    corners = [
        world_to_camera_view(scenes["branch"], scenes["branch"].camera, annotation.matrix_world @ Vector(v))
        for v in annotation.bound_box
    ]
    assert all(0.025 < c.x < 0.975 and 0.025 < c.y < 0.975 and c.z > 0 for c in corners), corners
    result.update(
        native_atoms=len(ids),
        protein_shells=len(groups),
        branch_angstrom=measured,
        font=annotation.data.font.name,
    )
elif "test" in scenes:
    assert set(scenes) == {"test", "wall", "ct-section"}
    info = json.loads((out / "inventory.json").read_text())
    arr = np.load(out / "geometry.npz")
    shell = next(o for o in scenes["test"].objects if o.get("specimen_id") == "GE_GAD_2")
    assert shell in list(scenes["wall"].objects) and np.array_equal(coords(shell.data), arr["shell_v"])
    ob = next(o for o in scenes["ct-section"].objects if "source_z_index" in o)
    material = ob.data.materials[0]
    im = next(n.image for n in material.node_tree.nodes if n.type == "TEX_IMAGE")
    assert im.packed_file and im.colorspace_settings.name == "Non-Color"
    shape = tuple(info["source_shape_zyx"])
    raw = np.memmap(out.parent / "sources/GE_GAD_2.raw", dtype=">u2", mode="r", shape=shape)
    height, width = shape[1:]
    assert tuple(im.size) == (width, height)
    rgba = np.empty(width * height * 4, "f4")
    im.pixels.foreach_get(rgba)
    difference = float(
        np.max(np.abs(rgba.reshape(height, width, 4)[:, :, 0] * 65535.0 - raw[ob["source_z_index"]]))
    )
    assert difference < 0.025, difference
    pts = coords(ob.data)
    z = (info["origin_um_xyz"][2] + ob["source_z_index"] * ob["pixel_size_um"]) / 100.0
    assert np.allclose(pts[:, 2], z) and np.allclose(
        np.ptp(pts, axis=0)[:2], np.array([width, height]) * ob["pixel_size_um"] / 100.0
    )
    h = scenes["test"]
    sample = coords(shell.data)[::81]
    points = np.array([list(world_to_camera_view(h, h.camera, Vector(v))) for v in sample])
    assert np.min(points[:, :2]) > 0.015 and np.max(points[:, :2]) < 0.985
    result.update(
        shell_faces=len(shell.data.polygons),
        shell_vertices=len(shell.data.vertices),
        specimen=info["specimen"],
        packed_slice_intensity_max_error=difference,
        whole_camera_bounds=[points[:, :2].min(0).tolist(), points[:, :2].max(0).tolist()],
    )
else:
    assert set(scenes) == {"working-head", "ribbon"}
    for s in scenes.values():
        sections = sorted((o for o in s.objects if "section_order" in o), key=lambda o: o["section_order"])
        assert [o["section_order"] for o in sections] == list(range(1, 9))
    assert all(o in list(scenes["working-head"].objects) for o in sections)
    water = scenes["ribbon"].objects["Water surface"]
    knife = scenes["ribbon"].objects["Diamond cutting edge"]
    edge = np.array(knife["edge_xyz"]).reshape(2, 3)
    assert np.allclose(edge[:, 1], 0) and np.allclose(edge[:, 2], water["nominal_z"])
    bpy.context.window.scene = scenes["ribbon"]
    bpy.context.view_layer.update()
    w = water.evaluated_get(bpy.context.evaluated_depsgraph_get())
    clear = []
    for sec in sections:
        for co in coords(sec.data):
            above = Vector(co) + Vector((0, 0, 0.2))
            hit, point, _, _ = w.ray_cast(above, Vector((0, 0, -1)))
            assert hit
            clear.append(float(co[2] - point.z))
    assert min(clear) > 0 and max(clear) < 0.0016, clear
    gap = float(edge[0, 1] - max(coords(scenes["ribbon"].objects["Trimmed block"].data)[:, 1]))
    assert 0 < gap < 0.04
    result.update(
        sections=len(sections),
        max_water_clearance=max(clear),
        min_water_clearance=min(clear),
        block_to_edge=gap,
    )
out.mkdir(parents=True, exist_ok=True)
for stem, s in scenes.items():
    if a.render:
        bpy.context.window.scene = s
        s.render.filepath = str(out / (stem + ".png"))
        bpy.ops.render.render(write_still=True, scene=s.name)
result["rendered"] = bool(a.render)
(out / "saved-check.json").write_text(json.dumps(result, indent=2))
print(json.dumps(result))
