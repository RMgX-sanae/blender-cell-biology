"""Open an existing Blender file, check its scene structure and render its named views."""

import argparse
import json
import sys
from pathlib import Path
import bpy
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--output", type=Path, required=True)
p.add_argument("--project", choices=("antenna", "diatom", "pipette"), required=True)
p.add_argument("--only", nargs="+")
p.add_argument("--check-only", action="store_true")
opts = p.parse_args(sys.argv[sys.argv.index("--") + 1 :])
out = opts.output.resolve()
out.mkdir(parents=True, exist_ok=True)
scenes = {s["delivery_stem"]: s for s in bpy.data.scenes if "delivery_stem" in s}
expected = {
    "antenna": ("antenna", "intact", "pigments", "binding_pair"),
    "diatom": ("frustule", "valve", "areolae"),
    "pipette": ("exploded", "assembled", "piston"),
}[opts.project]
assert set(scenes) == set(expected) and all(s.camera for s in scenes.values())
if opts.only and not set(opts.only) <= set(expected):
    p.error("Unknown named view")
details = {}
if opts.project == "antenna":
    s = scenes["intact"]
    objects = [o for o in s.objects if o.get("source") == "6KGX assembly 1"]
    ids = [o["label_asym"] for o in objects]
    assert len(set(ids)) == len(ids) == 2304
    assert sum(int(o["atoms"]) for o in objects) == 1011574
    assert (
        sum("source_host_label_asym" in o for o in scenes["binding_pair"].objects)
        == scenes["binding_pair"]["bound_bilins"]
    )
    bilins = sum(o.get("entity", 0) > 25 for o in scenes["pigments"].objects)
    assert bilins == 1598
    inv = json.loads((out / "prepared/inventory.json").read_text())
    raw = np.load(out / "prepared/meshes.npz")["vertices"]
    rr = {r["label_asym"]: r for r in inv["parts"]}
    maxerr = 0.0
    for obj in objects:
        row = rr[obj["label_asym"]]
        count = row["vertex_count"]
        assert len(obj.data.vertices) == count and len(obj.data.polygons) == row["face_count"]
        assert np.allclose(np.asarray(obj.matrix_world), np.eye(4), atol=1e-7)
        pos = np.empty(count * 3, dtype="f")
        obj.data.vertices.foreach_get("co", pos)
        maxerr = max(
            maxerr,
            float(
                np.abs(pos.reshape(-1, 3) - raw[row["vertex_start"] : row["vertex_start"] + count]).max()
            ),
        )
    assert maxerr < 1e-5
    details["all_envelope_vertices_max_error_nm"] = maxerr
    details.update(
        source_objects=len(objects),
        unique_source_ids=len(set(ids)),
        bilins=bilins,
        local_bilins=scenes["binding_pair"]["bound_bilins"],
    )
elif opts.project == "diatom":
    layers = {o["layer"] for o in scenes["areolae"].objects if "layer" in o}
    assert layers == {"foramen", "wall", "cribrum", "cribellum"}
    shell = scenes["frustule"]
    body = [o for o in shell.objects if o.type == "EMPTY" and o.animation_data]
    assert len(body) == 4
    poses = []
    for f in (1, 60):
        shell.frame_set(f)
        poses.append({o.name: list(o.location) for o in body})
    assert all(np.linalg.norm(np.array(x) - poses[1][k]) > 0.8 for k, x in poses[0].items())
    data = {o["layer"]: o.data for o in scenes["valve"].objects if "layer" in o}
    shell_meshes = [o for o in shell.objects if "layer" in o]
    assert all(o.data == data[o["layer"]] for o in shell_meshes) and len(shell_meshes) == 2 * len(data)
    from mathutils import Vector
    from mathutils.bvhtree import BVHTree

    macro = scenes["areolae"]
    bpy.context.window.scene = macro
    bpy.context.view_layer.update()
    mesh_by_layer = {o["layer"]: o for o in macro.objects if "layer" in o}

    def connected_sets(mesh):
        n = len(mesh.vertices)
        parent = list(range(n))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        edges = np.empty(len(mesh.edges) * 2, dtype="i")
        mesh.edges.foreach_get("vertices", edges)
        for i, j in edges.reshape(-1, 2):
            parent[find(int(i))] = find(int(j))
        groups = {}
        for i in range(n):
            groups.setdefault(find(i), []).append(i)
        return list(groups.values())

    plate = mesh_by_layer["foramen"]
    sections = connected_sets(plate.data)
    fine = connected_sets(mesh_by_layer["cribellum"].data)
    assert len(sections) == 7 and len(fine) == 343
    bvh = BVHTree.FromObject(plate, bpy.context.evaluated_depsgraph_get())
    open_rays = wall_rays = 0
    for ids in sections:
        xy = np.array([plate.data.vertices[i].co[:2] for i in ids])
        center = (xy.min(0) + xy.max(0)) / 2

        def hit(point):
            return bvh.ray_cast(Vector((*point, 3)), Vector((0, 0, -1)), 8)[0] is not None

        open_rays += not hit(center)
        radius = 0.4 * min(np.ptp(xy, axis=0))
        wall_rays += any(
            hit(center + [radius * np.cos(t), radius * np.sin(t)]) for t in np.arange(8) * np.pi / 4
        )
    assert open_rays == wall_rays == 7
    details.update(
        open_center_rays=int(open_rays),
        surrounding_plate_hits=int(wall_rays),
        fine_pore_mesh_components=len(fine),
    )
    details.update(
        separated_layers=sorted(layers), placed_valve_meshes=len(shell_meshes), registered_poses=poses
    )
else:
    pts = {o["assembly_group"] for o in scenes["assembled"].objects if "assembly_group" in o}
    assert len(pts) == 11
    refs = {o.data for o in scenes["assembled"].objects if o.get("assembly_group")}
    assert all(o.data in refs for o in scenes["exploded"].objects if "assembly_group" in o)
    fonts = [f for f in bpy.data.fonts if f.packed_file]
    assert (
        len(fonts) >= 2
        and sum(
            o.type == "FONT" and o.data.body in ("活塞与复位弹簧", "密封导向件与密封圈", "退吸头杆")
            for o in scenes["exploded"].objects
        )
        == 3
    )
    assembly = scenes["assembled"]

    def radial_section(name, z):
        obj = next(o for o in assembly.objects if o.name.startswith(name))
        pos = np.array([list(v.co) for v in obj.data.vertices])
        edge = np.array([list(e.vertices) for e in obj.data.edges])
        a = pos[edge[:, 0]]
        b = pos[edge[:, 1]]
        crossed = ((a[:, 2] - z) * (b[:, 2] - z) <= 0) & (np.abs(b[:, 2] - a[:, 2]) > 1e-7)
        a, b = a[crossed], b[crossed]
        xyz = a + (b - a) * ((z - a[:, 2]) / (b[:, 2] - a[:, 2]))[:, None]
        r = np.linalg.norm(xyz[:, :2], axis=1)
        return float(r.min()), float(r.max())

    piston = radial_section("Continuous polished piston", -1.46)[1]
    lip = radial_section("Annular seal guide", -1.46)[0]
    ring = radial_section("Elastomer O-ring", -1.71001)
    bore = radial_section("Tapered air channel", -1.71)[0]
    tip_inner = radial_section("Fine open polypropylene tip", -5.60)[0]
    cone = radial_section("Tapered air channel", -5.60)[1]
    assert (
        abs(piston - lip) < 0.006
        and abs(ring[0] - piston) < 0.006
        and -0.003 < bore - ring[1] < 0.02
        and abs(tip_inner - cone) < 0.01
    )
    details["model_radial_fit"] = {
        "piston_radius": piston,
        "seal_lip_radius": lip,
        "o_ring_inner_outer": ring,
        "shaft_bore_at_ring": bore,
        "tip_inner_radius": tip_inner,
        "holder_radius_at_tip": cone,
        "scope": "illustrative model in design units",
    }
    details.update(assembled_groups=len(pts), packed_fonts=[f.name for f in fonts], shared_parts=True)
result = {}
for stem in expected:
    s = scenes[stem]
    if opts.check_only or (opts.only and stem not in opts.only):
        continue
    bpy.context.window.scene = s
    bpy.context.view_layer.update()
    assert s.render.engine == "CYCLES" and s.render.resolution_percentage == 100
    s.render.filepath = str(out / (stem + ".png"))
    bpy.ops.render.render(write_still=True, scene=s.name)
    p = out / (stem + ".png")
    assert p.is_file() and p.stat().st_size > 15000
    result[stem] = {
        "scene": s.name,
        "width": s.render.resolution_x,
        "height": s.render.resolution_y,
        "samples": s.cycles.samples,
        "bytes": p.stat().st_size,
    }
(out / ("geometry_check.json" if opts.check_only else "saved_scene_check.json")).write_text(
    json.dumps(
        dict(project=opts.project, blender=bpy.app.version_string, checked=details, outputs=result),
        ensure_ascii=False,
        indent=2,
    )
)
