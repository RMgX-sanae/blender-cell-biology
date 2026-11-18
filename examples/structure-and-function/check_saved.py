"""Reopen native meshes, check biological or optical relationships and render actual saved scenes."""

import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Vector
from mathutils.bvhtree import BVHTree


def vertices(obj):
    p = np.empty(len(obj.data.vertices) * 3, dtype=np.float32)
    obj.data.vertices.foreach_get("co", p)
    return p.reshape(-1, 3)


def source_mesh(obj, v, f):
    assert np.array_equal(vertices(obj), v)
    a = np.empty(len(obj.data.loops), dtype="i")
    obj.data.loops.foreach_get("vertex_index", a)
    assert np.array_equal(a, f.ravel()) and np.array_equal(np.array(obj.matrix_world), np.eye(4))


def protein(project, main, detail):
    meta = json.loads((project / "inventory.json").read_text())
    with np.load(project / "prepared.npz") as archive:
        matched = {}
        for obj in main.objects:
            if "fu" in obj:
                k = obj["label"] + "_" + obj["fu"]
                source_mesh(obj, archive[k + "_v"], archive[k + "_f"])
                matched[k] = obj
    assert set(matched) == {r["key"] for r in meta["parts"]} and len(matched) == 80
    labels = Counter(o["label"] for o in detail.objects if "fu" in o)
    assert len(labels) == 6 and set(labels.values()) == {8}
    return dict(
        complete_protomers=10,
        domain_meshes=len(matched),
        source_atoms=meta["source_atoms"],
        triangles=meta["triangles"],
        visible_protomers_in_window=sorted(labels),
        mesh_save_deviation_nm=0,
    )


def mesh_parts(obj):
    xyz = vertices(obj)
    parent = list(range(len(xyz)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    edges = np.empty(2 * len(obj.data.edges), dtype="i")
    obj.data.edges.foreach_get("vertices", edges)
    for a, b in edges.reshape(-1, 2):
        parent[root(a)] = root(b)
    parts = {}
    for k in range(len(xyz)):
        parts.setdefault(root(k), []).append(k)
    return xyz, list(parts.values())


def tissue(project, main, detail, field):
    inner = [o for o in main.objects if o.name.startswith("Inner hair bundle")]
    outer = [o for o in main.objects if o.name.startswith("Outer hair bundle")]
    assert (len(inner), len(outer)) == (6, 18)
    dg = bpy.context.evaluated_depsgraph_get()
    roots = []
    link_short = []
    link_tall = []
    signed_ends = []
    for obj in inner + outer:
        (col,) = obj.users_collection
        cap = next(o for o in col.objects if "apical mosaic" in o.name)
        cap = cap.evaluated_get(dg)
        xyz, parts = mesh_parts(obj)
        shafts = [ids for ids in parts if len(ids) > 20]
        links = [ids for ids in parts if len(ids) <= 20]
        ncol = 27 if obj in inner else 23
        assert len(shafts) == 3 * ncol and len(links) == 2 * ncol
        mesh = obj.data
        ranks = [
            BVHTree.FromPolygons(xyz, [list(p.vertices) for p in mesh.polygons if p.material_index == r])
            for r in range(3)
        ]
        for ids in shafts:
            points = xyz[ids]
            low = points[:, 2].min()
            root = points[points[:, 2] < low + 0.00015].mean(0)
            hit, point, _, _ = cap.ray_cast(Vector(root) + Vector((0, 0, 2)), (0, 0, -1))
            assert hit, (obj.name, "outside its apical cell surface", root)
            roots.append(float(point.z - root[2]))
        for ids in links:
            points = xyz[ids]
            p = points[points[:, 2] < points[:, 2].min() + 0.0001].mean(0)
            q = points[points[:, 2] > points[:, 2].max() - 0.0001].mean(0)
            choices = [(ranks[r].find_nearest(p), ranks[r + 1].find_nearest(q)) for r in (0, 1)]
            a, b = min(choices, key=lambda pair: sum(hit[3] for hit in pair))
            link_short.append(a[3])
            link_tall.append(b[3])
            signed_ends.extend([(Vector(p) - a[0]).dot(a[1]), (Vector(q) - b[0]).dot(b[1])])
    assert min(roots) >= -0.002 and max(roots) < 0.26, (min(roots), max(roots))
    assert max(link_short) < 0.040 and max(link_tall) < 0.040, (max(link_short), max(link_tall))
    assert max(signed_ends) <= 0.002, max(signed_ends)
    sample_links = sum(o["illustrated_tip_links"] for o in inner + outer)
    for px, z in [(0.2, 10), (0.2, 14), (-1.3, 19)]:
        for collection in ("Sensory cells", "Pillars and Deiters cells"):
            for o in bpy.data.collections[collection].all_objects:
                if o.type == "MESH":
                    assert not o.evaluated_get(dg).ray_cast((px, -25, z), (0, 1, 0), distance=60)[0], (
                        px,
                        z,
                        o.name,
                    )
    assert all(o in list(main.objects) for o in detail.objects if "stereocilia_per_rank" in o)
    info = json.loads((project / "field.json").read_text())
    got = {}
    with np.load(project / "field.npz") as arc:
        for obj in field.objects:
            if "source_label" in obj:
                key = obj["source_label"]
                source_mesh(obj, arc[f"v{key}"], arc[f"f{key}"])
                got[key] = obj["source_voxels"]
    assert got == {r["label"]: r["voxels"] for r in info["bundles"]} and len(got) == 40
    return dict(
        inner_cells=len(inner),
        outer_cells=len(outer),
        illustrated_stereocilia=len(roots),
        root_insertion_range_um=[min(roots), max(roots)],
        tip_links=sample_links,
        maximum_link_end_distance_um=[max(link_short), max(link_tall)],
        maximum_signed_link_end_distance_um=max(signed_ends),
        longitudinal_tunnel_rays_clear=3,
        annotated_P21_bundles=len(got),
        annotated_voxel_count=sum(got.values()),
        raw_pitch_um=info["source_xyz_um"],
    )


def scanner(main, detail):
    rotor = main.objects["Common spindle axis"]
    mask = next(o for o in main.objects if o.name.startswith("Nipkow array"))
    lens = next(o for o in main.objects if o.name.startswith("Microlens array"))
    vertices_lens = vertices(lens)
    # Each independent terminal lenslet face marks its actual convex apex.
    centers = np.array(
        [vertices_lens[list(p.vertices)].mean(0) for p in lens.data.polygons if len(p.vertices) == 16]
    )
    assert len(centers) == 648 and set((mask.parent, lens.parent)) == {rotor}
    z = main["pinhole_plane_z"]

    def test_angle(theta):
        rotor.rotation_euler.z = theta
        bpy.context.view_layer.update()
        obj = mask.evaluated_get(bpy.context.evaluated_depsgraph_get())
        inv = obj.matrix_world.inverted()
        for x, y, _ in centers:
            radial = np.array([x, y]) / np.linalg.norm([x, y])
            # A vertical line through the upper mesh must also cross an opening in the lower mesh.
            pts = [(x, y, z + 0.14), (x + 0.064 * radial[0], y + 0.064 * radial[1], z + 0.14)]
            hits = [
                obj.ray_cast(inv @ (lens.matrix_world @ Vector(p)), (0, 0, -1), distance=0.25)[0]
                for p in pts
            ]
            assert hits == [False, True], (theta, x, y, hits)

    test_angle(0)
    test_angle(0.47)
    rotor.rotation_euler.z = 0
    bpy.context.view_layer.update()
    mirror = main.objects["Fixed dichroic mirror"]
    normal = mirror.matrix_world.to_quaternion() @ Vector((0, 0, 1))
    vertical = Vector((0, 0, 1))
    direction = vertical - 2 * vertical.dot(normal) * normal
    assert (direction - Vector((1, 0, 0))).length < 0.00001 and mirror.parent is None
    for i in range(3):
        incoming = main.objects[f"{i} returning fluorescence through pinhole"]
        reflected = main.objects[f"{i} reflected fluorescence"]
        hit = Vector(incoming["path_to"])
        assert (hit - Vector(reflected["path_from"])).length < 1e-5
        assert abs((hit - mirror.location).dot(normal)) < 1e-5
        assert ((Vector(reflected["path_to"]) - hit).normalized() - direction).length < 1e-5
        assert np.linalg.norm(centers[:, :2] - np.array(hit)[:2], axis=1).min() < 1e-5
        bead = main.objects[f"{i} fluorescent bead"]
        source = main.objects[f"{i} fluorescence from specimen"]
        assert (Vector(source["path_from"]) - Vector(bead["path_from"])).length < 1e-5
        train = [
            source,
            main.objects[f"{i} fluorescence in objective"],
            incoming,
            reflected,
            main.objects[f"{i} emission relay"],
        ]
        for first, second in zip(train, train[1:]):
            assert (Vector(first["path_to"]) - Vector(second["path_from"])).length < 1e-5
        sensor = main.objects["Camera sensor plane"].evaluated_get(bpy.context.evaluated_depsgraph_get())
        local = sensor.matrix_world.inverted()
        start = local @ Vector(train[-1]["path_from"])
        end = local @ Vector(train[-1]["path_to"])
        assert sensor.ray_cast(start, (end - start).normalized(), distance=(end - start).length)[0]
    return dict(
        open_pinhole_lens_pairs=len(centers),
        checked_rotor_angles_rad=[0, 0.47],
        adjacent_material_hits_at_each_angle=648,
        return_rays_to_sensor=3,
        reflected_direction=list(direction),
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--project", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--render", action="store_true")
    a = p.parse_args(sys.argv[sys.argv.index("--") + 1 :])
    scenes = sorted(bpy.data.scenes, key=lambda s: s.name)
    bpy.context.window.scene = scenes[0]
    for s in scenes:
        assert s.camera and s.compositing_node_group
        for branch in s["detail_compositor_branches"].values():
            assert all(n.scene == s for n in branch.nodes if n.type == "R_LAYERS")
    kind = a.project.name
    if kind == "01_hemocyanin":
        facts = protein(a.project, *scenes)
    elif kind == "02_corti":
        facts = tissue(a.project, *scenes)
    elif kind == "03_confocal":
        facts = scanner(*scenes)
    else:
        raise ValueError(kind)
    a.out.mkdir(parents=True, exist_ok=True)
    info = dict(project=kind, blender=bpy.app.version_string, scenes=[s.name for s in scenes], **facts)
    (a.out / "saved-check.json").write_text(json.dumps(info, indent=2))
    if a.render:
        for s in scenes:
            bpy.context.window.scene = s
            s.render.filepath = str((a.out / (s["delivery_stem"] + ".png")).resolve())
            bpy.ops.render.render(write_still=True, scene=s.name)


if __name__ == "__main__":
    main()
