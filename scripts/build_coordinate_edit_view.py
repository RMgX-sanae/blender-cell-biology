"""Display a reviewed local coordinate edit in two editable, equally projected Blender panels."""

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Matrix

sys.path.insert(0, str(Path(__file__).resolve().parent))
from detail_scene import collection, linear, material, path, save_and_render, studio
from layered_molecular_appearance import apply_fonts, shade
import coordinate_edit_animation as animation


def read_selection(candidate, receipt, view):
    if sha256(candidate.read_bytes()).hexdigest() != receipt["candidate_sha256"]:
        raise ValueError("Coordinate file differs from the audited candidate")
    ledger = receipt.get("ledger", receipt)
    if (ledger["source_sha256"], ledger["assembly_id"]) != (
        view["source_sha256"],
        str(view["assembly_id"]),
    ):
        raise ValueError("View belongs to a different source revision or assembly")
    selected = {k: tuple(spec["id"]) for k, spec in view["sites"].items()}
    wanted = set(selected.values())
    if not wanted or len(wanted) != len(selected) or any(len(k) != 7 for k in wanted):
        raise ValueError("Select unique explicit assembly atom identities")
    new = {tuple(e["id"]) for e in ledger["edits"] if e["origin"] == "modeled addition"}
    with np.load(candidate, allow_pickle=False) as f:
        ids = [tuple(row) for row in f["atom_ids"].tolist()]
        if len(ids[0]) == 6 and ledger.get("identity_only"):
            ids = [(k[0], "1", *k[1:]) for k in ids]
        source_ids = [k for k in ids if k not in new]
        if (
            len(source_ids) != ledger["source_heavy_atoms"]
            or len(ids) != ledger["model_heavy_atoms"]
            or not wanted <= set(ids)
        ):
            raise ValueError("View selection or archive inventory disagrees with the audit")
        panels = []
        for keys, coords in ((source_ids, f["source_A"]), (ids, f["model_A"])):
            rows = {key: pos for key, pos in zip(keys, coords, strict=True) if key in wanted}
            panels.append({name: rows[k] for name, k in selected.items() if k in rows})
    if not panels[0]:
        raise ValueError("Include at least one deposited site in the comparison")
    return panels


def projection(view):
    def unit(x):
        v = np.asarray(x, dtype=float)
        if v.shape != (3,) or not np.isfinite(v).all() or np.linalg.norm(v) < 1e-7:
            raise ValueError("Choose finite, non-parallel view and up directions")
        return v / np.linalg.norm(v)

    toward = unit(view["toward_camera"])
    right = unit(np.cross(unit(view["up"]), toward))
    return np.stack([right, np.cross(toward, right), toward])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("candidate", "audit", "view", "output", "font-zh", "font-en"):
        p.add_argument("--" + arg, type=Path, required=True)
    p.add_argument(
        "--profile",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "assets/style-structure-materials.json",
    )
    p.add_argument("--engine", choices=["cycles", "eevee"], default="cycles")
    p.add_argument("--compositor", choices=["clean", "reference", "monochrome"], default="clean")
    p.add_argument("--width", type=int, default=1600)
    p.add_argument("--samples", type=int, default=32)
    p.add_argument("--render", action="store_true")
    p.add_argument("--motion-plan", type=Path)
    p.add_argument("--motion-audit", type=Path)
    p.add_argument("--snapshot-frame", type=int)
    a = p.parse_args(sys.argv[sys.argv.index("--") + 1 :])
    if a.output.exists():
        p.error("Choose a fresh output directory")
    a.frames, a.fps, a.frame, a.animate = 1, 24, 1, False
    view = json.loads(a.view.read_text())
    receipt = json.loads(a.audit.read_text())
    panels = read_selection(a.candidate, receipt, view)
    if bool(a.motion_plan) != bool(a.motion_audit):
        p.error("Supply both the motion plan and its audit")
    if a.snapshot_frame is not None and not a.motion_plan:
        p.error("Choose a snapshot frame for a motion input")
    model = {tuple(view["sites"][name]["id"]): coord for name, coord in panels[1].items()}
    site_ids = list(model)
    track = None
    if a.motion_plan:
        frames, track = animation.read_track(a.motion_plan, a.motion_audit, a.audit, site_ids, model)
        a.frames = math.ceil(frames[-1])
        a.frame = a.snapshot_frame if a.snapshot_frame is not None else a.frames
        if not math.floor(frames[0]) <= a.frame <= a.frames:
            p.error("Snapshot outside the checked motion")
    basis = projection(view)
    sites, roles = view["sites"], view["roles"]
    if len(view["titles"]) != 2 or len(set(view["titles"])) != 2 or not all(view["titles"]):
        raise ValueError("Supply distinct state titles")
    if any(not np.isfinite(r["radius_A"]) or r["radius_A"] <= 0 for r in roles.values()):
        raise ValueError("Display radii must be finite and positive")
    if any(
        len(pair) != 2 or len(set(pair)) != 2 or not set(pair) <= sites.keys() for pair in view["bonds"]
    ):
        raise ValueError("Each declared display bond must name two different selected sites")
    if len({frozenset(pair) for pair in view["bonds"]}) != len(view["bonds"]):
        raise ValueError("Duplicate display bond")
    points = np.array([point for panel in panels for point in panel.values()]) @ basis.T
    radii = np.array([roles[sites[k]["role"]]["radius_A"] for panel in panels for k in panel])[:, None]
    low = (points - radii).min(axis=0)
    high = (points + radii).max(axis=0)
    if track is not None:
        vr = np.array([roles[sites[n]["role"]]["radius_A"] for n in panels[1]])[:, None]
        projected = track @ basis.T
        low = np.minimum(low, (projected - vr).min(axis=(0, 1)))
        high = np.maximum(high, (projected + vr).max(axis=(0, 1)))
    middle = (high + low) / 2
    w, h, d = high - low
    aspect = 16 / 9
    width = max(2.32 * w, 1.57 * h * aspect)
    height = width / aspect
    z_front = d / 2 + 1
    bpy.ops.wm.read_factory_settings(use_empty=True)
    profile = json.loads(a.profile.read_text())
    # Use the selected profile's lighting and materials with the common framing computed above.
    profile["camera"] = {}
    studio(
        profile,
        a,
        dict(position=(0, 0, max(30, d * 3 + 10)), target=(0, 0, 0), scale=width, light_scale=width / 22),
    )
    scene = bpy.context.scene
    scene.render.resolution_y = round(a.width / aspect)
    models = [collection(n) for n in ("01 Deposited selection", "02 Edited selection")]
    labels = collection("03 Panel labels")
    mats = {}
    spheres = {}
    for name, spec in roles.items():
        mat = mats[name] = shade(
            profile,
            name,
            spec["color"],
            spec.get("surface_role", name),
            profile.get("roughness", 0.48),
            profile.get("bump", 0.005),
        )
        bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=3, radius=1)
        template = bpy.context.object
        for face in template.data.polygons:
            face.use_smooth = True
        template.data.materials.append(mat)
        spheres[name] = template.data
        bpy.data.objects.remove(template, do_unlink=True)
    ink_color = (
        "172B37" if np.dot(linear(profile["background"]), [0.2126, 0.7152, 0.0722]) > 0.45 else "F1F5E9"
    )
    ink = material("Text and divider", ink_color, 0.6, 0)
    # Use a flat ink so that changing spot highlights never obscures teaching labels.
    tree = ink.node_tree
    emission = tree.nodes.new("ShaderNodeEmission")
    emission.inputs[0].default_value = (*linear(ink_color), 1)
    tree.links.new(emission.outputs[0], tree.nodes.get("Material Output").inputs[0])

    def text(body, x, y, size):
        curve = bpy.data.curves.new(body, "FONT")
        curve.body = body
        curve.size, curve.align_x, curve.align_y = size, "CENTER", "CENTER"
        obj = bpy.data.objects.new(body, curve)
        labels.objects.link(obj)
        obj.location = (x, y, z_front)
        curve.materials.append(ink)

    for i, (panel, group, title) in enumerate(zip(panels, models, view["titles"], strict=True)):
        side = (-1 if i == 0 else 1) * width / 4
        parent = bpy.data.objects.new(group.name + " frame", None)
        group.objects.link(parent)
        matrix = np.eye(4)
        matrix[:3, :3] = basis
        matrix[:3, 3] = [side, 0, 0] - middle
        parent.matrix_world = Matrix(matrix.tolist())
        for name, position in panel.items():
            role = sites[name]["role"]
            obj = bpy.data.objects.new(name, spheres[role])
            group.objects.link(obj)
            obj.parent = parent
            obj.location = position
            obj.scale = [roles[role]["radius_A"]] * 3
            obj["atom_id"] = json.dumps(sites[name]["id"])
            obj["comparison_state"] = i
        for left, right in view["bonds"]:
            if left not in panel or right not in panel:
                continue
            midpoint = (panel[left] + panel[right]) / 2
            for start, end in ((left, right), (right, left)):
                r = sites[start]["role"]
                radius = 0.32 * min(roles[r]["radius_A"], roles[sites[end]["role"]]["radius_A"])
                obj = path(start + " — " + end, [panel[start], midpoint], radius, mats[r], group)
                obj.parent = parent
                obj["bond_sites"] = json.dumps([sites[start]["id"], sites[end]["id"]])
                obj["comparison_state"] = i
        text(title, side, 0.414 * height, 0.062 * height)
    path(
        "Panel divider",
        [(0, -0.32 * height, z_front), (0, 0.33 * height, z_front)],
        width * 0.001,
        ink,
        labels,
    )
    legend = view.get("legend", [])
    for i, item in enumerate(legend):
        x = width * ((i + 0.5) / len(legend) - 0.5) * 0.70
        text(item["text"], x, -0.417 * height, 0.046 * height)
        swatch = bpy.data.objects.new(item["role"] + " legend", spheres[item["role"]])
        labels.objects.link(swatch)
        swatch.location = (x, -0.354 * height, z_front)
        swatch.scale = [height * 0.012] * 3
    zh, en = [bpy.data.fonts.load(str(f), check_existing=True) for f in (a.font_zh, a.font_en)]
    apply_fonts(scene, zh, en)
    zh.pack()
    en.pack()
    scene["coordinate_candidate_sha256"] = receipt["candidate_sha256"]
    scene["coordinate_audit_sha256"] = sha256(a.audit.read_bytes()).hexdigest()
    scene["coordinate_view_sha256"] = sha256(a.view.read_bytes()).hexdigest()
    scene["appearance_profile_sha256"] = sha256(a.profile.read_bytes()).hexdigest()
    scene["display_length_unit"] = "angstrom"
    scene["display_site_count"] = json.dumps([len(panel) for panel in panels])
    if track is not None:
        scene.frame_start = math.floor(frames[0])
        scene["coordinate_motion_plan_sha256"] = sha256(a.motion_plan.read_bytes()).hexdigest()
        scene["coordinate_motion_audit_sha256"] = sha256(a.motion_audit.read_bytes()).hexdigest()
        animation.bind(scene, site_ids, frames, track)
    save_and_render(
        a,
        dict(
            view=view,
            source_sha256=view["source_sha256"],
            audit_sha256=scene["coordinate_audit_sha256"],
            candidate_sha256=receipt["candidate_sha256"],
            shown_atoms=[len(p) for p in panels],
            view_basis=basis.tolist(),
            display_length_unit="angstrom",
            scope="Selected sites in a common orthographic projection; display bonds are the project's reviewed graph.",
        ),
    )


if __name__ == "__main__":
    main()
