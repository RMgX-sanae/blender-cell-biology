"""Reopen a paired view and verify the displayed edits against the serialized candidate and source ledger."""

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import sys

import bpy
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Vector
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from detail_scene import linear
import coordinate_edit_animation as animation


def check(condition, message):
    if not condition:
        raise ValueError(message)


def close(actual, expected, label, atol=4e-5):
    check(np.allclose(actual, expected, atol=atol, rtol=0), label)


def audit(scene, args):
    camera = scene.camera
    check(camera and camera.data.type == "ORTHO", "Require a common orthographic camera")
    inputs = {
        name: sha256(path.read_bytes()).hexdigest()
        for name, path in (
            ("coordinate_candidate_sha256", args.candidate),
            ("coordinate_audit_sha256", args.audit),
            ("coordinate_view_sha256", args.view),
            ("appearance_profile_sha256", args.profile),
        )
    }
    check(all(scene.get(k) == v for k, v in inputs.items()), "Saved scene input fingerprints differ")
    receipt = json.loads(args.audit.read_text())
    ledger = receipt.get("ledger", receipt)
    view = json.loads(args.view.read_text())
    profile = json.loads(args.profile.read_text())
    check(
        inputs["coordinate_candidate_sha256"] == receipt["candidate_sha256"]
        and (ledger["source_sha256"], ledger["assembly_id"])
        == (view["source_sha256"], str(view["assembly_id"])),
        "Source or candidate identity differs",
    )
    ids = {k: tuple(v["id"]) for k, v in view["sites"].items()}
    roles = {ids[k]: v["role"] for k, v in view["sites"].items()}
    edits = {tuple(e["id"]): e for e in ledger["edits"]}
    with np.load(args.candidate, allow_pickle=False) as f:
        wanted = set(roles)
        keys = f["atom_ids"].tolist()
        if len(keys[0]) == 6 and ledger.get("identity_only"):
            keys = [[k[0], "1", *k[1:]] for k in keys]
        model = {tuple(k): p for k, p in zip(keys, f["model_A"], strict=True) if tuple(k) in wanted}
    check(model.keys() == roles.keys(), "Requested site absent from model")
    # The independent ledger gives the original coordinates of every changed site.
    source = {
        k: np.asarray(edits[k]["source_A"]) if k in edits else p
        for k, p in model.items()
        if k not in edits or edits[k]["source_A"] is not None
    }
    bonds = [(ids[a], ids[b]) for a, b in view["bonds"]]
    motion_plan = getattr(args, "motion_plan", None)
    motion_audit = getattr(args, "motion_audit", None)
    is_animated = "coordinate_motion_plan_sha256" in scene
    check(
        bool(motion_plan) == bool(motion_audit) == is_animated,
        "Supply the saved motion plan and review for animated views",
    )
    track = None
    if is_animated:
        check(
            scene["coordinate_motion_plan_sha256"] == sha256(motion_plan.read_bytes()).hexdigest()
            and scene["coordinate_motion_audit_sha256"] == sha256(motion_audit.read_bytes()).hexdigest(),
            "Saved motion fingerprints differ",
        )
        site_ids = list(ids.values())
        frames, track = animation.read_track(motion_plan, motion_audit, args.audit, site_ids, model)
        check(
            scene.frame_start == int(np.floor(frames[0])) and scene.frame_end == int(np.ceil(frames[-1])),
            "Saved frame range differs from the checked stages",
        )
        animation.set_frame(scene, frames[0])
    toward = np.asarray(view["toward_camera"], float)
    toward /= np.linalg.norm(toward)
    up = np.asarray(view["up"], float)
    up -= toward * np.dot(up, toward)
    up /= np.linalg.norm(up)
    basis = np.stack([np.cross(up, toward), up, toward])
    shifts = []
    counts = []
    materials = set()
    for i, expected in enumerate((source, model)):
        objects = [o for o in scene.objects if o.get("comparison_state") == i]
        atoms = {}
        half_bonds = {}
        parents = set()
        for obj in objects:
            check(
                obj.parent
                and not obj.parent.hide_render
                and not obj.hide_render
                and all(not c.hide_render for c in obj.users_collection),
                "Selected object hidden or unparented",
            )
            parents.add(obj.parent)
            close(
                obj.matrix_parent_inverse, np.eye(4), "Object no longer shares its panel coordinate frame"
            )
            if "atom_id" in obj:
                k = tuple(json.loads(obj["atom_id"]))
                check(
                    obj.type == "MESH" and k in expected and k not in atoms,
                    "Incorrect displayed atom identity",
                )
                atoms[k] = obj
                close(obj.location, expected[k], "Displayed atom moved from the audited position")
                rad = view["roles"][roles[k]]["radius_A"]
                close(obj.scale, [rad] * 3, "Incorrect or stretched displayed atom")
                close(
                    obj.matrix_basis.to_translation(),
                    expected[k],
                    "Effective displayed atom position differs",
                )
                close(
                    np.linalg.svd(np.asarray(obj.matrix_basis)[:3, :3], compute_uv=False),
                    [rad] * 3,
                    "Effective atom display scale differs",
                )
                unit = np.array([v.co[:] for v in obj.data.vertices])
                close(np.linalg.norm(unit, axis=1), 1, "Deformed atom mesh")
                p = world_to_camera_view(scene, camera, obj.matrix_world.translation)
                w = camera.data.ortho_scale
                aspect = (
                    scene.render.resolution_x
                    * scene.render.pixel_aspect_x
                    / (scene.render.resolution_y * scene.render.pixel_aspect_y)
                )
                check(
                    p.z > 0
                    and i / 2 + rad / w <= p.x <= (i + 1) / 2 - rad / w
                    and rad * aspect / w <= p.y <= 1 - rad * aspect / w,
                    "Displayed atom outside its panel",
                )
                mat = obj.active_material
                spec = view["roles"][roles[k]]
                mat_role = spec.get("surface_role", roles[k])
                check(mat and mat.get("surface_role") == mat_role, "Displayed role material differs")
                close(mat.diffuse_color[:3], linear(spec["color"]), "Material identity color differs")
                shader = mat.node_tree.nodes.get("Principled BSDF")
                close(
                    shader.inputs["Base Color"].default_value[:3],
                    linear(spec["color"]),
                    "Shader identity color differs",
                )
                params = profile.get("surface_roles", {}).get(mat_role, {})
                for key, socket, default in (
                    ("roughness", "Roughness", profile.get("roughness", 0.48)),
                    ("subsurface_weight", "Subsurface Weight", 0.025),
                ):
                    close(
                        shader.inputs[socket].default_value,
                        params.get(key, default),
                        "Profile shader setting differs",
                    )
                if "bump" in params and params["bump"] > 0:
                    normal = shader.inputs["Normal"]
                    check(
                        normal.is_linked and normal.links[0].from_node.type == "BUMP",
                        "Missing authored bump shading",
                    )
                    close(
                        normal.links[0].from_node.inputs["Distance"].default_value,
                        params["bump"],
                        "Profile bump distance differs",
                    )
                materials.add(roles[k])
            elif "bond_sites" in obj:
                k, j = map(tuple, json.loads(obj["bond_sites"]))
                check(
                    obj.type == "CURVE" and (k, j) not in half_bonds and k in expected and j in expected,
                    "Incorrect displayed bond identity",
                )
                half_bonds[k, j] = obj
                spline = obj.data.splines[0]
                check(len(obj.data.splines) == 1 and len(spline.points) == 2, "Changed bond topology")
                close(
                    [p.co[:3] for p in spline.points],
                    [expected[k], (expected[k] + expected[j]) / 2],
                    "Displayed bond endpoint differs",
                )
                close(obj.matrix_basis, np.eye(4), "Transformed displayed bond")
            else:
                raise ValueError("Missing displayed source identity")
        expected_edges = {
            (l, r) for a, b in bonds if a in expected and b in expected for l, r in ((a, b), (b, a))
        }
        check(
            atoms.keys() == expected.keys() and half_bonds.keys() == expected_edges,
            "Incomplete displayed inventory",
        )
        check(len(parents) == 1, "State components do not share one parent")
        parent = parents.pop()
        check(set(parent.children) == set(objects), "Unidentified geometry attached to state")
        matrix = np.asarray(parent.matrix_world)
        close(matrix[:3, :3], basis, "Reflected, stretched or rotated panel")
        for (left, right), obj in half_bonds.items():
            check(obj.active_material == atoms[left].active_material, "Bond color differs from endpoint")
        shifts.append(matrix[:3, 3])
        counts.append({"atoms": len(atoms), "bonds": len(half_bonds) // 2})
    close(
        shifts[1] - shifts[0],
        [camera.data.ortho_scale / 2, 0, 0],
        "Panels differ in depth, scale or common framing",
    )
    close(
        np.asarray(camera.matrix_world)[:3, :3], np.eye(3), "Camera no longer uses the reviewed projection"
    )
    labels = [o for o in scene.objects if o.type == "FONT"]
    check(
        Counter(o.data.body for o in labels)
        == Counter([*view["titles"], *(v["text"] for v in view.get("legend", []))]),
        "Saved labels differ from the selected view",
    )
    fonts = [sha256(p.read_bytes()).hexdigest() for p in (args.font_zh, args.font_en)]
    for o in labels:
        check(
            not o.hide_render
            and all(
                f.packed_file and sha256(bytes(f.packed_file.data)).hexdigest() == d
                for f, d in zip((o.data.font, o.data.font_bold), fonts)
            ),
            "Incorrect or unpacked panel font",
        )
        corners = [world_to_camera_view(scene, camera, o.matrix_world @ Vector(p)) for p in o.bound_box]
        check(all(p.z > 0 and 0 < p.x < 1 and 0 < p.y < 1 for p in corners), "Clipped panel label")
    report = dict(
        panels=counts,
        roles=sorted(materials),
        matched_fonts=True,
        common_orthographic_scale=camera.data.ortho_scale,
        source_sha256=ledger["source_sha256"],
        fingerprints=inputs,
        scope="Saved local geometry, framing, role shading and labels. Inspect the actual rendered image separately.",
    )
    if track is not None:
        report["motion"] = animation.verify(scene, site_ids, source, frames, track)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("candidate", "audit", "view", "output", "font-zh", "font-en"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument(
        "--profile",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "assets/style-structure-materials.json",
    )
    p.add_argument("--motion-plan", type=Path)
    p.add_argument("--motion-audit", type=Path)
    a = p.parse_args(sys.argv[sys.argv.index("--") + 1 :])
    r = audit(bpy.context.scene, a)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(r, indent=2) + "\n")
    print(r["panels"])


if __name__ == "__main__":
    main()
