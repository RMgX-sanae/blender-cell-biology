"""Bind a checked sparse Cartesian path to the visible atoms and bond endpoints."""

from hashlib import sha256
import json
import math

import bpy
import numpy as np

from coordinate_animation import _curves


def read_track(plan_path, audit_path, coordinate_audit_path, site_ids, model):
    receipt = json.loads(audit_path.read_text())
    plan = json.loads(plan_path.read_text())
    coordinate = json.loads(coordinate_audit_path.read_text())
    if (
        receipt["plan_sha256"] != sha256(plan_path.read_bytes()).hexdigest()
        or receipt["coordinate_audit_sha256"] != sha256(coordinate_audit_path.read_bytes()).hexdigest()
        or receipt["candidate_sha256"] != coordinate["candidate_sha256"]
    ):
        raise ValueError("Motion, plan and coordinate receipts belong to different inputs")
    columns = {k: i for i, k in enumerate(site_ids)}
    state = np.asarray([model[k] for k in site_ids], dtype=float)
    times = []
    paths = []
    for spec, verified in zip(plan["stages"], receipt["stages"], strict=True):
        path = plan_path.parent / spec["file"]
        if (
            spec["id"] != verified["id"]
            or sha256(path.read_bytes()).hexdigest() != verified["motion_sha256"]
        ):
            raise ValueError("A reviewed motion stage was replaced")
        with np.load(path, allow_pickle=False) as f:
            ids = list(map(tuple, f["atom_ids"].tolist()))
            positions = f["positions_A"]
            frames = f["frames"] if "frames" in f else np.arange(len(positions)) + verified["frames"][0]
        use = [(i, columns[k]) for i, k in enumerate(ids) if k in columns]
        block = np.broadcast_to(state, (len(frames), *state.shape)).copy()
        if use:
            source, display = np.asarray(use).T
            block[:, display] = positions[:, source]
        # Frame zero of a later block is already present as the preceding end.
        skip = 1 if times else 0
        times.extend(frames[skip:].tolist())
        paths.append(block[skip:])
        state = block[-1]
    frames = np.asarray(times, dtype=float)
    if not len(frames) or np.any(np.diff(frames.astype("f4")) <= 0):
        raise ValueError("Motion frames are not distinct in Blender storage")
    return frames, np.concatenate(paths)


def set_frame(scene, time):
    whole = math.floor(time)
    scene.frame_set(whole, subframe=float(time - whole))


def uses_linear_samples(data, frames):
    curves = _curves(data)
    if not curves:
        raise ValueError("Missing native animation curves")
    for curve in curves:
        pts = curve.keyframe_points
        if (
            curve.modifiers
            or curve.extrapolation != "CONSTANT"
            or len(pts) != len(frames)
            or any(p.interpolation != "LINEAR" for p in pts)
            or not np.array_equal(np.array([p.co.x for p in pts], dtype="f4"), frames.astype("f4"))
        ):
            raise ValueError("Coordinate curves do not preserve the checked sample times or interpolation")


def bind(scene, site_ids, frames, points):
    columns = {k: i for i, k in enumerate(site_ids)}

    def key(data, values, point_index=None):
        if np.array_equal(values, np.broadcast_to(values[0], values.shape)):
            return
        accessor = "location" if point_index is None else f"splines[0].points[{point_index}].co"
        for t, v in zip(frames, values):
            if point_index is None:
                data.location = v
            else:
                data.splines[0].points[point_index].co = (*v, 1.0)
            data.keyframe_insert(data_path=accessor, frame=float(t))

    for obj in scene.objects:
        if obj.get("comparison_state") != 1:
            continue
        if "atom_id" in obj:
            values = points[:, columns[tuple(json.loads(obj["atom_id"]))]]
            data = obj
            key(data, values)
        else:
            a, b = (columns[tuple(k)] for k in json.loads(obj["bond_sites"]))
            data = obj.data
            key(data, points[:, a], 0)
            key(data, (points[:, a] + points[:, b]) / 2, 1)
        if data.animation_data:
            for curve in _curves(data):
                curve.extrapolation = "CONSTANT"
                for p in curve.keyframe_points:
                    p.interpolation = "LINEAR"
                curve.update()
            uses_linear_samples(data, frames)


def verify(scene, site_ids, source, frames, points):
    columns = {k: i for i, k in enumerate(site_ids)}
    atoms = []
    bonds = []
    camera = scene.camera
    if camera.animation_data or camera.constraints or camera.data.animation_data:
        raise ValueError("The matched-view camera must remain fixed during the motion")
    for obj in scene.objects:
        panel = obj.get("comparison_state")
        if panel is None:
            continue
        if obj.parent.animation_data or obj.parent.constraints:
            raise ValueError("The common panel frame must remain fixed throughout the motion")
        if obj.constraints or obj.modifiers:
            raise ValueError("Selected coordinate geometry includes an unreviewed deformation")
        if "atom_id" in obj:
            k = tuple(json.loads(obj["atom_id"]))
            atoms.append((obj, panel, k))
            animation = obj
        else:
            a, b = map(tuple, json.loads(obj["bond_sites"]))
            bonds.append((obj, panel, a, b))
            animation = obj.data
        if animation.animation_data:
            uses_linear_samples(animation, frames)
    # Inspect evaluated object transforms and curve geometry after loading the native file.
    times = np.empty(2 * len(frames) - 1)
    times[::2] = frames
    times[1::2] = (frames[1:] + frames[:-1]) / 2
    worst = 0.0
    world_worst = 0.0

    def error(actual, expected):
        residual = np.abs(np.asarray(actual) - np.asarray(expected))
        if not np.isfinite(residual).all():
            raise ValueError("Nonfinite evaluated Blender coordinates")
        return float(residual.max())

    original = scene.frame_current + scene.frame_subframe
    try:
        for j, t in enumerate(times):
            set_frame(scene, t)
            deps = bpy.context.evaluated_depsgraph_get()
            n = j // 2
            pose = points[n] if not j % 2 else (points[n] + points[n + 1]) / 2

            def co(panel, k):
                return pose[columns[k]] if panel == 1 else source[k]

            for obj, panel, k in atoms:
                ev = obj.evaluated_get(deps)
                local = np.asarray(ev.matrix_basis)[:3, 3]
                expected = co(panel, k)
                worst = max(worst, error(local, expected))
                parent = np.asarray(obj.parent.matrix_world)
                world_worst = max(
                    world_worst,
                    error(np.asarray(ev.matrix_world)[:3, 3], parent[:3, :3] @ expected + parent[:3, 3]),
                )
            for obj, panel, a, b in bonds:
                ev = obj.evaluated_get(deps)
                coords = np.asarray([p.co[:3] for p in ev.data.splines[0].points])
                target = np.stack([co(panel, a), (co(panel, a) + co(panel, b)) / 2])
                worst = max(worst, error(coords, target))
                if error(ev.matrix_world, obj.parent.matrix_world) > 4e-5:
                    raise ValueError("A bond moved away from its molecular coordinate frame")
    finally:
        set_frame(scene, original)
    if max(worst, world_worst) > 5e-5:
        raise ValueError(
            f"Saved motion differs from checked Cartesian samples by {max(worst, world_worst):.6g} A"
        )
    return dict(
        samples=len(frames),
        knots_and_midpoints=len(times),
        max_local_coordinate_error_A=worst,
        max_world_atom_coordinate_error_A=world_worst,
    )
