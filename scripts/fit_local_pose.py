"""Fit selected atoms of a reviewed coordinate archive to project-defined geometric restraints."""
import argparse
from hashlib import sha256
from io import BytesIO
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from scipy.optimize import least_squares

from reaction_pose_geometry import signed_height


ARITY = {"distances": 2, "min_distances": 2, "angles": 3, "signed_heights": 4, "dihedrals": 4}
COORDINATE_FIELDS = ("atom_ids", "elements", "source_A", "model_A")


def values(kind, points):
    """One measurement per row, in Å or degrees."""
    a, b = points[:, 0], points[:, 1]
    if kind in ("distances", "min_distances"):
        return np.linalg.norm(a - b, axis=1)
    if kind == "signed_heights":
        return np.array([signed_height(row) for row in points])
    c = points[:, 2]
    u, v = a - b, c - b
    lu, lv = np.linalg.norm(u, axis=1), np.linalg.norm(v, axis=1)
    if np.any(lu == 0) or np.any(lv == 0):
        raise ValueError(f"Coincident anchors in {kind}")
    if kind == "angles":
        return np.degrees(np.arctan2(np.linalg.norm(np.cross(u, v), axis=1), (u * v).sum(axis=1)))
    # Signed rotation from the first arm to the last around b -> c, as in Gemmi.
    axis = v / lv[:, None]
    w = points[:, 3] - c
    arm = u - (u * axis).sum(axis=1)[:, None] * axis
    last = w - (w * axis).sum(axis=1)[:, None] * axis
    if (np.any(np.linalg.norm(arm, axis=1) <= 1e-6 * lu)
            or np.any(np.linalg.norm(last, axis=1) <= 1e-6 * np.linalg.norm(w, axis=1))):
        raise ValueError("Cannot orient a dihedral with collinear or coincident anchors")
    return np.degrees(np.arctan2((np.cross(arm, last) * axis).sum(axis=1), (arm * last).sum(axis=1)))


def number(v):
    if type(v) not in (int, float) or not np.isfinite(v):
        raise ValueError("Restraint values must be finite numbers")
    return float(v)


def fit(arrays, plan):
    if (type(plan["version"]) is not int or plan["version"] != 1
            or not isinstance(plan["basis"], str) or not plan["basis"].strip()):
        raise ValueError("Supply a version-1 fit plan with the basis for the selected sites and restraints")
    ids, full = arrays["atom_ids"], arrays["model_A"]
    if ids.ndim != 2 or ids.shape[1] != 7 or ids.dtype.kind != "U":
        raise ValueError("Expected seven string fields per assembly atom ID")
    if (full.shape != (len(ids), 3) or full.dtype.kind != "f" or full.dtype.itemsize not in (4, 8)
            or not np.isfinite(full).all() or arrays["elements"].shape != (len(ids),)):
        raise ValueError("Expected finite float32 or float64 model coordinates and one element per atom")
    source = arrays["source_A"]
    if (source.ndim != 2 or source.shape[1] != 3 or len(source) > len(ids)
            or source.dtype.kind != "f" or not np.isfinite(source).all()):
        raise ValueError("Expected finite deposited coordinates in source_A")
    lookup = {tuple(k): i for i, k in enumerate(ids.tolist())}
    sites = plan["sites"]
    if len(lookup) != len(ids) or not isinstance(sites, dict) or not sites:
        raise ValueError("Expected a unique source inventory and named sites")
    rows = []
    for key in sites.values():
        if not isinstance(key, list) or not all(isinstance(k, str) for k in key) or tuple(key) not in lookup:
            raise ValueError(f"Unknown assembly atom ID: {key}")
        rows.append(lookup[tuple(key)])
    names = {name: i for i, name in enumerate(sites)}
    moving = plan["mobile"]
    if (len(set(rows)) != len(rows) or not isinstance(moving, list) or not moving
            or len(set(moving)) != len(moving) or not set(moving) <= names.keys()):
        raise ValueError("Give each physical site one name and select distinct named mobile sites")
    mobile = np.array([names[n] for n in moving])
    start = full[rows].astype(float)
    tether = number(plan["position_sigma_A"])
    budget = plan.get("max_nfev", 150)
    if tether <= 0 or type(budget) is not int or budget < 1:
        raise ValueError("Use a positive position_sigma_A and integer evaluation budget")
    groups = {}
    used = set()
    for kind, rules in plan["restraints"].items():
        if kind not in ARITY or not isinstance(rules, list):
            raise ValueError(f"Unknown restraint group or invalid rule list: {kind}")
        if not rules:
            continue
        indices, targets, sigmas = [], [], []
        for rule in rules:
            atoms = rule["atoms"]
            if (not isinstance(atoms, list) or len(atoms) != ARITY[kind] or len(set(atoms)) != len(atoms)
                    or not set(atoms) <= names.keys()):
                raise ValueError(f"Expected {ARITY[kind]} distinct named sites in {kind}")
            indices.append([names[a] for a in atoms])
            used.update(atoms)
            targets.append(rule["target"])
            sigmas.append(number(rule["sigma"]))
        indices = np.array(indices)
        ref = values(kind, start[indices])
        target = np.array([v if t == "start" else number(t) for t, v in zip(targets, ref)])
        sigma = np.array(sigmas)
        if (np.any(sigma <= 0) or (kind in ("distances", "min_distances") and np.any(target <= 0))
                or (kind == "angles" and np.any((target < 0) | (target > 180)))):
            raise ValueError(f"Invalid {kind} target or non-positive sigma")
        groups[kind] = (indices, target, sigma)
    if not set(moving) <= used:
        raise ValueError("Every mobile site needs a geometric restraint")

    def errors(coords):
        result = {}
        for kind, (indices, target, sigma) in groups.items():
            delta = values(kind, coords[indices]) - target
            if kind == "min_distances":
                delta = np.minimum(delta, 0)
            elif kind == "dihedrals":
                delta = (delta + 180) % 360 - 180
            result[kind] = delta / sigma
        return result

    def residual(vector):
        coords = start.copy()
        coords[mobile] = vector.reshape(-1, 3)
        return np.concatenate([*errors(coords).values(), ((coords[mobile] - start[mobile]) / tether).ravel()])

    result = least_squares(residual, start[mobile].ravel(), max_nfev=budget, ftol=2e-8, xtol=2e-8, gtol=2e-8)
    if not result.success:
        raise ValueError(f"Local fit did not converge: {result.message}")
    fitted = full.copy()
    fitted[np.asarray(rows)[mobile]] = result.x.reshape(-1, 3)
    stored = fitted[rows].astype(float)
    report = dict(status="fitted_candidate", evaluations=result.nfev, mobile_sites=len(mobile),
                  maximum_displacement_A=float(np.linalg.norm(stored[mobile] - start[mobile], axis=1).max()),
                  peak_scaled_residual={k: float(np.max(abs(v))) for k, v in errors(stored).items()})
    return fitted, report


def fit_archive(start, plan_path, output):
    source_bytes = start.read_bytes()
    source_hash = sha256(source_bytes).hexdigest()
    plan_bytes = plan_path.read_bytes()
    plan = json.loads(plan_bytes)
    if plan["start_sha256"] != source_hash:
        raise ValueError("Starting coordinate archive differs from the reviewed fit plan")
    if output.exists():
        raise FileExistsError(f"Choose a fresh output path: {output}")
    with np.load(BytesIO(source_bytes), allow_pickle=False) as src:
        arrays = {k: src[k] for k in COORDINATE_FIELDS}
    arrays["model_A"], report = fit(arrays, plan)
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent) as tmp:
        ready = Path(tmp) / "candidate.npz"
        np.savez_compressed(ready, **arrays)
        report.update(start_sha256=source_hash, plan_sha256=sha256(plan_bytes).hexdigest(),
                      candidate_sha256=sha256(ready.read_bytes()).hexdigest())
        os.link(ready, output)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("start", "plan", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(fit_archive(args.start, args.plan, args.output), indent=2))
