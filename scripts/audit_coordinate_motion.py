"""Check joined sparse molecular motion in an audited assembly with unchanged covalent topology.

Distances apply to piecewise-linear Cartesian coordinates between the supplied
samples, including the interior of each segment. Frames are presentation time.
"""

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


def positive(x, name, *, zero=False):
    if (
        isinstance(x, bool)
        or not isinstance(x, (int, float))
        or not np.isfinite(x)
        or (x < 0 if zero else x <= 0)
    ):
        raise ValueError("Invalid " + name)
    return float(x)


def closest_on_segments(start, end):
    delta = end - start
    denominator = np.einsum("ij,ij->i", delta, delta)
    t = np.divide(
        -np.einsum("ij,ij->i", start, delta),
        denominator,
        out=np.zeros_like(denominator),
        where=denominator > 0,
    )
    t = np.clip(t, 0, 1)
    return np.linalg.norm(start + t[:, None] * delta, axis=1), t


def read_model(candidate, receipt, plan):
    ledger = receipt.get("ledger", receipt)
    digest = sha256(candidate.read_bytes()).hexdigest()
    if (
        digest != receipt["candidate_sha256"]
        or digest != plan["candidate_sha256"]
        or (ledger["source_sha256"], ledger["assembly_id"])
        != (plan["source_sha256"], str(plan["assembly_id"]))
    ):
        raise ValueError("Motion plan belongs to a different audited candidate or assembly")
    with np.load(candidate, allow_pickle=False) as f:
        keys = [tuple(k) for k in f["atom_ids"].tolist()]
        if keys and len(keys[0]) == 6 and ledger.get("identity_only"):
            keys = [(k[0], "1", *k[1:]) for k in keys]
        xyz = f["model_A"]
        elements = f["elements"].tolist()
    if len(keys) != ledger["model_heavy_atoms"] or any(len(k) != 7 for k in keys):
        raise ValueError("Incomplete assembly identity table")
    return keys, xyz, elements


def sample_file(path, moves, expected_start):
    with np.load(path, allow_pickle=False) as f:
        ids = list(map(tuple, f["atom_ids"].tolist()))
        if len(ids) != len(moves) or set(ids) != set(moves):
            raise ValueError("Motion file differs from its declared moving-site inventory")
        positions = f["positions_A"]
        frames = f["frames"] if "frames" in f else np.arange(len(positions)) + expected_start
    if (
        positions.dtype.kind != "f"
        or positions.dtype.itemsize not in (4, 8)
        or positions.ndim != 3
        or positions.shape[1:] != (len(ids), 3)
        or len(positions) < 2
        or not np.isfinite(positions).all()
    ):
        raise ValueError("Invalid motion coordinate samples")
    if (
        frames.shape != (len(positions),)
        or not np.isfinite(frames).all()
        or np.any(np.diff(frames) <= 0)
        or frames[0] != expected_start
    ):
        raise ValueError("Motion frames must increase and join on the common boundary frame")
    return ids, positions, frames.astype(float)


def check_segments(
    keys, elements, state, mi, points, frames, edges, rest, exempt, base, floors, bond_limit
):
    local = np.full(len(keys), -1, dtype=int)
    local[mi] = np.arange(len(mi))
    fixed = np.flatnonzero(local < 0)
    tree = cKDTree(state[fixed])
    active = np.any(local[edges] >= 0, axis=1)
    edges, rest = edges[active], rest[active]
    max_floor = max([base, *[v for table in floors.values() for v in table.values()]])

    def positions(ii, sample):
        data = state[ii].copy()
        moving = local[ii] >= 0
        data[moving] = sample[local[ii[moving]]]
        return data

    worst = None
    bond_error = 0.0
    for step, (a, b) in enumerate(zip(points, points[1:])):
        if len(edges):
            first = positions(edges[:, 0], a) - positions(edges[:, 1], a)
            last = positions(edges[:, 0], b) - positions(edges[:, 1], b)
            nearest, _ = closest_on_segments(first, last)
            farthest = np.maximum(np.linalg.norm(first, axis=1), np.linalg.norm(last, axis=1))
            err = np.maximum(rest - nearest, farthest - rest)
            j = int(err.argmax())
            bond_error = max(bond_error, float(err[j]))
            if err[j] > bond_limit + 1e-9:
                raise ValueError(
                    f"Bond length changed by {err[j]:.6g} A at segment {step} for {keys[edges[j, 0]]}, {keys[edges[j, 1]]}"
                )
        middle = (a + b) / 2
        travel = np.linalg.norm(b - a, axis=1) / 2
        # Search an extra 0.5 A for a useful clearance measurement even on passing paths.
        radius = max_floor + 0.5
        neighborhoods = tree.query_ball_point(middle, radius + travel)
        pairs = [(i, int(fixed[j])) for i, js in enumerate(neighborhoods) for j in js]
        close = cKDTree(middle).query_pairs(radius + 2 * travel.max(), output_type="ndarray")
        pairs.extend((int(i), int(mi[j])) for i, j in close)
        eligible = [(i, j) for i, j in pairs if (min(mi[i], j), max(mi[i], j)) not in exempt]
        if not eligible:
            continue
        i, j = np.asarray(eligible, dtype=int).T
        required = np.array(
            [
                max(
                    base,
                    floors.get(keys[m][-2], {}).get(elements[n], 0),
                    floors.get(keys[n][-2], {}).get(elements[m], 0),
                )
                for m, n in zip(mi[i], j)
            ]
        )
        dist, t = closest_on_segments(a[i] - positions(j, a), b[i] - positions(j, b))
        margin = dist - required
        k = int(margin.argmin())
        minimum = float(margin[k])
        if worst is None or minimum < worst["margin_A"]:
            worst = dict(
                margin_A=minimum,
                distance_A=float(dist[k]),
                minimum_A=float(required[k]),
                frame=float(frames[step] + t[k] * (frames[step + 1] - frames[step])),
                ids=[keys[mi[i[k]]], keys[j[k]]],
            )
        if minimum < -1e-9:
            raise ValueError(
                f"Nonbonded collision at frame {worst['frame']:.6g}: {worst['ids']}, {worst['distance_A']:.6g} A < {worst['minimum_A']:.6g} A"
            )
    return dict(
        maximum_bond_length_change_A=bond_error,
        closest_screened_pair=worst,
        fixed_heavy_atoms=len(fixed),
        checked_covalent_edges=len(edges),
    )


def audit(candidate, receipt, plan, folder):
    keys, raw, elements = read_model(candidate, receipt, plan)
    index = {k: i for i, k in enumerate(keys)}

    def get(key):
        key = tuple(key)
        if key not in index:
            raise ValueError("Unknown assembly atom identity: " + str(key))
        return index[key]

    adj = defaultdict(set)
    edges = []
    for pair in plan["bonds"]:
        i, j = map(get, pair)
        if i == j or j in adj[i]:
            raise ValueError("Duplicate or self covalent edge")
        edges.append((i, j))
        adj[i].add(j)
        adj[j].add(i)
    edges = np.asarray(edges, dtype=int).reshape(-1, 2)
    state = raw.astype(float)
    rest = np.linalg.norm(state[edges[:, 0]] - state[edges[:, 1]], axis=1)
    exempt = {
        tuple(sorted((i, j))) for i, near in adj.items() for n in near for j in ({n} | adj[n]) if i != j
    }
    base = positive(plan["minimum_nonbonded_A"], "nonbonded floor")
    floors = {
        code: {el.upper(): positive(v, "component contact floor") for el, v in values.items()}
        for code, values in plan.get("component_floors_A", {}).items()
    }
    bond_limit = positive(plan["max_bond_change_A"], "bond tolerance", zero=True)
    frame = positive(plan.get("first_frame", 1), "first frame", zero=True)
    prior_float32 = np.full(len(keys), raw.dtype.itemsize == 4)
    records = []
    end_velocity = {}
    used = set()
    for stage in plan["stages"]:
        name = stage["id"]
        if not name or name in used:
            raise ValueError("Stage IDs must be distinct")
        used.add(name)
        moves = [tuple(k) for k in stage["moves"]]
        if not moves or len(set(moves)) != len(moves):
            raise ValueError("Missing or duplicate declared moving sites")
        src = folder / stage["file"]
        digest = sha256(src.read_bytes()).hexdigest()
        if "sha256" in stage and digest != stage["sha256"]:
            raise ValueError("Reviewed motion hash differs: " + name)
        ids, positions, frames = sample_file(src, moves, frame)
        mi = np.array([get(k) for k in ids])
        expected = state[mi]
        coarse = prior_float32[mi] | (positions.dtype.itemsize == 4)
        start = positions[0].astype(float)
        mismatch = np.any(abs(start - expected) > 1e-8, axis=1)
        mismatch[coarse] = np.any(start[coarse].astype("f4") != expected[coarse].astype("f4"), axis=1)
        distances = np.linalg.norm(start - expected, axis=1)
        if mismatch.any():
            j = int(np.argmax(np.where(mismatch, distances, -1)))
            raise ValueError(
                f"Boundary jump in {name}: {int(mismatch.sum())} sites; {ids[j]} moves {distances[j]:.6g} A at entry"
            )
        points = positions.astype(float, copy=False)
        steps = np.diff(points, axis=0)
        max_step = float(np.linalg.norm(steps, axis=2).max())
        if max_step > positive(stage["max_step_A"], "per-sample motion bound") + 1e-9:
            raise ValueError("Motion exceeds its declared sample displacement: " + name)
        velocity = dict(zip(mi, steps[0] / (frames[1] - frames[0])))
        dv = max(
            (
                float(np.linalg.norm(velocity.get(i, 0.0) - end_velocity.get(i, 0.0)))
                for i in set(velocity) | end_velocity.keys()
            ),
            default=0.0,
        )
        if (
            "max_entry_velocity_change_A_per_frame" in stage
            and dv
            > positive(
                stage["max_entry_velocity_change_A_per_frame"], "entry velocity tolerance", zero=True
            )
            + 1e-9
        ):
            raise ValueError("Velocity discontinuity at entry: " + name)
        measured = check_segments(
            keys, elements, state, mi, points, frames, edges, rest, exempt, base, floors, bond_limit
        )
        records.append(
            dict(
                id=name,
                motion_sha256=digest,
                moving_sites=len(mi),
                samples=len(frames),
                frames=[float(frames[0]), float(frames[-1])],
                max_start_error_A=float(distances.max()),
                max_sample_displacement_A=max_step,
                entry_velocity_change_A_per_frame=dv,
                **measured,
            )
        )
        end_velocity = dict(zip(mi, steps[-1] / (frames[-1] - frames[-2])))
        state[mi] = points[-1]
        prior_float32[mi] = positions.dtype.itemsize == 4
        frame = float(frames[-1])
    if not records:
        raise ValueError("Supply at least one coordinate motion stage")
    return dict(
        candidate_sha256=plan["candidate_sha256"],
        source_sha256=plan["source_sha256"],
        assembly_id=str(plan["assembly_id"]),
        assembly_heavy_atoms=len(keys),
        stages=records,
        topology="Constant reviewed covalent graph; bonded and 1-3 neighbors excluded from steric screen.",
        interpolation="Piecewise-linear Cartesian coordinates; frames are presentation time.",
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("candidate", "audit", "plan", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    a = p.parse_args()
    plan = json.loads(a.plan.read_text())
    result = audit(a.candidate, json.loads(a.audit.read_text()), plan, a.plan.parent)
    result["plan_sha256"] = sha256(a.plan.read_bytes()).hexdigest()
    result["coordinate_audit_sha256"] = sha256(a.audit.read_bytes()).hexdigest()
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"{len(result['stages'])} stages; {result['assembly_heavy_atoms']} assembly atoms")


if __name__ == "__main__":
    main()
