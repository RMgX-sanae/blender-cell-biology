"""Prepare the source-qualified, pre-attack ester reorientation for 9NO7.

The driving and compensating torsions are a teaching hypothesis. No new
covalent bond is formed. Calculate and screen motion in the complete assembly.
"""

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import gemmi
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from audit_coordinate_motion import audit as audit_motion
from audit_termination_substrate import audit_directory, references
from event_timing import phase_progress
from prepare_phosphotransfer_motion import smooth


def build_coordinates(candidate, source_path, receipt, mode, frames):
    with np.load(candidate / "coordinates.npz", allow_pickle=False) as f:
        ids = list(map(tuple, f["atom_ids"].tolist()))
        if ids and len(ids[0]) == 6 and receipt["ledger"]["identity_only"]:
            ids = [(k[0], "1", *k[1:]) for k in ids]
        xyz = f["model_A"]
        els = f["elements"]
    idx = {k: i for i, k in enumerate(ids)}
    s = ("CB", "1", "x", "76", "", "A")
    c3 = (*s, "C3'")
    o3 = (*s, "O3'")
    o2 = (*s, "O2'")
    pep = [k for k in ids if k[:3] == ("DB", "1", "z")]
    nodes = [c3, o3, *pep]
    lookup = {k: i for i, k in enumerate(nodes)}

    def site(num, name):
        matches = [k for k in pep if int(k[3]) == num and k[-1] == name and k[4] == ""]
        if len(matches) != 1:
            raise ValueError("Missing or ambiguous peptide site")
        return matches[0]

    block = gemmi.cif.read(str(source_path)).sole_block()
    ccd = block.get_mmcif_category("_chem_comp_bond.")
    orders = {
        (c, frozenset((x, y))): order
        for c, x, y, order in zip(ccd["comp_id"], ccd["atom_id_1"], ccd["atom_id_2"], ccd["value_order"])
    }
    graph = defaultdict(set)
    bonds = []

    def bond(a, b):
        i, j = lookup[a], lookup[b]
        if j in graph[i]:
            raise ValueError("Duplicate bond")
        graph[i].add(j)
        graph[j].add(i)
        bonds.append((i, j))

    for pref in dict.fromkeys(k[:-1] for k in pep):
        for code, x, y in zip(ccd["comp_id"], ccd["atom_id_1"], ccd["atom_id_2"]):
            a, b = (*pref, x), (*pref, y)
            if code == pref[-1] and a in lookup and b in lookup:
                bond(a, b)
    atoms = block.get_mmcif_category("_atom_site.")
    numbering = {
        (n, seq)
        for label, n, seq in zip(atoms["label_asym_id"], atoms["auth_seq_id"], atoms["label_seq_id"])
        if label == "DB"
    }
    if numbering != {(str(n), str(n)) for n in range(3, 8)}:
        raise ValueError("Peptide sequence or the deposited extent changed")
    for n in (3, 4, 5, 6):
        bond(site(n, "C"), site(n + 1, "N"))
    bond(c3, o3)
    bond(o3, site(7, "C"))
    conn = block.get_mmcif_category("_struct_conn.")
    external_links = [
        name
        for name, t, a, b in zip(
            conn["id"], conn["conn_type_id"], conn["ptnr1_label_asym_id"], conn["ptnr2_label_asym_id"]
        )
        if t == "covale" and "DB" in (a, b)
    ]
    if external_links != [receipt["source_connections"]["ester"]]:
        raise ValueError("Unexpected peptide crosslink in the source")
    real = np.array([idx[k] for k in nodes])
    start = xyz[real]
    for i, j in bonds:
        length = np.linalg.norm(start[i] - start[j])
        if not 1.1 < length < 1.95:
            raise ValueError(f"Unqualified covalent edge {nodes[i]}, {nodes[j]}")
    exempt = {frozenset((real[i], real[j])) for i, j in bonds}
    for i, near in graph.items():
        for j in near:
            for k in graph[j]:
                if i != k:
                    exempt.add(frozenset((real[i], real[k])))

    def joint(a, b):
        i, j = lookup[a], lookup[b]
        if a[:-1] != b[:-1] or orders.get((a[-2], frozenset((a[-1], b[-1])))) != "sing":
            raise ValueError("The chosen torsion must be a reviewed intracomponent single bond")
        seen = {i}
        todo = [j]
        members = []
        while todo:
            k = todo.pop()
            if k in seen:
                continue
            seen.add(k)
            members.append(k)
            todo.extend(graph[k] - seen)
        if set(graph[i]).intersection(members) != {j}:
            raise ValueError("The selected torsion would open a covalent ring")
        return i, j, np.asarray([k for k in sorted(members) if k != j], int)

    joints = [joint(c3, o3)]
    if set(joints[0][2]) != set(range(2, len(nodes))):
        raise ValueError("Disconnected deposited peptide atoms would be left behind")
    for n in (7, 6, 5, 4):
        joints.extend([joint(site(n, "C"), site(n, "CA")), joint(site(n, "CA"), site(n, "N"))])

    def pose(q):
        pts = start.copy()
        for (i, j, down), deg in zip(joints, q):
            if not len(down) or deg == 0:
                continue
            axis = pts[j] - pts[i]
            axis /= np.linalg.norm(axis)
            pts[down] = Rotation.from_rotvec(np.radians(deg) * axis).apply(pts[down] - pts[j]) + pts[j]
        return pts

    mobile = real[2:]
    fixed = np.flatnonzero(~np.isin(np.arange(len(ids)), mobile))
    tree = cKDTree(xyz[fixed])
    pairs = []
    for i, js in enumerate(tree.query_ball_point(start[2:], 6.5), 2):
        for j in js:
            fj = fixed[j]
            if frozenset((real[i], fj)) not in exempt:
                el = els[real[i]].upper()
                floor = (
                    2.38
                    if ids[fj][-2:] != ("HOH", "O")
                    else 2.90
                    if el in ("C", "P", "S")
                    else 2.50
                    if el in ("N", "O")
                    else 2.38
                )
                pairs.append((i, fj, floor))
    who = np.asarray([a for a, _, _ in pairs])
    around = xyz[[b for _, b, _ in pairs]]
    floors = np.asarray([v for _, _, v in pairs])
    intra = np.array(
        [
            (i, j)
            for i in range(2, len(nodes))
            for j in range(2, i)
            if frozenset((real[i], real[j])) not in exempt
        ]
    )
    carbon = lookup[site(7, "C")]
    ox = lookup[site(7, "O")]

    def residual(tail, deg):
        pts = pose([deg, *tail])
        m = np.linalg.norm(pts[who] - around, axis=1) - floors
        n = np.linalg.norm(pts[intra[:, 0]] - pts[intra[:, 1]], axis=1) - 2.38
        delta = pts[2:] - start[2:]
        return np.r_[
            0.13 * delta.ravel(),
            0.011 * tail,
            28 * np.minimum(m - 0.04, 0),
            28 * np.minimum(n - 0.04, 0),
            4 * np.maximum(np.linalg.norm(delta, axis=1) - 1.65, 0),
        ]

    # Drive the C3'-O3' torsion. Adjust only permitted backbone torsions
    # on the connected peptide to avoid forcing its remote end through the tunnel.
    targets = np.linspace(0, -18, 19)
    q = np.zeros(len(joints) - 1)
    knots = []
    for deg in targets:
        if deg:
            fit = least_squares(
                residual, q, args=(deg,), bounds=(-30, 30), max_nfev=180, ftol=1e-8, gtol=1e-8, xtol=1e-8
            )
            if not fit.success:
                raise ValueError("Torsion fit did not converge at the reviewed driving angle")
            q = fit.x
        knots.append([deg, *q])
    knots = np.asarray(knots)
    time = np.linspace(0, 1, frames)
    interval = (0.0, 1.0) if mode == "continuous" else (0.15, 0.80)
    fraction = smooth(phase_progress(time, {"approach": interval})["approach"])
    angles = np.stack(
        [np.interp(fraction, np.linspace(0, 1, len(knots)), column) for column in knots.T], axis=1
    )
    samples = np.stack([pose(row)[2:] for row in angles])
    if not np.array_equal(samples[0], xyz[mobile]):
        raise ValueError("First pose lost source coordinate precision")
    excursions = np.linalg.norm(samples - xyz[mobile], axis=2)
    if excursions.max() > 1.65:
        raise ValueError("Permitted local displacement exceeded")
    carbonyl = samples[:, carbon - 2, :]
    toward = xyz[idx[o2]] - carbonyl
    oxygen = samples[:, ox - 2, :] - carbonyl
    approach = np.linalg.norm(toward, axis=1)
    angle = np.degrees(
        np.arccos(
            np.clip(np.sum(toward * oxygen, axis=1) / (approach * np.linalg.norm(oxygen, axis=1)), -1, 1)
        )
    )
    if not (np.all(np.diff(approach) <= 1e-6) and 2.70 <= approach[-1] <= 2.80 and 105 <= angle[-1] <= 115):
        raise ValueError("Proposed pre-attack geometry was not reached")
    geometry = dict(
        mode=mode,
        driving_bond=[c3, o3],
        assembly_source_sha256=receipt["ledger"]["source_sha256"],
        compensated_torsions=[[nodes[i], nodes[j]] for i, j, _ in joints[1:]],
        fitted_torsions_deg=knots.tolist(),
        frames=frames,
        presentation_interval=interval,
        distance_O2_to_carbonyl_C_A=[float(approach[0]), float(approach[-1])],
        angle_O2_C_O_deg=[float(angle[0]), float(angle[-1])],
        maximum_peptide_displacement_A=float(excursions.max()),
        evidence="9NO7 supplies the ester and peptide. 9MTP contributes the upstream local sugar reference. Figure 6 of doi:10.1126/science.ads9030 proposes the chemical route. The driving torsion, compensation, endpoint and time are modeled.",
        scope="A finite open approach from the reconstructed substrate. The ester remains intact; the water and its chemical identity are unchanged throughout this stage.",
    )
    return pep, samples, [[nodes[i], nodes[j]] for i, j in bonds], geometry


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("candidate", "sources", "out"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--mode", choices=("continuous", "teaching"), default="continuous")
    p.add_argument("--frames", type=int, default=49)
    args = p.parse_args()
    if args.frames < 2 or args.out.exists():
        p.error("Use at least two samples and a fresh output directory")
    receipt_path = args.candidate / "audit.json"
    receipt = json.loads(receipt_path.read_text())
    actual = audit_directory(references(args.sources), args.candidate)
    if json.loads(json.dumps(actual)) != receipt:
        raise ValueError("Candidate state no longer matches its independent source audit")
    pep, samples, bonds, geometry = build_coordinates(
        args.candidate, args.sources / "9NO7.cif", receipt, args.mode, args.frames
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=args.out.parent) as tmp:
        stage = Path(tmp) / "stage"
        stage.mkdir()
        motion = stage / "approach.npz"
        np.savez_compressed(
            motion,
            atom_ids=np.array(pep),
            positions_A=samples,
            frames=np.arange(1, args.frames + 1, dtype=float),
        )
        plan = dict(
            source_sha256=receipt["ledger"]["source_sha256"],
            assembly_id="1",
            candidate_sha256=receipt["candidate_sha256"],
            minimum_nonbonded_A=2.38,
            component_floors_A={"HOH": {"C": 2.90, "P": 2.90, "S": 2.90, "N": 2.50, "O": 2.50}},
            max_bond_change_A=0.01,
            bonds=bonds,
            stages=[
                dict(
                    id="pre_attack_approach",
                    moves=pep,
                    file=motion.name,
                    sha256=sha256(motion.read_bytes()).hexdigest(),
                    max_step_A=0.12,
                    max_entry_velocity_change_A_per_frame=0.005,
                )
            ],
        )
        plan_path = stage / "stages.json"
        plan_path.write_text(json.dumps(plan, indent=2) + "\n")
        audit = audit_motion(args.candidate / "coordinates.npz", receipt, plan, stage)
        audit["plan_sha256"] = sha256(plan_path.read_bytes()).hexdigest()
        audit["coordinate_audit_sha256"] = sha256(receipt_path.read_bytes()).hexdigest()
        (stage / "motion-audit.json").write_text(json.dumps(audit, indent=2) + "\n")
        (stage / "geometry.json").write_text(json.dumps(geometry, indent=2) + "\n")
        stage.rename(args.out)
    print(geometry["distance_O2_to_carbonyl_C_A"], geometry["angle_O2_C_O_deg"])


if __name__ == "__main__":
    main()
