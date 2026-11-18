"""Restore the terminal ribose in 9NO7 using explicitly mapped local 9MTP geometry.

All coordinates are in angstroms. No bound water is assigned a catalytic role.
"""

import argparse
import json
from tempfile import TemporaryDirectory
from pathlib import Path

import gemmi
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

from audit_coordinate_edits import FIELDS
from audit_site_links import inspect
from audit_termination_substrate import audit_directory, references
from prepare_binding_candidates import rigid_fit_points, transform
from prepare_phosphotransfer import ccd_file
from source_assembly import SourceAssembly

RING = ["O4'", "C1'", "C2'", "C3'", "C4'"]
CHIRAL = {
    "C1'": ("O4'", "C2'", "N9"),
    "C2'": ("C1'", "C3'", "O2'"),
    "C3'": ("C2'", "C4'", "O3'"),
    "C4'": ("O4'", "C3'", "C5'"),
}
TERMINAL = ("CB", "x", 76, "", "A")
PEPTIDE_END = ("DB", "z", 7, "", "CYS")
WATER = ("GQB", "A", 6983, "", "HOH", "O")
METAL = ("BQB", "x", 106, "", "MG", "MG")


def kept(atom):
    return atom.occ > 0 and not atom.element.is_hydrogen and atom.altloc in ("\x00", "A")


def site_atoms(residue):
    atoms = {a.name: np.array(list(a.pos)) for a in residue if kept(a)}
    if len(atoms) != sum(kept(a) for a in residue):
        raise ValueError("Ambiguous retained sites within a residue")
    return atoms


def local_residue(structure, label, auth, number, component):
    group = structure.groups[label]
    if group["chain"] != auth or (label, "1") not in structure.operations:
        raise ValueError("Unexpected assembly chain")
    if not np.allclose(structure.operations[label, "1"], np.eye(4), atol=1e-8, rtol=0):
        raise ValueError("This residue recipe uses coordinates from the identity operation")
    matches = [
        r
        for r in group["residues"]
        if r.seqid.num == number and not r.seqid.icode.strip() and r.name == component
    ]
    if len(matches) != 1:
        raise ValueError("Missing or ambiguous reviewed residue")
    return site_atoms(matches[0])


def assembly_atoms(structure):
    """Return all deposited positive-occupancy heavy atoms of this identity assembly."""
    if set(structure.operations) != {(label, "1") for label in structure.groups} or any(
        not np.allclose(m, np.eye(4), rtol=0, atol=1e-8) for m in structure.operations.values()
    ):
        raise ValueError("A different assembly requires explicit instance transforms and atom mapping")
    keys, elements, positions, polymer = [], [], [], []
    for label, group in structure.groups.items():
        for residue in group["residues"]:
            prefix = (label, group["chain"], residue.seqid.num, residue.seqid.icode.strip(), residue.name)
            for atom in residue:
                if kept(atom):
                    keys.append((*prefix, atom.name))
                    elements.append(atom.element.name)
                    positions.append(list(atom.pos))
                    polymer.append(residue.entity_type == gemmi.EntityType.Polymer)
    if len(set(keys)) != len(keys):
        raise ValueError("Ambiguous retained assembly atom identities")
    return keys, np.array(positions), elements, np.array(polymer)


def ring_torsion(points):
    a, b, c, d = points
    axis = (c - b) / np.linalg.norm(c - b)
    first = (a - b) - np.dot(a - b, axis) * axis
    last = (d - c) - np.dot(d - c, axis) * axis
    return np.arctan2(np.dot(np.cross(axis, first), last), np.dot(first, last))


def angles(points, triples):
    left = points[triples[:, 0]] - points[triples[:, 1]]
    right = points[triples[:, 2]] - points[triples[:, 1]]
    return np.arccos(
        np.clip(
            np.sum(left * right, axis=1) / (np.linalg.norm(left, axis=1) * np.linalg.norm(right, axis=1)),
            -1,
            1,
        )
    )


def fit_ribose(
    source,
    template,
    bonds,
    nearby,
    neighbors,
    neighbor_elements,
    carbonyl,
    peptide_ca,
    peptide_o,
    *,
    move_water,
):
    """Fit a native 2-prime hydroxyl and sugar pucker while preserving fixed covalent context."""
    if "O2'" in source or set(template) != set(source) | {"O2'"}:
        raise ValueError("Require one expressly missing hydroxyl and a complete mapped local template")
    names = list(source) + ["O2'"]
    index = {name: i for i, name in enumerate(names)}
    anchors = ["C1'", "O4'", "C4'", "C5'", "N9"]
    matrix, rmsd = rigid_fit_points(
        np.array([template[n] for n in anchors]), np.array([source[n] for n in anchors])
    )
    if rmsd > 0.35:
        raise ValueError("Local template anchor fit exceeds 0.35 angstrom")
    target = transform(np.array([template[n] for n in names]), matrix)
    raw = np.array([source[n] if n in source else target[index[n]] for n in names])
    mobile = np.array([index[n] for n in [*RING, "O2'", "O3'"]])
    edges = np.array([(index[a], index[b]) for a, b in bonds if a in index and b in index], dtype=int)
    adjacency = [set() for _ in names]
    for i, j in edges:
        adjacency[i].add(j)
        adjacency[j].add(i)
    linked = np.array([e for e in edges if any(i in mobile for i in e)])
    base_lengths = np.linalg.norm(raw[linked[:, 0]] - raw[linked[:, 1]], axis=1)
    added_bond = np.any(linked == index["O2'"], axis=1)
    if added_bond.sum() != 1 or index["C2'"] not in linked[added_bond][0]:
        raise ValueError("The CCD hydroxyl must be bonded only to C2-prime")
    base_lengths[added_bond] = np.linalg.norm(template["C2'"] - template["O2'"])
    triples = np.array(
        [
            (a, b, c)
            for b, near in enumerate(adjacency)
            for a in sorted(near)
            for c in sorted(near)
            if a < c and any(i in mobile for i in (a, b, c))
        ]
    )
    target_angles = angles(target, triples)
    quartets = np.array([[index[RING[(i + j) % 5]] for j in range(4)] for i in range(5)])
    target_torsions = np.array([ring_torsion(target[q]) for q in quartets])
    chiral_centers = np.array([index[n] for n in CHIRAL])
    chiral_neighbors = np.array([[index[n] for n in around] for around in CHIRAL.values()])

    def chiral_volumes(p):
        return np.linalg.det(p[chiral_neighbors] - p[chiral_centers, None])

    reference_volumes = chiral_volumes(target)
    self_pairs = np.array(
        [
            (a, b)
            for a in range(len(names))
            for b in range(a)
            if (a in mobile or b in mobile) and b not in adjacency[a] and not adjacency[a] & adjacency[b]
        ]
    )
    if neighbors.count(WATER) != 1:
        raise ValueError("Retain the reviewed deposited solvent site exactly once")
    water_index = neighbors.index(WATER)
    water_source = nearby[water_index].copy()
    nearby = np.delete(nearby, water_index, axis=0)
    neighbors = [k for k in neighbors if k != WATER]
    neighbor_elements = np.delete(neighbor_elements, water_index)
    water_partner_elements = np.concatenate(([n[0] for n in names], neighbor_elements))
    # Conservative contact screens for a modeled solvent oxygen; no H atoms or energies are inferred.
    water_floors = np.array(
        [
            2.90 if el in ("C", "P", "S") else 2.50 if el in ("N", "O") else 2.38
            for el in water_partner_elements
        ]
    )
    external_mask = np.ones((len(mobile), len(nearby)), dtype=bool)
    # These four pairs are the ester bond or its 1-3 neighbors across two residues.
    for row, i in enumerate(mobile):
        for col, key in enumerate(neighbors):
            if key[:5] == PEPTIDE_END and (
                (names[i] == "O3'" and key[5] in ("C", "CA", "O")) or (names[i] == "C3'" and key[5] == "C")
            ):
                external_mask[row, col] = False
    ester_context = np.stack((carbonyl, peptide_ca, peptide_o))
    ester_distances = np.linalg.norm(ester_context - source["O3'"], axis=1)
    normal = np.cross(peptide_ca - carbonyl, peptide_o - carbonyl)
    normal /= np.linalg.norm(normal)
    source_ester_height = float(normal @ (source["O3'"] - carbonyl))
    ester_angle_indices = np.array([[0, 1, 2]])

    def ester_angle(p3, po):
        return angles(np.stack((p3, po, carbonyl)), ester_angle_indices)[0]

    source_ester_angle = ester_angle(source["C3'"], source["O3'"])

    def positions(values):
        p = raw.copy()
        p[mobile] = values[: len(mobile) * 3].reshape(-1, 3)
        return p, values[-3:] if move_water else water_source

    def distances(p, water):
        return (
            np.linalg.norm(p[mobile, None] - nearby[None], axis=2)[external_mask],
            np.linalg.norm(p[self_pairs[:, 0]] - p[self_pairs[:, 1]], axis=1),
            np.linalg.norm(np.vstack((p, nearby)) - water, axis=1),
        )

    def residual(values):
        p, water = positions(values)
        external, internal, solvent = distances(p, water)
        delta = np.array([ring_torsion(p[q]) for q in quartets]) - target_torsions
        wrapped = np.arctan2(np.sin(delta), np.cos(delta))
        return np.r_[
            55 * (np.linalg.norm(p[linked[:, 0]] - p[linked[:, 1]], axis=1) - base_lengths),
            8 * (angles(p, triples) - target_angles),
            12 * wrapped,
            45 * np.maximum(np.abs(wrapped) - np.deg2rad(4.7), 0),
            5 * (chiral_volumes(p) - reference_volumes),
            35 * np.minimum(external - 2.48, 0),
            35 * np.minimum(internal - 2.48, 0),
            42 * np.minimum(solvent - water_floors - 0.025, 0),
            40 * (np.linalg.norm(ester_context - p[index["O3'"]], axis=1) - ester_distances),
            55 * (normal @ (p[index["O3'"]] - carbonyl) - source_ester_height),
            5 * (ester_angle(p[index["C3'"]], p[index["O3'"]]) - source_ester_angle),
            42 * np.maximum(abs(angles(p, triples) - target_angles) - np.deg2rad(7.5), 0),
            0.2 * (p[mobile] - target[mobile]).ravel(),
            (water - water_source) * 1.5,
            45 * max(0, np.linalg.norm(water - water_source) - 0.975),
        ]

    initial = np.r_[target[mobile].ravel(), water_source] if move_water else target[mobile].ravel()
    solution = least_squares(residual, initial, max_nfev=600, ftol=1e-10, xtol=1e-10, gtol=1e-10)
    if not solution.success:
        raise ValueError("Local substrate fit did not converge")
    fitted, water = positions(solution.x)
    ext, local, solvent = distances(fitted, water)
    new_torsions = np.array([ring_torsion(fitted[q]) for q in quartets])
    delta = new_torsions - target_torsions
    signed = {}
    for center, around in CHIRAL.items():
        idx = [index[k] for k in around]
        volume = lambda a: float(np.linalg.det(a[idx] - a[index[center]]))
        signed[center] = [volume(target), volume(fitted)]
    report = {
        "anchor_rmsd_A": rmsd,
        "proper_transform": matrix.tolist(),
        "bond_error_A": float(
            np.max(abs(np.linalg.norm(fitted[linked[:, 0]] - fitted[linked[:, 1]], axis=1) - base_lengths))
        ),
        "angle_error_deg": float(np.max(abs(np.degrees(angles(fitted, triples) - target_angles)))),
        "ring_torsion_error_deg": float(np.max(abs(np.degrees(np.arctan2(np.sin(delta), np.cos(delta)))))),
        "source_ring_torsions_deg": np.degrees([ring_torsion(raw[q]) for q in quartets]).tolist(),
        "template_ring_torsions_deg": np.degrees(target_torsions).tolist(),
        "model_ring_torsions_deg": np.degrees(new_torsions).tolist(),
        "signed_chiral_volumes_A3": signed,
        "minimum_fixed_nonbonded_A": float(ext.min()),
        "minimum_internal_nonbonded_A": float(local.min()),
        "ester_context_error_A": float(
            max(abs(np.linalg.norm(ester_context - fitted[index["O3'"]], axis=1) - ester_distances))
        ),
        "max_source_displacement_A": float(np.linalg.norm(fitted[:-1] - raw[:-1], axis=1).max()),
        "mobile_names": [names[i] for i in mobile],
    }
    new_height = float(normal @ (fitted[index["O3'"]] - carbonyl))
    report.update(
        ester_source_height_A=source_ester_height,
        ester_model_height_A=new_height,
        ester_height_change_A=abs(new_height - source_ester_height),
        ester_angle_error_deg=float(
            abs(np.degrees(ester_angle(fitted[index["C3'"]], fitted[index["O3'"]]) - source_ester_angle))
        ),
    )
    report.update(
        water_fixed=not move_water,
        water_displacement_A=float(np.linalg.norm(water - water_source)),
        water_minimum_nonbonded_A=float(solvent.min()),
        water_to_added_O2_A=float(np.linalg.norm(water - fitted[index["O2'"]])),
        water_minimum_clearance_over_element_screen_A=float((solvent - water_floors).min()),
    )
    return dict(zip(names, fitted)), water, report


def fit_checks(m):
    chiral = all(
        ref * actual > 0 and 0.7 < abs(actual / ref) < 1.3
        for ref, actual in m["signed_chiral_volumes_A3"].values()
    )
    tests = {
        "local_fit": m["anchor_rmsd_A"] <= 0.35,
        "bond_lengths": m["bond_error_A"] <= 0.03,
        "bond_angles": m["angle_error_deg"] <= 8.0,
        "pucker": m["ring_torsion_error_deg"] <= 5.0,
        "handedness": chiral,
        "unbonded_contacts": min(
            m["minimum_fixed_nonbonded_A"],
            m["minimum_internal_nonbonded_A"],
            m["water_minimum_nonbonded_A"],
        )
        >= 2.38,
        "ester_context": m["ester_context_error_A"] <= 0.04,
        "ester_plane": m["ester_height_change_A"] <= 0.04,
        "ester_angle": m["ester_angle_error_deg"] <= 8.0,
        "local_displacement": m["max_source_displacement_A"] <= 1.0,
        "water_displacement": m["water_displacement_A"] <= 1.0,
        "water_contacts_by_element": m["water_minimum_clearance_over_element_screen_A"] >= 0,
    }
    return tests


def prepare(sources, cases):
    qualified = {
        c["id"]: inspect(sources, c) for c in cases if c["id"] in ("rf1-deoxy-ester", "rf1-analogue")
    }
    if len(qualified) != 2 or not all(q["passed"] for q in qualified.values()):
        raise ValueError("Verify both source site inventories before modeling")
    original = SourceAssembly(sources / "9NO7.cif", "1")
    donor = SourceAssembly(sources / "9MTP.cif", "1")
    raw = local_residue(original, "CB", "x", 76, "A")
    analog = local_residue(donor, "CB", "1x", 76, "8AN")
    mapped = dict(analog)
    for dest, foreign in [("OP1", "O1P"), ("OP2", "O2P"), ("O3'", "N3'")]:
        mapped[dest] = mapped.pop(foreign)
    peptide = local_residue(original, "DB", "z", 7, "CYS")
    keys, source_xyz, elements, polymer = assembly_atoms(original)
    selected = [
        i
        for i in cKDTree(source_xyz).query_ball_point(raw["C2'"], 12.0)
        if keys[i][:5] != ("CB", "x", 76, "", "A")
    ]
    selected.sort()
    ccd, _ = ccd_file(sources, "A")
    graph = ccd.get_mmcif_category("_chem_comp_bond.")
    fit_args = (
        raw,
        mapped,
        list(zip(graph["atom_id_1"], graph["atom_id_2"])),
        source_xyz[selected],
        [keys[i] for i in selected],
        [elements[i] for i in selected],
        peptide["C"],
        peptide["CA"],
        peptide["O"],
    )
    fitted, water, metrics = fit_ribose(*fit_args, move_water=False)
    fixed_baseline = {"metrics": dict(metrics), "checks": fit_checks(metrics)}
    if not all(fixed_baseline["checks"].values()):
        fitted, water, metrics = fit_ribose(*fit_args, move_water=True)
    tests = fit_checks(metrics)
    if not all(tests.values()):
        raise ValueError("Rejected candidate " + str(tests) + " " + str(metrics))
    model_xyz = np.vstack((source_xyz, fitted["O2'"]))
    for i, k in enumerate(keys):
        if k[:5] == TERMINAL:
            model_xyz[i] = fitted[k[5]]
        if k == WATER:
            model_xyz[i] = water
    added = ("CB", "x", 76, "", "A", "O2'")
    changed = np.flatnonzero(np.any(source_xyz != model_xyz[:-1], axis=1))
    edits = [
        {
            "id": list(keys[i]),
            "element": elements[i],
            "origin": "modeled reposition",
            "source_A": source_xyz[i].tolist(),
            "model_A": model_xyz[i].tolist(),
        }
        for i in changed
    ]
    edits.append(
        {
            "id": list(added),
            "element": "O",
            "origin": "modeled addition",
            "source_A": None,
            "model_A": model_xyz[-1].tolist(),
        }
    )
    expected_keys = {(*TERMINAL, n) for n in metrics["mobile_names"]}
    if not metrics["water_fixed"]:
        expected_keys.add(WATER)
    if {tuple(e["id"]) for e in edits} != expected_keys:
        raise ValueError("Unexpected moved atom inventory")
    mobile_points = [p for n, p in fitted.items() if n in metrics["mobile_names"]] + [water]
    if not all(np.linalg.norm(p - raw["C2'"]) < 12.0 - 2.48 for p in mobile_points):
        raise ValueError("Local contact window would miss a possible neighbor")
    metals = [k for k in keys if k[:5] == ("BQB", "x", 106, "", "MG")]
    waters = [k for k in keys if k[:5] == ("GQB", "A", 6983, "", "HOH")]
    if len(metals) != 1 or len(waters) != 1:
        raise ValueError("Deposited active-site companion missing")
    report = {
        "sources": {
            name: {"pdb": v["source"]["id"], "sha256": v["source"]["sha256"], "passed": v["passed"]}
            for name, v in qualified.items()
        },
        "ccd_sources": {
            name: {code: record["sha256"] for code, record in v["ccd_sources"].items()}
            for name, v in qualified.items()
        },
        "row_map": "atom_ids and elements describe model_A; identities of source_A use the first len(source_A) rows. The appended O2-prime exists only in model_A.",
        "source_heavy_atoms": len(keys),
        "model_heavy_atoms": len(model_xyz),
        "unchanged_atoms": len(keys) - len(changed),
        "source_polymer_atoms": int(polymer.sum()),
        "retained_nonpolymer_atoms": int((~polymer).sum()),
        "identity_fields": list(FIELDS),
        "assembly": "9NO7 assembly 1, identity operation; one retained blank/A, positive-occupancy heavy atom per site",
        "edits": edits,
        "metrics": metrics,
        "checks": tests,
        "fixed_water_baseline": fixed_baseline,
        "nearby_source_atoms": len(selected),
        "retained_companions": {
            "phosphate_coordinating_magnesium": list(metals[0]),
            "nearby_deposited_water": list(waters[0]),
        },
        "template_use": "9MTP 8AN provides local ring, hydroxyl and adjacent bond directions; its amide N3-prime is used only as a direction for target O3-prime. No 9MTP crosslinks or partners are transferred.",
        "scope": "A candidate native ribose at one deposited deoxy site. Edits enumerate every moved or introduced source site. The ester-linked peptide, all metal sites and the rest of the ribosome remain at their experimental coordinates. Water 6983 is retained with its original identity. Its possible catalytic role and a chemical path are open.",
    }

    # This recipe established identity operation 1 for the complete assembly above.
    def instance_id(key):
        return [str(key[0]), "1", *map(str, key[1:])]

    for row in report["edits"]:
        row["id"] = instance_id(row["id"])
    report["retained_companions"] = {k: instance_id(v) for k, v in report["retained_companions"].items()}
    return [instance_id(k) for k in [*keys, added]], source_xyz, elements, model_xyz, report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sources", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    if a.out.exists():
        p.error("Choose a fresh output directory")
    cases = Path(__file__).resolve().parents[1] / "assets/site-link-examples.json"
    ids, source, elems, model, report = prepare(a.sources, json.loads(cases.read_text())["cases"])
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=a.out.parent) as tmp:
        staged = Path(tmp) / "candidate"
        staged.mkdir()
        (staged / "residue.json").write_text(json.dumps(report, indent=2) + "\n")
        np.savez_compressed(
            staged / "coordinates.npz",
            source_A=source,
            model_A=model,
            atom_ids=np.asarray(ids, dtype=str),
            elements=np.array(elems + ["O"]),
        )
        verified = audit_directory(references(a.sources), staged)
        (staged / "audit.json").write_text(json.dumps(verified, indent=2) + "\n")
        staged.rename(a.out)
    print(
        {
            k: report[k]
            for k in (
                "source_heavy_atoms",
                "model_heavy_atoms",
                "retained_nonpolymer_atoms",
                "unchanged_atoms",
            )
        },
        flush=True,
    )
    print("independent geometry", verified["checks"], flush=True)


if __name__ == "__main__":
    main()
