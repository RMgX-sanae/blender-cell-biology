"""Independently audit the static A76 candidate, including its cross-residue ester and solvent."""

import argparse
from collections import defaultdict
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import gemmi
import numpy as np

from audit_coordinate_edits import FIELDS, atom_key, compare, read_source

SUGAR = ("CB", "1", "x", "76", "", "A")
DONOR = ("CB", "1", "1x", "76", "", "8AN")
CYS = ("DB", "1", "z", "7", "", "CYS")
WATER = ("GQB", "1", "A", "6983", "", "HOH", "O")
MG = ("BQB", "1", "x", "106", "", "MG", "MG")
RING = ("O4'", "C1'", "C2'", "C3'", "C4'")
HYDROXYL = (*SUGAR, "O2'")
MOBILE = {*RING, "O3'", "O2'"}
CENTRES = {
    "C1'": ("O4'", "C2'", "N9"),
    "C2'": ("C1'", "C3'", "O2'"),
    "C3'": ("C2'", "C4'", "O3'"),
    "C4'": ("O4'", "C3'", "C5'"),
}
DIRECT = {"O1P": "OP1", "O2P": "OP2", "N3'": "O3'"}


def connection(block, category, left, right):
    c = block.get_mmcif_category("_struct_conn.")
    matches = []
    for i, kind in enumerate(c["conn_type_id"]):
        if kind != category:
            continue
        ends = []
        for j in (1, 2):
            prefix = f"ptnr{j}_"
            fields = [c[prefix + f][i] for f in ["label_asym_id", "auth_asym_id", "auth_seq_id"]]
            ends.append(
                (
                    *fields,
                    c[f"pdbx_ptnr{j}_PDB_ins_code"][i] or "",
                    c[prefix + "label_comp_id"][i],
                    c[prefix + "label_atom_id"][i],
                )
            )
        if set(ends) == {left[:1] + left[2:], right[:1] + right[2:]}:
            if any(c[f"ptnr{j}_symmetry"][i] != "1_555" for j in (1, 2)):
                raise ValueError("Cross-operation connection")
            matches.append(c["id"][i])
    if len(matches) != 1:
        raise ValueError("Missing or ambiguous reviewed " + category + " connection")
    return matches[0]


def references(sources):
    config = Path(__file__).resolve().parents[1] / "assets/site-link-examples.json"
    cases = {
        c["id"]: c
        for c in json.loads(config.read_text())["cases"]
        if c["id"] in ("rf1-deoxy-ester", "rf1-analogue")
    }
    inputs = {
        case: read_source(sources / (c["source"]["pdb"] + ".cif"), "1", c["source"]["sha256"])
        for case, c in cases.items()
    }
    original = inputs["rf1-deoxy-ester"]
    donor = inputs["rf1-analogue"]
    for spec in cases.values():
        for code, digest in spec["ccd_sha256"].items():
            if sha256((sources / (code + ".ccd.cif")).read_bytes()).hexdigest() != digest:
                raise ValueError("CCD hash mismatch: " + code)
    block = gemmi.cif.read(str(sources / "9NO7.cif")).sole_block()
    analogue = gemmi.cif.read(str(sources / "9MTP.cif")).sole_block()
    links = dict(
        ester=connection(block, "covale", (*SUGAR, "O3'"), (*CYS, "C")),
        analogue_amide=connection(
            analogue, "covale", (*DONOR, "N3'"), ("DB", "1", "1z", "7", "", "CYS", "C")
        ),
        phosphate_metal=connection(block, "metalc", (*SUGAR, "OP2"), MG),
    )

    def residue(src, prefix):
        return {k[-1]: p for k, p in zip(src.keys, src.xyz) if k[:-1] == prefix}

    dep = residue(original, SUGAR)
    tmp = {DIRECT.get(k, k): v for k, v in residue(donor, DONOR).items()}
    a = gemmi.cif.read(str(sources / "A.ccd.cif")).sole_block()
    at = a.get_mmcif_category("_chem_comp_atom.")
    canonical = {k for k, el in zip(at["atom_id"], at["type_symbol"]) if el.upper() not in ("H", "D")}
    if set(dep) != canonical - {"OP3", "O2'"} or set(tmp) != canonical - {"OP3"}:
        raise ValueError("Unexpected source or donor residue inventory")
    graph = defaultdict(set)

    def edge(u, v):
        graph[u].add(v)
        graph[v].add(u)

    bonds = a.get_mmcif_category("_chem_comp_bond.")
    for u, v in zip(bonds["atom_id_1"], bonds["atom_id_2"]):
        if u in tmp and v in tmp:
            edge((*SUGAR, u), (*SUGAR, v))
    edge((*SUGAR, "O3'"), (*CYS, "C"))
    edge((*CYS, "C"), (*CYS, "CA"))
    edge((*CYS, "C"), (*CYS, "O"))
    policy = dict(
        assembly_id="1",
        source_sha256=original.digest,
        moves=[dict(id=(*SUGAR, n), max_displacement_A=1.0) for n in sorted(MOBILE - {"O2'"})]
        + [dict(id=WATER, max_displacement_A=1.0)],
        additions=[dict(id=HYDROXYL, element="O")],
        preserve_planes=[
            dict(
                id="peptide_ester",
                plane=[(*CYS, n) for n in ("C", "CA", "O")],
                site=(*SUGAR, "O3'"),
                max_height_change_A=0.04,
            )
        ],
    )
    return dict(
        source=original,
        source_sugar=dep,
        template=tmp,
        graph=graph,
        policy=policy,
        cases=cases,
        connections=links,
    )


def angle(a, b, c):
    return float(np.degrees(gemmi.calculate_angle(*(gemmi.Position(*p) for p in (a, b, c)))))


def torsion(points):
    return float(np.degrees(gemmi.calculate_dihedral(*(gemmi.Position(*p) for p in points))))


def geometry(ref, ids, xyz, els):
    ix = {k: i for i, k in enumerate(ids)}
    get = lambda key: xyz[ix[key]]
    model = {k[-1]: xyz[i] for k, i in ix.items() if k[:-1] == SUGAR}
    raw = ref["source_sugar"]
    template = ref["template"]
    graph = ref["graph"]
    lengths = []
    angles = []
    for n in model:
        neighbors = [k[-1] for k in graph[(*SUGAR, n)] if k[:-1] == SUGAR]
        for b in neighbors:
            if b <= n or not MOBILE.intersection((n, b)):
                continue
            standard = template if "O2'" in (n, b) else raw
            lengths.append(
                abs(np.linalg.norm(model[n] - model[b]) - np.linalg.norm(standard[n] - standard[b]))
            )
        for a, c in combinations(neighbors, 2):
            if MOBILE.intersection((a, n, c)):
                angles.append(
                    abs(angle(model[a], model[n], model[c]) - angle(template[a], template[n], template[c]))
                )
    if not lengths or not angles:
        raise ValueError("Missing sugar covalent constraints")
    torsions = []
    for i in range(5):
        q = [RING[(i + j) % 5] for j in range(4)]
        d = torsion([model[k] for k in q]) - torsion([template[k] for k in q])
        torsions.append(abs((d + 180) % 360 - 180))
    ratios = {}
    for n, triple in CENTRES.items():

        def signed(data):
            a, b, c = [data[v] - data[n] for v in triple]
            return np.dot(np.cross(a, b), c)

        ratios[n] = float(signed(model) / signed(template))
    carbon, ca, o = [get((*CYS, k)) for k in ("C", "CA", "O")]
    normal = np.cross(ca - carbon, o - carbon)
    normal /= np.linalg.norm(normal)
    heights = [float(normal @ (d["O3'"] - carbon)) for d in [raw, model]]
    est = max(
        abs(np.linalg.norm(raw["O3'"] - v) - np.linalg.norm(model["O3'"] - v)) for v in (carbon, ca, o)
    )
    # All model atoms are searched. Exclude only the same atom and covalent 1–2 / 1–3 neighbors.
    closest = []
    other_els = np.char.upper(np.asarray(els, dtype=str))
    for k in [(*SUGAR, n) for n in sorted(MOBILE)] + [WATER]:
        ds = np.linalg.norm(xyz - get(k), axis=1)
        near = graph.get(k, set())
        excluded = {k} | near
        for bonded in near:
            excluded = excluded | graph[bonded]
        ds[[ix[n] for n in excluded]] = np.inf
        i = int(ds.argmin())
        closest.append(dict(atom=list(k), neighbor=list(ids[i]), distance_A=float(ds[i])))
        if k == WATER:
            floor = np.where(
                np.isin(other_els, ["C", "P", "S"]),
                2.90,
                np.where(np.isin(other_els, ["N", "O"]), 2.50, 2.38),
            )
            clearance = float((ds - floor).min())
    idx = ref["source"].keys.index(WATER)
    wsrc = ref["source"].xyz[idx]
    m = dict(
        bond_error_A=float(max(lengths)),
        angle_error_deg=float(max(angles)),
        ring_torsion_error_deg=float(max(torsions)),
        ester_context_error_A=float(est),
        ester_height_change_A=abs(heights[1] - heights[0]),
        ester_source_height_A=heights[0],
        ester_model_height_A=heights[1],
        ester_angle_error_deg=abs(
            angle(raw["C3'"], raw["O3'"], carbon) - angle(model["C3'"], model["O3'"], carbon)
        ),
        water_displacement_A=float(np.linalg.norm(wsrc - get(WATER))),
        water_to_added_O2_A=float(np.linalg.norm(get(WATER) - get(HYDROXYL))),
        water_minimum_clearance_over_element_screen_A=clearance,
    )
    checks = {
        "bond_lengths": m["bond_error_A"] <= 0.03,
        "bond_angles": m["angle_error_deg"] <= 8.0,
        "ring_pucker": m["ring_torsion_error_deg"] <= 5.0,
        "chirality": all(0.7 < v < 1.3 for v in ratios.values()),
        "ester_distances": est <= 0.04,
        "ester_plane": m["ester_height_change_A"] <= 0.04,
        "ester_angle": m["ester_angle_error_deg"] <= 8.0,
        "nonbonded": min(row["distance_A"] for row in closest) >= 2.38,
        "solvent_elements": clearance >= 0,
    }
    context = dict(
        water_to_carbonyl_source_A=float(np.linalg.norm(wsrc - carbon)),
        water_to_carbonyl_model_A=float(np.linalg.norm(get(WATER) - carbon)),
        metal_to_OP2_A=float(np.linalg.norm(get(MG) - get((*SUGAR, "OP2")))),
        chiral_volume_ratios=ratios,
        nearest_nonbonded=closest,
    )
    return m, {key: bool(value) for key, value in checks.items()}, context


def audit(ref, arrays, report):
    src = ref["source"]
    ledger = compare(src, arrays, ref["policy"])
    field_count = arrays["atom_ids"].shape[1]
    fields = FIELDS if field_count == len(FIELDS) else (FIELDS[0], *FIELDS[2:])
    if tuple(report["identity_fields"]) != fields:
        raise ValueError("Reported identity columns disagree with archive")
    companions = report["retained_companions"]
    if (
        atom_key(companions["phosphate_coordinating_magnesium"], src) != MG
        or atom_key(companions["nearby_deposited_water"], src) != WATER
    ):
        raise ValueError("Reported companion differs from its source identity")
    moves = {tuple(row["id"]): row for row in ledger["edits"]}
    for k in ("source_heavy_atoms", "model_heavy_atoms", "unchanged_atoms", "retained_nonpolymer_atoms"):
        if ledger[k] != report[k]:
            raise ValueError("Reported inventory differs: " + k)
    if int(src.polymer.sum()) != report["source_polymer_atoms"]:
        raise ValueError("Polymer count differs")
    rows = {atom_key(r["id"], src): r for r in report["edits"]}
    if moves.keys() != rows.keys() or len(rows) != len(report["edits"]):
        raise ValueError("Incomplete recorded edit identities")
    for k, row in rows.items():
        for field in ("element", "origin", "source_A", "model_A"):
            if row[field] != moves[k][field]:
                raise ValueError("Edit record disagrees with coordinates: " + field)
    for case, spec in ref["cases"].items():
        if (
            report["sources"][case]["sha256"] != spec["source"]["sha256"]
            or report["sources"][case]["pdb"] != spec["source"]["pdb"]
            or report["ccd_sources"][case] != spec["ccd_sha256"]
        ):
            raise ValueError("Source or CCD provenance mismatch")
    ids = [atom_key(k, src) for k in arrays["atom_ids"]]
    m, checks, extra = geometry(ref, ids, arrays["model_A"], arrays["elements"])
    failed = [k for k, v in checks.items() if not v]
    if failed:
        raise ValueError("Independent geometry rejected: " + ", ".join(failed))
    for k, v in m.items():
        if not np.isclose(report["metrics"][k], v, atol=1e-7, rtol=0):
            raise ValueError("Generator metric disagrees with independent geometry: " + k)
    anchors = ["C1'", "O4'", "C4'", "C5'", "N9"]
    donor = np.array([ref["template"][k] for k in anchors])
    dep = np.array([ref["source_sugar"][k] for k in anchors])
    matrix = np.asarray(report["metrics"]["proper_transform"])
    if matrix.shape != (4, 4) or not np.isfinite(matrix).all():
        raise ValueError("Invalid donor registration")
    rot = matrix[:3, :3]
    residual = donor @ rot.T + matrix[:3, 3] - dep
    rms = float(np.sqrt(np.mean(np.sum(residual**2, axis=1))))
    if not (
        np.allclose(matrix[3], [0, 0, 0, 1], rtol=0, atol=1e-8)
        and np.allclose(rot.T @ rot, np.eye(3), atol=1e-8)
        and np.isclose(np.linalg.det(rot), 1)
        and rms <= 0.35
        and np.isclose(rms, report["metrics"]["anchor_rmsd_A"], atol=1e-8)
    ):
        raise ValueError("Invalid donor registration")
    return dict(
        ledger=ledger,
        source_connections=ref["connections"],
        metrics=m,
        checks=checks,
        context=extra,
        scope="Static teaching candidate. Solvent position, reconstructed hydroxyl and local rearrangements are modeled; no catalytic water or reaction trajectory was determined.",
    )


def audit_directory(ref, directory):
    candidate = directory / "coordinates.npz"
    recorded = directory / "residue.json"
    with np.load(candidate, allow_pickle=False) as arrays:
        result = audit(ref, arrays, json.loads(recorded.read_text()))
    result["candidate_sha256"] = sha256(candidate.read_bytes()).hexdigest()
    result["residue_report_sha256"] = sha256(recorded.read_bytes()).hexdigest()
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sources", type=Path, required=True)
    p.add_argument("--candidate", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    ref = references(a.sources)
    result = audit_directory(ref, a.candidate)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["metrics"]))


if __name__ == "__main__":
    main()
