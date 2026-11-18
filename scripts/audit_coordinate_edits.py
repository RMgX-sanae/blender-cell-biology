"""Compare a sparse edit with every retained atom of a reviewed mmCIF assembly.

Use Angstrom coordinates in the original assembly frame. The policy file is a
separate input, reviewed before fitting. Atom identities include the operation,
both chain IDs, author number, insertion code, component and atom name.
"""

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path

import gemmi
import numpy as np

from source_assembly import operation_products

FIELDS = ("label_asym_id", "operator", "auth_asym_id", "auth_seq_id", "insertion_code", "component", "atom")


@dataclass
class SourceAtoms:
    keys: list
    xyz: np.ndarray
    elements: list
    polymer: np.ndarray
    operators: set
    identity_only: bool
    digest: str
    assembly_id: str


def read_source(path, assembly_id, expected_sha256):
    digest = sha256(path.read_bytes()).hexdigest()
    if digest != expected_sha256:
        raise ValueError("Source SHA256 differs from the reviewed revision")
    block = gemmi.cif.read(str(path)).sole_block()
    op = block.get_mmcif_category("_pdbx_struct_oper_list.")
    matrices = {}
    for n, name in enumerate(op["id"]):
        m = np.eye(4)
        for i in range(3):
            m[i, :3] = [float(op[f"matrix[{i + 1}][{j + 1}]"][n]) for j in range(3)]
            m[i, 3] = float(op[f"vector[{i + 1}]"][n])
        if name in matrices:
            raise ValueError("Duplicate operation")
        matrices[name] = m
    gen = block.get_mmcif_category("_pdbx_struct_assembly_gen.")
    operations = {}
    for a, expr, asym in zip(gen["assembly_id"], gen["oper_expression"], gen["asym_id_list"]):
        if a != str(assembly_id):
            continue
        for name, matrix in operation_products(expr, matrices):
            for label in asym.split(","):
                if (label, name) in operations:
                    raise ValueError("Duplicate assembly instance")
                operations[label, name] = matrix
    if not operations:
        raise ValueError("Missing assembly")
    entity = block.get_mmcif_category("_entity.")
    kind = dict(zip(entity["id"], entity["type"]))
    cat = block.get_mmcif_category("_atom_site.")
    if set(cat["pdbx_PDB_model_num"]) != {"1"}:
        raise ValueError("Select one deposited model before auditing")
    # Group actual atom_site rows once, including all waters, ions and other nonpolymers.
    rows = {}
    for i, (label, el, alt, occ) in enumerate(
        zip(cat["label_asym_id"], cat["type_symbol"], cat["label_alt_id"], cat["occupancy"])
    ):
        if el.upper() not in ("H", "D") and alt in (None, False, "A") and float(occ) > 0:
            rows.setdefault(label, []).append(i)
    raw = np.column_stack([cat["Cartn_" + c] for c in "xyz"]).astype(float)
    keys, els, poly, positions = [], [], [], []
    for (label, name), matrix in operations.items():
        indices = rows.get(label, [])
        if not indices:
            continue
        positions.append(raw[indices] @ matrix[:3, :3].T + matrix[:3, 3])
        for i in indices:
            keys.append(
                (
                    label,
                    name,
                    cat["auth_asym_id"][i],
                    cat["auth_seq_id"][i],
                    cat["pdbx_PDB_ins_code"][i] or "",
                    cat["label_comp_id"][i],
                    cat["label_atom_id"][i],
                )
            )
            els.append(cat["type_symbol"][i].upper())
            poly.append(kind[cat["label_entity_id"][i]] == "polymer")
    if not keys or len(set(keys)) != len(keys):
        raise ValueError("Missing or ambiguous source identities")
    xyz = np.concatenate(positions)
    if not np.isfinite(xyz).all():
        raise ValueError("Nonfinite raw coordinates")
    identity_only = all(name == "1" and np.array_equal(m, np.eye(4)) for (_, name), m in operations.items())
    return SourceAtoms(
        keys, xyz, els, np.asarray(poly), set(operations), identity_only, digest, str(assembly_id)
    )


def atom_key(row, source):
    key = tuple(str(v) for v in row)
    if len(key) == 6 and source.identity_only:
        key = (key[0], "1", *key[1:])
    if len(key) != len(FIELDS):
        raise ValueError("Require an unambiguous assembly operation in each atom ID")
    if any(key[i] in ("", ".", "?") for i in (0, 1, 2, 3, 5, 6)):
        raise ValueError("Incomplete assembly atom identity")
    return key


def matches_source(stored, target, identity_only):
    expected = target.astype(stored.dtype)
    if identity_only or stored.dtype.itemsize == 4:
        return stored == expected
    # Independent float64 implementations of a matrix product may round its final sum differently.
    return np.isclose(stored, expected, rtol=0, atol=1e-9)


def signed_height(points):
    origin, a, b, site = points
    u, v = a - origin, b - origin
    normal = np.cross(u, v)
    length = np.linalg.norm(normal)
    if length <= 1e-6 * np.linalg.norm(u) * np.linalg.norm(v):
        raise ValueError("Cannot orient a plane from coincident or collinear anchors")
    return float((site - origin) @ (normal / length))


def compare(source, arrays, policy):
    if policy["source_sha256"] != source.digest or str(policy["assembly_id"]) != source.assembly_id:
        raise ValueError("Edit policy names a different source revision or assembly")
    ids = [atom_key(k, source) for k in arrays["atom_ids"]]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate candidate atom ID")
    indices = {k: i for i, k in enumerate(ids)}
    original = {k: i for i, k in enumerate(source.keys)}
    allowed = {}
    for row in policy["moves"]:
        key = atom_key(row["id"], source)
        limit = row["max_displacement_A"]
        if (
            key not in original
            or key in allowed
            or isinstance(limit, bool)
            or not np.isfinite(limit)
            or limit <= 0
        ):
            raise ValueError("Invalid permitted move or displacement")
        allowed[key] = limit
    owners = {k[:-1] for k in original}
    added = {}
    for row in policy["additions"]:
        key = atom_key(row["id"], source)
        el = row["element"].upper()
        if (
            key[:-1] not in owners
            or key in original
            or key in added
            or gemmi.Element(el).atomic_number <= 1
        ):
            raise ValueError("Invalid added heavy atom or source residue instance")
        added[key] = el
    if set(ids) != original.keys() | added.keys():
        raise ValueError("Candidate atom inventory differs from reviewed source plus additions")
    elems = [str(el).upper() for el in arrays["elements"]]
    if (
        len(elems) != len(ids)
        or any(elems[indices[k]] != el for k, el in zip(source.keys, source.elements))
        or any(elems[indices[k]] != el for k, el in added.items())
    ):
        raise ValueError("Element differs from its reviewed atom identity")
    reference = np.asarray(arrays["source_A"])
    model = np.asarray(arrays["model_A"])
    if (
        reference.dtype.kind != "f"
        or reference.dtype.itemsize < 4
        or model.dtype.kind != "f"
        or model.dtype.itemsize < 4
        or reference.shape != source.xyz.shape
        or model.shape != (len(ids), 3)
        or not np.isfinite(reference).all()
        or not np.isfinite(model).all()
    ):
        raise ValueError("Invalid coordinate shape or nonfinite positions")
    # A paired archive may permute rows. source_A follows the original identities in atom_ids order.
    raw_ids = [k for k in ids if k in original]
    raw_order = np.asarray([original[k] for k in raw_ids])
    paired = np.asarray([indices[k] for k in raw_ids])
    target = source.xyz[raw_order]
    if not matches_source(reference, target, source.identity_only).all():
        raise ValueError("Recorded source positions disagree with raw mmCIF identities")
    changed = ~matches_source(model[paired], target, source.identity_only).all(axis=1)
    delta = np.linalg.norm(model[paired] - target, axis=1)
    # Separate numerical roundoff of a symmetry transform from a modeled displacement.
    edit_rows = []
    for j in np.flatnonzero(changed):
        key = raw_ids[j]
        if key not in allowed:
            raise ValueError("Unapproved source atom moved: " + str(key))
        if delta[j] > allowed[key]:
            raise ValueError("Displacement exceeds the reviewed bound: " + str(key))
        edit_rows.append(
            dict(
                id=list(key),
                element=elems[paired[j]],
                source_A=target[j].tolist(),
                model_A=model[paired[j]].tolist(),
                displacement_A=float(delta[j]),
                origin="modeled reposition",
            )
        )
    for key, el in added.items():
        edit_rows.append(
            dict(
                id=list(key),
                element=el,
                source_A=None,
                model_A=model[indices[key]].tolist(),
                displacement_A=None,
                origin="modeled addition",
            )
        )
    planes = {}
    for rule in policy.get("preserve_planes", []):
        name, limit = rule["id"], rule["max_height_change_A"]
        keys = [atom_key(k, source) for k in [*rule["plane"], rule["site"]]]
        if (
            not isinstance(name, str)
            or not name
            or name in planes
            or len(keys) != 4
            or len(set(keys)) != 4
            or not set(keys) <= original.keys()
            or isinstance(limit, bool)
            or not np.isfinite(limit)
            or limit < 0
        ):
            raise ValueError("Invalid reviewed plane, original sites or tolerance")
        raw = signed_height(source.xyz[[original[k] for k in keys]])
        edited = signed_height(model[[indices[k] for k in keys]].astype(float))
        if abs(edited - raw) > limit + 1e-9:
            raise ValueError("Oriented plane changed beyond the reviewed bound: " + name)
        planes[name] = dict(
            plane=[list(k) for k in keys[:3]],
            site=list(keys[3]),
            source_signed_height_A=raw,
            model_signed_height_A=edited,
            change_A=abs(edited - raw),
            max_height_change_A=limit,
        )
    return dict(
        source_sha256=source.digest,
        assembly_id=source.assembly_id,
        identity_only=source.identity_only,
        identity_fields=FIELDS,
        source_heavy_atoms=len(original),
        retained_nonpolymer_atoms=int((~source.polymer).sum()),
        model_heavy_atoms=len(ids),
        unchanged_atoms=len(original) - int(np.count_nonzero(changed)),
        edits=edit_rows,
        preserved_planes=planes,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--candidate", type=Path, required=True)
    p.add_argument("--policy", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    policy = json.loads(a.policy.read_text())
    source = read_source(a.source, policy["assembly_id"], policy["source_sha256"])
    with np.load(a.candidate, allow_pickle=False) as arrays:
        report = compare(source, arrays, policy)
    report["candidate_sha256"] = sha256(a.candidate.read_bytes()).hexdigest()
    report["policy_sha256"] = sha256(a.policy.read_bytes()).hexdigest()
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"{report['source_heavy_atoms']} source heavy atoms, {len(report['edits'])} recorded edits")


if __name__ == "__main__":
    main()
