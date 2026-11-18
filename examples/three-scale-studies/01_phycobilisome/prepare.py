"""Extract all deposited proteins and bilins of 6KGX in one right-handed nm frame."""

from pathlib import Path
import collections
import hashlib
import json
import os
import sys
import tempfile
import urllib.request

import gemmi
import numpy as np

SKILL = Path(os.environ.get("CELL_BIOLOGY_SKILL", Path.home() / ".codex/skills/blender-cell-biology"))
sys.path.insert(0, str(SKILL / "scripts"))
from prepare_molecular_asset import surface

SOURCE_DIR = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("SHOWCASE_PROJECT", SOURCE_DIR))
URL = "https://files.rcsb.org/download/6KGX.cif.gz"
SHA = "35f420989955fe6f316340f7f35f45cc283f5c3566a43e289def2d2f610d5b44"
src = ROOT / "sources/6KGX.cif.gz"
src.parent.mkdir(exist_ok=True)
if not src.exists() or hashlib.sha256(src.read_bytes()).hexdigest() != SHA:
    with urllib.request.urlopen(URL, timeout=150) as response:
        coordinate_file = response.read()
    if hashlib.sha256(coordinate_file).hexdigest() != SHA:
        raise ValueError("Downloaded coordinates differ from the reviewed 6KGX file")
    with tempfile.TemporaryDirectory(dir=src.parent) as staging:
        verified_file = Path(staging) / src.name
        verified_file.write_bytes(coordinate_file)
        verified_file.replace(src)
block = gemmi.cif.read(str(src)).sole_block()
model = gemmi.make_structure_from_block(block)[0]
assembly = list(
    block.find("_pdbx_struct_assembly_gen.", ["assembly_id", "oper_expression", "asym_id_list"])
)
operators = list(
    block.find(
        "_pdbx_struct_oper_list.",
        [
            "id",
            "type",
            "matrix[1][1]",
            "matrix[2][2]",
            "matrix[3][3]",
            "vector[1]",
            "vector[2]",
            "vector[3]",
        ],
    )
)
assert len(assembly) == len(operators) == 1 and assembly[0][1] == "1" and "identity" in operators[0][1]
descriptions = {
    row[0]: gemmi.cif.as_string(row[1]) for row in block.find("_entity.", ["id", "pdbx_description"])
}
parts = []
for chain in model:
    for span in chain.subchains():
        r = span[0]
        eid = int(r.entity_id)
        coords = np.array(
            [
                list(at.pos)
                for res in span
                for at in res
                if at.occ > 0 and not at.element.is_hydrogen and at.altloc in ("\x00", "A")
            ],
            dtype=float,
        )
        if not len(coords):
            continue
        role = (
            "PE"
            if eid in (1, 2)
            else "PC"
            if eid in (9, 10)
            else "APC"
            if 17 <= eid <= 20
            else "linker"
            if eid <= 25
            else r.name
        )
        parts.append(
            (
                {
                    "label_asym": r.subchain,
                    "author_chain": chain.name,
                    "entity": eid,
                    "name": descriptions[str(eid)],
                    "role": role,
                    "first_auth_res": str(r.seqid),
                    "resolved_residues": len(span),
                    "atoms": len(coords),
                },
                coords,
            )
        )
expected_asym = set(gemmi.cif.as_string(assembly[0][2]).split(","))
assert {m["label_asym"] for m, _ in parts} == expected_asym and len(parts) == 2304
counts = collections.Counter(m["role"] for m, _ in parts)
assert counts == dict(PE=528, PC=72, APC=46, linker=60, PEB=1430, PUB=48, CYC=120)
points = np.concatenate([pts for _, pts in parts])
origin = (points.max(0) + points.min(0)) / 2
c = np.cov(points[:, :2].T)
v = np.linalg.eigh(c)[1][:, -1]
if v[0] < 0:
    v = -v
basis = np.array([[v[0], -v[1], 0], [v[1], v[0], 0], [0, 0, 1]])
assert np.isclose(np.linalg.det(basis), 1)
out = ROOT / "prepared"
out.mkdir(exist_ok=True)
# A single ordered mesh array avoids filenames colliding for mixed-case mmCIF IDs.
verts = []
faces = []
records = []
voff = foff = 0
for i, (meta, pts) in enumerate(parts):
    local = (pts - origin) @ basis
    ligand = meta["entity"] > 25
    vs, fs = surface(local, 1.00 if ligand else 2.05, 1.05 if ligand else 2.15, 0.75 if ligand else 1.65)
    vs *= 0.1
    verts.append(vs)
    faces.append(fs.astype("int32"))
    meta.update(
        vertex_start=voff,
        vertex_count=len(vs),
        face_start=foff,
        face_count=len(fs),
        center_nm=(local.mean(0) * 0.1).tolist(),
    )
    records.append(meta)
    voff += len(vs)
    foff += len(fs)
vs = np.vstack(verts)
fs = np.vstack(faces)
np.savez_compressed(out / "meshes.npz", vertices=vs, faces=fs)
report = {
    "entry": "6KGX",
    "source_url": URL,
    "sha256": SHA,
    "assembly": "1, identity operation; all deposited subchains",
    "origin_A": origin.tolist(),
    "basis_columns": basis.tolist(),
    "units": "nm, 10 angstroms per unit",
    "scope": "Complete modeled positive-occupancy heavy atom sets, blank or A conformer. Gaussian illustration surfaces, not density or absolute atomic radii.",
    "counts": dict(counts),
    "heavy_atoms": len(points),
    "bounds_nm": [vs.min(0).tolist(), vs.max(0).tolist()],
    "parts": records,
}
(out / "inventory.json").write_text(json.dumps(report, ensure_ascii=False))
print(
    {k: report[k] for k in ("counts", "heavy_atoms", "bounds_nm")},
    "vertices",
    len(vs),
    "triangles",
    len(fs),
)
