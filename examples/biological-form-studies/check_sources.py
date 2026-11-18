"""Validate atom identities/positions and the raw CT slice using independent source reads."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import gemmi
import numpy as np
from PIL import Image

ap = argparse.ArgumentParser()
ap.add_argument("project", type=Path)
a = ap.parse_args()
root = a.project.resolve()
p = root / "01_spliceosome"
info = json.loads((p / "inventory.json").read_text())
cif = root / "sources/7B9V.cif.gz"
block = gemmi.cif.read(str(cif)).sole_block()
tab = block.get_mmcif_category("_atom_site.")
atoms = {}
for i, id_ in enumerate(tab["id"]):
    if (
        tab["type_symbol"][i] in ("H", "D")
        or float(tab["occupancy"][i]) <= 0
        or str(tab["label_alt_id"][i]) not in ("False", "None", "", "A", ".")
    ):
        continue
    atoms[int(id_)] = (
        np.array([float(tab[f"Cartn_{ax}"][i]) for ax in "xyz"]),
        tuple(
            tab[k][i]
            for k in ("label_asym_id", "auth_asym_id", "auth_seq_id", "label_comp_id", "label_atom_id")
        ),
    )
arr = np.load(p / "geometry.npz")
ids = arr["atom_site_ids"]
assert set(ids) == set(atoms)
src = np.stack([atoms[id_][0] for id_ in ids])
err = float(np.max(np.abs(src - (arr["xyz_nm"].astype(float) * 10.0 + info["origin_A"]))))
assert err < 0.00002
pair = [int(ids[k]) for k in info["substrand_link"]]
assert [atoms[k][1] for k in pair] == [("L", "I", "70", "A", "O2'"), ("L", "I", "1", "G", "P")]
dist = float(np.linalg.norm(atoms[pair[0]][0] - atoms[pair[1]][0]))
assert abs(dist - info["branch_distance_A"]) < 1e-8
p = root / "02_foraminifer"
ct = json.loads((p / "inventory.json").read_text())
raw = root / "sources/GE_GAD_2.raw"
with raw.open("rb") as f:
    digest = hashlib.file_digest(f, "sha256").hexdigest()
assert digest == "c8874ba38ca0740a7c9b78008a17224f19ecd637446de9cb75bb0c70284c0906"
a3 = np.memmap(raw, ">u2", "r", shape=(384, 765, 738))
plane = np.array(Image.open(p / "registered-slice.png"))
assert np.array_equal(plane[::-1], a3[192])
report = dict(
    pdb_sha256=hashlib.sha256(gzip.open(cif, "rb").read()).hexdigest(),
    native_atom_count=len(ids),
    native_positions_max_error_angstrom=err,
    branch_atom_site_ids=pair,
    branch_angstrom=dist,
    raw_scan_sha256=digest,
    raw_voxel_um=ct["source_voxel_um"],
    matched_native_slice_pixels=plane.size,
)
(root / "source-check.json").write_text(json.dumps(report, indent=2))
print(json.dumps(report))
