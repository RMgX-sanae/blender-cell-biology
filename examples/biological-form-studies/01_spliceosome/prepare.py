"""Prepare 7B9V assembly 1 for a whole C-complex view and a registered branch-site view."""

import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import urllib.request

import gemmi
import numpy as np

HERE = Path(__file__).resolve().parent
SKILL = Path(os.environ.get("CELL_BIOLOGY_SKILL", Path.home() / ".codex/skills/blender-cell-biology"))
sys.path.insert(0, str(SKILL / "scripts"))
from prepare_molecular_asset import surface

SHA = "d993bd9c0bbcf881ecafff1a831c3944c2680b014fd655727fe272fc25efce86"
URL = "https://files.rcsb.org/download/7B9V.cif.gz"


def main():
    source = HERE.parent / "sources/7B9V.cif.gz"
    source.parent.mkdir(exist_ok=True)
    if not source.exists():
        with urllib.request.urlopen(URL, timeout=120) as response:
            packed = response.read()
        if hashlib.sha256(gzip.decompress(packed)).hexdigest() != SHA:
            raise ValueError("7B9V coordinate revision changed")
        source.write_bytes(packed)
    if hashlib.sha256(gzip.open(source, "rb").read()).hexdigest() != SHA:
        raise ValueError("Unreviewed mmCIF coordinate table")
    block = gemmi.cif.read(str(source)).sole_block()
    structure = gemmi.make_structure_from_block(block)
    (assy,) = structure.assemblies
    (gen,) = assy.generators
    (op,) = gen.operators
    assert assy.name == "1" and op.transform.is_identity()
    entity = block.get_mmcif_category("_entity.")
    names = dict(zip(entity["id"], entity["pdbx_description"]))
    category = block.get_mmcif_category("_struct_asym.")
    parts = dict(zip(category["id"], category["entity_id"]))
    assert set(gen.subchains) == set(parts)
    raw, serials, groups, records = [], [], [], []
    for chain in structure[0]:
        for sub in chain.subchains():
            label = sub[0].subchain
            group = {
                "label": label,
                "author": chain.name,
                "entity": parts[label],
                "name": names[parts[label]],
                "start": len(raw),
            }
            for r in sub:
                atoms = {}
                for a in r:
                    if a.element.is_hydrogen or a.altloc not in ("\0", "A") or a.occ <= 0:
                        continue
                    atoms[a.name] = len(raw)
                    raw.append(list(a.pos))
                    serials.append(a.serial)
                if label in ("A", "B", "C", "H", "L") or r.entity_type != gemmi.EntityType.Polymer:
                    records.append(dict(label=label, auth=r.seqid.num, comp=r.name, atoms=atoms))
            group["end"] = len(raw)
            group["polymer"] = sub[0].entity_type == gemmi.EntityType.Polymer
            groups.append(group)
    raw = np.asarray(raw)
    assert len(serials) == len(set(serials)) and len(groups) == len(gen.subchains)
    origin = raw.mean(0)
    index = {(r["label"], r["auth"]): r for r in records}
    anchor = index["L", 1]["atoms"]["P"]
    target = raw[anchor]
    linkage = [index["L", 70]["atoms"]["O2'"], anchor]
    dist = np.linalg.norm(raw[linkage[0]] - target)
    assert index["L", 70]["comp"] == "A" and index["L", 1]["comp"] == "G" and 1.5 < dist < 1.7
    coords = (raw - origin) / 10.0
    arrays = {"xyz_nm": coords.astype("f4"), "atom_site_ids": np.asarray(serials, "i4")}
    for g in groups:
        if not g["polymer"] or g["label"] in ("A", "B", "C", "H", "L"):
            continue
        xyz = raw[g["start"] : g["end"]] - origin
        v, f = surface(xyz, 1.20, 1.40, 0.78)
        tri = v[f] - v.mean(0)
        if np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() < 0:
            f = f[:, [0, 2, 1]]
        arrays["v_" + g["label"]] = v / 10.0
        arrays["f_" + g["label"]] = f
    np.savez_compressed(HERE / "geometry.npz", **arrays)
    info = dict(
        pdb="7B9V",
        assembly=1,
        source_sha256=SHA,
        origin_A=origin.tolist(),
        unit="nm",
        groups=groups,
        substrand_link=linkage,
        branch_distance_A=float(dist),
        records=records,
    )
    (HERE / "inventory.json").write_text(json.dumps(info, separators=(",", ":")))
    print(f"{len(groups)} subchains, {len(raw)} mapped source atoms, branch O2'-P {dist:.4f} Å")


if __name__ == "__main__":
    main()
