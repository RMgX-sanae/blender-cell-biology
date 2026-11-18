"""Audit prepared structures against public input records, independently of Blender."""

from collections import Counter, defaultdict
import io
import json
from pathlib import Path
import struct
import sys
import xml.etree.ElementTree as ET
import zlib

import gemmi
import numpy as np
from PIL import Image


def protein(root):
    meta = json.loads((root / "01_hemocyanin/inventory.json").read_text())
    block = gemmi.cif.read(str(root / "sources/6R83.cif.gz")).sole_block()
    (gen,) = block.find("_pdbx_struct_assembly_gen.", ["assembly_id", "oper_expression", "asym_id_list"])
    assert list(gen) == ["1", "1", "A,B,C,D,E,F,G,H,I,J"]
    # Independently count atom-site rows by label sequence intervals, not by author numbering.
    edges = np.array([411, 826, 1240, 1656, 2071, 2492, 2907, 3314])
    names = ("a", "b", "c", "d", "d-star", "e", "f", "g")
    count = Counter()
    sums = defaultdict(lambda: np.zeros(3))
    fields = [
        "label_asym_id",
        "auth_asym_id",
        "label_seq_id",
        "auth_seq_id",
        "type_symbol",
        "label_alt_id",
        "occupancy",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
    ]
    for row in block.find("_atom_site.", fields):
        if row[4] in ("H", "D") or row[5] not in (".", "A") or float(row[6]) <= 0:
            continue
        fu = names[np.searchsorted(edges, int(row[2]))]
        key = (row[0], row[1], fu)
        count[key] += 1
        sums[key] += np.array(list(row)[7:], float)
    errors = []
    for r in meta["parts"]:
        key = r["label"], r["author"], r["fu"]
        assert count[key] == r["atoms"]
        errors.append(np.linalg.norm((sums[key] / count[key] - meta["origin_A"]) / 10 - r["center_nm"]))
    assert len(count) == 80 and max(errors) < 1e-8
    return dict(
        source_atom_rows=sum(count.values()), source_domains=len(count), max_centroid_error_nm=max(errors)
    )


def field(root):
    src = root / "sources"
    data = (src / "p21-labels.bin").read_bytes()
    # Use the signed ZIP data descriptors as separate end markers and verify the per-member CRC.
    ends = []
    start = 0
    while (hit := data.find(b"PK\x07\x08", start)) >= 0:
        ends.append(hit)
        start = hit + 16
    counts = Counter()
    shape = None
    first = 0
    for idx, end in enumerate(ends):
        _, crc, length, plain = struct.unpack_from("<4I", data, end)
        _, _, _, method, _, _, _, _, _, n, extra = struct.unpack_from("<4s5H3I2H", data, first)
        filename = data[first + 30 : first + 30 + n].decode()
        assert (
            filename.endswith(f"_{idx + 36:04d}.png")
            and method == 8
            and first + 30 + n + extra == end - length
        )
        png = zlib.decompress(data[end - length : end], -15)
        assert zlib.crc32(png) == crc and len(png) == plain
        image = np.asarray(Image.open(io.BytesIO(png)))
        shape = image.shape
        for val, size in enumerate(np.bincount(image.ravel())):
            if val and size:
                counts[val] += int(size)
        first = end + 16
    assert len(ends) == 21 and first == len(data)
    source = json.loads((root / "02_corti/field.json").read_text())
    assert source["annotation_shape_zyx"] == [len(ends), *shape]
    assert dict(counts) == {r["label"]: r["voxels"] for r in source["bundles"]}
    at = ET.parse(src / "p21-calibration.xml")
    um = {x.get("Id"): float(x.findtext("Value")) * 1e6 for x in at.findall(".//Scaling/Items/Distance")}
    assert um == source["source_xyz_um"] and len(counts) == 40
    return dict(
        manually_labeled_instances=len(counts),
        nonzero_voxels=sum(counts.values()),
        original_czi_xyz_um=um,
        source_png_planes=len(ends),
    )


if __name__ == "__main__":
    root = Path(sys.argv[1]).resolve()
    report = {"6R83": protein(root), "P21_manual_labels": field(root)}
    (root / "source-check.json").write_text(json.dumps(report, indent=2))
