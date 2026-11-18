"""Build source-registered functional-unit surfaces of squid hemocyanin 6R83."""

from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
import tempfile
import sys
from urllib.request import urlopen

import gemmi
import numpy as np

HERE = Path(__file__).resolve().parent
SKILL = Path(os.environ.get("CELL_BIOLOGY_SKILL", Path.home() / ".codex/skills/blender-cell-biology"))
sys.path.insert(0, str(SKILL / "scripts"))
from prepare_molecular_asset import surface

URL = "https://files.rcsb.org/download/6R83.cif.gz"
DIGEST = "8b1916892c2ff8474d26885555cdb2826ba6207bf93b2518ee7ff6823667687c"
DOMAINS = ("a", "b", "c", "d", "d-star", "e", "f", "g")


def source_file():
    path = HERE.parent / "sources/6R83.cif.gz"
    path.parent.mkdir(exist_ok=True)
    if not path.is_file():
        with tempfile.TemporaryDirectory(dir=path.parent) as tmp:
            download = Path(tmp) / path.name
            with urlopen(URL, timeout=120) as remote:
                download.write_bytes(remote.read())
            if hashlib.sha256(download.read_bytes()).hexdigest() != DIGEST:
                raise ValueError("Download does not match the reviewed 6R83 revision")
            download.replace(path)
    if hashlib.sha256(path.read_bytes()).hexdigest() != DIGEST:
        raise ValueError("Recheck this 6R83 source before rebuilding")
    return path


def main():
    block = gemmi.cif.read(str(source_file())).sole_block()
    structure = gemmi.make_structure_from_block(block)
    (gen,) = structure.assemblies[0].generators
    (op,) = gen.operators
    if not op.transform.is_identity() or set(gen.subchains) != set("ABCDEFGHIJ"):
        raise ValueError("Unexpected source biological assembly")
    coords = defaultdict(list)
    labels = defaultdict(list)
    for chain in structure[0]:
        for res in chain:
            # Author numbering starts a new 3314-position block for each modeled FU.
            # label_seq is the continuous sequence numbering; keep both explicitly.
            fu = DOMAINS[(res.seqid.num - 1) // 3314]
            key = res.subchain, chain.name, fu
            labels[key].append((res.label_seq, res.seqid.num))
            for atom in res:
                if atom.element.is_hydrogen or atom.occ <= 0 or atom.altloc not in ("\0", "A"):
                    continue
                coords[key].append(list(atom.pos))
    assert len(coords) == 80 and all(len(v) > 1000 for v in coords.values())
    center = np.concatenate(list(map(np.array, coords.values()))).mean(axis=0)
    data, parts = {}, []
    for (label, author, fu), points in coords.items():
        key = f"{label}_{fu}"
        xyz = np.asarray(points) - center
        vertices, faces = surface(xyz, 1.34, 1.64, 0.79)
        tri = vertices[faces] - vertices.mean(0)
        if np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() < 0:
            faces = faces[:, [0, 2, 1]]
        data[key + "_v"] = (vertices / 10).astype(np.float32)
        data[key + "_f"] = faces
        ids = labels[(label, author, fu)]
        parts.append(
            dict(
                key=key,
                label=label,
                author=author,
                fu=fu,
                atoms=len(points),
                residues=len(ids),
                label_seq=[ids[0][0], ids[-1][0]],
                auth_seq=[ids[0][1], ids[-1][1]],
                center_nm=(xyz.mean(0) / 10).tolist(),
            )
        )
    np.savez_compressed(HERE / "prepared.npz", **data)
    info = dict(
        pdb="6R83",
        assembly="1",
        sha256=DIGEST,
        origin_A=center.tolist(),
        units="nm",
        source_atoms=sum(x["atoms"] for x in parts),
        triangles=sum(len(data[x["key"] + "_f"]) for x in parts),
        description="Gaussian display envelopes of deposited coordinates; colors follow functional-unit identity.",
        parts=parts,
    )
    (HERE / "inventory.json").write_text(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
