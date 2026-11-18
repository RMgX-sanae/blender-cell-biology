"""Build a calcite display isosurface and its registered source slice from specimen GE_GAD_2."""

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from skimage.measure import marching_cubes
from PIL import Image

HERE = Path(__file__).resolve().parent
SHAPE = (384, 765, 738)
PITCH_UM = 1.6827
LEVEL = 30250.0
SHA = "c8874ba38ca0740a7c9b78008a17224f19ecd637446de9cb75bb0c70284c0906"


def largest(mask):
    labels, count = ndi.label(mask, structure=ndi.generate_binary_structure(3, 2))
    counts = np.bincount(labels.ravel())
    counts[0] = 0
    which = counts.argmax()
    return labels == which, int(counts[which]), count


def surface(field, level, pitch, origin):
    v, f, _, _ = marching_cubes(field, level, allow_degenerate=False)
    # RAW indices are Z,Y,X. Blender uses source X,Y,Z. One display unit = 100 micrometres.
    v = ((v[:, ::-1] * pitch + origin) / 100.0).astype("f4")
    tri = v[f] - v.mean(0)
    if np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() < 0:
        f = f[:, [0, 2, 1]]
    return v, np.asarray(f, dtype="i4")


def main():
    path = HERE.parent / "sources/GE_GAD_2.raw"
    if path.stat().st_size != np.prod(SHAPE) * 2:
        raise ValueError("Incomplete GE_GAD_2 big-endian scan")
    with path.open("rb") as f:
        digest = hashlib.file_digest(f, "sha256").hexdigest()
    if digest != SHA:
        raise ValueError("GE_GAD_2 scan hash differs from the reviewed specimen")
    a = np.memmap(path, dtype=">u2", shape=SHAPE, mode="r")
    # Average each acquired 2x2x2 block. Last unmatched Y row is outside the specimen.
    reduced = np.empty((192, 382, 369), dtype="f4")
    for j in range(192):
        reduced[j] = a[j * 2 : j * 2 + 2, :764].reshape(2, 382, 2, 369, 2).mean((0, 2, 4))
    signal = ndi.gaussian_filter(reduced, 0.45)
    del reduced
    shell, ns, ncomps = largest(signal > LEVEL)
    kept = ndi.binary_dilation(shell, iterations=1)
    signal[~kept] = LEVEL - 1200
    del kept
    zero = -np.array([737, 764, 383]) * PITCH_UM / 2
    origin2 = zero + 0.5 * PITCH_UM
    v, f = surface(signal, LEVEL, 2 * PITCH_UM, origin2)
    del signal
    arrays = {"shell_v": v, "shell_f": f}
    mid = SHAPE[0] // 2
    part = np.asarray(a[mid])
    # Image files start at the upper left; Blender UV (0,0) is at the bottom left.
    Image.fromarray(part[::-1].astype("uint16")).save(HERE / "registered-slice.png")
    np.savez_compressed(HERE / "geometry.npz", **arrays)
    d = dict(
        specimen="GE_GAD_2",
        taxon="Globigerinella adamsi",
        data_doi="10.1594/PANGAEA.949585",
        raw_sha256=digest,
        source_shape_zyx=SHAPE,
        source_voxel_um=PITCH_UM,
        threshold=LEVEL,
        unit_um=100.0,
        origin_um_xyz=zero.tolist(),
        raw_z_index=mid,
        selected_component_voxels=ns,
        raw_components=ncomps,
        triangles=len(arrays["shell_f"]),
    )
    (HERE / "inventory.json").write_text(json.dumps(d, indent=2))
    print("scan", digest, "vertices", len(v), "faces", len(arrays["shell_f"]), flush=True)


if __name__ == "__main__":
    main()
