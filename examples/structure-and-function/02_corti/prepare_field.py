"""Retrieve one annotated P21 field and its own CZI calibration from the deposited ZIP."""

import hashlib
import io
import json
from pathlib import Path
import struct
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET
import zlib

import numpy as np
from PIL import Image
from scipy.ndimage import find_objects, gaussian_filter
from skimage.measure import marching_cubes

HERE = Path(__file__).resolve().parent
URL = "https://zenodo.org/api/records/17905676/files/VASCilia_Data.zip/content"
SPECIMEN = "Litter19Mouse3APEX-Eps8++Cdh23correctedcontro"
RANGES = {
    "p21-labels.bin": (
        217304858,
        372335,
        "398a521b2254f636fcfe32c3de0a0b9271e7a84ba605b001dfdd08aca3efca89",
    ),
    "p21-calibration.xml": (
        39368172766,
        28192,
        "3162f0344d3d0426bb6c23473c69e7d870a722ee3b32d541a44cf448915ca915",
    ),
}


def asset(name):
    start, length, sha = RANGES[name]
    path = HERE.parent / "sources" / name
    path.parent.mkdir(exist_ok=True)
    if path.is_file():
        raw = path.read_bytes()
    else:
        last = start + length - 1
        request = Request(URL + f"?range={start}-{last}", headers={"Range": f"bytes={start}-{last}"})
        with urlopen(request, timeout=100) as response:
            if (
                response.status != 206
                or response.headers["Content-Range"] != f"bytes {start}-{last}/43307707254"
            ):
                raise ValueError("Server did not return the requested archive slice")
            raw = response.read(length + 1)
    if len(raw) != length or hashlib.sha256(raw).hexdigest() != sha:
        raise ValueError("Source range changed: " + name)
    if not path.is_file():
        tmp = path.with_suffix(path.suffix + ".part")
        tmp.write_bytes(raw)
        tmp.replace(path)
    return raw


def annotations(raw):
    cursor = 0
    images = []
    for z in range(36, 57):
        head = struct.unpack_from("<4s5H3I2H", raw, cursor)
        if head[0] != b"PK\x03\x04" or head[3] != 8:
            raise ValueError("Expected a deflated local PNG member")
        name_start = cursor + 30
        member = raw[name_start : name_start + head[9]].decode()
        if member != f"VASCilia_3D_DataSet/P21/{SPECIMEN}/new_assignment_obj/{SPECIMEN}_{z:04d}.png":
            raise ValueError("Specimen or stack order changed")
        start = name_start + head[9] + head[10]
        stream = zlib.decompressobj(-15)
        png = stream.decompress(raw[start:])
        end = len(raw) - len(stream.unused_data)
        signature, crc, compressed, uncompressed = struct.unpack_from("<4I", raw, end)
        if (signature, crc, compressed, uncompressed) != (
            0x08074B50,
            zlib.crc32(png),
            end - start,
            len(png),
        ):
            raise ValueError("PNG member failed its ZIP descriptor")
        images.append(np.asarray(Image.open(io.BytesIO(png))))
        cursor = end + 16
    if cursor != len(raw):
        raise ValueError("Unexpected bytes after the final label plane")
    return np.stack(images)


def main():
    xml = ET.fromstring(asset("p21-calibration.xml"))
    pitch = {
        el.get("Id"): float(el.findtext("Value")) * 1e6 for el in xml.findall(".//Scaling/Items/Distance")
    }
    volume = annotations(asset("p21-labels.bin"))
    xy = [int(xml.findtext("Metadata/Information/Image/Size" + axis)) for axis in ("X", "Y")]
    if list(reversed(volume.shape[1:])) != xy or not np.isclose(pitch["X"], pitch["Y"]):
        raise ValueError("Annotation shape and original CZI grid differ")
    surfaces = {}
    records = []
    for label, region in enumerate(find_objects(volume), 1):
        if region is None:
            continue
        mask = volume[region] == label
        start = np.array([sl.start for sl in region])
        stop = np.array([sl.stop for sl in region])
        center = start + np.argwhere(mask).mean(0)
        # Subsample the display mesh in XY only. Source counts remain on the original grid.
        value = np.pad(mask[:, ::2, ::2], 2).astype(np.float32)
        value = gaussian_filter(value, (0.50, 0.83, 0.83))
        v, f, _, _ = marching_cubes(value, 0.5)
        k, y, x = ((v - 2) * [1, 2, 2] + start).T
        v = np.column_stack(
            ((x - (xy[0] - 1) / 2) * pitch["X"], ((xy[1] - 1) / 2 - y) * pitch["Y"], (k - 10) * pitch["Z"])
        )
        tri = v[f] - v.mean(0)
        if np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() < 0:
            f = f[:, [0, 2, 1]]
        surfaces[f"v{label}"] = v.astype(np.float32)
        surfaces[f"f{label}"] = f.astype(np.int32)
        records.append(
            dict(
                label=label,
                voxels=int(mask.sum()),
                start_zyx=start.tolist(),
                stop_zyx=stop.tolist(),
                centroid_zyx=center.tolist(),
                touch_xy=bool(np.any(start[1:] == 0) or np.any(stop[1:] == xy)),
                touch_z=bool(start[0] == 0 or stop[0] == volume.shape[0]),
            )
        )
    # Labels identify individual bundles. Row assignment is an illustrated reading of this aligned field.
    order = sorted(records, key=lambda d: d["centroid_zyx"][1])
    gaps = np.diff([r["centroid_zyx"][1] for r in order])
    cuts = sorted(np.argsort(gaps)[-3:] + 1)
    for row, subset in zip(
        ("outer 3", "outer 2", "outer 1", "inner"), np.split(np.array(order, dtype=object), cuts)
    ):
        for record in subset:
            record["display_row"] = row
    np.savez_compressed(HERE / "field.npz", **surfaces)
    (HERE / "field.json").write_text(
        json.dumps(
            dict(
                specimen=SPECIMEN,
                accession="10.5281/zenodo.17905676",
                slice_ids=[36, 56],
                source_xyz_um=pitch,
                annotation_shape_zyx=volume.shape,
                annotated_voxels=int((volume != 0).sum()),
                mesh="Half-voxel contour after 2x XY subsampling and a mild display filter; original voxel counts retained.",
                bundles=records,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
