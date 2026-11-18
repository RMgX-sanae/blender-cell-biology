"""Mouse mid-cochlear sensory epithelium and a registered outer-hair-cell bundle."""

from math import cos, sin, pi, tau
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import studio as st
from views import bio_material, loft, view


def hair_bundle(name, center, rotation, inner, col, mats, phase):
    x0, y0, z0 = center
    c, s = cos(rotation), sin(rotation)

    def point(x, y, z):
        return (x0 + c * x - s * y, y0 + s * x + c * y, z0 + z)

    verts, faces, material_ids = [], [], []

    def rings(levels, n, mat):
        start = len(verts)
        for cx, cy, z, radius, slope in levels:
            for k in range(n):
                th = tau * k / n
                verts.append(point(cx + radius * cos(th), cy + radius * sin(th), z + slope * cos(th)))
        for j in range(len(levels) - 1):
            for k in range(n):
                a = start + j * n + k
                b = start + j * n + (k + 1) % n
                faces.append((a, b, b + n, a + n))
        faces.extend([tuple(range(start + n - 1, start - 1, -1)), tuple(range(len(verts) - n, len(verts)))])
        material_ids.extend([mat] * (n * (len(levels) - 1) + 2))

    ncol = 27 if inner else 23
    seeds = []
    for row in range(3):
        layer = []
        for k in range(ncol):
            u = 2 * k / (ncol - 1) - 1
            y = (2.96 if inner else 2.22) * u
            sweep = 0.29 * u * u if inner else 1.90 * (u * u + 0.012) ** 0.55
            x = (0.44 if inner else 0.90) - sweep + (row - 1) * (0.35 if inner else 0.375)
            if inner:
                x += 0.09 * sin(6.8 * u + phase + row * 0.33)
            y += 0.018 * sin(11 * u + phase)
            h = (1.08, 1.72, 2.42)[row] * (1 - 0.10 * u * u) * (1 + 0.035 * sin(k * 0.72 + phase))
            r = (0.092 if not inner else 0.106) * (1 + 0.06 * cos(4 * u + row))
            lean = 0.10 * h
            levels = [
                (x, y, -0.155, r * 0.56, 0),
                (x + 0.02, y, 0.11, r * 0.66, 0),
                (x + 0.035, y, 0.25, r * 0.96, 0),
                (x + 0.55 * lean, y, 0.55 * h, r, 0),
                (x + lean, y, h - 0.14, r * 0.99, 0.035),
                (x + lean + 0.015, y, h - 0.051, r * 0.87, 0.05),
                (x + lean + 0.024, y, h - 0.012, r * 0.57, 0.027),
                (x + lean + 0.024, y, h, r * 0.20, 0.007),
            ]
            rings(levels, 12, row)
            layer.append((x, y, h, lean, r))
        seeds.append(layer)
    for row in range(2):
        for short, tall in zip(seeds[row], seeds[row + 1]):
            x, y, h, lean, r = short
            tx, ty, th, tl, tr = tall
            zz = min(th - 0.17, h + 0.24)
            # A schematic oblique tip link from the shorter tip to the taller shaft.
            xa = x + lean + 0.043
            xb = tx + tl * (zz / th) - tr * 0.77
            rings([(xa, y, h - 0.025, 0.017, 0), (xb, ty, zz, 0.017, 0)], 7, 3)
    obj = st.mesh(name, verts, faces, col, mats[0])
    for mat in mats[1:]:
        obj.data.materials.append(mat)
    obj.data.polygons.foreach_set("material_index", np.asarray(material_ids, dtype="i"))
    obj["stereocilia_per_rank"] = ncol
    obj["rank_count"] = 3
    obj["illustrated_tip_links"] = 2 * ncol
    obj["evidence"] = "Illustrated size, count and positions; adult mouse bundle profiles after SEM Fig. 4"
    return obj


def outer_cell(parent, k, j, y, mats):
    phase = 0.77 * k + 1.26 * j
    x = 7.2 + j * 7.9 + 0.22 * sin(phase)
    y += (0, 2.5, -1.6)[j] + 0.18 * cos(phase)
    z = 31.5 + j * 0.93 + 0.38 * sin(0.65 * k)
    bx, by, bz = x + 5.8, y + 3.15, 13.3 + 0.4 * j
    cell = st.collection(f"OHC {j + 1} / station {k + 1}", parent)
    sections = []
    for t, r, ry in [
        (0, 0.55, 0.65),
        (0.035, 1.47, 1.48),
        (0.09, 2.30, 2.24),
        (0.20, 2.69, 2.55),
        (0.44, 2.80, 2.57),
        (0.70, 2.73, 2.55),
        (0.91, 2.67, 2.56),
        (0.984, 2.63, 2.61),
        (1, 2.56, 2.50),
    ]:
        sections.append(
            (bx + (x - bx) * t + 0.22 * sin(pi * t), by + (y - by) * t, bz + (z - bz) * t, r, ry)
        )
    body = loft("OHC soma", sections, cell, mats["ohc"], phase=phase, ripple=0.012)
    body["row"] = j + 1
    loft(
        "Subapical cuticular plate",
        [(x, y, z - 0.28, 2.53, 2.55), (x, y, z - 0.04, 2.59, 2.61)],
        cell,
        mats["plate"],
        phase=phase,
        ripple=0.014,
    )
    hair_bundle(
        "Outer hair bundle", (x, y, z + 0.14), 0.055 * sin(phase), False, cell, mats["cilia"], phase
    )
    return cell, (x, y, z), (bx, by, bz)


def apical_mosaic(sites):
    """Bound adjacent apical cell ends by a weighted Voronoi mosaic in the reticular lamina."""
    for i, (x, y, z, weight, col, mat, name) in enumerate(sites):
        seed = np.array([x, y])
        angle = 0.13 * sin(i * 0.71)
        polygon = [
            seed + np.array([5.4 * cos(tau * k / 32 + angle), 4.6 * sin(tau * k / 32 + angle)])
            for k in range(32)
        ]
        for j, (u, v, _, other, *_) in enumerate(sites):
            if i == j or (x - u) ** 2 + (y - v) ** 2 > 150:
                continue
            axis = np.array([u, v]) - seed
            bound = (u * u + v * v - x * x - y * y + weight - other) / 2
            new = []
            for p, q in zip(polygon, polygon[1:] + polygon[:1]):
                f0, f1 = p @ axis - bound, q @ axis - bound
                if f0 <= 0:
                    new.append(p)
                if (f0 < 0) != (f1 < 0):
                    new.append(p + (q - p) * (f0 / (f0 - f1)))
            polygon = new
        if len(polygon) < 3:
            continue
        rounded = []
        for p, q in zip(polygon, polygon[1:] + polygon[:1]):
            rounded.extend((0.90 * p + 0.10 * q, 0.10 * p + 0.90 * q))
        pts = []
        faces = []
        n = len(rounded)
        for t, dz in [
            (0.53, -1.9),
            (0.84, -0.82),
            (0.985, -0.22),
            (0.994, -0.10),
            (0.969, -0.008),
            (0.76, 0.014),
            (0.45, 0.039),
        ]:
            for p in rounded:
                a, b = seed + t * (p - seed)
                pts.append(
                    (a, b, z + dz + 0.013 * (a - x) + 0.018 * sin(0.8 * (a - x) + 0.4 * (b - y) + i) * t)
                )
        for row in range(6):
            for k in range(n):
                a = row * n + k
                b = row * n + (k + 1) % n
                faces.append((a, b, b + n, a + n))
        faces.append(tuple(range(6 * n, 7 * n)))
        ob = st.mesh(name + " apical mosaic", pts, faces, col, mat)
        ob["morphology"] = "Illustrated apical cell boundary in the reticular lamina"


def main():
    a = st.options()
    st.init()
    tissue = st.scene("01  Pillar tunnel and sensory epithelium", a)
    close = st.scene("02  The three ranks of one outer hair bundle", a)
    field = st.scene("03  P21 apical field from manual annotations", a)
    # Each collection represents an intact cell or a family of supporting structures.
    hair = st.collection("Sensory cells", tissue.collection)
    support = st.collection("Pillars and Deiters cells", tissue.collection)
    base = st.collection("Basilar membrane", tissue.collection)
    mats = {
        k: bio_material(k, c, rough=r)
        for k, c, r in [
            ("ohc", "C77B61", 0.44),
            ("ihc", "A95D65", 0.44),
            ("plate", "AF795C", 0.45),
            ("dc", "427F80", 0.45),
            ("pillar", "669A9C", 0.43),
            ("head", "699998", 0.43),
            ("bm", "426167", 0.53),
        ]
    }
    mats["cilia"] = [
        bio_material(n, c, rough=0.37, sss=0.07, bump=0.009, scale=7.2)
        for n, c in [
            ("short stereocilia", "CEAD78"),
            ("middle stereocilia", "DEBF88"),
            ("tall stereocilia", "EFD6A0"),
            ("tip links", "B7EBDB"),
        ]
    ]
    # Subtle membrane roughness and pigmentation at two length scales; no additional anatomy is implied.
    for family, mat in list(mats.items())[:-1] + [("shaft", m) for m in mats["cilia"][:3]]:
        nodes = mat.node_tree.nodes
        links = mat.node_tree.links
        bsdf = nodes.get("Principled BSDF")
        base_color = np.asarray(bsdf.inputs["Base Color"].default_value[:3])
        tex = nodes.new("ShaderNodeTexCoord")
        noise = nodes.new("ShaderNodeTexNoise")
        noise.inputs["Scale"].default_value = 7.5 if family == "shaft" else 0.85
        noise.inputs["Detail"].default_value = 3.2
        noise.inputs["Roughness"].default_value = 0.60
        links.new(tex.outputs["Object"], noise.inputs["Vector"])
        ramp = nodes.new("ShaderNodeValToRGB")
        ramp.label = "Fine local color variation"
        ramp.color_ramp.elements[0].position = 0.20
        ramp.color_ramp.elements[0].color = (*(0.78 * base_color), 1)
        ramp.color_ramp.elements[1].position = 0.80
        ramp.color_ramp.elements[1].color = (*np.minimum(1, 1.09 * base_color), 1)
        links.new(noise.outputs["Fac"], ramp.inputs["Fac"])
        links.new(ramp.outputs[0], bsdf.inputs["Base Color"])
    count = 6
    # The basilar membrane is a gently relieved finite ribbon beneath this sample.
    v, f = [], []
    for h in (-0.95, -0.12):
        for yy in np.linspace(-21, 23, 11):
            for xx in np.linspace(-17.5, 35.7, 15):
                v.append((xx, yy, h + 0.085 * sin(0.13 * xx + 0.055 * yy)))
    n = 11 * 15
    for j in range(10):
        for i in range(14):
            t = j * 15 + i
            f.extend([(t + 15, t + 16, t + 1, t), (t + n, t + n + 1, t + n + 16, t + n + 15)])
    edge = (
        list(range(15))
        + [j * 15 + 14 for j in range(1, 11)]
        + list(range(163, 149, -1))
        + [j * 15 for j in range(9, 0, -1)]
    )
    f.extend((i, k, k + n, i + n) for i, k in zip(edge, edge[1:] + edge[:1]))
    st.mesh("Basilar membrane segment", v, f, base, mats["bm"])
    focus = None
    sites = []
    for k in range(count):
        y = (k - (count - 1) / 2) * 7.15
        # Medial and lateral pillar shafts meet at the top, leaving an open tunnel below.
        for name, foot, top, shaft, bottom_y, head_y in [
            ("Inner pillar", -6.9, -2.15, 1.17, 2.58, 2.87),
            ("Outer pillar", 10.2, -0.55, 1.12, 2.20, 2.75),
        ]:
            values = []
            for t, r, ry in [
                (0, 2.8, bottom_y),
                (0.05, 2.75, bottom_y),
                (0.18, 1.62, 1.46),
                (0.43, shaft, 1.20),
                (0.76, shaft * 0.85, 1.1),
                (0.93, 1.18, 1.43),
                (1, 1.68, head_y),
            ]:
                values.append((foot + (top - foot) * t + 0.36 * sin(pi * t), y, 30.7 * t, r, ry))
            loft(name, values, support, mats["pillar"], phase=k * 0.7, ripple=0.023)
        loft(
            "Joined pillar heads",
            [
                (-1.6, y, 30.45, 2.44, 2.91),
                (-1.5, y, 30.92, 3.05, 3.02),
                (-1.35, y, 31.25, 3.10, 3.1),
                (-1.35, y, 31.39, 2.93, 3.01),
            ],
            support,
            mats["head"],
        )
        sites.append((-1.5, y, 31.31, 4.7, support, mats["head"], "Pillar"))
        inner = st.collection(f"IHC / station {k + 1}", hair)
        ix, iy, iz = -9.15, y + 0.31 * cos(k), 30.65 + 0.21 * sin(k)
        loft(
            "Flask-shaped IHC",
            [
                (ix - 2.8, iy, 7.6, 0.9, 1.15),
                (ix - 2.6, iy, 8.7, 2.35, 2.25),
                (ix - 2.5, iy, 11.8, 3.83, 2.86),
                (ix - 2.25, iy, 16.8, 4.5, 3.03),
                (ix - 1.4, iy, 22.4, 3.75, 2.93),
                (ix, iy, 28.35, 2.81, 2.88),
                (ix, iy, iz, 2.66, 3.38),
            ],
            inner,
            mats["ihc"],
            ripple=0.024,
            phase=k,
        )
        loft(
            "Inner subapical cuticular plate",
            [(ix, iy, iz - 0.29, 2.45, 3.14), (ix, iy, iz - 0.12, 2.46, 3.16)],
            inner,
            mats["plate"],
            phase=k,
            ripple=0.015,
        )
        hair_bundle("Inner hair bundle", (ix, iy, iz + 0.10), -0.02, True, inner, mats["cilia"], 0.66 * k)
        sites.append((ix, iy, iz + 0.085, 16.0, inner, mats["plate"], "IHC"))
        # Inner phalangeal cells cup the hair-cell bases and run beside the IHC.
        loft(
            "Inner phalangeal cell",
            [
                (ix - 1.5, iy + 1.6, 0, 2.86, 2.68),
                (ix - 1.9, iy + 1.8, 2.2, 3.08, 2.59),
                (ix - 2.4, iy + 1.9, 6.2, 3.1, 2.65),
                (ix - 2.8, iy + 2.35, 9.0, 2.15, 1.79),
                (ix - 2.8, iy + 3.04, 15.6, 1.28, 0.82),
                (ix - 1.7, iy + 3.45, 22.7, 1.19, 0.74),
                (ix - 0.7, iy + 3.6, iz - 0.29, 1.38, 0.77),
            ],
            support,
            mats["dc"],
            ripple=0.031,
            phase=k,
        )
        sites.append((ix - 0.7, iy + 3.6, iz - 0.14, 2.0, support, mats["dc"], "Inner phalangeal"))
        for j in range(3):
            cell, (x, cy, z), (bx, by, bz) = outer_cell(hair, k, j, y, mats)
            if k == 2 and j == 1:
                close.collection.children.link(cell)
                focus = (x, cy, z)
            elif (k, j) in {(2, 2), (3, 2)}:
                close.collection.children.link(cell)
            sites.append((x, cy, z + 0.125, 8.8, cell, mats["plate"], "OHC"))
            # Soma seated in its supporting cup. Phalangeal extension leans the opposite way.
            loft(
                f"Deiters cell body {j + 1}",
                [
                    (bx + 2.1, by - 0.7, 0.02, 2.63, 2.35),
                    (bx + 2, by - 0.6, 1.3, 2.99, 2.68),
                    (bx + 1.4, by - 0.45, 5.5, 3.17, 2.85),
                    (bx + 0.5, by - 0.20, 10.9, 2.56, 2.50),
                    (bx, by, bz - 0.50, 1.97, 1.90),
                    (bx, by, bz + 0.09, 1.26, 1.34),
                ],
                support,
                mats["dc"],
                phase=k + j * 0.9,
                ripple=0.033,
            )
            xx = bx + 1.5
            loft(
                "Deiters phalangeal process",
                [
                    (xx, by, 10.2, 1.48, 1.38),
                    (xx + 0.8, by + 0.65, 13, 1.48, 1.29),
                    (xx + 0.63, by + 2.0, 17.6, 0.85, 0.92),
                    (x + 4.22, cy + 6.8, 24.4, 0.67, 0.83),
                    (x + 3.26, cy + 8.9, z - 0.65, 0.89, 1.25),
                    (x + 3.05, cy + 9.27, z + 0.02, 1.21, 2.48),
                ],
                support,
                mats["dc"],
                phase=k + j,
                ripple=0.028,
            )
            sites.append((x + 3.05, cy + 9.27, z + 0.035, 3.0, support, mats["dc"], "Deiters"))
    apical_mosaic(sites)
    close.collection.children.link(support)
    tissue["scope"] = (
        "Illustrated adult mouse mid-cochlear sensory epithelium; open field by omission of tectorial membrane and distant tissue."
    )
    close["scope"] = (
        "Same complete OHC and same mesh as tissue view. Heights, tips and link widths are illustrative."
    )
    view(
        tissue,
        (8.0, 1.3, 15.6),
        (0.36, -1.0, 0.51),
        96,
        ("293D43", "101D26"),
        roll=-0.06,
        fstop=0.47,
        light=0.49,
    )
    x, y, z = focus
    view(
        close,
        (x - 0.03, y, z + 1.44),
        (-0.96, -0.20, 0.76),
        8.42,
        ("37505A", "14232F"),
        roll=0.045,
        fstop=0.48,
        light=0.43,
    )
    source = json.loads((HERE / "field.json").read_text())
    surfaces = np.load(HERE / "field.npz")
    row_colors = {"outer 3": "8CCFC0", "outer 2": "68B9B6", "outer 1": "4FA2A8", "inner": "E6A76C"}
    for row, color in row_colors.items():
        col = st.collection(row, field.collection)
        mat = bio_material("Annotated " + row, color, rough=0.48, bump=0)
        for r in [v for v in source["bundles"] if v["display_row"] == row]:
            label = r["label"]
            o = st.tri_mesh(
                f"Manual bundle {label:02d}", surfaces[f"v{label}"], surfaces[f"f{label}"], col, mat
            )
            o["source_label"] = label
            o["source_voxels"] = r["voxels"]
            o["display_row"] = row
            o["touches_xy_boundary"] = r["touch_xy"]
            o["touches_z_boundary"] = r["touch_z"]
    field["source"] = source["specimen"]
    field["source_doi"] = source["accession"]
    field["scope"] = (
        "Measured annotation locations in a P21 apical stack, distinct from the adult mid-cochlear teaching geometry in scenes 01 and 02."
    )
    view(field, (0.1, 2.84, 0), (0.06, -0.55, 1), 85.8, ("233A45", "0D1822"), roll=0.008, light=0.48)
    st.finish(a, [(tissue, "organ-of-corti"), (close, "hair-bundle"), (field, "measured-field")])


if __name__ == "__main__":
    main()
