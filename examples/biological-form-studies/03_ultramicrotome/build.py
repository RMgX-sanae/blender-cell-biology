"""An editable ultramicrotome workstation detail. The two cameras share one physical assembly."""

from math import cos, sin, tau
from pathlib import Path
import sys

import numpy as np
from mathutils import Vector

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import studio as st
from views import view


def finish_edges(obj, r):
    if r:
        mod = obj.modifiers.new("Machined edge highlights", "BEVEL")
        mod.width = r
        mod.segments = 4
        obj.modifiers.new("Flat faces", "WEIGHTED_NORMAL")
    return obj


def block(name, center, dims, mat, col, r=0.04):
    x, y, z = np.asarray(dims) / 2
    p = [(a * x, b * y, c * z) for c in (-1, 1) for b in (-1, 1) for a in (-1, 1)]
    f = [(0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4), (2, 6, 7, 3), (0, 4, 6, 2), (1, 3, 7, 5)]
    ob = st.mesh(name, p, f, col, mat, smooth=False)
    ob.location = center
    return finish_edges(ob, r)


def turned(name, center, axis, profile, mat, col, knurl=False, n=96):
    verts = []
    for z, r in profile:
        for j in range(n):
            t = tau * j / n
            rr = r * (1 - 0.022 * (j % 2)) if knurl else r
            verts.append((rr * cos(t), rr * sin(t), z))
    faces = []
    for i in range(len(profile) - 1):
        for j in range(n):
            k = (j + 1) % n
            faces.append((i * n + j, i * n + k, (i + 1) * n + k, (i + 1) * n + j))
    faces.extend([tuple(reversed(range(n))), tuple(range(len(verts) - n, len(verts)))])
    ob = st.mesh(name, verts, faces, col, mat, smooth=not knurl)
    ob.location = center
    ob.rotation_euler = Vector(axis).to_track_quat("Z", "Y").to_euler()
    return finish_edges(ob, 0.012 if knurl else 0)


def xprofile(name, left, right, yz, mat, col, bevel=0):
    verts = [(x, y, z) for x in (left, right) for y, z in yz]
    n = len(yz)
    f = [tuple(reversed(range(n))), tuple(range(n, 2 * n))]
    f.extend((i, (i + 1) % n, n + (i + 1) % n, n + i) for i in range(n))
    return finish_edges(st.mesh(name, verts, f, col, mat, smooth=False), bevel)


def body_ring(name, bottom, top, corners, mat, col, bevel):
    n = len(corners)
    v = [(x, y, z) for z in (bottom, top) for x, y in corners]
    f = [tuple(reversed(range(n))), tuple(range(n, 2 * n))] + [
        (i, (i + 1) % n, n + (i + 1) % n, n + i) for i in range(n)
    ]
    return finish_edges(st.mesh(name, v, f, col, mat, smooth=False), bevel)


def optical(name, color, transmission, ior):
    m = st.shader(name, color, rough=0.16, coat=0.38)
    p = m.node_tree.nodes.get("Principled BSDF")
    p.inputs["Transmission Weight"].default_value = transmission
    p.inputs["IOR"].default_value = ior
    return m


def water_height(x, y):
    s = max(0.0, min(1.0, (y - 0.001) / 4.36))
    u = x / (0.865 + 0.63 * s**0.71)
    return 1.37 + 0.043 * abs(u) ** 12 + 0.023 * s**22 + 0.004 * sin(2 * y + 0.7 * u) * s * (1 - s)


def main():
    a = st.options()
    st.init()
    full = st.scene("01  Specimen to water boat", a)
    detail = st.scene("02  Section ribbon", a)
    assembly = st.collection("Working head", full.collection)
    detail.collection.children.link(assembly)
    slide = st.collection("Knife stage and specimen arm", assembly)
    bath = st.collection("Blue water boat and diamond insert", assembly)
    spec = st.collection("Trimmed embedded specimen", assembly)
    ribbon = st.collection("Serial sections", assembly)
    full.cycles.max_bounces = detail.cycles.max_bounces = 9
    lacquer = st.shader("Warm enamel", "CCCEC1", rough=0.32, coat=0.26)
    graphite = st.shader("Anthracite stage", "344845", rough=0.32, metal=0.52, coat=0.14)
    rail = st.shader("Ground stainless steel", "ACB9BD", rough=0.22, metal=0.83)
    satin = st.shader("Satin titanium", "777F83", rough=0.30, metal=0.75)
    dark = st.shader("Recesses and gaskets", "202C31", rough=0.46, metal=0.18)
    blue = st.shader("Blue anodized boat", "22689F", rough=0.27, metal=0.57, coat=0.24)
    edgeblue = st.shader("Blue rim", "448CC5", rough=0.25, metal=0.58)
    gold = st.shader("Diamond mounting insert", "A6854E", rough=0.24, metal=0.79)
    glass = optical("Diamond facet", "D9F7F1", 0.67, 2.41)
    water = optical("Water meniscus", "E4F7FF", 0.98, 1.333)
    water.node_tree.nodes.get("Principled BSDF").inputs["Roughness"].default_value = 0.048
    resin = optical("Embedding resin", "C99A63", 0.72, 1.53)
    tissue = st.shader("Embedded specimen cross section", "96523D", rough=0.53)
    block("Common cast base", (0, -0.56, -1.79), (7.25, 12.5, 0.86), graphite, slide, 0.24)
    block("Housing foot", (0, -5.65, -0.89), (4.85, 2.23, 0.95), lacquer, slide, 0.20)
    block("Floating stage bearing", (0, 2.12, -1.35), (5.8, 6.70, 0.29), satin, slide, 0.11)
    block("Knife stage support", (0, 2.0, -1.0), (6.3, 7.1, 0.89), graphite, slide, 0.19)
    block("Ground dovetail", (0, 2.03, -0.46), (5.35, 6.3, 0.22), rail, slide, 0.065)
    block("Sliding knife carrier", (0, 2.03, -0.21), (5.6, 6.45, 0.48), graphite, slide, 0.14)
    for x in [-2.25, 2.25]:
        block("Recessed rail", (x, 2.1, 0.039), (0.24, 5.46, 0.06), dark, slide, 0.024)
    for x, y in [(-2.2, -0.55), (2.2, -0.55), (-2.2, 4.65), (2.2, 4.65)]:
        turned(
            "Countersunk fastener",
            (x, y, 0.027),
            (0, 0, 1),
            [(0, 0.2), (0.034, 0.2), (0.047, 0.17)],
            satin,
            slide,
        )
        block("Recessed drive slot", (x, y, 0.076), (0.19, 0.038, 0.008), dark, slide, 0.01)
    turned(
        "Lateral fine adjustment",
        (2.78, 2.60, -0.16),
        (1, 0, 0),
        [(0, 0.20), (0.18, 0.2), (0.25, 0.39), (0.85, 0.39), (0.91, 0.34)],
        satin,
        slide,
        True,
    )
    block("Specimen side housing", (0, -5.70, 1.65), (4.9, 2.2, 5.0), lacquer, slide, 0.44)
    turned(
        "Specimen arm socket",
        (0, -4.65, 1.56),
        (0, 1, 0),
        [(0, 1.52), (0.18, 1.57), (0.30, 1.42)],
        graphite,
        slide,
    )
    for name, y, pr, mat, kn in [
        ("Arm", -4.6, [(0, 1.10), (0.47, 1.11), (0.69, 0.95), (1.04, 0.92)], satin, False),
        ("Collet lock", -3.63, [(0, 0.96), (0.07, 1.02), (0.58, 1.02), (0.66, 0.92)], rail, True),
        ("Chuck", -3.02, [(0, 0.92), (0.12, 0.96), (0.32, 0.91), (1.18, 0.59), (1.34, 0.56)], rail, False),
        ("Collet nose", -1.74, [(0, 0.55), (0.14, 0.57), (0.39, 0.42)], dark, False),
    ]:
        turned(name, (0, y, 1.56), (0, 1, 0), pr, mat, slide, kn)
    for i in range(4):
        ang = tau * i / 4
        ob = block(
            "Collet split",
            (0.515 * cos(ang), -1.68, 1.56 + 0.515 * sin(ang)),
            (0.037, 0.34, 0.065),
            dark,
            slide,
            0.007,
        )
        ob.rotation_euler[1] = -ang
    # A trapezoidal front face is carried through the ribbon so section width follows the block.
    verts = []
    for y, w, lo, hi in [(-1.52, 0.42, 1.19, 1.99), (-0.68, 0.44, 1.17, 1.92), (-0.028, 0.34, 1.372, 1.76)]:
        verts.extend([(-w * 0.85, y, lo), (w * 0.85, y, lo), (w, y, hi), (-w, y, hi)])
    faces = [(0, 3, 2, 1), (8, 9, 10, 11)] + [
        (j * 4 + k, j * 4 + (k + 1) % 4, (j + 1) * 4 + (k + 1) % 4, (j + 1) * 4 + k)
        for j in range(2)
        for k in range(4)
    ]
    obj = finish_edges(st.mesh("Trimmed block", verts, faces, spec, resin, smooth=False), 0.016)
    obj["face_width"] = 0.68
    shape = [
        (-0.14, 1.48),
        (-0.06, 1.435),
        (0.075, 1.448),
        (0.15, 1.52),
        (0.113, 1.615),
        (0.02, 1.66),
        (-0.103, 1.64),
        (-0.165, 1.55),
    ]
    st.mesh(
        "Tissue within the block",
        [(x, -0.031, z) for x, z in shape],
        [tuple(range(len(shape)))],
        spec,
        tissue,
        smooth=False,
    )
    corners = [
        (-1.25, -0.14),
        (1.25, -0.14),
        (1.40, 0.18),
        (1.78, 4.26),
        (1.66, 4.55),
        (-1.66, 4.55),
        (-1.78, 4.26),
        (-1.4, 0.18),
    ]
    body_ring("Blue boat base", 0.09, 0.78, corners, blue, bath, 0.12)
    left = [
        (-1.25, -0.14),
        (-0.86, -0.14),
        (-1.02, 0.28),
        (-1.40, 4.28),
        (-1.68, 4.48),
        (-1.78, 4.26),
        (-1.40, 0.18),
    ]
    for sign in (1, -1):
        body_ring("Raised bath side", 0.60, 1.53, [(x * sign, y) for x, y in left], edgeblue, bath, 0.07)
    body_ring(
        "Rear lip",
        0.65,
        1.54,
        [(-1.68, 4.48), (-1.40, 4.23), (1.40, 4.23), (1.70, 4.49), (1.59, 4.58), (-1.59, 4.58)],
        edgeblue,
        bath,
        0.055,
    )
    xprofile(
        "Clamped support for diamond",
        -0.91,
        0.91,
        [(-0.13, 0.74), (0.61, 0.75), (0.49, 1.18), (0.19, 1.32), (-0.065, 1.15)],
        gold,
        bath,
        0.024,
    )
    knife = xprofile(
        "Diamond cutting edge",
        -0.87,
        0.87,
        [(0, 1.37), (0.53, 1.285), (0.34, 1.012), (0.16, 1.18)],
        glass,
        bath,
    )
    knife["edge_xyz"] = [-0.87, 0, 1.37, 0.87, 0, 1.37]
    for sign in (-1, 1):
        xprofile(
            "Diamond clamp cheek",
            sign * 0.90,
            sign * 0.69,
            [(0.02, 1.37), (0.37, 1.40), (0.65, 1.25), (0.61, 1.12)],
            gold,
            bath,
            0.014,
        )
    # A height field closes against both side walls, ends at the lip and wets the knife edge.
    wv = []
    wf = []
    ny, nx = 74, 40
    for j in range(ny):
        s = j / (ny - 1)
        y = 0.001 + 4.36 * s
        width = 0.865 + 0.63 * s**0.71
        for i in range(nx):
            u = 2 * i / (nx - 1) - 1
            z = water_height(width * u, y)
            wv.append((width * u, y, z))
            if i and j:
                k = j * nx + i
                wf.append((k - nx - 1, k - nx, k, k - 1))
    w = st.mesh("Water surface", wv, wf, bath, water)
    w["nominal_z"] = 1.37
    solid = w.modifiers.new("Depth of the water bath", "SOLIDIFY")
    solid.thickness = 0.53
    solid.offset = -1
    film = st.shader("Floating section film", "D9D2AA", rough=0.19, metal=0.34, coat=0.40)
    p = film.node_tree.nodes.get("Principled BSDF")
    if "Thin Film Thickness" in p.inputs:
        p.inputs["Thin Film Thickness"].default_value = 355
    nodes = film.node_tree.nodes
    links = film.node_tree.links
    mix = nodes.new("ShaderNodeMixShader")
    mix.inputs[0].default_value = 0.30
    links.new(nodes.new("ShaderNodeBsdfTransparent").outputs[0], mix.inputs[1])
    links.new(p.outputs[0], mix.inputs[2])
    links.new(mix.outputs[0], nodes.get("Material Output").inputs["Surface"])
    section_mats = []
    for i, color in enumerate(["CDCDB8", "CABF9C", "B5C0BE", "CFC9AD"]):
        m = film.copy()
        m.name = "Pale film %d" % i
        m.node_tree.nodes.get("Principled BSDF").inputs["Base Color"].default_value = (*st.linear(color), 1)
        section_mats.append(m)
    # Contiguous edges form one floating ribbon. Spot outlines refer to the same block face.
    rows = []
    nsec = 8
    for j in range(nsec + 1):
        y = 0.095 + j * 0.402
        x = -0.065 + 0.075 * y * y + 0.075 * sin(1.75 * y)
        tangent = np.array([0.15 * y + 0.13125 * cos(1.75 * y), 1.0])
        tangent /= np.linalg.norm(tangent)
        cross = np.array([tangent[1], -tangent[0]])
        rows.append((np.array([x, y]), cross))
    for j in range(nsec):
        (p0, u0), (p1, u1) = rows[j : j + 2]
        xy = [p0 - u0 * 0.30, p0 + u0 * 0.30, p1 + u1 * 0.316, p1 - u1 * 0.316]
        sec = st.mesh(
            "Section %02d" % (j + 1),
            [(*v, water_height(*v) + 0.0009) for v in xy],
            [(0, 1, 2, 3)],
            ribbon,
            section_mats[j % 4],
            smooth=False,
        )
        sec["section_order"] = j + 1
    view(full, (0, -0.43, 0.35), (1.08, 1.02, 0.94), 17.85, ("354445", "151F28"), roll=-0.10, light=0.64)
    view(detail, (0.08, 0.93, 1.35), (0.42, 0.84, 1.18), 5.85, ("324247", "12212A"), roll=-0.18, light=0.62)
    for sc in [full, detail]:
        lights = [o for o in sc.objects if o.type == "LIGHT"]
        for o in lights:
            if o.name.startswith("Key"):
                o.data.shape = "RECTANGLE"
                o.data.size *= 0.57
                o.data.size_y = o.data.size * 0.31
            if o.name.startswith("Soft fill"):
                o.data.energy *= 0.48
        sc.view_settings.look = "AgX - Medium High Contrast"
    full["scope"] = detail["scope"] = (
        "Generic external working head after published instrument photographs; display proportions and film colors are illustrative."
    )
    st.finish(a, [(full, "working-head"), (detail, "ribbon")])


if __name__ == "__main__":
    main()
