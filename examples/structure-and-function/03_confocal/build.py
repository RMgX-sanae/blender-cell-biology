"""A dual-disk confocal scanner with real apertures, registered lenslets and an optical return path."""

from math import atan2, cos, pi, sin, tau
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Vector

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import studio as st
from views import view


def turned(name, profile, col, mat, *, loc=(0, 0, 0), axis=(0, 0, 1), sides=96):
    """Closed lathe of a radial/axial section. Hollow profiles retain their bore."""
    verts = [
        (r * cos(tau * i / sides), r * sin(tau * i / sides), z) for r, z in profile for i in range(sides)
    ]
    faces = []
    for j in range(len(profile)):
        a, b = j * sides, ((j + 1) % len(profile)) * sides
        faces.extend((a + i, a + (i + 1) % sides, b + (i + 1) % sides, b + i) for i in range(sides))
    ob = st.mesh(name, verts, faces, col, mat)
    ob.location = loc
    ob.rotation_euler = Vector(axis).to_track_quat("Z", "Y").to_euler()
    return ob


def collar(name, r, z, w, inner, col, mat, xy=(0, 0), axis=(0, 0, 1), center=None):
    edge = min(0.040, w * 0.17, (r - inner) * 0.20)
    return turned(
        name,
        [
            (inner, -w / 2),
            (r - edge, -w / 2),
            (r, -w / 2 + edge),
            (r, w / 2 - edge),
            (r - edge, w / 2),
            (inner, w / 2),
        ],
        col,
        mat,
        loc=center if center is not None else (*xy, z),
        axis=axis,
    )


def box(name, loc, dimensions, col, mat, bevel=0.06):
    sx, sy, sz = np.array(dimensions) / 2
    verts = [(x * sx, y * sy, z * sz) for z in (-1, 1) for y in (-1, 1) for x in (-1, 1)]
    faces = [(0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4), (2, 6, 7, 3), (0, 4, 6, 2), (1, 3, 7, 5)]
    ob = st.mesh(name, verts, faces, col, mat, smooth=False)
    ob.location = loc
    if bevel:
        m = ob.modifiers.new("Machined edge radius", "BEVEL")
        m.width = bevel
        m.segments = 3
        ob.modifiers.new("Flat machined faces", "WEIGHTED_NORMAL")
    return ob


def ray(name, p, q, r1, r2, col, mat):
    start, end = Vector(p), Vector(q)
    d = end - start
    # Tube with a small end face at each end; a changing radius illustrates convergence.
    obj = turned(
        name,
        [(0.00008, 0), (r1, 0), (r2, d.length), (0.00008, d.length)],
        col,
        mat,
        loc=p,
        axis=d,
        sides=12,
    )
    obj["path_from"] = list(p)
    obj["path_to"] = list(q)
    return obj


def optical_glass(name, color, transmission, opacity):
    m = st.shader(name, color, rough=0.18, coat=0.28)
    n = m.node_tree.nodes
    l = m.node_tree.links
    p = n.get("Principled BSDF")
    p.inputs["Transmission Weight"].default_value = transmission
    p.inputs["IOR"].default_value = 1.48
    mix = n.new("ShaderNodeMixShader")
    mix.inputs[0].default_value = opacity
    l.new(n.new("ShaderNodeBsdfTransparent").outputs[0], mix.inputs[1])
    l.new(p.outputs[0], mix.inputs[2])
    l.new(mix.outputs[0], n.get("Material Output").inputs["Surface"])
    return m


def light_material(name, color, factor, strength):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    n = m.node_tree.nodes
    l = m.node_tree.links
    emit = n.new("ShaderNodeEmission")
    emit.inputs[0].default_value = (*st.linear(color), 1)
    emit.inputs[1].default_value = strength
    transparent = n.new("ShaderNodeBsdfTransparent")
    mix = n.new("ShaderNodeMixShader")
    mix.inputs[0].default_value = factor
    l.new(transparent.outputs[0], mix.inputs[1])
    l.new(emit.outputs[0], mix.inputs[2])
    l.new(mix.outputs[0], n.get("Material Output").inputs["Surface"])
    return m


def arrays(col, mats, rotor):
    nangle, nr = 36, 18
    lo, hi = 1.23, 4.07
    dr = (hi - lo) / nr
    zpin, zlens = 0.36, 3.42
    v, f, topverts, topfaces = [], [], [], []
    centers = []
    bounds = [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]
    border = []
    for a, b in zip(bounds, bounds[1:] + bounds[:1]):
        for q in (0, 0.25, 0.5, 0.75):
            border.append(tuple(u + (w - u) * q for u, w in zip(a, b)))
    ns = len(border)

    def xy(i, j, u, v):
        r = lo + (i + 0.5 + u) * dr
        theta = tau / nangle * (j + v + 0.38 * (i + 0.5 + u))
        return r * cos(theta), r * sin(theta)

    for i in range(nr):
        for j in range(nangle):
            cx, cy = xy(i, j, 0, 0)
            centers.append((cx, cy))
            cell = [xy(i, j, u, w) for u, w in border]
            angles = [atan2(y - cy, x - cx) for x, y in cell]
            start = len(v)
            for ring in ((0, 0), (0.048, 0), (0.048, -0.070), (0, -0.070)):
                radius, dz = ring
                for a, (xx, yy) in zip(angles, cell):
                    x, y = (cx + radius * cos(a), cy + radius * sin(a)) if radius else (xx, yy)
                    v.append((x, y, zpin + dz))
            for q in range(ns):
                p = (q + 1) % ns
                f.extend(
                    [
                        (start + q, start + p, start + ns + p, start + ns + q),
                        (start + ns + q, start + ns + p, start + 2 * ns + p, start + 2 * ns + q),
                        (start + 3 * ns + p, start + 3 * ns + q, start + 2 * ns + q, start + 2 * ns + p),
                    ]
                )
            # Sample the lenslet differently from the aperture. Both have exactly the same optical axis.
            start = len(topverts)
            for rad, bump in [(0.073, 0), (0.070, 0.020), (0.055, 0.038), (0.032, 0.049), (0.009, 0.053)]:
                topverts.extend(
                    (cx + rad * cos(tau * k / 16), cy + rad * sin(tau * k / 16), zlens + bump)
                    for k in range(16)
                )
            for t in range(4):
                for k in range(16):
                    q = start + t * 16 + k
                    p = start + t * 16 + (k + 1) % 16
                    topfaces.append((q, p, p + 16, q + 16))
            topfaces.append(tuple(range(start + 64, start + 80)))
    mask = st.mesh("Nipkow array   648 open pinholes", v, f, col, mats["mask"], smooth=False)
    weld = mask.modifiers.new("Continuous perforated foil", "WELD")
    weld.merge_threshold = 0.00001
    mask["aperture_radius"] = 0.048
    mask["sampled_spirals"] = nangle
    mask["sites"] = len(centers)
    lens = st.mesh("Microlens array   648 convex lenslets", topverts, topfaces, col, mats["lens"])
    lens["sites"] = len(centers)
    # Rings join the two open annulus edges and add the central mounting webs.
    for z, mat, name in [
        (zpin, mats["mask"], "Pinhole carrier"),
        (zlens, mats["substrate"], "Microlens substrate"),
    ]:
        collar(name + " hub web", lo + 0.003, z - 0.035, 0.070, 0.34, col, mat)
        collar(name + " outer lip", 4.21, z - 0.035, 0.073, hi - 0.001, col, mats["rim"])
    collar("Optical-glass annulus", hi, zlens - 0.04, 0.080, lo, col, mats["substrate"])
    for ob in col.objects:
        if ob != rotor:
            ob.parent = rotor
    return np.array(centers), zpin, zlens


def main():
    a = st.options()
    st.init()
    sc = st.scene("01  Dual spinning disk fluorescence scanner", a)
    detail = st.scene("02  Lenslets pinholes and fixed dichroic", a)
    base = st.collection("Scanner assembly", sc.collection)
    detail.collection.children.link(base)
    disk = st.collection("Two mechanically coupled rotors", base)
    mech = st.collection("Spindle motor and bearings", base)
    optics = st.collection("Fixed illumination and objective optics", base)
    ret = st.collection("Dichroic emission relay and camera", base)
    rays = st.collection("Representative light paths at the saved rotor angle", base)
    rotor = st.group("Common spindle axis", disk)
    M = {
        n: st.shader(n, c, rough=r, metal=m, coat=0.16)
        for n, c, r, m in [
            ("mask", "416178", 0.30, 0.81),
            ("rim", "A0B2BB", 0.22, 0.88),
            ("carbon", "273744", 0.30, 0.68),
            ("brass", "B5834C", 0.27, 0.76),
            ("steel", "A9BBBB", 0.22, 0.85),
            ("housing", "495C67", 0.34, 0.73),
            ("sensor", "5B606F", 0.27, 0.51),
            ("slide", "8FAEB6", 0.40, 0.36),
        ]
    }
    M["lens"] = optical_glass("Raised fused-silica lenslets", "95DFED", 0.45, 0.83)
    M["substrate"] = optical_glass("Carrier optical glass", "78B3C6", 0.58, 0.49)
    M["mirror"] = optical_glass("Dichroic coating", "B396D2", 0.32, 0.72)
    M["filter"] = optical_glass("Emission barrier filter", "96BB92", 0.38, 0.75)
    M["exc"] = light_material("Excitation light", "4FE0FC", 0.58, 1.05)
    M["em"] = light_material("Fluorescence return", "FFB14F", 0.79, 1.35)
    M["exc_cone"] = light_material("Converging beam envelope", "62DCFC", 0.12, 0.85)
    centers, zpin, zlens = arrays(disk, M, rotor)
    members = np.argsort(np.linalg.norm(centers - np.array([2.73, -0.83]), axis=1))[:3]
    cx, cy = centers[members].mean(0)
    xy = (cx, cy)
    shaft = turned(
        "Shared precision spindle",
        [(0.002, -1.62), (0.30, -1.62), (0.30, 3.67), (0.002, 3.67)],
        mech,
        M["steel"],
    )
    shaft.parent = rotor
    for z in (zpin, zlens):
        collar("Clamping hub", 0.79, z + 0.05, 0.22, 0.301, disk, M["carbon"]).parent = rotor
        collar("Machined clamp rim", 0.69, z + 0.15, 0.053, 0.309, disk, M["brass"]).parent = rotor
        for k in range(6):
            th = tau * k / 6
            pos = (0.56 * cos(th), 0.56 * sin(th), z + 0.196)
            turned(
                "Socket screw",
                [(0.030, -0.048), (0.086, -0.048), (0.086, 0), (0.032, 0)],
                disk,
                M["steel"],
                loc=pos,
                sides=6,
            ).parent = rotor
    collar("Stationary motor case", 0.96, -1.59, 2.34, 0.27, mech, M["housing"])
    for k in range(9):
        collar("Motor heat fin", 1.07, -2.43 + k * 0.185, 0.085, 0.83, mech, M["carbon"])
    for z in (-2.87, -0.43):
        collar("Motor bearing flange", 1.13, z, 0.25, 0.27, mech, M["steel"])
        collar("Bearing race", 0.54, z + 0.16, 0.12, 0.31, mech, M["brass"])
    # A stationary beam splitter in the active off-axis sector. With n=(-1,0,1)/sqrt(2),
    # returning +Z fluorescence turns into +X. Downward excitation transmits.
    zm = 1.89
    mirror = box("Fixed dichroic mirror", (cx, cy, zm), (1.74, 1.61, 0.052), ret, M["mirror"], 0.023)
    mirror.rotation_euler.y = -pi / 4
    for side in (-1, 1):
        rail = box(
            "Dichroic mount", (cx, cy + side * 0.841, zm), (1.93, 0.105, 0.085), ret, M["brass"], 0.025
        )
        rail.rotation_euler.y = -pi / 4
    for z, r, w, inner, mat in [
        (5.88, 0.80, 0.76, 0.54, "carbon"),
        (5.48, 0.72, 0.16, 0.54, "steel"),
        (6.31, 0.72, 0.14, 0.53, "brass"),
        (6.63, 0.46, 0.52, 0.24, "housing"),
        (6.99, 0.35, 0.16, 0.17, "steel"),
    ]:
        collar("Illumination collimator", r, z, w, inner, optics, M[mat], xy)
    collar("Collimating optical element", 0.58, 5.63, 0.12, 0.001, optics, M["lens"], xy)
    for z, r, w, inside, kind in [
        (-0.68, 0.82, 0.30, 0.58, "steel"),
        (-1.17, 0.94, 0.56, 0.61, "housing"),
        (-1.58, 0.97, 0.19, 0.58, "brass"),
        (-2.18, 0.84, 1.03, 0.54, "carbon"),
        (-2.81, 0.68, 0.32, 0.43, "steel"),
        (-3.11, 0.55, 0.36, 0.31, "housing"),
        (-3.35, 0.43, 0.16, 0.25, "brass"),
    ]:
        collar("Microscope objective and relay", r, z, w, inside, optics, M[kind], xy)
    for z, r in [(-0.62, 0.59), (-1.52, 0.59), (-3.38, 0.28)]:
        collar("Objective optical element", r, z, 0.10, 0.0008, optics, M["lens"], xy)
    for k in range(8):
        collar("Objective grip ring", 0.872, -2.52 + 0.106 * k, 0.037, 0.82, optics, M["housing"], xy)
    box("Specimen slide", (cx, cy, -4.12), (3.35, 1.70, 0.115), optics, M["substrate"])
    box("Cover glass", (cx, cy, -4.018), (1.50, 1.05, 0.043), optics, M["lens"], 0.015)
    for k in (-1, 1):
        box("Stage clip", (cx + 1.56 * k, cy, -4.05), (0.27, 0.65, 0.12), optics, M["steel"], 0.04)
    # The relay is open between fixed lens holders, as in an optical layout drawing.
    for x, r, th, ir, kind in [
        (5.08, 0.62, 0.18, 0.47, "brass"),
        (5.20, 0.54, 0.08, 0.001, "filter"),
        (6.02, 0.72, 0.48, 0.48, "housing"),
        (6.00, 0.51, 0.12, 0.0005, "lens"),
        (6.96, 0.68, 0.25, 0.43, "carbon"),
        (6.95, 0.48, 0.13, 0.0005, "lens"),
        (7.53, 0.73, 0.30, 0.51, "steel"),
    ]:
        collar("Emission filter or relay", r, 0, th, ir, ret, M[kind], axis=(1, 0, 0), center=(x, cy, zm))
    box("CCD camera body", (8.78, cy, zm), (1.60, 2.25, 2.12), ret, M["housing"], 0.15)
    collar(
        "Camera lens mount", 0.87, 0, 0.26, 0.54, ret, M["carbon"], axis=(1, 0, 0), center=(7.78, cy, zm)
    )
    box("Camera sensor plane", (7.91, cy, zm), (0.055, 0.83, 0.76), ret, M["sensor"], 0.015)
    for k in range(6):
        box(
            "Camera cooling ribs",
            (8.11 + k * 0.21, cy, zm + 1.04),
            (0.085, 1.83, 0.12),
            ret,
            M["carbon"],
            0.025,
        )
    mount = st.collection("Stationary optical bench and brackets", sc.collection)
    box("Dovetail base rail", (4.0, 1.39, -4.84), (12.1, 1.08, 0.43), mount, M["housing"], 0.10)
    box("Rail shoulder", (4, 1.40, -4.58), (11.8, 0.73, 0.15), mount, M["steel"], 0.045)
    for px in (-0.88, 8.75):
        box("Optical bench foot", (px, 1.4, -5.10), (2.02, 2.08, 0.25), mount, M["carbon"], 0.08)
    for px, py, top in [(0, 0.74, -2.89), (4.66, 1.10, 5.93), (8.78, 0.61, 0.89)]:
        box(
            "Fixed instrument support",
            (px, py, (-4.53 + top) / 2),
            (0.30, 0.39, top + 4.53),
            mount,
            M["housing"],
            0.056,
        )
        box("Post foot", (px, 1.04, -4.53), (0.70, 1.22, 0.22), mount, M["carbon"], 0.055)
    box("Collimator bracket", (3.73, 1.1, 5.84), (2.13, 0.28, 0.20), mount, M["housing"], 0.045)
    box("Collimator retaining arm", (cx, 0.38, 5.84), (0.24, 1.43, 0.20), mount, M["housing"], 0.04)
    box(
        "Dichroic perimeter arm",
        (3.88, cy + 0.86, zm + 0.45),
        (1.73, 0.18, 0.14),
        mount,
        M["housing"],
        0.035,
    )
    box("Dichroic arm bracket", (4.65, 0.59, zm + 0.45), (0.20, 1.18, 0.17), mount, M["housing"], 0.035)
    box("Camera cradle", (8.78, cy + 0.04, 0.79), (1.35, 2.18, 0.22), mount, M["carbon"], 0.06)
    for sign in (-1, 1):
        box(
            "Slide stage side support",
            (cx + sign * 1.54, cy + 0.43, -4.27),
            (0.29, 2.65, 0.20),
            mount,
            M["housing"],
            0.045,
        )
    box("Stage rear cross member", (cx, cy + 1.57, -4.34), (3.38, 0.43, 0.23), mount, M["carbon"], 0.07)
    for n, m in enumerate(members):
        x, y = centers[m]
        delta = x - cx
        ray(
            f"{n} excitation after collimation",
            (x, y, 5.59),
            (x, y, zlens + 0.056),
            0.030,
            0.030,
            rays,
            M["exc"],
        )
        ray(f"{n} lenslet focus envelope", (x, y, zlens), (x, y, zpin), 0.073, 0.010, rays, M["exc_cone"])
        ray(f"{n} chief excitation ray", (x, y, zlens), (x, y, zpin), 0.013, 0.013, rays, M["exc"])
        ray(f"{n} pinhole image to objective", (x, y, zpin), (x, y, -0.64), 0.016, 0.022, rays, M["exc"])
        sample = (cx - 0.40 * delta, cy - 0.40 * (y - cy), -3.986)
        lower = (cx + 0.36 * delta, cy + 0.36 * (y - cy), -3.35)
        ray(f"{n} objective internal chief ray", (x, y, -0.64), lower, 0.009, 0.010, rays, M["exc"])
        ray(f"{n} objective focal spot", lower, sample, 0.064, 0.010, rays, M["exc_cone"])
        ray(f"{n} fluorescence from specimen", sample, lower, 0.010, 0.025, rays, M["em"])
        ray(f"{n} fluorescence in objective", lower, (x, y, -0.59), 0.013, 0.011, rays, M["em"])
        ray(
            f"{n} fluorescent bead",
            sample,
            (sample[0], sample[1], sample[2] + 0.027),
            0.033,
            0.019,
            rays,
            M["em"],
        )
        # The representative return ray passes through exactly the same aperture.
        hit = (x, y, zm + delta)
        ray(f"{n} returning fluorescence through pinhole", (x, y, -0.59), hit, 0.011, 0.011, rays, M["em"])
        ray(f"{n} reflected fluorescence", hit, (6.0, y, hit[2]), 0.018, 0.018, rays, M["em"])
        ray(
            f"{n} emission relay",
            (6.0, y, hit[2]),
            (7.884, cy - 0.47 * (y - cy), zm - 0.47 * delta),
            0.015,
            0.008,
            rays,
            M["em"],
        )
    rotor["registered_sites"] = len(centers)
    rotor["scope"] = (
        "Two arrays remain registered under a shared rotation. Light paths describe the saved angle."
    )
    for s in (sc, detail):
        s["evidence"] = "Yokogawa CSU-X1 Fig. 5; optical schematic, reduced site count and designed spacing"
        s["reference_mirror_center"] = [cx, cy, zm]
        s["pinhole_plane_z"] = zpin
        s["microlens_plane_z"] = zlens
    view(sc, (2.52, 0.13, 0.90), (0.22, -1, 0.41), 21.1, ("2A3946", "101923"), roll=-0.06, light=0.60)
    view(
        detail,
        (cx + 0.08, cy + 0.20, 1.82),
        (-0.20, -1, 0.29),
        9.15,
        ("2D3E4C", "121E28"),
        roll=-0.065,
        light=0.65,
    )
    st.finish(a, [(sc, "scanner"), (detail, "optical-path")])


if __name__ == "__main__":
    main()
