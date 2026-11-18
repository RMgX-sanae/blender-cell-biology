"""Air-displacement teaching model. One common set of part meshes in three scenes."""

from pathlib import Path
import math
import os
import sys
import numpy as np
import bpy
from mathutils import Vector

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "src"))
import studio as st

a = st.options()
st.init()
exploded = st.scene("01  Exploded assembly", a)
closed = st.scene("02  Assembled pipette", a)
detail = st.scene("03  Piston and sealing stack", a)
shell = st.shader("Body | ivory polyamide", "DCE1D8", 0.34, coat=0.23)
navy = st.shader("Grip | ink-blue polymer", "183447", 0.38, coat=0.18)
dark = st.shader("Seal | graphite elastomer", "263542", 0.58)
char = st.shader("Window | blue-black enamel", "0C202F", 0.28, coat=0.2)
steel = st.shader("Piston | polished steel", "B8C7CB", 0.22, metal=0.86)
springmat = st.shader("Spring | tempered steel", "82949E", 0.26, metal=0.84)
amber = st.shader("Thumb | saffron polymer", "E5A52F", 0.29, coat=0.30)
amber2 = st.shader("Ring | dark honey", "AE6620", 0.34, metal=0.19)
pale = st.shader("Seal guide | engineering polymer", "EFE5C4", 0.37)
tipmat = st.shader("Tip | translucent polypropylene", "AECDC9", 0.22, coat=0.23)
bs = tipmat.node_tree.nodes.get("Principled BSDF")
bs.inputs["Transmission Weight"].default_value = 0.62
bs.inputs["IOR"].default_value = 1.47
tree = tipmat.node_tree
tr = tree.nodes.new("ShaderNodeBsdfTransparent")
mix = tree.nodes.new("ShaderNodeMixShader")
mix.inputs[0].default_value = 0.40
tree.links.new(bs.outputs[0], mix.inputs[1])
tree.links.new(tr.outputs[0], mix.inputs[2])
tree.links.new(mix.outputs[0], tree.nodes.get("Material Output").inputs["Surface"])
for mat in (steel, springmat):
    b = mat.node_tree.nodes.get("Principled BSDF")
    (b.inputs.get("Anisotropic") or b.inputs["Anisotropic IOR Level"]).default_value = 0.28
water = st.shader("Liquid | cyan tint", "53B7C2", 0.15, coat=0.22)
bs = water.node_tree.nodes.get("Principled BSDF")
bs.inputs["Transmission Weight"].default_value = 0.54
bs.inputs["IOR"].default_value = 1.333
ink = st.shader("Readout | ivory ink", "F3F0DB", 0.5, emission=0.11)
parts = {}


def part(k):
    col = st.collection(k, closed.collection)
    parts[k] = col
    return col


def lathe(name, prof, col, mat, n=80, knurl=0, grooves=0):
    v = []
    f = []
    for z, r in prof:
        for i in range(n):
            th = 2 * np.pi * i / n
            rr = r * (1 + knurl * np.cos(grooves * th)) if r > 0.55 else r
            v.append((rr * np.cos(th), rr * np.sin(th), z))
    for j in range(len(prof)):
        jj = (j + 1) % len(prof)
        for i in range(n):
            h = (i + 1) % n
            f.append((j * n + i, j * n + h, jj * n + h, jj * n + i))
    return st.mesh(name, v, f, col, mat)


def cuff(name, z0, z1, r, inner, col, mat):
    e = min(0.025, (z1 - z0) * 0.23, (r - inner) * 0.25)
    return lathe(
        name,
        [
            (z0, inner),
            (z0, r - e),
            (z0 + e, r),
            (z1 - e, r),
            (z1, r - e),
            (z1, inner),
            (z1 - e, inner - e * 0.2),
            (z0 + e, inner - e * 0.2),
        ],
        col,
        mat,
    )


def tube(name, pts, rad, col, mat):
    cr = bpy.data.curves.new(name, "CURVE")
    cr.dimensions = "3D"
    cr.resolution_u = 8
    cr.bevel_depth = rad
    cr.bevel_resolution = 4
    cr.use_fill_caps = True
    sp = cr.splines.new("POLY")
    sp.points.add(len(pts) - 1)
    for p, co in zip(sp.points, pts):
        p.co = (*co, 1)
    ob = bpy.data.objects.new(name, cr)
    col.objects.link(ob)
    cr.materials.append(mat)
    return ob


def panel(name, c, w, h, d, mat, col, br=0.07):
    # A bevelled rounded rectangular product insert; local y is the front.
    v = [(xx * w / 2, yy * d / 2, zz * h / 2) for zz in (-1, 1) for yy in (-1, 1) for xx in (-1, 1)]
    faces = [(0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4), (1, 3, 7, 5), (3, 2, 6, 7), (2, 0, 4, 6)]
    ob = st.mesh(name, v, faces, col, mat, False)
    ob.location = c
    bev = ob.modifiers.new("Moulded edge", "BEVEL")
    bev.width = br
    bev.segments = 5
    ob.modifiers.new("Face normals", "WEIGHTED_NORMAL")
    return ob


fontdir = Path(os.environ["SHOWCASE_FONTS"])
enfont = bpy.data.fonts.load(str(fontdir / "Times-New-Roman.ttf"))
enfont.pack()


def text_obj(name, body, loc, size, col, mat):
    d = bpy.data.curves.new(name, "FONT")
    d.body = body
    d.font = enfont
    d.align_x = "CENTER"
    d.align_y = "CENTER"
    d.size = size
    d.extrude = 0.001
    ob = bpy.data.objects.new(name, d)
    col.objects.link(ob)
    ob.location = loc
    ob.rotation_euler = (np.pi / 2, 0, 0)
    d.materials.append(mat)
    return ob


body = part("01  Handle and volume window")
# A varying superelliptical shell with a real lower and upper opening.
levels = [
    (0.60, 0.66, 0.60, 0),
    (0.66, 0.71, 0.65, 0),
    (0.85, 0.76, 0.70, 0),
    (1.35, 0.77, 0.71, 0),
    (1.60, 0.84, 0.73, 0),
    (2.0, 0.96, 0.83, 0),
    (2.65, 1.10, 0.97, -0.02),
    (3.15, 1.17, 1.04, -0.03),
    (4.55, 1.20, 1.07, -0.04),
    (6, 1.24, 1.09, -0.06),
    (7.05, 1.29, 1.14, -0.07),
    (7.54, 1.32, 1.16, -0.07),
    (8.12, 1.26, 1.10, -0.06),
    (8.59, 1.11, 0.94, 0),
    (8.89, 0.91, 0.81, 0.02),
    (9.08, 0.74, 0.70, 0.03),
    (9.18, 0.68, 0.63, 0.03),
    (9.20, 0.44, 0.42, 0.03),
    (9.10, 0.42, 0.40, 0.03),
    (0.61, 0.41, 0.40, 0),
]
v = []
f = []
idx = []
n = 96
for z, rx, ry, dx in levels:
    for i in range(n):
        th = 2 * np.pi * i / n
        cx = np.cos(th)
        cy = np.sin(th)
        v.append((dx + rx * np.sign(cx) * abs(cx) ** 0.78, ry * np.sign(cy) * abs(cy) ** 0.78, z))
for j in range(len(levels)):
    z = levels[j][0]
    for i in range(n):
        ii = (i + 1) % n
        jj = (j + 1) % len(levels)
        f.append((j * n + i, j * n + ii, jj * n + ii, jj * n + i))
        th = 2 * np.pi * (i + 0.5) / n
        idx.append(int(j <= 15 and ((np.sin(th) > -0.50 and 1.95 <= z <= 8.15) or z < 1.60)))
handle = st.mesh("Ergonomic hollow housing", v, f, body, shell)
handle.data.materials.append(navy)
handle.data.polygons.foreach_set("material_index", np.array(idx, dtype="i"))
# Sculpted finger rest across the back/side with tapered, rounded ends.
v = []
f = []
n = 48
for x, cy, z, ry, rz in [
    (-0.58, 0.22, 7.95, 0.70, 0.50),
    (-1.03, 0.20, 7.99, 0.73, 0.49),
    (-1.47, 0.18, 8.03, 0.70, 0.42),
    (-1.95, 0.15, 8.14, 0.62, 0.31),
    (-2.32, 0.13, 8.24, 0.50, 0.22),
    (-2.52, 0.12, 8.25, 0.32, 0.145),
    (-2.63, 0.12, 8.25, 0.035, 0.021),
]:
    for i in range(n):
        t = 2 * np.pi * i / n
        v.append((x, cy + ry * np.cos(t), z + rz * np.sin(t)))
for j in range(6):
    for i in range(n):
        ii = (i + 1) % n
        f.append((j * n + i, j * n + ii, (j + 1) * n + ii, (j + 1) * n + i))
f.append(tuple(reversed(range(n))))
f.append(tuple(6 * n + i for i in range(n)))
st.mesh("Tapered finger rest", v, f, body, navy)
panel("Recessed bezel", (-0.04, -1.055, 5.54), 0.86, 2.06, 0.14, navy, body, 0.068)
panel("Mechanical counter glass", (-0.04, -1.145, 5.56), 0.60, 1.52, 0.065, char, body, 0.03)
for z, t in [(6.04, "1"), (5.56, "0"), (5.08, "0")]:
    text_obj("Counter " + str(z), t, (-0.04, -1.19, z), 0.465, body, ink)
for z in (5.81, 5.32):
    panel("Counter drum separation", (-0.04, -1.186, z), 0.47, 0.018, 0.012, steel, body, 0.003)
text_obj("Microlitre unit", "µL", (-0.04, -1.113, 4.21), 0.28, body, navy)
for i in range(4):
    z = 3.17 + i * 0.18
    tube(
        "Side grip fluting %d" % (i + 1),
        [(0.74, -0.87, z), (1.0, -0.63, z + 0.04), (1.13, -0.27, z + 0.07), (1.17, 0.12, z + 0.11)],
        0.026,
        body,
        navy,
    )
thumb = part("02  Thumb plunger and setting dial")
cuff("Piston actuation spindle", 8.6, 10.94, 0.30, 0.04, thumb, steel)
cuff("Dial spindle cover", 9.18, 9.92, 0.48, 0.34, thumb, navy)
lathe(
    "Volume adjustment knurl",
    [(9.78, 0.32), (9.78, 0.56), (9.84, 0.61), (10.20, 0.61), (10.26, 0.55), (10.26, 0.32)],
    thumb,
    navy,
    144,
    0.024,
    36,
)
cuff("Dial brass index", 10.16, 10.22, 0.617, 0.52, thumb, amber)
lathe(
    "Concave thumb button",
    [
        (10.52, 0.09),
        (10.52, 0.55),
        (10.60, 0.71),
        (10.76, 0.78),
        (11.04, 0.78),
        (11.13, 0.70),
        (11.16, 0.55),
        (11.09, 0.19),
        (11.085, 0.01),
        (10.54, 0.01),
    ],
    thumb,
    amber,
)
ej = part("03  Independent tip-ejector lever")
panel("Eject key linkage", (0.86, -0.08, 8.90), 0.39, 0.95, 0.32, steel, ej, 0.12)
panel("Eject key", (0.92, -0.05, 9.40), 0.72, 0.33, 0.91, navy, ej, 0.13)
panel("Eject key contact pad", (0.93, -0.07, 9.56), 0.56, 0.13, 0.71, amber, ej, 0.052)
nut = part("04  Connecting nut")
lathe(
    "Hollow connecting nut",
    [
        (0.30, 0.47),
        (0.30, 0.67),
        (0.37, 0.82),
        (0.50, 0.86),
        (1.34, 0.84),
        (1.46, 0.76),
        (1.49, 0.53),
        (1.40, 0.49),
    ],
    nut,
    navy,
    180,
    0.028,
    36,
)
cuff("Nut retaining lip", 0.37, 0.43, 0.835, 0.62, nut, steel)
pist = part("05  Piston and return spring")
lathe(
    "Continuous polished piston",
    [
        (-4.50, 0.018),
        (-4.50, 0.143),
        (-4.46, 0.164),
        (-1.30, 0.165),
        (-1.22, 0.22),
        (-1.14, 0.23),
        (-0.98, 0.25),
        (-0.90, 0.21),
        (1.60, 0.21),
        (1.72, 0.16),
        (1.72, 0.012),
    ],
    pist,
    steel,
)
for z in (-1.18, 1.05):
    cuff("Spring seat", z, z + 0.15, 0.46, 0.185, pist, steel)
turn = np.linspace(0, 2 * np.pi * 9.6, 385)
tube(
    "Helical compression spring",
    [(0.335 * np.cos(t), 0.335 * np.sin(t), -1.01 + 2.07 * t / turn[-1]) for t in turn],
    0.058,
    pist,
    springmat,
)
cuff("Spring carrier collar", 1.18, 1.47, 0.36, 0.195, pist, navy)
seal = part("06  Piston seal guide")
lathe(
    "Annular seal guide",
    [
        (-1.59, 0.164),
        (-1.59, 0.39),
        (-1.53, 0.425),
        (-1.37, 0.425),
        (-1.34, 0.376),
        (-1.26, 0.376),
        (-1.24, 0.25),
        (-1.27, 0.164),
    ],
    seal,
    pale,
)
oring = part("07  O-ring")
ang = np.linspace(0, 2 * np.pi, 25)[:-1]
lathe(
    "Elastomer O-ring", [(-1.71 + 0.125 * np.sin(t), 0.290 + 0.125 * np.cos(t)) for t in ang], oring, dark
)
shaft = part("08  Hollow tip holder")
lathe(
    "Tapered air channel and seal seat",
    [
        (-6.68, 0.084),
        (-6.68, 0.154),
        (-6.61, 0.187),
        (-6.00, 0.236),
        (-5.74, 0.262),
        (-5.71, 0.315),
        (-5.59, 0.321),
        (-5.55, 0.289),
        (-4.72, 0.31),
        (-2.49, 0.425),
        (-2.10, 0.515),
        (-1.74, 0.60),
        (-0.92, 0.63),
        (-0.68, 0.69),
        (0.13, 0.69),
        (0.43, 0.615),
        (0.53, 0.58),
        (0.53, 0.45),
        (-1.25, 0.45),
        (-1.79, 0.414),
        (-1.87, 0.27),
        (-4.64, 0.245),
        (-5.83, 0.138),
        (-6.64, 0.083),
    ],
    shaft,
    shell,
)
for z in (-0.63, -0.50, -0.37):
    cuff("Tip-holder engagement ridge", z, z + 0.052, 0.705, 0.58, shaft, pale)
tip = part("09  Disposable hollow tip")
lathe(
    "Fine open polypropylene tip",
    [
        (-11.17, 0.027),
        (-11.17, 0.038),
        (-10.97, 0.064),
        (-9.0, 0.146),
        (-7.43, 0.245),
        (-6.25, 0.336),
        (-5.65, 0.368),
        (-5.60, 0.39),
        (-5.52, 0.39),
        (-5.51, 0.324),
        (-6.22, 0.293),
        (-7.52, 0.211),
        (-9.03, 0.126),
        (-10.99, 0.036),
    ],
    tip,
    tipmat,
)
ejector = part("10  External tip ejector")
# A stamped rail is offset from the pneumatic axis. Only the distal sleeve is concentric.
tube(
    "Ejection rail",
    [
        (0.72, 0.12, 1.98),
        (0.72, 0.12, 1.28),
        (0.76, 0.12, 0.55),
        (0.74, 0.12, -0.55),
        (0.66, 0.105, -1.25),
        (0.55, 0.08, -2.91),
        (0.43, 0.04, -4.39),
        (0.36, 0.025, -5.15),
    ],
    0.095,
    ejector,
    steel,
)
lathe(
    "Tip-ejector sleeve",
    [
        (-5.45, 0.316),
        (-5.45, 0.407),
        (-5.40, 0.43),
        (-5.22, 0.44),
        (-5.08, 0.417),
        (-5.06, 0.335),
        (-5.44, 0.316),
    ],
    ejector,
    steel,
)
panel("Ejector rail shoe", (0.69, 0.12, 1.59), 0.31, 1.02, 0.24, navy, ejector, 0.05)
fluid = part("11  Illustrative aspirate")
lathe(
    "Liquid meniscus inside the distal tip",
    [
        (-10.96, 0.01),
        (-10.96, 0.034),
        (-9.65, 0.095),
        (-8.66, 0.145),
        (-8.62, 0.142),
        (-8.67, 0.10),
        (-8.69, 0.008),
    ],
    fluid,
    water,
)
# Reuse the actual meshes across views, keeping the parts independently selectable.
defs = list(parts.items())
closed_root = st.group("Instrument axis | assembled", closed.collection)
closed_root.rotation_euler = (0, math.radians(57), math.radians(-3))
for label, c in defs:
    for o in c.objects:
        o.parent = closed_root
        o["assembly_group"] = label
shift = [
    (0, 0, 5.1),
    (0, 0, 7.05),
    (0.43, 0, 5.5),
    (0, 0, 2.28),
    (0, 0, 0),
    (0, 0, -4.55),
    (0, 0, -5.43),
    (0, 0, -8.25),
    (0, 0, -11.82),
    (3.10, -0.25, -6.4),
    None,
]


def place(target, sel, moves, angle):
    root = st.group("Presentation axis | " + target.name, target.collection)
    root.rotation_euler = (0, math.radians(angle), math.radians(-3))
    for i in sel:
        label, original = defs[i]
        col = st.collection(label + " | " + target.name[:2], target.collection)
        g = st.group(label, col)
        g.parent = root
        g.location = moves[i]
        for obj in original.objects:
            o = obj.copy()
            o.data = obj.data
            col.objects.link(o)
            o.parent = g
    return root


boom = place(exploded, range(10), shift, 61)
local = {3: (0, 0, 2.0), 4: (0, 0, 0), 5: (0, 0, -4.13), 6: (0, 0, -4.88), 7: (0, 0, -7.42)}
place(detail, (3, 4, 5, 6, 7), local, 65)
for s, loc, t, span, colors in [
    (exploded, (-7, -53, 15), (-2.2, 0, -1.15), 43, ("37525D", "0F1B27")),
    (closed, (-8, -36, 14), (-0.25, 0, 0), 24, ("405863", "101A29")),
    (detail, (-4, -27, 10), (-4.95, 0, -2.4), 19.8, ("344B52", "101A21")),
]:
    c = st.camera(s, loc, t, span, perspective=False)
    st.studio(s, c, t, span, colors, world_color="A3B8C4", intensity=0.54)
    s["design"] = (
        "Air displacement. Functional layout after Gilson PIPETMAN Classic P200 parts drawings; exterior, materials and dimensions are illustrative."
    )
# Mark the three small working groups directly on the exploded plate.
bpy.context.window.scene = exploded
bpy.context.view_layer.update()
legend = st.collection("Functional labels", exploded.collection)
fontzh = bpy.data.fonts.load(str(fontdir / "Songti-SC-Regular.ttf"))
fontzh.pack()
labelmat = bpy.data.materials.new("Legend in warm ivory")
labelmat.use_nodes = True
nt = labelmat.node_tree
nt.nodes.clear()
em = nt.nodes.new("ShaderNodeEmission")
em.inputs[0].default_value = (*st.linear("EFEEE0"), 1)
em.inputs[1].default_value = 1.15
nt.links.new(em.outputs[0], nt.nodes.new("ShaderNodeOutputMaterial").inputs[0])
wire = st.shader("Leader | pale teal", "9CB8BA", 0.64, emission=0.65)
cam = exploded.camera
q = cam.rotation_euler.to_quaternion()
right = q @ Vector((1, 0, 0))
up = q @ Vector((0, 1, 0))
front = q @ Vector((0, 0, 1))
target = Vector((-2.2, 0, -1.15))
span = 43
h = span * 2 / 3


def textpoint(uv):
    return target + front * 6 + right * ((uv[0] - 0.5) * span) + up * ((uv[1] - 0.5) * h)


def project(part_index, pt):
    world = boom.matrix_world @ (Vector(pt) + Vector(shift[part_index]))
    delta = world - target
    return (0.5 + delta.dot(right) / span, 0.5 + delta.dot(up) / h)


def callout(body, start, line, anchors):
    d = bpy.data.curves.new(body, "FONT")
    d.font = fontzh
    d.body = body
    d.align_x = "LEFT"
    d.align_y = "BOTTOM_BASELINE"
    d.size = 1.37
    obj = bpy.data.objects.new(body, d)
    legend.objects.link(obj)
    obj.location = textpoint(start)
    obj.rotation_euler = cam.rotation_euler
    d.materials.append(labelmat)
    pts = [textpoint(p) for p in line]
    tube("Leader " + body, pts, 0.024, legend, wire)
    for anchor in anchors:
        end = textpoint(anchor)
        tube("Pointer " + body, [pts[-1], end], 0.024, legend, wire)
        radius = 0.066
        vv = [
            end + right * (radius * np.cos(t)) + up * (radius * np.sin(t))
            for t in np.linspace(0, 2 * np.pi, 16, endpoint=False)
        ]
        st.mesh("Index " + body, vv, [tuple(range(16))], legend, labelmat, False)


callout("活塞与复位弹簧", (0.29, 0.717), [(0.293, 0.693), (0.482, 0.693)], [project(4, (0, -0.34, -0.08))])
callout(
    "密封导向件与密封圈",
    (0.075, 0.508),
    [(0.08, 0.49), (0.349, 0.49)],
    [project(5, (0, -0.26, -1.44)), project(6, (0, -0.26, -1.71))],
)
callout("退吸头杆", (0.59, 0.236), [(0.59, 0.219), (0.51, 0.219)], [project(9, (0.56, -0.03, -2.55))])
st.finish(a, [(exploded, "exploded"), (closed, "assembled"), (detail, "piston")])
