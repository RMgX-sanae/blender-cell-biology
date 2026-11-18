"""Assembly 1 from 6KGX, with a second independently selectable pigment view."""

import json
from pathlib import Path
import os
import sys
import numpy as np

SOURCE_DIR = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("SHOWCASE_PROJECT", SOURCE_DIR))
sys.path.insert(0, str(SOURCE_DIR.parent / "src"))
import studio as st

a = st.options()
st.init()
main = st.scene("01  Light-harvesting antenna", a)
complete = st.scene("02  Intact assembly", a)
secondary = st.scene("03  Bilin array", a)
local = st.scene("04  Phycoerythrin binding pair", a)
m = json.loads((ROOT / "prepared/inventory.json").read_text())
arr = np.load(ROOT / "prepared/meshes.npz")
verts = arr["vertices"]
faces = arr["faces"]
profile = json.loads((st.SKILL / "assets/style-structure-materials.json").read_text())["surface_roles"]
view = json.loads((ROOT / "prepared/views.json").read_text())
views = np.load(ROOT / "prepared/views.npz")
protein = st.collection("Assembly 1 | 706 proteins", complete.collection)
exposed = st.collection("Core view | whole selected surfaces", main.collection)
cartoons = st.collection("Core view | whole-chain cartoons", main.collection)
selected_pig = st.collection("Core view | bilins of displayed proteins", main.collection)
pig = st.collection("Assembly 1 | 1598 bilins", complete.collection)
secondary.collection.children.link(pig)
roles = {k: st.collection(k, protein if k in ("PE", "PC", "APC", "linker") else pig) for k in m["counts"]}
source_colors = {
    1: "A83158",
    2: "CC5868",
    9: "525AA8",
    10: "657EBE",
    17: "097E88",
    18: "25AA9E",
    19: "D1A055",
    20: "62B2A6",
    26: "DE962E",
    27: "F6D681",
    28: "40CBD1",
}
materials = {}
for eid in range(1, 29):
    role = "protein" if eid <= 25 else "cofactor"
    base = dict(profile[role])
    base.update(
        bump=0.015 if eid <= 25 else 0.003, roughness=0.41 if eid <= 25 else 0.29, subsurface_weight=0.018
    )
    col = source_colors.get(eid, "BC995B")
    mat = st.material("Entity %02d  %s" % (eid, col), col, settings=base)
    if eid > 25:
        bs = mat.node_tree.nodes.get("Principled BSDF")
        bs.inputs["Emission Color"].default_value = (*st.linear(col), 1)
        bs.inputs["Emission Strength"].default_value = 0.13
    materials[eid] = mat
for r in m["parts"]:
    v = verts[r["vertex_start"] : r["vertex_start"] + r["vertex_count"]]
    f = faces[r["face_start"] : r["face_start"] + r["face_count"]]
    # Check orientation per disconnected envelope in a local frame, for correct exterior lighting.
    t = v[f] - v.mean(0)
    volume = np.einsum("ij,ij->i", t[:, 0], np.cross(t[:, 1], t[:, 2])).sum()
    if volume < 0:
        f = f[:, [0, 2, 1]]
    obj = st.tri_mesh(
        r["role"] + " | " + r["label_asym"] + " | " + r["author_chain"],
        v,
        f,
        roles[r["role"]],
        materials[r["entity"]],
    )
    for key in ("label_asym", "author_chain", "entity", "role", "atoms", "first_auth_res"):
        obj[key] = r[key]
    obj["source"] = "6KGX assembly 1"
    label = r["label_asym"]
    if view["styles"].get(label) == "surface":
        exposed.objects.link(obj)
    if label in view["host"] and view["styles"][view["host"][label]] != "hidden":
        selected_pig.objects.link(obj)
rows = {r["label_asym"]: r for r in m["parts"]}
v_all = views["v"]
f_all = views["f"]
ribbon_paint = {1: "D27182", 2: "E89975", 9: "8F83CE", 10: "86AACB"}
ribbon_mats = {
    k: st.shader("Backbone ribbon " + str(k), v, 0.37, coat=0.12) for k, v in ribbon_paint.items()
}
for label, vs, nv, fs, nf in view["ribbons"]:
    if view["styles"][label] != "ribbon":
        continue
    eid = rows[label]["entity"]
    ob = st.tri_mesh(
        "Backbone | " + label, v_all[vs : vs + nv], f_all[fs : fs + nf], cartoons, ribbon_mats[eid]
    )
    ob["source"] = "6KGX, C-alpha path and deposited helix ranges"
    ob["label_asym"] = label
# Full molecular coordinates are unaltered; membrane below the core is contextual scenery.
env = st.collection("Context | stromal face", main.collection)
complete.collection.children.link(env)
vv = []
ff = []
nx, ny = 91, 57
for j in range(ny):
    y = -80 + j * 6.5
    for i in range(nx):
        x = -140 + i * 3.2
        vv.append(
            (
                x,
                y,
                -20.95
                + 0.25 * np.sin(x * 0.059)
                + 0.32 * np.sin(y * 0.071 + x * 0.019)
                + 0.14 * np.cos(x * 0.15 + y * 0.052),
            )
        )
        if i and j:
            k = j * nx + i
            ff.append((k - nx - 1, k - nx, k, k - 1))
mem = st.material(
    "Thylakoid surface | contextual",
    "0A252C",
    settings={
        "roughness": 0.81,
        "specular": 0.09,
        "bump": 0.13,
        "texture_scale": 1.05,
        "bump_strength": 0.22,
        "subsurface_weight": 0.01,
    },
)
st.mesh("Stromal membrane surface | illustrative placement", vv, ff, env, mem)
distance_mat = st.shader("Distant antenna | optically subdued", "112D37", 0.79)
context = st.tri_mesh(
    "Neighbouring antenna | coarse source envelope", views["env_v"], views["env_f"], env, distance_mat
)
context.location = (-73, 107, -0.65)
context.rotation_euler.z = -0.70
other = context.copy()
other.data = context.data
env.objects.link(other)
other.location = (78, 142, -0.65)
other.rotation_euler.z = 1.05
for s, pal, e in [
    (main, ("2B414B", "0B1927"), 0.35),
    (complete, ("2B414B", "0B1927"), 0.35),
    (secondary, ("203A53", "091321"), 0.40),
]:
    # Slight elevation shows both rod tiers and the lower core, retaining a comparable projection.
    target = (0, 0, -0.65)
    c = st.camera(s, (9, -122, 15), target, 76, perspective=True, fstop=0.13 if s != secondary else None)
    st.studio(s, c, target, 72, pal, intensity=e)
    s["entry"] = "6KGX"
    s["protein_subunits"] = 706
    s["bilin_chromophores"] = 1598
    s["coordinates"] = "1 BU = 1 nm; source assembly intact, shared rigid frame"
# Near-field model of an actual contacting pair; all coordinates stay in the same frame.
from detail_scene import Parts
from mathutils import Matrix, Vector

pair = json.loads((ROOT / "prepared/pair.json").read_text())
cg = st.collection("Adjacent source proteins and their bound bilins", local.collection)
allv = []
for k, color in zip(pair["chains"], ["C45C82", "23868D"]):
    label, vs, nv, fs, nf = next(r for r in view["ribbons"] if r[0] == k)
    points = v_all[vs : vs + nv]
    allv.append(points)
    ob = st.tri_mesh(
        "PE backbone | " + label,
        points,
        f_all[fs : fs + nf],
        cg,
        st.shader("Local protein " + label, color, 0.36, coat=0.13),
    )
    ob["label_asym"] = label
    ob["author_chain"] = rows[label]["author_chain"]
atomcols = [
    st.shader("Bilin carbon", "EBA93B", 0.28, metal=0.10),
    st.shader("Bilin nitrogen", "658DA9", 0.32),
    st.shader("Bilin oxygen", "D66345", 0.31),
    st.shader("Bilin heteroatom", "D4B68D", 0.35),
]
for l, r in pair["pigments"].items():
    g = Parts()
    for at in r["atoms"].values():
        g.ball(at["pos"], 0.069, {"C": 0, "N": 1, "O": 2}.get(at["el"], 3), 14, 9)
    for aa, bb in r["bonds"]:
        g.tube([r["atoms"][aa]["pos"], r["atoms"][bb]["pos"]], 0.031, 0, 9)
    ob = g.object(r["ccd"] + " | " + l, atomcols, cg)
    ob["label_asym"] = l
    ob["source_host_label_asym"] = r["protein"]
coords = np.vstack(
    allv + [np.array([a["pos"] for g in pair["pigments"].values() for a in g["atoms"].values()])]
)
t = (coords.max(0) + coords.min(0)) / 2
normal = np.linalg.eigh(np.cov(coords.T))[1][:, 0]
if np.dot(normal, t) < 0:
    normal = -normal
x = np.array(rows["FO"]["center_nm"]) - np.array(rows["EO"]["center_nm"])
x -= np.dot(x, normal) * normal
x /= np.linalg.norm(x)
y = np.cross(normal, x)
theta = np.deg2rad(14)
right = x * np.cos(theta) - y * np.sin(theta)
up = x * np.sin(theta) + y * np.cos(theta)
proj = np.stack([(coords - t) @ right, (coords - t) @ up], axis=1)
span = max(np.ptp(proj[:, 0]), np.ptp(proj[:, 1]) * 1.5) * 1.22
t += right * (proj[:, 0].max() + proj[:, 0].min()) / 2 + up * (proj[:, 1].max() + proj[:, 1].min()) / 2
c = st.camera(local, t + normal * 18, t, span)
c.rotation_euler = Matrix((Vector(right), Vector(up), Vector(normal))).transposed().to_euler()
st.studio(local, c, t, span, ("283F50", "0A1829"), intensity=0.34)
local["selected_deposited_subchains"] = "EO and FO"
local["bound_bilins"] = len(pair["pigments"])
st.finish(a, [(main, "antenna"), (complete, "intact"), (secondary, "pigments"), (local, "binding_pair")])
