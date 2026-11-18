"""Whole proximal valve and a registered exploded window, with open mesh pores."""

from pathlib import Path
import os
import sys
import math
import numpy as np

SOURCE_DIR = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("SHOWCASE_PROJECT", SOURCE_DIR))
sys.path.insert(0, str(SOURCE_DIR.parent / "src"))
import studio as st

a = st.options()
st.init()
shellscene = st.scene("01  Two valves and girdle", a)
main = st.scene("02  Proximal valve", a)
macro = st.scene("03  Areola in layers", a)
pal = {
    "foramen": "B1CDBA",
    "wall": "83AAA3",
    "cribrum": "A4BFBA",
    "cribellum": "D6A45E",
    "mantle": "C9DFCE",
    "rimoportula": "E2E4C7",
}
names = {
    "foramen": "Proximal plate with foramina",
    "wall": "Areola chamber walls",
    "cribrum": "Distal perforated cribrum",
    "cribellum": "Fine distal cribella",
    "mantle": "Valve margin",
    "rimoportula": "Marginal rimoportulae",
}
mats = {}
for k, color in pal.items():
    mat = st.material(
        names[k] + " | hydrated silica",
        color,
        settings={
            "roughness": 0.31,
            "specular": 0.38,
            "subsurface_weight": 0.08,
            "subsurface_radius": [0.09, 0.07, 0.04],
            "subsurface_scale": 0.12,
            "bump": 0.0045,
            "bump_strength": 0.13,
            "texture_scale": 17.5,
        },
    )
    bs = mat.node_tree.nodes.get("Principled BSDF")
    bs.inputs["Coat Weight"].default_value = 0.21
    bs.inputs["Coat Roughness"].default_value = 0.27
    mats[k] = mat
for s, fname in [(main, "valve"), (macro, "window")]:
    geo = st.collection("Silica | " + fname, s.collection)
    arrays = np.load(ROOT / "prepared" / (fname + ".npz"))
    for k in pal:
        if k + "_v" not in arrays:
            continue
        ob = st.tri_mesh(names[k], arrays[k + "_v"], arrays[k + "_f"], geo, mats[k])
        ob["layer"] = k
        if s == macro:
            shifts = {
                "foramen": (0.95, 0, 1.23),
                "wall": (0.08, 0, 0.09),
                "cribrum": (-0.90, 0, -0.65),
                "cribellum": (-0.90, 0, -0.65),
            }
            ob.location = shifts[k]
            ob["exploded_translation_um"] = list(shifts[k])
    s["species_reference"] = "Coscinodiscus granii"
    s["evidence"] = "Literature-constrained constructed valve, schematic dimensions and optical materials."
c = st.camera(main, (18, -74, 55), (0, 0, -1.15), 56, perspective=True)
st.studio(main, c, (0, 0, -1.15), 56, ("214048", "07141C"), intensity=0.29)
for light in main.objects:
    if light.type == "LIGHT" and light.data.name.startswith("Warm edge"):
        light.data.color = st.linear("FFAB59")
        light.data.energy *= 2.4
    if light.type == "LIGHT" and light.data.name.startswith("Soft fill"):
        light.data.energy *= 0.58
c = st.camera(macro, (1.65, -11, 5.7), (0, 0, 0.40), 6.9, perspective=True)
st.studio(macro, c, (0, 0, 0.40), 6.9, ("264449", "0B1721"), intensity=0.35)
# The frustule's two thecae keep their separate complete valves. The narrower lower
# half nests in the broad upper half. A slight valve tilt makes the girdle wedge-shaped.
base = next(c for c in main.collection.children if c.name.startswith("Silica"))
rim = np.load(ROOT / "prepared/rim.npy")
ang = np.arctan2(rim[:, 1], rim[:, 0])
radius = np.linalg.norm(rim, axis=1)
theta = np.arange(480) * 2 * np.pi / 480


def perimeter(sign):
    interpolated = np.interp(sign * theta, ang, radius, period=2 * np.pi)
    harmonics = np.fft.rfft(interpolated)
    harmonics[11:] = 0
    return np.fft.irfft(harmonics, n=len(theta))


base_upper = perimeter(-1) + 0.20
base_lower = 0.96 * (perimeter(1) + 0.20)


def animate_group(name, children, at, apart, rotation=(0, 0, 0), scale=1):
    c = st.collection(name, shellscene.collection)
    g = st.group(name, c)
    for obj in children:
        ob = obj.copy()
        ob.data = obj.data
        c.objects.link(ob)
        ob.parent = g
    g.rotation_euler = rotation
    g.scale = (scale,) * 3
    g.location = at
    g.keyframe_insert(data_path="location", frame=1)
    g.location = np.array(at) + apart
    g.keyframe_insert(data_path="location", frame=60)
    return g


animate_group("Upper epivalve", list(base.objects), (0, 0, 2.8), (1.0, 0, 4.7), (math.pi, -0.072, 0))
animate_group("Lower hypovalve", list(base.objects), (0, 0, -2.4), (-0.8, 0, -4.0), scale=0.96)
bands = st.material(
    "Translucent girdle silica",
    "85AAA6",
    settings={
        "roughness": 0.32,
        "specular": 0.34,
        "subsurface_weight": 0.095,
        "subsurface_scale": 0.16,
        "bump": 0.009,
        "texture_scale": 10,
    },
)
# The band rims follow the wedged perimeter; grooves are shallow surface relief.
for idx, (r, z0, z1, shift) in enumerate(
    [(base_upper, -0.10, 2.91, 1.70), (base_lower, -2.60, 0.29, -1.44)]
):
    name = ("Epicingulum", "Hypocingulum")[idx]
    c = st.collection(name, shellscene.collection)
    g = st.group(name, c)
    profile = [
        (0, 0),
        (0.025, 0.072),
        (0.07, 0.11),
        (0.13, 0.06),
        (0.20, 0.095),
        (0.28, 0.09),
        (0.32, 0.043),
        (0.38, 0.085),
        (0.60, 0.09),
        (0.73, 0.066),
        (0.77, 0.043),
        (0.82, 0.076),
        (0.93, 0.085),
        (0.997, 0.055),
        (1, -0.07),
        (0.92, -0.09),
        (0.04, -0.10),
        (0, -0.06),
    ]
    cv = []
    cf = []
    n = len(theta)
    for t, expand in profile:
        rr = r + expand
        xx = rr * np.cos(theta)
        yy = rr * np.sin(theta)
        z = z0 + (z1 - z0) * t
        weight = (0.38 + 0.62 * t) if idx == 0 else (0.36 * t)
        zz = z + 0.072 * xx * weight + 0.13 * np.sin(2 * theta + 0.5) * weight
        cv.extend(zip(xx, yy, zz))
    for j in range(len(profile)):
        for k in range(n):
            h = (k + 1) % n
            m = ((j + 1) % len(profile)) * n
            cf.append((j * n + k, j * n + h, m + h, m + k))
    ob = st.mesh(name + " | siliceous girdle band", cv, cf, c, bands)
    ob.parent = g
    g.location = (0, 0, 0)
    g.keyframe_insert(data_path="location", frame=1)
    g.location = ((0.25 if idx == 0 else -0.25), 0, shift)
    g.keyframe_insert(data_path="location", frame=60)
shellscene.frame_start = 1
shellscene.frame_end = 60
shellscene.frame_set(60)
t = (0, 0, 0)
cam = st.camera(shellscene, (18, -82, 28), t, 57, perspective=True)
st.studio(shellscene, cam, t, 57, ("29444B", "0B1724"), intensity=0.29)
shellscene["model"] = (
    "Two illustrative complementary intact valves and two girdle bands. The wedge, openings and relative ordering are source-informed; dimensions, placement and scaling are model parameters. Key 1 assembled, key 60 exploded."
)
st.finish(a, [(shellscene, "frustule"), (main, "valve"), (macro, "areolae")])
