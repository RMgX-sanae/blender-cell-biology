"""Two registered views of the 7B9V yeast C complex and its branched RNA."""

from pathlib import Path
import json
import os
import sys

import bpy
import numpy as np
from mathutils import Vector

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import studio as st
from views import view, bio_material
from detail_scene import Parts

RNA = {
    "A": ("U2", "63BDB8"),
    "B": ("U5", "87B9D9"),
    "C": ("U6", "3AADB3"),
    "H": ("5-prime exon", "DB8191"),
    "L": ("Branched intron", "F3B65D"),
}
SUGAR = ["C1'", "C2'", "C3'", "C4'", "O4'"]
BACK = ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]
SIX = ["N1", "C2", "N3", "C4", "C5", "C6"]
FIVE = ["C4", "C5", "N7", "C8", "N9"]


def filled_ring(name, rings, vertices, faces):
    for ring in rings:
        a = np.asarray(ring)
        normal = np.cross(a[1] - a[0], a[2] - a[0])
        normal /= np.linalg.norm(normal)
        k = len(vertices)
        n = len(a)
        vertices.extend(a + normal * 0.020)
        vertices.extend(a - normal * 0.020)
        faces.extend([tuple(range(k, k + n)), tuple(reversed(range(k + n, k + 2 * n)))])
        for j in range(n):
            i = (j + 1) % n
            faces.append((k + j, k + n + j, k + n + i, k + i))


def draw_rna(label, records, xyz, col, backmat, basemat, local=False):
    frames = Parts()
    plates = []
    polygons = []
    previous = None
    for r in records:
        v = {k: xyz[i] for k, i in r["atoms"].items()}
        radius = 0.105 if not local else 0.067
        for path in [BACK, SUGAR + SUGAR[:1], ["C1'", "N9" if r["comp"] in ("A", "G") else "N1"]]:
            # Partial source residues at a modeled boundary may omit terminal atoms.
            for l, h in zip(path, path[1:]):
                if l in v and h in v:
                    frames.tube([v[l], v[h]], radius if path == BACK else radius * 0.58, sides=9)
        for atom in set(BACK + SUGAR):
            if atom in v:
                frames.ball(v[atom], radius if atom in BACK else radius * 0.62, seg=8, rings=5)
        if previous is not None and r["auth"] == previous[0] + 1 and "P" in v:
            if previous[1] is not None and np.linalg.norm(previous[1] - v["P"]) < 0.22:
                frames.tube([previous[1], v["P"]], radius, sides=9)
        previous = (r["auth"], v.get("O3'"))
        rings = [SIX] + ([FIVE] if r["comp"] in ("A", "G") else [])
        filled_ring(
            label, [[v[k] for k in ring] for ring in rings if all(k in v for k in ring)], plates, polygons
        )
    name = RNA[label][0]
    o = frames.object(name + " | measured backbone", [backmat], col)
    o["label_asym_id"] = label
    o["auth_seq_ids"] = ",".join(str(r["auth"]) for r in records)
    if len(plates):
        bases = st.mesh(name + " | planar bases", plates, polygons, col, basemat, smooth=False)
        bases["label_asym_id"] = label


def atomistic_branch(records, xyz, col, link_indices):
    # Known RNA covalent graph; the two isolated resolved intron segments have no interpolated bridge.
    selectors = {"L": {1, 2, 69, 70, 71}, "H": {-1}}
    node_colors = {
        "intron": "D6944C",
        "base": "BC925C",
        "exon": "B88D94",
        "O": "BF755F",
        "N": "698F9F",
        "P": "E5B263",
    }
    mats = {k: st.shader("Branch " + k, v, rough=0.34) for k, v in node_colors.items()}
    for label, selection in selectors.items():
        previous = None
        for r in [r for r in records if r["label"] == label and r["auth"] in selection]:
            ids = r["atoms"]
            comp = r["comp"]
            v = {k: xyz[i] for k, i in ids.items()}
            pieces = Parts()
            keys = list(mats)

            def color(k):
                if label == "H":
                    return keys.index("exon")
                return keys.index(k[0] if k[0] in ("O", "N", "P") else "intron" if "'" in k else "base")

            for k, i in ids.items():
                if i not in link_indices:
                    radius = 0.061 if k == "P" else 0.045 if label == "H" else 0.053
                    pieces.ball(v[k], radius, color(k), 14, 9)
            paths = [BACK, SUGAR + SUGAR[:1], ["P", "OP1"], ["P", "OP2"], ["C2'", "O2'"], SIX + SIX[:1]]
            if comp in ("A", "G"):
                paths += [
                    FIVE + FIVE[:1],
                    ["C1'", "N9"],
                    ["C6", "N6" if comp == "A" else "O6"],
                    ["C2", "N2"],
                ]
            else:
                paths += [["C1'", "N1"], ["C4", "N4" if comp == "C" else "O4"], ["C2", "O2"]]
            seen = set()
            for path in paths:
                for a, b in zip(path, path[1:]):
                    pair = tuple(sorted((a, b)))
                    if a not in v or b not in v or pair in seen:
                        continue
                    seen.add(pair)
                    mid = (v[a] + v[b]) * 0.5
                    pieces.tube([v[a], mid], 0.024, color(a), sides=12)
                    pieces.tube([mid, v[b]], 0.024, color(b), sides=12)
            if (
                previous
                and previous[0] + 1 == r["auth"]
                and "P" in v
                and np.linalg.norm(previous[1] - v["P"]) < 0.22
            ):
                pieces.tube([previous[1], v["P"]], 0.025, color("P"), sides=12)
            previous = (r["auth"], v["O3'"]) if "O3'" in v else None
            obj = pieces.object(f"{label} | {comp}{r['auth']}", list(mats.values()), col)
            obj["source_label"] = label
            obj["auth_seq_id"] = r["auth"]


def branch_caption(scene, cam, point):
    q = cam.rotation_euler.to_quaternion()
    right = q @ Vector((1, 0, 0))
    up = q @ Vector((0, 1, 0))
    front = q @ Vector((0, 0, 1))
    material = st.shader("Bond annotation", "FAE7AA", rough=0.38, emission=0.6)
    # The annotation occupies camera-facing space. Its leader projects onto the actual bond.
    pt = Vector(point)
    anchor = pt + front * 0.71
    end = anchor + right * 0.47 + up * 0.235
    line = Parts()
    line.tube([anchor + right * 0.05 + up * 0.035, end, end + right * 0.095], 0.007, sides=8)
    line.object("Leader to branch", [material], scene.collection)
    td = bpy.data.curves.new("Branch label", "FONT")
    td.body = "2′–5′"
    td.size = 0.139
    td.extrude = 0.0002
    font = os.environ.get("CELL_FONT_EN")
    if not font:
        raise ValueError("Set CELL_FONT_EN to Times New Roman Regular for the branch label")
    td.font = bpy.data.fonts.load(font, check_existing=True)
    td.font.pack()
    ob = bpy.data.objects.new("2′–5′ phosphodiester", td)
    scene.collection.objects.link(ob)
    ob.rotation_euler = cam.rotation_euler.copy()
    ob.location = end + right * 0.125 - up * 0.032
    td.materials.append(material)


def coordinates(s, xyz, ids):
    data = bpy.data.meshes.new("Native atom coordinates in nm")
    data.vertices.add(len(xyz))
    data.vertices.foreach_set("co", np.asarray(xyz, dtype="f4").ravel())
    att = data.attributes.new("atom_site_id", "INT", "POINT")
    att.data.foreach_set("value", ids)
    data.update()
    ob = bpy.data.objects.new("All deposited heavy atoms | select by mmCIF atom_site.id", data)
    col = st.collection("Source coordinates", s.collection)
    col.objects.link(ob)
    col.hide_render = True
    col.hide_viewport = True


def main():
    a = st.options()
    st.init()
    whole = st.scene("01  Assembled spliceosome", a)
    core = st.scene("02  Ordered RNA and docking factors", a)
    close = st.scene("03  The RNA branch", a)
    info = json.loads((HERE / "inventory.json").read_text())
    z = np.load(HERE / "geometry.npz")
    xyz = z["xyz_nm"]
    records = info["records"]
    ilink = info["substrand_link"]
    target = xyz[ilink].mean(0)
    strands = {k: bio_material(v[0], v[1], rough=0.36, sss=0.07, bump=0.004) for k, v in RNA.items()}
    bases = {k: st.shader(v[0] + " bases", v[1], rough=0.28, coat=0.23) for k, v in RNA.items()}
    body = [
        bio_material("Protein envelope " + str(i), v, rough=0.46, bump=0.008)
        for i, v in enumerate(["385970", "437C7A", "608A93", "4D6779", "547E88", "4C5F73"])
    ]
    prp8 = bio_material("Prp8 scaffold", "886D86", rough=0.42, bump=0.012)
    factors = bio_material("Branch docking factors", "B98B79", rough=0.40, bump=0.008)
    ions = st.shader("Magnesium", "9BD996", rough=0.28, metal=0.14, coat=0.25)
    ligand = st.shader("Other deposited cofactors", "97A9A1", rough=0.48)
    linkmat = st.shader("Branch phosphodiester", "FFF1A0", rough=0.25, coat=0.18, emission=0.12)
    c = st.collection("Protein scaffold by deposited chain", whole.collection)
    for i, g in enumerate(info["groups"]):
        k = g["label"]
        if "v_" + k not in z:
            continue
        mat = prp8 if k == "D" else factors if k in {"G", "I", "J"} else body[i % len(body)]
        o = st.tri_mesh(k + " | " + g["name"], z["v_" + k], z["f_" + k], c, mat)
        if k in {"G", "I", "J"}:
            core.collection.objects.link(o)
        o["label_asym_id"] = k
        o["auth_asym_id"] = g["author"]
        o["entity_id"] = g["entity"]
    coordinates(whole, xyz, z["atom_site_ids"])
    for scene, local in [(whole, False), (core, False), (close, True)]:
        rc = st.collection("RNA in its deposited position", scene.collection)
        if local:
            atomistic_branch(records, xyz, rc, ilink)
        else:
            for k in RNA:
                rs = [r for r in records if r["label"] == k]
                draw_rna(k, rs, xyz, rc, strands[k], bases[k])
        marked = Parts()
        marked.tube(xyz[ilink], 0.041 if local else 0.085, sides=18)
        for i in ilink:
            marked.ball(xyz[i], 0.061 if local else 0.12, seg=18, rings=10)
        o = marked.object("A70 O2' to G1 P | 2'-5' branch", [linkmat], rc)
        o["distance_angstrom"] = float(np.linalg.norm(xyz[ilink[0]] - xyz[ilink[1]]) * 10)
        o["endpoint_indices"] = ilink
        lc = st.collection("Deposited ions and ligands", scene.collection)
        for r in records:
            if r["label"] in RNA:
                continue
            ids = list(r["atoms"].values())
            if scene != whole and r["label"] not in ("YA", "ZA"):
                continue
            balls = Parts()
            for i in ids:
                balls.ball(
                    xyz[i], (0.085 if local else 0.124) if r["comp"] == "MG" else 0.14, seg=12, rings=8
                )
            ob = balls.object(r["label"] + " | " + r["comp"], [ions if r["comp"] == "MG" else ligand], lc)
            ob["source_label"] = r["label"]
    # Views use the same coordinates. The close view removes whole occluding proteins.
    view(whole, (0, 0, 0), (-0.74, -0.48, 0.85), 49.8, ("303543", "0C1522"), roll=0.90, light=0.56)
    view(core, (-0.1, 1.7, -1.3), (-0.40, 0.80, 0.53), 40.4, ("293944", "0A1723"), roll=-0.58, light=0.60)
    direction = Vector((-0.36, 0.39, 0.84)).normalized()
    q = (-direction).to_track_quat("-Z", "Y")
    rot = q.to_euler()
    rot.rotate_axis("Z", 0.23)
    q = rot.to_quaternion()
    right = np.asarray(q @ Vector((1, 0, 0)))
    up = np.asarray(q @ Vector((0, 1, 0)))
    selected = [
        i
        for r in records
        if (
            (r["label"] == "L" and r["auth"] in (1, 2, 69, 70, 71))
            or (r["label"] == "H" and r["auth"] == -1)
            or r["label"] in ("YA", "ZA")
        )
        for i in r["atoms"].values()
    ]
    screen = (xyz[selected] - target) @ np.stack([right, up], axis=1)
    lo = screen.min(0)
    hi = screen.max(0)
    hi = np.maximum(hi, [1.03, 0.44])
    mid = (hi + lo) / 2
    span = max(hi[0] - lo[0], 1.5 * (hi[1] - lo[1])) * 1.26
    cam = view(
        close,
        target + mid[0] * right + mid[1] * up,
        direction,
        span,
        ("293B41", "101C29"),
        roll=0.23,
        light=0.52,
    )
    branch_caption(close, cam, target)
    for scene in (whole, core, close):
        scene.view_settings.look = "AgX - Medium High Contrast"
    close["display"] = "Intron G1/U2 and 69/70/71, terminal 5-prime exon residue and Mg 201/202"
    core["display"] = (
        "All ordered U2, U5, U6, exon and intron RNA; complete deposited Yju2, Cwc25 and Isy1; two local Mg ions"
    )
    whole["source"] = close["source"] = (
        "7B9V biological assembly 1; 1 native unit = 1 nm. PDB atom-site ID maps to coordinates."
    )
    st.finish(a, [(whole, "assembly"), (core, "rna-core"), (close, "branch")])


if __name__ == "__main__":
    main()
