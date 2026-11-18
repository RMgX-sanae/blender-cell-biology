"""Backbone cartoons for intact foreground subunits and an uncut cellular context LOD."""

from pathlib import Path
import collections
import json
import os
import sys
import gemmi
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.ndimage import uniform_filter1d
from scipy.spatial import cKDTree

SOURCE_DIR = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("SHOWCASE_PROJECT", SOURCE_DIR))
out = ROOT / "prepared"
sys.path.insert(
    0,
    str(
        Path(os.environ.get("CELL_BIOLOGY_SKILL", Path.home() / ".codex/skills/blender-cell-biology"))
        / "scripts"
    ),
)
from prepare_molecular_asset import surface

m = json.loads((out / "inventory.json").read_text())
rows = {r["label_asym"]: r for r in m["parts"]}
origin = np.array(m["origin_A"])
basis = np.array(m["basis_columns"])
b = gemmi.cif.read(str(ROOT / "sources/6KGX.cif.gz")).sole_block()
model = gemmi.make_structure_from_block(b)[0]
helices = collections.defaultdict(set)
for r in b.find(
    "_struct_conf.", ["conf_type_id", "beg_label_asym_id", "beg_label_seq_id", "end_label_seq_id"]
):
    if r[0].startswith("HELX"):
        helices[r[1]].update(range(int(r[2]), int(r[3]) + 1))


def normal(v):
    return v / np.maximum(np.linalg.norm(v, axis=-1, keepdims=True), 1e-9)


def cartoon(pos, flags):
    v = []
    f = []
    ends = np.r_[0, np.where(np.linalg.norm(np.diff(pos, axis=0), axis=1) > 0.48)[0] + 1, len(pos)]
    for a, b in zip(ends, ends[1:]):
        q = pos[a:b]
        fl = flags[a:b]
        if len(q) < 3:
            continue
        s = np.arange(len(q))
        t = np.linspace(0, len(q) - 1, (len(q) - 1) * 6 + 1)
        cs = CubicSpline(s, q, axis=0)
        xyz = cs(t)
        tan = normal(cs(t, 1))
        axial = CubicSpline(s, uniform_filter1d(q, size=5, axis=0, mode="nearest"), axis=0)(t, 1)
        wide = axial - (axial * tan).sum(1)[:, None] * tan
        h = np.interp(t, s, fl)
        # Parallel transport keeps a band continuous when the projected helix axis
        # becomes ambiguous in a turn. Flip the target, never the ring indexing.
        framed = np.empty_like(wide)
        for i in range(len(wide)):
            axis = tan[i]
            if i:
                carry = framed[i - 1] - axis * np.dot(axis, framed[i - 1])
            else:
                carry = np.cross(axis, (0, 0, 1) if abs(axis[2]) < 0.85 else (1, 0, 0))
            if np.linalg.norm(carry) < 1e-7:
                carry = np.cross(axis, (0, 0, 1) if abs(axis[2]) < 0.85 else (1, 0, 0))
            carry /= np.linalg.norm(carry)
            target = wide[i]
            twist = 0.0
            if h[i] > 0.01 and np.linalg.norm(target) > 0.001:
                target = target / np.linalg.norm(target)
                if np.dot(target, carry) < 0:
                    target = -target
                twist = (
                    np.clip(
                        np.arctan2(np.dot(axis, np.cross(carry, target)), np.dot(carry, target)),
                        -0.34,
                        0.34,
                    )
                    * h[i]
                )
            framed[i] = carry * np.cos(twist) + np.cross(axis, carry) * np.sin(twist)
        wide = framed
        side = normal(np.cross(tan, wide))
        rw = 0.046 + 0.077 * h
        rn = 0.046 - 0.021 * h
        n = 10
        ang = np.arange(n) * 2 * np.pi / n
        ring = (
            xyz[:, None, :]
            + wide[:, None, :] * rw[:, None, None] * np.cos(ang)[None, :, None]
            + side[:, None, :] * rn[:, None, None] * np.sin(ang)[None, :, None]
        )
        k = sum(len(x) for x in v)
        v.append(ring.reshape(-1, 3))
        i = np.arange((len(xyz) - 1) * n)
        j = i // n * n + (i + 1) % n
        f.append(np.vstack((np.stack((i, j, j + n), axis=1), np.stack((i, j + n, i + n), axis=1))) + k)
        near = k + len(xyz) * n
        far = near + 1
        v.append(np.stack([xyz[0], xyz[-1]]))
        end = k + (len(xyz) - 1) * n
        f.extend(
            [
                np.array([(near, k + (u + 1) % n, k + u) for u in range(n)]),
                np.array([(far, end + u, end + (u + 1) % n) for u in range(n)]),
            ]
        )
    return np.vstack(v).astype("f"), np.vstack(f).astype("i")


styles = {}
host = {}
ca_all = []
id_all = []
ids = []
verts = []
faces = []
records = []
vo = fo = 0
for c in model:
    for span in c.subchains():
        label = span[0].subchain
        row = rows[label]
        if row["entity"] > 25:
            continue
        x, y, z = row["center_nm"]
        role = row["role"]
        if role == "PE":
            rep = (
                "surface"
                if y >= -0.6 or (x > 10 and y >= -6.7)
                else ("ribbon" if y >= -7.2 and x < -5 and -2 < z < 11 else "hidden")
            )
        elif role == "PC":
            rep = "ribbon" if y < -0.6 and x < 0 else ("hidden" if y < -5.3 else "surface")
        elif role == "linker":
            rep = "surface" if y >= -2.0 or row["entity"] in (21, 22) else "hidden"
        else:
            rep = "surface"
        styles[label] = rep
        ars = []
        fl = []
        for res in span:
            at = res.find_atom("CA", "*")
            if at:
                ars.append((np.array(list(at.pos)) - origin) @ basis * 0.1)
                fl.append(int(res.label_seq in helices[label]))
        if not ars:
            raise ValueError(label)
        ars = np.array(ars)
        ca_all.append(ars)
        id_all.extend([label] * len(ars))
        if rep == "ribbon" or label in ("EO", "FO"):
            v, f = cartoon(ars, np.array(fl))
            records.append([label, vo, len(v), fo, len(f)])
            vo += len(v)
            fo += len(f)
            verts.append(v)
            faces.append(f)
links = collections.defaultdict(list)
for r in b.find("_struct_conn.", ["conn_type_id", "ptnr1_label_asym_id", "ptnr2_label_asym_id"]):
    if not r[0].startswith("covale"):
        continue
    aa, bb = r[1], r[2]
    if (aa in styles) != (bb in styles):
        lig, prot = (bb, aa) if aa in styles else (aa, bb)
        links[lig].append(prot)
ca = np.vstack(ca_all)
tree = cKDTree(ca)
for label, row in rows.items():
    if label in styles:
        continue
    host[label] = links[label][0] if links[label] else id_all[tree.query(row["center_nm"])[1]]
env_v, env_f = surface(ca * 10, 4.4, 4.9, 0.8)
tv = env_v[env_f] - env_v.mean(0)
if np.einsum("ij,ij->i", tv[:, 0], np.cross(tv[:, 1], tv[:, 2])).sum() < 0:
    env_f = env_f[:, [0, 2, 1]]
np.savez_compressed(
    out / "views.npz", v=np.vstack(verts), f=np.vstack(faces), env_v=env_v * 0.1, env_f=env_f
)
(out / "views.json").write_text(
    json.dumps(
        dict(
            styles=styles,
            host=host,
            ribbons=records,
            counts=collections.Counter(styles.values()),
            display_host="Deposited covalent partner when present, otherwise nearest C-alpha for display grouping only. No new bond inferred.",
            from_connections=len(links),
        )
    )
)
print("display", collections.Counter(styles.values()), "bilins with source links", len(links))

# Two specific adjacent PE subunits, with the highest-contact beta partner for EO.
sp = {u[0].subchain: u for c in model for u in c.subchains()}


def coord_atoms(span):
    return np.array(
        [
            list(a.pos)
            for res in span
            for a in res
            if not a.element.is_hydrogen and a.occ > 0 and a.altloc in ("\x00", "A")
        ]
    )


alpha = coord_atoms(sp["EO"])
near = cKDTree(alpha)
scores = {
    r["label_asym"]: int(near.count_neighbors(cKDTree(coord_atoms(sp[r["label_asym"]])), 4.2))
    for r in rows.values()
    if r["entity"] == 2 and np.linalg.norm(np.array(r["center_nm"]) - rows["EO"]["center_nm"]) < 8
}
assert max(scores, key=scores.get) == "FO" and scores["FO"] == 219
selected = ("EO", "FO")
ligands = [l for l, p in host.items() if p in selected]
ccd = collections.defaultdict(list)
for r in b.find("_chem_comp_bond.", ["comp_id", "atom_id_1", "atom_id_2"]):
    ccd[r[0]].append([r[1], r[2]])
local = {}
for l in ligands:
    residue = sp[l][0]
    atoms = {
        a.name: dict(pos=((np.array(list(a.pos)) - origin) @ basis * 0.1).tolist(), el=a.element.name)
        for a in residue
        if not a.element.is_hydrogen and a.occ > 0 and a.altloc in ("\x00", "A")
    }
    bonds = [list(pair) for pair in ccd[residue.name] if all(k in atoms for k in pair)]
    assert len(bonds) >= len(atoms) - 1
    local[l] = dict(ccd=residue.name, protein=host[l], atoms=atoms, bonds=bonds)
(out / "pair.json").write_text(
    json.dumps(
        dict(
            chains=list(selected),
            contact_pairs_A4_2=scores["FO"],
            method="Highest number of heavy atom pairs within 4.2 A between EO and neighboring deposited PE beta chains. Residue numbers, positions, chromophores and bond names from 6KGX.",
            pigments=local,
        )
    )
)
