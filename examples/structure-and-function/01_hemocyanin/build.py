"""The complete decamer and a window made by hiding four whole protomers."""

from pathlib import Path
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import studio as st
from views import bio_material, view


def main():
    a = st.options()
    st.init()
    total = st.scene("01  Whole oxygen carrier", a)
    inside = st.scene("02  Inside the asymmetric collar", a)
    info = json.loads((HERE / "inventory.json").read_text())
    archive = np.load(HERE / "prepared.npz")
    colors = {
        "shell1": "609DA1",
        "shell2": "487982",
        "shell3": "84B9B3",
        "shell4": "629195",
        "g": "EFB65F",
        "d-star": "EF7869",
    }
    mats = {n: bio_material(n, c, rough=0.44, bump=0.012, scale=6.0) for n, c in colors.items()}
    # Prefer the intact protomers furthest from the second camera as an open window.
    eye = np.array([0.36, -1.0, 0.22])
    proj = {
        p: np.mean([r["center_nm"] for r in info["parts"] if r["label"] == p], axis=0) @ eye
        for p in "ABCDEFGHIJ"
    }
    retained = set(sorted(proj, key=proj.get)[:6])
    for i, label in enumerate("ABCDEFGHIJ"):
        parent = st.collection(f"Protomer {label}   author {i + 1}a", total.collection)
        if label in retained:
            inside.collection.children.link(parent)
        for r in (r for r in info["parts"] if r["label"] == label):
            fu = r["fu"]
            material = mats[fu] if fu in ("g", "d-star") else mats["shell" + str(1 + (i % 4))]
            obj = st.tri_mesh(
                f"{label}  FU-{fu}", archive[r["key"] + "_v"], archive[r["key"] + "_f"], parent, material
            )
            for k in ("label", "author", "fu", "atoms", "residues", "label_seq", "auth_seq"):
                obj[k] = r[k]
    for s, pose in [
        (total, ((0, 0, 0), (0.37, -0.65, 0.83), 52.8, -0.14)),
        (inside, ((-0.3, 0.8, -0.5), eye, 41.5, -0.06)),
    ]:
        target, direction, span, roll = pose
        view(s, target, direction, span, ("263A48", "0C1822"), roll=roll, light=0.55)
        s["source"] = "6R83 assembly 1; 1 Blender unit = 1 nm; fixed fit, no new molecular coordinates"
        s["palette"] = "gold: FU-g; coral: FU-d*; blue-green: cylinder-wall FUs"
    inside["visible_complete_protomers"] = ",".join(sorted(retained))
    st.finish(a, [(total, "decamer"), (inside, "inner-collar")])


if __name__ == "__main__":
    main()
