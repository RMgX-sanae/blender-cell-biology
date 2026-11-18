"""Exercise sparse edits using the existing assembly fixture, including two copies of each atom."""

import copy
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from audit_coordinate_edits import compare, read_source
from test_source_assembly import fixture


class CoordinateEdits(unittest.TestCase):
    def setUp(self):
        tmp = TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = Path(tmp.name) / "source.cif"

        def fractional(block):
            block.find_values("_atom_site.Cartn_x")[0] = "0.321052931"

        digest = fixture(path, mutate=fractional)
        self.src = src = read_source(path, "1", digest)
        self.ion = ("Z", "1", "A", "99", "B", "CA", "CA")
        self.ix = src.keys.index(self.ion)
        self.added = ("X", "2", "A", "11", "A", "ALA", "N")
        self.plane = [k for k in src.keys if k[:2] == ("X", "1")]
        self.policy = {
            "source_sha256": digest,
            "assembly_id": "1",
            "moves": [{"id": self.ion, "max_displacement_A": 0.25}],
            "additions": [{"id": self.added, "element": "N"}],
            "preserve_planes": [
                {"id": "signed", "plane": self.plane, "site": self.ion, "max_height_change_A": 0.10}
            ],
        }
        model = np.vstack((src.xyz, [0.25, 1.5, 0.1]))
        model[self.ix, 2] += 0.05
        self.arrays = dict(
            atom_ids=np.array([*src.keys, self.added], dtype="U16"),
            elements=np.array([*src.elements, "N"], dtype="U2"),
            source_A=src.xyz.copy(),
            model_A=model,
        )

    def test_assembly_instances_reordering_and_precision(self):
        src = self.src
        np.testing.assert_allclose(
            src.xyz[src.keys.index(("X", "2", "A", "11", "A", "ALA", "CA"))], [0, 2, 0]
        )
        for dtype in (np.float32, np.float64):
            a = self.arrays
            # The insertion is now first; the full row order and the source-only order both change.
            shuffled = {k: v[::-1].copy() for k, v in a.items()}
            for k in ("source_A", "model_A"):
                shuffled[k] = shuffled[k].astype(dtype)
            r = compare(src, shuffled, self.policy)
            self.assertEqual(r["unchanged_atoms"], len(src.keys) - 1)
            self.assertEqual(r["retained_nonpolymer_atoms"], 2)
            self.assertEqual({tuple(e["id"]) for e in r["edits"]}, {self.ion, self.added})
            self.assertAlmostEqual(r["preserved_planes"]["signed"]["change_A"], 0.05, places=5)
        comp = copy.deepcopy(self.arrays)
        comp["atom_ids"] = np.delete(comp["atom_ids"], 1, axis=1)
        with self.assertRaisesRegex(ValueError, "assembly operation"):
            compare(src, comp, self.policy)

    def test_identity_and_inventory(self):
        for column, value in ((0, "Y"), (1, "2"), (2, "B"), (3, "98"), (4, ""), (5, "MG"), (6, "C")):
            a = copy.deepcopy(self.arrays)
            a["atom_ids"][self.ix, column] = value
            with self.subTest(column=column), self.assertRaisesRegex(ValueError, "inventory|Duplicate"):
                compare(self.src, a, self.policy)
        for k in ("label", "element", "remove", "parent", "source"):
            a, p = copy.deepcopy(self.arrays), copy.deepcopy(self.policy)
            if k == "label":
                a["atom_ids"][0] = a["atom_ids"][1]
            if k == "element":
                a["elements"][self.ix] = "MG"
            if k == "remove":
                a = {n: np.delete(v, self.ix, axis=0) for n, v in a.items()}
            if k == "parent":
                p["additions"][0]["id"] = (*self.added[:3], "50", *self.added[4:])
            if k == "source":
                p["assembly_id"] = "4"
            with (
                self.subTest(k=k),
                self.assertRaisesRegex(
                    ValueError, "Duplicate|Element|inventory|residue instance|different source"
                ),
            ):
                compare(self.src, a, p)

    def test_wrong_placements_and_rows(self):
        for problem in ("source", "nan", "bound", "wrongcopy", "ulp", "ulp_near_zero"):
            a = copy.deepcopy(self.arrays)
            if problem == "source":
                a["source_A"][[0, 1]] = a["source_A"][[1, 0]]
            if problem == "nan":
                a["model_A"][-1, 0] = np.nan
            if problem == "bound":
                a["model_A"][self.ix, 0] += 5
            if problem == "wrongcopy":
                i = self.src.keys.index((self.ion[0], "2", *self.ion[2:]))
                a["model_A"][i, 0] += 0.1
            if problem in ("ulp", "ulp_near_zero"):
                a["source_A"] = a["source_A"].astype("f4")
                v = a["model_A"] = a["model_A"].astype("f4")
                column = 0 if problem == "ulp" else 2
                v[0, column] = np.nextafter(v[0, column], np.float32(np.inf))
            with (
                self.subTest(problem=problem),
                self.assertRaisesRegex(ValueError, "source positions|nonfinite|exceeds|Unapproved"),
            ):
                compare(self.src, a, self.policy)

    def test_reflected_site_preserves_distances_to_each_anchor(self):
        p = copy.deepcopy(self.policy)
        p["moves"][0]["max_displacement_A"] = 2.5
        a = self.arrays
        anchors = self.src.xyz[[self.src.keys.index(k) for k in self.plane]]
        before = self.src.xyz[self.ix]
        mirrored = a["model_A"][self.ix] = before * [1, 1, -1]
        np.testing.assert_allclose(
            np.linalg.norm(before - anchors, axis=1), np.linalg.norm(mirrored - anchors, axis=1)
        )
        with self.assertRaisesRegex(ValueError, "Oriented plane"):
            compare(self.src, a, p)

    def test_degenerate_plane_and_invented_source_site(self):
        for problem in ("duplicate", "collinear", "added", "nan"):
            p = copy.deepcopy(self.policy)
            rule = p["preserve_planes"][0]
            if problem == "duplicate":
                rule["plane"][0] = rule["plane"][1]
            if problem == "collinear":
                rule["plane"][2] = ("Y", "1", "B", "10", "", "ALA", "CA")
            if problem == "added":
                rule["site"] = self.added
            if problem == "nan":
                rule["max_height_change_A"] = float("nan")
            with self.subTest(problem=problem), self.assertRaisesRegex(ValueError, "plane"):
                compare(self.src, self.arrays, p)


if __name__ == "__main__":
    unittest.main()
