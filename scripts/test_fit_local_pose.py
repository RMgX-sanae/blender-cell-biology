"""Local fits must preserve assembly identity, handedness conventions and intact output on failure."""
import copy
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from fit_local_pose import fit, fit_archive, values
from reaction_pose_geometry import measure_nm


class LocalPose(unittest.TestCase):
    def setUp(self):
        self.sites = {n: ["A", "1", "a", "2", "", "CYS", n] for n in ("N", "C", "CB", "CA")}
        copy_ca = [*self.sites["CA"]]
        copy_ca[1] = "2"
        # The last site is a second assembly copy with the same author atom name.
        coords = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [.25, .25, -.5], [2, 3, 9]], np.float32)
        self.arrays = dict(atom_ids=np.array([*self.sites.values(), copy_ca]), elements=np.array(list("NCCCC")),
                           source_A=coords.copy(), model_A=coords)
        self.plan = dict(version=1, basis="technical reverse-hand control", sites=self.sites, mobile=["CA"],
                         position_sigma_A=10, restraints={
                             "signed_heights": [dict(atoms=list(self.sites), target=.5, sigma=.005)],
                             "distances": [dict(atoms=[n, "CA"], target="start", sigma=.01) for n in ("N", "C", "CB")]})

    def test_fit_by_physical_identity_after_reordering(self):
        data = {k: a[[4, 3, 2, 1, 0]].copy() for k, a in self.arrays.items()}
        output, report = fit(data, self.plan)
        np.testing.assert_array_equal(output[[0, 2, 3, 4]], data["model_A"][[0, 2, 3, 4]])
        np.testing.assert_allclose(output[1], [.25, .25, .5], atol=1e-5)
        self.assertEqual(output.dtype, np.float32)
        self.assertEqual(report["mobile_sites"], 1)
        self.assertEqual(data["model_A"][1, 2], -.5)

    def test_orientation_and_short_rotation_across_180(self):
        def coords(deg):
            r = np.radians(deg)
            return np.array([[1., 0, 0], [0, 0, 0], [0, 0, 1], [np.cos(r), np.sin(r), 1]])
        q = coords(179.)
        for k, n, old in (("angles", 3, "angles"), ("signed_heights", 4, "signed_heights")):
            self.assertAlmostEqual(values(k, q[None, :n])[0], measure_nm(old, q[:n] / 10))
        self.assertAlmostEqual(values("dihedrals", q[None])[0], 179.)
        arr = dict(self.arrays, model_A=np.vstack([q, [3., 4., 5.]]))
        plan = dict(self.plan, restraints={"dihedrals": [dict(atoms=list(self.sites), target=-179., sigma=.03)]})
        end, _ = fit(arr, plan)
        self.assertAlmostEqual(values("dihedrals", end[None, :4])[0], -179., places=4)
        self.assertLess(np.linalg.norm(end[3] - q[3]), .04)

    def test_distance_floor_only_pushes_when_violated(self):
        rule = dict(atoms=["N", "CA"], target=.5, sigma=.001)
        plan = dict(self.plan, restraints={"min_distances": [rule]})
        end, _ = fit(self.arrays, plan)
        np.testing.assert_array_equal(end, self.arrays["model_A"])
        rule["target"] = 1.2
        end, _ = fit(self.arrays, plan)
        self.assertAlmostEqual(np.linalg.norm(end[3] - end[0]), 1.2, places=5)

    def test_invalid_plans_and_degenerate_geometry(self):
        for change in (dict(mobile=["CA", "CA"]), dict(mobile=["absent"]), dict(position_sigma_A=0),
                       dict(restraints={"dihedral": []}), dict(sites={**self.sites, "copy": self.sites["CA"]})):
            with self.subTest(change=change), self.assertRaises(ValueError):
                fit(self.arrays, dict(self.plan, **change))
        for field, value in (("target", float("nan")), ("sigma", -2), ("atoms", ["N", "N"])):
            plan = copy.deepcopy(self.plan)
            plan["restraints"]["distances"][0][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                fit(self.arrays, plan)
        points = np.array([[[0., 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 1]]])
        for kind in ("signed_heights", "dihedrals"):
            with self.subTest(kind=kind), self.assertRaisesRegex(ValueError, "collinear"):
                values(kind, points)

    def test_output_is_bound_to_start_and_published_only_after_convergence(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            start, plan, output = (root / name for name in ("start.npz", "plan.json", "end.npz"))
            np.savez_compressed(start, **self.arrays)
            spec = dict(self.plan, start_sha256="wrong")
            plan.write_text(json.dumps(spec))
            with self.assertRaisesRegex(ValueError, "differs"):
                fit_archive(start, plan, output)
            spec.update(start_sha256=sha256(start.read_bytes()).hexdigest(), max_nfev=1)
            plan.write_text(json.dumps(spec))
            with self.assertRaisesRegex(ValueError, "did not converge"):
                fit_archive(start, plan, output)
            self.assertFalse(output.exists())
            del spec["max_nfev"]
            plan.write_text(json.dumps(spec))
            receipt = fit_archive(start, plan, output)
            with np.load(output, allow_pickle=False) as saved:
                for field in ("source_A", "atom_ids", "elements"):
                    np.testing.assert_array_equal(saved[field], self.arrays[field])
            with self.assertRaises(FileExistsError):
                fit_archive(start, plan, output)
            self.assertEqual(receipt["candidate_sha256"], sha256(output.read_bytes()).hexdigest())


if __name__ == "__main__":
    unittest.main()
