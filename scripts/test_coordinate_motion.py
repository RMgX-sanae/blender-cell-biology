"""Boundary identity, carried atoms and collisions hidden between accepted samples."""

from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from audit_coordinate_motion import audit


class MotionHandoff(unittest.TestCase):
    def setUp(self):
        tmp = TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.path = Path(tmp.name)
        self.candidate = self.path / "candidate.npz"
        # Two copies of the same author residue, one deposited water and a non-collinear fixed site.
        self.keys = [("X", op, "A", "17", "A", "GLY", "CA") for op in ("1", "2")]
        self.keys += [("W", "2", "A", "19", "", "HOH", "O"), ("Z", "1", "B", "8", "", "ALA", "N")]
        self.xyz = np.array([[-2, 0, 0], [2, 0, 0], [0, 0, 0], [0, 8, 0]], dtype=float)
        np.savez(
            self.candidate, model_A=self.xyz, atom_ids=np.asarray(self.keys), elements=["C", "C", "O", "N"]
        )
        digest = sha256(self.candidate.read_bytes()).hexdigest()
        self.receipt = dict(
            source_sha256="verified-raw", assembly_id="1", candidate_sha256=digest, model_heavy_atoms=4
        )
        self.plan = dict(
            source_sha256="verified-raw",
            assembly_id="1",
            candidate_sha256=digest,
            minimum_nonbonded_A=0.15,
            max_bond_change_A=0.01,
            bonds=[],
            stages=[],
        )

    def stage(self, name, rows, positions, frames=(1, 2)):
        p = self.path / (name + ".npz")
        np.savez(
            p,
            positions_A=np.asarray(positions, dtype=float),
            atom_ids=[self.keys[i] for i in rows],
            frames=frames,
        )
        self.plan["stages"].append(
            dict(id=name, file=p.name, moves=[self.keys[i] for i in reversed(rows)], max_step_A=10)
        )

    def run_plan(self):
        return audit(self.candidate, self.receipt, self.plan, self.path)

    def test_reordering_and_carry_across_an_inactive_stage(self):
        a, b, w, c = self.xyz
        self.stage("move", [1, 0], [[b, a], [b + [0, 1, 0], a + [0, -1, 0]]])
        self.stage("other", [3], [[c], [c + [0, 0, 1]]], (2, 6))
        self.stage("resume", [0], [[a + [0, -1, 0]], [a + [0, -2, 0]]], (6, 9))
        out = self.run_plan()
        self.assertEqual([s["max_start_error_A"] for s in out["stages"]], [0] * 3)
        self.assertEqual(out["stages"][-1]["frames"], [6, 9])
        p = self.path / "resume.npz"
        np.savez(p, positions_A=np.array([[a], [a + [0, -2, 0]]]), atom_ids=[self.keys[0]], frames=[6, 9])
        with self.assertRaisesRegex(ValueError, "Boundary jump"):
            self.run_plan()

    def test_stationary_obstacle_at_interior_only(self):
        a = self.xyz[0]
        self.stage("cross", [0], [[a], [[1.0, 0, 0]]])
        with self.assertRaisesRegex(ValueError, "Nonbonded collision at frame 1.66667"):
            self.run_plan()

    def test_moving_pair_crossing_and_bond_collapse(self):
        a, b = np.array([[-2.0, 2, 0], [2.0, 2, 0]])
        self.stage("rise", [0, 1], [[self.xyz[0], self.xyz[1]], [a, b]])
        self.stage("cross", [0, 1], [[a, b], [b, a]], (2, 3))
        with self.assertRaisesRegex(ValueError, "Nonbonded collision at frame 2.5"):
            self.run_plan()
        self.plan["bonds"] = [[self.keys[0], self.keys[1]]]
        with self.assertRaisesRegex(ValueError, "Bond length changed by 4"):
            self.run_plan()

    def test_contact_floor_and_global_frames(self):
        a = self.xyz[0]
        self.stage("near_water", [0], [[a], [a + [0.2, 0, 0]]])
        self.assertEqual(len(self.run_plan()["stages"]), 1)
        self.plan["component_floors_A"] = {"HOH": {"C": 1.9}}
        with self.assertRaisesRegex(ValueError, "Nonbonded collision"):
            self.run_plan()
        self.plan["component_floors_A"] = {}
        self.stage("misaligned", [0], [[a + [0.2, 0, 0]], [a]], (3, 4))
        with self.assertRaisesRegex(ValueError, "common boundary frame"):
            self.run_plan()

    def test_identity_hash_and_velocity_constraints(self):
        a = self.xyz[0]
        self.stage("motion", [0], [[a], [a + [0, 0.2, 0]]])
        stage = self.plan["stages"][0]
        stage["moves"] = [self.keys[1]]
        with self.assertRaisesRegex(ValueError, "inventory"):
            self.run_plan()
        stage["moves"] = [self.keys[0]]
        stage["max_entry_velocity_change_A_per_frame"] = 0.1
        with self.assertRaisesRegex(ValueError, "Velocity discontinuity"):
            self.run_plan()
        stage["max_entry_velocity_change_A_per_frame"] = 0.3
        stage["max_step_A"] = 0.05
        with self.assertRaisesRegex(ValueError, "sample displacement"):
            self.run_plan()
        stage["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "motion hash"):
            self.run_plan()
        self.receipt["candidate_sha256"] = "1" * 64
        with self.assertRaisesRegex(ValueError, "different audited candidate"):
            self.run_plan()


if __name__ == "__main__":
    unittest.main()
