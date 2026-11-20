#!/usr/bin/env python3
"""Regression tests of the Linear Strain Element mapping in real SOFA scenes.

Both scenes are run headless with runSofa and compared with the independent console solver of
src/liegroups/Tests/unit/test_LinearStrainElementCantilever.cpp (same stiffness matrix, same kernel,
Gauss-Newton):

  * tutorials/advanced/lse/lse_cantilever.py : 4-element cantilever, F = 1 N (EI = 0.2 N.m^2,
    EA = 1000 N, nu = 0.3), tip at (0.612151294, 0, 0.715351765);
  * tutorials/advanced/lse/lse_ring.py : closed ring with 16 elements under diametral compression,
    diameter changes 1.00385 (vertical) and 1.00125 (horizontal) times the thin-ring analytic values;
  * tutorials/advanced/lse/lse_gravity.py : cantilever under gravity with lumped nodal mass; the damped implicit
    dynamics must settle on the static equilibrium (extension of the reference formulation, which is static).

The tests are skipped when runSofa cannot be found. Point them to a SOFA build with the SOFA_ROOT
environment variable (the directory containing bin/runSofa), or put runSofa on the PATH.
"""

import os
import re
import shutil
import subprocess
import unittest
from pathlib import Path

SCENES = Path(__file__).resolve().parents[2] / "tutorials" / "advanced" / "lse"


def find_run_sofa():
    root = os.environ.get("SOFA_ROOT")
    if root and (Path(root) / "bin" / "runSofa").exists():
        return str(Path(root) / "bin" / "runSofa")
    return shutil.which("runSofa")


def run_scene(name, iterations, env=None):
    run_sofa = find_run_sofa()
    if run_sofa is None:
        raise unittest.SkipTest("runSofa not found (set SOFA_ROOT to a SOFA build)")
    result = subprocess.run(
        [run_sofa, "-g", "batch", "-n", str(iterations), "-l", "SofaPython3", str(SCENES / name)],
        capture_output=True,
        text=True,
        timeout=600,
        env={**os.environ, **(env or {})},
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output[-2000:]
    assert "[ERROR]" not in output, output[-2000:]
    return output


class LinearStrainElementSceneTest(unittest.TestCase):
    def test_cantilever_tip_matches_reference_solver(self):
        # tip of the 4-element rod at F = 1 N; SOFA Newton stops at residualThreshold = 1e-12
        expected = (0.612151294, 0.0, 0.715351765)
        output = run_scene("lse_cantilever.py", 14)
        tips = re.findall(r"LSE_TIP step=(\d+) F=([\d.]+) x=(-?[\d.]+) y=(-?[\d.]+) z=(-?[\d.]+)", output)
        self.assertGreaterEqual(len(tips), 12, "the scene did not print the tip position")
        step, force, x, y, z = tips[-1]
        self.assertAlmostEqual(float(force), 1.0)
        for value, ref, axis in zip((x, y, z), expected, "xyz"):
            self.assertAlmostEqual(float(value), ref, delta=2e-6, msg=f"tip {axis} at step {step}")

    def test_closed_ring_matches_reference_solver_and_theory(self):
        output = run_scene("lse_ring.py", 3)
        reports = re.findall(r"LSE_RING step=(\d+) N=(\d+) vertical_ratio=([\d.]+) horizontal_ratio=([\d.]+)", output)
        self.assertGreaterEqual(len(reports), 3, "the scene did not print the diameter changes")
        _, n, vertical, horizontal = reports[-1]
        self.assertEqual(int(n), 16)
        # same stiffness matrix as the console solver
        self.assertAlmostEqual(float(vertical), 1.00385, delta=2e-4)
        self.assertAlmostEqual(float(horizontal), 1.00125, delta=2e-4)
        # and within 0.5 % of the thin-ring theory (Timoshenko)
        self.assertAlmostEqual(float(vertical), 1.0, delta=5e-3)
        self.assertAlmostEqual(float(horizontal), 1.0, delta=5e-3)

    def test_gravity_dynamics_settles_on_static_equilibrium(self):
        pattern = r"LSE_GRAVITY step=(\d+) x=(-?[\d.]+) y=(-?[\d.]+) z=(-?[\d.]+)"
        static = re.findall(pattern, run_scene("lse_gravity.py", 14, {"LSE_GRAVITY_MODE": "static"}))
        dynamic = re.findall(pattern, run_scene("lse_gravity.py", 3000, {"LSE_GRAVITY_MODE": "dynamic"}))
        self.assertEqual(len(static), 14)
        self.assertEqual(len(dynamic), 3000)
        # measured: static (0.94807457, -0.29650403), dynamic after 3000 steps (0.94807452, -0.29650414)
        for axis, s, d in zip("xyz", static[-1][1:], dynamic[-1][1:]):
            self.assertAlmostEqual(float(d), float(s), delta=1e-5, msg=f"tip {axis}")
        self.assertAlmostEqual(float(static[-1][3]), -0.29650403, delta=1e-5)


if __name__ == "__main__":
    unittest.main()
