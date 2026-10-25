#!/usr/bin/env python3
"""
Regression test: BeamHookeLawForceField (Vec6d) must use the per-section Young's modulus.

A cantilever (clamped at frame 0, along X, under gravity) is built from two halves, stiff then soft
(Young's modulus ratio 50).  At the junction the bending moment is continuous, so the bending strain
must jump by about the ratio of the moduli (a little less, since the moment decreases towards the free
end: ((L - s_stiff) / (L - s_soft))^2 = 0.82 for the sections next to the junction, so about 41).

Before the fix, the Vec6d specialisation ignored ``variantSections`` / ``youngModulusList`` and every
section kept the common ``youngModulus``: the strain was continuous (ratio about 0.8).

Run it with the Python that SOFA is built against, with the SofaPython3 modules on the path:

    SOFA_ROOT=<sofa build> PYTHONPATH=<sofa build>/lib/python3/site-packages python3 test_beam_hooke_variant_sections.py

``COSSERAT_PLUGIN`` (default ``Cosserat``) can be the full path of a libCosserat to test.
"""
import os
import unittest

try:
    import Sofa
    import Sofa.Core
    import Sofa.Simulation
    SOFA_AVAILABLE = True
except ImportError as error:                                    # the test is skipped without SofaPython3
    print(f"Warning: could not import SOFA modules: {error}")
    SOFA_AVAILABLE = False

SECTIONS = 20
SECTION_LENGTH = 8.0                    # mm
E_STIFF, E_SOFT = 4.0e5, 8.0e3          # kPa (units: mm, kg, s)
POISSON = 0.4
STEPS = 320                             # 1.6 s at dt = 5 ms: settled under gravity

PLUGINS = ("Sofa.Component.Constraint.Projective", "Sofa.Component.LinearSolver.Direct",
           "Sofa.Component.Mapping.NonLinear", "Sofa.Component.Mass",
           "Sofa.Component.IntegrationScheme.Backward", "Sofa.Component.StateContainer")


def settled_bending_strains(**hooke_options):
    """Bending strain (about z) of every section of the cantilever once it has settled."""
    root = Sofa.Core.Node("root")
    root.addObject("RequiredPlugin", name="cosserat", pluginName=os.environ.get("COSSERAT_PLUGIN", "Cosserat"))
    for name in PLUGINS:
        root.addObject("RequiredPlugin", name="p_" + name.replace(".", "_"), pluginName=name)
    root.addObject("DefaultAnimationLoop")
    root.gravity = [0.0, -9810.0, 0.0]
    root.dt = 0.005

    length = SECTIONS * SECTION_LENGTH
    curv_abs = [length * i / SECTIONS for i in range(SECTIONS + 1)]
    solver = root.addChild("solver")
    solver.addObject("EulerImplicitIntegrationScheme", name="euler", firstOrder=False,
                     rayleighMass=0.5, rayleighStiffness=0.02)
    solver.addObject("SparseLDLSolver", name="linear_solver", template="CompressedRowSparseMatrixd")

    frames = solver.addChild("frames")
    frames_mo = frames.addObject("MechanicalObject", template="Rigid3d", name="FramesMO",
                                 position=[[s, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0] for s in curv_abs])
    frames.addObject("FixedProjectiveConstraint", name="clamp", indices=[0])
    frames.addObject("UniformMass", totalMass=0.02)

    bending = frames.addChild("bending")
    strain_mo = bending.addObject("MechanicalObject", template="Vec6d", name="strain",
                                  position=[[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]] * SECTIONS)
    options = dict(crossSectionShape="circular", length=[SECTION_LENGTH] * SECTIONS, radius=4.5,
                   innerRadius=3.0, youngModulus=E_STIFF, poissonRatio=POISSON)
    options.update(hooke_options)
    bending.addObject("BeamHookeLawForceField", **options)
    bending.addObject("Frames2StrainCosseratMapping", name="mapping", curv_abs_input=curv_abs,
                      curv_abs_output=curv_abs, input=frames_mo.getLinkPath(),
                      output=strain_mo.getLinkPath(), radius=4.5, debug=0)

    Sofa.Simulation.init(root)
    for _ in range(STEPS):
        Sofa.Simulation.animate(root, root.dt.value)
    return [float(s[2]) for s in strain_mo.position.value]


@unittest.skipUnless(SOFA_AVAILABLE, "SofaPython3 is not available")
class TestBeamHookeVariantSections(unittest.TestCase):
    half = SECTIONS // 2

    def junction_ratio(self, strains):
        return strains[self.half] / strains[self.half - 1]

    def test_variant_sections_use_their_own_young_modulus(self):
        strains = settled_bending_strains(
            variantSections=True, youngModulusList=[E_STIFF] * self.half + [E_SOFT] * self.half,
            poissonRatioList=[POISSON] * SECTIONS)
        ratio = self.junction_ratio(strains)
        expected = (E_STIFF / E_SOFT) * ((SECTIONS - self.half - 0.5) / (SECTIONS - self.half + 0.5)) ** 2
        self.assertAlmostEqual(ratio, expected, delta=0.1 * expected)

    def test_uniform_beam_is_unchanged(self):
        strains = settled_bending_strains()
        self.assertLess(self.junction_ratio(strains), 1.0)       # same stiffness: the moment only decreases

    def test_lists_of_the_wrong_size_fall_back_to_the_common_stiffness(self):
        strains = settled_bending_strains(variantSections=True, youngModulusList=[E_STIFF] * 3,
                                          poissonRatioList=[POISSON] * 3)
        self.assertLess(self.junction_ratio(strains), 1.0)       # an error is printed, nothing explodes
        self.assertEqual(strains, settled_bending_strains())     # and it is the same beam as the uniform one


if __name__ == "__main__":
    unittest.main()
