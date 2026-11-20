# -*- coding: utf-8 -*-
"""Linear Strain Element cantilever under gravity: static equilibrium or damped implicit dynamics.

Extension of the reference formulation (Xun & Tamadazte, Comput. Struct. 325 (2026) 108181), which is
static: mass is lumped on the nodal Rigid3 poses (UniformMass, so the mass matrix is native on the parent
degrees of freedom) and the rod is integrated either

  * statically (default for LSE_GRAVITY_MODE=static): Newton equilibrium with gravity ramped over several steps, or
  * dynamically (LSE_GRAVITY_MODE=dynamic): implicit Euler with Rayleigh damping, which settles on the same
    equilibrium.

The two modes must agree at rest; tests/unit/test_lse_scenes.py checks it. The tip position is printed at every step.
Run:  LSE_GRAVITY_MODE=dynamic runSofa -g batch -n 1500 lse_gravity.py
"""

import math
import os
import sys

import Sofa

# Make the maintained cosserat package win over any installed legacy duplicate.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))

from cosserat.linear_strain import add_linear_strain_rod  # noqa: E402

LENGTH = 1.0  # m
BENDING_STIFFNESS = 0.2  # EI, N.m^2
AXIAL_STIFFNESS = 1000.0  # EA, N
POISSON_RATIO = 0.3
TOTAL_MASS = 0.05  # kg
GRAVITY = -9.81  # m/s^2 along z
NB_ELEMENTS = 8
STATIC_RAMP_STEPS = 10
DT_DYNAMIC = 0.01

RADIUS = math.sqrt(4.0 * BENDING_STIFFNESS / AXIAL_STIFFNESS)  # circular section giving EI and EA
YOUNG_MODULUS = AXIAL_STIFFNESS / (math.pi * RADIUS**2)

PLUGINS = [
    "Sofa.Component.AnimationLoop",
    "Sofa.Component.Constraint.Projective",
    "Sofa.Component.IntegrationScheme.Backward",
    "Sofa.Component.LinearSolver.Direct",
    "Sofa.Component.Mapping.Linear",
    "Sofa.Component.Mass",
    "Sofa.Component.SolidMechanics.Spring",
    "Sofa.Component.StateContainer",
    "Cosserat",
]


class TipReport(Sofa.Core.Controller):
    """Ramps the gravity in static mode and prints the tip position after every step."""

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.root = kwargs["root"]
        self.nodes = kwargs["nodes_mo"]
        self.tip = kwargs["tip_index"]
        self.ramp_steps = kwargs["ramp_steps"]  # 0: no ramp (dynamic)
        self.step = 0

    def onAnimateBeginEvent(self, event):
        self.step += 1
        if self.ramp_steps > 0:
            scale = min(self.step, self.ramp_steps) / self.ramp_steps
            self.root.gravity = [0.0, 0.0, GRAVITY * scale]

    def onAnimateEndEvent(self, event):
        tip = self.nodes.position.value[self.tip]
        print(f"LSE_GRAVITY step={self.step} x={tip[0]:.8f} y={tip[1]:.8f} z={tip[2]:.8f}")


def createScene(root, n_elements=NB_ELEMENTS):
    mode = os.environ.get("LSE_GRAVITY_MODE", "dynamic")
    if mode not in ("static", "dynamic"):
        raise ValueError("LSE_GRAVITY_MODE must be 'static' or 'dynamic'")

    root.addObject("RequiredPlugin", pluginName=PLUGINS)
    root.addObject("DefaultAnimationLoop")

    solver_node = root.addChild("solverNode")
    if mode == "static":
        root.gravity = [0.0, 0.0, 0.0]
        root.dt = 1.0
        solver_node.addObject(
            "StaticEquilibriumIntegrationScheme",
            maxNbIterationsNewton=100,
            maxNbIterationsLineSearch=10,
            residualThreshold=1e-12,
            alwaysAdvanceNewton=True,
        )
        ramp_steps = STATIC_RAMP_STEPS
    else:
        root.gravity = [0.0, 0.0, GRAVITY]
        root.dt = DT_DYNAMIC
        solver_node.addObject("EulerImplicitSolver", rayleighStiffness=0.1, rayleighMass=0.1)
        ramp_steps = 0
    solver_node.addObject("SparseLDLSolver", name="linearSolver", template="CompressedRowSparseMatrixd")

    h = LENGTH / n_elements
    node_poses = [[i * h, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0] for i in range(n_elements + 1)]
    elements = [[i, i + 1] for i in range(n_elements)]
    rod = add_linear_strain_rod(
        solver_node,
        node_poses,
        elements,
        young_modulus=YOUNG_MODULUS,
        poisson_ratio=POISSON_RATIO,
        radius=RADIUS,
    )
    # the mass is native on the parent poses (lumped, equal on every node)
    rod.node.addObject("UniformMass", name="mass", totalMass=TOTAL_MASS)
    rod.node.addObject("FixedProjectiveConstraint", name="clamp", indices=[0])

    root.addObject(
        TipReport(name="report", root=root, nodes_mo=rod.nodes_mo, tip_index=n_elements, ramp_steps=ramp_steps)
    )
    return root
