# -*- coding: utf-8 -*-
"""Static cantilever with Linear Strain Elements (LSE).

Large-deflection cantilever of Xun & Tamadazte, "Geometrically explicit Cosserat-rod modeling with
piecewise linear strain for complex rod systems", Computers & Structures 325 (2026) 108181, §6.1:
L = 1 m, bending stiffness EI = 0.2 N.m^2, axial stiffness EA = 1000 N, a dead force F along z at the tip.

The rod has N elements. Each element has two Rigid3 nodal poses and one strain slope; the mean strain is
rebuilt by LinearStrainElementMapping (xi_bar = A^-1 Log(g_a^-1 g_b), A = h I - h^3/12 ad_beta) and the elastic
energy is two standard Hooke laws (see cosserat.linear_strain). With 4 elements the tip is within 0.1 % of a
fine mesh (the constant strain element, beta = 0, is off by 0.5 to 1.7 %).

The load is ramped over several time steps so that every Newton solve starts close to the equilibrium.
Run:  runSofa -g batch -n 12 lse_cantilever.py     (prints the tip position at every step)
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
NB_ELEMENTS = 4
FORCE = 1.0  # N, along +z
LOAD_STEPS = 10

# circular section giving EI and EA: I = pi r^4 / 4, A = pi r^2  =>  r^2 = 4 EI / EA
RADIUS = math.sqrt(4.0 * BENDING_STIFFNESS / AXIAL_STIFFNESS)
YOUNG_MODULUS = AXIAL_STIFFNESS / (math.pi * RADIUS**2)

PLUGINS = [
    "Sofa.Component.AnimationLoop",
    "Sofa.Component.Constraint.Projective",
    "Sofa.Component.IntegrationScheme.Backward",
    "Sofa.Component.LinearSolver.Direct",
    "Sofa.Component.Mapping.Linear",
    "Sofa.Component.MechanicalLoad",
    "Sofa.Component.SolidMechanics.Spring",
    "Sofa.Component.StateContainer",
    "Cosserat",
]


class LoadRamp(Sofa.Core.Controller):
    """Ramps the tip force and prints the tip position once the solver has run."""

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.force = kwargs["force_object"]
        self.nodes = kwargs["nodes_mo"]
        self.f_max = kwargs["f_max"]
        self.steps = kwargs["steps"]
        self.step = 0
        self.tip_index = kwargs["tip_index"]

    def onAnimateBeginEvent(self, event):
        self.step += 1
        scale = min(self.step, self.steps) / self.steps
        with self.force.forces.writeable() as f:
            f[0][2] = self.f_max * scale

    def onAnimateEndEvent(self, event):
        tip = self.nodes.position.value[self.tip_index]
        print(
            f"LSE_TIP step={self.step} F={self.f_max * min(self.step, self.steps) / self.steps:.4f} "
            f"x={tip[0]:.8f} y={tip[1]:.8f} z={tip[2]:.8f}"
        )


def createScene(root, n_elements=NB_ELEMENTS, force=FORCE, load_steps=LOAD_STEPS):
    root.addObject("RequiredPlugin", pluginName=PLUGINS)
    root.gravity = [0.0, 0.0, 0.0]
    root.dt = 1.0
    root.addObject("DefaultAnimationLoop")

    solver_node = root.addChild("solverNode")
    solver_node.addObject(
        "StaticEquilibriumIntegrationScheme",
        maxNbIterationsNewton=100,
        maxNbIterationsLineSearch=10,
        residualThreshold=1e-12,
        alwaysAdvanceNewton=True,
    )
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

    rod.node.addObject("FixedProjectiveConstraint", name="clamp", indices=[0])
    load = rod.node.addObject(
        "ConstantForceField", name="tipLoad", template="Rigid3d", indices=[n_elements], forces=[[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
    )

    root.addObject(
        LoadRamp(
            name="loadRamp",
            force_object=load,
            nodes_mo=rod.nodes_mo,
            f_max=force,
            steps=load_steps,
            tip_index=n_elements,
        )
    )
    return root
