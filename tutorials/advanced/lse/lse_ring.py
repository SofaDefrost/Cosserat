# -*- coding: utf-8 -*-
"""Closed-loop rod network with Linear Strain Elements: thin ring under diametral compression.

Elements are node pairs, so a closed loop needs no special treatment: the last element simply joins node N-1
to node 0 (a chain-based Cosserat mapping cannot do this without kinematic-closure constraints). This scene
checks the loop against the classical thin-ring solution (Timoshenko): two opposite forces P along a diameter
shorten that diameter by (pi/4 - 2/pi) P R^3 / EI = 0.1488 P R^3 / EI and lengthen the perpendicular one by
(2/pi - 1/2) P R^3 / EI = 0.1366 P R^3 / EI.

Node k sits at angle -pi/2 + 2 pi k / N and its frame has x along the counter-clockwise tangent. Node 0 (bottom) is
clamped; the force acts downwards on the top node N/2. The rest strain of each element is the mean strain of its
initial poses (curvature 1/R), computed by the mapping.

Run:  runSofa -g batch -n 6 lse_ring.py     (prints the two diameter changes relative to the analytic values)
"""

import math
import os
import sys

import Sofa

# Make the maintained cosserat package win over any installed legacy duplicate.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))

from cosserat.linear_strain import add_linear_strain_rod  # noqa: E402

RADIUS_RING = 1.0  # m
BENDING_STIFFNESS = 0.2  # EI, N.m^2
AXIAL_STIFFNESS = 1000.0  # EA, N
POISSON_RATIO = 0.3
NB_ELEMENTS = 16  # multiple of 4
FORCE = 0.001  # N, P R^2 / EI = 0.005: linear regime

RADIUS = math.sqrt(4.0 * BENDING_STIFFNESS / AXIAL_STIFFNESS)  # circular section giving EI and EA
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


class RingReport(Sofa.Core.Controller):
    """Prints the diameter changes relative to the thin-ring analytic values after every step."""

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.nodes = kwargs["nodes_mo"]
        self.n = kwargs["n_elements"]
        self.force = kwargs["force"]
        self.initial = [list(p[:3]) for p in self.nodes.position.value]
        self.step = 0
        r, ei = RADIUS_RING, BENDING_STIFFNESS
        self.analytic_vertical = (math.pi / 4.0 - 2.0 / math.pi) * self.force * r**3 / ei
        self.analytic_horizontal = (2.0 / math.pi - 0.5) * self.force * r**3 / ei

    def onAnimateEndEvent(self, event):
        self.step += 1
        x = self.nodes.position.value
        n = self.n
        vertical = self.initial[n // 2][1] - x[n // 2][1]
        horizontal = (x[n // 4][0] - x[3 * n // 4][0]) - (self.initial[n // 4][0] - self.initial[3 * n // 4][0])
        print(
            f"LSE_RING step={self.step} N={n} "
            f"vertical_ratio={vertical / self.analytic_vertical:.6f} "
            f"horizontal_ratio={horizontal / self.analytic_horizontal:.6f}"
        )


def createScene(root, n_elements=NB_ELEMENTS, force=FORCE):
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

    # frames: x along the counter-clockwise tangent, rotation about z by theta + pi/2
    node_poses = []
    for k in range(n_elements):
        theta = -0.5 * math.pi + 2.0 * math.pi * k / n_elements
        alpha = theta + 0.5 * math.pi
        node_poses.append(
            [
                RADIUS_RING * math.cos(theta),
                RADIUS_RING * math.sin(theta),
                0.0,
                0.0,
                0.0,
                math.sin(0.5 * alpha),
                math.cos(0.5 * alpha),
            ]
        )
    elements = [[k, (k + 1) % n_elements] for k in range(n_elements)]  # the last element closes the loop
    arc = 2.0 * math.pi * RADIUS_RING / n_elements

    rod = add_linear_strain_rod(
        solver_node,
        node_poses,
        elements,
        young_modulus=YOUNG_MODULUS,
        poisson_ratio=POISSON_RATIO,
        radius=RADIUS,
        lengths=[arc] * n_elements,
        name="ring",
    )
    rod.node.addObject("FixedProjectiveConstraint", name="clamp", indices=[0])
    rod.node.addObject(
        "ConstantForceField",
        name="load",
        template="Rigid3d",
        indices=[n_elements // 2],
        forces=[[0.0, -force, 0.0, 0.0, 0.0, 0.0]],
    )

    root.addObject(RingReport(name="report", nodes_mo=rod.nodes_mo, n_elements=n_elements, force=force))
    return root
