# -*- coding: utf-8 -*-
"""
Tutorial 07: Modeling Initially Curved Beams
=============================================

This tutorial demonstrates how to create a Cosserat beam with an initial
(stress-free) curvature. This is essential for modeling objects that are
naturally curved, like hooks, arches, or pre-bent robotic arms.

Key concepts:
- Defining a non-zero resting curvature.
- Programmatically generating `bendingState` to create complex shapes.
"""

import os
import sys
import numpy as np

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

from introduction_and_setup import (_add_cosserat_frame, _add_cosserat_state,
                                    _add_rigid_base, add_mini_header)

def createScene(root_node):
    """Create a scene with an initially curved Cosserat beam."""
    add_mini_header(root_node)
    root_node.gravity = [0, -9.81, 0] # Add some gravity to see it deform

    # --- Solver ---
    solver_node = root_node.addChild("solver")
    solver_node.addObject("EulerImplicitSolver", rayleighStiffness="0.0", rayleighMass="0.0", vdamping=0.1)
    solver_node.addObject("SparseLDLSolver", name="solver")

    # --- Beam Geometry ---
    beam_geometry_params = BeamGeometryParameters(
        beam_length=30.0,
        nb_section=30,
        nb_frames=30,
    )
    beam_geometry = CosseratGeometry(beam_geometry_params)

    # --- Define Initial Curvature ---
    # To create a curved beam, we define a non-zero resting curvature.
    # Here, we'll create a simple circular arch by applying a constant
    # curvature around the Z-axis for each section.
    num_sections = beam_geometry.get_number_of_sections()
    # The total curvature will define a semi-circle (pi radians)
    total_angle = np.pi
    curvature_per_section = total_angle / beam_geometry.get_beam_length()

    custom_bending_states = [[0, 0, curvature_per_section] for _ in range(num_sections)]

    # --- Create the Beam ---
    base_node = _add_rigid_base(solver_node)
    bending_node = _add_cosserat_state(solver_node, beam_geometry,
                                       custom_bending_states=custom_bending_states)
    frame_node = _add_cosserat_frame(
        base_node, bending_node, beam_geometry, beam_mass=5.0
    )

    print("✨ Created an initially curved beam (arch).")

    return root_node
