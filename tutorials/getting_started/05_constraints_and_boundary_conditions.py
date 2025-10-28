# -*- coding: utf-8 -*-
"""
Tutorial 05: Constraints and Boundary Conditions
================================================

This tutorial demonstrates how to apply constraints to a Cosserat beam.
We will create a simple "bridge" by fixing both ends of the beam,
showing how it deforms under gravity.

Key concepts:
- Applying constraints to specific degrees of freedom.
- Using `FixedConstraint` to lock a point in space.
"""

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from python.cosserat import BeamGeometryParameters, CosseratGeometry

from introduction_and_setup import (_add_cosserat_frame, _add_cosserat_state,
                                    _add_rigid_base, add_mini_header)

v_damping_param: float = 0.4  # Damping parameter for dynamics

def createScene(root_node):
    """Create a Cosserat beam scene with constraints."""
    # Configure scene with time integration
    add_mini_header(root_node)

    # Add gravity
    root_node.gravity = [0, -9.81, 0]

    # Configure time integration and solver
    solver_node = root_node.addChild("solver")
    solver_node.addObject(
        "EulerImplicitSolver",
        rayleighStiffness="0.0",
        rayleighMass="0.0",
        vdamping=v_damping_param,
    )
    solver_node.addObject("SparseLDLSolver", name="solver")

    # Define beam geometry
    beam_geometry_params = BeamGeometryParameters(
        beam_length=40.0,
        nb_section=40,
        nb_frames=40,
    )
    beam_geometry = CosseratGeometry(beam_geometry_params)

    # Create the beam nodes
    base_node = _add_rigid_base(solver_node)
    bending_node = _add_cosserat_state(solver_node, beam_geometry)
    frame_node = _add_cosserat_frame(
        base_node, bending_node, beam_geometry, beam_mass=10.0
    )

    # --- CONSTRAINT ---
    # Fix the tip of the beam to create a bridge
    tip_frame_index = beam_geometry.get_number_of_frames() -1

    # Add a FixedConstraint to the last frame of the beam.
    # This will lock its position and orientation.
    frame_node.addObject(
        "FixedConstraint",
        name="bridgeConstraint",
        indices=[tip_frame_index], # Index of the frame to constrain
    )

    print("✨ Created a beam bridge by fixing both ends.")
    print(f"   - Tip frame index constrained: {tip_frame_index}")

    return root_node

