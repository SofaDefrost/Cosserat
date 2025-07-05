# -*- coding: utf-8 -*-
"""
Tutorial 02: Cosserat Beam with Forces
=====================================

This tutorial builds on Tutorial 01 by adding:
- Gravity forces
- Applied forces at the beam tip
- Mass distribution
- Solver configuration for dynamic simulation

Key improvements over manual approach:
- CosseratGeometry handles all geometry calculations
- Easy to modify beam parameters
- Clean, readable code structure
"""

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

from _00_introduction_and_setup import (_add_cosserat_frame, _add_cosserat_state,
                                    _add_rigid_base, add_mini_header)

v_damping_param: float = 8.e-1  # Damping parameter for dynamics

def createScene(root_node):
    """Create a Cosserat beam scene with forces and dynamics."""
    # Configure scene with time integration
    add_mini_header(root_node)

    # Add gravity
    root_node.gravity = [0, -9.81, 0]  # Add gravity!

    # Configure time integration and solver
    solver_node = root_node.addChild("solver_1")

    solver_node.addObject(
        "EulerImplicitSolver",
        firstOrder="0",
        rayleighStiffness="0.0",
        rayleighMass="0.0",
        vdamping=v_damping_param,  # Damping parameter for dynamics
    )
    solver_node.addObject("SparseLDLSolver", name="solver")

    # === NEW APPROACH: Use CosseratGeometry with more sections for smoother dynamics ===
    beam_geometry_params = BeamGeometryParameters(
        beam_length=30.0,  # Same beam length
        nb_section=30,  # 30 sections for good physics resolution
        nb_frames=30,  # 30 frames for smooth visualization
    )

    # Create geometry object
    beam_geometry = CosseratGeometry(beam_geometry_params)

    print(f"🚀 Created dynamic beam with:")
    print(f"   - Length: {beam_geometry.get_beam_length()}")
    print(f"   - Sections: {beam_geometry.get_number_of_sections()}")
    print(f"   - Frames: {beam_geometry.get_number_of_frames()}")
    print(f"   - Mass will be distributed across frames")

    # Create rigid base
    base_node = _add_rigid_base(solver_node)

    # Create bending states with a curve (last section has more bending)
    custom_bending_states = []
    for i in range(beam_geometry.get_number_of_sections()):
        custom_bending_states.append([0, 0.0, 0.0])

    # Create cosserat state using geometry
    bending_node = _add_cosserat_state(solver_node, beam_geometry, node_name="cosserat_states",
                                       custom_bending_states=custom_bending_states)

    # Create cosserat frame with mass (important for dynamics!)
    frame_node = _add_cosserat_frame(
        base_node, bending_node, beam_geometry, beam_mass=5.0
    )


    # # -------------------------------------------------
    # # === ADD SECOND BEAM WITH DIFFERENT PARAMETERS===
    # # -------------------------------------------------

    # Configure time integration and solver
    solver_node2 = root_node.addChild("solver_2")

    solver_node2.addObject(
        "EulerImplicitSolver",
        firstOrder="0",
        rayleighStiffness="0.0",
        rayleighMass="0.0",
        name="euler_solver2",
        vdamping=v_damping_param
    )
    solver_node2.addObject("SparseLDLSolver", name="solver2")

    # # Define second beam geometry parameters
    beam_geometry_params2 = BeamGeometryParameters(
        beam_length=30.0,  # Same beam length
        nb_section=15,  # 30 sections for good physics resolution
        nb_frames=15,  # 3 frames for smooth visualization
    )

    # # Create second geometry object
    # This beam has fewer sections (15) than the first one (30).
    # This will make the simulation faster, but less accurate.
    # It's a trade-off between performance and physical fidelity.
    beam_geometry2 = CosseratGeometry(beam_geometry_params2)
    print(f"🚀 Created second dynamic beam with:")
    print(f"   - Length: {beam_geometry2.get_beam_length()}")
    print(f"   - Sections: {beam_geometry2.get_number_of_sections()}")
    print(f"   - Frames: {beam_geometry2.get_number_of_frames()}")
    print(f"   - Mass will be distributed across frames")

    # # Create rigid base for second beam
    base_node2 = _add_rigid_base(solver_node2, node_name="rigid_base2")

    # # Create cosserat state for the second beam
    # # -------------------------------------------------
    custom_bending_states2 = []
    for i in range(beam_geometry2.get_number_of_sections()):
        custom_bending_states.append([0, 0.0, 0.0])

    bending_node2 = _add_cosserat_state(
        solver_node2, beam_geometry2, node_name="cosserat_states2",
        custom_bending_states=custom_bending_states2
    )

    # # Create cosserat frame for the second beam
    _add_cosserat_frame(
        base_node2, bending_node2, beam_geometry2, node_name="cosserat_in_Sofa_frame_node2",
        beam_mass=5.0
    )


    # # === ADD FORCES ===
    # # Add a force at the tip of the beam
    # tip_frame_index = beam_geometry.get_number_of_frames()  # Last frame
    # applied_force = [-10.0, 0.0, 0.0, 0, 0, 0]  # Force in -X direction
    #
    # frame_node.addObject(
    #     "ConstantForceField",
    #     name="tipForce",
    #     indices=[tip_frame_index],
    #     forces=[applied_force],
    #     showArrowSize=1,
    # )
    #
    # print(f"💪 Applied force {applied_force[:3]} at frame {tip_frame_index}")

    return root_node
