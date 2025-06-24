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

from examples.advanced.tuto_4 import force_null

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

from tutorial_00_basic_beam import (_add_cosserat_frame, _add_cosserat_state,
                                    _add_rigid_base, add_mini_header)
from force_controller import ForceController

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


    # === ADD FORCES ===
    # Add a force at the tip of the beam
    # this constance force is used only in the case we are doing force_type 1 or 2
    force_null = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # No force initially

    const_force_node = frame_node.addObject('ConstantForceField', name='constForce', showArrowSize=1.e-8,
                                                 indices=beam_geometry.nb_frames, forces=force_null)

    # The effector is used only when force_type is 3
    # create a rigid body to control the end effector of the beam
    tip_controller = root_node.addChild('tip_controller')
    controller_state = tip_controller.addObject('MechanicalObject', template='Rigid3d', name="controlEndEffector",
                                                showObjectScale=0.3, position=[beam_geometry.beam_length, 0, 0, 0, 0, 0, 1],
                                                showObject=True)

    return root_node
