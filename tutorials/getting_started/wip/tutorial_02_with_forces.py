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
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'python'))
# The helpers of the former tutorial_01_basic_beam now live in using_lie_algebra/introduction_and_setup.py
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'using_lie_algebra'))

from cosserat import BeamGeometryParameters, CosseratGeometry

from introduction_and_setup import (_add_cosserat_frame, _add_cosserat_state,
                                    _add_rigid_base)

stiffness_param: float = 1.e10
beam_radius: float = 1.
vdamping_param: float = 1.e-1


def createScene(root_node):
    """Create a Cosserat beam scene with forces and dynamics."""
    # Load required plugins
    root_node.addObject('RequiredPlugin', name='Sofa.Component.LinearSolver.Direct')
    root_node.addObject('RequiredPlugin', name='Sofa.Component.Mass')
    root_node.addObject('RequiredPlugin', name='Sofa.Component.IntegrationScheme.Backward')
    root_node.addObject('RequiredPlugin', name='Sofa.Component.SolidMechanics.Spring')
    root_node.addObject('RequiredPlugin', name='Sofa.Component.SolidMechanics.FEM.Elastic')
    root_node.addObject('RequiredPlugin', name='Sofa.Component.StateContainer')
    root_node.addObject('RequiredPlugin', name='Sofa.Component.MechanicalLoad')  # ConstantForceField
    root_node.addObject('RequiredPlugin', name='Sofa.Component.Mapping.NonLinear')
    root_node.addObject('RequiredPlugin', name='Sofa.Component.Visual')
    root_node.addObject("RequiredPlugin", name='Cosserat')

    # Configure scene with dynamics
    root_node.addObject(
        "VisualStyle",
        displayFlags="showBehaviorModels showCollisionModels showMechanicalMappings",
    )
    root_node.gravity = [0, -9.81, 0]  # Add gravity!
    # vdamping no longer exists (SOFA >= 26.06): its damping is added to rayleighMass of EulerImplicitIntegrationScheme.
    solver_node = root_node.addObject(
        "EulerImplicitIntegrationScheme",
        firstOrder="0",
        rayleighStiffness="0.0",
        rayleighMass=vdamping_param,  # Damping parameter for dynamics
    )
    # solver_node.setAttribute("vdamping", 0.02)  # Set time step for dynamics
    root_node.addObject("SparseLDLSolver", name="solver")

    # === NEW APPROACH: Use CosseratGeometry with more sections for smoother dynamics ===
    beam_geometry_params = BeamGeometryParameters(
        beam_length=30.0,   # Same beam length
        nb_section=6,       # 6 sections for good physics resolution
        nb_frames=32        # 32 frames for smooth visualization
    )

    # Create geometry object
    beam_geometry = CosseratGeometry(beam_geometry_params)

    print(f"🚀 Created dynamic beam with:")
    print(f"   - Length: {beam_geometry.get_beam_length()}")
    print(f"   - Sections: {beam_geometry.get_number_of_sections()}")
    print(f"   - Frames: {beam_geometry.get_number_of_frames()}")
    print(f"   - Mass will be distributed across frames")

    # Create rigid base
    base_node = _add_rigid_base(root_node)

    # Create bending states with a curve (last section has more bending)
    custom_bending_states = []
    for i in range(beam_geometry.get_number_of_sections()):
        if i == beam_geometry.get_number_of_sections() - 1:
            # Last section has more bending
            custom_bending_states.append([0, 0.0, 0.2])
        else:
            # Other sections have slight bending
            custom_bending_states.append([0, 0.0, 0.1])

    # Create cosserat state using geometry
    bending_node = _add_cosserat_state(root_node, beam_geometry, node_name="cosserat_state", custom_bending_states=custom_bending_states)

    # Create cosserat frame with mass (important for dynamics!)
    frame_node = _add_cosserat_frame(base_node, bending_node, beam_geometry, beam_mass=5.0)

    # === ADD FORCES ===
    # Add a force at the tip of the beam
    tip_frame_index = beam_geometry.get_number_of_frames()  # Last frame
    applied_force = [-10.0, 0.0, 0.0, 0, 0, 0]  # Force in -X direction

    frame_node.addObject(
        'ConstantForceField',
        name='tipForce',
        indices=[tip_frame_index],
        forces=[applied_force],
        showArrowSize=1e-1,
    )

    print(f"💪 Applied force {applied_force[:3]} at frame {tip_frame_index}")

    return root_node
