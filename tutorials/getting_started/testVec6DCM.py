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

_beam_radius = 0.5
_beam_length = 30.
_nb_section = 32

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

from force_controller import ForceController
from introduction_and_setup import (_add_cosserat_frame, _add_cosserat_state,
                                    _add_rigid_base, add_mini_header)

v_damping_param: float = 8.e-1  # Damping parameter for dynamics

def createScene(root_node):
    """Create a Cosserat beam scene with forces and dynamics."""
    # Configure scene with time integration
    add_mini_header(root_node)
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.LinearSolver.Direct') # Needed to use components [SparseLDLSolver]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.ODESolver.Backward') # Needed to use components [EulerImplicitSolver]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.MechanicalLoad') # Needed to use components [ConstantForceField] 
    
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
    solver_node.addObject("SparseLDLSolver", template="CompressedRowSparseMatrixMat3x3d", name="solver")

    beam_geometry_params = BeamGeometryParameters(
        beam_length=10.0,
        nb_section=8,
        nb_frames=8,
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
        custom_bending_states.append([0, 0, i*0.1, 0, 0, 0])

    # Create cosserat state using geometry
    bending_node = solver_node.addChild("bending_node")
    bending_node.addObject(
        "MechanicalObject",
        template="Vec6d",
        name="cosserat_state",
        position=custom_bending_states,
    )
    bending_node.addObject(
        "BeamHookeLawForceField",
        crossSectionShape="circular",
        length=beam_geometry.section_lengths,  # Use geometry data
        radius=2.0,
        youngModulus=1.0e3,
        poissonRatio=0.4,
    )

    # Create cosserat frame with mass (important for dynamics!)
    frame_node = base_node.addChild("frame_node")
    bending_node.addChild(frame_node)
    frames_mo = frame_node.addObject(
        "MechanicalObject",
        template="Rigid3d",
        name="FramesMO",
        position=beam_geometry.frames,  # Use geometry data
        showIndices=1,
        showObject=1,
        showObjectScale=0.8,
    )

    frame_node.addObject("UniformMass", totalMass=0.0)

    frame_node.addObject(
        "Strain2FramesCosseratMapping",
        curv_abs_input=beam_geometry.curv_abs_sections,  # Use geometry data
        curv_abs_output=beam_geometry.curv_abs_frames,  # Use geometry data
        name="cosseratMapping",
        input1=bending_node.cosserat_state.getLinkPath(),
        input2=base_node.cosserat_base_mo.getLinkPath(),
        output=frames_mo.getLinkPath(),
        debug=0,
        radius=1.0,
    )


#     # === ADD FORCES ===
#     # Add a force at the tip of the beam
#     # this constance force is used only in the case we are doing force_type 1 or 2
#     force_null = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # No force initially
#    # force_active = [0.0, 30.0, 0.0, 0.0, 0.0, 0.0]  # Pull the rod upward

#     const_force_node = frame_node.addObject('ConstantForceField', name='constForce', showArrowSize=1.e-8,
#                                                  indices=[beam_geometry.get_number_of_frames()], forces=force_null)

#     # The effector is used only when force_type is 3
#     # create a rigid body to control the end effector of the beam
#     tip_controller = root_node.addChild('tip_controller')
#     controller_state = tip_controller.addObject('MechanicalObject', template='Rigid3d', name="controlEndEffector",
#                                                 showObjectScale=0.3, position=[beam_geometry.get_beam_length(), 0, 0, 0, 0, 0, 1],
#                                                 showObject=True)

#     # add the controller to animate the force
#     root_node.addObject(ForceController(
#         name="ForceController",
#         forceNode=const_force_node,     # ConstantForceField
#         frame_node=frame_node,          # node containing FramesMO
#         force_type=2,                   # Change 1, 2 or 3 to test all force type
#         tip_controller=controller_state,# a MechanicalObject used to control the beam's tip (for force_type 3)
#         geoParams=beam_geometry_params  # geometric params
#     ))

    return root_node
