# -*- coding: utf-8 -*-
"""
Tutorial 02: Cosserat Beam with Forces
=====================================

This tutorial builds on Tutorial 01 by adding:
- Gravity forces
- Applied forces at the beam tip
- Mass distribution
- Solver configuration for dynamic simulation

"""

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

from basic_functions import (_add_cosserat_frame, _add_cosserat_state,
                                    _add_rigid_base, add_mini_header)

v_damping_param: float =  8e-1  # Damping parameter for dynamics
beam_mass: float = 5
nb_section: int = 10
beam_length: int = 4
stiffness_param: float = 1e10

radius=0.5
youngModulus=1e3
poissonRatio=0.4

def createScene(root):
    """Create a Cosserat beam scene with forces and dynamics."""
    # Configure scene with time integration
    add_mini_header(root)

    # Add gravity
    root.gravity = [0, -9.81, 0]  # Add gravity!

    # Configure time integration and solver

    ## solver node
    solver = root.addChild("solver_node")
    solver.addObject("EulerImplicitSolver", firstOrder="0", rayleighMass=0.1, rayleighStiffness=0.1, vdamping=v_damping_param)
    solver.addObject("CGLinearSolver", iterations=30, tolerance=1e-10, threshold=1e-10)

    beam_geometry_params = BeamGeometryParameters(
        beam_length=beam_length,
        nb_section=nb_section,
        nb_frames=nb_section,
    )    
        
    beam_geometry = CosseratGeometry(beam_geometry_params)


    ## base node
    rigid_base = _add_rigid_base(solver, node_name="rigid_base")
    
    ## frame node
    frame_node = _add_cosserat_frame(solver, beam_geometry, beam_mass=beam_mass, constraint=True)

    ## bending node
    strain_node = _add_cosserat_state(rigid_base, frame_node, beam_geometry)


    # ## Beam 2 (red beam)

    # rigid_base2 = solver.addChild("rigid_base2")
    # rigid_base2.addObject(
    #     "MechanicalObject",
    #     template="Rigid3d",
    #     name="cosserat_base_mo2",
    #     position=[0, 0, 0, 0, 0, 0, 1],
    #     showObject=True,
    #     showObjectScale="0.1",
    # )

    # rigid_base2.addObject(
    #     "RestShapeSpringsForceField",
    #     name="spring",
    #     stiffness=stiffness_param,
    #     angularStiffness=stiffness_param,
    #     external_points="0",
    #     mstate="@cosserat_base_mo2",
    #     points="0",
    #     template="Rigid3d",
    # )

    # frame_node2 = solver.addChild("frame_node2")

    # frames_mo2 = frame_node2.addObject(
    #     "MechanicalObject",
    #     template="Rigid3d",
    #     name="FramesMO2",
    #     position=beam_geometry.frames,  # Use geometry data
    #     showIndices=1,
    #     showObject=1,
    #     showObjectScale=0.8,
    # )

    # frame_node2.addObject("UniformMass", totalMass=beam_mass)
    # # Utiliser un RestShapeSpringsForceField à la place du FixedProjectiveConstraint

    # frame_node2.addObject(
    #     "RestShapeSpringsForceField",
    #     name="spring",
    #     stiffness=stiffness_param,
    #     angularStiffness=stiffness_param,
    #     external_points="0",
    #     mstate="@FramesMO2",
    #     points="0",
    #     template="Rigid3d",
    # )

    # strain_node2 = rigid_base2.addChild("strain_node2")
    # frame_node2.addChild(strain_node2)

    # bending_states = []
    # for i in range(beam_geometry.params.nb_section):
    #     bending_states.append([0, 0, 0, 1, 0, 0]) # strain à l'équilibre 

    # strain_node2.addObject(
    #     "MechanicalObject",
    #     template="Vec6d",
    #     name="cosserat_state",
    #     position=bending_states,
    # )
    # strain_node2.addObject(
    #     "BeamHookeLawForceField",
    #     crossSectionShape="circular",
    #     length=beam_geometry.section_lengths,  # Use geometry data
    #     radius=radius,
    #     youngModulus=youngModulus,
    #     poissonRatio=poissonRatio,
    # )

    # strain_node2.addObject(
    # "Frames2StrainCosseratMapping", 
    # curv_abs_input=beam_geometry.curv_abs_sections, 
    # curv_abs_output=beam_geometry.curv_abs_frames, 
    # name="cosseratMapping", 
    # input1=frame_node2.getLinkPath(), 
    # input2=rigid_base2.cosserat_base_mo2.getLinkPath(), 
    # output=strain_node2.cosserat_state.getLinkPath(), 
    # debug=0,
    # radius=radius, 
    # color=[1., 0., 0., 0.5],
    # )


    return root

