# -*- coding: utf-8 -*-
"""
Comparison between Frames2StrainCosseratMapping and Strain2RigidCosseratMappin
"""

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

from basic_functions import (_add_cosserat_frame, _add_cosserat_state,
                                    _add_rigid_base, add_mini_header)

stiffness_param: float = 1e10
v_damping_param: float =  8e-1  # Damping parameter for dynamics
beam_length :float = 5
nb_section : int = 10
beam_mass: float = 1.0


def createScene(root):
    add_mini_header(root)

    # Add gravity
    root.gravity = [0, -9.81, 0]  # Add gravity!

    # Configure time integration and solver
    solver = root.addChild("solver")

    solver.addObject(
        "EulerImplicitSolver",
        firstOrder="0",
        rayleighStiffness="0.1",
        rayleighMass="0.1",
        vdamping=v_damping_param,  # Damping parameter for dynamics
    )
    solver.addObject("CGLinearSolver", iterations=100, tolerance=1e-9, threshold=1e-9)

    #### Beam 1 (Green beam using Frames2StrainCosseratMapping) ####
 
    beam_geometry_params = BeamGeometryParameters(
        beam_length=beam_length,  # Total beam length
        nb_section=nb_section,
        nb_frames=nb_section,
    )    
        
    beam_geometry = CosseratGeometry(beam_geometry_params)
    
    ## base node
    rigid_base = _add_rigid_base(solver, node_name="rigid_base")
    
    ## frame node
    frame_node = _add_cosserat_frame(solver, beam_geometry, beam_mass=beam_mass)

    ## bending node
    strain_node = _add_cosserat_state(rigid_base, frame_node, beam_geometry)
   

    #### Beam 2 (Red beam using Strain2RigidCosseratMapping) ####
    
    #Using same beam parameters

    base_node2 = solver.addChild("base_node2")
    base_node2.addObject("MechanicalObject", template="Rigid3d", name="base_mo2", 
                        position=[0., 0., 0., 0., 0., 0., 1.], showIndices="1", showObject="1", showObjectScale=0.1)
    
    base_node2.addObject("RestShapeSpringsForceField", name="spring2", stiffness=stiffness_param, angularStiffness=stiffness_param,
                         external_points="0", mstate="@base_mo2", points="0", template="Rigid3d")
        
    custom_bending_states = [[0., 0., 0, 0, 0, 0] for _ in range(nb_section)]
    
    bending_node = solver.addChild("bending_node")
    bending_node.addObject("MechanicalObject", template="Vec6d", position=custom_bending_states, name="bending_state")
    bending_node.addObject("BeamHookeLawForceField",
                           crossSectionShape="circular", 
                           length=beam_geometry.section_lengths, 
                           radius=0.5, 
                           youngModulus=1.0e3, 
                           poissonRatio=0.4)

    frame_node2 = base_node2.addChild("frame_node")
    bending_node.addChild(frame_node2)

    frames_mo2 = frame_node2.addObject(
        "MechanicalObject",
        template="Rigid3d",
        name="FramesMO2",
        position=beam_geometry.frames,  # Use geometry data
        showIndices=1,
        showObject=1,
        showObjectScale=0.8,
    )

    frame_node2.addObject("UniformMass", totalMass=beam_mass)

    frame_node2.addObject(
        "Strain2RigidCosseratMapping",
        curv_abs_input=beam_geometry.curv_abs_sections,
        curv_abs_output=beam_geometry.curv_abs_frames,
        name="cosseratMapping2",
        input1=bending_node.bending_state.getLinkPath(),
        input2=base_node2.base_mo2.getLinkPath(),
        output=frames_mo2.getLinkPath(),
        debug=0,
        radius=0.5,
        color=[1.0, 0.0, 0.0, 0.5], #red
    )


    return root