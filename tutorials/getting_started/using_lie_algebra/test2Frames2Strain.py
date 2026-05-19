# -*- coding: utf-8 -*-
"""
Test de la fonction applyJ de Frames2StrainCosseratMapping
"""

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

from introduction_and_setup import (_add_cosserat_frame, _add_cosserat_state,
                                    _add_rigid_base, add_mini_header)

v_damping_param: float =  8e-1  # Damping parameter for dynamics
beam_length :float = 15.0
nb_section : int = 3
section_length: float = beam_length/float(nb_section)
stiffness_param: float = 1.e10

def createScene(root_node):
    """Create a Cosserat beam scene with forces and dynamics."""
    # Configure scene with time integration
    add_mini_header(root_node)
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.LinearSolver.Direct') # Needed to use components [SparseLDLSolver]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.ODESolver.Backward') # Needed to use components [EulerImplicitSolver]  

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

    # === NEW APPROACH: Use CosseratGeometry with more sections for smoother dynamics ===
    beam_geometry_params = BeamGeometryParameters(
        beam_length=beam_length,
        nb_section=nb_section,
        nb_frames=nb_section,
    )

    # Create geometry object
    beam_geometry = CosseratGeometry(beam_geometry_params)

    # Create rigid base
    base_node = _add_rigid_base(solver_node)
    # Create bending states with a curve (last section has more bending)
    # We start with a slightly bent beam to better visualize the effect of gravity.
    # custom_bending_states = []
    # for i in range(beam_geometry.get_number_of_sections()):
    #     custom_bending_states.append([0, 0.0, 0.1])

    # # Create cosserat state using geometry
    # bending_node = _add_cosserat_state(solver_node, beam_geometry, node_name="cosserat_states",
    #                                    custom_bending_states=custom_bending_states)

    # # Create cosserat frame with mass (important for dynamics!)
    # # The mass is distributed uniformly along the beam. Without mass, the beam
    # # would not be affected by gravity.
    # frame_node = _add_cosserat_frame(
    #     base_node, bending_node, beam_geometry, beam_mass=5.0
    # )

    ### Tige 2: (same base_node)
    section_curv_abs = [i*section_length for i in range(nb_section+1)]
    frame_curv_abs = section_curv_abs

    frame_positions = [[0., 0., 0., 0., 0., 0., 1.],
                       [4.79426, 1.22417, 0., 0., 0., 0.247404, 0.968912],
                       [8.41471, 4.59698, 0., 0., 0., 0.479426, 0.877583],
                       [9.97495, 9.29263, 0., 0., 0., 0.681639, 0.731689]]
    

    frame_node2 = solver_node.addChild("frame_node")
    base_node.addChild(frame_node2)

    frames_mo2 = frame_node2.addObject("MechanicalObject", template="Rigid3d",
                                    name="FramesMO", position=frame_positions, 
                                    showIndices="1", showObject="1", showObjectScale=0.8)

    frame_node2.addObject("RestShapeSpringsForceField", stiffness=stiffness_param, angularStiffness=stiffness_param,
                         external_points="0", points="0", template="Rigid3d")
    
    frame_node2.addObject("UniformMass", totalMass=5.0)
    ##Strain (sortie du mapping)
    #*******strain state**********
    strain_node = frame_node2.addChild("strain_node")
    # frame_node2.addChild(strain_node)
    
    strain_node.addObject("MechanicalObject", template="Vec3d", name="strain_state")
    
    strain_node.addObject("BeamHookeLawForceField",
                           crossSectionShape="circular", 
                           length=[section_length]*nb_section, 
                           radius=0.5, 
                           youngModulus=1.0e4, 
                           poissonRatio=0.38)
    
    strain_node.addObject("Frames2StrainCosseratMapping", curv_abs_input=section_curv_abs, 
                        curv_abs_output=frame_curv_abs, name="cosseratMapping2", 
                        input1=frames_mo2.getLinkPath(), 
                        input2=base_node.cosserat_base_mo.getLinkPath(), 
                        output=strain_node.strain_state.getLinkPath(), 
                        debug=0,
                        radius=0.5,
                        color=[0., 1., 0., 0.5])   
 

    return root_node
