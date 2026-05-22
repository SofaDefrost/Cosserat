# -*- coding: utf-8 -*-
"""
Test de la fonction applyJ de Frames2StrainCosseratMapping
"""

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

from introduction_and_setup import (_add_cosserat_frame_v2, _add_cosserat_state,
                                    _add_rigid_base, add_mini_header)

v_damping_param: float =  8e-1  # Damping parameter for dynamics
beam_length :float = 15.0
nb_section : int = 2
nb_frame : int = 5
dist_frame: float = beam_length/float(nb_frame)
section_length: float = beam_length/float(nb_section)
stiffness_param: float = 1.e10

def createScene(root_node):
    """Create a Cosserat beam scene with forces and dynamics."""
    # Configure scene with time integration
    add_mini_header(root_node)
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.LinearSolver.Direct') # Needed to use components [SparseLDLSolver]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.ODESolver.Backward') # Needed to use components [EulerImplicitSolver]  
    root_node.addObject('RequiredPlugin', name='Sofa.Component.LinearSolver.Iterative') # Needed to use components [CGLinearSolver]  
    root_node.addObject('RequiredPlugin', name='Sofa.Component.MechanicalLoad') # Needed to use components [ConstantForceField]  
   
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
    solver_node.addObject("CGLinearSolver", iterations=25, tolerance=1e-7, threshold=1e-7)

    # # 1. rigid base
    base_node = solver_node.addChild("rigid_base")
    base_node.addObject("MechanicalObject", template="Rigid3d", name="cosserat_base_mo", 
                        position=[0., 0., 0., 0., 0., 0., 1.], showIndices="1", showObject="1", showObjectScale=0.1)
    base_node.addObject("RestShapeSpringsForceField", name="spring", stiffness=stiffness_param, angularStiffness=stiffness_param,
                         external_points="0", mstate="@cosserat_base_mo", points="0", template="Rigid3d")
    
    # 2. Frame node
    frame_positions = [[i*dist_frame, 0., 0., 0., 0., 0., 1.] for i in range(nb_frame+1)]    

    frame_node = base_node.addChild("frame_node")
    # base_node.addChild(frame_node)
    frames_mo = frame_node.addObject("MechanicalObject", template="Rigid3d",
                                    name="FramesMO", position=frame_positions, 
                                    showIndices="1", showObject="1", showObjectScale=0.8)
    frame_node.addObject("UniformMass", totalMass=0.0)

    # # 3. Strain node
    strain_node = base_node.addChild("strain_node")
    frame_node.addChild(strain_node)
    strain_node.addObject("MechanicalObject", template="Vec3d", name="strain_state", position=[[0., 0., 0.1] for _ in range(nb_section)])
        
    strain_node.addObject("BeamHookeLawForceField",
                           crossSectionShape="circular", 
                           length=[section_length]*nb_section, 
                           radius=0.5, 
                           youngModulus=1.0e4, 
                           poissonRatio=0.38)
    
    section_curv_abs = [i*section_length for i in range(nb_section+1)]
    frame_curv_abs = section_curv_abs
    strain_node.addObject("Frames2StrainCosseratMapping", curv_abs_input=section_curv_abs, 
                        curv_abs_output=frame_curv_abs, name="cosseratMapping", 
                        input1=frames_mo.getLinkPath(), 
                        input2=base_node.cosserat_base_mo.getLinkPath(), 
                        output=strain_node.strain_state.getLinkPath(), 
                        debug=0,
                        radius=0.5,
                        color=[0., 1., 0., 0.5])   


    return root_node