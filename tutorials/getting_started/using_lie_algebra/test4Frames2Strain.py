# -*- coding: utf-8 -*-
"""
Test de la fonction applyJT de Frames2StrainCosseratMapping
"""

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

from introduction_and_setup import (_add_cosserat_frame_v2, _add_cosserat_state,
                                    _add_rigid_base, add_mini_header)

v_damping_param: float =  8e-1  # Damping parameter for dynamics
beam_length :float = 5
nb_section : int = 5
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
    root_node.addObject('RequiredPlugin', name='Sofa.Component.Constraint.Projective') # Needed to use components [FixedProjectiveConstraint]
    
    # Add gravity
    root_node.gravity = [0, -9.81, 0]  # Add gravity!

    root_node.dt = 0.01
    # Configure time integration and solver
    solver_node = root_node.addChild("solver")

    solver_node.addObject(
        "EulerImplicitSolver",
        firstOrder="0",
        rayleighStiffness="0.1",
        rayleighMass="0.1",
        vdamping=v_damping_param,  # Damping parameter for dynamics
    )
    solver_node.addObject("CGLinearSolver", iterations=1000, tolerance=1e-10, threshold=1e-10)

    # 1. rigid base (Pas utiliser du tout/ seulement considérer comme second entrée du Frames2StrainCosseratMapping)
    base_node = solver_node.addChild("rigid_base")
    base_node.addObject("MechanicalObject", template="Rigid3d", name="cosserat_base_mo", 
                        position=[0., 0., 0., 0., 0., 0., 1.], showIndices="1", showObject="1", showObjectScale=0.1)
    
    base_node.addObject("RestShapeSpringsForceField", name="spring", stiffness=stiffness_param, angularStiffness=stiffness_param,
                         external_points="0", mstate="@cosserat_base_mo", points="0", template="Rigid3d")
    
    # 2. Frame node
    frame_positions = [[i*section_length, 0., 0., 0., 0., 0., 1.] for i in range(nb_section+1)] 
      
    frame_node = solver_node.addChild("frame_node")
    frames_mo = frame_node.addObject("MechanicalObject", template="Rigid3d",
                                    name="FramesMO", position=frame_positions, 
                                    showIndices="1", showObject="1", showObjectScale=0.8)
    frame_node.addObject("FixedProjectiveConstraint", indices="0")
    frame_node.addObject("UniformMass", totalMass=1.0)
    # frame_node.addObject("ConstantForceField", indices="5", forces="0 0 0 0 0 1")
    # 3. Strain node
    strain_node = frame_node.addChild("strain_node")
    strain_node.addObject("MechanicalObject", template="Vec6d", name="strain_state", position=[[0., 0., 0., 1., 0, 0] for _ in range(nb_section)])
        
    strain_node.addObject("BeamHookeLawForceField",
                           crossSectionShape="circular", 
                           length=[section_length]*nb_section, 
                           radius=0.5, 
                           youngModulus=1.0e2, 
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
 


 ########## Tige 2

    # base_node2 = solver_node.addChild("base_node2")
    # base_node2.addObject("MechanicalObject", template="Rigid3d", name="cosserat_base_mo2", 
    #                     position=[0., 0., 0., 0., 0., 0., 1.], showIndices="1", showObject="1", showObjectScale=0.1)
    
    # base_node2.addObject("RestShapeSpringsForceField", name="spring2", stiffness=stiffness_param, angularStiffness=stiffness_param,
    #                      external_points="0", mstate="@cosserat_base_mo2", points="0", template="Rigid3d")
        
    # custom_bending_states = [[0., 0., 0., 1, 0, 0.] for _ in range(nb_section)]
    
    # bending_node = solver_node.addChild("bending_node")
    # bending_node.addObject("MechanicalObject", template="Vec6d", position=custom_bending_states, name="cosserat_state")
    # bending_node.addObject("BeamHookeLawForceField",
    #                        crossSectionShape="circular", 
    #                        length=[section_length]*nb_section, 
    #                        radius=0.5, 
    #                        youngModulus=1.0e2, 
    #                        poissonRatio=0.38)
    # frame_node2 = base_node2.addChild("frame_node")
    # bending_node.addChild(frame_node2)

    # frames_mo2 = frame_node2.addObject(
    #     "MechanicalObject",
    #     template="Rigid3d",
    #     name="FramesMO2",
    #     position=[[i*section_length, 0., 0, 0, 0, 0, 1] for i in range(nb_section+1)],  # Use geometry data
    #     showIndices=1,
    #     showObject=1,
    #     showObjectScale=0.8,
    # )

    # section_curv_abs = [i*section_length for i in range(nb_section+1)]
    # frame_curv_abs = section_curv_abs

    # frame_node2.addObject("UniformMass", totalMass=1.0)
    # # frame_node2.addObject("ConstantForceField", indices="5", forces="0 0 0 0 0 1")

    # frame_node2.addObject(
    #     "Strain2RigidCosseratMapping",
    #     curv_abs_input=section_curv_abs,
    #     curv_abs_output=frame_curv_abs,
    #     name="cosseratMapping",
    #     input1=bending_node.cosserat_state.getLinkPath(),
    #     input2=base_node2.cosserat_base_mo2.getLinkPath(),
    #     output=frames_mo2.getLinkPath(),
    #     debug=0,
    #     radius=0.5,
    #     color=[1.0, 0.0, 0.0, 0.5], #red
    # )


    return root_node