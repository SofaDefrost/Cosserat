# -*- coding: utf-8 -*-
"""
Test de la fonction applyJ de Frames2StrainCosseratMapping
"""

from introduction_and_setup import add_mini_header

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

    # Create rigid base
    base_node = solver_node.addChild("rigid_base")
    base_node.addObject("MechanicalObject", template="Rigid3d", name="cosserat_base_mo", 
                        position=[0., 0., 0., 0., 0., 0., 1.], showIndices="1", showObject="1", showObjectScale=0.1)
    base_node.addObject("RestShapeSpringsForceField", name="spring", stiffness=stiffness_param, angularStiffness=stiffness_param,
                         external_points="0", mstate="@cosserat_base_mo", points="0", template="Rigid3d")
    
    # frame_positions = [[i*section_length, 0., 0., 0., 0., 0., 1.] for i in range(nb_section+1)]
    frame_positions = [[0., 0., 0., 0., 0., 0., 1.],
                    [4.79426, 1.22417, 0., 0., 0., 0.247404, 0.968912],
                    [8.41471, 4.59698, 0., 0., 0., 0.479426, 0.877583],
                    [9.97495, 9.29263, 0., 0., 0., 0.681639, 0.731689]]    

    frame_node = solver_node.addChild("frame_node")
    frames_mo = frame_node.addObject("MechanicalObject", template="Rigid3d",
                                    name="FramesMO", position=frame_positions, 
                                    showIndices="1", showObject="1", showObjectScale=0.8)
    frame_node.addObject("UniformMass", totalMass=5.0)

    # frame_node.addObject("ConstantForceField", forces=[[0, 9.81, 0, 0, 0, 0]] * (nb_section+1))

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



# frame_positions = [[0., 0., 0., 0., 0., 0., 1.],
#                    [4.79426, 1.22417, 0., 0., 0., 0.247404, 0.968912],
#                    [8.41471, 4.59698, 0., 0., 0., 0.479426, 0.877583],
#                    [9.97495, 9.29263, 0., 0., 0., 0.681639, 0.731689]]