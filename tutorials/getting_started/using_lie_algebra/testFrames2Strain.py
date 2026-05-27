
stiffness_param: float = 1.e10
beam_radius: float = 1.0
beam_mass: float = 1.
beam_length: float = 5.
nb_section: int = 3
section_length: float = beam_length/float(nb_section)

# nb_frame = nb_section
# dist_frame = beam_length/nb_frame

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))


from cosserat import BeamGeometryParameters, CosseratGeometry

from introduction_and_setup import (_add_cosserat_frame, _add_cosserat_state,
                                    _add_rigid_base)


def createScene(root):
    root.gravity = [0., -9.81, 0.]
    
    #*********required plugins
    root.addObject("RequiredPlugin", name="Sofa.Component.StateContainer")
    root.addObject("RequiredPlugin", name="Sofa.Component.Visual")    
    root.addObject("RequiredPlugin", name="Sofa.Component.Mass")    
    root.addObject("RequiredPlugin", name="Sofa.Component.SolidMechanics.Spring")    
    root.addObject("RequiredPlugin", name="Cosserat")
    root.addObject("RequiredPlugin", name="Sofa.Component.ODESolver.Backward")
    root.addObject("RequiredPlugin", name="Sofa.Component.LinearSolver.Direct")
    root.addObject('RequiredPlugin', name='Sofa.Component.Constraint.Lagrangian.Correction')   
    root.addObject("VisualStyle", displayFlags="showBehaviorModels showMechanicalMappings showCollisionModels")
    root.addObject('RequiredPlugin', name='Sofa.Component.MechanicalLoad') # Needed to use components [ConstantForceField]  
    root.addObject('RequiredPlugin', name='Sofa.Component.Constraint.Projective') # Needed to use components [FixedProjectiveConstraint]  
   
    #********solver
    solver = root.addChild("solver_node")
    solver.addObject("EulerImplicitSolver", rayleighMass=0.1, rayleighStiffness=0.1)
    solver.addObject("SparseLDLSolver", template="CompressedRowSparseMatrixd")
    solver.addObject("GenericConstraintCorrection")
   

    #*****rigid base****** (OK ici)
    rigid_base = solver.addChild("rigid_base")
    rigid_base.addObject("MechanicalObject", template="Rigid3d", name="cosserat_base_mo", 
                        position=[0., 0., 0., 0., 0., 0., 1.], showIndices="1", showObject="1", showObjectScale=0.1)
    rigid_base.addObject("RestShapeSpringsForceField", name="spring", stiffness=stiffness_param, angularStiffness=stiffness_param,
                         external_points="0", mstate="@cosserat_base_mo", points="0", template="Rigid3d")
    
    
    ##inverser la logique strain->frames du Strain2Rigid ou du DiscreteCosserat
    #*******frames***********
    section_curv_abs = [i*section_length for i in range(nb_section+1)]
    frame_curv_abs = section_curv_abs

    frame_positions = [[0., 0., 0., 0., 0., 0., 1.],
                       [1.65896, 0.138568, 0, 0, 0, 0.0832369, 0.99653],
                       [3.22661, 0.681371, 0, 0, 0, 0.247404, 0.968912],
                       [4.43343, 1.80564, 0, 0, 0, 0.479426, 0.877583]]    


    frame_node = solver.addChild("frame_node")

    frames_mo = frame_node.addObject("MechanicalObject", template="Rigid3d",
                                    name="Frames_MO", position=frame_positions, 
                                    showIndices="1", showObject="1", showObjectScale=0.8)

    frame_node.addObject("UniformMass", totalMass=beam_mass)

    ##Strain (sortie du mapping)
    #*******strain state**********
    strain_node = rigid_base.addChild("strain_node")
    frame_node.addChild(strain_node)
    
    strain_node.addObject("MechanicalObject", template="Vec6d", name="cosserat_state", position=[0, 0, 0, 0, 0, 0])
    
    strain_node.addObject("BeamHookeLawForceField",
                           crossSectionShape="circular", 
                           length=[section_length]*nb_section, 
                           radius=beam_radius, 
                           youngModulus=1.0e4, 
                           poissonRatio=0.38)


    strain_node.addObject("Frames2StrainCosseratMapping", curv_abs_input=section_curv_abs, 
                        curv_abs_output=frame_curv_abs, name="cosseratMapping", 
                        input1=frames_mo.getLinkPath(), 
                        input2=rigid_base.cosserat_base_mo.getLinkPath(), 
                        output=strain_node.cosserat_state.getLinkPath(), 
                        debug=0,
                        radius=beam_radius)
    

    ################### Tige 2 ####################
    # # Define beam geometry parameters
    # base_node = _add_rigid_base(root, node_name="rigid_base")

    # custom_bending_states = [[0., 0., 0.1, 1, 0, 0.], [0., 0., 0.2, 1, 0, 0.], [0., 0., 0.3, 1, 0, 0.]]
    
    # bending_node = solver.addChild("bending_node")
    # bending_node.addObject("MechanicalObject", template="Vec6d", position=custom_bending_states, name="cosserat_state")
    # bending_node.addObject("BeamHookeLawForceField",
    #                        crossSectionShape="circular", 
    #                        length=[section_length]*nb_section, 
    #                        radius=0.5, 
    #                        youngModulus=1.0e4, 
    #                        poissonRatio=0.38)
    # frame_node = base_node.addChild("frame_node")
    # bending_node.addChild(frame_node)

    # frames_mo = frame_node.addObject(
    #     "MechanicalObject",
    #     template="Rigid3d",
    #     name="FramesMO",
    #     position=[[i*section_length, 0., 0, 0, 0, 0, 1] for i in range(nb_section+1)],  # Use geometry data
    #     showIndices=1,
    #     showObject=1,
    #     showObjectScale=0.8,
    # )

    # section_curv_abs = [i*section_length for i in range(nb_section+1)]
    # frame_curv_abs = section_curv_abs

    # frame_node.addObject("FixedProjectiveConstraint", indices="0")
    # frame_node.addObject("UniformMass", totalMass=1.0)
    # frame_node.addObject(
    #     "Strain2RigidCosseratMapping",
    #     curv_abs_input=section_curv_abs,
    #     curv_abs_output=frame_curv_abs,
    #     name="cosseratMapping",
    #     input1=bending_node.cosserat_state.getLinkPath(),
    #     input2=base_node.cosserat_base_mo.getLinkPath(),
    #     output=frames_mo.getLinkPath(),
    #     debug=0,
    #     radius=0.5,
    #     color=[1.0, 0.0, 0.0, 0.5], #red
    # )

    return root


# frame_positions = [[0., 0., 0., 0., 0., 0., 1.],
#                        [0.998334, 0.0499583, 0., 0., 0., 0.0499792, 0.99875],
#                        [1.98669, 0.199334, 0., 0., 0., 0.0998334, 0.995004],
#                        [2.9552, 0.446635, 0., 0., 0., 0.149438, 0.988771]]    
