
stiffness_param: float = 1.e10
beam_radius: float = 1.0
beam_mass: float = 0.
beam_length: float = 30
nb_section: int = 3
section_length: float = beam_length/float(nb_section)


import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))


from cosserat import BeamGeometryParameters, CosseratGeometry

from introduction_and_setup import (_add_cosserat_frame, _add_cosserat_state,
                                    _add_rigid_base)


def createScene(root):
    root.gravity = [0., 0., 0.]
    
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
                       [8.41471, 4.59698, 0., 0., 0., 0.479426, 0.877583],
                       [9.09297, 14.1615, 0., 0., 0., 0.841471, 0.540302],
                       [1.4112, 19.8999, 0., 0., 0., 0.997495, 0.0707372]]

    frame_node = solver.addChild("frame_node")

    frames_mo = frame_node.addObject("MechanicalObject", template="Rigid3d",
                                    name="Frames_MO", position=frame_positions, 
                                    showIndices="1", showObject="1", showObjectScale=0.8)

    frame_node.addObject("UniformMass", totalMass=beam_mass)

    ##Strain (sortie du mapping)
    #*******strain state**********
    strain_node = rigid_base.addChild("strain_node")
    frame_node.addChild(strain_node)
    
    strain_node.addObject("MechanicalObject", template="Vec3d", name="cosserat_state", 
                                position=[[0., 0., 0.]*nb_section])
    
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
    # Define beam geometry parameters
    beam_geometry_params = BeamGeometryParameters(
        beam_length=beam_length,  # Total beam length
        nb_section=nb_section,  # Number of sections for physics
        nb_frames=nb_section,  # Number of frames for visualization
    )

    # Create geometry object - this automatically calculates all the geometry!
    beam_geometry = CosseratGeometry(beam_geometry_params)

    base_node = _add_rigid_base(root, node_name="rigid_base")
    custom_bending_states = [
        [0.0, 0., 0.1] for _ in range(beam_geometry_params.nb_section)
    ]

    # Create cosserat state using the geometry object
    bending_node = _add_cosserat_state(root, beam_geometry, node_name="c_state_1",
                                       custom_bending_states=custom_bending_states)
    _add_cosserat_frame(base_node, bending_node, beam_geometry, node_name="frame_1", beam_mass=0.0)
    



    return root


