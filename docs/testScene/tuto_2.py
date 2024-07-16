# -*- coding: utf-8 -*-

stiffness_param = 1.0e10
beam_radius = 1.0


def _add_rigid_base(p_node):
    rigid_base_node = p_node.addChild("rigid_base")
    rigid_base_node.addObject(
        "MechanicalObject",
        template="Rigid3d",
        name="cosserat_base_mo",
        position="0 0 0  0 0 0. 1",
        showObject=1,
        showObjectScale="0.1",
    )
    rigid_base_node.addObject(
        "RestShapeSpringsForceField",
        name="spring",
        stiffness=stiffness_param,
        angularStiffness=stiffness_param,
        external_points="0",
        mstate="@cosserat_base_mo",
        points="0",
        template="Rigid3d",
    )
    return rigid_base_node


def _add_cosserat_state(p_node, bending_states, list_sections_length, _radius=2.0):
    cosserat_coordinate_node = p_node.addChild("cosseratCoordinate")
    print(f" ===> bendind state : {bending_states}")
    cosserat_coordinate_node.addObject(
        "MechanicalObject",
        template="Vec3d",
        name="cosserat_state",
        position=bending_states,
    )
    testNode = cosserat_coordinate_node.addObject(
        "BeamHookeLawForceField",
        crossSectionShape="circular",
        length=list_sections_length,
        radius=_radius,
        youngModulus=1.0e4,
        poissonRatio=0.4,
    )
    print(f" the dire of node is : {dir(testNode)}")
    return cosserat_coordinate_node


def _add_cosserat_frame(
    p_node,
    _bending_node,
    framesF,
    _section_curv_abs,
    _frame_curv_abs,
    _radius,
    _beam_mass=0.0,
):
    cosserat_in_Sofa_frame_node = p_node.addChild("cosserat_in_Sofa_frame_node")

    _bending_node.addChild(cosserat_in_Sofa_frame_node)
    frames_mo = cosserat_in_Sofa_frame_node.addObject(
        "MechanicalObject",
        template="Rigid3d",
        name="FramesMO",
        position=framesF,
        showIndices=0.0,
        showObject=0,
        showObjectScale=0.8,
    )

    cosserat_in_Sofa_frame_node.addObject("UniformMass", totalMass=_beam_mass)

    cosserat_in_Sofa_frame_node.addObject(
        "DiscreteCosseratMapping",
        curv_abs_input=_section_curv_abs,
        curv_abs_output=_frame_curv_abs,
        name="cosseratMapping",
        input1=_bending_node.cosserat_state.getLinkPath(),
        input2=p_node.cosserat_base_mo.getLinkPath(),
        output=frames_mo.getLinkPath(),
        debug=0,
        radius=_radius,
    )
    return cosserat_in_Sofa_frame_node


def createScene(root_node):
    root_node.addObject('RequiredPlugin', name='Sofa.Component.LinearSolver.Direct') # Needed to use components [SparseLDLSolver]  
    root_node.addObject('RequiredPlugin', name='Sofa.Component.Mass') # Needed to use components [UniformMass]  
    root_node.addObject('RequiredPlugin', name='Sofa.Component.ODESolver.Backward') # Needed to use components [EulerImplicitSolver]  
    root_node.addObject('RequiredPlugin', name='Sofa.Component.SolidMechanics.Spring') # Needed to use components [RestShapeSpringsForceField]  
    root_node.addObject('RequiredPlugin', name='Sofa.Component.StateContainer') # Needed to use components [MechanicalObject]  
    root_node.addObject('RequiredPlugin', name='Sofa.Component.Visual') # Needed to use components [VisualStyle]  
  
    root_node.addObject(
        "VisualStyle",
        displayFlags="showBehaviorModels showCollisionModels showMechanicalMappings",
    )
    root_node.gravity = [0, 0.0, 0]
    root_node.addObject(
        "EulerImplicitSolver",
        firstOrder="0",
        rayleighStiffness="0.0",
        rayleighMass="0.0",
    )
    root_node.addObject("SparseLDLSolver", name="solver")

    # Add rigid base
    base_node = _add_rigid_base(root_node)

    # build beam geometry
    nb_sections = 6
    beam_length = 30
    length_s = beam_length / float(nb_sections)
    bending_states = []
    list_sections_length = []
    temp = 0.0  # where to start the base position
    section_curv_abs = [0.0]  # section/segment curve abscissa

    for i in range(nb_sections):
        bending_states.append([0, 0.0, 0.2])  # torsion, y_bending, z_bending
        list_sections_length.append((((i + 1) * length_s) - i * length_s))
        temp += list_sections_length[i]
        section_curv_abs.append(temp)
    bending_states[nb_sections - 1] = [0, 0.0, 0.2]

    # call add cosserat state and force field
    bending_node = _add_cosserat_state(root_node, bending_states, list_sections_length)

    # comment : ???
    nb_frames = 32
    length_f = beam_length / float(nb_frames)
    cosserat_G_frames = []
    frames_curv_abs = []
    cable_position_f = []  # need sometimes for drawing segment
    x, y, z = 0, 0, 0

    for i in range(nb_frames):
        sol = i * length_f
        cosserat_G_frames.append([sol + x, y, z, 0, 0, 0, 1])
        cable_position_f.append([sol + x, y, z])
        frames_curv_abs.append(sol + x)

    cosserat_G_frames.append([beam_length + x, y, z, 0, 0, 0, 1])
    cable_position_f.append([beam_length + x, y, z])
    frames_curv_abs.append(beam_length + x)

    _add_cosserat_frame(
        base_node,
        bending_node,
        cosserat_G_frames,
        section_curv_abs,
        frames_curv_abs,
        beam_radius,
    )

    return root_node
