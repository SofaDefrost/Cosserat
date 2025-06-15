# -*- coding: utf-8 -*-

from tutorial_01_basic_beam import _add_rigid_base, _add_cosserat_state, _add_cosserat_frame

stiffness_param: float = 1.e10
beam_radius: float = 1.


def createScene(root_node):
    root_node.addObject('RequiredPlugin', name='Sofa.Component.LinearSolver.Direct') # Needed to use components [SparseLDLSolver]  
    root_node.addObject('RequiredPlugin', name='Sofa.Component.Mass') # Needed to use components [UniformMass]  
    root_node.addObject('RequiredPlugin', name='Sofa.Component.ODESolver.Backward') # Needed to use components [EulerImplicitSolver]  
    root_node.addObject('RequiredPlugin', name='Sofa.Component.SolidMechanics.Spring') # Needed to use components [RestShapeSpringsForceField]  
    root_node.addObject('RequiredPlugin', name='Sofa.Component.StateContainer') # Needed to use components [MechanicalObject]  
    root_node.addObject('RequiredPlugin', name='Sofa.Component.Visual') # Needed to use components [VisualStyle]
    root_node.addObject("RequiredPlugin", name='Cosserat')

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
