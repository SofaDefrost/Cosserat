# -*- coding: utf-8 -*-
from useful.header import addHeader, addSolverNode, addVisual

stiffness_param = 1.e10
beam_radius = 1.

nb_sections = 6
nb_frames = 12
beam_length = 30


def _add_rigid_base(p_node):
    rigid_base_node = p_node.addChild('rigid_base')
    rigid_base_node.addObject('MechanicalObject', template='Rigid3d', name="cosserat_base_mo",
                              position="0 0 0  0 0 0. 1",
                              showObject=1, showObjectScale='0.1')
    rigid_base_node.addObject('RestShapeSpringsForceField', name='spring', stiffness=stiffness_param,
                              angularStiffness=stiffness_param, external_points="0", mstate="@cosserat_base_mo",
                              points="0", template="Rigid3d", activeDirections=[0,1,1,1,1,1,1])
    return rigid_base_node


def _add_cosserat_state(p_node, bending_states, list_sections_length, _radius=2.):
    cosserat_coordinate_node = p_node.addChild('cosseratCoordinate')
    print(f' ===> bendind state : {bending_states}')
    cosserat_coordinate_node.addObject('MechanicalObject', template='Vec6d', name='cosserat_state',
                                       position=bending_states)
    cosserat_coordinate_node.addObject('BeamHookeLawForceField', crossSectionShape='circular',
                                       length=list_sections_length, radius=2., youngModulus=1.e3,
                                       poissonRatio=0.4)
    return cosserat_coordinate_node


def _add_cosserat_frame(p_node, _bending_node, framesF, _section_curv_abs, _frame_curv_abs, _radius, _beam_mass=5):
    cosserat_in_Sofa_frame_node = p_node.addChild('cosserat_in_Sofa_frame_node')

    _bending_node.addChild(cosserat_in_Sofa_frame_node)
    frames_mo = cosserat_in_Sofa_frame_node.addObject('MechanicalObject', template='Rigid3d',
                                                      name="FramesMO", position=framesF, showIndices=0., showObject=0,
                                                      showObjectScale=0.8)

    cosserat_in_Sofa_frame_node.addObject('UniformMass', totalMass=_beam_mass)

    cosserat_in_Sofa_frame_node.addObject('DiscreteCosseratMapping', curv_abs_input=_section_curv_abs,
                                          curv_abs_output=_frame_curv_abs, name='cosseratMapping',
                                          input1=_bending_node.cosserat_state.getLinkPath(),
                                          input2=p_node.cosserat_base_mo.getLinkPath(),
                                          output=frames_mo.getLinkPath(), debug=1, radius=_radius)
    return cosserat_in_Sofa_frame_node


def createScene(root_node):
    root_node.addObject('VisualStyle', displayFlags='showBehaviorModels showCollisionModels showMechanicalMappings')
#    root_node.addObject('EulerImplicitSolver', firstOrder="0", rayleighStiffness="0.0", rayleighMass='0.0')
#    root_node.addObject('SparseLDLSolver', name='solver')

    addHeader(root_node)
    # root_node.gravity = [0, 0, 0]
    root_node.gravity = [0, -9.81, 0.]

    solver_node = addSolverNode(root_node, name="solver_node")

    # Add rigid base
    base_node = _add_rigid_base(solver_node)

    # build beam geometry

    length_s = beam_length/float(nb_sections)
    bending_states = []
    list_sections_length = []
    temp = 0.  # where to start the base position
    section_curv_abs = [0.]  # section/segment curve abscissa

    for i in range(nb_sections):
        bending_states.append([0, 0., 0., 0, 0., 0.])  # torsion, y-bending, z-bending, elongation, y-shear and z-shear
        list_sections_length.append((((i + 1) * length_s) - i * length_s))
        temp += list_sections_length[i]
        section_curv_abs.append(temp)
#    bending_states[nb_sections-1] = [0, 0.0, 0.3, 0, 0., 0.]
    bending_states[nb_sections-1] = [1., 0., 0., 0., 0., 0.]

    # call add cosserat state and force field
    bending_node = _add_cosserat_state(solver_node, bending_states, list_sections_length)

    # comment : ???

    length_f = beam_length/float(nb_frames)
    cosserat_G_frames = []
    frames_curv_abs = []
    cable_position_f = []  # need sometimes for drawing segment
    x, y, z = 0, 0, 0

    for i in range(nb_frames+1):
        sol = i * length_f
        cosserat_G_frames.append([sol + x, y, z, 0, 0, 0, 1])
        cable_position_f.append([sol + x, y, z])
        frames_curv_abs.append(sol + x)

    cosserat_G_frames[nb_frames] = [beam_length + x, y, z, 0, 0, 0, 1]
    cable_position_f[nb_frames] = [beam_length + x, y, z]
    frames_curv_abs[nb_frames] = beam_length + x

    cosserat_frames = _add_cosserat_frame(base_node, bending_node, cosserat_G_frames, section_curv_abs, frames_curv_abs,
                        beam_radius)

    ## Add a force component to test the stretch
    applied_force =[-1.e1, 0.,  0, 0, 0, 0]
    cosserat_frames.addObject('ConstantForceField', name='constForce', showArrowSize=1e10, indices=nb_frames, forces=applied_force)

    return root_node
