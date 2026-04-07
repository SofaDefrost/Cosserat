# -*- coding: utf-8 -*-
from useful.header import addHeader, addSolverNode, addVisual
from math import pi

stiffness_param = 1.e10

nb_sections = 3
nb_frames = 3


youngModulus: float = 1.205e11
poissonRatio: float = 0.499
beam_radius: float = 0.004
density: float = 7.850e3
beam_length: float = 2

totalMass = density * beam_length * beam_radius * beam_radius * pi


def _add_rigid_base(p_node, _name='rigid_base'):
    rigid_base_node = p_node.addChild(_name)
    rigid_base_node.addObject('MechanicalObject', template='Rigid3d', name="cosserat_base_mo",
                              position="0 0 0  0 0 0. 1",
                              showObject=0, showObjectScale='0.')
    rigid_base_node.addObject('RestShapeSpringsForceField', name='spring', stiffness=stiffness_param,
                              angularStiffness=stiffness_param, external_points="0", mstate="@cosserat_base_mo",
                              points="0", template="Rigid3d", activeDirections=[1,1,1,1,1,1,1])
    return rigid_base_node


def _add_cosserat_state(p_node, bending_states, list_sections_length, _radius=beam_radius,_template='Vec6d'):
    cosserat_coordinate_node = p_node.addChild('cosseratCoordinate')
    print(f' ===> bendind state : {bending_states}')
    cosserat_coordinate_node.addObject('MechanicalObject', template=_template, name='cosserat_state',
                                       position=bending_states)
    cosserat_coordinate_node.addObject('BeamHookeLawForceField', crossSectionShape='circular',
                                       length=list_sections_length, radius=_radius, youngModulus=youngModulus,
                                       poissonRatio=poissonRatio)
    return cosserat_coordinate_node


def _add_cosserat_frame(p_node, _bending_node, framesF, _section_curv_abs, _frame_curv_abs, _radius=beam_radius, _beam_mass=totalMass):
    cosserat_in_Sofa_frame_node = p_node.addChild('cosserat_in_Sofa_frame_node')

    _bending_node.addChild(cosserat_in_Sofa_frame_node)
    frames_mo = cosserat_in_Sofa_frame_node.addObject('MechanicalObject', template='Rigid3d',
                                                      name="FramesMO", position=framesF, showIndices=0., showObject=0,
                                                      showObjectScale=0.0)

    cosserat_in_Sofa_frame_node.addObject('UniformMass', totalMass=_beam_mass)

    cosserat_in_Sofa_frame_node.addObject('DiscreteCosseratMapping', 
    					  curv_abs_input=_section_curv_abs,
                                          curv_abs_output=_frame_curv_abs, 
                                          name='cosseratMapping',
                                          input1=_bending_node.cosserat_state.getLinkPath(),
                                          input2=p_node.cosserat_base_mo.getLinkPath(),
                                          output=frames_mo.getLinkPath(), printLog=True, radius=_radius)
    return cosserat_in_Sofa_frame_node


def createScene(root_node):
    root_node.addObject('VisualStyle', displayFlags='hideBehaviorModels hideCollisionModels hideMechanicalMappings')
#    root_node.addObject('EulerImplicitSolver', firstOrder="0", rayleighStiffness="0.0", rayleighMass='0.0')
#    root_node.addObject('SparseLDLSolver', name='solver')

    addHeader(root_node)
    # root_node.gravity = [0, 0, 0]
    root_node.gravity = [0, -9.81, 0.]

    solver_node_1 = addSolverNode(root_node, name="solver_node_1")
    solver_node_2 = addSolverNode(root_node, name="solver_node_2")

    # Add rigid base
    base_node = _add_rigid_base(solver_node_1)
    base_node_2 = _add_rigid_base(solver_node_2,_name='rigid_base2')

    # build beam geometry

    length_s = beam_length/float(nb_sections)
    bending_states = []
    bending_states_2 = []
    list_sections_length = []
    temp = 0.  # where to start the base position
    section_curv_abs = [0.]  # section/segment curve abscissa

    for i in range(nb_sections):
        bending_states.append([0, 0., 0., 0, 0., 0.])  # torsion, y-bending, z-bending, elongation, y-shear and z-shear
        bending_states_2.append([0, 0., 0.])  # torsion, y-bending
        list_sections_length.append((((i + 1) * length_s) - i * length_s))
        temp += list_sections_length[i]
        section_curv_abs.append(temp)
    bending_states[nb_sections-1] = [0, 0.0, 0.0, 1, 0., 0.]
    bending_states_2[nb_sections-1] = [0, 0.0, 0.0]
#    bending_states[nb_sections-1] = [1., 0., 0., 0., 0., 0.]

    # call add cosserat state and force field
    bending_node = _add_cosserat_state(solver_node_1, bending_states, list_sections_length)
    bending_node_2 = _add_cosserat_state(solver_node_2, bending_states_2, list_sections_length,_template='Vec3d')

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
    cosserat_frames_2 = _add_cosserat_frame(base_node_2, bending_node_2, cosserat_G_frames, section_curv_abs, frames_curv_abs,
                                            beam_radius)

    ## Add a force component to test the stretch
    #applied_force =[-1.e1, 0.,  0, 0, 0, 0]
    #cosserat_frames.addObject('ConstantForceField', name='constForce', showArrowSize=1e10, indices=nb_frames, forces=applied_force)

    return root_node
