# -*- coding: utf-8 -*-

import Sofa

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
    cosserat_coordinate_node.addObject(
        "MechanicalObject",
        template="Vec3d",
        name="cosserat_state",
        position=bending_states,
    )
    cosserat_coordinate_node.addObject(
        "BeamHookeLawForceField",
        crossSectionShape="circular",
        length=list_sections_length,
        radius=2.0,
        youngModulus=1.0e4,
        poissonRatio=0.4,
    )
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
        showIndices=1,
        showObject=1,
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
    root_node.addObject("RequiredPlugin", name='Sofa.Component.Mass')
    root_node.addObject("RequiredPlugin", name='Sofa.Component.SolidMechanics.Spring')
    root_node.addObject("RequiredPlugin", name='Sofa.Component.StateContainer')
    root_node.addObject("RequiredPlugin", name='Sofa.Component.Visual')
    
    root_node.addObject(
        "VisualStyle",
        displayFlags="showBehaviorModels showCollisionModels showMechanicalMappings",
    )
    root_node.gravity = [0, 0.0, 0]
    #
    base_node = _add_rigid_base(root_node)

    #
    cos_nul_state = [0.0, 0.1, 0.1]  # torsion, y_bending, z_bending
    bending_states = [cos_nul_state, cos_nul_state, cos_nul_state]
    list_sections_length = [10, 10, 10]  # SI units

    bending_node = _add_cosserat_state(root_node, bending_states, list_sections_length)

    section_curv_abs = [0, 10, 20, 30]
    frames_curv_abs = [0, 10, 20, 30]
    cosserat_G_frames = [
        [0.0, 0, 0, 0, 0, 0, 1],
        [10.0, 0, 0, 0, 0, 0, 1],
        [20.0, 0, 0, 0, 0, 0, 1],
        [30.0, 0, 0, 0, 0, 0, 1],
    ]
    _add_cosserat_frame(
        base_node,
        bending_node,
        cosserat_G_frames,
        section_curv_abs,
        frames_curv_abs,
        beam_radius,
    )

    return root_node
