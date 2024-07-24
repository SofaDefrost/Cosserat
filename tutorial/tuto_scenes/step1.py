# -*- coding: utf-8 -*-

import Sofa
from useful.params import Parameters

from math import sin, cos, sqrt, pi

# import os
# path = os.path.dirname(os.path.abspath(__file__))+'/mesh/'
#
#
# _tension = 0.0
# class TensionComputing(Sofa.PythonScriptController):
#     def initGraph(self, node):
#         self.tension = 500
#         self.node = node;
#         self.BeamHookeLawForce = self.node.getObject('BeamHookeLawForce')
#
#     def onBeginAnimationStep(self, dt):
#         self.tension = self.tension + 8000 * dt;
#         self.BeamHookeLawForce.findData('tension').value = self.tension

stiffness_param = 1.e10
beam_radius = 1.
# params = Parameters(beamGeoParams=BeamGeometryParameters(init_pos=[0, 0, 0]))


def _add_rigid_base(p_node):
    rigid_base_node = p_node.addChild('rigid_base')
    rigid_base_node.addObject('MechanicalObject', template='Rigid3d', name="cosserat_base_mo",
                              position="0 0 0  0 0 0. 1",
                              showObject=1, showObjectScale='0.1')
    rigid_base_node.addObject('RestShapeSpringsForceField', name='spring', stiffness=stiffness_param,
                              angularStiffness=stiffness_param, external_points="0", mstate="@cosserat_base_mo",
                              points="0", template="Rigid3d")
    return rigid_base_node


def _add_cosserat_state(p_node, bending_states, list_sections_length, _radius=2.):
    cosserat_coordinate_node = p_node.addChild('cosseratCoordinate')
    cosserat_coordinate_node.addObject('MechanicalObject', template='Vec3d', name='cosserat_state',
                                       position=bending_states)
    cosserat_coordinate_node.addObject('BeamHookeLawForceField', crossSectionShape='circular',
                                       length=list_sections_length, radius=2., youngModulus=1.e4,
                                       poissonRatio=0.4)
    return cosserat_coordinate_node


def _add_cosserat_frame(p_node, _bending_node, framesF, _section_curv_abs, _frame_curv_abs, _radius, _beam_mass=0.0):
    cosserat_in_Sofa_frame_node = p_node.addChild('cosserat_in_Sofa_frame_node')

    _bending_node.addChild(cosserat_in_Sofa_frame_node)
    frames_mo = cosserat_in_Sofa_frame_node.addObject('MechanicalObject', template='Rigid3d',
                                                      name="FramesMO", position=framesF, showIndices=1, showObject=1,
                                                      showObjectScale=0.8)

    cosserat_in_Sofa_frame_node.addObject('UniformMass', totalMass=_beam_mass)

    cosserat_in_Sofa_frame_node.addObject('DiscreteCosseratMapping', curv_abs_input=_section_curv_abs,
                                          curv_abs_output=_frame_curv_abs, name='cosseratMapping',
                                          input1=_bending_node.cosserat_state.getLinkPath(),
                                          input2=p_node.cosserat_base_mo.getLinkPath(),
                                          output=frames_mo.getLinkPath(), debug=0, radius=_radius)
    return cosserat_in_Sofa_frame_node



def createScene(root_node):
    #
    base_node = _add_rigid_base(root_node)

    #
    cos_nul_state = [0.0, 0.0, 0.0]  # torsion, y_bending, z_bending
    bending_states = [cos_nul_state, cos_nul_state, cos_nul_state]
    list_sections_length = [10, 10, 10]
    bending_node = _add_cosserat_state(root_node, cos_nul_state, list_sections_length)

    section_curv_abs = [0, 10, 20, 30]
    frames_curv_abs = [0., 5, 10, 15, 20, 25, 30]
    cosserat_G_frames = [[0., 0, 0,  0, 0, 0, 1], [5., 0, 0, 0, 0, 0, 1], [10., 0, 0,  0, 0, 0, 1],
                         [15, 0, 0, 0, 0, 0, 1], [20., 0, 0,  0, 0, 0, 1], [25., 0, 0,  0, 0, 0, 1],
                         [30., 0, 0,  0, 0, 0, 1]]

    _add_cosserat_frame(base_node, bending_node, cosserat_G_frames, section_curv_abs, frames_curv_abs,
                        beam_radius)

    return root_node

    ###############
    ## Rate of angular Deformation  (2 sections)
    ###############

    distance1 = [0.0, 0.2, 0.0]
    distance2 = [0.0, 0.2, 0.0]
    distance3 = [0.0, 0.2, 0.0]
    _distance = [distance1, distance2, distance3]

    ddistance1 = [0.0, 0.0, 0.0]
    ddistance2 = [0.0, 0.0, 0.0]
    ddistance3 = [0.0, 0.0, 0.0]
    _ddistance = [ddistance1, ddistance2, ddistance3]

    rateAngularDeformNode = rootNode.createChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.createObject('MechanicalObject', template='Vec3d',
                                                             name='rateAngularDeformMO', position=pos,
                                                             velocity='0 0 0 0 0 0 0 0 0', length='10 10 10', )
    # (2 series of 3 angles for 2 sections. we suppose that the lenght is 10 for each)
    # BeamHookeLawForce = rateAngularDeformNode.createObject('CosseratInternalActuation',
    # crossSectionShape='circular', length='10 10 10', radius='0.5', youngModulus='5e6')
    BeamHookeLawForce = rateAngularDeformNode.createObject('CosseratInternalActuation', name="BeamHookeLawForce",
                                                           crossSectionShape='circular', length='10 10 10',
                                                           radius='0.5',
                                                           youngModulus='1e6', distance=_distance, ddistance=_ddistance,
                                                           tension=_tension)
    rateAngularDeformNode.createObject('PythonScriptController', classname="TensionComputing")

    ##############
    ## Frames
    ##############
    # the node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.createChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)
    framesMO = mappedFrameNode.createObject('MechanicalObject', template='Rigid3d', name="FramesMO",
                                            position="0.5 0 0  0 0 0 1  5 0 0 0 0 0 1   10 0 0  0 0 0 1   15 0 0 0 0 0 1  20 0 0  0 0 0 1 25 0 0  0 0 0 1 30 0 0  0 0 0 1",
                                            showObject='1', showObjectScale='1')

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    # one output: FramesMO

    inputMO = rateAngularDeformMO.getLinkPath()  # + " " + RigidBaseMO.getLinkPath()
    # inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()
    # TODO:
    mappedFrameNode.createObject('DiscretCosseratMapping', curv_abs_input='0 10 20 30',
                                 curv_abs_output='0.5 5 10 15 20 25 30', input1=inputMO, input2=inputMO_rigid,
                                 output=outputMO, debug='0')

    #### CylinderGridTop
    CylinderCollision = mappedFrameNode.createChild('CylinderCollision')

    # CylinderCollision.createObject('MeshSTLLoader', filename=path+'trunk.stl', name='loader', rotation='0 90 0', scale='0.155')
    CylinderCollision.createObject('CylinderGridTopology', name="loader", nx="8", ny="8", nz="20", length="30",
                                   radius="0.5", axis="1 0 0")
    CylinderCollision.createObject('Mesh', src='@loader')
    CylinderCollision.createObject('MechanicalObject', template='Vec3d')
    CylinderCollision.createObject('Triangle')
    CylinderCollision.createObject('SkinningMapping', nbRef='2')

    # rootNode.createObject('BilateralInteractionConstraint', template='Rigid3d', object2='@rigidBase/MappedFrames/FramesMO', object1='@targetPos/target', first_point='0', second_point='6')

    return rootNode
