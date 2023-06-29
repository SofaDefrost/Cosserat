# -*- coding: utf-8 -*-
"""Basic scene using Cosserat in SofaPython3.

Based on the work done with SofaPython. See POEMapping.py
"""

__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com, yinoussa.adagolodjo@inria.fr"
__version__ = "1.0.0"
__copyright__ = "(c) 2020,Inria"
__date__ = "March 17 2021"

import Sofa
import os
import sys
from stlib3.scene import MainHeader
from splib3.numerics import Quat

sys.path.append('../')
from createFemRegularGrid import createFemCube
from usefulFunctions import BuildCosseratGeometry, AddPointProcess

pluginNameList = 'SofaConstraint SofaDeformable SofaImplicitOdeSolver SofaMeshCollision SofaPreconditioner' \
                 ' SofaGeneralTopology SofaOpenglVisual SofaGeneralRigid SoftRobots SofaSparseSolver' \
                 ' Cosserat SofaBoundaryCondition'


class Animation(Sofa.Core.Controller):

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.rigidBaseMO = args[0]
        self.rateAngularDeformMO = args[1]

        self.rate = 0.2
        self.angularRate = 0.02
        return

    def onKeypressedEvent(self, event):
        key = event['key']
        # ######## Rate angular #########
        if key == "I":  # +
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[5][1] += self.angularRate

        if key == "K":  # -
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[5][1] -= self.angularRate

        if ord(key) == 18:  # left
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][0] -= self.rate

        if ord(key) == 20:  # right
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][0] += self.rate

        if ord(key) == 21:  # down
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][1] -= self.rate

        if ord(key) == 19:  # up
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][1] += self.rate


def createScene(rootNode):

    rootNode.addObject('RequiredPlugin', pluginName=pluginNameList, printLog='0')

    rootNode.addObject('VisualStyle', displayFlags='showBehaviorModels hideCollisionModels hideBoundingCollisionModels '
                                                   'showForceFields hideInteractionForceFields showWireframe')
    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('GenericConstraintSolver', tolerance="1e-20", maxIterations="500", printLog="0")

    gravity = [0, 0, 0]
    rootNode.gravity.value = gravity
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.addObject('OglSceneFrame', style="Arrows", alignment="TopRight")
    # ###############
    # New adds to use the sliding Actuator
    ###############
    cableNode = rootNode.addChild('cableNode')
    cableNode.addObject('EulerImplicitSolver', firstOrder="0", rayleighStiffness="0.1", rayleighMass='0.1')
    cableNode.addObject('SparseLUSolver', name='solver')
    cableNode.addObject('GenericConstraintCorrection')

    # ###############
    # RigidBase
    ###############
    rigidBaseNode = cableNode.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode.addObject('MechanicalObject', template='Rigid3d', name="RigidBaseMO",
                                          position="0 0 0  0 0 0 1", translation="-40. 0. 0.", showObject='1',
                                          showObjectScale='5.')
    rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="50000",
                            angularStiffness="50000", external_points="0", mstate="@RigidBaseMO", points="0",
                            template="Rigid3d")

    #############################################
    # Rate of angular Deformation  (2 sections)
    #############################################
    cosserat_config = {'init_pos': [0., 0., 0.], 'tot_length': 80, 'nbSectionS': 8,
                       'nbFramesF': 16, 'buildCollisionModel': 0, 'beamMass': 0.22}
    [positionS, curv_abs_inputS, longeurS, framesF, curv_abs_outputF, cable_positionF] = \
        BuildCosseratGeometry(cosserat_config)

    rateAngularDeformNode = cableNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject(
        'MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=positionS, showIndices="0")
    rateAngularDeformNode.addObject(
        'BeamHookeLawForceField', crossSectionShape='circular', length=longeurS, radius='2.0', youngModulus='1.e12')
    ################################
    # Animation (to move the dofs) #
    ################################
    rootNode.addObject(Animation(RigidBaseMO, rateAngularDeformMO))

    ##############
    #   Frames   #
    ##############

    # the node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.addChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)
    framesMO = mappedFrameNode.addObject(
        'MechanicalObject', template='Rigid3d', name="FramesMO", position=framesF, showObject='1', showObjectScale='1')

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=curv_abs_inputS,
                              curv_abs_output=curv_abs_outputF, input1=inputMO, input2=inputMO_rigid,
                              output=outputMO, debug='0', max=6.e-2, deformationAxis=1, nonColored="0", radius=5)

    needleCollision = mappedFrameNode.addChild('needleCollision')
    edgeList = "0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16"
    needleCollision.addObject('EdgeSetTopologyContainer', name="Container", position=cable_positionF, edges=edgeList)
    needleCollisionMO = needleCollision.addObject('MechanicalObject', name="cablePos", position=cable_positionF,
                                                  showObject="1", showIndices="0")
    needleCollision.addObject('IdentityMapping')

    # Create FEM Node
    femPos = [" 40.0 0 0 45 0 0 50 0 0 55 0 0 60 0 0 65 0 0  70 0 0  80 0 0 "]
    edgeTrajectory = "0 1 1 2 2 3 3 4 4 5 5 6 6 7 "
    cubeNode = createFemCube(rootNode)
    gelNode = cubeNode.getChild('gelNode')

    Trajectory = gelNode.addChild('Trajectory')
    Trajectory.addObject('VisualStyle', displayFlags='showCollisionModels')
    Trajectory.addObject('EdgeSetTopologyContainer', name="Container", position=femPos, edges=edgeTrajectory)
    Trajectory.addObject('EdgeSetTopologyModifier', name='Modifier')
    Trajectory.addObject('MechanicalObject', name="pointsInFEM", position=femPos, showIndices="1")
    Trajectory.addObject('BarycentricMapping')

    constraintPoints = gelNode.addChild('constraintPoints')
    inputConstraintPointsMO = constraintPoints.addObject('MechanicalObject', name="pointsInFEM", position=femPos,
                                                         showIndices="1")
    constraintPoints.addObject('BarycentricMapping')

    mappedPointsNode = needleCollision.addChild('MappedPoints')
    constraintPoints.addChild(mappedPointsNode)
    mappedPoints = mappedPointsNode.addObject('MechanicalObject', template='Vec3d', name="FramesMO")

    inputCableMO = needleCollisionMO.getLinkPath()
    inputFEMCableMO = inputConstraintPointsMO.getLinkPath()
    outputPointMO = mappedPoints.getLinkPath()

    mappedPointsNode.addObject('CosseratNeedleSlidingConstraint', name="QPConstraint")

    diffMapping = mappedPointsNode.addObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO,
                                             lastPointIsFixed=0, input2=inputCableMO, output=outputPointMO,
                                             direction="@../../FramesMO.position")
    mappedPointsNode.addObject(AddPointProcess(inputConstraintPointsMO, needleCollisionMO, diffMapping))

    return rootNode


def main():
    import SofaRuntime
    import Sofa.Gui
    SofaRuntime.importPlugin("SofaOpenglVisual")
    root = Sofa.Core.Node("root")
    print("0. ====> ")
    createScene(root)
    print("1. ====> ")

    Sofa.Gui.GUIManager.Init("myscene", "qglviewer")
    Sofa.Gui.GUIManager.createGUI(root, __file__)
    Sofa.Gui.GUIManager.SetDimension(1080, 1080)
    Sofa.Gui.GUIManager.MainLoop(root)
    Sofa.Gui.GUIManager.closeGUI()

    print("End of simulation.")


if __name__ == '__main__':
    main()
