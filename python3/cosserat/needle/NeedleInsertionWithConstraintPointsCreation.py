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
sys.path.append('../')
from createFemRegularGrid import createFemCube
from usefulFunctions import BuildCosseratGeometry

from stlib3.scene import MainHeader
from splib3.numerics import Quat


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

        # ######## Reste rigid position #########
        if key == "+":  # up
            with self.rigidBaseMO.rest_position.writeable() as posA:
                qOld = Quat()
                for i in range(0, 4):
                    qOld[i] = posA[0][i+3]
                qNew = Quat.createFromEuler([0., self.angularRate, 0.], 'ryxz')
                qNew.normalize()
                qNew.rotateFromQuat(qOld)
                for i in range(0, 4):
                    posA[0][i+3] = qNew[i]

        if key == "-":  # down
            with self.rigidBaseMO.rest_position.writeable() as posA:
                qOld = Quat()
                for i in range(0, 4):
                    qOld[i] = posA[0][i+3]

                qNew = Quat.createFromEuler([0., -self.angularRate,  0.], 'ryxz')
                qNew.normalize()
                qNew.rotateFromQuat(qOld)
                for i in range(0, 4):
                    posA[0][i+3] = qNew[i]

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

    MainHeader(rootNode, plugins=["SoftRobots", "SofaSparseSolver", "SofaPreconditioner",
                                  "SofaOpenglVisual", "Cosserat", "BeamAdapter", "SofaDeformable",
                                  "SofaImplicitOdeSolver", 'SofaEngine', 'SofaMeshCollision', 'SofaSimpleFem',
                                  'SofaConstraint', 'SofaTopologyMapping'],
               repositoryPaths=[os.getcwd()])

    # rootNode.addObject('VisualStyle', displayFlags='showVisualModels showInteractionForceFields')
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
                                          position="0 0 0  0 0 0 1", showObject='1', showObjectScale='5.')
    rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="50000",
                            angularStiffness="50000", external_points="0", mstate="@RigidBaseMO", points="0",
                            template="Rigid3d")

    ###############
    # Rate of angular Deformation  (2 sections)
    ###############
    [positionS, curv_abs_inputS, longeurS, framesF, curv_abs_outputF, cable_positionF] = \
        BuildCosseratGeometry(nbSection=8, nbFrames=16, totalLength=80)

    rateAngularDeformNode = cableNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject(
        'MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=positionS, showIndices="1")
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

    slidingPoint = mappedFrameNode.addChild('slidingPoint')
    edgeList = "0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16"
    slidingPoint.addObject('EdgeSetTopologyContainer', name="Container", position=cable_positionF, edges=edgeList)
    slidingPointMO = slidingPoint.addObject('MechanicalObject', name="cablePos", position=cable_positionF,
                                            showObject="1", showIndices="0")
    slidingPoint.addObject('NeedleGeometry', name='Needle')  ####  NeedleConstriantPlugin
    slidingPoint.addObject('IdentityMapping')

    # Create FEM Node
    femPos = [" 41.0 0 0 45 0 0 50 0 0 55 0 0 60 0 0 60 0 0  70 0 0  80 0 0 "]
    cubeNode = createFemCube(rootNode)
    gelNode = cubeNode.getChild('gelNode')

    Trajectory = gelNode.addChild('Trajectory')
    Trajectory.addObject('EdgeSetTopologyContainer', name="Container")
    Trajectory.addObject('EdgeSetTopologyModifier', name='Modifier')
    inputFEMCable = Trajectory.addObject('MechanicalObject', name="pointsInFEM", position=femPos, showIndices="1")
    Trajectory.addObject('NeedleTrajectoryGeometry',     name="trajectory", entryDist="1", constraintDist="3",
                            clearTrajectory="true")
    Trajectory.addObject('BarycentricMapping')

    mappedPointsNode = slidingPoint.addChild('MappedPoints')
    Trajectory.addChild(mappedPointsNode)
    mappedPoints = mappedPointsNode.addObject('MechanicalObject', template='Vec3d', position=femPos, name="FramesMO")

    inputCableMO = slidingPointMO.getLinkPath()
    inputFEMCableMO = inputFEMCable.getLinkPath()
    outputPointMO = mappedPoints.getLinkPath()

    mappedPointsNode.addObject('CosseratNeedleSlidingConstraint', name="QPConstraint")
    mappedPointsNode.addObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO, lastPointIsFixed=0,
                               input2=inputCableMO, output=outputPointMO, direction="@../../FramesMO.position")

    return rootNode

