# -*- coding: utf-8 -*-
"""Basic scene using Cosserat in SofaPython3.

Based on the work done with SofaPython. See POEMapping.py
"""

from splib3.numerics import Quat
import sys

sys.path.append('../')

__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com, yinoussa.adagolodjo@inria.fr"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "March 8 2021"

import Sofa
from cosserat.cosseratObject import Cosserat
from cosserat.createFemRegularGrid import createFemCubeWithParams
from cosserat.usefulFunctions import pluginList
from params import NeedleParameters


params = NeedleParameters()

needleGeometryConfig = {'init_pos': [0., 0., 0.], 'tot_length': params.Geometry.totalLength,
                        'nbSectionS': params.Geometry.nbSections, 'nbFramesF': params.Geometry.nbFrames,
                        'buildCollisionModel': 1, 'beamMass': params.Physics.mass}


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
                for i in range(4):
                    qOld[i] = posA[0][i + 3]
                qNew = Quat.createFromEuler([0., self.angularRate, 0.], 'ryxz')
                qNew.normalize()
                qNew.rotateFromQuat(qOld)
                for i in range(4):
                    posA[0][i + 3] = qNew[i]

        if key == "-":  # down
            with self.rigidBaseMO.rest_position.writeable() as posA:
                qOld = Quat()
                for i in range(4):
                    qOld[i] = posA[0][i + 3]

                qNew = Quat.createFromEuler(
                    [0., -self.angularRate, 0.], 'ryxz')
                qNew.normalize()
                qNew.rotateFromQuat(qOld)
                for i in range(0, 4):
                    posA[0][i + 3] = qNew[i]

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
    rootNode.addObject(
        'RequiredPlugin', pluginName=pluginList, printLog='0')

    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideCollisionModels '
                                                   'hideBoundingCollisionModels hireForceFields '
                                                   'hideInteractionForceFields hideWireframe showMechanicalMappings')
    rootNode.addObject('DefaultPipeline')
    rootNode.addObject("DefaultVisualManagerLoop")
    rootNode.addObject('RuleBasedContactManager', responseParams='mu=0.1', response='FrictionContactConstraint')
    rootNode.addObject('BruteForceBroadPhase')
    rootNode.addObject('BVHNarrowPhase')
    # rootNode.addObject('LocalMinDistance', alarmDistance=1.0, contactDistance=0.01)
    rootNode.addObject('LocalMinDistance', name="Proximity", alarmDistance=2.,
                       contactDistance=params.contact.contactDistance, coneFactor=params.contact.coneFactor,
                       angleCone=0.1)

    rootNode.addObject('FreeMotionAnimationLoop')
    # rootNode.addObject('CollisionPipeline', verbose="0")
    # rootNode.addObject('BruteForceDetection', name="N2")
    # rootNode.addObject('DefaultContactManager',
    #                    response="FrictionContactConstraint", responseParams=params.contact.dataMu)
    # rootNode.addObject('LocalMinDistance', name="Proximity", alarmDistance=params.contact.alarmDistance,
    #                    contactDistance=params.contact.contactDistance, coneFactor=params.contact.coneFactor,
    #                    angleCone=params.contact.angleCone)

    rootNode.addObject('GenericConstraintSolver', tolerance="1e-20", maxIterations="500", printLog="0")

    gravity = [0, 0, 0]
    rootNode.gravity.value = gravity
    rootNode.addObject('BackgroundSetting', color='1 1 1')
    rootNode.addObject('OglSceneFrame', style="Arrows", alignment="TopRight")
    # ###############
    # New adds to use the sliding Actuator
    ###############
    solverNode = rootNode.addChild('solverNode')
    solverNode.addObject('EulerImplicitSolver', rayleighStiffness=params.Physics.rayleighStiffness)
    solverNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")
    solverNode.addObject('GenericConstraintCorrection')

    needle = solverNode.addChild(
        Cosserat(parent=solverNode, cosseratGeometry=needleGeometryConfig, radius=params.Geometry.radius,
                 name="needle", youngModulus=params.Physics.youngModulus, poissonRatio=params.Physics.poissonRatio,
                 rayleighStiffness=params.Physics.rayleighStiffness))
    collisionModel = needle.addPointCollisionModel()
    slidingPoint = needle.addSlidingPoints()

    # Create FEM Node
    # TODO: Where we handle Sliding constraints,
    #  we need to be able added or remove dynamically
    # TODO: @Paul is in charge of creating this in python
    #
    femPos = []
    cubeNode = createFemCubeWithParams(rootNode, params.FemParams)
    gelNode = cubeNode.getChild('gelNode')
    femPoints = gelNode.addChild('femPoints')
    inputFEMCable = femPoints.addObject(
        'MechanicalObject', name="pointsInFEM", position=femPos, showIndices="1")
    femPoints.addObject('BarycentricMapping')

    mappedPointsNode = slidingPoint.addChild('MappedPoints')
    femPoints.addChild(mappedPointsNode)
    mappedPoints = mappedPointsNode.addObject(
        'MechanicalObject', template='Vec3d', position=femPos, name="FramesMO")

    inputCableMO = slidingPoint.slidingPointMO.getLinkPath()
    inputFEMCableMO = inputFEMCable.getLinkPath()
    outputPointMO = mappedPoints.getLinkPath()

    rootNode.addObject(Animation(needle.rigidBaseNode.RigidBaseMO, needle.cosseratCoordinateNode.cosseratCoordinateMO))

    mappedPointsNode.addObject(
        'CosseratNeedleSlidingConstraint', name="QPConstraint")
    mappedPointsNode.addObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO, lastPointIsFixed=0,
                               input2=inputCableMO, output=outputPointMO, direction="@../../FramesMO.position")
    return rootNode

