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
        self.contactListener = args[2]
        self.generic = args[3]
        self.entryPoint = []
        self.threshold = 4.
        self.needleCollisionModel = args[4]
        self.constraintPoints = args[5]
        self.inside = False

        self.rate = 0.2
        self.angularRate = 0.02
        return

    def onAnimateEndEvent(self, event):
        if self.contactListener.getContactPoints() and not self.inside:
            vec = self.contactListener.getContactPoints()[0][1]
            tip = [vec[0], vec[1], vec[2]]
            # print(f' ====> The tip is : {tip}')
            # print(f' ====> The force is : {self.generic.constraintForces[0]}')

            if self.generic.constraintForces and self.generic.constraintForces[0] > self.threshold:
                # @info 1. Save the entryPoint
                self.entryPoint = tip
                print(f' ====> The entryPoint is : {self.entryPoint}')
                # @info 2. deactivate the contact constraint
                # todo: add the code to deactivate the contact constraint
                self.needleCollisionModel.findData('activated').value = 0
                # @info 3. Add entryPoint point as the first sliding constraint
                # todo: add the code according to #3
                # with self.constraintPoints.position.writeable() as pos:
                #     print(f'pos is : {dir(self.constraintPoints.position)}')
                #     pos = self.entryPoint
                self.inside = True

            elif self.generic.constraintForces[0] > self.threshold:
                print("Please activate computeConstraintForces data field inside the GenericConstraint component")
        elif self.inside:
            # 4. todo: add new constraint point inside the volume if needed.
            # todo: This depends on the choice of algorithm
            #  expl1: one can compare tip position related to the last constraint position inside the volume and
            #  when this > than the constraintDistance we add new constraint point
            # addNewConstraintPoint()

            # 5. todo: If the user is pulling out the needle and the needle tip is behind is before the entryPoint,
            # todo: activated the contact constraint.
            # 5.1 self.inside=False
            # 5.2 self.needleCollisionModel.findData('activated').value = 1
            pass

    def onKeypressedEvent(self, event):
        key = event['key']
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
                for i in range(4):
                    posA[0][i + 3] = qNew[i]

        if ord(key) == 18:  # left
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][0] -= self.rate

        if ord(key) == 20:  # right
            print(f' ====> contactListener : {self.contactListener.getContactPoints()}')
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
    rootNode.addObject('LocalMinDistance', name="Proximity", alarmDistance=0.5,
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

    generic = rootNode.addObject('GenericConstraintSolver', tolerance="1e-20", maxIterations="500",
                                 computeConstraintForces=1, printLog="0")

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
    needleCollisionModel = needle.addPointCollisionModel()
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

    conttactL = rootNode.addObject('ContactListener', name="contactListener",
                                   collisionModel1=cubeNode.gelNode.surfaceNode.surface.getLinkPath(),
                                   collisionModel2=needleCollisionModel.pointColli.getLinkPath())
    rootNode.addObject(Animation(needle.rigidBaseNode.RigidBaseMO, needle.cosseratCoordinateNode.cosseratCoordinateMO,
                                 conttactL, generic, needleCollisionModel,inputFEMCable))
    mappedPointsNode.addObject(
        'CosseratNeedleSlidingConstraint', name="QPConstraint")
    mappedPointsNode.addObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO, lastPointIsFixed=0,
                               input2=inputCableMO, output=outputPointMO, direction="@../../FramesMO.position")
    return rootNode
