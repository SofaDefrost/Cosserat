# -*- coding: utf-8 -*-
"""Basic scene using Cosserat in SofaPython3.

Based on the work done with SofaPython. See POEMapping.py
"""

from params import NeedleParameters
from cosserat.usefulFunctions import pluginList
from cosserat.createFemRegularGrid import createFemCubeWithParams
from cosserat.cosseratObject import Cosserat
import Sofa
from splib3.numerics import Quat
import sys

sys.path.append('../')

__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com, yinoussa.adagolodjo@inria.fr"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "March 8 2021"


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
        self.threshold = 1.
        self.needleCollisionModel = args[4]
        self.constraintPoints = args[5]
        self.inside = False

        self.rate = 0.2
        self.angularRate = 0.02
        self.tipForce = [0., 0., 0]
        return

    def onAnimateEndEvent(self, event):
        if self.contactListener.getContactPoints() and not self.inside:
            vec = self.contactListener.getContactPoints()[0][1]
            tip = [vec[0], vec[1], vec[2]]

            if self.generic.constraintForces and self.generic.constraintForces[0] > self.threshold:
                # @info 1. Save the entryPoint
                self.entryPoint = tip

                Force = self.generic.constraintForces
                self.tipForce = [Force[0], Force[1], Force[2]]

                # @info 2. deactivate the contact constraint
                self.needleCollisionModel.findData('activated').value = 0

                # @info 3. Add entryPoint point as the first constraint point in FEM
                with self.constraintPoints.position.writeable() as pos:
                    pos.resize(1)
                    pos[0] = self.entryPoint
                    print(f' ====> The tip is : {tip}')
                    print(f' ====> The entryPoint is : {self.entryPoint}')
                    print(f' ====> The constraintPoints is : {pos[0]}')

                self.inside = True

            elif self.tipForce[0] > self.threshold:
                print(
                    "Please activate computeConstraintForces data field inside the GenericConstraint component")
        elif self.inside:
            # 4. todo: add new constraint points inside the volume if needed.
            # todo: This depends on the choice of the algorithm
            #  expl1: one can compare tip position related to the last constraint position inside the volume and
            #  when this > than the constraintDistance we add new constraint point
            # addNewConstraintPoint()

            # 5. todo: If the user is pulling out the needle and the needle tip is behind is before the entryPoint,
            # todo: activated contact constraint.
            # 5.1 self.inside=False
            # 5.2 self.needleCollisionModel.findData('activated').value = 1
            pass

    def onKeypressedEvent(self, event):
        key = event['key']
        if key == "k":  # -
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[5][1] -= self.angularRate

        # ######## Reste rigid position #########
        elif key == "+":  # up
            with self.rigidBaseMO.rest_position.writeable() as posA:
                qOld = Quat()
                for i in range(4):
                    qOld[i] = posA[0][i + 3]
                qNew = Quat.createFromEuler([0., self.angularRate, 0.], 'ryxz')
                qNew.normalize()
                qNew.rotateFromQuat(qOld)
                for i in range(4):
                    posA[0][i + 3] = qNew[i]

        elif key == "-":  # down
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
        elif ord(key) == 20:  # right
            print(
                f' ====> contactListener : {self.contactListener.getContactPoints()}')
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][0] += self.rate

        elif ord(key) == 21:  # down
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][1] -= self.rate

        elif ord(key) == 19:  # up
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
    rootNode.addObject('RuleBasedContactManager',
                       responseParams='mu=0.1', response='FrictionContactConstraint')
    rootNode.addObject('BruteForceBroadPhase')
    rootNode.addObject('BVHNarrowPhase')
    # rootNode.addObject('LocalMinDistance', alarmDistance=1.0, contactDistance=0.01)
    rootNode.addObject('LocalMinDistance', name="Proximity", alarmDistance=0.5,
                       contactDistance=params.contact.contactDistance,
                       coneFactor=params.contact.coneFactor, angleCone=0.1)

    rootNode.addObject('FreeMotionAnimationLoop')
    generic = rootNode.addObject('GenericConstraintSolver', tolerance="1e-20",
                                 maxIterations="500", computeConstraintForces=1, printLog="0")

    gravity = [0, 0, 0]
    rootNode.gravity.value = gravity
    rootNode.addObject('BackgroundSetting', color='1 1 1')
    rootNode.addObject('OglSceneFrame', style="Arrows", alignment="TopRight")
    # ###############
    # New adds to use the sliding Actuator
    ###############
    solverNode = rootNode.addChild('solverNode')
    solverNode.addObject('EulerImplicitSolver',
                         rayleighStiffness=params.Physics.rayleighStiffness)
    solverNode.addObject('SparseLDLSolver', name='solver',
                         template="CompressedRowSparseMatrixd")
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
                                 conttactL, generic, needleCollisionModel, inputFEMCable))
    mappedPointsNode.addObject(
        'CosseratNeedleSlidingConstraint', name="QPConstraint")
    mappedPointsNode.addObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO, lastPointIsFixed=0,
                               input2=inputCableMO, output=outputPointMO, direction="@../../FramesMO.position")
    return rootNode
