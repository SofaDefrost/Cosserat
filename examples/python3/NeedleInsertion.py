# -*- coding: utf-8 -*-
"""Basic scene using Cosserat in SofaPython3.

Based on the work done with SofaPython. See POEMapping.py
"""

from cosserat.needle.needleController import Animation
from cosserat.needle.params import NeedleParameters, GeometryParams, PhysicsParams, FemParams, ContactParams
from cosserat.usefulFunctions import pluginList
from cosserat.createFemRegularGrid import createFemCubeWithParams
from cosserat.cosseratObject import Cosserat
from cosserat.utils import addConstraintPoint
import sys

# params = NeedleParameters()

sys.path.append('../')

nbFrames = GeometryParams.nbFrames
needleGeometryConfig = {'init_pos': [0., 0., 0.], 'tot_length': GeometryParams.totalLength,
                        'nbSectionS': GeometryParams.nbSections, 'nbFramesF': nbFrames,
                        'buildCollisionModel': 1, 'beamMass': PhysicsParams.mass}


def createScene(rootNode):

    rootNode.addObject('RequiredPlugin', name='plugins', pluginName=[pluginList])

    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideCollisionModels '
                                                   'hideBoundingCollisionModels hideForceFields '
                                                   'hideInteractionForceFields hideWireframe showMechanicalMappings')
    rootNode.addObject('CollisionPipeline')
    rootNode.addObject("DefaultVisualManagerLoop")
    rootNode.addObject('RuleBasedContactManager',
                       responseParams='mu=0.1', response='FrictionContactConstraint')
    rootNode.addObject('BruteForceBroadPhase')
    rootNode.addObject('BVHNarrowPhase')
    # rootNode.addObject('LocalMinDistance', alarmDistance=1.0, contactDistance=0.01)
    rootNode.addObject('LocalMinDistance', name="Proximity", alarmDistance=0.5,
                       contactDistance=ContactParams.contactDistance,
                       coneFactor=ContactParams.coneFactor, angleCone=0.1)

    rootNode.addObject('FreeMotionAnimationLoop')
    generic = rootNode.addObject('GenericConstraintSolver', tolerance="1e-20",
                                 maxIterations="500", computeConstraintForces=1, printLog="0")

    gravity = [0, 0, 0]
    rootNode.gravity.value = gravity
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    # rootNode.addObject('OglSceneFrame', style="Arrows", alignment="TopRight")
    # ###############
    # New adds to use the sliding Actuator
    ###############
    solverNode = rootNode.addChild('solverNode')
    solverNode.addObject('EulerImplicitSolver',
                         rayleighStiffness=PhysicsParams.rayleighStiffness)
    solverNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")
    solverNode.addObject('GenericConstraintCorrection')

    needle = solverNode.addChild(
        Cosserat(parent=solverNode, cosseratGeometry=needleGeometryConfig, radius=GeometryParams.radius,
                 name="needle", youngModulus=PhysicsParams.youngModulus, poissonRatio=PhysicsParams.poissonRatio,
                 rayleighStiffness=PhysicsParams.rayleighStiffness))
    needleCollisionModel = needle.addPointCollisionModel("needleCollision")

    # These state is mapped on the needle and used to compute the distance between the needle and the
    # FEM constraint points
    slidingPoint = needle.addSlidingPoints()

    # -----------------
    # Start the volume definition
    # -----------------
    cubeNode = createFemCubeWithParams(rootNode, FemParams)
    gelNode = cubeNode.getChild('gelNode')
    # FEM constraint points
    constraintPointNode = addConstraintPoint(
        gelNode, slidingPoint.getLinkPath())

    # @info : This is the constraint point that will be used to compute the distance between the needle and the volume
    conttactL = rootNode.addObject('ContactListener', name="contactListener",
                                   collisionModel1=cubeNode.gelNode.surfaceNode.surface.getLinkPath(),
                                   collisionModel2=needleCollisionModel.collisionStats.getLinkPath())

    # These stats will represents the distance between the contraint point in the volume and
    # their projection on the needle
    # It 's also important to say that the x direction is not taken into account
    distanceStatsNode = slidingPoint.addChild('distanceStatsNode')
    constraintPointNode.addChild(distanceStatsNode)
    constraintPoinMo = distanceStatsNode.addObject('MechanicalObject', name="distanceStats", template="Vec3d",
                                                   position=[], listening="1", showObject="1", showObjectScale="0.1")
    inputVolumeMo = constraintPointNode.constraintPointsMo.getLinkPath()
    inputNeedleMo = slidingPoint.slidingPointMO.getLinkPath()
    outputDistanceMo = distanceStatsNode.distanceStats.getLinkPath()

    # ---------------------------------------------------
    # @info: Start controller node
    rootNode.addObject(Animation(needle, conttactL, generic,
                       constraintPointNode, rootNode, constraintPoinMo))

    distanceStatsNode.addObject(
        'CosseratNeedleSlidingConstraint', name="computeDistanceComponent")
    distanceStatsNode.addObject('DifferenceMultiMapping', name="pointsMulti", input1=inputVolumeMo, lastPointIsFixed=0,
                                input2=inputNeedleMo, output=outputDistanceMo, direction="@../../FramesMO.position")
