# -*- coding: utf-8 -*-
"""Needle insertion into a deformable FEM gel, using the Cosserat plugin.

Rewritten against the *maintained* ``cosserat`` package (python/cosserat):
    - beam           -> cosserat.beam.CosseratBase (+ params dataclasses)
    - FEM gel        -> cosserat.fem.createFemCubeWithParams
    - constraints    -> cosserat.utils.addConstraintPoint (uses PointsManager)
    - controller     -> cosserat.needle.needleController.Animation

The legacy version depended on a stale duplicate ``examples/python3/cosserat``
package (deprecated SofaEngine plugin, broken ``useful.utils`` imports); this
port removes that dependency.

NOTE: the fine collision/constraint wiring (marked "VALIDATE") should be checked
with runSofa; it could not be validated at authoring time.
"""

import os
import sys

# Make the maintained cosserat package win over any installed legacy duplicate.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))

from cosserat.beam import CosseratBase
from cosserat.fem import createFemCubeWithParams
from cosserat.needle.needleController import Animation
from cosserat.needle.params import (
    ContactParams,
    FemParams,
    GeometryParams,
    PhysicsParams,
)
from cosserat.params import (
    BeamGeometryParameters,
    BeamPhysicsParametersNoInertia,
    Parameters,
    SimulationParameters,
)
from cosserat.utils import addConstraintPoint

# Plugins: the "Sofa.Component" meta-plugin pulls in every standard component,
# which keeps this list robust; plus GL rendering and the Cosserat plugin.
pluginList = [
    "Sofa.Component",
    "Sofa.GL.Component.Rendering3D",
    "Cosserat",
]


def _build_beam_params():
    geo = BeamGeometryParameters(
        beam_length=GeometryParams.totalLength,
        nb_section=GeometryParams.nbSections,
        nb_frames=GeometryParams.nbFrames,
        build_collision_model=1,
    )
    phys = BeamPhysicsParametersNoInertia(
        young_modulus=PhysicsParams.youngModulus,
        poisson_ratio=PhysicsParams.poissonRatio,
        beam_mass=PhysicsParams.mass,
        beam_radius=GeometryParams.radius,
        beam_length=GeometryParams.totalLength,
    )
    simu = SimulationParameters(rayleigh_stiffness=PhysicsParams.rayleighStiffness)
    return Parameters(beam_geo_params=geo, beam_physics_params=phys, simu_params=simu)


def createScene(rootNode):
    rootNode.addObject("RequiredPlugin", name="plugins", pluginName=pluginList)

    rootNode.addObject(
        "VisualStyle",
        displayFlags="showVisualModels showBehaviorModels hideCollisionModels "
        "hideBoundingCollisionModels hideForceFields hideInteractionForceFields "
        "hideWireframe showMechanicalMappings",
    )
    rootNode.addObject("CollisionPipeline")
    rootNode.addObject("DefaultVisualManagerLoop")
    rootNode.addObject(
        "RuleBasedContactManager",
        responseParams="mu=0.1",
        response="FrictionContactConstraint",
    )
    rootNode.addObject("BruteForceBroadPhase")
    rootNode.addObject("BVHNarrowPhase")
    rootNode.addObject(
        "LocalMinDistance",
        name="Proximity",
        alarmDistance=0.5,
        contactDistance=ContactParams.contactDistance,
        coneFactor=ContactParams.coneFactor,
        angleCone=0.1,
    )

    rootNode.addObject("FreeMotionAnimationLoop")
    generic = rootNode.addObject(
        "BlockGaussSeidelConstraintSolver",
        tolerance="1e-20",
        maxIterations="500",
        computeConstraintForces=1,
        printLog="0",
    )

    rootNode.gravity.value = [0, 0, 0]
    rootNode.addObject("BackgroundSetting", color="0 0.168627 0.211765")

    solverNode = rootNode.addChild("solverNode")
    solverNode.addObject(
        "EulerImplicitSolver", rayleighStiffness=PhysicsParams.rayleighStiffness
    )
    solverNode.addObject(
        "SparseLDLSolver", name="solver", template="CompressedRowSparseMatrixd"
    )
    solverNode.addObject("GenericConstraintCorrection")

    # ----- The needle (Cosserat beam) -----
    needle = solverNode.addChild(
        CosseratBase(parent=solverNode, params=_build_beam_params(), name="needle")
    )
    # point collision model on the needle frames (node "needleCollision")
    needleCollisionModel = needle.add_point_collision_model("needleCollision")
    # sliding points used to measure needle <-> FEM constraint-point distance
    slidingPoint = needle.add_sliding_points()

    # ----- Base controller -----
    # A free-standing rigid MechanicalObject the user drives with the keyboard;
    # the needle base follows it through a stiff spring (external_rest_shape).
    baseControllerNode = rootNode.addChild("baseController")
    baseControlMO = baseControllerNode.addObject(
        "MechanicalObject",
        template="Rigid3d",
        name="controlMO",
        position=[[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]],
        showObject=True,
        showObjectScale=2,
    )
    needle.rigid_base_node.addObject(
        "RestShapeSpringsForceField",
        name="baseSpring",
        stiffness=1e8,
        angularStiffness=1e8,
        external_rest_shape=baseControlMO.getLinkPath(),
        points=0,
        external_points=0,
        template="Rigid3d",
    )

    # ----- The deformable FEM gel volume -----
    cubeNode = createFemCubeWithParams(rootNode, FemParams)
    gelNode = cubeNode.getChild("gelNode")

    # FEM constraint points (holds the PointsManager that adds/removes points).
    # PointsManager.beamPath is a plain string path (Data<string>): it must NOT
    # start with '@' (which SOFA would parse as a Data link), so strip the
    # leading '@' returned by getLinkPath().
    beam_path = slidingPoint.slidingPointMO.getLinkPath().lstrip("@")
    constraintPointNode = addConstraintPoint(gelNode, beam_path)

    # VALIDATE: ContactListener needs two *collision models*. On the needle side
    # that is the PointCollisionModel "beamColMod" created by addPointsCollision
    # (collisionStats is the MechanicalObject, not a collision model).
    conttactL = rootNode.addObject(
        "ContactListener",
        name="contactListener",
        collisionModel1=cubeNode.gelNode.surfaceNode.surface.getLinkPath(),
        collisionModel2=needleCollisionModel.beamColMod.getLinkPath(),
    )

    # Distance stats between the FEM constraint points and their needle projection
    distanceStatsNode = slidingPoint.addChild("distanceStatsNode")
    constraintPointNode.addChild(distanceStatsNode)
    constraintPoinMo = distanceStatsNode.addObject(
        "MechanicalObject",
        name="distanceStats",
        template="Vec3d",
        position=[],
        listening="1",
        showObject="1",
        showObjectScale="0.1",
    )
    inputVolumeMo = constraintPointNode.constraintPointsMo.getLinkPath()
    inputNeedleMo = slidingPoint.slidingPointMO.getLinkPath()
    outputDistanceMo = distanceStatsNode.distanceStats.getLinkPath()

    # Insertion controller (drives PointsManager add/remove during insertion,
    # and moves the base control MO with the arrow keys)
    rootNode.addObject(
        Animation(
            needle, conttactL, generic, constraintPointNode, rootNode, constraintPoinMo,
            control_mo=baseControlMO,
        )
    )

    distanceStatsNode.addObject(
        "CosseratNeedleSlidingConstraint", name="computeDistanceComponent"
    )
    distanceStatsNode.addObject(
        "DifferenceMultiMapping",
        name="pointsMulti",
        input1=inputVolumeMo,
        lastPointIsFixed=0,
        input2=inputNeedleMo,
        output=outputDistanceMo,
        direction="@../../FramesMO.position",
    )

    return rootNode
