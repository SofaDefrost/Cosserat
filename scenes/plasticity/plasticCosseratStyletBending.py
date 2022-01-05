# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 2021

@author: ckrewcun
"""

# -*- coding: utf-8 -*-
import Sofa
import sys
import os
from CosseratBeamCreator import *  # generateRegularSectionsAndFrames
from CosseratBeamControllers import *  # InsertionController, ColorMapController

# ----- Beam topology ----- #

# /!\ Anatomic mesh in cm /!\

totLength = 55.00
nbSections = 15
nbFrames = 95
beamCreatorResult = generateRegularSectionsAndFrames(totLength, nbSections, nbFrames)

beamSectionLengths = beamCreatorResult['sectionLengths']
beamCurvAbsInput = beamCreatorResult['curvAbsInput']
beamSectionDoFs = beamCreatorResult['sectionDoFs']
beamFrameRigidDoFs = beamCreatorResult['frameRigidDoFs']
beamFrame3DDoFs = beamCreatorResult['frame3DDoFs']
beamFrameEdges = beamCreatorResult['frameEdges']
beamCurvAbsOutput = beamCreatorResult['curvAbsOutput']

# ----- Section geometry ----- #

beamCrossSectionShape = 'circular'
sectionRadius = 3.5e-2  # in cm

# ----- Anatomy mesh file ----- #

anatomyMeshFile = 'Meshes/cylinder_varying_diameter_mesh.stl'

# ----- BeamHookeLawForceField lists ----- #

## Random params
poissonRatio = 0.3
youngModulus = 2.00e4 # in cm
beamYieldStress = 2.20e2
beamPlasticModulus = 34628.0
beamBendingMoment = 3.0e2
## Params in SI:
# youngModulus = 2.00e11 # in cm
# beamYieldStress = 2.20e8
# beamPlasticModulus = 34628588874.0
## Params in cm
# youngModulus = 2.00e9 # in cm
# beamYieldStress = 2.20e6
# beamPlasticModulus = 346285888.7
# beamBendingMoment = 3.0e6
beamYoungModulusList = [youngModulus]*nbSections
beamPoissonRatioList = [poissonRatio]*nbSections

# ----- Miscellaneous parameters ----- #

beamInsertionRate = 0.02
beamInsertionDirection = np.array([1.0, 0., 0.])
forceEstimationSpringStiffness = 1e8
bendingIndices = nbSections-1

# Friction parameter contact
frictionCoefficient = 0.1
tolerance = 1e-4
frictionCoefficientParam = "mu="+str(frictionCoefficient)+"&tol="+str(tolerance)

pluginNameList = 'SofaConstraint SofaDeformable SofaImplicitOdeSolver SofaMeshCollision SofaPreconditioner' \
                 ' SofaGeneralTopology SofaOpenglVisual SofaGeneralRigid SoftRobots SofaSparseSolver' \
                 ' CosseratPlugin SofaBoundaryCondition SofaGeneralAnimationLoop'  # BeamAdapter


def createScene(rootNode):
    rootNode.addObject('RequiredPlugin', pluginName=pluginNameList, printLog='0')
    rootNode.addObject('VisualStyle', displayFlags='showBehaviorModels hideCollisionModels '
                                                   'hideVisualModels'
                                                   'hideBoundingCollisionModels showForceFields '
                                                   'showInteractionForceFields hideWireframe')
    rootNode.addObject('DefaultVisualManagerLoop')
    rootNode.findData('gravity').value = [9.810, 0., 0.]
    rootNode.findData('dt').value = 0.01

    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('DefaultPipeline', verbose="0")
    rootNode.addObject('BruteForceBroadPhase', name="BroadPhase")
    rootNode.addObject('BVHNarrowPhase', name="NarrowPhase")
    rootNode.addObject('DefaultContactManager', response="FrictionContact", responseParams=frictionCoefficientParam)
    rootNode.addObject('LocalMinDistance', name="Proximity", alarmDistance=1.6*sectionRadius, contactDistance=sectionRadius)
    rootNode.addObject('GenericConstraintSolver', tolerance=1e-10, maxIterations=1e3)
    rootNode.addObject('Camera', position="-35 0 280", lookAt="0 0 0")

    # -------------------------------------#
    # -----       Control point      ----- #
    # -------------------------------------#

    controlPointNode = rootNode.addChild('controlPointNode')
    controlPointMO = controlPointNode.addObject('MechanicalObject', template='Rigid3d', name="controlPointMO",
                                                position="0 0 0  0 0 0. 1", showObject='1', showObjectScale='0.1')

    # -------------------------------------#
    # -----          Solver          ----- #
    # -------------------------------------#

    cosseratBeamNode = rootNode.addChild('cosseratBeamNode')
    cosseratBeamNode.addObject('EulerImplicitSolver', printLog="false", rayleighStiffness="0.043",
                                rayleighMass="0.01")
    cosseratBeamNode.addObject('SparseLDLSolver', name='solverImplant', template='CompressedRowSparseMatrixd')
    cosseratBeamNode.addObject('GenericConstraintCorrection')

    # -------------------------------------#
    # -----        Rigid base        ----- # (beam base that is not part of the beam)
    # -------------------------------------#

    rigidBaseNode = cosseratBeamNode.addChild('rigidBaseNode')

    RigidBaseMO = rigidBaseNode.addObject('MechanicalObject', template='Rigid3d', name="RigidBaseMO",
                                          position="0 0 0  0 0 0. 1", showObject='1', showObjectScale='0.001')
    rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="1.e6", angularStiffness="1.e6",
                            external_rest_shape="@../../controlPointNode/controlPointMO", external_points="0",
                            mstate="@RigidBaseMO", points="0", template="Rigid3d")

    # ----------------------------------------#
    # ----- Rate of angular deformation ----- #
    # ----------------------------------------#

    rateAngularDeformNode = cosseratBeamNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject('MechanicalObject', template='Vec3d',
                                                          position=beamSectionDoFs, name='rateAngularDeformMO')
    rateAngularDeformNode.addObject('BeamPlasticLawForceField', name="beamHookeLaw",
                                    crossSectionShape=beamCrossSectionShape,
                                    radius=sectionRadius, youngModulus='2e4', varianteSections="true",
                                    length=beamCreatorResult['sectionLengths'], poissonRatioList=beamPoissonRatioList,
                                    youngModulusList=beamYoungModulusList,
                                    initialYieldStress=beamYieldStress, plasticModulus=beamPlasticModulus)

    rateAngularDeformNode.addObject('ConstantForceField', name='bendingMoment', indices=bendingIndices,
                                    forces=[0., 0., 0.])
    fixedSections = list(range(nbSections - 2))
    rateAngularDeformNode.addObject('FixedConstraint', name="fixation", indices=fixedSections)

    # -------------------------------------#
    # -----          Frames          ----- #
    # -------------------------------------#

    # the node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.addChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)
    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d', name="FramesMO",
                                         position=beamFrameRigidDoFs,
                                         showObject='1', showObjectScale='0.1')
    mappedFrameNode.addObject('UniformMass', totalMass="0.022", showAxisSizeFactor='0.')

    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()
    mappedFrameNode.addObject('DiscreteCosseratMapping', name="mapping", input1=inputMO, input2=inputMO_rigid,
                              output=outputMO, debug='0', nonColored=False,
                              curv_abs_input=beamCurvAbsInput, curv_abs_output=beamCurvAbsOutput,
                              forcefield='@../../rateAngularDeform/beamHookeLaw')
    cosseratBeamNode.addObject('MechanicalMatrixMapper', template='Vec3,Rigid3', object1=inputMO, object2=inputMO_rigid,
                                nodeToParse=mappedFrameNode.getLinkPath())

    # -------------------------------------#
    # -----   Beam collision model   ----- #
    # -------------------------------------#

    CollisInstrumentCombined = mappedFrameNode.addChild('CollisInstrumentCombined')
    CollisInstrumentCombined.addObject('EdgeSetTopologyContainer', name="collisEdgeSet", position=beamCreatorResult['frame3DDoFs'],
                                       edges=beamCreatorResult['frameEdges'])
    CollisInstrumentCombined.addObject('EdgeSetTopologyModifier', name="colliseEdgeModifier")
    CollisInstrumentCombined.addObject('MechanicalObject', name="CollisionDOFs")
    CollisInstrumentCombined.addObject('LineCollisionModel', bothSide="1", group='2', proximity="0.01")
    CollisInstrumentCombined.addObject('PointCollisionModel', bothSide="1", group='2')
    CollisInstrumentCombined.addObject('IdentityMapping', name="mapping")

    # -------------------------------------#
    # -----           Tube           ----- #
    # -------------------------------------#

    anatomyRefNode = rootNode.addChild('anatomyRefNode')
    anatomyRefNode.addObject('MeshSTLLoader', name='loader', filename=anatomyMeshFile)

    anatomyRefNode.addObject('MechanicalObject', template='Vec3d', name='anatomyRefDoFs', position="@loader.position")
    # anatomyRefNode.addObject('OglModel', src="@loader", color="red")

    anatomyNode = rootNode.addChild('anatomyNode')
    anatomyNode.addObject('VisualStyle', displayFlags='showCollisionModels showWireframe')
    anatomyNode.addObject('EulerImplicitSolver', name='odesolver', firstOrder="0", rayleighMass="0.1",
                          rayleighStiffness="0.1")
    anatomyNode.addObject('ShewchukPCGLinearSolver', name='linearSolver', iterations='500', tolerance='1.0e-8',
                          preconditioners="precond")
    anatomyNode.addObject('SparseLDLSolver', name='precond', template='CompressedRowSparseMatrixMat3x3d,')
    anatomyNode.addObject('MeshTopology', src='@../anatomyRefNode/loader')
    anatomyNode.addObject('MechanicalObject', name='anatomyDoFs', template='Vec3d', showIndices='false',
                          showIndicesScale='4e-5', rx='0', printLog="0", position="@../anatomyRefNode/loader.position")
    anatomyNode.addObject('TriangleCollisionModel', group='1', bothSide="1")
    anatomyNode.addObject('LineCollisionModel', group='1', bothSide="1")
    anatomyNode.addObject('PointCollisionModel', group='1', bothSide="1")
    # anatomyNode.addObject('OglModel', src="@../anatomyRefNode/loader", color="green")

    anatomyNode.addObject('RestShapeSpringsForceField', stiffness=forceEstimationSpringStiffness,
                          external_rest_shape="@../anatomyRefNode/anatomyRefDoFs",
                          mstate="@anatomyDoFs")

    anatomyNode.addObject('GenericConstraintCorrection', solverName='precond')

    # -------------------------------------#
    # -----    Python controllers    ----- #
    # -------------------------------------#

    switchFrameDistance = 0.75 * totLength / nbFrames
    incrementAngle = 5.0  # in degrees
    rootNode.addObject(InsertionController(rootNode=rootNode, insertionRate=beamInsertionRate,
                                           insertionDirection=beamInsertionDirection, name="InsertionController"))
    # rootNode.addObject(InteractiveInsertionController(rootNode=rootNode, initFrameId=nbFrames - 1,
    #                                                   insertionRate=beamInsertionRate,
    #                                                   insertionDirection=beamInsertionDirection,
    #                                                   switchFrameDistance=switchFrameDistance,
    #                                                   incrementAngle=incrementAngle,
    #                                                   name="InsertionController"))

    rootNode.addObject(BendingController(rootNode=rootNode, cosseratNode=rateAngularDeformNode,
                                         frameNode=mappedFrameNode, bendingMoment=beamBendingMoment,
                                         fixedIndices=fixedSections, momentAxis='y', name="BendingController"))
