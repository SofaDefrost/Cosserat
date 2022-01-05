# -*- coding: utf-8 -*-
"""
Created on Thu Nov 4 2021

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

totLength = 20.0
nbSections = 5
nbFrames = 34
beamCreatorResult = generateRegularSectionsAndFrames(totLength, nbSections, nbFrames)

beamSectionLengths = beamCreatorResult['sectionLengths']
beamCurvAbsInput = beamCreatorResult['curvAbsInput']
beamSectionDoFs = beamCreatorResult['sectionDoFs']
beamFrameRigidDoFs = beamCreatorResult['frameRigidDoFs']
beamFrame3DDoFs = beamCreatorResult['frame3DDoFs']
beamFrameEdges = beamCreatorResult['frameEdges']
beamCurvAbsOutput = beamCreatorResult['curvAbsOutput']

totalMass = 0.022

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

    controlPointPos = np.array([1., 0., 0., 0., 0., 0., 1.])
    controlPointNode = rootNode.addChild('controlPointNode')
    controlPointNode.addObject('MechanicalObject', template='Rigid3d', name="controlPointMO",
                               position=controlPointPos, showObject='1', showObjectScale='0.5')

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
    # rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="1.e6", angularStiffness="1.e6",
    #                         external_rest_shape="@../../controlPointNode/controlPointMO", external_points="0",
    #                         mstate="@RigidBaseMO", points="0", template="Rigid3d")

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

    # rateAngularDeformNode.addObject('ConstantForceField', name='bendingMoment', indices=bendingIndices,
    #                                 forces=[0., 0., 0.])
    # fixedSections = list(range(nbSections - 2))
    # rateAngularDeformNode.addObject('FixedConstraint', name="fixation", indices=fixedSections)

    # -------------------------------------#
    # -----          Frames          ----- #
    # -------------------------------------#

    # the node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.addChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)
    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d', name="FramesMO",
                                         position=beamFrameRigidDoFs,
                                         showObject='1', showObjectScale='0.1')
    mappedFrameNode.addObject('UniformMass', totalMass=totalMass, showAxisSizeFactor='0.')

    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()
    mappedFrameNode.addObject('DiscreteCosseratMapping', name="mapping", input1=inputMO, input2=inputMO_rigid,
                              output=outputMO, debug='0', nonColored=False,
                              curv_abs_input=beamCurvAbsInput, curv_abs_output=beamCurvAbsOutput,
                              forcefield='@../../rateAngularDeform/beamHookeLaw')
    cosseratBeamNode.addObject('MechanicalMatrixMapper', name='mechanicalMatrixMapper', template='Vec3,Rigid3',
                               object1=inputMO, object2=inputMO_rigid, nodeToParse=mappedFrameNode.getLinkPath())

    # Second frame nodes, used by the smooth controller
    mappedFrameNode2 = rigidBaseNode.addChild('MappedFrames2')
    rateAngularDeformNode.addChild(mappedFrameNode2)
    frameInitId = 5
    frameInitPos = beamFrameRigidDoFs[frameInitId]x
    frames2MO = mappedFrameNode2.addObject('MechanicalObject', template='Rigid3d', name="Frames2MO",
                                           position=frameInitPos,
                                           showObject='1', showObjectScale='0.9')
    outputMO2 = frames2MO.getLinkPath()
    mappedFrameNode2.addObject('DiscreteCosseratMapping', name="controlMapping", input1=inputMO, input2=inputMO_rigid,
                               output=outputMO2, debug='0', nonColored=False,
                               curv_abs_input=beamCurvAbsInput, curv_abs_output=frameInitPos[0],
                               forcefield='@../../rateAngularDeform/beamHookeLaw')
    # TO DO : is the second MechanicalMatrixMapper correct ?
    cosseratBeamNode.addObject('MechanicalMatrixMapper', name='controlMechanicalMatrixMapper', template='Vec3,Rigid3',
                               object1=inputMO, object2=inputMO_rigid, nodeToParse=mappedFrameNode2.getLinkPath())
    mappedFrameNode2.addObject('UniformMass', totalMass=totalMass/nbFrames, showAxisSizeFactor='0.')
    mappedFrameNode2.addObject('RestShapeSpringsForceField', name='controlSpring', stiffness="1.e6",
                               angularStiffness="1.e6", external_rest_shape="@../../../controlPointNode/controlPointMO",
                               external_points=0, mstate="@FramesMO", points=0, template="Rigid3d", drawSpring=1)

    # -------------------------------------#
    # -----    Python controllers    ----- #
    # -------------------------------------#

    incrementAngle = 5.0  # in degrees
    incrementDistance = 1.0  # in cm in this scene
    # rootNode.addObject(InsertionController(rootNode=rootNode, insertionRate=beamInsertionRate,
    #                                        insertionDirection=beamInsertionDirection, name="InsertionController"))
    rootNode.addObject(SmoothInsertionController(rootNode=rootNode, incrementDistance=incrementDistance,
                                                 incrementAngle=incrementAngle,
                                                 name="InsertionController"))

    # rootNode.addObject(BendingController(rootNode=rootNode, cosseratNode=rateAngularDeformNode,
    #                                      frameNode=mappedFrameNode, bendingMoment=beamBendingMoment,
    #                                      fixedIndices=fixedSections, momentAxis='y', name="BendingController"))
