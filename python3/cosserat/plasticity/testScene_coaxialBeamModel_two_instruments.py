# -*- coding: utf-8 -*-
"""
Created on April 22 2022

@author: ckrewcun
"""

import Sofa
import os
import numpy as np
from CosseratNavigationController import CombinedInstrumentsController

# constants
GRAVITY = 0.0 # 9810
TOT_MASS = 0.1
DENSITY = 0.02
DT = 1e-2

# -------------------------------------#
# -----        SOFA scene        ----- #
# -------------------------------------#

pluginNameList = 'SofaPython3 CosseratPlugin' \
                 ' Sofa.Component.MechanicalLoad' \
                 ' Sofa.Component.ODESolver.Backward' \
                 ' Sofa.Component.SolidMechanics.Spring' \
                 ' Sofa.Component.StateContainer' \
                 ' Sofa.Component.Visual ' \
                 ' Sofa.Component.Constraint.Projective ' \
                 ' Sofa.Component.LinearSolver.Direct'
                 # ' Sofa.Component.LinearSolver.Iterative'

visualFlagList = 'showVisualModels showBehaviorModels showCollisionModels' \
                 ' hideBoundingCollisionModels hideForceFields' \
                 ' hideInteractionForceFields hideWireframe' \
                 ' showMechanicalMappings'

def createScene(rootNode):

    rootNode.addObject('RequiredPlugin', pluginName=pluginNameList, printLog='0')
    rootNode.addObject('VisualStyle', displayFlags=visualFlagList)

    rootNode.addObject('DefaultVisualManagerLoop')
    rootNode.findData('dt').value = DT
    rootNode.findData('gravity').value = [0., 0., -GRAVITY]

    rootNode.addObject('DefaultAnimationLoop')

    # -------------------------------------------------------------------- #
    # -----                      Beam parameters                     ----- #
    # -------------------------------------------------------------------- #

    # --- Instrument0 --- #

    # Define: the total length, number of beams, and number of frames
    totalLength0 = 15

    nbBeams0 = 5
    oneBeamLength = totalLength0 / nbBeams0

    # nbFramesMax = 14
    # distBetweenFrames = totalLength0 / nbFrames
    nbFrames0 = 14

    # --- Instrument1 --- #

    # Define: the total length, number of beams, and number of frames
    totalLength1 = 20

    nbBeams1 = 5
    oneBeamLength = totalLength1 / nbBeams1

    # distBetweenFrames = totalLength1 / nbFrames
    nbFrames1 = 17

    # --- Common --- #
    nbMaxInstrumentBeams = max(nbBeams0, nbBeams1)

    # -------------------------------------------------------------------- #
    # -----                   First beam components                  ----- #
    # -------------------------------------------------------------------- #

    # Redefining the number of beams 'in stock' in order to handle all
    # discretisation scenarios
    nbBeams0PlusStock = nbBeams0 + nbMaxInstrumentBeams

    # ----- Rigid base ----- #

    instrument0Node = rootNode.addChild('Instrument0')
    instrument0Node.addObject('EulerImplicitSolver', rayleighStiffness="0.", rayleighMass='0.')
    # instrument0Node.addObject('CGLinearSolver', name="solver",
    #                    iterations='100', tolerance='1e-5', threshold='1e-5')
    instrument0Node.addObject('SparseLUSolver',
                              template='CompressedRowSparseMatrixd',
                              printLog="false")

    rigidBaseNode0 = instrument0Node.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode0.addObject('MechanicalObject', template='Rigid3d',
                                           name="RigidBaseMO", position=[0., 0., 0., 0, 0, 0, 1],
                                           showObject=True,
                                           showObjectScale=2.)
    rigidBaseNode0.addObject('RestShapeSpringsForceField', name='spring',
                             stiffness="5.e8", angularStiffness="5.e8",
                             external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")

    # ----- Rate of angular deformation ----- #
    # Define the length of each beam in a list, the positions of each beam

    beamStrainDoFs = []
    beamLengths = []
    sum = 0.
    beamCurvAbscissa = []
    beamCurvAbscissa.append(0.0)

    # EXPERIMENTAL: insertion navigation
    # Initially defining beams with 0 length, at 0
    # For each beam, we:
    #     - Add a Vec3 in beamStrainDoFs, with each strain component set at 0.
    #     - Add 0 in the length data field for the Cosserat mapping
    #     - Add the last curvilinear abscissa (total length) in the beams
    # curvilinear abcissas in the Cosserat mapping.
    # This is equivalent to saying that the beams are 0-length,
    # 0-strain beams, added at the beginning of the instrument.
    for i in range(0, nbBeams0PlusStock):
        beamStrainDoFs.append([0, 0, 0])
        beamLengths.append(0)
        beamCurvAbscissa.append(0)

    # for i in range(nbBeams0PlusStock):
    #     beamStrainDoFs.append([0, 0, 0])
    #     beamLengths.append(oneBeamLength)
    #     sum += beamLengths[i+nbStockBeams]
    #     beamCurvAbscissa.append(sum)
    # beamCurvAbscissa[nbBeams0PlusStock+nbStockBeams] = totalLength0

    # Define angular rate which is the torsion(x) and bending (y, z) of each section
    rateAngularDeformNode0 = instrument0Node.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode0.addObject('MechanicalObject',
                                                           template='Vec3d',
                                                           name='rateAngularDeformMO',
                                                           position=beamStrainDoFs,
                                                           showIndices=0,
                                                           rest_position=[0.0, 0.0, 0.0])

    beamCrossSectionShape='circular'
    sectionRadius = 0.5
    poissonRatio = 0.42
    beamPoissonRatioList = [poissonRatio]*(nbBeams0PlusStock)
    youngModulus = 5.0e6
    beamYoungModulusList = [youngModulus]*(nbBeams0PlusStock)
    yieldStress = 5.0e4
    yieldStressList = [yieldStress]*(nbBeams0PlusStock)
    plasticModulus = 2.0e5
    plasticModulusList = [plasticModulus]*(nbBeams0PlusStock)
    hardeningCoeff = 0.5
    hardeningCoefficientList = [hardeningCoeff]*(nbBeams0PlusStock)
    ### Plastic FF version
    rateAngularDeformNode0.addObject('BeamPlasticLawForceField', name="beamForceField",
                                     crossSectionShape=beamCrossSectionShape,
                                     radius=sectionRadius, variantSections="true",
                                     length=beamLengths, poissonRatioList=beamPoissonRatioList,
                                     youngModulusList=beamYoungModulusList,
                                     initialYieldStresses=yieldStressList,
                                     plasticModuli=plasticModulusList,
                                     mixedHardeningCoefficients= hardeningCoefficientList)
    ### Elastic FF version
    # rateAngularDeformNode.addObject('BeamHookeLawForceField', name="beamForceField",
    #                                 crossSectionShape=beamCrossSectionShape,
    #                                 radius=sectionRadius, variantSections="true",
    #                                 length=beamLengths, poissonRatioList=beamPoissonRatioList,
    #                                 youngModulusList=beamYoungModulusList)

    beamBendingMoment = 1.0e5
    bendingForces = np.array([0, beamBendingMoment, beamBendingMoment])
    # momentIndices = range(1, nbBeams0PlusStock)
    momentIndices = [nbBeams0PlusStock-1]
    # rateAngularDeformNode.addObject('ConstantForceField', name='Moment',
    #                                 indices=momentIndices,
    #                                 forces=bendingForces)

    # EXPERIMENTAL: navigation simulation
    # Adding constraints on the additional beams which are not meant to be
    # simulated at the beginning
    # fixedIndices = list(range(nbBeams0PlusStock, nbBeams0PlusStock+nbStockBeams))
    fixedIndices = list(range(0, nbBeams0PlusStock))
    rateAngularDeformNode0.addObject('FixedConstraint', name='FixedConstraintOnStock',
                                    indices=fixedIndices)

    # ----- Frames ----- #

    # Define local frames related to each section and parameters framesCurvAbscissa
    frames6DDoFs = []
    frames3DDoFs = []
    framesCurvAbscissa = []

    # EXPERIMENTAL: navigation simulation
    # At the beginning of the simulation, when no beam element is actually
    # simulated, all rigid frames are set at 0
    for i in range(nbFrames0):
        frames3DDoFs.append([0, 0, 0])
        frames6DDoFs.append([0, 0, 0,  0, 0, 0, 1])
        framesCurvAbscissa.append(0)

    frames6DDoFs.append([0, 0, 0, 0, 0, 0, 1])
    framesCurvAbscissa.append(0)

    # The node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode0.addChild('MappedFrames')
    rateAngularDeformNode0.addChild(mappedFrameNode)

    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                         name="FramesMO", position=frames6DDoFs,
                                         showObject=True, showObjectScale=1)

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=beamCurvAbscissa,
                              curv_abs_output=framesCurvAbscissa, input1=inputMO, input2=inputMO_rigid,
                              output=outputMO, forcefield='@../../rateAngularDeform/beamForceField',
                              drawBeamSegments=True, nonColored=False, debug=0, printLog=False)



    # -------------------------------------------------------------------- #
    # -----                  Second beam components                  ----- #
    # -------------------------------------------------------------------- #

    # Redefining the number of beams 'in stock' in order to handle all
    # discretisation scenarios
    nbBeams1PlusStock = nbBeams1 + nbMaxInstrumentBeams

    # ----- Rigid base ----- #

    instrument1Node = rootNode.addChild('Instrument1')
    instrument1Node.addObject('EulerImplicitSolver', rayleighStiffness="0.", rayleighMass='0.')
    # instrument1Node.addObject('CGLinearSolver', name="solver",
    #                    iterations='100', tolerance='1e-5', threshold='1e-5')
    instrument1Node.addObject('SparseLUSolver',
                              template='CompressedRowSparseMatrixd',
                              printLog="false")

    rigidBaseNode1 = instrument1Node.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode1.addObject('MechanicalObject', template='Rigid3d',
                                           name="RigidBaseMO", position=[0., 0., 3., 0, 0, 0, 1],
                                           showObject=True,
                                           showObjectScale=2.)
    rigidBaseNode1.addObject('RestShapeSpringsForceField', name='spring',
                             stiffness="5.e8", angularStiffness="5.e8",
                             external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")

    # ----- Rate of angular deformation ----- #
    # Define the length of each beam in a list, the positions of each beam

    beamStrainDoFs = []
    beamLengths = []
    sum = 0.
    beamCurvAbscissa = []
    beamCurvAbscissa.append(0.0)

    # EXPERIMENTAL: insertion navigation
    # Initially defining beams with 0 length, at 0
    # For each beam, we:
    #     - Add a Vec3 in beamStrainDoFs, with each strain component set at 0.
    #     - Add 0 in the length data field for the Cosserat mapping
    #     - Add the last curvilinear abscissa (total length) in the beams
    # curvilinear abcissas in the Cosserat mapping.
    # This is equivalent to saying that the beams are 0-length,
    # 0-strain beams, added at the beginning of the instrument.
    for i in range(0, nbBeams1PlusStock):
        beamStrainDoFs.append([0, 0, 0])
        beamLengths.append(0)
        beamCurvAbscissa.append(0)

    # for i in range(nbBeams1PlusStock):
    #     beamStrainDoFs.append([0, 0, 0])
    #     beamLengths.append(oneBeamLength)
    #     sum += beamLengths[i+nbStockBeams]
    #     beamCurvAbscissa.append(sum)
    # beamCurvAbscissa[nbBeams1PlusStock+nbStockBeams] = totalLength1

    # Define angular rate which is the torsion(x) and bending (y, z) of each section
    rateAngularDeformNode0 = instrument1Node.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode0.addObject('MechanicalObject',
                                                           template='Vec3d',
                                                           name='rateAngularDeformMO',
                                                           position=beamStrainDoFs,
                                                           showIndices=0,
                                                           rest_position=[0.0, 0.0, 0.0])

    beamCrossSectionShape='circular'
    sectionRadius = 0.4
    poissonRatio = 0.45
    beamPoissonRatioList = [poissonRatio]*(nbBeams1PlusStock)
    youngModulus = 5.0e6
    beamYoungModulusList = [youngModulus]*(nbBeams1PlusStock)
    yieldStress = 5.0e4
    yieldStressList = [yieldStress]*(nbBeams1PlusStock)
    plasticModulus = 2.0e5
    plasticModulusList = [plasticModulus]*(nbBeams1PlusStock)
    hardeningCoeff = 0.5
    hardeningCoefficientList = [hardeningCoeff]*(nbBeams1PlusStock)
    ### Plastic FF version
    rateAngularDeformNode0.addObject('BeamPlasticLawForceField', name="beamForceField",
                                     crossSectionShape=beamCrossSectionShape,
                                     radius=sectionRadius, variantSections="true",
                                     length=beamLengths, poissonRatioList=beamPoissonRatioList,
                                     youngModulusList=beamYoungModulusList,
                                     initialYieldStresses=yieldStressList,
                                     plasticModuli=plasticModulusList,
                                     mixedHardeningCoefficients= hardeningCoefficientList)
    ### Elastic FF version
    # rateAngularDeformNode.addObject('BeamHookeLawForceField', name="beamForceField",
    #                                 crossSectionShape=beamCrossSectionShape,
    #                                 radius=sectionRadius, variantSections="true",
    #                                 length=beamLengths, poissonRatioList=beamPoissonRatioList,
    #                                 youngModulusList=beamYoungModulusList)

    beamBendingMoment = 1.0e5
    bendingForces = np.array([0, beamBendingMoment, beamBendingMoment])
    # momentIndices = range(1, nbBeams1PlusStock)
    momentIndices = [nbBeams1PlusStock-1]
    # rateAngularDeformNode.addObject('ConstantForceField', name='Moment',
    #                                 indices=momentIndices,
    #                                 forces=bendingForces)

    # EXPERIMENTAL: navigation simulation
    # Adding constraints on the additional beams which are not meant to be
    # simulated at the beginning
    # fixedIndices = list(range(nbBeams1PlusStock, nbBeams1PlusStock+nbStockBeams))
    fixedIndices = list(range(0, nbBeams1PlusStock))
    rateAngularDeformNode0.addObject('FixedConstraint', name='FixedConstraintOnStock',
                                    indices=fixedIndices)

    # rateAngularDeformNode0.addObject('FixedConstraint', name='FixedConstraint',
    #                                 indices=fixedIndices)

    # ----- Frames ----- #

    # Define local frames related to each section and parameters framesCurvAbscissa
    frames6DDoFs = []
    frames3DDoFs = []
    framesCurvAbscissa = []

    # EXPERIMENTAL: navigation simulation
    # At the beginning of the simulation, when no beam element is actually
    # simulated, all rigid frames are set at 0
    for i in range(nbFrames1):
        frames3DDoFs.append([0, 0, 3])
        frames6DDoFs.append([0, 0, 3,  0, 0, 0, 1])
        framesCurvAbscissa.append(0)

    frames6DDoFs.append([0, 0, 3, 0, 0, 0, 1])
    framesCurvAbscissa.append(0)

    # The node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode1.addChild('MappedFrames')
    rateAngularDeformNode0.addChild(mappedFrameNode)

    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                         name="FramesMO", position=frames6DDoFs,
                                         showObject=True, showObjectScale=1)

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=beamCurvAbscissa,
                              curv_abs_output=framesCurvAbscissa, input1=inputMO, input2=inputMO_rigid,
                              output=outputMO, forcefield='@../../rateAngularDeform/beamForceField',
                              drawBeamSegments=True, nonColored=False, debug=0, printLog=False)




    # -------------------------------------------------------------------- #
    # -----                    Python controllers                    ----- #
    # -------------------------------------------------------------------- #

    # rootNode.addObject(BendingController(name="BendingController",
    #                                      bendingMoment=beamBendingMoment))
    nbInstruments=2
    instrumentBeamNumbers=[nbBeams0, nbBeams1]
    instrumentFrameNumbers=[nbFrames0, nbFrames1]
    incrementDistance=0.1
    incrementDirection = np.array([1., 0., 0.])
    isInstrumentStraightVect=[True, True]
    curvAbsTolerance= 1.0e-4
    instrumentLengths=[totalLength0, totalLength1]

    rootNode.addObject(CombinedInstrumentsController(
                            name="NavigationController",
                            rootNode=rootNode,
                            nbInstruments=nbInstruments,
                            instrumentBeamNumberVect=instrumentBeamNumbers,
                            instrumentFrameNumberVect=instrumentFrameNumbers,
                            incrementDistance=incrementDistance,
                            incrementDirection=incrementDirection,
                            isInstrumentStraightVect=isInstrumentStraightVect,
                            curvAbsTolerance=curvAbsTolerance,
                            instrumentLengths=instrumentLengths))

    return rootNode
