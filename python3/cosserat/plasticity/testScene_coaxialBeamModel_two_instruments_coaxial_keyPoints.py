# -*- coding: utf-8 -*-
"""
Created on January 20 2023

@author: ckrewcun
"""

import Sofa
import os
import numpy as np
from CosseratNavigationController import CombinedInstrumentsController
from instrument import Instrument # defining a class for instrument characteristics

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
                 ' Sofa.Component.LinearSolver.Direct' \
                 ' Sofa.Component.AnimationLoop' \
                 ' Sofa.Component.Constraint.Lagrangian.Correction' \
                 ' Sofa.Component.Constraint.Lagrangian.Solver' \
                 ' Sofa.Component.Mass'

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

    rootNode.addObject('FreeMotionAnimationLoop')

    # --- Constraint handling --- #
    rootNode.addObject('GenericConstraintSolver', tolerance=1e-8,
                       maxIterations=2e3, printLog=False)


    # --- Collisions --- #
    # rootNode.addObject('DefaultPipeline', verbose="0")
    # rootNode.addObject('BruteForceBroadPhase', name="BroadPhase")
    # rootNode.addObject('BVHNarrowPhase', name="NarrowPhase")
    # rootNode.addObject('DefaultContactManager', response="FrictionContactConstraint",
    #                    responseParams=frictionCoefficientParam)
    # rootNode.addObject('LocalMinDistance', name="Proximity",
    #                    alarmDistance=100.*sectionRadius,
    #                    contactDistance=10.*sectionRadius, angleCone=0.1)


    # -------------------------------------------------------------------- #
    # -----                      Beam parameters                     ----- #
    # -------------------------------------------------------------------- #

    # --- Instrument0 --- #

    # Define: the total length, number of beams, and number of frames
    totalLength0 = 25

    nbBeams0_0 = 3
    nbBeams0_1 = 2
    nbBeams0 = nbBeams0_0 + nbBeams0_1

    # nbFramesMax = 14
    # distBetweenFrames = totalLength0 / nbFrames
    nbFrames0 = 30

    # --- Instrument1 --- #

    # Define: the total length, number of beams, and number of frames
    totalLength1 = 30

    nbBeams1_0 = 4
    nbBeams1_1 = 1
    nbBeams1 = nbBeams1_0 + nbBeams1_1

    # distBetweenFrames = totalLength1 / nbFrames
    nbFrames1 = 35

    # --- Common --- #
    nbMaxInstrumentBeams = max(nbBeams0, nbBeams1)
    totalMass = 0.022



    # -------------------------------------#
    # -----      Control points      ----- #
    # -------------------------------------#

    controlPointPos = np.array([0., 0., 0., 0., 0., 0., 1.])
    # Instrument 0
    controlPointNode0 = rootNode.addChild('controlPointNode0')
    controlPointNode0.addObject('MechanicalObject', template='Rigid3d', name="controlPointMO",
                               position=controlPointPos, showObject='1', showObjectScale='0.5')

    # Instrument1
    controlPointNode1 = rootNode.addChild('controlPointNode1')
    controlPointNode1.addObject('MechanicalObject', template='Rigid3d', name="controlPointMO",
                               position=controlPointPos, showObject='1', showObjectScale='0.5')


    # -------------------------------------------------------------------- #
    # -----                        Solver node                       ----- #
    # -------------------------------------------------------------------- #

    solverNode = rootNode.addChild('solverNode')

    solverNode.addObject('EulerImplicitSolver', rayleighStiffness="0.", rayleighMass='0.')
    # solverNode.addObject('CGLinearSolver', name="solver",
    #                      iterations='2000', tolerance='1e-8', threshold='1e-12')
    solverNode.addObject('SparseLUSolver',
                         template='CompressedRowSparseMatrixd',
                         printLog="false")
    solverNode.addObject('GenericConstraintCorrection', printLog=False,
                         linearSolver="@SparseLUSolver",
                         ODESolver="@EulerImplicitSolver")

    # -------------------------------------------------------------------- #
    # -----                   First beam components                  ----- #
    # -------------------------------------------------------------------- #

    # Redefining the number of beams 'in stock' in order to handle all
    # discretisation scenarios
    nbBeams0PlusStock = nbBeams0 + nbMaxInstrumentBeams

    # ----- Rigid base ----- #

    instrument0Node = solverNode.addChild('Instrument0')

    rigidBaseNode0 = instrument0Node.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode0.addObject('MechanicalObject', template='Rigid3d',
                                           name="RigidBaseMO", position=[0., 0., 0., 0, 0, 0, 1],
                                           showObject=False,
                                           showObjectScale=2.)
    # rigidBaseNode0.addObject('RestShapeSpringsForceField', name='spring',
    #                          stiffness="5.e8", angularStiffness="5.e8",
    #                          external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")
    rigidBaseNode0.addObject('RestShapeSpringsForceField', name='controlSpring',
                             stiffness="5.e8", angularStiffness="1.0e5",
                             external_rest_shape="@../../../controlPointNode0/controlPointMO",
                             external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")


    # ----- Rate of angular deformation ----- #
    # Define the length of each beam in a list, the positions of each beam

    beamStrainDoFs0 = []
    beamLengths = []
    sum = 0.
    beamCurvAbscissa = []
    beamCurvAbscissa.append(0.0)

    instrument0ConstantRestStrain0 = [0., 0., 0.]
    instrument0ConstantRestStrain1 = [0., 0., 0.05]

    # EXPERIMENTAL: insertion navigation
    # Initially defining beams with 0 length, at 0
    # For each beam, we:
    #     - Add a Vec3 in beamStrainDoFs0, with each strain component set at 0.
    #     - Add 0 in the length data field for the Cosserat mapping
    #     - Add the last curvilinear abscissa (total length) in the beams
    # curvilinear abcissas in the Cosserat mapping.
    # This is equivalent to saying that the beams are 0-length,
    # 0-strain beams, added at the beginning of the instrument.
    for i in range(0, nbBeams0PlusStock-nbBeams0_1):
        beamStrainDoFs0.append(instrument0ConstantRestStrain0)
        beamLengths.append(0)
        beamCurvAbscissa.append(0)

    for i in range(nbBeams0PlusStock-nbBeams0_1, nbBeams0PlusStock):
        beamStrainDoFs0.append(instrument0ConstantRestStrain1)
        beamLengths.append(0)
        beamCurvAbscissa.append(0)

    # Define angular rate which is the torsion(x) and bending (y, z) of each section
    rateAngularDeformNode0 = instrument0Node.addChild('rateAngularDeform')
    rateAngularDeformMO0 = rateAngularDeformNode0.addObject('MechanicalObject',
                                                            template='Vec3d',
                                                            name='rateAngularDeformMO',
                                                            position=beamStrainDoFs0,
                                                            showIndices=0,
                                                            rest_position=beamStrainDoFs0)

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
    # rateAngularDeformNode0.addObject('BeamPlasticLawForceField', name="beamForceField",
    #                                  crossSectionShape=beamCrossSectionShape,
    #                                  radius=sectionRadius, variantSections="true",
    #                                  length=beamLengths, poissonRatioList=beamPoissonRatioList,
    #                                  youngModulusList=beamYoungModulusList,
    #                                  initialYieldStresses=yieldStressList,
    #                                  plasticModuli=plasticModulusList,
    #                                  mixedHardeningCoefficients= hardeningCoefficientList)
    ### Elastic FF version
    rateAngularDeformNode0.addObject('BeamHookeLawForceField', name="beamForceField",
                                     crossSectionShape=beamCrossSectionShape,
                                     radius=sectionRadius, variantSections="true",
                                     length=beamLengths, poissonRatioList=beamPoissonRatioList,
                                     youngModulusList=beamYoungModulusList, innerRadius=0.4)

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
                                         showObject=False, showObjectScale=1)

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO0.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=beamCurvAbscissa,
                              curv_abs_output=framesCurvAbscissa, input1=inputMO, input2=inputMO_rigid,
                              output=outputMO, forcefield='@../../rateAngularDeform/beamForceField',
                              drawBeamSegments=True, nonColored=False, debug=0, printLog=False)

    # Second node of mapped frames, to apply 'constraints' on the coaxial beam segments
    coaxialFrameNode0 = rigidBaseNode0.addChild('coaxialSegmentFrames')
    rateAngularDeformNode0.addChild(coaxialFrameNode0)

    # By default, the nbIntermediateConstraintFrames parameter of the
    # CosseratNavigationController is 0, meaning that we don't add intermediate
    # coaxial frames on each coaxial beam segments (between the beam extremities)
    # to further constrain the segments. In such scenario, we need at most
    # nbBeams0PlusStock+1 coaxial frames during the simulation.
    # /!\ Here, we will pass the nbIntermediateConstraintFrames parameter with
    # a value of 1, meaning that we need at most 2*nbBeams0PlusStock + 1 coaxial
    # frames
    nbIntermediateConstraintFrames = 1
    nbCoaxialFrames0 = (nbIntermediateConstraintFrames+1)*nbBeams0PlusStock + 1
    coaxialFrames0InitPos = [[0., 0., 0., 0., 0., 0., 1.]]*nbCoaxialFrames0
    coaxialFrame0CurvAbscissa = [0.]*nbCoaxialFrames0

    coaxialFramesMO0 = coaxialFrameNode0.addObject('MechanicalObject', template='Rigid3d',
                                                   name="coaxialFramesMO",
                                                   position=coaxialFrames0InitPos,
                                                   showObject=True, showObjectScale=1)

    coaxialFrameNode0.addObject('UniformMass', totalMass=totalMass,
                                name="UniformMass0", showAxisSizeFactor='0.')

    coaxialFrameNode0.addObject('DiscreteCosseratMapping',
                               name='CoaxialCosseratMapping',
                               curv_abs_input=beamCurvAbscissa,
                               curv_abs_output=coaxialFrame0CurvAbscissa,
                               input1=rateAngularDeformMO0.getLinkPath(),
                               input2=RigidBaseMO.getLinkPath(),
                               output=coaxialFramesMO0.getLinkPath(),
                               forcefield='@../../rateAngularDeform/beamForceField',
                               drawBeamSegments=False, nonColored=False,
                               debug=0, printLog=False)



    # -------------------------------------------------------------------- #
    # -----                  Second beam components                  ----- #
    # -------------------------------------------------------------------- #

    # Redefining the number of beams 'in stock' in order to handle all
    # discretisation scenarios
    nbBeams1PlusStock = nbBeams1 + nbMaxInstrumentBeams

    # ----- Rigid base ----- #

    instrument1Node = solverNode.addChild('Instrument1')

    rigidBaseNode1 = instrument1Node.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode1.addObject('MechanicalObject', template='Rigid3d',
                                           name="RigidBaseMO", position=[0., 0., 0., 0, 0, 0, 1],
                                           showObject=False,
                                           showObjectScale=2.)
    # rigidBaseNode1.addObject('RestShapeSpringsForceField', name='spring',
    #                          stiffness="5.e8", angularStiffness="5.e8",
    #                          external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")
    rigidBaseNode1.addObject('RestShapeSpringsForceField', name='controlSpring',
                             stiffness="5.e8", angularStiffness="1.0e5",
                             external_rest_shape="@../../../controlPointNode1/controlPointMO",
                             external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")

    # ----- Rate of angular deformation ----- #
    # Define the length of each beam in a list, the positions of each beam

    beamStrainDoFs1 = []
    beamLengths = []
    sum = 0.
    beamCurvAbscissa = []
    beamCurvAbscissa.append(0.0)

    # instrument1ConstantRestStrain = [0., 0.1, 0.]
    instrument1ConstantRestStrain0 = [0., 0., 0.]
    instrument1ConstantRestStrain1 = [0., 0.1, 0.]

    # EXPERIMENTAL: insertion navigation
    # Initially defining beams with 0 length, at 0
    # For each beam, we:
    #     - Add a Vec3 in beamStrainDoFs1, with each strain component set at 0.
    #     - Add 0 in the length data field for the Cosserat mapping
    #     - Add the last curvilinear abscissa (total length) in the beams
    # curvilinear abcissas in the Cosserat mapping.
    # This is equivalent to saying that the beams are 0-length,
    # 0-strain beams, added at the beginning of the instrument.
    for i in range(0, nbBeams1PlusStock-nbBeams1_1):
        beamStrainDoFs1.append(instrument1ConstantRestStrain0)
        beamLengths.append(0)
        beamCurvAbscissa.append(0)

    for i in range(nbBeams1PlusStock-nbBeams1_1, nbBeams1PlusStock):
        beamStrainDoFs1.append(instrument1ConstantRestStrain1)
        beamLengths.append(0)
        beamCurvAbscissa.append(0)

    # Define angular rate which is the torsion(x) and bending (y, z) of each section
    rateAngularDeformNode1 = instrument1Node.addChild('rateAngularDeform')
    rateAngularDeformMO1 = rateAngularDeformNode1.addObject('MechanicalObject',
                                                            template='Vec3d',
                                                            name='rateAngularDeformMO',
                                                            position=beamStrainDoFs1,
                                                            showIndices=0,
                                                            rest_position=beamStrainDoFs1)

    beamCrossSectionShape='circular'
    sectionRadius = 0.3
    poissonRatio = 0.45
    beamPoissonRatioList = [poissonRatio]*(nbBeams1PlusStock)
    youngModulus = 5.0e10
    beamYoungModulusList = [youngModulus]*(nbBeams1PlusStock)
    yieldStress = 5.0e4
    yieldStressList = [yieldStress]*(nbBeams1PlusStock)
    plasticModulus = 2.0e5
    plasticModulusList = [plasticModulus]*(nbBeams1PlusStock)
    hardeningCoeff = 0.5
    hardeningCoefficientList = [hardeningCoeff]*(nbBeams1PlusStock)
    ### Plastic FF version
    # rateAngularDeformNode1.addObject('BeamPlasticLawForceField', name="beamForceField",
    #                                  crossSectionShape=beamCrossSectionShape,
    #                                  radius=sectionRadius, variantSections="true",
    #                                  length=beamLengths, poissonRatioList=beamPoissonRatioList,
    #                                  youngModulusList=beamYoungModulusList,
    #                                  initialYieldStresses=yieldStressList,
    #                                  plasticModuli=plasticModulusList,
    #                                  mixedHardeningCoefficients= hardeningCoefficientList)
    ### Elastic FF version
    rateAngularDeformNode1.addObject('BeamHookeLawForceField', name="beamForceField",
                                     crossSectionShape=beamCrossSectionShape,
                                     radius=sectionRadius, variantSections="true",
                                     length=beamLengths, poissonRatioList=beamPoissonRatioList,
                                     youngModulusList=beamYoungModulusList)

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
    rateAngularDeformNode1.addObject('FixedConstraint', name='FixedConstraintOnStock',
                                     indices=fixedIndices)

    # At the beginning of the simulation, no instrument is deployed
    constraintSpringStiffness = 5.0e8
    constraintSpringDamping = 0.


    # ----- Frames ----- #

    # Define local frames related to each section and parameters framesCurvAbscissa
    frames6DDoFs = []
    frames3DDoFs = []
    framesCurvAbscissa = []

    # EXPERIMENTAL: navigation simulation
    # At the beginning of the simulation, when no beam element is actually
    # simulated, all rigid frames are set at 0
    for i in range(nbFrames1):
        frames3DDoFs.append([0, 0, 0])
        frames6DDoFs.append([0, 0, 0, 0, 0, 0, 1])
        framesCurvAbscissa.append(0)

    frames6DDoFs.append([0, 0, 0, 0, 0, 0, 1])
    framesCurvAbscissa.append(0)

    # The node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode1.addChild('MappedFrames')
    rateAngularDeformNode1.addChild(mappedFrameNode)

    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                         name="FramesMO", position=frames6DDoFs,
                                         showObject=False, showObjectScale=1)

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO1.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=beamCurvAbscissa,
                              curv_abs_output=framesCurvAbscissa, input1=inputMO, input2=inputMO_rigid,
                              output=outputMO, forcefield='@../../rateAngularDeform/beamForceField',
                              drawBeamSegments=True, nonColored=False, debug=0, printLog=False)

    # Second node of mapped frames, to apply 'constraints' on the coaxial beam segments
    coaxialFrameNode1 = rigidBaseNode1.addChild('coaxialSegmentFrames')
    rateAngularDeformNode1.addChild(coaxialFrameNode1)

    # Cf comment on nbCoaxialFrames0 definition for details on the number of coaxial
    # frames
    nbCoaxialFrames1 = (nbIntermediateConstraintFrames+1)*nbBeams1PlusStock + 1
    coaxialFrames1InitPos = [[0., 0., 0., 0., 0., 0., 1.]]*nbCoaxialFrames1
    coaxialFrame1CurvAbscissa = [0.]*nbCoaxialFrames1

    coaxialFramesMO1 = coaxialFrameNode1.addObject('MechanicalObject', template='Rigid3d',
                                                  name="coaxialFramesMO",
                                                  position=coaxialFrames1InitPos,
                                                  showObject=True, showObjectScale=1)

    coaxialFrameNode1.addObject('UniformMass', totalMass=totalMass,
                                name="UniformMass1", showAxisSizeFactor='0.')

    coaxialFrameNode1.addObject('DiscreteCosseratMapping',
                               name='CoaxialCosseratMapping',
                               curv_abs_input=beamCurvAbscissa,
                               curv_abs_output=coaxialFrame1CurvAbscissa,
                               input1=rateAngularDeformMO1.getLinkPath(),
                               input2=RigidBaseMO.getLinkPath(),
                               output=coaxialFramesMO1.getLinkPath(),
                               forcefield='@../../rateAngularDeform/beamForceField',
                               drawBeamSegments=False, nonColored=False,
                               debug=0, printLog=False)

    ## Difference mapping node

    constraintWith0Node = coaxialFrameNode1.addChild('constraint1With0')
    coaxialFrameNode0.addChild(constraintWith0Node)

    rigidDiffMO = constraintWith0Node.addObject('MechanicalObject', template='Rigid3d',
                                                 name="rigidDiffMO", position=[0., 0., 0., 0., 0., 0., 1.],
                                                 showObject=False, showObjectScale=1)

    constraintWith0Node.addObject('RigidDistanceMapping', name='coaxialFramesDistanceMapping',
                                  input1=coaxialFramesMO1.getLinkPath(),
                                  input2=coaxialFramesMO0.getLinkPath(),
                                  output=rigidDiffMO.getLinkPath(),
                                  first_point=[], second_point=[])

    constraintWith0Node.addObject('CosseratNeedleSlidingConstraint',
                                  name='constraintMappingConstraint',
                                  template="Rigid3d",
                                  useDirections=np.array([0, 1, 1, 0, 0, 0]))


    # -------------------------------------------------------------------- #
    # -----                    Python controllers                    ----- #
    # -------------------------------------------------------------------- #

    nbInstruments=2
    instrument0KeyPoints = [15]
    instrument1KeyPoints = [25]
    nbBeamDistribution0 = [nbBeams0_0, nbBeams0_1]
    nbBeamDistribution1 = [nbBeams1_0, nbBeams1_1]
    instrumentFrameNumbers=[nbFrames0, nbFrames1]
    incrementDistance=0.1
    incrementAngle=5.0
    incrementDirection = np.array([1., 0., 0.])
    curvAbsTolerance= 1.0e-4

    instrument0 = Instrument(instrumentNode=instrument0Node,
                             totalLength=totalLength0,
                             keyPoints=instrument0KeyPoints,
                             nbBeamDistribution=nbBeamDistribution0,
                             restStrain=[instrument0ConstantRestStrain0, instrument0ConstantRestStrain1],
                             curvAbsTolerance=curvAbsTolerance)
    instrument1 = Instrument(instrumentNode=instrument1Node,
                             totalLength=totalLength1,
                             keyPoints=instrument1KeyPoints,
                             nbBeamDistribution=nbBeamDistribution1,
                             restStrain=[instrument1ConstantRestStrain0, instrument1ConstantRestStrain1],
                             curvAbsTolerance=curvAbsTolerance)
    instrumentList = [instrument0, instrument1]

    rootNode.addObject(CombinedInstrumentsController(
                            name="NavigationController",
                            rootNode=rootNode,
                            solverNode=solverNode,
                            nbInstruments=nbInstruments,
                            instrumentFrameNumberVect=instrumentFrameNumbers,
                            incrementDistance=incrementDistance,
                            incrementAngle=incrementAngle,
                            incrementDirection=incrementDirection,
                            instrumentList=instrumentList,
                            curvAbsTolerance=curvAbsTolerance,
                            nbIntermediateConstraintFrames=nbIntermediateConstraintFrames))

    return rootNode
