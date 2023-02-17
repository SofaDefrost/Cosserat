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
                 ' Sofa.Component.Mass' \
                 ' Sofa.Component.Collision.Geometry' \
                 ' Sofa.Component.Mapping.Linear' \
                 ' Sofa.Component.Topology.Container.Dynamic' \
                 ' Sofa.Component.Collision.Detection.Algorithm' \
                 ' Sofa.Component.Collision.Detection.Intersection' \
                 ' Sofa.Component.Collision.Response.Contact' \
                 ' Sofa.Component.Engine.Transform' \
                 ' Sofa.Component.IO.Mesh' \
                 ' Sofa.Component.Topology.Container.Constant'

visualFlagList = 'showVisualModels showBehaviorModels showCollisionModels' \
                 ' hideBoundingCollisionModels hideForceFields' \
                 ' showInteractionForceFields hideWireframe' \
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

    # --- Collision handling --- #

    frictionCoefficient = 0.1
    frictionCoefficientParam = "mu="+str(frictionCoefficient)
    largestSectionRadius = (6.5e-1/6.0)
    alarmDistance = 10.0*largestSectionRadius
    contactDistance = 1.0*largestSectionRadius


    rootNode.addObject('DefaultPipeline', verbose="0")
    rootNode.addObject('BruteForceBroadPhase', name="BroadPhase")
    rootNode.addObject('BVHNarrowPhase', name="NarrowPhase")
    rootNode.addObject('DefaultContactManager', response="FrictionContactConstraint",
                       responseParams=frictionCoefficientParam)
    rootNode.addObject('LocalMinDistance', name="Proximity",
                       alarmDistance=alarmDistance,
                       contactDistance=contactDistance, angleCone=0.1)


    # -------------------------------------------------------------------- #
    # -----                      Beam parameters                     ----- #
    # -------------------------------------------------------------------- #

    # --- Instrument0 --- #

    # Define: the total length, number of beams, and number of frames
    totalLength0 = 50

    nbBeams0_0 = 60
    nbBeams0_1 = 10
    nbBeams0 = nbBeams0_0 + nbBeams0_1

    nbFrames0 = 100

    # --- Instrument1 --- #

    # Define: the total length, number of beams, and number of frames
    totalLength1 = 50

    nbBeams1_0 = 60
    nbBeams1_1 = 15
    nbBeams1 = nbBeams1_0 + nbBeams1_1

    nbFrames1 = 100

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
    instrument0ConstantRestStrain1 = [0., 0., 0.]

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
    sectionRadius = (6.5e-1/6.0) # 6.5Fr diameter = (6.5/3) mm diameter = (6.5/6) mm radius = (6.5/6)e-1 cm radius (NB: 0.1083cm)
    leadInnerRadius = 3.5e-2 / 2.0 # radius of the stylet = 0.35 mm diameter, in cm
    poissonRatio = 0.3
    beamPoissonRatioList = [poissonRatio]*(nbBeams0PlusStock)
    youngModulus = 1.15e5 # in kg.cm-1.s-2 (1 Pa = 1 kg.m-1.s-2 = 0.01 kg.cm-1.s-2)
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
                                     youngModulusList=beamYoungModulusList,
                                     innerRadius=leadInnerRadius)

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
    frameEdges = []
    framesCurvAbscissa = []

    # EXPERIMENTAL: navigation simulation
    # At the beginning of the simulation, when no beam element is actually
    # simulated, all rigid frames are set at 0
    for i in range(nbFrames0):
        frames3DDoFs.append([0, 0, 0])
        frameEdges.append(i)
        frameEdges.append(i+1)
        frames6DDoFs.append([0, 0, 0,  0, 0, 0, 1])
        framesCurvAbscissa.append(0)

    frames3DDoFs.append([0, 0, 0])
    frames6DDoFs.append([0, 0, 0, 0, 0, 0, 1])
    framesCurvAbscissa.append(0)

    # The node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode0.addChild('MappedFrames')
    rateAngularDeformNode0.addChild(mappedFrameNode)

    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                         name="FramesMO", position=frames6DDoFs,
                                         showObject=True, showObjectScale=0.2)

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO0.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=beamCurvAbscissa,
                              curv_abs_output=framesCurvAbscissa, input1=inputMO, input2=inputMO_rigid,
                              output=outputMO, forcefield='@../../rateAngularDeform/beamForceField',
                              drawBeamSegments=True, nonColored=False, debug=0, printLog=False)

    # Collision model

    instrumentCollisionNode0 = mappedFrameNode.addChild('InstrumentCollisionNode')
    instrumentCollisionNode0.addObject('EdgeSetTopologyContainer', name="collisEdgeSet",
                                       position=frames3DDoFs, edges=frameEdges)
    instrumentCollisionNode0.addObject('EdgeSetTopologyModifier', name="colliseEdgeModifier")
    instrumentCollisionNode0.addObject('MechanicalObject', name="CollisionDOFs")
    instrumentCollisionNode0.addObject('LineCollisionModel', bothSide="1", group='1', proximity="0.01")
    instrumentCollisionNode0.addObject('PointCollisionModel', bothSide="1", group='1')
    instrumentCollisionNode0.addObject('IdentityMapping', name="mapping")


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
                                                   showObject=False, showObjectScale=1)

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
    sectionRadius = leadInnerRadius
    poissonRatio = 0.3
    beamPoissonRatioList = [poissonRatio]*(nbBeams1PlusStock)
    youngModulus = 1.5e9 # 150 GPa = 1.5e11 Pa = 1.5e9 kg.cm-1.s-2 (1 Pa = 1 kg.m-1.s-2 = 0.01 kg.cm-1.s-2)
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
    frameEdges = []
    framesCurvAbscissa = []

    # EXPERIMENTAL: navigation simulation
    # At the beginning of the simulation, when no beam element is actually
    # simulated, all rigid frames are set at 0
    for i in range(nbFrames1):
        frames3DDoFs.append([0, 0, 0])
        frameEdges.append(i)
        frameEdges.append(i+1)
        frames6DDoFs.append([0, 0, 0, 0, 0, 0, 1])
        framesCurvAbscissa.append(0)

    frames3DDoFs.append([0, 0, 0])
    frames6DDoFs.append([0, 0, 0, 0, 0, 0, 1])
    framesCurvAbscissa.append(0)

    # The node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode1.addChild('MappedFrames')
    rateAngularDeformNode1.addChild(mappedFrameNode)

    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                         name="FramesMO", position=frames6DDoFs,
                                         showObject=True, showObjectScale=0.2)

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO1.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=beamCurvAbscissa,
                              curv_abs_output=framesCurvAbscissa, input1=inputMO, input2=inputMO_rigid,
                              output=outputMO, forcefield='@../../rateAngularDeform/beamForceField',
                              drawBeamSegments=True, nonColored=False, debug=0, printLog=False)

    # Collision model

    instrumentCollisionNode1 = mappedFrameNode.addChild('InstrumentCollisionNode')
    instrumentCollisionNode1.addObject('EdgeSetTopologyContainer', name="collisEdgeSet",
                                       position=frames3DDoFs, edges=frameEdges)
    instrumentCollisionNode1.addObject('EdgeSetTopologyModifier', name="colliseEdgeModifier")
    instrumentCollisionNode1.addObject('MechanicalObject', name="CollisionDOFs")
    instrumentCollisionNode1.addObject('LineCollisionModel', bothSide="1", group='1', proximity="0.01")
    instrumentCollisionNode1.addObject('PointCollisionModel', bothSide="1", group='1')
    instrumentCollisionNode1.addObject('IdentityMapping', name="mapping")

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
                                                  showObject=False, showObjectScale=1)

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
                                  first_point=[], second_point=[],
                                  newVersionOfFrameComputation=True)

    constraintWith0Node.addObject('CosseratNeedleSlidingConstraint',
                                  name='constraintMappingConstraint',
                                  template="Rigid3d",
                                  useDirections=np.array([0, 1, 1, 0, 0, 0]))



    # -------------------------------------#
    # -----         Anatomy          ----- #
    # -------------------------------------#

    anatomyMeshFile = 'Meshes/nurbs_tortuous_cylinder_mesh.stl'
    forceEstimationSpringStiffness = 1e8

    anatomyRefNode = rootNode.addChild('anatomyRefNode')
    anatomyRefNode.addObject('MeshSTLLoader', name='loader', filename=anatomyMeshFile,
                             flipNormals=True)

    meshTranslation = np.array([-6.5e1, -1.0, 0.])
    meshRotation = np.array([0., 0., 0.])
    scaleFactor = 1.0 # Original mesh is in meters, scene is in mm
    meshScale = np.array([scaleFactor, scaleFactor, scaleFactor])
    anatomyRefNode.addObject("TransformEngine", name="loaderEngine", template="Vec3d",
                             translation=meshTranslation, rotation=meshRotation,
                             scale=meshScale,
                             input_position="@loader.position")

    anatomyRefNode.addObject('MechanicalObject', template='Vec3d', name='anatomyRefDoFs',
                             position="@./loaderEngine.output_position")

    anatomyNode = rootNode.addChild('anatomyNode')

    anatomyNode.addObject('MeshSTLLoader', name='loader', filename=anatomyMeshFile)

    # anatomyNode.addObject('VisualStyle', displayFlags='showBehaviorModels showCollisionModels showForceFields '
    #                                                   'showWireframe')
    anatomyNode.addObject('EulerImplicitSolver', name='odesolver', firstOrder="0",
                          rayleighMass="0.", rayleighStiffness="0.5")
    # anatomyNode.addObject('ShewchukPCGLinearSolver', name='linearSolver', iterations='500', tolerance='1.0e-8',
    #                       preconditioners="precond")
    anatomyNode.addObject('SparseLUSolver', name='precond', template='CompressedRowSparseMatrixd')
    # anatomyNode.addObject('LinearSolverConstraintCorrection', solverName='solver')
    anatomyNode.addObject('MechanicalObject', name='anatomyDoFs', template='Vec3d',
                          showIndices='false', showIndicesScale='4e-5', rx='0',
                          position="@../anatomyRefNode/loaderEngine.output_position",
                          showObject=False, showObjectScale=5.0, showColor='red',
                          printLog=0)
    anatomyNode.addObject('MeshTopology', name="anatomyTopo",
                          position="@anatomyDoFs.position",
                          triangles="@../anatomyRefNode/loader.triangles",
                          drawTriangles='false', drawEdges='false')
    anatomyNode.addObject('TriangleCollisionModel', group='2', bothSide="0",
                          topology='@anatomyTopo')
    anatomyNode.addObject('LineCollisionModel', group='2', bothSide="0")
    anatomyNode.addObject('PointCollisionModel', group='2', bothSide="0")
    # anatomyNode.addObject('OglModel', color="green",
    #                       position="@./loaderEngine.output_position",
    #                       triangles="@./loader.triangles")

    # firstTargetMeshMO = "@../anatomyLoaderNode/anatomyLoader0Node/frame0AnatomyDoFs"
    anatomyNode.addObject('RestShapeSpringsForceField', stiffness=forceEstimationSpringStiffness,
                          external_rest_shape="@../anatomyRefNode/anatomyRefDoFs",
                          mstate="@anatomyDoFs", drawSpring=True)

    anatomyNode.addObject('GenericConstraintCorrection')



    # -------------------------------------------------------------------- #
    # -----                    Python controllers                    ----- #
    # -------------------------------------------------------------------- #

    nbInstruments=2
    instrument0KeyPoints = [42]
    instrument1KeyPoints = [40]
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
