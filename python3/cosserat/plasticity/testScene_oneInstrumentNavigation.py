# -*- coding: utf-8 -*-
"""
Created on May 12 2022

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

    #----------------------------------------------#
    # -----           Cosserat beam          ----- #
    #----------------------------------------------#

    beamNode = rootNode.addChild('beamNode')

    #----- Solvers -----#

    beamNode.addObject('EulerImplicitSolver', rayleighStiffness="1.2", rayleighMass='1.1')
    # beamNode.addObject('CGLinearSolver', name="solver",
    #                    iterations='100', tolerance='1e-5', threshold='1e-5')
    beamNode.addObject('SparseLUSolver',
                               template='CompressedRowSparseMatrixd',
                               printLog="false")

    #----- Component 1 = rigid base of the beams -----#

    rigidBaseNode= beamNode.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode.addObject('MechanicalObject', template='Rigid3d',
                                          name="RigidBaseMO", position=[0., 0., 0., 0, 0, 0, 1], showObject=1,
                                          showObjectScale=2.)
    rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="5.e8", angularStiffness="5.e8",
                            external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")

    #----- Beam property definition -----#

    # /!\ In addition to the rigid base defined above, a Cosserat model generally
    # requires 4 components :
    #  - A first MechanicalObject, composed of 3D elements, representing the
    # 'strain' degrees of freedom (not equivalent to classic 3D position DoFs)
    #  - A ForceField, specific to these strain degrees of freedom, implementing
    # the mechanical behaviour of the Cosserat model.
    #  - A second MechanicalObject, composed of 6D rigid frames, used to
    # represent the position and orientation of the beams in the global frame.
    # These frames are computed from the strain DoFs of the Cosserat model.
    #  -  A dedicated mapping (e.g. DiscreteCosseratMapping), implementing the
    # computation to switch from the strain DoFs to the rigid representation

    # Before creating these 4 components, we define a set of properties which
    # they depend on.

    # Parameters for the Cosserat strain DoFs
    beamStrainDoFs = []
    beamLengths = []
    beamCurvAbscissa = []
    # This data, which will be passed as 'curv_abs_input' to the
    # DiscreteCosseratMapping, contains the curvilinear abscissas of all the beam
    # extremities in the Cosserat model. Therefore it's first element is always 0,
    # and its size is nbBeams+1
    beamCurvAbscissa.append(0.)

    # EXPERIMENTAL: insertion navigation
    # In this example, the Cosserat model is used to represent slender interventional
    # instruments, which are pushed or pulled from an insertion point. At the
    # beginning of the scene, the beams are defined with 0-length and 0-strain,
    # as the instruments are fully retracted. These quantities are meant to be
    # modified by a controller during the simulation of navigation.
    # For each beam, we have:
    #     - Add a Vec3 in beamStrainDoFs, with each strain component set at 0.
    #     - Add 0 in the length data field for the Cosserat mapping
    #     - Add the last curvilinear abscissa (total length) in the beams
    # curvilinear abcissas in the Cosserat mapping.
    # This is equivalent to saying that these beams are 0-length,
    # 0-strain beams, added at the end of the instrument.
    nbBeams = 3
    for i in range(0, nbBeams):
        beamStrainDoFs.append([0., 0., 0.])
        beamLengths.append(0.)
        beamCurvAbscissa.append(0.)

    # Parameters for the forcefield (mechanical behaviour)
    beamCrossSectionShape='circular'
    sectionRadius = 0.5
    poissonRatio = 0.45
    beamPoissonRatioList = [poissonRatio]*nbBeams
    youngModulus = 5.0e6
    beamYoungModulusList = [youngModulus]*nbBeams
    yieldStress = 5.0e4
    yieldStressList = [yieldStress]*nbBeams
    plasticModulus = 2.0e5
    plasticModulusList = [plasticModulus]*nbBeams
    hardeningCoeff = 0.5
    hardeningCoefficientList = [hardeningCoeff]*nbBeams

    # Parameters for the 6D rigid frames

    # EXPERIMENTAL: insertion navigation
    # Following the same idea as for the Cosserat strain DoFs, we initially
    # define all the desired frames at the same point, with a neutral orientation.
    # These positions and orientations will be automatically updated by the
    # Cosserat mapping when the navigation controller will update the frames
    # curvilinear abscissas in the mapping component.
    nbFrames = nbBeams*2 + 1
    frames6DDoFs = []
    frames3DDoFs = []
    framesCurvAbscissa = []
    for i in range(nbFrames):
        frames3DDoFs.append([0., 0., 0.])
        frames6DDoFs.append([0., 0., 0., 0., 0., 0., 1.])
        framesCurvAbscissa.append(0.)

    frames6DDoFs.append([0., 0., 0., 0., 0., 0., 1.])
    framesCurvAbscissa.append(0.)

    #----- Components 2 and 3 = Cosserat strain DoFs and forcefield -----#

    rateAngularDeformNode = beamNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject('MechanicalObject',
                                                          template='Vec3d',
                                                          name='rateAngularDeformMO',
                                                          position=beamStrainDoFs,
                                                          showIndices=0,
                                                          rest_position=[0.0, 0.0, 0.0])


    rateAngularDeformNode.addObject('BeamPlasticLawForceField', name="beamForceField",
                                    crossSectionShape=beamCrossSectionShape,
                                    radius=sectionRadius, varianteSections="true",
                                    length=beamLengths, poissonRatioList=beamPoissonRatioList,
                                    youngModulusList=beamYoungModulusList,
                                    initialYieldStresses=yieldStressList,
                                    plasticModuli=plasticModulusList,
                                    mixedHardeningCoefficients= hardeningCoefficientList)

    # EXPERIMENTAL: navigation simulation
    # Adding constraints on the additional beams which are not meant to be
    # simulated yieldStress
    fixedIndices = list(range(0, nbBeams))
    rateAngularDeformNode.addObject('FixedConstraint', name='FixedConstraint',
                                    indices=fixedIndices)

    # ----- Component 4 = 6D rigid frames for position and orientation ----- #

    # The node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.addChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)

    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                         name="FramesMO", position=frames6DDoFs,
                                         showObject=1, showObjectScale=1)


    # ----- Component 5 = Cosserat mapping ----- #

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO, and one
    # output: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=beamCurvAbscissa,
                              curv_abs_output=framesCurvAbscissa, input1=inputMO, input2=inputMO_rigid,
                              output=outputMO, forcefield='@../../rateAngularDeform/beamForceField',
                              nonColored=False, debug=0)

    # ----- Python controller ----- #

    # Navigation Parameters
    nbInstruments = 1
    instrumentBeamDensityVect = np.array([20])
    incrementDistance = 0.1
    isInstrumentStraightVect = [True]
    curvAbsTolerance = 1.0e-5
    instrumentLengths = np.array([15.])

    rootNode.addObject(CombinedInstrumentsController(name="CombinedInstrumentsController",
                                                     rootNode=rootNode,
                                                     nbInstruments=nbInstruments,
                                                     instrumentBeamDensityVect=instrumentBeamDensityVect,
                                                     incrementDistance=incrementDistance,
                                                     isInstrumentStraightVect=isInstrumentStraightVect,
                                                     curvAbsTolerance=curvAbsTolerance,
                                                     instrumentLengths=instrumentLengths))

    return rootNode
