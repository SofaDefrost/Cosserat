# -*- coding: utf-8 -*-
"""
Created on April 22 2022

@author: ckrewcun
"""

import Sofa
import os
import numpy as np

# constants
GRAVITY = 0.0 # 9810
TOT_MASS = 0.1
DENSITY = 0.02
DT = 1e-2

# -------------------------------------#
# -----     Python controller    ----- #
# -------------------------------------#

class BendingController(Sofa.Core.Controller):

    def __init__(self, bendingMoment, *args, **kwargs):
        # These are needed (and the normal way to override from a python class)
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.totalTime = 0.0
        self.nbIterations = 0
        self.bendingMoment = bendingMoment
        self.activated = True

    def onAnimateBeginEvent(self, event):  # called at each begin of animation step
        if self.nbIterations == 2030:
            self.triggerBendingMoment()

        self.totalTime = self.totalTime + DT
        self.nbIterations = self.nbIterations + 1

    def onAnimateEndEvent(self, event):  # called at each end of animation step
        pass

    def onKeypressedEvent(self, event):
        # Press L key triggers the creation of new objects in the scene
        if event['key'] == 'F':
            self.triggerBendingMoment()

    def triggerBendingMoment(self):
        root = self.getContext()
        with root.beamNode.rateAngularDeform.Moment.forces.writeable() as moment:
            if (self.activated):
                # Bending moment is currently activated, it should be deactivated
                print("Bending moment deactivated")
                moment[0] = [0.0, 0.0, 0.0]
                self.activated = False
            else:
                # self.activated = False => bending moment should be reactivated
                print("Bending moment activated")
                moment[0] = [0.0, self.bendingMoment, self.bendingMoment]
                self.activated = True

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

    # ----- Plastic Cosserat beam parameters ----- #

    # Define: the total length of the beam
    totalLength = 15

    # Define: the number of section, the total length and the length of each beam.
    nbBeams = 3
    oneBeamLength = totalLength / nbBeams

    # Define: the number of frame and the length between each frame.
    nbFrames = 12
    distBetweenFrames = totalLength / nbFrames

    # ----- Rigid base ----- #

    beamNode = rootNode.addChild('beamNode')
    beamNode.addObject('EulerImplicitSolver', rayleighStiffness="1.2", rayleighMass='1.1')
    # beamNode.addObject('CGLinearSolver', name="solver",
    #                    iterations='100', tolerance='1e-5', threshold='1e-5')
    beamNode.addObject('SparseLUSolver',
                               template='CompressedRowSparseMatrixd',
                               printLog="false")

    rigidBaseNode= beamNode.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode.addObject('MechanicalObject', template='Rigid3d',
                                          name="RigidBaseMO", position=[0., 0., 0., 0, 0, 0, 1], showObject=1,
                                          showObjectScale=2.)
    rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="5.e8", angularStiffness="5.e8",
                            external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")

    # ----- Rate of angular deformation ----- #
    # Define the length of each beam in a list, the positions of each beam



    beamStrainDoFs = []
    beamLengths = []
    sum = 0.
    beamCurvAbscissa = []
    beamCurvAbscissa.append(0.0)

    # EXPERIMENTAL: insertion navigation
    # Increasing the number of DoFs and the corresponding parameters.
    # For each additional beam, we:
    #     - Add a Vec3 in beamStrainDoFs, with each strain component set at 0.
    #     - Add 0 in the length data field for the Cosserat mapping
    #     - Add the last curvilinear abscissa (total length) in the beams
    # curvilinear abcissas in the Cosserat mapping.
    # This is equivalent to saying that the additional beams are 0-length,
    # 0-strain beams, added at the end of the instrument (in case of navigation)
    nbStockBeams = 1
    for i in range(0, nbStockBeams):
        beamStrainDoFs.append([0, 0, 0])
        beamLengths.append(0)
        beamCurvAbscissa.append(0)

    for i in range(nbBeams):
        beamStrainDoFs.append([0, 0, 0])
        beamLengths.append(oneBeamLength)
        sum += beamLengths[i+nbStockBeams]
        beamCurvAbscissa.append(sum)
    beamCurvAbscissa[nbBeams+nbStockBeams] = totalLength

    # Define angular rate which is the torsion(x) and bending (y, z) of each section
    rateAngularDeformNode = beamNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject('MechanicalObject',
                                                          template='Vec3d',
                                                          name='rateAngularDeformMO',
                                                          position=beamStrainDoFs,
                                                          showIndices=0,
                                                          rest_position=[0.0, 0.0, 0.0])

    beamCrossSectionShape='circular'
    sectionRadius = 0.5
    poissonRatio = 0.45
    beamPoissonRatioList = [poissonRatio]*(nbBeams+1)
    youngModulus = 5.0e6
    beamYoungModulusList = [youngModulus]*(nbBeams+1)
    yieldStress = 5.0e4
    yieldStressList = [yieldStress]*nbBeams
    plasticModulus = 2.0e5
    plasticModulusList = [plasticModulus]*nbBeams
    hardeningCoeff = 0.5
    hardeningCoefficientList = [hardeningCoeff]*nbBeams
    rateAngularDeformNode.addObject('BeamPlasticLawForceField', name="beamForceField",
                                    crossSectionShape=beamCrossSectionShape,
                                    radius=sectionRadius, varianteSections="true",
                                    length=beamLengths, poissonRatioList=beamPoissonRatioList,
                                    youngModulusList=beamYoungModulusList,
                                    initialYieldStresses=yieldStressList,
                                    plasticModuli=plasticModulusList,
                                    mixedHardeningCoefficients= hardeningCoefficientList)

    beamBendingMoment = 1.0e5
    bendingForces = np.array([0, beamBendingMoment, beamBendingMoment])
    # momentIndices = range(1, nbBeams)
    momentIndices = [nbBeams-1]
    # rateAngularDeformNode.addObject('ConstantForceField', name='Moment',
    #                                 indices=momentIndices,
    #                                 forces=bendingForces)

    # EXPERIMENTAL: navigation simulation
    # Adding constraints on the additional beams which are not meant to be
    # simulated yieldStress
    # fixedIndices = list(range(nbBeams, nbBeams+nbStockBeams))
    fixedIndices = list(range(0, nbStockBeams))
    rateAngularDeformNode.addObject('FixedConstraint', name='FixedConstraint',
                                    indices=fixedIndices)

    # ----- Frames ----- #

    # Define local frames related to each section and parameters framesCurvAbscissa
    frames6DDoFs = []
    frames3DDoFs = []
    framesCurvAbscissa = []
    for i in range(nbFrames):
        sol = i * distBetweenFrames
        frames3DDoFs.append([sol, 0, 0])
        frames6DDoFs.append([sol, 0, 0,  0, 0, 0, 1])
        framesCurvAbscissa.append(sol)

    frames6DDoFs.append([totalLength, 0, 0, 0, 0, 0, 1])
    framesCurvAbscissa.append(totalLength)

    # The node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.addChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)

    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                         name="FramesMO", position=frames6DDoFs,
                                         showObject=1, showObjectScale=1)

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=beamCurvAbscissa,
                              curv_abs_output=framesCurvAbscissa, input1=inputMO, input2=inputMO_rigid,
                              output=outputMO, forcefield='@../../rateAngularDeform/beamForceField',
                              nonColored=False, debug=0)

    mappedFrameNode.addObject('ConstantForceField', name='Moment',
                              indices=nbFrames-4,
                              forces=np.array([0, 0, 0, 0, 0, 8e4]))

    # ----- Python controller ----- #

    rootNode.addObject(BendingController(name="BendingController",
                                         bendingMoment=beamBendingMoment))

    return rootNode