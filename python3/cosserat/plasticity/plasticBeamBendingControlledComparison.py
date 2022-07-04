# -*- coding: utf-8 -*-
"""
Created on Thu Oct 5 2021

@author: ckrewcun
Most of the code is duplicated from python3/cosserat/scenes/testScenes/testCosseratScene
"""

import Sofa
import os

# constants
GRAVITY = 9810
TOT_MASS = 0.1
DENSITY = 0.02
DT = 1e-2

# -------------------------------------#
# -----     Python controller    ----- #
# -------------------------------------#

class BendingController(Sofa.Core.Controller):

    def __init__(self, *args, **kwargs):
        # These are needed (and the normal way to override from a python class)
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.totalTime = 0.0
        self.nbIterations = 0
        self.bendingMoment = 2e8
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
            with root.elaBeamNode.rateAngularDeform.Moment.forces.writeable() as elasticMoment:
                if (self.activated):
                    # Bending moment is currently activated, it should be deactivated
                    print("Bending moment deactivated")
                    moment[0] = [0.0, 0.0, 0.0]
                    elasticMoment[0] = [0.0, 0.0, 0.0]
                    self.activated = False
                else:
                    # self.activated = False => bending moment should be reactivated
                    print("Bending moment activated")
                    moment[0] = [0.0, self.bendingMoment, self.bendingMoment]
                    elasticMoment[0] = [0.0, self.bendingMoment, self.bendingMoment]
                    self.activated = True

# -------------------------------------#
# -----        SOFA scene        ----- #
# -------------------------------------#

pluginNameList = 'Sofa.Component.LinearSolver.Direct' \
                 ' Sofa.Component.MechanicalLoad' \
                 ' Sofa.Component.ODESolver.Backward' \
                 ' Sofa.Component.SolidMechanics.Spring' \
                 ' Sofa.Component.StateContainer' \
                 ' Sofa.Component.Visual' \
                 ' SofaPython3 CosseratPlugin'

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

    # ----- Python controller ----- #

    rootNode.addObject(BendingController(name="BendingController"))

    # -------------------------------------#
    # -----  Plastic Cosserat beam   ----- #
    # -------------------------------------#

    # ----- Cosserat beam parameters ----- #

    # Define: the total length of the beam
    totalLength = 100.0

    # Define: the number of section, the total length and the length of each beam.
    nbBeams = 1
    oneBeamLength = totalLength / nbBeams

    # Define: the number of frame and the length between each frame.
    nbFrames = 20
    distBetweenFrames = totalLength / nbFrames

    # ----- Rigid base ----- #

    beamNode = rootNode.addChild('beamNode')
    beamNode.addObject('EulerImplicitSolver', rayleighStiffness="1.2", rayleighMass='1.1')
    beamNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")

    rigidBaseNode= beamNode.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode.addObject('MechanicalObject', template='Rigid3d',
                                          name="RigidBaseMO", position=[0., 0., 0., 0, 0, 0, 1],
                                          showObject=1,
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
    for i in range(nbBeams):
        beamStrainDoFs.append([0, 0, 0])
        beamLengths.append(oneBeamLength)
        sum += beamLengths[i]
        beamCurvAbscissa.append(sum)
    beamCurvAbscissa[nbBeams] = totalLength

    # Define angular rate which is the torsion(x) and bending (y, z) of each section
    rateAngularDeformNode = beamNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject('MechanicalObject',
                                                          template='Vec3d',
                                                          name='rateAngularDeformMO',
                                                          position=beamStrainDoFs,
                                                          showIndices=0,
                                                          rest_position=[0.0, 0.0, 0.0])

    beamCrossSectionShape='circular'
    sectionRadius = 2.0
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
    rateAngularDeformNode.addObject('BeamPlasticLawForceField', name="beamForceField",
                                    crossSectionShape=beamCrossSectionShape,
                                    radius=sectionRadius, variantSections="true",
                                    length=beamLengths, poissonRatioList=beamPoissonRatioList,
                                    youngModulusList=beamYoungModulusList,
                                    initialYieldStresses=yieldStressList,
                                    plasticModuli=plasticModulusList,
                                    mixedHardeningCoefficients= hardeningCoefficientList)

    rateAngularDeformNode.addObject('ConstantForceField', name='Moment', indices="0",
                                    forces="0 2e8 2e8")

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

    # -------------------------------------#
    # -----  Elastic Cosserat beam   ----- #
    # -------------------------------------#

    # ----- Cosserat beam parameters ----- #

    # Same parameters as for the plastic beam :
    #   totalLength = 100.0
    #   nbBeams = 1
    #   oneBeamLength = totalLength / nbBeams
    #   nbFrames = 20
    #   distBetweenFrames = totalLength / nbFrames

    # ----- Rigid base ----- #

    elaBeamNode = rootNode.addChild('elaBeamNode')
    elaBeamNode.addObject('EulerImplicitSolver', rayleighStiffness="1.2", rayleighMass='1.1')
    elaBeamNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")

    rigidBaseNode= elaBeamNode.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode.addObject('MechanicalObject', template='Rigid3d',
                                          name="RigidBaseMO", position=[0., 0., 10., 0, 0, 0, 1],
                                          showObject=1,
                                          showObjectScale=2.)
    rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="5.e8", angularStiffness="5.e8",
                            external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")

    # ----- Rate of angular deformation ----- #
    # Define: the length of each beam in a list, the positions of each beam

    # beamStrainDoFs, beamLengths, beamCurvAbscissa are the same as for the plastic beam

    rateAngularDeformNode = elaBeamNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject('MechanicalObject',
                                                          template='Vec3d',
                                                          name='rateAngularDeformMO',
                                                          position=beamStrainDoFs,
                                                          showIndices=0,
                                                          rest_position=[0.0, 0.0, 0.0])

    # Elastic behaviour law
    rateAngularDeformNode.addObject('BeamHookeLawForceField',
                                    name='hookeBeamForceField',
                                    crossSectionShape='circular', length=beamLengths,
                                    radius=2., youngModulus=5e6)
    rateAngularDeformNode.addObject('ConstantForceField', name='Moment', indices="0",
                                    forces="0 2e8 2e8")

    # ----- Frames ----- #

    # Define local frames related to each section and parameters curv_abs_outputF
    frames6DDoFs = []
    frames3DDoFs = []
    framesCurvAbscissa = []
    for i in range(nbFrames):
        sol = i * distBetweenFrames
        frames3DDoFs.append([sol, 0, 10])
        frames6DDoFs.append([sol, 0, 10,  0, 0, 0, 1])
        framesCurvAbscissa.append(sol)

    frames6DDoFs.append([totalLength, 0, 10, 0, 0, 0, 1])
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
                              output=outputMO, forcefield='@../../rateAngularDeform/hookeBeamForceField',
                              nonColored=False, debug=0)


    return rootNode
