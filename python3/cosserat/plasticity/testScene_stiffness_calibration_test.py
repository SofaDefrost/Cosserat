# -*- coding: utf-8 -*-
"""
Created on Februrary 15 2023

@author: ckrewcun
"""

import Sofa
import os
import numpy as np
from CosseratNavigationController import CombinedInstrumentsController
from CosseratBeamCreator import *  # generateRegularSectionsAndFrames

# ----- Controller definition ----- #

# Controller to apply a force on the tip of Cosserat beam instrument model:
#   - Ctrl + 0, Ctrl + 1 : switch the currently controlled instrument
#   - Ctrl + F : Activate/deactivate force on the tip of the instrument
#
# Required structure of the scene :
# * rootNode
#
#   * controlPointNode0
#       controlPointMO (Rigid)
#   * cosseratBeamNode
#       MechanicalMatrixMapper
#       * rigidBaseNode
#           RigidBaseMO (Rigid)
#           * MappedFrames
#               FramesMO (Rigid)
#               controlSpring (RestShapeSpringsForceField)
#               mapping (DiscreteCosseratMapping)
#       * rateAngularDeform
#           rateAngularDeformMO (Cosserat strains)
#           bendingMoment (ConstantForceField)
#           fixation (FixedConstraint)
class ForceController(Sofa.Core.Controller):

    def __init__(self, rootNode, solverNode, instrumentNodeList, forceList,
                 *args, **kwargs):

        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.rootNode = rootNode
        self.solverNode = solverNode
        self.instrumentNodeList = instrumentNodeList

        if len(instrumentNodeList) != len(forceList):
            raise ValueError("[ForceController]: parameters instrumentNodeList "
                            "and forceList should have the same size (one Vec6 "
                            "should be provided in forceList for each instrument "
                            "in instrumentNodeList. ")

        self.forceList = forceList

        self.totalTime = 0.0
        self.nbIterations = 0

        self.activatedForceList = [False]*len(instrumentNodeList)
        self.currentInstrumentId = 0


    def onAnimateBeginEvent(self, event):  # called at each begin of animation step

        self.totalTime = self.totalTime + self.rootNode.findData('dt').value
        self.nbIterations = self.nbIterations + 1

    def onAnimateEndEvent(self, event):  # called at each end of animation step
        pass

    def onKeypressedEvent(self, event):

        if event['key'] == 'F':
            if (self.activatedForceList[self.currentInstrumentId]):
                print("turning off the force on instrument {}".format(self.currentInstrumentId))
                self.deactivateForce(self.currentInstrumentId)
            else:
                print("turning on the force on instrument {}".format(self.currentInstrumentId))
                self.activateForce(self.currentInstrumentId)
            self.activatedForceList[self.currentInstrumentId] = not self.activatedForceList[self.currentInstrumentId]

        if event['key'] == '0':
            self.currentInstrumentId = 0
            print("Currently controlled: instrument 0")

        if event['key'] == '1':
            if len(self.instrumentNodeList) <= 1:
                warnings.warn("Instrument number 1 doesn't exist (only one instrument (0) is available).".format(self.nbInstruments))
            else:
                self.currentInstrumentId = 1
                print("Currently controlled: instrument 1")


    def activateForce(self, instrumentId):
        # Retrieve the instrument node
        instrumentNode = self.instrumentNodeList[instrumentId]

        rigidBaseNode = instrumentNode.getChild('rigidBase')
        if rigidBaseNode is None:
            raise NameError("[ForceController]: Node \'rigidBase\' "
                            "not found. Your scene should contain a node named "
                            "\'rigidBase\' among the children of the current "
                            "instrument ({}) node, "
                            "where the base and rigid frames of the "
                            "Cosserat model are defined".format(self.currentInstrumentId))

        mappedFramesNode = rigidBaseNode.getChild('mappedFrames')
        if mappedFramesNode is None:
            raise NameError("[ForceController]: Node \'mappedFrames\' "
                            "not found. The \'rigidBase\' node should have a child "
                            "node called \'mappedFrames\', in which the Cosserat "
                            "mapping output frames should be defined.")

        # Set force value
        mappedFramesNode.externalForce.force = self.forceList[instrumentId]

    def deactivateForce(self, instrumentId):
        # Retrieve the instrument node
        instrumentNode = self.instrumentNodeList[instrumentId]

        rigidBaseNode = instrumentNode.getChild('rigidBase')
        if rigidBaseNode is None:
            raise NameError("[ForceController]: Node \'rigidBase\' "
                            "not found. Your scene should contain a node named "
                            "\'rigidBase\' among the children of the current "
                            "instrument ({}) node, "
                            "where the base and rigid frames of the "
                            "Cosserat model are defined".format(self.currentInstrumentId))

        mappedFramesNode = rigidBaseNode.getChild('mappedFrames')
        if mappedFramesNode is None:
            raise NameError("[ForceController]: Node \'mappedFrames\' "
                            "not found. The \'rigidBase\' node should have a child "
                            "node called \'mappedFrames\', in which the Cosserat "
                            "mapping output frames should be defined.")

        # Set force value
        mappedFramesNode.externalForce.force = np.array([0., 0., 0., 0., 0., 0.])




# ----- Constants ----- #

GRAVITY = 0. # 981
TOT_MASS = 0.1
DENSITY = 0.02
DT = 1e-2


# ----- Instrument parameters ----- #

instrument0TotLength = 10.0 # in cm
instrument0NbBeams = 3
instrument0NbFrames = 10
beamCreatorResult = generateRegularSectionsAndFrames(instrument0TotLength,
                                                     instrument0NbBeams,
                                                     instrument0NbFrames)

instrument0BeamLengths = beamCreatorResult['sectionLengths']
instrument0BeamCurvAbsInput = beamCreatorResult['curvAbsInput']
instrument0BeamStrainDoFs = beamCreatorResult['sectionDoFs']
instrument0BeamFrameRigidDoFs = beamCreatorResult['frameRigidDoFs']
instrument0BeamFrame3DDoFs = beamCreatorResult['frame3DDoFs']
instrument0BeamFrameEdges = beamCreatorResult['frameEdges']
instrument0BeamCurvAbsOutput = beamCreatorResult['curvAbsOutput']



instrument1TotLength = 10.0 # in cm
instrument1NbBeams = 3
instrument1NbFrames = 10
beamCreatorResult = generateRegularSectionsAndFrames(instrument1TotLength,
                                                     instrument1NbBeams,
                                                     instrument1NbFrames)

instrument1BeamLengths = beamCreatorResult['sectionLengths']
instrument1BeamCurvAbsInput = beamCreatorResult['curvAbsInput']
instrument1BeamStrainDoFs = beamCreatorResult['sectionDoFs']
instrument1BeamFrameRigidDoFs = beamCreatorResult['frameRigidDoFs']
instrument1BeamFrame3DDoFs = beamCreatorResult['frame3DDoFs']
instrument1BeamFrameEdges = beamCreatorResult['frameEdges']
instrument1BeamCurvAbsOutput = beamCreatorResult['curvAbsOutput']

# Common
totalMass = 0.022



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
    rootNode.findData('gravity').value = [0., -GRAVITY, 0.]

    rootNode.addObject('DefaultAnimationLoop')


    # -------------------------------------#
    # -----      Control points      ----- #
    # -------------------------------------#


    # Instrument 0
    controlPointPos0 = np.array([0., 0., 0., 0., 0., 0., 1.])
    controlPointNode0 = rootNode.addChild('controlPointNode0')
    controlPointNode0.addObject('MechanicalObject', template='Rigid3d', name="controlPointMO",
                               position=controlPointPos0, showObject='1', showObjectScale='0.5')

    # Instrument1
    controlPointPos1 = np.array([0., 0., 3., 0., 0., 0., 1.])
    controlPointNode1 = rootNode.addChild('controlPointNode1')
    controlPointNode1.addObject('MechanicalObject', template='Rigid3d', name="controlPointMO",
                               position=controlPointPos1, showObject='1', showObjectScale='0.5')


    # -------------------------------------------------------------------- #
    # -----                        Solver node                       ----- #
    # -------------------------------------------------------------------- #

    solverNode = rootNode.addChild('solverNode')

    stiffnessDamping = 0.5
    massDamping = 0.
    solverNode.addObject('EulerImplicitSolver',
                         rayleighStiffness=stiffnessDamping,
                         rayleighMass=massDamping)
    # solverNode.addObject('CGLinearSolver', name="solver",
    #                      iterations='2000', tolerance='1e-8', threshold='1e-12')
    solverNode.addObject('SparseLUSolver',
                         template='CompressedRowSparseMatrixd',
                         printLog="false")

    # -------------------------------------------------------------------- #
    # -----                First instrument components               ----- #
    # -------------------------------------------------------------------- #

    # ----- Rigid base ----- #

    instrument0Node = solverNode.addChild('Instrument0')

    rigidBaseNode0 = instrument0Node.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode0.addObject('MechanicalObject', template='Rigid3d',
                                           name="RigidBaseMO",
                                           position=controlPointPos0,
                                           showObject=False, showObjectScale=2.)
    rigidBaseNode0.addObject('RestShapeSpringsForceField', name='controlSpring',
                             stiffness="5.e8", angularStiffness="1.0e5",
                             external_rest_shape="@../../../controlPointNode0/controlPointMO",
                             external_points="0", mstate="@RigidBaseMO",
                             points="0", template="Rigid3d")


    # ----- Rate of angular deformation ----- #

    rateAngularDeformNode0 = instrument0Node.addChild('rateAngularDeform')
    rateAngularDeformMO0 = rateAngularDeformNode0.addObject('MechanicalObject',
                                                            template='Vec3d',
                                                            name='rateAngularDeformMO',
                                                            position=instrument0BeamStrainDoFs,
                                                            showIndices=0,
                                                            rest_position=instrument0BeamStrainDoFs)

    beamCrossSectionShape='circular'
    sectionRadius = (6.5e-1/6.0) # 6.5Fr diameter = (6.5/3) mm diameter = (6.5/6) mm radius = (6.5/6)e-1 cm radius (NB: 0.1083cm)
    leadInnerRadius = 3.5e-2 / 2.0 # radius of the stylet = 0.35 mm diameter, in cm
    poissonRatio = 0.3
    beamPoissonRatioList = [poissonRatio]*(instrument0NbBeams)
    youngModulus = 1.15e5 # in kg.cm-1.s-2 (1 Pa = 1 kg.m-1.s-2 = 0.01 kg.cm-1.s-2)
    beamYoungModulusList = [youngModulus]*(instrument0NbBeams)
    yieldStress = 5.0e4
    yieldStressList = [yieldStress]*(instrument0NbBeams)
    plasticModulus = 2.0e5
    plasticModulusList = [plasticModulus]*(instrument0NbBeams)
    hardeningCoeff = 0.5
    hardeningCoefficientList = [hardeningCoeff]*(instrument0NbBeams)
    ### Plastic FF version
    # rateAngularDeformNode0.addObject('BeamPlasticLawForceField', name="beamForceField",
    #                                  crossSectionShape=beamCrossSectionShape,
    #                                  radius=sectionRadius, variantSections="true",
    #                                  length=instrument0BeamLengths, poissonRatioList=beamPoissonRatioList,
    #                                  youngModulusList=beamYoungModulusList,
    #                                  initialYieldStresses=yieldStressList,
    #                                  plasticModuli=plasticModulusList,
    #                                  mixedHardeningCoefficients= hardeningCoefficientList)
    ### Elastic FF version
    rateAngularDeformNode0.addObject('BeamHookeLawForceField', name="beamForceField",
                                     crossSectionShape=beamCrossSectionShape,
                                     radius=sectionRadius, variantSections="true",
                                     length=instrument0BeamLengths, poissonRatioList=beamPoissonRatioList,
                                     youngModulusList=beamYoungModulusList,
                                     innerRadius=leadInnerRadius)

    # ----- Frames ----- #

    # The node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode0.addChild('mappedFrames')
    rateAngularDeformNode0.addChild(mappedFrameNode)

    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                         name="FramesMO",
                                         position=instrument0BeamFrameRigidDoFs,
                                         showObject=True, showObjectScale=0.2)

    mappedFrameNode.addObject('UniformMass', totalMass=totalMass,
                              name="UniformMass0", showAxisSizeFactor='0.')

    mappedFrameNode.addObject('ConstantForceField', name="externalForce",
                              indices=instrument0NbFrames-1,
                              force=np.array([0., 0., 0., 0., 0., 0.]),
                              showArrowSize=1.0)

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO0.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping',
                              curv_abs_input=instrument0BeamCurvAbsInput,
                              curv_abs_output=instrument0BeamCurvAbsOutput,
                              input1=inputMO, input2=inputMO_rigid,
                              output=outputMO,
                              forcefield='@../../rateAngularDeform/beamForceField',
                              drawBeamSegments=True, nonColored=False, debug=0,
                              printLog=False)




    # -------------------------------------------------------------------- #
    # -----               Second instrument components               ----- #
    # -------------------------------------------------------------------- #


    # ----- Rigid base ----- #

    instrument1Node = solverNode.addChild('Instrument1')

    rigidBaseNode1 = instrument1Node.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode1.addObject('MechanicalObject', template='Rigid3d',
                                           name="RigidBaseMO",
                                           position=controlPointPos1,
                                           showObject=False, showObjectScale=2.)
    rigidBaseNode1.addObject('RestShapeSpringsForceField', name='controlSpring',
                             stiffness="5.e8", angularStiffness="1.0e5",
                             external_rest_shape="@../../../controlPointNode1/controlPointMO",
                             external_points="0", mstate="@RigidBaseMO",
                             points="0", template="Rigid3d")

    # ----- Rate of angular deformation ----- #

    # Define angular rate which is the torsion(x) and bending (y, z) of each section
    rateAngularDeformNode1 = instrument1Node.addChild('rateAngularDeform')
    rateAngularDeformMO1 = rateAngularDeformNode1.addObject('MechanicalObject',
                                                            template='Vec3d',
                                                            name='rateAngularDeformMO',
                                                            position=instrument1BeamStrainDoFs,
                                                            showIndices=0,
                                                            rest_position=instrument1BeamStrainDoFs)

    beamCrossSectionShape='circular'
    sectionRadius = leadInnerRadius
    poissonRatio = 0.3
    beamPoissonRatioList = [poissonRatio]*(instrument1NbBeams)
    youngModulus = 1.5e9 # 150 GPa = 1.5e11 Pa = 1.5e9 kg.cm-1.s-2 (1 Pa = 1 kg.m-1.s-2 = 0.01 kg.cm-1.s-2)
    beamYoungModulusList = [youngModulus]*(instrument1NbBeams)
    yieldStress = 5.0e4
    yieldStressList = [yieldStress]*(instrument1NbBeams)
    plasticModulus = 2.0e5
    plasticModulusList = [plasticModulus]*(instrument1NbBeams)
    hardeningCoeff = 0.5
    hardeningCoefficientList = [hardeningCoeff]*(instrument1NbBeams)
    ### Plastic FF version
    # rateAngularDeformNode1.addObject('BeamPlasticLawForceField', name="beamForceField",
    #                                  crossSectionShape=beamCrossSectionShape,
    #                                  radius=sectionRadius, variantSections="true",
    #                                  length=instrument1BeamLengths, poissonRatioList=beamPoissonRatioList,
    #                                  youngModulusList=beamYoungModulusList,
    #                                  initialYieldStresses=yieldStressList,
    #                                  plasticModuli=plasticModulusList,
    #                                  mixedHardeningCoefficients= hardeningCoefficientList)
    ### Elastic FF version
    rateAngularDeformNode1.addObject('BeamHookeLawForceField', name="beamForceField",
                                     crossSectionShape=beamCrossSectionShape,
                                     radius=sectionRadius, variantSections="true",
                                     length=instrument1BeamLengths, poissonRatioList=beamPoissonRatioList,
                                     youngModulusList=beamYoungModulusList)

    # ----- Frames ----- #

    # The node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode1.addChild('mappedFrames')
    rateAngularDeformNode1.addChild(mappedFrameNode)

    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                         name="FramesMO",
                                         position=instrument1BeamFrameRigidDoFs,
                                         showObject=True, showObjectScale=0.2)

    mappedFrameNode.addObject('UniformMass', totalMass=totalMass,
                              name="UniformMass1", showAxisSizeFactor='0.')

    mappedFrameNode.addObject('ConstantForceField', name="externalForce",
                              indices=instrument1NbFrames-1,
                              force=np.array([0., 0., 0., 0., 0., 0.]),
                              showArrowSize=1.0)

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO1.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping',
                              curv_abs_input=instrument1BeamCurvAbsInput,
                              curv_abs_output=instrument1BeamCurvAbsOutput,
                              input1=inputMO, input2=inputMO_rigid,
                              output=outputMO,
                              forcefield='@../../rateAngularDeform/beamForceField',
                              drawBeamSegments=True, nonColored=False, debug=0,
                              printLog=False)


    # -------------------------------------------------------------------- #
    # -----                    Python controllers                    ----- #
    # -------------------------------------------------------------------- #

    sceneInstrumentNodeList = [instrument0Node, instrument1Node]
    forceIntensity = 2.0e-1 # 1e-3N in cm.kg.s-2 (1 cm.kg.s-2 = 1e-2 m.kg.s-2)
    forceDirection = np.array([0., -1., 0., 0., 0., 0.])
    instrument0Force = forceIntensity * forceDirection
    forceIntensity = 1.0 # 1e-3N in cm.kg.s-2 (1 cm.kg.s-2 = 1e-2 m.kg.s-2)
    forceDirection = np.array([0., -1., 0., 0., 0., 0.])
    instrument1Force = forceIntensity * forceDirection
    sceneForceList = [instrument0Force, instrument1Force]

    rootNode.addObject(ForceController(name="ForceController",
                                       rootNode=rootNode,
                                       solverNode=solverNode,
                                       instrumentNodeList=sceneInstrumentNodeList,
                                       forceList=sceneForceList))

    return rootNode
