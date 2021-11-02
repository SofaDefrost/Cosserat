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

class PythonController(Sofa.Core.Controller):

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
        with root.cableNode.rateAngularDeform.Moment.forces.writeable() as moment:
            with root.elaCableNode.rateAngularDeform.Moment.forces.writeable() as elasticMoment:
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

def createScene(rootNode):

    rootNode.addObject('RequiredPlugin', name='SoftRobots')
    rootNode.addObject('RequiredPlugin', name='BeamAdapter')
    rootNode.addObject('RequiredPlugin', name='SofaBoundaryCondition')
    rootNode.addObject('RequiredPlugin', name='SofaPython3')
    rootNode.addObject('RequiredPlugin', name='SofaSparseSolver')
    rootNode.addObject('RequiredPlugin', name='SofaOpenglVisual')
    rootNode.addObject('RequiredPlugin', name='SofaConstraint')
    rootNode.addObject('RequiredPlugin', name='SofaLoader')
    rootNode.addObject('RequiredPlugin', name='SofaImplicitOdeSolver')
    rootNode.addObject('RequiredPlugin', name='SofaMeshCollision')
    rootNode.addObject('RequiredPlugin', name='SofaRigid')
    rootNode.addObject('RequiredPlugin', name='CosseratPlugin')
    rootNode.addObject('RequiredPlugin', name='SofaDeformable')
    rootNode.addObject('RequiredPlugin', name='SofaGeneralLinearSolver')
    rootNode.addObject('RequiredPlugin', name='SofaGeneralRigid')

    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels'
                                                   'hideBoundingCollisionModels hideForceFields'
                                                   'hideInteractionForceFields hideWireframe')
    rootNode.addObject('DefaultVisualManagerLoop')
    rootNode.findData('dt').value = DT
    rootNode.findData('gravity').value = [0., 0., -GRAVITY]

    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('GenericConstraintSolver', tolerance=1e-6, maxIterations=1000)

    # ----- Python controller ----- #

    rootNode.addObject(PythonController(name="PythonController"))

    # -------------------------------------#
    # -----  Plastic Cosserat beam   ----- #
    # -------------------------------------#

    # ----- Cosserat beam parameters ----- #

    # Define the total length of the beam
    tot_length = 100.0

    # Define the number of section, the total length and the length of each beam.
    nbSectionS = 1
    lengthS = tot_length / nbSectionS

    # Define the number of frame and the length between each frame.
    nbFramesF = 20
    lengthF = tot_length /nbFramesF

    points = []
    position = []
    lines = []

    position3D = []
    for i in range(nbFramesF):
        sol = i * lengthF
        points.append(i)
        position.append([sol, 0, 0, 0, 0, 0, 1])
        position3D.append([sol, 0, 0])
        if i != nbFramesF-1:
            lines += [i, i+1]

    # ----- Rigid base ----- #

    cableNode = rootNode.addChild('cableNode')
    cableNode.addObject('EulerImplicitSolver', rayleighStiffness="1.2", rayleighMass='1.1')
    cableNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixdouble")
    cableNode.addObject('GenericConstraintCorrection')

    rigidBaseNode= cableNode.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode.addObject('MechanicalObject', template='Rigid3d',
                                          name="RigidBaseMO", position=[0., 0., 0., 0, 0, 0, 1], showObject=1,
                                          showObjectScale=2.)
    rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="5.e8", angularStiffness="5.e8",
                            external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")

    # ----- Rate of angular deformation ----- #
    # Define the length of each beam in a list, the positions of each beam

    positionS = []
    longeurS = []
    sum = 0.
    curv_abs_inputS = []
    curv_abs_inputS.append(0.0)
    for i in range(nbSectionS):
        positionS.append([0, 0, 0])
        longeurS.append(lengthS)
        sum += longeurS[i]
        curv_abs_inputS.append(sum)

    curv_abs_inputS[nbSectionS] = tot_length

    # Define angular rate which is the torsion(x) and bending (y, z) of each section
    rateAngularDeformNode = cableNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject('MechanicalObject',
                                    template='Vec3d', name='rateAngularDeformMO',
                                    position=positionS, showIndices=0, rest_position=[0.0, 0.0, 0.0])
    rateAngularDeformNode.addObject('BeamPlasticLawForceField',
                                    crossSectionShape='circular', length=longeurS,
                                    radius=2., youngModulus=5e6,
                                    initialYieldStress=5e4, plasticModulus=2e5)
    rateAngularDeformNode.addObject('ConstantForceField', name='Moment', indices="0", forces="0 2e8 2e8")

    # ----- Frames ----- #

    # Define local frames related to each section and parameters curv_abs_outputF
    framesF = []
    curv_abs_outputF = []
    for i in range(nbFramesF):
        sol = i * lengthF
        framesF.append([sol, 0, 0,  0, 0, 0, 1])
        curv_abs_outputF.append(sol)

    framesF.append([tot_length, 0, 0, 0, 0, 0, 1])
    curv_abs_outputF.append(tot_length)

    # The node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.addChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)

    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                         name="FramesMO", position=framesF,
                                         showObject=1, showObjectScale=1)

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input= curv_abs_inputS,
                              curv_abs_output=curv_abs_outputF, input1=inputMO, input2=inputMO_rigid,
                              output=outputMO, forcefield='@../../rateAngularDeform/BeamPlasticLawForceField',
                              nonColored=False, debug=0)

    # ----- Mapping ----- #

    cosCollisionPoints = mappedFrameNode.addChild('cosCollisionPoints')
    cosCollisionPoints.addObject('MechanicalObject', name="cosColliPoints", template="Vec3d",
                                 position=position3D)
    cosCollisionPoints.addObject('SkinningMapping', nbRef='2')

    # -------------------------------------#
    # -----  Elastic Cosserat beam   ----- #
    # -------------------------------------#

    # ----- Cosserat beam parameters ----- #

    # Define the total length of the beam
    tot_length = 100.0

    # Define the number of section, the total length and the length of each beam.
    nbSectionS = 1
    lengthS = tot_length / nbSectionS

    # Define the number of frame and the length between each frame.
    nbFramesF = 20
    lengthF = tot_length /nbFramesF

    points = []
    position = []
    lines = []

    position3D = []
    for i in range(nbFramesF):
        sol = i * lengthF
        points.append(i)
        position.append([sol, 0, 1, 0, 0, 0, 1])
        position3D.append([sol, 0, 1])
        if i != nbFramesF-1:
            lines += [i, i+1]

    # ----- Rigid base ----- #

    elaCableNode = rootNode.addChild('elaCableNode')
    elaCableNode.addObject('EulerImplicitSolver', rayleighStiffness="1.2", rayleighMass='1.1')
    elaCableNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixdouble")
    elaCableNode.addObject('GenericConstraintCorrection')

    rigidBaseNode= elaCableNode.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode.addObject('MechanicalObject', template='Rigid3d',
                                          name="RigidBaseMO", position=[0., 0., 10., 0, 0, 0, 1], showObject=1,
                                          showObjectScale=2.)
    rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="5.e8", angularStiffness="5.e8",
                            external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")

    # ----- Rate of angular deformation ----- #
    # Define: the length of each beam in a list, the positions of each beam

    positionS = []
    longeurS = []
    sum = 0.
    curv_abs_inputS = []
    curv_abs_inputS.append(0.0)
    for i in range(nbSectionS):
        positionS.append([0, 0, 0])
        longeurS.append(lengthS)
        sum += longeurS[i]
        curv_abs_inputS.append(sum)

    curv_abs_inputS[nbSectionS] = tot_length

    # Define angular rate which is the torsion(x) and bending (y, z) of each section
    rateAngularDeformNode = elaCableNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject('MechanicalObject',
                                    template='Vec3d', name='rateAngularDeformMO',
                                    position=positionS, showIndices=0, rest_position=[0.0, 0.0, 0.0])
    rateAngularDeformNode.addObject('BeamHookeLawForceField',
                                    crossSectionShape='circular', length=longeurS,
                                    radius=2., youngModulus=5e6)
    rateAngularDeformNode.addObject('ConstantForceField', name='Moment', indices="0", forces="0 2e8 2e8")

    # ----- Frames ----- #

    # Define local frames related to each section and parameters curv_abs_outputF
    framesF = []
    curv_abs_outputF = []
    for i in range(nbFramesF):
        sol = i * lengthF
        framesF.append([sol, 0, 1,  0, 0, 0, 1])
        curv_abs_outputF.append(sol)

    framesF.append([tot_length, 0, 1, 0, 0, 0, 1])
    curv_abs_outputF.append(tot_length)

    # The node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.addChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)

    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                         name="FramesMO", position=framesF,
                                         showObject=1, showObjectScale=1)

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input= curv_abs_inputS,
                                 curv_abs_output=curv_abs_outputF, input1=inputMO, input2=inputMO_rigid,
                                 output=outputMO, debug=0)

    # ----- Mapping ----- #

    cosCollisionPoints = mappedFrameNode.addChild('cosCollisionPoints')
    cosCollisionPoints.addObject('MechanicalObject', name="cosColliPoints", template="Vec3d",
                                 position=position3D)
    cosCollisionPoints.addObject('SkinningMapping', nbRef='2')


    return rootNode
