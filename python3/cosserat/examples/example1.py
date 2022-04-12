# -*- coding: utf-8 -*-
"""
    Cosserat class in SofaPython3.
"""

__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "Jan, 17 2022"

import Sofa
from cosserat.cosseratObject import Cosserat
from cosserat.nonLinearCosserat import NonLinearCosserat as nonCosserat
from cosserat.usefulFunctions import buildEdges, pluginList, BuildCosseratGeometry
from math import sqrt

# @todo ================ Unit: N, m, Kg, Pa  ================
LegendrePolyOrder = 5  # P_{k_2}=P_{k_3}
initialStrain = [[0., 0., 0], [0., 0., 0], [0., 0., 0], [0., 0., 0], [0., 0., 0]]

YM = 1.0e8
PR = 0.
rayleighStiffness = 1.e-3  # Nope
firstOrder = 1

EI = 1.e2
coeff = 1

F1 = [0., 0., 0., 0., (coeff*1.)/sqrt(2), (coeff*1.)/sqrt(2)]  # N

Rb = 0.01  # beam radius in m
length = 1  # in m
nbSection = 5  #
deltaT = 0.02  # s

nonLinearConfig = {'init_pos': [0., 0., 0.], 'tot_length': length, 'nbSectionS': nbSection,
                   'nbFramesF': 15, 'buildCollisionModel': 0, 'beamMass': 0.}


class ForceController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.frames = kwargs['cosseratFrames']
        self.forceNode = kwargs['forceNode']
        self.size = nonLinearConfig['nbFramesF']
        self.applyForce = True
        self.forceCoeff = coeff
        # self.cosseratGeometry = kwargs['cosseratGeometry']

    def onAnimateEndEvent(self, event):
        if self.applyForce:
            with self.forceNode.force.writeable() as force:
                vec = [0., 0., 0., 0., (self.forceCoeff * 1.) / sqrt(2), (self.forceCoeff * 1.) / sqrt(2)]
                for i, v in enumerate(vec):
                    force[i] = v
                # print(f' The new force: {force}')

    def onKeypressedEvent(self, event):
        key = event['key']
        if key == "+":
            self.forceCoeff += 1
            print(f' The new force coeff is : {self.forceCoeff}')
        elif key == "-":
            self.forceCoeff -= 1
            print(f' The new force coeff is : {self.forceCoeff}')


def createScene(rootNode):
    rootNode.addObject('RequiredPlugin', name='plugins', pluginName=[pluginList,
                                                                     ['SofaEngine', 'SofaLoader', 'SofaSimpleFem',
                                                                      'SofaExporter']])
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideCollisionModels '
                                                   'hideBoundingCollisionModels hireForceFields '
                                                   'hideInteractionForceFields hideWireframe')
    rootNode.findData('dt').value = deltaT
    # rootNode.findData('gravity').value = [0., -9.81, 0.]
    rootNode.findData('gravity').value = [0., 0., 0.]
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    # rootNode.addObject('FreeMotionAnimationLoop')
    # rootNode.addObject('GenericConstraintSolver', tolerance=1e-5, maxIterations=5e2)
    rootNode.addObject('Camera', position="-35 0 280", lookAt="0 0 0")

    solverNode = rootNode.addChild('solverNode')
    solverNode.addObject('EulerImplicitSolver', rayleighStiffness=rayleighStiffness, rayleighMass='0.',
                         firstOrder=firstOrder)
    solverNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")
    # solverNode.addObject('SparseLUSolver', name='solver', template="CompressedRowSparseMatrixd")
    # solverNode.addObject('CGLinearSolver', tolerance=1.e-12, iterations=1000, threshold=1.e-18)

    needCollisionModel = 0  # use this if the collision model if the beam will interact with another object
    nonLinearCosserat = solverNode.addChild(
        nonCosserat(parent=solverNode, cosseratGeometry=nonLinearConfig, useCollisionModel=needCollisionModel,
                    name="cosserat", radius=Rb, youngModulus=YM, legendreControlPoints=initialStrain, poissonRatio=PR,
                    order=LegendrePolyOrder, rayleighStiffness=rayleighStiffness))
    cosseratNode = nonLinearCosserat.legendreControlPointsNode
    cosseratNode.addObject('MechanicalMatrixMapper', template='Vec3,Vec3',
                           object1=cosseratNode.getLinkPath(),
                           object2=cosseratNode.getLinkPath(),
                           name='cosseratCoordinateNodeMapper',
                           nodeToParse=nonLinearCosserat.cosseratCoordinateNode.getLinkPath())

    beamFrame = nonLinearCosserat.cosseratFrame

    constForce = beamFrame.addObject('ConstantForceField', name='constForce', showArrowSize=1.e-8,
                        indices=nonLinearConfig['nbFramesF'], force=F1)

    nonLinearCosserat = solverNode.addObject(
        ForceController(parent=solverNode, cosseratFrames=beamFrame.FramesMO, forceNode=constForce))


    # # solverNode2 = rootNode.addChild('solverNode2')
    # # solverNode2.addObject('EulerImplicitSolver', rayleighStiffness="0.2", rayleighMass='0.1')
    # # solverNode2.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")
    # # solverNode2.addObject('GenericConstraintCorrection')
    # # cosserat2 = solverNode2.addChild(Cosserat(parent=solverNode2, cosseratGeometry=linearConfig,
    # #                                           useCollisionModel=needCollisionModel, name="cosserat2", radius=0.1))
    #
    # beamFrame2 = cosserat2.cosseratFrame
    # beamFrame2.addObject('ConstantForceField', name='constForce', showArrowSize=0, indices=12,
    #                      force=[0., 0., 0., 0., 450., 0.])

    return rootNode