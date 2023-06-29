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
initialStrain = [[0., 0., 0], [0., 0., 0],
                 [0., 0., 0], [0., 0., 0], [0., 0., 0]]

YM = 1.0e8
PR = 0.
rayleighStiffness = 1.e-3  # Nope
firstOrder = 1

EI = 1.e2
coeff = 1

F1 = [0., 0., 0., 0., (coeff*1.)/sqrt(2), (coeff*1.)/sqrt(2)]  # N

# Inertia parameter
Rb = 0.02/2.  # beam radius in m
length = 1.  # in m
nbSection = 5  #
deltaT = 0.02  # s

inertialParams = {'GI': 1.5708, 'GA': 3.1416e4,
                  'EI': 0.7854, 'EA': 3.1416e4, 'L': length,  'Rb': Rb}
nonLinearConfig = {'init_pos': [0., 0., 0.], 'tot_length': length, 'nbSectionS': nbSection,
                   'nbFramesF': 30, 'buildCollisionModel': 0, 'beamMass': 0.}


class ForceController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.frames = kwargs['cosseratFrames']
        self.forceNode = kwargs['forceNode']
        self.size = nonLinearConfig['nbFramesF']
        self.applyForce = True
        self.forceCoeff = coeff

    def onAnimateEndEvent(self, event):
        if self.applyForce:
            with self.forceNode.force.writeable() as force:
                vec = [0., 0., 0., 0., (self.forceCoeff * 1.) /
                       sqrt(2), (self.forceCoeff * 1.) / sqrt(2)]
                for i, v in enumerate(vec):
                    force[i] = v

    def onKeypressedEvent(self, event):
        key = event['key']
        if key == "+":
            self.forceCoeff += 1
            print(f' The new force coeff is : {self.forceCoeff}')
        elif key == "-":
            self.forceCoeff -= 1
            print(f' The new force coeff is : {self.forceCoeff}')


def createScene(rootNode):
    rootNode.addObject('RequiredPlugin', name='plugins', pluginName=[pluginList])
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideCollisionModels '
                                                   'hideBoundingCollisionModels hideForceFields '
                                                   'hideInteractionForceFields hideWireframe showMechanicalMappings')
    rootNode.dt.value = deltaT
    rootNode.gravity.value = [0., 0., 0.]

    rootNode.addObject('Camera', position="-35 0 280", lookAt="0 0 0")
    rootNode.addObject('DefaultAnimationLoop')

    solverNode = rootNode.addChild('solverNode')
    solverNode.addObject('EulerImplicitSolver', rayleighStiffness=rayleighStiffness, rayleighMass='0.',
                         firstOrder=firstOrder)
    solverNode.addObject('CGLinearSolver', tolerance=1.e-12, iterations=1000, threshold=1.e-18)

    # use this if the collision model if the beam will interact with another object
    needCollisionModel = 0
    PCS_Cosserat = solverNode.addChild(
        Cosserat(parent=solverNode, cosseratGeometry=nonLinearConfig, inertialParams=inertialParams, radius=Rb,
                 useCollisionModel=needCollisionModel, name="cosserat", youngModulus=YM, poissonRatio=PR,
                 rayleighStiffness=rayleighStiffness))
    
    beamFrame = PCS_Cosserat.cosseratFrame
    constForce = beamFrame.addObject('ConstantForceField', name='constForce', showArrowSize=1.e-8,
                                     indices=nonLinearConfig['nbFramesF'], force=F1)

    solverNode.addObject(ForceController(
        parent=solverNode, cosseratFrames=beamFrame.FramesMO, forceNode=constForce))

    return rootNode
