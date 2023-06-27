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
from splib3.numerics import Quat

LegendrePolyOrder = 5
# initialStrain = [[0., 0., 0], [0., 0., 0], [0., 0., 0]]
initialStrain = [[0., 0., 0], [0., 0., 0], [0., 0., 0], [0., 0., 0], [0., 0., 0]]

YM = 4.015e8
rayleighStiffness = 0.2  # Nope
forceCoeff = 100
F1 = [0., forceCoeff, 0., 0., 0., 0.]  # Nope
Rb = 0.57/2.  # @todo ==> beam radius in m
length = 100.  # @todo ==>  beam length in m
nbSection = 30  # P_{k_2}=P_{k_3}
nbFrames = 60
firstOrder = 1

inertialParams = {'GI': 3.5e7, 'GA': 1.61e8, 'EI': 3.5e7, 'EA': 1.61e8, 'L': length,  'Rb': Rb}
nonLinearConfig = {'init_pos': [0., 0., 0.], 'tot_length': length, 'nbSectionS': nbSection,
                   'nbFramesF': nbFrames, 'buildCollisionModel': 0, 'beamMass': 0.}


class ForceController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.frames = kwargs['cosseratFrames']
        self.forceNode = kwargs['forceNode']
        self.size = nonLinearConfig['nbFramesF']
        self.applyForce = True
        self.forceCoeff = forceCoeff
        # self.cosseratGeometry = kwargs['cosseratGeometry']


        if self.applyForce:
            position = self.frames.position[self.size]  # get the last rigid of the cosserat frame
            orientation = Quat(position[3], position[4], position[5], position[6])  # get the orientation
            # Get the force direction in order to remain orthogonal to the last section of beam
            with self.forceNode.force.writeable() as force:
                vec = orientation.rotate([0., self.forceCoeff, 0.])
                # print(f' The new vec is : {vec}')
                for count in range(3):
                    force[count] = vec[count]
            if self.forceCoeff < 13.1e4:
                self.forceCoeff += 100
            else:
                print(f' The new force coeff is : {self.forceCoeff}')

    def onKeypressedEvent(self, event):
        key = event['key']
        if key == "+":
            self.forceCoeff += 10
            print(f' The new force coeff is : {self.forceCoeff}')
        elif key == "-":
            self.forceCoeff -= 0.2
            print(f' The new force coeff is : {self.forceCoeff}')


def createScene(rootNode):
    rootNode.addObject('RequiredPlugin', name='plugins', pluginName=[pluginList])
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideCollisionModels '
                                                   'hideBoundingCollisionModels hideForceFields '
                                                   'hideInteractionForceFields hideWireframe showMechanicalMappings')
    rootNode.dt.value = 0.02
    rootNode.gravity.value = [0., 0., 0.]
    rootNode.addObject('Camera', position="-35 0 280", lookAt="0 0 0")
    rootNode.addObject('DefaultAnimationLoop')

    solverNode = rootNode.addChild('solverNode')
    solverNode.addObject('EulerImplicitSolver', rayleighStiffness=0, rayleighMass='0.',
                         firstOrder=firstOrder)
    solverNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")

    needCollisionModel = 0  # use this if the collision model if the beam will interact with another object
    PCS_Cosserat = solverNode.addChild(
        Cosserat(parent=solverNode, cosseratGeometry=nonLinearConfig, useCollisionModel=needCollisionModel,
                 inertialParams=inertialParams, name="cosserat", radius=Rb, youngModulus=YM))

    beamFrame = PCS_Cosserat.cosseratFrame

    constForce = beamFrame.addObject('ConstantForceField', name='constForce', showArrowSize=0.0,
                        indices=nonLinearConfig['nbFramesF'], force=F1)

    solverNode.addObject(ForceController(parent=solverNode, cosseratFrames=beamFrame.FramesMO, forceNode=constForce))

    return rootNode
