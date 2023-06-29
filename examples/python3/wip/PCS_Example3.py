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
# from cosserat.nonLinearCosserat import NonLinearCosserat as nonCosserat
from cosserat.usefulFunctions import buildEdges, pluginList, BuildCosseratGeometry
from math import sqrt

# @todo ================ Unit: N, m, Kg, Pa  ================
from splib3.numerics import Quat

LegendrePolyOrder = 5
# initialStrain = [[0., 0., 0], [0., 0., 0], [0., 0., 0]]
initialStrain = [[0., 0., 0], [0., 0., 0],
                 [0., 0., 0], [0., 0., 0], [0., 0., 0]]
YM = 4.015e8
rayleighStiffness = 0.2  # Nope
forceCoeff = 0.05
F1 = [0., forceCoeff, 0., 0., 0., 0.]  # Nope
Rb = 0.01/2  # @todo ==> 0.57/2. # beam radius in m
length = 1  # @todo ==>  100 # in m
nbSection = 5  # P_{k_2}=P_{k_3}

nonLinearConfig = {'init_pos': [0., 0., 0.], 'tot_length': length, 'nbSectionS': nbSection,
                   'nbFramesF': 15, 'buildCollisionModel': 0, 'beamMass': 0.}

inertialParams = {'GI': 1800, 'GA': 1e4, 'EI': 1000, 'EA': 5.137e7, 'L': length,  'Rb': Rb}
class ForceController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.frames = kwargs['cosseratFrames']
        self.forceNode = kwargs['forceNode']
        self.size = nonLinearConfig['nbFramesF']
        self.applyForce = True
        self.forceCoeff = forceCoeff
        # self.cosseratGeometry = kwargs['cosseratGeometry']

    def onAnimateEndEvent(self, event):
        if self.applyForce:
            # get the last rigid of the cosserat frame
            position = self.frames.position[self.size]
            # get the orientation
            orientation = Quat(
                position[3], position[4], position[5], position[6])
            # Get the force direction in order to remain orthogonal to the last section of beam
            with self.forceNode.force.writeable() as force:
                vec = orientation.rotate([0., self.forceCoeff, 0.])
                # print(f' The new vec is : {vec}')
                for count in range(3):
                    force[count] = vec[count]

    def onKeypressedEvent(self, event):
        key = event['key']
        if key == "+":
            self.forceCoeff += 0.2
            print(f' The new force coeff is : {self.forceCoeff}')
        elif key == "-":
            self.forceCoeff -= 0.2
            print(f' The new force coeff is : {self.forceCoeff}')


def createScene(rootNode):
    rootNode.addObject('RequiredPlugin', name='plugins', pluginName=[pluginList])
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideCollisionModels '
                                                   'hideBoundingCollisionModels hideForceFields '
                                                   'hideInteractionForceFields hideWireframe')
    rootNode.dt.value = 0.02
    rootNode.gravity.value = [0., 0., 0.]
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.addObject('Camera', position="-35 0 280", lookAt="0 0 0")
    rootNode.addObject('DefaultAnimationLoop')

    solverNode = rootNode.addChild('solverNode')
    solverNode.addObject('EulerImplicitSolver',
                         rayleighStiffness="0.2", rayleighMass='0.')
    solverNode.addObject('EigenSimplicialLDLT', name='solver',
                         template="CompressedRowSparseMatrixMat3x3d")

    # use this if the collision model if the beam will interact with another object
    needCollisionModel = 0
    PCS_Cosserat = solverNode.addChild(
        Cosserat(parent=solverNode, cosseratGeometry=nonLinearConfig, useCollisionModel=needCollisionModel,
                 inertialParams=inertialParams, name="cosserat", radius=Rb, youngModulus=YM))

    beamFrame = nonLinearCosserat.cosseratFrame

    constForce = beamFrame.addObject('ConstantForceField', name='constForce', showArrowSize=0.02,
                                     indices=nonLinearConfig['nbFramesF'], force=F1)

    nonLinearCosserat = solverNode.addObject(
        ForceController(parent=solverNode, cosseratFrames=beamFrame.FramesMO, forceNode=constForce))

    # ----------------

    return rootNode
