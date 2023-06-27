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
from splib3.numerics import Quat
from math import pi

# @todo ================ Unit: N, m, Kg, Pa  ================
from splib3.numerics import Quat

LegendrePolyOrder = 5
# initialStrain = [[0., 0., 0], [0., 0., 0], [0., 0., 0]]
initialStrain = [[0., 0., 0], [0., 0., 0], [0., 0., 0], [0., 0., 0], [0., 0., 0]]

YM = 4.015e8
rayleighStiffness = 0.2  # Nope
forceCoeff = 1.
F1 = [0., 0, 0., 0., 0., 0.]  # Nope
Rb = 0.0176/2.  # @todo ==> beam radius in m
length = 20.  # @todo ==> beam length in m
nbSection = 20  # P_{k_2}=P_{k_3}
nbFrames = 60
firstOrder = 1

inertialParams = {'GI': 1800, 'GA': 1e4, 'EI': 1000, 'EA': 5.137e7, 'L': length,  'Rb': Rb}
nonLinearConfig = {'init_pos': [0., 0., 0.], 'tot_length': length, 'nbSectionS': nbSection,
                   'nbFramesF': nbFrames, 'buildCollisionModel': 0, 'beamMass': 0.}


class ForceController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.frames = kwargs['cosseratFrames']
        self.controllerNode = kwargs['controller']
        self.size = nonLinearConfig['nbFramesF']
        self.applyForce = True
        self.forceCoeff = forceCoeff
        self.phase1 = True
        self.phase2 = False
        # self.cosseratGeometry = kwargs['cosseratGeometry']

    def onAnimateEndEvent(self, event):
        if self.applyForce:
            if self.phase1:
                if self.forceCoeff < 2*pi:
                    with self.controllerNode.position.writeable() as position:
                        # vec = Quat()
                        frames = self.frames.position[self.size]
                        vec = Quat(frames[3], frames[4], frames[5], frames[6])  # get the orientation
                        # print(f' The forceCoeff is : {self.forceCoeff}')
                        print(f' The vec before : {vec}')
                        # get the orientation
                        vec.rotateFromEuler([0.01, 0., 0.])
                        vec.normalize()
                        print(f' The vec after : {vec}')
                        for i, v in enumerate(vec):
                            position[0][i+3] = v
                        self.forceCoeff += 0.01
                else:
                    self.phase2 = True

            if self.phase2:
                with self.controllerNode.position.writeable() as position:
                    if position[0][0] > 10.:
                        position[0][0] -= 0.01
                    print(f' The new end effector is : {position[0]}')


    def onKeypressedEvent(self, event):
        key = event['key']
        if key == "+":
            self.phase1 = False
            self.phase2 = True
            print(' The phase 2 is activated !')
        elif key == "-":
            self.phase1 = True
            self.phase2 = False
            print(' The phase 1 is activated !')


def createScene(rootNode):
    rootNode.addObject('RequiredPlugin', name='plugins', pluginName=[pluginList])
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideCollisionModels '
                                                   'hideBoundingCollisionModels hideForceFields '
                                                   'hideInteractionForceFields hideWireframe showMechanicalMappings')
    rootNode.dt.value = 0.02
    rootNode.gravity.value = [0., 0., 0.]
    rootNode.addObject('BackgroundSetting', color='1 1 1')

    solverNode = rootNode.addChild('solverNode')
    solverNode.addObject('EulerImplicitSolver', rayleighStiffness=0, rayleighMass='0.',
                         firstOrder=firstOrder)
    solverNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")

    needCollisionModel = 0  # use this if the collision model if the beam will interact with another object
    nonLinearCosserat = solverNode.addChild(
        nonCosserat(parent=solverNode, cosseratGeometry=nonLinearConfig, useCollisionModel=needCollisionModel,
                    name="cosserat", radius=Rb, youngModulus=YM, legendreControlPoints=initialStrain,
                    order=LegendrePolyOrder, inertialParams=inertialParams, showObject='1',
                    activatedMMM=True))
    beamFrame = nonLinearCosserat.cosseratFrame

    # @todo attach the end effector of the beam with a control point
    EndEffectorColler = rootNode.addChild('EndEffectorController')
    controlMo = EndEffectorColler.addObject('MechanicalObject', template='Rigid3d', name="controlEndEffector",
                            showObjectScale=0.3, position=[length, 0, 0, 0, 0, 0, 1], showObject=True)

    beamFrame.addObject('RestShapeSpringsForceField', name='spring',
                        stiffness=1e8, angularStiffness=1e8, external_points=0,
                        external_rest_shape=controlMo.getLinkPath(), points=nbFrames, template="Rigid3d")

    cosseratNode = nonLinearCosserat.legendreControlPointsNode
   
    constForce = beamFrame.addObject('ConstantForceField', name='constForce', showArrowSize=1.e-9,
                        indices=nonLinearConfig['nbFramesF'], force=F1)

    nonLinearCosserat = solverNode.addObject(
        ForceController(parent=solverNode, cosseratFrames=beamFrame.FramesMO, controller=controlMo))

    return rootNode
