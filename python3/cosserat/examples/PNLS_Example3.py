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
forceCoeff = 0.
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
                if self.forceCoeff < pi:
                    with self.controllerNode.position.writeable() as position:
                        print(f' The forceCoeff is : {self.forceCoeff}')
                        print(f' The control point is : {position[0]}')

                        # orientation = Quat(position[0][3], position[0][4], position[0][5], position[0][6])  # get the orientation
                        vec = Quat()  # get the orientation
                        vec.rotateFromEuler([self.forceCoeff, 0., 0.])
                        vec.normalize()
                        print(f' The control vec is : {vec}')
                        for i, v in enumerate(vec):
                            position[0][i+3] = v
                        self.forceCoeff += 0.1
                else:
                    self.phase2 = True

            if self.phase2:
                with self.controllerNode.position.writeable() as position:
                    if position[0][0] > 10.:
                        position[0][0] -= 0.1
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
    rootNode.addObject('RequiredPlugin', name='plugins', pluginName=[pluginList,
                                                                     ['SofaEngine', 'SofaLoader', 'SofaSimpleFem',
                                                                      'SofaExporter']])
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideCollisionModels '
                                                   'hideBoundingCollisionModels hireForceFields '
                                                   'hideInteractionForceFields hideWireframe showMechanicalMappings')
    rootNode.findData('dt').value = 0.02
    rootNode.findData('gravity').value = [0., -9.81, 0.]
    # rootNode.findData('gravity').value = [0., 0., 0.]
    rootNode.addObject('BackgroundSetting', color='1 1 1')
    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('GenericConstraintSolver', tolerance=1e-5, maxIterations=5e2)
    # rootNode.addObject('Camera', position="-35 0 280", lookAt="0 0 0")

    solverNode = rootNode.addChild('solverNode')
    # solverNode.addObject('EulerImplicitSolver', rayleighStiffness="0.", rayleighMass='0.')
    solverNode.addObject('EulerImplicitSolver', rayleighStiffness=0, rayleighMass='0.',
                         firstOrder=firstOrder)
    solverNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")
    # solverNode.addObject('EigenSimplicialLDLT', name='solver', template="CompressedRowSparseMatrixMat3x3d" )

    # solverNode.addObject('CGLinearSolver', tolerance=1.e-12, iterations=1000, threshold=1.e-18)

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
                            showObjectScale=0.2, position=[length, 0, 0, 0, 0, 0, 1])

    beamFrame.addObject('RestShapeSpringsForceField', name='spring',
                        stiffness=1e10, angularStiffness=1e10, external_points=0,
                        external_rest_shape=controlMo.getLinkPath(), points=nbFrames, template="Rigid3d")

    cosseratNode = nonLinearCosserat.legendreControlPointsNode
    cosseratNode.addObject('MechanicalMatrixMapper', template='Vec3,Vec3',
                           object1=cosseratNode.getLinkPath(),
                           object2=cosseratNode.getLinkPath(),
                           name='cosseratCoordinateNodeMapper',
                           nodeToParse=nonLinearCosserat.cosseratCoordinateNode.getLinkPath())

    constForce = beamFrame.addObject('ConstantForceField', name='constForce', showArrowSize=1.e-9,
                        indices=nonLinearConfig['nbFramesF'], force=F1)

    nonLinearCosserat = solverNode.addObject(
        ForceController(parent=solverNode, cosseratFrames=beamFrame.FramesMO, controller=controlMo))

    return rootNode