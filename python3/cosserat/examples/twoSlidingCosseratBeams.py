# -*- coding: utf-8 -*-
"""Basic scene using Cosserat in SofaPython3.

Based on the work done with SofaPython. See POEMapping.pyscn.
"""

__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com, yinoussa.adagolodjo@inria.fr"
__version__ = "1.0.0"
__copyright__ = "(c) 2020,Inria"
__date__ = "March 8 2020"

import SofaRuntime
import Sofa
import os
from splib3.numerics import Quat
from cosserat.cosseratObject import Cosserat
from cosserat.usefulFunctions import pluginList

# from stlib3.scene import Node
path = os.path.dirname(os.path.abspath(__file__)) + '/../mesh/'


class Animation(Sofa.Core.Controller):
    """
        Implements the AnimationManager as a PythonScriptController
    """

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.rigidBaseMO = args[0]
        self.rateAngularDeformMO = args[1]

        self.rate = 0.2
        self.angularRate = 0.02

    def applyRotation(self, _qOld, posA):
        qNew = Quat.createFromEuler([0., self.angularRate, 0.], 'ryxz')
        qNew.normalize()
        qNew.rotateFromQuat(_qOld)
        for i in range(4):
            posA[0][i + 3] = qNew[i]

    def onKeypressedEvent(self, event):
        key = event['key']
        if ord(key) == 19:  # up
            with self.rigidBaseMO.rest_position.writeable() as posA:
                qOld = Quat()
                for i in range(4):
                    qOld[i] = posA[0][i + 3]
                self.applyRotation(qOld, posA)

        if ord(key) == 21:  # down
            with self.rigidBaseMO.rest_position.writeable() as posA:
                qOld = Quat()
                for i in range(4):
                    qOld[i] = posA[0][i + 3]
                self.applyRotation(qOld, posA)

        if ord(key) == 18:  # left
            with self.rigidBaseMO.rest_position.writeable() as posA:
                print(f' The position at beginning is : {posA}')
                posA[0][0] -= self.rate
                print(f' The position at the end is : {posA}')
        if ord(key) == 20:  # right
            with self.rigidBaseMO.rest_position.writeable() as posA:
                print(f' The position at beginning is : {posA}')
                posA[0][0] += self.rate
                print(f' The position at the end is : {posA}')


def createScene(rootNode):
    rootNode.addObject('RequiredPlugin', name='plugins', pluginName=[pluginList,
                                                                     ['SofaEngine', 'SofaLoader', 'SofaSimpleFem',
                                                                      'SofaExporter']])
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideCollisionModels '
                                                   'hideBoundingCollisionModels hireForceFields '
                                                   'hideInteractionForceFields hideWireframe')
    rootNode.findData('dt').value = 0.01
    rootNode.findData('gravity').value = [0., 0, 0.]
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('GenericConstraintSolver', tolerance=1e-5, maxIterations=5e2)
    rootNode.addObject('Camera', position="-35 0 280", lookAt="0 0 0")

    solverNode = rootNode.addChild('solverNode')
    solverNode.addObject('EulerImplicitSolver', rayleighStiffness="0.2", rayleighMass='0.1')
    solverNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")
    # solverNode.addObject('SparseLUSolver', name='solver')
    solverNode.addObject('GenericConstraintCorrection')

    cosserat_config0 = {'init_pos': [0., 2., 0.], 'tot_length': 20, 'nbSectionS': 20,
                        'nbFramesF': 40, 'beamMass': 0.22}
    cosserat_config1 = {'init_pos': [0., 0, 0.], 'tot_length': 20, 'nbSectionS': 20,
                        'nbFramesF': 40, 'beamMass': 0.22}

    cosserat0 = solverNode.addChild(
        Cosserat(parent=solverNode, cosseratGeometry=cosserat_config0, name="cosseratBeam0", radius=0.2, showObject='1',
                 attachingToLink='1', position=[[0., 2, 0., 0, 0, 0, 1]]))

    # create slidingPoint
    FemPos = [[i * 4, 0., 0.] for i in range(6)]
    femPoints = cosserat0.cosseratFrame.addChild('femPoints')
    inputFEMCable = femPoints.addObject('MechanicalObject', name="pointsInFEM", position=FemPos, showObject="1",
                                        showIndices="1")
    femPoints.addObject('SphereCollisionModel', radius=0.2)
    femPoints.addObject('SkinningMapping', nbRef='1')

    cosserat1 = solverNode.addChild(
        Cosserat(parent=solverNode, cosseratGeometry=cosserat_config1, name="cosseratBeam1", radius=0.2, showObject='1',
                 attachingToLink='1',position=[[0., -2, 0., 0, 0, 0, 1]]))

    cable_position = [[i*2., 2., 0.] for i in range(11)]

    slidingPoint = cosserat1.cosseratFrame.addChild('slidingPoint')
    slidingPointMO = slidingPoint.addObject('MechanicalObject', name="cablePos",
                                            position=cable_position, showObject="1", showIndices="0")
    slidingPoint.addObject('SphereCollisionModel', radius=0.3)
    slidingPoint.addObject('SkinningMapping', nbRef='1')

    mappedPointsNode = slidingPoint.addChild('MappedPoints')
    femPoints.addChild(mappedPointsNode)
    mappedPoints = mappedPointsNode.addObject('MechanicalObject', template='Vec3d', position=cable_position,
                                              name="FramesMO", showObject='1', showObjectScale='1')
    inputCableMO = slidingPointMO.getLinkPath()
    inputFEMCableMO = inputFEMCable.getLinkPath()
    outputPointMO = mappedPoints.getLinkPath()

    controller = solverNode.addObject(
        Animation(cosserat1.rigidBaseNode.RigidBaseMO, cosserat1.cosseratCoordinateNode.cosseratCoordinateMO))

    mappedPointsNode.addObject('QPSlidingConstraint', name="QPConstraint")
    # mappedPointsNode.addObject('CosseratSlidingConstraint', name="QPConstraint")
    mappedPointsNode.addObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO,
                               input2=inputCableMO, output=outputPointMO, direction="@../../FramesMO.position")

    # solverNode.addObject('LinearSolverConstraintCorrection')
    return rootNode
