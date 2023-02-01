#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__authors__ = "Yinoussa:Younes"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021, Inria"
__date__ = "September 2 2021"

"""
    Create the scene with the
    Units: mm, kg, s.

    In this version(2), the goal is to have a rigid between the two or more Cosserat beams.
    Means, instead of having a constraint between two Cosserat beams we have two
    [RigidDistanceMapping].
    {Cosserat1} <-> [RigidDistanceMapping] {RigidElement} [RigidDistanceMapping] <-> {Cosserat2}
"""

import Sofa
from splib3.numerics import Vec3, Quat
from math import pi
import numpy as np

# constants
GRAVITY = 0
TOT_MASS = 0.1
YOUNG_MODULUS = 7e6
DENSITY = 0.02
TOTAL_LENGTH = 4
NB_FRAMES = 6
NB_SECTIONS = 4
cosserat_config1 = {'init_pos': [0., 0., 0.], 'tot_length': TOTAL_LENGTH, 'nbSectionS': NB_SECTIONS,
                    'nbFramesF': NB_FRAMES}
cosserat_config2 = {'init_pos': [0., 0., 0.], 'tot_length': TOTAL_LENGTH, 'nbSectionS': NB_SECTIONS,
                    'nbFramesF': NB_FRAMES}

POISSON_RATIO = 0.4
LENGTH_Y = 1.
LENGTH_Z = 1.
config_material = {'youngModulus': YOUNG_MODULUS, 'poissonRatio': POISSON_RATIO, 'length_Y': LENGTH_Y,
                   'length_Z': LENGTH_Z, 'shape': 'circular'}

is_constrained = True
is_CG = 0

pluginNameList = 'CosseratPlugin' \
                 ' Sofa.Component.AnimationLoop' \
                 ' Sofa.Component.Constraint.Lagrangian.Correction' \
                 ' Sofa.Component.Constraint.Lagrangian.Solver' \
                 ' Sofa.Component.LinearSolver.Direct' \
                 ' Sofa.Component.Mapping.MappedMatrix' \
                 ' Sofa.Component.ODESolver.Backward' \
                 ' Sofa.Component.SolidMechanics.Spring' \
                 ' Sofa.Component.StateContainer' \
                 ' Sofa.Component.Visual'


class Animation(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.rigidBaseMO = args[0].RigidBaseMO
        self.rate = 0.1

    def onKeypressedEvent(self, event):
        key = event['key']
        if ord(key) == 18:  # left
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][0] -= self.rate
        elif ord(key) == 20:  # right
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][0] += self.rate
        elif ord(key) == 19:  # up
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][1] += self.rate
        elif ord(key) == 21:  # down
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][1] -= self.rate
        elif key == "i" or key == "I":  # page up
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][2] += self.rate
        elif key == "k" or key == "K":  # page down
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][2] -= self.rate
        elif key == "+":  # rotate aroud x
            angle = Vec3(1, 0, 0)
            with self.rigidBaseMO.rest_position.writeable() as posA:
                self._extracted_from_onKeypressedEvent_17(posA, angle)
        elif key == "-":  # rotate aroud x
            angle = Vec3(0, 1, 0)
            with self.rigidBaseMO.rest_position.writeable() as posA:
                self._extracted_from_onKeypressedEvent_17(posA, angle)
                # TODO Rename this here and in `onKeypressedEvent`

    def _extracted_from_onKeypressedEvent_17(self, posA, angle):
        # QInit = Quat.createFromEuler([0.5, 0, 0], 'ryxz', inDegree=True)
        QInit = Quat.createFromAxisAngle(angle, pi / 2.)
        QInit.normalize()
        qTemp = Quat()

        for i in range(4):
            qTemp[i] = posA[0][3 + i]

        rot = QInit.rotateFromQuat(qTemp)
        rot.normalize()
        print(f'rot: {rot}')
        posA[0][3:7] = rot


def createScene(rootNode):
    rootNode.addObject('RequiredPlugin', pluginName=pluginNameList, printLog='0')
    rootNode.addObject('VisualStyle',
                       displayFlags='showVisualModels showBehaviorModels showCollisionModels '
                                    'hideBoundingCollisionModels hideForceFields hideInteractionForceFields '
                                    'hideWireframe')
    rootNode.findData('dt').value = 0.01
    rootNode.findData('gravity').value = [0., 0., -GRAVITY]
    if is_constrained:
        rootNode.addObject('FreeMotionAnimationLoop')
        rootNode.addObject('GenericConstraintSolver', tolerance="1e-20", maxIterations="500", printLog="0")

    solverNode = rootNode.addChild('solverNode')
    solverNode.addObject('EulerImplicitSolver', firstOrder=0, rayleighStiffness=0.01, rayleighMass=0.)
    if is_CG:
        solverNode.addObject('CGLinearSolver', name='solver')
    else:
        solverNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")

    if is_constrained:
        solverNode.addObject('GenericConstraintCorrection')

    [x, y, z] = cosserat_config2['init_pos']
    # orientation = [0.7071, 0, 0.7071,0]
    orientation = [0, 0, 0, 1]
    rigidBaseNode1 = solverNode.addChild('rigidBase1')
    RigidBaseMO = rigidBaseNode1.addObject('MechanicalObject', template='Rigid3d',
                                           name="RigidBaseMO", position=[x, y, z] + orientation, showObject=1,
                                           showObjectScale=0.5, rest_position=[0., 0., 0., 0., 0., 0., 1.])
    rigidBaseNode1.addObject('RestShapeSpringsForceField', name='spring1', stiffness=2e1, angularStiffness=1e1,
                             external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")

    # distanceNode.addObject('RestShapeSpringsForceField', name='spring2', stiffness="1e5", angularStiffness="1e5",
    #                      external_points="0", mstate="@FramesMO", points="0", template="Rigid3d")

    rigidBaseNode2 = solverNode.addChild('rigidBase2')
    rigidBaseNode2.addObject('MechanicalObject', name="rigidState", template='Rigid3d', position="2 0. 0 0 0 0 1",
                             showObject=True, showObjectScale=1.0)
    # position="2 4 1 0.7071 0 0.7071 0",

    rigidBaseNode2.addObject('RestShapeSpringsForceField', name='rigid3', stiffness=0.2e1, template="Rigid3d",
                             angularStiffness=0.2e1, external_points=0, mstate="@rigidState", points=0)

    # rigidBase2Child = rigidBaseNode2.addChild('rigidBase2Child')
    # interRigidChild.addChild(distanceNode)  # add double sub_node
    # rigidBase2Child.addObject('MechanicalObject', name="interRigidChildMo", template='Rigid3d',
    #                          position="1. 0.7 0 0 0 0 1", showObject=True, showObjectScale=0.2)


    distanceNode = rigidBaseNode1.addChild('distanceNode')
    distanceMo = distanceNode.addObject('MechanicalObject', template='Rigid3d', name="distanceMo",
                                        position=[1., 0., 0., 0., 0., 0., 1.], showObject=1, showObjectScale=0.1)
    distanceMapping = True
    if distanceMapping:
        distanceNode.addObject('RigidDistanceMapping', input1=rigidBaseNode2.getLinkPath(),
                               input2=rigidBaseNode1.getLinkPath(), newVersionOfFrameComputation='1',
                               output=distanceMo.getLinkPath(), first_point=[0], second_point=[0], name='distanceMap1')
        solverNode.addObject('MechanicalMatrixMapper', template='Rigid3d,Rigid3d', object1=rigidBaseNode1.getLinkPath(),
                             object2=rigidBaseNode2.getLinkPath(), name='mapper1',
                             nodeToParse=distanceNode.getLinkPath())
        distanceNode.addObject('CosseratNeedleSlidingConstraint', name='constraintMappingConstraint',
                               template="Rigid3d", useDirections=np.array([0, 1, 1, 0, 0, 0]))
    else:
        distanceNode.addObject('RigidRigidMapping', name="interRigidMap", globalToLocalCoords='0', template="Rigid3d,Rigid3d")



    solverNode.addObject(Animation(rigidBaseNode1))
    return rootNode