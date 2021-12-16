# -*- coding: utf-8 -*-
"""
    Cosserat class in SofaPython3.
"""

__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "October, 26 2021"

from dataclasses import dataclass
import Sofa
from usefulFunctions import buildEdges, pluginList, BuildCosseratGeometry
from cosserat.cosseratObject import Cosserat

cosserat_config = {'init_pos': [0., 0., 0.], 'tot_length': 1, 'nbSectionS': 3,
                   'nbFramesF': 12, 'buildCollisionModel': 1, 'beamMass': 0.22}


class Animation(Sofa.Core.Controller):

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.bending = kwargs.get('bending', None)
        return

    def onKeypressedEvent(self, event):
        key = event['key']
        # ######## bending #########
        if key == "P":  # +
            with self.bending.position.writeable() as bend:
                print(f'{bend}')


initialStrain = [[0., 0., -0.55157735], [0., 0., -0.38899746], [0., 0., -0.23071231], [0., 0., -0.07628332]]


def createScene(rootNode):
    rootNode.addObject('RequiredPlugin', name='plugins', pluginName=[pluginList,
                                                                     ['SofaEngine', 'SofaLoader', 'SofaSimpleFem',
                                                                      'SofaExporter']])
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideCollisionModels '
                                                   'hideBoundingCollisionModels hireForceFields '
                                                   'hideInteractionForceFields hideWireframe')
    rootNode.findData('dt').value = 0.01
    rootNode.findData('gravity').value = [0., -9.81, 0.]
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('GenericConstraintSolver', tolerance=1e-5, maxIterations=5e2)
    rootNode.addObject('Camera', position="-35 0 280", lookAt="0 0 0")

    solverNode = rootNode.addChild('solverNode')
    solverNode.addObject('EulerImplicitSolver', rayleighStiffness="0.2", rayleighMass='0.1')
    solverNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")
    solverNode.addObject('GenericConstraintCorrection')

    needCollisionModel = 0  # use this if the collision model if the beam will interact with another object
    cosserat = solverNode.addChild(
        Cosserat(parent=solverNode, cosseratGeometry=cosserat_config, useCollisionModel=needCollisionModel,
                 name="cosserat", radius=0.1))

    beamFrame = cosserat.cosseratFrame
    beamFrame.addObject('ConstantForceField', name='constForce', showArrowSize=1.e-2, indices=12,
                        force=[0., -100., 0., 0., 0., 0.])


# def createScene(rootNode):
#     rootNode.addObject('RequiredPlugin', name='plugins', pluginName=[pluginList,
#                                                                      ['SofaEngine', 'SofaLoader', 'SofaSimpleFem',
#                                                                       'SofaExporter']])
#     rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideCollisionModels '
#                                                    'hideBoundingCollisionModels hireForceFields '
#                                                    'hideInteractionForceFields hideWireframe')
#     rootNode.findData('dt').value = 0.01
#     rootNode.findData('gravity').value = [0., -9.81, 0.]
#     rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
#     rootNode.addObject('FreeMotionAnimationLoop')
#     rootNode.addObject('GenericConstraintSolver', tolerance=1e-5, maxIterations=5e2)
#     rootNode.addObject('Camera', position="-35 0 280", lookAt="0 0 0")
#
#     solverNode = rootNode.addChild('solverNode')
#     solverNode.addObject('EulerImplicitSolver', rayleighStiffness="0.2", rayleighMass='0.1')
#     solverNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")
#     solverNode.addObject('GenericConstraintCorrection')
#
#     needCollisionModel = 0  # use this if the collision model if the beam will interact with another object
#     cosserat = solverNode.addChild(
#         Cosserat(parent=solverNode, cosseratGeometry=cosserat_config, useCollisionModel=needCollisionModel,
#                  name="cosserat", radius=0.1))
#
#     cosseratStrain = solverNode.addChild("cosseratStrain")
#     cosseratStrain.addObject("MechanicalObject", position=initialStrain)
#     # attach force at the beam tip,
#     # we can attach this force to non mechanical node thanks to the MechanicalMatrixMapper component
#     beamFrame = cosserat.cosseratFrame
#     # beamFrame.addObject('ConstantForceField', name='constForce', showArrowSize=1.e-2, indices=12,
#     #                     force=[0., -50., 0., 0., 0., 0.])
#
    ################################
    # Animation (to move the dofs) #
    ################################
    bendingVector = cosserat.cosseratCoordinateNode.cosseratCoordinateMO
    rootNode.addObject(Animation(bending=bendingVector))
#
#     return rootNode
