# -*- coding: utf-8 -*-
"""
    Cosserat class in SofaPython3.
"""

__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "October, 26 2021"

cosserat_config = {'init_pos': [0., 0., 0.], 'tot_length': 1, 'nbSectionS': 3,
                   'nbFramesF': 12, 'buildCollisionModel': 1, 'beamMass': 0.22}

import Sofa
from cosserat.cosseratObject import Cosserat
from cosserat.usefulFunctions import buildEdges, pluginList, BuildCosseratGeometry
from math import sqrt

# @todo ================ Unit: N, m, Kg, Pa  ================
#LegendrePolyOrder = 3
#initialStrain = [[0., 0., 0], [0., 0., 0], [0., 0., 0]]
coeff = 0.3
YM = 1.0e8
PR = 0.
rayleighStiffness = 0.2  # Nope
F1 = [0., 0., 0., 0., (coeff*1.)/sqrt(2), (coeff*1.)/sqrt(2)]  # Nope
Rb = 0.01/2. # beam radius in m
length = 1  # in m
nbSection = 5  # P_{k_2}=P_{k_3}
deltaT = 0.02  # s
nbFrames = 30
beamMass = 0

beamConfig = {'init_pos': [0., 0., 0.], 'tot_length': length, 'nbSectionS': nbSection,
                   'nbFramesF': nbFrames, 'buildCollisionModel': 0, 'beamMass': beamMass}


def createScene(rootNode):
    rootNode.addObject('RequiredPlugin', name='plugins', pluginName=[pluginList,
                                                                     ['SofaEngine', 'SofaLoader', 'SofaSimpleFem',
                                                                      'SofaExporter']])
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideCollisionModels '
                                                   'hideBoundingCollisionModels hireForceFields '
                                                   'hideInteractionForceFields hideWireframe')
    rootNode.findData('dt').value = deltaT
    rootNode.findData('gravity').value = [0., 0, 0.]
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('GenericConstraintSolver', tolerance=1e-5, maxIterations=5e2)
    rootNode.addObject('Camera', position="-35 0 280", lookAt="0 0 0")

    solverNode = rootNode.addChild('solverNode')
    solverNode.addObject('EulerImplicitSolver', rayleighStiffness=rayleighStiffness, rayleighMass='0.', firstOrder=1)
    solverNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")
    solverNode.addObject('GenericConstraintCorrection')

    needCollisionModel = 0  # use this if the collision model if the beam will interact with another object
    cosserat = solverNode.addChild(
        Cosserat(parent=solverNode, cosseratGeometry=beamConfig, useCollisionModel=beamConfig['buildCollisionModel'],
                 name="cosserat", radius=Rb, youngModulus=YM, poissonRatio=PR))

    # attach force at the beam tip,
    # we can attach this force to non mechanical node thanks to the MechanicalMatrixMapper component
    beamFrame = cosserat.cosseratFrame
    beamFrame.addObject('ConstantForceField', name='constForce', showArrowSize=1.e-2, indices=nbFrames,
                        force=F1)

