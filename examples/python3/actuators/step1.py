# -*- coding: utf-8 -*-
"""
    Scene robot ISIR.
"""
__authors__ = "Yinoussa"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "July, 20 2023"

from useful.header import addHeader, addSolverNode
from useful.params import BeamPhysicsParameters, BeamGeometryParameters, SimulationParameters
from useful.params import Parameters
from cosserat.CosseratBase import CosseratBase
from cable import PullingCable
from utils import createRigidDisk, create_cable_points, FingerController
from splib3.loaders import loadPointListFromFile

import os

"""  @info
    To construct this robot, our initial step involves building its central stem. 
    In this regard, we've opted for employing the Cosserat model for constructing this rod. Our approach entails 
    dividing this rod into 14 identical segments, each of which terminates with a rigid element. 
    Subsequently, we will affix a disc to each of these rigid elements.
    Initially, the beam is affected by gravity, oriented along the -y axis. 
    Later on, we will adjust the direction of gravity to match the real-world scenario where the robot is oriented 
    downward, aligning with the gravitational force.
    To simplify scene construction, we've also activated constraints from the beginning using the 'isConstrained' 
    keyword set to 'true'.
"""

# @todo
geoParams = BeamGeometryParameters(init_pos=[0., 0., 0.], beamLength=65.5, show_frames_object=1,
                                   nbSection=14, nbFrames=14, buildCollisionModel=0)
physicsParams = BeamPhysicsParameters(beamMass=5., youngModulus=1.0e8, poissonRatio=0.38, beamRadius=6.2e-1 / 2.,
                                      beamLength=65.5)
simuParams = SimulationParameters()
Params = Parameters(beamGeoParams=geoParams, beamPhysicsParams=physicsParams, simuParams=simuParams)

rotation = [0.0, 0.0, 0.0]
translation = [0.0, 0.0, 0.0]
pullPointLocation = [0.0, 0.0, 0.0]
youngModulus = 18000
valueType = 'position'


def createScene(rootNode):
    addHeader(rootNode, isConstrained=True)
    rootNode.gravity = [0, -9.81, 0.]

    stemNode = addSolverNode(rootNode, name="stemNode", isConstrained=True)

    stemNode.addChild(CosseratBase(parent=stemNode, params=Params))

    return rootNode
