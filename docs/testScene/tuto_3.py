# -*- coding: utf-8 -*-
"""
    Cosserat class in SofaPython3.
"""

__authors__ = "adagolodjo"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "October, 26 2021"

from useful.header import addHeader, addSolverNode, addVisual
from useful.params import BeamPhysicsParameters, BeamGeometryParameters, SimulationParameters
from useful.params import Parameters
from cosserat.CosseratBase import CosseratBase

geoParams = BeamGeometryParameters(init_pos=[0., 0., 0.], beamLength=30., showFramesObject=1, nbSection=2, nbFrames=12, buildCollisionModel=0)
physicsParams = BeamPhysicsParameters(beamMass=0.3, youngModulus=1.0e3, poissonRatio=0.38, beamRadius=1., beamLength=30)
simuParams = SimulationParameters()
Params = Parameters(beamGeoParams=geoParams, beamPhysicsParams=physicsParams, simuParams=simuParams)


def createScene(root_node):
    addHeader(root_node)
    root_node.gravity = [0, -9.81, 0.]

    solver_node = addSolverNode(root_node, name="solver_node")

    # create cosserat Beam
    solver_node.addChild(CosseratBase(parent=solver_node, params=Params))

    return root_node
