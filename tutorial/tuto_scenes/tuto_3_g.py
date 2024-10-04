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

nb_segment = 32
geoParams = BeamGeometryParameters(init_pos=[0., 0., 0.], beamLength=0.650, show_frames_object=1, nbSection=nb_segment,
                                   nbFrames=nb_segment, buildCollisionModel=0)
physicsParams = BeamPhysicsParameters(beamMass=0.081, youngModulus=0.8e9, poissonRatio=0.35,
                                      beamRadius=2.5e-3, beamLength=0.650)
simuParams = SimulationParameters(rayleighMass=0.0, rayleighStiffness=0.0, firstOrder=False)

Params = Parameters(beamGeoParams=geoParams, beamPhysicsParams=physicsParams, simuParams=simuParams)


def createScene(root_node):
    addHeader(root_node)
    root_node.gravity = [0, -9.810, 0.]

    solver_node = addSolverNode(root_node, name="solver_node")

    # create cosserat Beam
    solver_node.addChild(CosseratBase(parent=solver_node, params=Params))

    return root_node
