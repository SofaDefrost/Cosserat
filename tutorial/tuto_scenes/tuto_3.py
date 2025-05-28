# -*- coding: utf-8 -*-
"""
    Cosserat class in SofaPython3.
"""

__authors__ = "adagolodjo"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "October, 26 2021"

from cosserat.CosseratBase import CosseratBase
from useful.header import addHeader, addSolverNode
from useful.params import (BeamGeometryParameters,
                           BeamPhysicsParametersNoInertia, Parameters)

geoParams = BeamGeometryParameters(beam_length=30.,  nb_section=32, nb_frames=32, build_collision_model=0)
physicsParams = BeamPhysicsParametersNoInertia(beam_mass=0.3, young_modulus=1.0e3, poisson_ratio=0.38, beam_radius=1.,
                                               beam_length=30)
Params = Parameters(beam_geo_params=geoParams, beam_physics_params=physicsParams)


def createScene(root_node):
    addHeader(root_node)
    root_node.gravity = [0, -9.81, 0.]
    solver_node = addSolverNode(root_node, name="solver_node")

    # create cosserat Beam
    beam = solver_node.addChild(CosseratBase(parent=solver_node, params=Params))
    # Attach beam base using a spring force field
    beam.rigidBaseNode.addObject(
            "RestShapeSpringsForceField",
            name="spring",
            stiffness=1e8,
            angularStiffness=1.0e8,
            external_points=0,
            points=0,
            template="Rigid3d"
        )

    return root_node
