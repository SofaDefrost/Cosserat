# -*- coding: utf-8 -*-
"""
    Cosserat class in SofaPython3.
"""

__authors__ = "adagolodjo"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "October, 09 2024"

from useful.header import addHeader, addSolverNode
from useful.params import Parameters
from useful.params import BeamPhysicsParametersNoInertia, BeamGeometryParameters, \
    ContactParameters
from cosserat.CosseratBase import CosseratBase
from softrobots.actuators import PullingCable
import Sofa

beam_length = 1.
geoParams = BeamGeometryParameters(beam_length=beam_length, nb_section=12, nb_frames=12, build_collision_model=0)
physicsParams = BeamPhysicsParametersNoInertia(beam_mass=0.03, young_modulus=1.0e5, poisson_ratio=0.4, beam_radius=0.04,
                                               beam_length=beam_length)
contactParams = ContactParameters()
Params = Parameters(beam_geo_params=geoParams, beam_physics_params=physicsParams, contact_params=contactParams)


def add_mecha_points_with_skinng_maps(node_name, parent_node, positions, _show=True):
    node = parent_node.addChild(node_name)
    meca_node = node.addObject('MechanicalObject', position=positions, template="Vec3d")
    node.addObject('SkinningMapping', nbRef="1", name="skinning_mapping")
    show_mecha_visual(meca_node, show=_show)
    return node


def show_mecha_visual(node, show=True):
    node.showObject = show
    node.showIndices = show


class FingerController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, args, kwargs)
        self.cable = args[0]

    def onKeypressedEvent(self, event):
        displacement = self.cable.CableConstraint.value[0]
        if event["key"] == "+":
            displacement += 1.

        elif event["key"] == "-":
            displacement -= 1.
            if displacement < 0:
                displacement = 0
        self.cable.CableConstraint.value = [displacement]


def createScene(root_node):
    addHeader(root_node, is_constrained=1, is_contact=1, contact_params=Params)
    root_node.gravity = [0, -9.81, 0.]
    solver_node = addSolverNode(root_node, name="solver_node")

    solver_node.addObject('GenericConstraintCorrection')
    # create cosserat Beam
    beam = solver_node.addChild(CosseratBase(parent=solver_node, beam_params=Params))
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

    # add points inside the beam
    frame_node = beam.rigidBaseNode.cosseratInSofaFrameNode
    cable_position = [[0, 0, 0.02], [0.2, 0, 0.02], [0.4, 0, 0.02], [0.6, 0, 0.02], [0.8, 0, 0.02], [1, 0, 0.02]]
    add_mecha_points_with_skinng_maps(node_name='mapped_poind', parent_node=frame_node, positions=cable_position)

    return root_node

    PullingCable(frame_node, cableGeometry=cable_position, name="cable")

    return root_node
