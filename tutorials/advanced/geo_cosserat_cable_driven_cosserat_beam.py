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
import Sofa

# Define the beam parameters
beam_length = 1.
geoParams = BeamGeometryParameters(beam_length=beam_length, nb_section=12, nb_frames=12, build_collision_model=0)
physicsParams = BeamPhysicsParametersNoInertia(beam_mass=0.03, young_modulus=1.0e5, poisson_ratio=0.4, beam_radius=0.04,
                                               beam_length=beam_length)
contactParams = ContactParameters()
Params = Parameters(beam_geo_params=geoParams, beam_physics_params=physicsParams, contact_params=contactParams)

# Define Cable parameters
cable_length = 1.2
cable_geo_params = BeamGeometryParameters(beam_length=cable_length, nb_section=16, nb_frames=16,
                                          build_collision_model=1)
cable_physics_params = BeamPhysicsParametersNoInertia(beam_mass=0.03, young_modulus=1.0e7, poisson_ratio=0.4,
                                                      beam_radius=0.005, beam_length=cable_length)
cable_params = Parameters(beam_geo_params=cable_geo_params, beam_physics_params=cable_physics_params)

def add_mecha_points_with_skinng_maps(node_name, parent_node, positions, _show=True):
    node = parent_node.addChild(node_name)
    meca_node = node.addObject('MechanicalObject', position=positions, template="Vec3d", name=node_name+"_mo")
    node.addObject('SkinningMapping', nbRef="1", name="skinning_mapping")
    return meca_node


def show_mecha_visual(node, show=True):
    node.showObject = show
    node.showIndices = show


class CableController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, args, kwargs)
        self.cable = args[0]
        self.move = True

    def onAnimateEndEvent(self, event):
        if self.move:
            self.pull_cable()

    def pull_cable(self):
        with self.cable.rest_position.writeable() as pos:
            pos[0][0] -= 0.01

    def onKeypressedEvent(self, event):
        if event["key"] == "+":
            self.move = True
        elif event["key"] == "-":
            self.move -= False


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
    sliding_position = [[0, 0, 0.02], [0.2, 0, 0.02], [0.4, 0, 0.02], [0.6, 0, 0.02], [0.8, 0, 0.02], [1, 0, 0.02]]
    sliding_mo = add_mecha_points_with_skinng_maps(node_name='mapped_poind', parent_node=frame_node, positions=sliding_position)

    # Add cable to the scene
    cable_solver = addSolverNode(root_node, name="cable_solver")
    cable_node = cable_solver.addChild(CosseratBase(parent=cable_solver, translation=[-0.2, 0, 0.02],
                                                    beam_params=cable_params, name="cable"))
    cable_solver.addObject('GenericConstraintCorrection')
    cable_frames_node = cable_node.cosseratFrame

    #  This creates a new node in the scene. This node is appended to the finger's node.
    cable_state_node = cable_frames_node.addChild('cable_state_node')

    # This creates a MechanicalObject, a component holding the degree of freedom of our
    # mechanical modelling. In the case of a cable it is a set of positions specifying
    # the points where the cable is passing by.
    cable_state = cable_state_node.addObject('MechanicalObject', name="cablePos",
                                             position=cable_node.frames3D)
    show_mecha_visual(cable_state, show=True)
    cable_state_node.addObject('IdentityMapping')

    """ These positions are in fact the distance between the 3D points mapped inside the Beam and the cable points"""
    distance_node = cable_state_node.addChild('distance_node')
    beam.cosseratFrame.addChild(distance_node)

    distance = distance_node.addObject('MechanicalObject', template='Vec3d', position=sliding_position,
                                       name="distancePointsMO", showObject='1', showObjectScale='1')

    """The controller of the cable is added to the scene"""
    # cable_state_node.addObject(CableController(cable_node.rigidBaseNode.RigidBaseMO))

    inputCableMO = cable_state.getLinkPath()
    sliding_points = sliding_mo.getLinkPath()
    outputPointMO = distance.getLinkPath()
    """ This constraint is used to compute the distance between the cable and the fem points"""
    distance_node.addObject('QPSlidingConstraint', name="QPConstraint")
    distance_node.addObject('DifferenceMultiMapping', name="pointsMulti", input1=sliding_points, indices="5",
                            input2=inputCableMO, output=outputPointMO, direction="@../../FramesMO.position")

    return root_node