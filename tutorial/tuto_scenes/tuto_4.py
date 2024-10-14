# -*- coding: utf-8 -*-
"""
    Cosserat class in SofaPython3.
"""

__authors__ = "adagolodjo"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "October, 26 2021"

from useful.header import addHeader, addSolverNode
from useful.params import BeamPhysicsParametersNoInertia, BeamGeometryParameters, SimulationParameters
from useful.params import Parameters
from cosserat.CosseratBase import CosseratBase
from math import sqrt
from splib3.numerics import Quat
import Sofa
from math import pi

_beam_radius = 0.5
_beam_length = 30.
_nb_section = 32
force_null = [0., 0., 0., 0., 0., 0.]  # N
geoParams = BeamGeometryParameters(beam_length=30.,
                                   nb_section=_nb_section, nb_frames=_nb_section, build_collision_model=0)
physicsParams = BeamPhysicsParametersNoInertia(beam_mass=1., young_modulus=1.0e4, poisson_ratio=0.38,
                                               beam_radius=_beam_radius, beam_length=30)
simuParams = SimulationParameters()
Params = Parameters(beam_geo_params=geoParams, beam_physics_params=physicsParams, simu_params=simuParams)


class ForceController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.forceNode = kwargs['forceNode']
        self.frames = kwargs['frame_node'].FramesMO
        self.force_type = kwargs['force_type']
        self.tip_controller = kwargs['tip_controller']

        self.size = geoParams.nb_frames
        self.applyForce = True
        self.forceCoeff = 0.
        self.theta = 0.1
        self.incremental = 0.01

    def onAnimateEndEvent(self, event):
        if self.applyForce:
            self.forceCoeff += self.incremental
        else:
            self.forceCoeff -= self.incremental

        # choose the type of force 
        if self.force_type == 1:
            # print('inside force type 1')
            self.incremental = 0.1
            self.compute_force()
        elif self.force_type == 2:
            self.incremental = 0.4
            self.compute_orthogonal_force()
        elif self.force_type == 3:
            self.rotate_force()

    def compute_force(self):
        with self.forceNode.forces.writeable() as force:
            vec = [0., 0., 0., 0., self.forceCoeff / sqrt(2), self.forceCoeff / sqrt(2)]
            for i, v in enumerate(vec):
                force[0][i] = v

    def compute_orthogonal_force(self):
        position = self.frames.position[self.size]  # get the last rigid of the cosserat frame
        orientation = Quat(position[3], position[4], position[5], position[6])  # get the orientation
        # Calculate the direction of the force in order to remain orthogonal to the x_axis
        # of the last frame of the beam.
        with self.forceNode.forces.writeable() as force:
            vec = orientation.rotate([0., self.forceCoeff * 5.e-2, 0.])
            # vec.normalize()
            # print(f' The new vec is : {vec}')
            for count in range(3):
                force[0][count] = vec[count]

    def rotate_force(self):
        if self.forceCoeff <= 100. * pi:
            with self.tip_controller.position.writeable() as position:
                last_frame = self.frames.position[self.size]
                vec = Quat(last_frame[3], last_frame[4], last_frame[5], last_frame[6])  # get the orientation

                vec.rotateFromEuler([self.theta, 0., 0.])  # apply rotation around x-axis
                vec.normalize()
                for i, v in enumerate(vec):
                    position[0][i + 3] = v

    def onKeypressedEvent(self, event):
        key = event['key']
        if key == "+":
            self.applyForce = True
        elif key == "-":
            self.applyForce = False


controller_type: int = 1


def createScene(root_node):
    addHeader(root_node)
    root_node.gravity = [0, 0., 0.]

    solver_node = addSolverNode(root_node, name="solver_node")

    # create cosserat Beam
    cosserat_beam = solver_node.addChild(CosseratBase(parent=solver_node, params=Params))
    cosserat_beam.rigidBaseNode.addObject(
        "RestShapeSpringsForceField",
        name="spring",
        stiffness=1e8,
        angularStiffness=1.0e8,
        external_points=0,
        points=0,
        template="Rigid3d"
    )
    cosserat_frames = cosserat_beam.cosseratFrame

    # this constance force is used only in the case we are doing force_type 1 or 2
    const_force_node = cosserat_frames.addObject('ConstantForceField', name='constForce', showArrowSize=1.e-8,
                                                 indices=geoParams.nb_frames, forces=force_null)

    # The effector is used only when force_type is 3
    # create a rigid body to control the end effector of the beam
    tip_controller = root_node.addChild('tip_controller')
    controller_state = tip_controller.addObject('MechanicalObject', template='Rigid3d', name="controlEndEffector",
                                                showObjectScale=0.3, position=[geoParams.beam_length, 0, 0, 0, 0, 0, 1],
                                                showObject=True)
    if controller_type == 3:
        cosserat_frames.addObject('RestShapeSpringsForceField', name='spring', stiffness=0., angularStiffness=1.e8,
                                  external_points=0, external_rest_shape=controller_state.getLinkPath(),
                                  points=geoParams.nb_frames, template="Rigid3d")

    solver_node.addObject(ForceController(forceNode=const_force_node, frame_node=cosserat_frames,
                                          force_type=controller_type, tip_controller=controller_state))

    return root_node
