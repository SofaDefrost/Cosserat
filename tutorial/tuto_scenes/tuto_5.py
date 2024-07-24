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
from math import sqrt
from splib3.numerics import Quat
import Sofa
import os
from math import pi
from useful.header import addHeader, addVisual, addSolverNode, Finger
from controller import FingerController
from numpy import array


path = f'{os.path.dirname(os.path.abspath(__file__))}/../../examples/python3/actuators/mesh/'

geoParams = BeamGeometryParameters(init_pos=[0., 0., 0.], beamLength=81., showFramesObject=1,
                                   nbSection=12, nbFrames=32, buildCollisionModel=0)
physicsParams = BeamPhysicsParameters(beamMass=1., youngModulus=5.e6, poissonRatio=0.4, beamRadius=0.5,
                                      beamLength=30)

simuParams = SimulationParameters()
Params = Parameters(beamGeoParams=geoParams, beamPhysicsParams=physicsParams, simuParams=simuParams)

force_null = [0., 0., 0., 0., 0., 0.]  # N
femPos = [[0.0, 0, 0], [15, 0, 0], [30, 0, 0], [45, 0, 0], [60, 0, 0], [66, 0, 0], [81, 0.0, 0.0]]


class ForceController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.forceNode = kwargs['forceNode']
        self.frames = kwargs['frame_node'].FramesMO
        self.force_type = kwargs['force_type']
        self.tip_controller = kwargs['tip_controller']

        self.size = geoParams.nbFrames
        self.applyForce = True
        self.forceCoeff = 0.
        self.theta = 0.1

    def onAnimateEndEvent(self, event):
        if self.applyForce:
            self.forceCoeff += 0.01
        else:
            self.forceCoeff -= 0.01

        # choose the type of force 
        if self.force_type == 1:
            self.compute_force()
        elif self.force_type == 2:
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
        # Calculate the direction of the force in order to remain orthogonal to the x axis of the last frame of the beam.
        with self.forceNode.forces.writeable() as force:
            vec = orientation.rotate([0., self.forceCoeff * 5.e-2, 0.])
            vec.normalize()
            # print(f' The new vec is : {vec}')
            for count in range(3):
                force[0][count] = vec[count]

    def rotate_force(self):
        if self.forceCoeff <= 100. * pi:
            with self.tip_controller.position.writeable() as position:
                last_frame = self.frames.position[self.size]
                vec = Quat(last_frame[3], last_frame[4], last_frame[5], last_frame[6])  # get the orientation

                vec.rotateFromEuler([self.theta, 0., 0.])  # apply rotation arround x-axis 
                vec.normalize()
                for i, v in enumerate(vec):
                    position[0][i + 3] = v

    def onKeypressedEvent(self, event):
        key = event['key']
        if key == "+":
            self.applyForce = True
        elif key == "-":
            self.applyForce = False


def createScene(root_node):
    addHeader(root_node, isConstrained=True)
    root_node.gravity = [0, 0., 0.]

    solver_node = addSolverNode(root_node, name="solver_node", isConstrained=True)

    # create cosserat Beam
    cosserat_beam = solver_node.addChild(CosseratBase(parent=solver_node, params=Params))
    cosserat_frames_node = cosserat_beam.cosseratFrame

    # Finger node
    femFingerNode = root_node.addChild('femFingerNode')
    """ Add FEM finger to the scene"""
    finger_node, fem_points_node = Finger(femFingerNode, name="Finger", rotation=array([0.0, 180.0, 0.0]),
                                       translation=array([-17.5, -12.5, 7.5]), path=path)

    #  This creates a new node in the scene. This node is appended to the finger's node.
    cable_state_node = cosserat_frames_node.addChild('cable_state_node')

    # This creates a MechanicalObject, a componant holding the degree of freedom of our
    # mechanical modelling. In the case of a cable it is a set of positions specifying
    # the points where the cable is passing by.
    cable_state = cable_state_node.addObject('MechanicalObject', name="cablePos",
                                             position=cosserat_beam.frames3D, showObject="1", showIndices="0")
    cable_state_node.addObject('IdentityMapping')

    """ These positions are in fact the distance between fem points and the cable points"""
    distance_node = cable_state_node.addChild('distance_node')
    fem_points_node.addChild(distance_node)
    distance = distance_node.addObject('MechanicalObject', template='Vec3d', position=femPos,
                                           name="distancePointsMO", showObject='1', showObjectScale='1')


    """The controller of the cable is added to the scene"""
    cable_state_node.addObject(FingerController(cosserat_beam.rigidBaseNode.RigidBaseMO,
                                         cosserat_beam.cosseratCoordinateNode.cosseratCoordinateMO))

    inputCableMO = cable_state.getLinkPath()
    inputFEMCableMO = fem_points_node.getLinkPath()
    outputPointMO = distance.getLinkPath()
    """ This constraint is used to compute the distance between the cable and the fem points"""
    distance_node.addObject('QPSlidingConstraint', name="QPConstraint")
    distance_node.addObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO, indices="5",
                                 input2=inputCableMO, output=outputPointMO, direction="@../../FramesMO.position")

    return root_node
    # # this constance force is used only in the case we are doing force_type 1 or 2
    # const_force_node = cosserat_frames.addObject('ConstantForceField', name='constForce', showArrowSize=1.e-8,
    #                                              indices=geoParams.nbFrames, forces=force_null)
    #
    # # The effector is used only when force_type is 3
    # # create a rigid body to control the end effector of the beam
    # tip_controller = root_node.addChild('tip_controller')
    # controller_state = tip_controller.addObject('MechanicalObject', template='Rigid3d', name="controlEndEffector",
    #                                             showObjectScale=0.3, position=[geoParams.beamLength, 0, 0, 0, 0, 0, 1],
    #                                             showObject=True)
    #
    # cosserat_frames.addObject('RestShapeSpringsForceField', name='spring', stiffness=1e8, angularStiffness=1e8,
    #                           external_points=0, external_rest_shape=controller_state.getLinkPath(),
    #                           points=geoParams.nbFrames, template="Rigid3d")
    #
    # solver_node.addObject(ForceController(forceNode=const_force_node, frame_node=cosserat_frames, force_type=3,
    #                                       tip_controller=controller_state))

    return root_node
