from math import pi, sqrt

import Sofa
from splib3.numerics import Quat


class ForceController(Sofa.Core.Controller):
    """
    A controller to apply different types of forces to the tip of the Cosserat beam.
    This controller is called at each simulation step (onAnimateEndEvent).
    """
    def __init__(self, *args, **kwargs):
        """
        Initializes the controller.

        Args:
            forceNode: The node containing the ConstantForceField to be controlled.
            frame_node: The node containing the beam's frames (MechanicalObject).
            force_type: An integer (1, 2, or 3) specifying the type of force to apply.
            tip_controller: A node with a MechanicalObject used to control the beam's tip (for force_type 3).
            geoParams: The BeamGeometryParameters object for the beam.
        """
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.forceNode = kwargs["forceNode"]
        self.frames = kwargs["frame_node"].FramesMO
        self.force_type = kwargs["force_type"]
        self.tip_controller = kwargs["tip_controller"]
        self.geoParams = kwargs["geoParams"]

        self.nb_frames = self.geoParams.nb_frames
        self.applyForce = True
        self.forceCoeff = 0.0
        self.theta = 0.1
        self.incremental = 0.01

    def onAnimateEndEvent(self, event):
        """
        Called at the end of each animation step.
        This method updates the force coefficient and calls the appropriate force function.
        """
        if self.applyForce:
            self.forceCoeff += self.incremental
        else:
            self.forceCoeff -= self.incremental

        # choose the type of force
        if self.force_type == 1:
            # Applies a constant torque.
            self.incremental = 0.1
            self.compute_force()
        elif self.force_type == 2:
            # Applies a force that is always orthogonal to the beam's tip.
            self.incremental = 0.4
            self.compute_orthogonal_force()
        elif self.force_type == 3:
            # Rotates a target frame that the beam's tip will follow.
            self.rotate_force()

    def compute_force(self):
        """Applies a constant torque to the beam's tip."""
        with self.forceNode.forces.writeable() as force:
            # This vector represents a torque around the Y and Z axes.
            vec = [
                0.0,
                0.0,
                0.0,
                0.0,
                self.forceCoeff / sqrt(2),
                self.forceCoeff / sqrt(2),
            ]
            for i, v in enumerate(vec):
                force[0][i] = v

    def compute_orthogonal_force(self):
        """Applies a force that is always orthogonal to the beam's tip."""
        position = self.frames.position[
            self.nb_frames
        ]  # get the last rigid of the cosserat frame
        orientation = Quat(
            position[3], position[4], position[5], position[6]
        )  # get the orientation
        # Calculate the direction of the force in order to remain orthogonal to the x_axis
        # of the last frame of the beam.
        with self.forceNode.forces.writeable() as force:
            # Rotate a vector [0, y, 0] by the tip's orientation to get the force in world coordinates.
            vec = orientation.rotate([0.0, self.forceCoeff * 5.0e-2, 0.0])
            for count in range(3):
                force[0][count] = vec[count]

    def rotate_force(self):
        """Rotates the target frame (tip_controller) that the beam's tip is constrained to follow."""
        if self.forceCoeff <= 100.0 * pi:
            with self.tip_controller.position.writeable() as position:
                # Get the orientation of the beam's tip
                last_frame = self.frames.position[self.nb_frames]
                vec = Quat(
                    last_frame[3], last_frame[4], last_frame[5], last_frame[6]
                )  # get the orientation

                vec.rotateFromEuler(
                    [self.theta, 0.0, 0.0]
                )  # apply rotation around x-axis
                vec.normalize()
                for i, v in enumerate(vec):
                    position[0][i + 3] = v

    def onKeypressedEvent(self, event):
        """
        Handles key presses to enable or disable the force application.
        Press '+' to enable/increase the force.
        Press '-' to disable/decrease the force.
        """
        key = event["key"]
        if key == "+":
            self.applyForce = True
        elif key == "-":
            self.applyForce = False
