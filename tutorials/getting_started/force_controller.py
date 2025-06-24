from math import pi, sqrt

import Sofa
from splib3.numerics import Quat


class ForceController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
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

        # print(f"💪 Applied force {force[:3]} at frame {tip_frame_index}")


    def compute_force(self):
        with self.forceNode.forces.writeable() as force:
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
        position = self.frames.position[
            self.nb_frames
        ]  # get the last rigid of the cosserat frame
        orientation = Quat(
            position[3], position[4], position[5], position[6]
        )  # get the orientation
        # Calculate the direction of the force in order to remain orthogonal to the x_axis
        # of the last frame of the beam.
        with self.forceNode.forces.writeable() as force:
            vec = orientation.rotate([0.0, self.forceCoeff * 5.0e-2, 0.0])
            # vec.normalize()
            # print(f' The new vec is : {vec}')
            for count in range(3):
                force[0][count] = vec[count]
            print(f"💪 Applied force {applied_force[:3]} at frame {tip_frame_index}")

    def rotate_force(self):
        if self.forceCoeff <= 100.0 * pi:
            with self.tip_controller.position.writeable() as position:
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
        key = event["key"]
        if key == "+":
            self.applyForce = True
        elif key == "-":
            self.applyForce = False
