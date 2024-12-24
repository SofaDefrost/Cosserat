import Sofa.Core
from splib3.numerics import Quat


class PullingCosseratCable(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.cable = kwargs['frame_node'].FramesMO
        self.tip_controller = kwargs['tip_controller']

    def onAnimateEndEvent(self, event):
        pass

    def pull_cable(self):
        with self.cable.restPosition.writeable() as pos:
            pos[0][0] -= self.rate

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

