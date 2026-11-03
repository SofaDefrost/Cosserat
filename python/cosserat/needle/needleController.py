"""Interactive needle-insertion controller.

Ported from the legacy ``examples/python3/cosserat/needle/needleController.py``
into the maintained ``cosserat`` package. Adapted to the ``CosseratBase`` prefab
API (``rigid_base_node`` / ``cosserat_coordinate_node`` / ``cosserat_frame``) and
to the maintained distance helpers in :mod:`cosserat.utils`.
"""

__authors__ = "adagolodjo"
__contact__ = "adagolodjo@protonmail.com"

import Sofa

from cosserat.needle.params import ConstraintsParams
from cosserat.utils import (
    computePositiveAlongXDistanceBetweenPoints,
    computeNegativeAlongXDistanceBetweenPoints,
)


class Animation(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        # Optional external control MechanicalObject the arrow keys drive; the
        # needle base follows it through a spring (external_rest_shape). Pop it
        # before forwarding kwargs to the SOFA base controller.
        self._control_mo = kwargs.pop("control_mo", None)
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        beam = args[0]
        self.rigidBaseMO = beam.rigid_base_node.RigidBaseMO
        # Data the keyboard moves: the external control MO's position if given,
        # otherwise fall back to driving the base's own rest_position.
        self._drive_data = (
            self._control_mo.position
            if self._control_mo is not None
            else self.rigidBaseMO.rest_position
        )
        # unused, kept for reference
        self.rateAngularDeformMO = beam.cosserat_coordinate_node.cosseratCoordinateMO

        self.needleSlidingState = beam.cosserat_frame.slidingPoint.slidingPointMO
        self.needleCollisionModel = beam.cosserat_frame.needleCollision

        self.contactListener = args[1]
        self.generic = args[2]
        self.entryPoint = []
        self.threshold = 3.0  # lambda_1

        self.constraintPointsNode = args[3]
        self.pointManager = self.constraintPointsNode.pointsManager

        self.constraintPts = self.constraintPointsNode.constraintPointsMo
        self.constraintPtsContainer = self.constraintPointsNode.constraintPtsContainer
        self.constraintPtsModifier = self.constraintPointsNode.constraintPtsModifier
        self.inside = False
        self.addNewPoint = False
        self.rootNode = args[4]

        self.rate = 0.2
        self.angularRate = 0.02
        self.tipForce = [0.0, 0.0, 0]
        return

    def onAnimateEndEvent(self, event):
        if self.contactListener.getContactPoints() and not self.inside:
            # Check if the contact force is large enough to go through the tissue
            if (
                self.generic.constraintForces
                and self.generic.constraintForces[0] > self.threshold
            ):
                # 1. Save the entryPoint and the contact force
                vec = self.contactListener.getContactPoints()[0][1]
                self.entryPoint = [vec[0], vec[1], vec[2]]
                self.tipForce = [
                    self.generic.constraintForces[0],
                    self.generic.constraintForces[1],
                    self.generic.constraintForces[2],
                ]

                # 2. Deactivate the contact constraint
                self.needleCollisionModel.findData("activated").value = 0

                # 3. Add the entry point as the first constraint point in the FEM
                self.pointManager.addNewPointToState()
                self.inside = True

            elif self.tipForce[0] > self.threshold:
                print(
                    "Please activate computeConstraintForces data field inside "
                    "the GenericConstraint component"
                )

        else:
            if self.inside:
                # 1. Check if the needle reached the distance to create/remove a point
                slidingPos = self.needleSlidingState.position.array()
                constraintPos = self.constraintPts.position.array()

                # Add a constraint point when going forwards
                if (
                    computePositiveAlongXDistanceBetweenPoints(slidingPos, constraintPos)
                    > ConstraintsParams.constraintDistance
                ):
                    self.pointManager.addNewPointToState()

                # Going backwards ...
                elif (
                    computeNegativeAlongXDistanceBetweenPoints(slidingPos, constraintPos)
                    > 0
                ):
                    # If last constraint, remove the entry point and get out of the gel
                    if len(constraintPos) == 1:
                        self.inside = False
                        self.needleCollisionModel.findData("activated").value = True
                        self.generic.computeConstraintForces.value = True
                        self.tipForce = [0.0, 0.0, 0]
                        self.pointManager.removeLastPointfromState()

                    # Remove the previous constraint point when going backwards
                    elif (
                        computeNegativeAlongXDistanceBetweenPoints(
                            slidingPos, constraintPos
                        )
                        > ConstraintsParams.constraintDistance
                    ):
                        self.pointManager.removeLastPointfromState()
                else:
                    pass

    def onKeypressedEvent(self, event):
        key = event["key"]
        if key == "M":
            self.pointManager.addNewPointToState()
        if key == "D":
            self.pointManager.removeLastPointfromState()
        if key == "B":
            self.rootNode.findData("animate").value = 1

        if ord(key) == 18:  # left
            with self._drive_data.writeable() as posA:
                posA[0][0] -= self.rate
        elif ord(key) == 20:  # right
            with self._drive_data.writeable() as posA:
                posA[0][0] += self.rate
        elif ord(key) == 21:  # down
            with self._drive_data.writeable() as posA:
                posA[0][1] -= self.rate
        elif ord(key) == 19:  # up
            with self._drive_data.writeable() as posA:
                posA[0][1] += self.rate

        if key == "I":
            with self._drive_data.writeable() as posA:
                posA[0][2] -= self.rate
        elif key == "K":
            with self._drive_data.writeable() as posA:
                posA[0][2] += self.rate
