__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com, yinoussa.adagolodjo@inria.fr"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "March 8 2021"

import Sofa
from splib3.numerics import Quat
import numpy as np
from params import ConstraintsParams
from cosserat.utils import computeDistanceBetweenPoints
import params


class Animation(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.rigidBaseMO = args[0].rigidBaseNode.RigidBaseMO
        self.rateAngularDeformMO = args[0].cosseratCoordinateNode.cosseratCoordinateMO

        self.needleSlidingState = args[0].cosseratFrame.slidingPoint.slidingPointMO
        self.needleCollisionModel = args[0].cosseratFrame.needleCollision

        self.contactListener = args[1]
        self.generic = args[2]
        self.entryPoint = []
        self.threshold = 3.

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
        self.tipForce = [0., 0., 0]

        return

    def onAnimateEndEvent(self, event):
        if self.contactListener.getContactPoints() and not self.inside:
            vec = self.contactListener.getContactPoints()[0][1]
            tip = [vec[0], vec[1], vec[2]]
            print(f'the contact force is : {self.generic.constraintForces[0]}')

            # @infi: check if the contact force is large enough to go through the tissue
            if self.generic.constraintForces and self.generic.constraintForces[0] > self.threshold:
                # @info 1. Save the entryPoint and contact force
                self.entryPoint = tip
                self.tipForce = [self.generic.constraintForces[0], self.generic.constraintForces[1],
                                 self.generic.constraintForces[2]]

                # @info 2. deactivate the contact constraint
                self.needleCollisionModel.findData('activated').value = 0

                # @info 3. Add entryPoint point as the first constraint point in FEM
                # self.constraintPtsModifier.addPoints(1, True)
                # self.addNewPoint = True
                self.inside = True
                print(
                    "=====> inside the volume, the first point is added to the state <=====")
                print("@Todo: call the binding of addNewPointToState()")

            elif self.tipForce[0] > self.threshold:
                print(
                    "Please activate computeConstraintForces data field inside the GenericConstraint component")
            # 5. todo: If the user is pulling out the needle and the needle tip is behind the entryPoint,
            # 5. todo: call the remove constraint function to remove those constraint points
            # todo: this can be handle with :
            # todo: 5.1 self.inside=False
            # todo: 5.2 self.needleCollisionModel.findData('activated').value = 1
        else:
            if self.inside:
                # @info 1. check if the needle is we reach the distance to create a new constraint point
                if computeDistanceBetweenPoints(self.constraintPts, self.needleSlidingState) > \
                        params.ConstraintsParams.constraintDistance:
                    self.addNewPoint = True
                    # Todo: call the binding of addNewPointToState()
                    # print("=====> Todo: call the binding of addNewPointToState() <=====")
                else:
                    pass
                    # print("=====> The distance between the tip and the last constraint is < constraint Distance <=====")

    def addConstraintPointToFem(self, point):
        with self.constraintPts.position.writeable() as pos:
            sz = len(pos)
            if sz != 0:
                print(f'2:  ====> The state size is : {sz}')
                print(f' ====> The tip is : {point}')
                print(f'1 ====> The state : {pos}')
                # pos[sz - 1][0] = point[0]
                # pos[sz - 1][1] = point[1]
                # pos[sz - 1][2] = point[2]

                print(f'2 ====> The state : {pos}')
                self.addNewPoint = False
            else:
                print(
                    "The state is empty, add a point in modifier before adding a point in the state")

    def onKeypressedEvent(self, event):
        key = event['key']
        if key == "A":  # -
            self.pointManager.addNewPoint
        if key == "D":  # +
            self.constraintPtsModifier.removePoints(np.array([0]), True)
        if key == "B":  # +
            self.rootNode.findData('animate').value = 1

        if ord(key) == 18:  # left
            print("======= > left")
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][0] -= self.rate
        elif ord(key) == 20:  # right
            print("======= > right")
            # print(
            #     f' ====> contactListener : {self.contactListener.getContactPoints()}')
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][0] += self.rate

        elif ord(key) == 21:  # down
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][1] -= self.rate

        elif ord(key) == 19:  # up
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][1] += self.rate

    # def onAnimateBeginEvent(self, event):
    #     if self.inside & self.addNewPoint:
    #         print("Inside the volume")
    #         self.addConstraintPointToFem(self.entryPoint)
    #
    #         # 4. todo: add new constraint points inside the volume if needed.
    #         # todo: This depends on the choice of the algorithm
    #         #  expl1: one can compare tip position related to the last constraint position inside the volume and
    #         #  when this > than the constraintDistance we add new constraint point
    #         # self.addPointToFem()
    #         # expl2: one can compare the tip position related to the entry point and when
    #         # this > than the constraintDistance
    #
    #
