__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com, yinoussa.adagolodjo@inria.fr"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "March 8 2021"

import Sofa

import numpy as np
import Sofa
import Cosserat
from cosserat.needle.params import ConstraintsParams
from cosserat.utils import computePositiveAlongXDistanceBetweenPoints, computeNegativeAlongXDistanceBetweenPoints


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
            # @Info: check if the contact force is large enough to go through the tissue
            if self.generic.constraintForces and self.generic.constraintForces[0] > self.threshold:
                # 1. @Info: Save the entryPoint and contact force
                vec = self.contactListener.getContactPoints()[0][1]
                self.entryPoint = [vec[0], vec[1], vec[2]]
                self.tipForce = [self.generic.constraintForces[0], self.generic.constraintForces[1],
                                 self.generic.constraintForces[2]]

                # 2. @Info: deactivate the contact constraint
                self.needleCollisionModel.findData('activated').value = 0

                # @Info 3. Add entryPoint point as the first constraint point in FEM
                self.pointManager.addNewPointToState()
                self.inside = True

            elif self.tipForce[0] > self.threshold:
                print("Please activate computeConstraintForces data field inside the GenericConstraint component")

        else:
            if self.inside:
                # @info 1. check if the needle reached the distance to create/remove a constraint point
                slidingPos = self.needleSlidingState.position.array()
                constraintPos = self.constraintPts.position.array()
                
                # Add constraint point when going forwards
                if computePositiveAlongXDistanceBetweenPoints(slidingPos, constraintPos) > ConstraintsParams.constraintDistance:
                    self.pointManager.addNewPointToState()

                # Going backwards ...
                elif computeNegativeAlongXDistanceBetweenPoints(slidingPos, constraintPos) > 0:

                    # If last constraint, remove the entry point and get out from the gel
                    if len(constraintPos) == 1 :
                        self.inside=False
                        self.needleCollisionModel.findData('activated').value = True
                        self.generic.computeConstraintForces.value = True
                        self.tipForce = [0., 0., 0]
                        self.pointManager.removeLastPointfromState()

                    # Remove previous constraint point when going backwards
                    elif computeNegativeAlongXDistanceBetweenPoints(slidingPos, constraintPos) > ConstraintsParams.constraintDistance:
                        self.pointManager.removeLastPointfromState()
                else:
                    pass

    def onKeypressedEvent(self, event):
        key = event['key']
        if key == "M": 
            self.pointManager.addNewPointToState()
        if key == "D": 
            self.pointManager.removeLastPointfromState()
        if key == "B":  
            self.rootNode.findData('animate').value = 1

        if ord(key) == 18:  # left
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][0] -= self.rate
        elif ord(key) == 20:  # right
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][0] += self.rate

        elif ord(key) == 21:  # down
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][1] -= self.rate

        elif ord(key) == 19:  # up
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][1] += self.rate
        
        if key == "I":  # 
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][2] -= self.rate
        elif key == "K":  # 
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][2] += self.rate
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
