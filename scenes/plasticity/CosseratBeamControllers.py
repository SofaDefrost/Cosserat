# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 2021

@author: ckrewcun
"""

# -*- coding: utf-8 -*-

"""
    Auxiliary methods for Cosserat beam topology generation
"""

import Sofa
import Sofa.constants.Key as Key
import numpy as np
import math
from pyquaternion import Quaternion as Quat


class InsertionController(Sofa.Core.Controller):

    def __init__(self, rootNode, insertionRate, insertionDirection, *args, **kwargs):
        # These are needed (and the normal way to override from a python class)
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.rootNode = rootNode
        self.controlPointNode = rootNode.getChild('controlPointNode')
        cosseratBeamNode = rootNode.getChild('cosseratBeamNode')
        self.rigidBaseNode = cosseratBeamNode.getChild('rigidBaseNode')
        self.insertionRate = insertionRate
        self.insertionDirection = np.array(insertionDirection / np.linalg.norm(insertionDirection))  # normalised
        self.totalTime = 0.0
        self.nbIterations = 0
        self.autoInsertion = True
        self.backwards = False

        # constructs a grid of indices to access only position DoFs of the rigid particle
        self.posDoFsIdGrid = np.ix_([0], [0, 1, 2])

    def onAnimateBeginEvent(self, event):  # called at each begin of animation step
        if self.autoInsertion:
            with self.controlPointNode.controlPointMO.position.writeable() as controlPointPos:
                if not self.backwards:
                    controlPointPos[self.posDoFsIdGrid] += self.insertionDirection * self.insertionRate
                else:
                    controlPointPos[self.posDoFsIdGrid] -= self.insertionDirection * self.insertionRate

        self.totalTime = self.totalTime + self.rootNode.findData('dt').value
        self.nbIterations = self.nbIterations + 1

    def onAnimateEndEvent(self, event):  # called at each end of animation step
        pass

    def onKeypressedEvent(self, event):
        # Pressing A key activates/deactivates the automatic insertion (in the X axis direction)
        if event['key'] == 'A':
            if self.autoInsertion:
                print("Automatic insertion stopped")
            else:
                print("Automatic insertion started")
            self.autoInsertion = not self.autoInsertion

        # Pressing R key inverse the insertion (switch between pushing and pulling)
        if event['key'] == 'C':
            if self.backwards:
                print("Now pushing the catheter")
            else:
                print("Now pulling the catheter")
            self.backwards = not self.backwards


class InteractiveInsertionController(Sofa.Core.Controller):

    def __init__(self, rootNode, initFrameId, insertionRate, insertionDirection, switchFrameDistance,
                 incrementAngle, *args, **kwargs):
        # These are needed (and the normal way to override from a python class)
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.rootNode = rootNode
        self.controlPointNode = rootNode.getChild('controlPointNode')
        cosseratBeamNode = rootNode.getChild('cosseratBeamNode')
        self.rigidBaseNode = cosseratBeamNode.getChild('rigidBaseNode')
        self.mappedFramesNode = self.rigidBaseNode.getChild('MappedFrames')

        self.insertionRate = insertionRate
        self.insertionDirection = np.array(insertionDirection / np.linalg.norm(insertionDirection))  # normalised
        self.currentFrameId = initFrameId
        self.switchFrameDistance = switchFrameDistance
        self.incrementAngle = incrementAngle

        # Computing the incremental quaternions for rotation
        # Taking the direction of insertion as rotation axis
        qw = math.cos(math.radians(self.incrementAngle)/2)
        plusQuat = self.insertionDirection * math.sin(math.radians(self.incrementAngle)/2)
        minusQuat = -plusQuat
        self.plusQuat = Quat(np.insert(plusQuat, 0, qw))
        self.minusQuat = Quat(np.insert(minusQuat, 0, qw))

        self.totalTime = 0.0
        self.nbIterations = 0
        self.nbFrames = len(self.mappedFramesNode.FramesMO.position)
        self.autoInsertion = False
        self.backwards = False

        self.controlPointRefPos = np.array(self.mappedFramesNode.FramesMO.position.value[self.currentFrameId])  # copy

        # constructs a grid of indices to access only position DoFs of the rigid particle
        self.posDoFsIdGrid = np.ix_([0], [0, 1, 2])
        self.posDoFsIdGridSingle = np.ix_([0, 1, 2])
        # constructs a grid of indices to access only orientation DoFs of the rigid particle
        self.quatDoFsIdGrid = np.ix_([0], [3, 4, 5, 6])

        # Placing the control rigid particle on the first frame to 'push'
        with self.controlPointNode.controlPointMO.position.writeable() as controlPointPos:
            with self.mappedFramesNode.FramesMO.position.writeable() as currentFramePos:
                controlPointPos[0] = currentFramePos[self.currentFrameId]
                print("Control point moved to new position : {}".format(controlPointPos[0]))
        # Attaching the RestShapeSpringsForceField to this frame
        with self.mappedFramesNode.controlSpring.points.writeable() as controlSpringPoints:
            controlSpringPoints = self.currentFrameId

    def onAnimateBeginEvent(self, event):  # called at each begin of animation step
        if self.autoInsertion:
            with self.controlPointNode.controlPointMO.position.writeable() as controlPointPos:
                if not self.backwards:
                    controlPointPos[self.posDoFsIdGrid] += self.insertionDirection * self.insertionRate
                else:
                    controlPointPos[self.posDoFsIdGrid] -= self.insertionDirection * self.insertionRate

        self.totalTime = self.totalTime + self.rootNode.findData('dt').value
        self.nbIterations = self.nbIterations + 1

    def onAnimateEndEvent(self, event):  # called at each end of animation step
        pass

    def onKeypressedEvent(self, event):
        # Pressing A key activates/deactivates the automatic insertion (in the X axis direction)
        if event['key'] == 'A':
            if self.autoInsertion:
                print("Automatic insertion stopped")
            else:
                print("Automatic insertion started")
            self.autoInsertion = not self.autoInsertion

        # Pressing R key inverse the insertion (switch between pushing and pulling)
        if event['key'] == 'C':
            if self.backwards:
                print("Now pushing the catheter")
            else:
                print("Now pulling the catheter")
            self.backwards = not self.backwards

        if event['key'] == Key.uparrow:  # Up arrow
            with self.controlPointNode.controlPointMO.position.writeable() as controlPointPos:
                controlPointPos[self.posDoFsIdGrid] += self.insertionDirection * self.insertionRate
                if self.currentFrameId != 0:
                    euclDist = controlPointPos[self.posDoFsIdGrid] - self.controlPointRefPos[self.posDoFsIdGridSingle]
                    if np.linalg.norm(euclDist) > self.switchFrameDistance:
                        # Moving the control point to the previous frame (as we are pushing)
                        self.currentFrameId -= 1
                        # Storing the new reference position by copy
                        self.controlPointRefPos = np.array(self.mappedFramesNode.FramesMO.position.value[self.currentFrameId])
                        # Updating the control point position and changing the frame to which the spring is attached
                        with self.mappedFramesNode.FramesMO.position.writeable() as currentFramePos:
                            controlPointPos[0] = currentFramePos[self.currentFrameId]
                            with self.mappedFramesNode.controlSpring.points.writeable() as attachedFrame:
                                attachedFrame[0] = self.currentFrameId
                # else : we are pushing on the extremity frame, so nothing is to be changed

        if event['key'] == Key.downarrow:  # Down arrow
            with self.controlPointNode.controlPointMO.position.writeable() as controlPointPos:
                controlPointPos[self.posDoFsIdGrid] -= self.insertionDirection * self.insertionRate
                if self.currentFrameId != self.nbFrames-1:
                    euclDist = controlPointPos[self.posDoFsIdGrid] - self.controlPointRefPos[self.posDoFsIdGridSingle]
                    if np.linalg.norm(euclDist) > self.switchFrameDistance:
                        # Moving the control point to the following frame (as we are pulling)
                        self.currentFrameId += 1
                        # Storing the new reference position by copy
                        self.controlPointRefPos = np.array(self.mappedFramesNode.FramesMO.position.value[self.currentFrameId])
                        # Updating the control point position and changing the frame to which the spring is attached
                        with self.mappedFramesNode.FramesMO.position.writeable() as currentFramePos:
                            controlPointPos[0] = currentFramePos[self.currentFrameId]
                            with self.mappedFramesNode.controlSpring.points.writeable() as attachedFrame:
                                attachedFrame[0] = self.currentFrameId
                # else : we are pulling on the extremity frame, so nothing is to be changed

        if event['key'] == Key.leftarrow:  # Left arrow
            # We apply a rotation of self.incrementAngle degrees around the insertion direction.
            # We have to convert the quaternion part of the control point DoFs to a pyquaternion.Quaternion object,
            # as the order of the quaternion coordinates in pyquaternion (w, x, y, z) is not the same as the Rigid3d
            # objects in Sofa (x, y, z, w)
            with self.controlPointNode.controlPointMO.position.writeable() as controlPointPos:
                controlPointSofaQuat = controlPointPos[self.quatDoFsIdGrid]  # np matrix of size 1x4
                qx = controlPointSofaQuat[0][0]
                qy = controlPointSofaQuat[0][1]
                qz = controlPointSofaQuat[0][2]
                qw = controlPointSofaQuat[0][3]
                controlPointQuat = Quat(qw, qx, qy, qz)
                newControlPointQuat = self.minusQuat * controlPointQuat
                controlPointPos[self.quatDoFsIdGrid] = np.array([[newControlPointQuat[1], newControlPointQuat[2],
                                                                newControlPointQuat[3], newControlPointQuat[0]]])

        if event['key'] == Key.rightarrow:  # Right arrow
            # Same as above, but with a quaternion corresponding to the opposite angle rotation
            with self.controlPointNode.controlPointMO.position.writeable() as controlPointPos:
                controlPointSofaQuat = controlPointPos[self.quatDoFsIdGrid]  # np matrix of size 1x4
                qx = controlPointSofaQuat[0][0]
                qy = controlPointSofaQuat[0][1]
                qz = controlPointSofaQuat[0][2]
                qw = controlPointSofaQuat[0][3]
                controlPointQuat = Quat(qw, qx, qy, qz)
                newControlPointQuat = self.plusQuat * controlPointQuat
                controlPointPos[self.quatDoFsIdGrid] = np.array([[newControlPointQuat[1], newControlPointQuat[2],
                                                                  newControlPointQuat[3], newControlPointQuat[0]]])


class ColorMapController(Sofa.Core.Controller):

    def __init__(self, rootNode, springStiffness, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.rootNode = rootNode
        self.springStiffness = springStiffness
        self.anatomyRefNode = self.rootNode.getChild('anatomyRefNode')
        self.anatomyNode = self.rootNode.getChild('anatomyNode')

        self.anatomyDoFs = self.anatomyNode.anatomyDoFs
        self.anatomyRefDoFs = self.anatomyRefNode.anatomyRefDoFs
        self.dataDisplayNode = self.anatomyNode.getChild('dataDisplayNode')
        self.DataDisplay = self.dataDisplayNode.DataDisplay

        with self.anatomyDoFs.position.writeable() as anatomyPos:
            self.nbNodes = len(anatomyPos)

        anatomyForce = []
        for i in range(self.nbNodes):
            anatomyForce += [0.0]
        with self.DataDisplay.pointData.writeable() as pointData:
            pointData = anatomyForce

    def onAnimateEndEvent(self, event):
        anatomyPos = self.anatomyDoFs.position.value
        anatomyRefPos = self.anatomyRefDoFs.position.value

        nodeDistance = []
        forceIntensity = []
        nbNodes = len(anatomyPos)
        for i in range(nbNodes):
            euclDist = np.linalg.norm(anatomyPos[i]-anatomyRefPos[i])
            nodeDistance += [euclDist]
            forceIntensity += [self.springStiffness*euclDist]
        self.DataDisplay.pointData = forceIntensity

class BendingController(Sofa.Core.Controller):

    def __init__(self, rootNode, cosseratNode, frameNode, bendingMoment, fixedIndices, momentAxis, *args, **kwargs):
        # These are needed (and the normal way to override from a python class)
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.rootNode = rootNode
        self.frameNode = frameNode
        self.previousAttachedFrameId = [0]  # not used until bending is activated
        self.cosseratNode = cosseratNode  # rateAngularDeform
        self.bendingMoment = bendingMoment
        self.fixedIndices = fixedIndices
        self.momentAxis = momentAxis
        self.totalTime = 0.0
        self.nbIterations = 0
        self.bendingActivated = False
        self.fixationActivated = False
        self.baseFixationActivated = False

    def onAnimateBeginEvent(self, event):  # called at each begin of animation step

        self.totalTime = self.totalTime + self.rootNode.findData('dt').value
        self.nbIterations = self.nbIterations + 1

    def onAnimateEndEvent(self, event):  # called at each end of animation step
        pass

    def onKeypressedEvent(self, event):
        if event['key'] == 'B':
            self.triggerBendingMoment()
        if event['key'] == 'F':
            self.triggerFixation()
        if event['key'] == 'T':
            self.triggerBaseFixation()

    def triggerBendingMoment(self):
        with self.cosseratNode.bendingMoment.forces.writeable() as moment:
            if self.bendingActivated:
                # Bending moment is currently activated, it should be deactivated
                print("Bending moment deactivated")
                for i in range(len(moment)):
                    moment[i] = [0.0, 0.0, 0.0]
                self.bendingActivated = False

            else:
                # self.bendingActivated = False => bending moment should be reactivated
                print("Bending moment activated")
                if self.momentAxis == 'x':
                    for i in range(len(moment)):
                        moment[i] = [self.bendingMoment, 0.0, 0.0]
                elif self.momentAxis == 'y':
                    for i in range(len(moment)):
                        moment[i] = [0.0, self.bendingMoment, 0.0]
                elif self.momentAxis == 'z':
                    for i in range(len(moment)):
                        moment[i] = [0.0, 0.0, self.bendingMoment]
                else:
                    print("Only a cardinal axis is expected to apply the bending moment around, please specify "
                          "'x', 'y', or 'z' as a value for parameter momentAxis. By default, the moment is applied "
                          "around the x axis")
                    moment[0] = [self.bendingMoment, 0.0, 0.0]
                self.bendingActivated = True

    def triggerFixation(self):

        if self.fixationActivated:
            # Fixation is currently activated, it should be deactivated
            print("Fixation deactivated")
            self.cosseratNode.fixation.indices = []
            self.fixationActivated = False

        else:
            # self.fixationActivated = False => fixation should be reactivated
            print("Fixation activated")
            self.cosseratNode.fixation.indices = self.fixedIndices
            self.fixationActivated = True

    def triggerBaseFixation(self):

        if self.baseFixationActivated:
            # Base fixation is currently activated, it should be deactivated
            print("Base fixation deactivated")
            self.baseFixationActivated = False

            # Moving the control point back to its original position, before bending
            self.frameNode.controlSpring.points = self.previousAttachedFrameId
            controlPointNode = self.rootNode.getChild('controlPointNode')
            with controlPointNode.controlPointMO.position.writeable() as controlPointPos:
                controlPointPos[0] = self.frameNode.FramesMO.position.value[self.previousAttachedFrameId[0]]
        else:
            # self.baseFixationActivated = False => Base fixation should be reactivated
            print("Base fixation activated")
            self.baseFixationActivated = True

            # Moving the control point to the beam base, in order to prevent motion during tip bending
            self.previousAttachedFrameId = np.array(self.frameNode.controlSpring.points.value)  # copy
            self.frameNode.controlSpring.points = [0]
            controlPointNode = self.rootNode.getChild('controlPointNode')
            with controlPointNode.controlPointMO.position.writeable() as controlPointPos:
                controlPointPos[0] = self.frameNode.FramesMO.position.value[0]
