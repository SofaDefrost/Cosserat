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
import numpy as np

class InsertionController(Sofa.Core.Controller):

    def __init__(self, rootNode, insertionRate, *args, **kwargs):
        # These are needed (and the normal way to override from a python class)
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.rootNode = rootNode
        self.controlPointNode = rootNode.getChild('controlPointNode')
        cosseratBeamNode = rootNode.getChild('cosseratBeamNode')
        self.rigidBaseNode = cosseratBeamNode.getChild('rigidBaseNode')
        self.insertionRate = insertionRate
        self.totalTime = 0.0
        self.nbIterations = 0
        self.autoInsertion = True
        self.backwards = False

    def onAnimateBeginEvent(self, event): # called at each begin of animation step
        if self.autoInsertion:
            with self.controlPointNode.controlPointMO.position.writeable() as controlPointPos:
                if (not self.backwards):
                    controlPointPos[0][0] += self.insertionRate[0]
                    controlPointPos[0][1] += self.insertionRate[1]
                    controlPointPos[0][2] += self.insertionRate[2]
                else:
                    controlPointPos[0][0] -= self.insertionRate[0]
                    controlPointPos[0][1] -= self.insertionRate[1]
                    controlPointPos[0][2] -= self.insertionRate[2]

        self.totalTime = self.totalTime + self.rootNode.findData('dt').value
        self.nbIterations = self.nbIterations + 1

    def onAnimateEndEvent(self, event): # called at each end of animation step
        pass

    def onKeypressedEvent(self, event):
        # Pressing A key activates/deactivates the automatic insertion (in the X axis direction)
        if event['key']=='A':
            if (self.autoInsertion):
                print("Automatic insertion stopped")
            else:
                print("Automatic insertion started")
            self.autoInsertion = not self.autoInsertion

        # Pressing R key inverse the insertion (switch between pushing and pulling)
        if event['key'] == 'C':
            if (self.backwards):
                print("Now pushing the catheter")
            else:
                print("Now pulling the catheter")
            self.backwards = not self.backwards


class InteractiveInsertionController(Sofa.Core.Controller):

    def __init__(self, rootNode, initFrame, *args, **kwargs):
        # These are needed (and the normal way to override from a python class)
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.rootNode = rootNode
        self.controlPointNode = rootNode.getChild('controlPointNode')
        cosseratBeamNode = rootNode.getChild('cosseratBeamNode')
        self.rigidBaseNode = cosseratBeamNode.getChild('rigidBaseNode')
        self.insertionRate = 25
        self.totalTime = 0.0
        self.nbIterations = 0
        self.autoInsertion = True
        self.backwards = False

    # Default Events *********************************************
    def onAnimateBeginEvent(self, event): # called at each begin of animation step
        if self.autoInsertion:
            with self.controlPointNode.controlPointMO.position.writeable() as controlPointPos:
                if (not self.backwards):
                    controlPointPos[0][0] += self.insertionRate[0]
                    controlPointPos[0][1] += self.insertionRate[1]
                    controlPointPos[0][2] += self.insertionRate[2]
                else:
                    controlPointPos[0][0] -= self.insertionRate[0]
                    controlPointPos[0][1] -= self.insertionRate[1]
                    controlPointPos[0][2] -= self.insertionRate[2]

        self.totalTime = self.totalTime + self.rootNode.findData('dt').value
        self.nbIterations = self.nbIterations + 1

    def onAnimateEndEvent(self, event): # called at each end of animation step
        pass

    def onKeypressedEvent(self, event):
        # Pressing A key activates/deactivates the automatic insertion (in the X axis direction)
        if event['key']=='A':
            if (self.autoInsertion):
                print("Automatic insertion stopped")
            else:
                print("Automatic insertion started")
            self.autoInsertion = not self.autoInsertion

        # Pressing R key inverse the insertion (switch between pushing and pulling)
        if event['key'] == 'C':
            if (self.backwards):
                print("Now pushing the catheter")
            else:
                print("Now pulling the catheter")
            self.backwards = not self.backwards


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

    def __init__(self, rootNode, cosseratNode, bendingMoment, fixedIndices, momentAxis, *args, **kwargs):
        # These are needed (and the normal way to override from a python class)
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.rootNode = rootNode
        self.cosseratNode = cosseratNode  # rateAngularDeform
        self.bendingMoment = bendingMoment
        self.fixedIndices = fixedIndices
        self.momentAxis = momentAxis
        self.totalTime = 0.0
        self.nbIterations = 0
        self.bendingActivated = False
        self.fixationActivated = True

    def onAnimateBeginEvent(self, event):  # called at each begin of animation step

        self.totalTime = self.totalTime + self.rootNode.findData('dt').value
        self.nbIterations = self.nbIterations + 1

    def onAnimateEndEvent(self, event):  # called at each end of animation step
        pass

    def onKeypressedEvent(self, event):
        # Pressing B key activates/deactivates the bending of the Cosserat beam
        if event['key'] == 'B':
            self.triggerBendingMoment()
        if event['key'] == 'F':
            self.triggerFixation();

    def triggerBendingMoment(self):
        with self.cosseratNode.bendingMoment.forces.writeable() as moment:
            if (self.bendingActivated):
                # Bending moment is currently activated, it should be deactivated
                print("Bending moment deactivated")
                for i in range(len(moment)):
                    moment[i] = [0.0, 0.0, 0.0]
                self.bendingActivated = False
            else:
                # self.bendingActivated = False => bending moment should be reactivated
                print("Bending moment activated")
                if (self.momentAxis == 'x'):
                    for i in range(len(moment)):
                        moment[i] = [self.bendingMoment, 0.0, 0.0]
                elif (self.momentAxis == 'y'):
                    for i in range(len(moment)):
                        moment[i] = [0.0, self.bendingMoment, 0.0]
                elif (self.momentAxis == 'z'):
                    for i in range(len(moment)):
                        moment[i] = [0.0, 0.0, self.bendingMoment]
                else:
                    print("Only a cardinal axis is expected to apply the bending moment around, please specify "
                          "'x', 'y', or 'z' as a value for parameter momentAxis. By default, the moement is applied "
                          "around the x axis")
                    moment[0] = [self.bendingMoment, 0.0, 0.0]
                self.bendingActivated = True

    def triggerFixation(self):
        # with self.cosseratNode.fixation.indices.value as fixedIndices:
        fixedIndices = self.cosseratNode.fixation.indices.value
        if (self.fixationActivated):
            # Fixation is currently activated, it should be deactivated
            print("Fixation deactivated")
            self.cosseratNode.fixation.indices = []
            self.fixationActivated = False
        else:
            # self.fixationActivated = False => fiwation should be reactivated
            print("Fixation activated")
            self.cosseratNode.fixation.indices = self.fixedIndices
            self.fixationActivated = True