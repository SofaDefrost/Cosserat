# -*- coding: utf-8 -*-
"""Basic scene using Cosserat in SofaPython3.

Based on the work done with SofaPython. See POEMapping.py
"""

__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com, yinoussa.adagolodjo@inria.fr"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "March 16 2021"

import Sofa
import numpy as np

def BuildCosseratGeometry(nbSection=8, nbFrames=12, totalLength=80):
    # Define: the number of section, the total length and the length of each beam.

    lengthS = totalLength / nbSection

    # Define: the length of each beam in a list, the positions of each beam
    # (flexion, torsion), the abs of each section
    positionS = []
    longeurS = []
    temp = 0.
    curv_abs_inputS = [0.0]
    for i in range(nbSection):
        positionS.append([0, 0, 0])
        longeurS.append((((i + 1) * lengthS) - i * lengthS))
        temp += longeurS[i]
        curv_abs_inputS.append(temp)
    curv_abs_inputS[nbSection] = totalLength

    # ###############################################
    # Define: the number of frame and the length between each frame.

    lengthF = totalLength / nbFrames

    # Define: the abs of each frame and the position of each frame.
    framesF = []
    curv_abs_outputF = []
    cable_positionF = []
    for i in range(nbFrames):
        sol = i * lengthF
        framesF.append([sol, 0, 0, 0, 0, 0, 1])
        cable_positionF.append([sol, 0, 0])
        curv_abs_outputF.append(sol)

    framesF.append([totalLength, 0, 0, 0, 0, 0, 1])
    cable_positionF.append([totalLength, 0, 0])
    curv_abs_outputF.append(totalLength)

    return [positionS, curv_abs_inputS, longeurS, framesF, curv_abs_outputF, cable_positionF]


class MoveTargetProcess(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.targetMO = args[0]
        self.targetIndex = 0
        self.targetList = [[85.0, 0,  0.35857], [85.0, 0.5,  0.35857], [86.0, 1.5,  0.35857], [87.0, 2.5,  0.35857]]
        pos = self.targetMO.findData('effectorGoal').value
        print("effectorGoal :", pos[0])


    def onKeypressedEvent(self, event):
        key = event['key']
        if key == "+":  # +
            # with self.targetMO.findData('effectorGoal').writeable() as posR:
            #     print("effectorGoal :", posR)
            #     posR[0][1] = 0.5

            with self.targetMO.position.writeable() as posA:
                print("0. Position target : ", posA)
                if self.targetIndex < len(self.targetList) - 1:
                    self.targetIndex += 1
                print("Target index : ", self.targetIndex)
                for i in range(3):
                    print("I :", i)
                    print("target is :", self.targetList[self.targetIndex][i])
                    posA[0][1] = self.targetList[self.targetIndex][1]
                print("1. Position target : ", posA)

        if key == "-":  # -
            with self.targetMO.rest_position.writeable() as posR:
                with self.targetMO.position.writeable() as posA:
                    if self.targetIndex > 0:
                        self.targetIndex -= 1
                    print("Target index : ", self.targetIndex)

                    for i in range(3):
                        print("I :", i)
                        posA[0][i] = self.targetList[self.targetIndex][i]
                        posR[0][i] = self.targetList[self.targetIndex][i]
                    print("Target index : ", self.targetIndex)
                    print("Position target : ", posA)


class AddPointProcess(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.inputConstraintPointsMO = args[0]
        self.needleCollisionMO = args[1]
        self.diffMapping = args[2]
        self.tipIndex = len(self.needleCollisionMO.position) - 1
        print("tip indice :", self.tipIndex)

    def onKeypressedEvent(self, event):
        key = event['key']
        if key == "+":  # +
            with self.inputConstraintPointsMO.position.writeable() as posA:
                print("0. List The position :", posA)

            self.inputConstraintPointsMO.position.value = np.array([[40, 0, 0], [45, 0, 0]])

            with self.inputConstraintPointsMO.position.writeable() as posA:
                print("1. List The position :", posA)
                # tip = self.needleCollisionMO.position[self.tipIndex]
                # print("Tip pose :", tip)
                # posA = [[40.0, 0., 0.], tip]
                # print("1. List The position :", posA)
                # posA.append(self.needleCollisionMO.position[self.tipIndex])

            indice = self.diffMapping.findData('indices')
            print("====> indices :", self.diffMapping.addPointProcess([40.0, 0., 0.]))

        if key == "-":  # -
            with self.inputConstraintPointsMO.rest_position.writeable() as posA:
                sz = len(posA)
                print("Tip pose :", posA[sz-1])

