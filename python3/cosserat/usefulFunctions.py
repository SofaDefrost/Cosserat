# -*- coding: utf-8 -*-
"""Basic scene using Cosserat in SofaPython3.

Based on the work done with SofaPython. See POEMapping.py
"""

__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com, yinoussa.adagolodjo@inria.fr"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "March 16 2021"


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