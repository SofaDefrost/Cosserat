#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Basic scene using Cosserat in SofaPython3.

Created on Tue Feb  2 16:21:56 2021

@author: younesssss
"""

__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com, yinoussa.adagolodjo@inria.fr"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "March 16 2021"

import Sofa
import numpy as np
import math
from splib3.numerics import Vec3, Quat

# from stlib3.scene import MainScene

lenSoftPart = 0.5
rateLength = [1.2, 0.5, 0.7, 0.5, 0.7, 0.5, 0.7, 0.5, 0.7, 0.5, 0.7, 0.5, 0.7, 0.5, 0.7, 0.5, 0.7, 0.5,
              0.7, 0.5, 0.7, 0.5, 0.7, 0.5, 0.7, 0.5, 0.7, 0.5, 0.7, 0.5, 0.7, 0.5, 0.7, 0.5, 0.7, 0.5,
              0.7, 0.5, 0.7, 0.5, 0.45]

draw_cylinder = 1
add_collision_point = 1

pluginList = ['Cosserat',
              'Sofa.Component.AnimationLoop',
              'Sofa.Component.LinearSolver.Iterative', 'Sofa.Component.LinearSolver.Direct', 'Sofa.Component.StateContainer',
              'Sofa.Component.SolidMechanics.FEM.Elastic', 'Sofa.Component.SolidMechanics.Spring', 'Sofa.Component.MechanicalLoad', 'Sofa.Component.Mass',
              'Sofa.Component.Topology.Container.Dynamic', 'Sofa.Component.Topology.Container.Grid',
              'Sofa.Component.ODESolver.Backward',
              'Sofa.Component.Mapping.Linear', 'Sofa.Component.Mapping.NonLinear', 'Sofa.Component.Topology.Mapping',
              'Sofa.Component.Collision.Geometry', 'Sofa.Component.Collision.Detection.Algorithm','Sofa.Component.Collision.Detection.Intersection','Sofa.Component.Collision.Response.Contact',
              'Sofa.Component.Constraint.Lagrangian.Correction','Sofa.Component.Constraint.Lagrangian.Solver',
              'Sofa.Component.Engine.Select',
              'Sofa.GL.Component.Rendering3D', 'Sofa.Component.Setting', 'Sofa.Component.SceneUtility', 'Sofa.Component.IO.Mesh', 'Sofa.Component.Visual'
              ]


PRig = 0.38  # poison ratio for the rigid part of the beam
PSoft = 0.43  # poison ratio for the soft part of the beam
poisonRatioList = [PSoft, PRig, PSoft, PRig, PSoft, PRig, PSoft, PRig, PSoft, PRig, PSoft, PRig, PSoft, PRig, PSoft,
                   PRig, PSoft, PRig, PSoft, PRig, PSoft, PRig, PSoft, PRig, PSoft, PRig, PSoft, PRig, PSoft, PRig,
                   PSoft, PRig, PSoft, PRig, PSoft, PRig, PSoft, PRig, PSoft, PRig, PSoft]
# YRig = 215.e11
YRig = 215.e9

temp = [8.1e4, YRig, 3.1e4, YRig, 2.5e4, YRig, 2.3e4, YRig, 2.4e4, YRig, 2.3e4, YRig, 1.95e4, YRig,
        1.93e4, YRig, 1.52e4, YRig, 1.03e4, YRig, 1.14e4, YRig, 1.07e4, YRig, 7.9e3, YRig, 7.3e3,
        YRig, 6.0e3, YRig, 5.0e3, YRig, 2.975e3, YRig, 1.77e3, YRig, 1.615e3, YRig, 9.25e2, YRig,
        9.25e2]

youngModulusList = [x * 4 for x in temp]

adds = 1.4e4
# youngModulusList = [8.1e4, YRig, 3.1e4, YRig, 2.5e4, YRig, 2.3e4, YRig, 2.4e4, YRig, 2.3e4, YRig, 1.95e4, YRig,
#                    2.e4, YRig, 2.3e4, YRig, 2.2e4, YRig, 1.14e4+adds, YRig, 1.07e4+adds, YRig, 7.9e3+adds, YRig,
#                    7.3e3+adds, YRig, 6.0e3+adds, YRig, 5.0e3+5.e3, YRig, 2.975e3, YRig, 1.615e3, YRig, 9.25e2,
#                    YRig, 9.25e2, YRig, 9.25e2]

"""
 The plasticity parameter table.
The first soft beam has no parameter, because it was not tested in the experimental part. But we set it to 0.36 like beam 1.
The rigid part has no parameter, but we set it to 100, in order to simplify the code
"""
rPlastic = 100.
plasticityTab = [0.36, rPlastic, 0.36, rPlastic, 0.35, rPlastic, 0.41, rPlastic, 0.34, rPlastic, 0.42, rPlastic, 0.48,
                 rPlastic, 0.45, rPlastic, 0.54, rPlastic, 0.17, rPlastic, 0.20, rPlastic, 0.2, rPlastic, 0.15,
                 rPlastic,
                 0.15, rPlastic, 0.12, rPlastic, 0.2, rPlastic, 0.3, rPlastic, 0.3, rPlastic, 0.3, rPlastic, 0.3,
                 rPlastic, 0.3]


def formated(value):
    return "%.8f" % value


def buildEdges(cable3DPos):
    """ This function is used to build edges required in the EdgeSetTopologyContainer component"""
    points_size = len(cable3DPos)
    edgeList = []
    for i in range(0, points_size - 1):
        edgeList.append(i)
        edgeList.append(i + 1)
    return edgeList


class CosseratComponent(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.rigidBaseNode = args[0]
        self.rateAngularDeformNode = args[1]
        self.frameNode = args[0].getChild('MappedFrames')
        self.cable3DPosNode = self.frameNode.getChild(
            'CollisInstrumentCombined')
        self.rigidBaseMO = self.rigidBaseNode.getChild('RigidBaseMO')
        self.dX = 0.0025
        self.angularRate = 15.
        self.length = 12.95
        # Cosserat frames parameters
        self.nbFrames = 55
        self.frames = []
        self.curvOut = []
        self.mappedPos = []
        self.cable3DPos = []
        self.edgeList = []
        # angular stress parameters, this also defines the numbers of section
        self.nbSection = 22
        self.rateAngular = []

        leng = 0.
        for i in range(len(rateLength)):
            leng = rateLength[i] + leng
        self.length = leng

        self.curvInput = [0.0]
        temp = 0.0
        for i in range(len(rateLength)):
            temp += rateLength[i]
            self.curvInput.append(temp)
            self.rateAngular.append([0, 0., 0.])
        counter = 0

        for i in range(len(self.curvInput) - 1):
            middle = (self.curvInput[i + 1] + self.curvInput[i]) / 2.
            self.frames.append([self.curvInput[i], 0., 0., 0., 0., 0., 1.])
            self.frames.append([middle, 0., 0., 0., 0., 0., 1.])

            self.curvOut.append(self.curvInput[i])
            self.curvOut.append(middle)

            self.cable3DPos.append([self.curvInput[i], 0., 0.])
            self.cable3DPos.append([middle, 0., 0.])
            counter += 2

        self.curvOut.append(self.length)
        self.cable3DPos.append([self.length, 0., 0.])
        # print("self.cable3DPos: ", self.cable3DPos)
        self.frames.append([self.length, 0., 0., 0., 0., 0., 1.])
        counter += 1
        self.nbFrames = counter
        self.initGraph()
        print("# ====> End CosseratComponent.__init__(self, *args, **kwargs)")

    def applyRotation(self, angle=None):
        if angle is None:
            angle = [0., 0., 0.]
        QInit = Quat.createFromEuler(angle, 'ryxz', inDegree=True)
        QInit.normalize()
        return QInit.rotate([self.dX, 0., 0.])

    def computeCosseratParameters(self):
        points_size = len(self.cable3DPos)
        for i in range(points_size - 1):
            self.edgeList.append(i)
            self.edgeList.append(i + 1)

    def createCylinderNode(self, length=0.5, abscissa=0.0, index=0):
        translation = [abscissa, 0., 0.]
        color = 'grey' if length == 0.5 else 'blue'
        CylinderFEMNode = self.frameNode.addChild(
            'CylinderFEMNode' + str(index))
        CylinderFEMNode.addObject('CylinderGridTopology', name="grid", nx="5", ny="5", nz="8", length=length,
                                  radius=0.25, axis="1 0 0")
        CylinderFEMNode.addObject('MeshTopology', src="@grid")
        CylinderFEMNode.addObject(
            'MechanicalObject', name="cylinder", template="Vec3d", translation=translation)
        CylinderFEMNode.addObject('SkinningMapping', nbRef='1')
        Visu = CylinderFEMNode.addChild('Visu')
        Visu.addObject('OglModel', name="Visual", color=color)
        Visu.addObject('IdentityMapping',
                       input="@../cylinder", output="@Visual")

    def initGraph(self):
        # @info compute cosserat parameters
        self.computeCosseratParameters()
        self.frameNode.FramesMO.position.value = self.frames
        self.frameNode.mapping.curv_abs_output.value = self.curvOut
        self.frameNode.mapping.curv_abs_input.value = self.curvInput
        self.rateAngularDeformNode.rateAngularDeformMO.position.value = self.rateAngular
        self.rateAngularDeformNode.beamHookeLaw.findData(
            'length').value = rateLength
        self.rateAngularDeformNode.beamHookeLaw.findData(
            'youngModululsList').value = youngModulusList
        self.rateAngularDeformNode.beamHookeLaw.findData(
            'poissonRatioList').value = poisonRatioList

        if draw_cylinder == 1:
            # Draw the cylinder over the beam
            rate = 0
            for i in range(len(rateLength)):
                self.createCylinderNode(rateLength[i], rate, i)
                rate += rateLength[i]

        if add_collision_point == 1:
            # creat collision points for the beam
            edgeList = self.cable3DPosNode.collisEdgeSet
            with self.cable3DPosNode.collisEdgeSet.edges.writeable() as edges:
                edges = self.edgeList
                print("self.edgeList :", self.edgeList)
            self.cable3DPosNode.collisEdgeSet.position.value = self.cable3DPos

    # def onKeypressedEvent(self, event):
    #     c = event['key']
    #     print("The computeCosseratParameters :", c)
    #     if c == 'L':  # up
    #         with self.frameNode.FramesMO.position.writeable() as posA:
    #             print("onKeypressedEvent: ==============> posA :", posA[7])
    #     ######### Rate angular #########
    #     if ord(c) == 19:  # up
    #         with self.rigidBaseMO.rest_position.writeable() as posA:
    #             posA[0][0] += self.ratePos[0]  #self.rate
    #             posA[0][1] += self.ratePos[1]  # self.rate
    #             posA[0][2] += self.ratePos[2]  # self.rate
    #             # //self.rigidBaseMO.findData('rest_position').value = posA
    #
    #     if ord(c) == 21:  # down
    #         with self.rigidBaseMO.rest_position.writeable() as posA:
    #             posA[0][0] -= self.ratePos[0]  # self.rate
    #             posA[0][1] -= self.ratePos[1]  # self.rate
    #             posA[0][2] -= self.ratePos[2]  # self.rate
    # self.rigidBaseMO.findData('rest_position').value = posA


k = 1e5


class forceControl(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.node = args[0]
        self.refNode = self.node.getChild('PlaneRef')
        self.cochleaNode = self.node.getChild('cochleaNode')

        self.MecaObject = self.cochleaNode.dofs
        self.RefObject = self.refNode.Ref
        self.colorMapNode = self.cochleaNode.getChild('ColorMap')
        self.colorMap = self.colorMapNode.ColorMap
        self.DataDisplay = self.colorMapNode.DataDisplay

        with self.MecaObject.position.writeable() as pos:
            self.nbNode = len(pos)

        forceZ = [0.0 for _ in range(self.nbNode)]
        with self.DataDisplay.pointData.writeable() as pos:
            pos = forceZ

    def onAnimateEndEvent(self, event):
        position = self.MecaObject.position
        posRef = self.RefObject.position

        distX = []
        forceX = []
        nbNode = len(position)
        # print("1 =========> The len is : ", nbNode)
        for i in range(nbNode):
            distX += [math.fabs(position[i][2] - posRef[i][2])]
            # print("> =========", forceX)
            forceX += [math.fabs(0. * (position[i][2] - posRef[i][2]))]
        self.colorMapNode.DataDisplay.pointData = forceX


class buildElasticityRateAngulare(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.mechanicalObjectMO = args[0]
        self.insertionIndex = args[1]
        self.elsaticityRate = float(args[2])
        self.build_plasticity = args[3]

    def onAnimateEndEvent(self, event):
        """@todo: no need to have this boolean.
        switch this to the choice of use table or not"""
        if self.build_plasticity:
            self.buildPlasticity()
        else:
            self.usePlasticityWithTable()

    """============================================================
        find plasticity parameter for each soft beam of the implant
    ============================================================"""

    def buildPlasticity(self):
        positions = self.mechanicalObjectMO.position.value
        with self.mechanicalObjectMO.rest_position.writeable() as restPos:
            for idx, pos in enumerate(positions):
                # print("=+++++> idx: {} and rate is : {}".format(idx,pos[self.insertionIndex]))
                if abs(pos[self.insertionIndex]) > self.elsaticityRate:
                    # print("=+++++> idx: {} and rate is : {}".format(idx,pos[self.insertionIndex]))
                    restPos[idx][self.insertionIndex] = pos[self.insertionIndex] - \
                        self.elsaticityRate
                    # print("The positions: {} restPos: {}".format(pos[self.insertionIndex], restPos[idx][self.insertionIndex]))

    """ ============================================================
            This plasticity function is base on the tabe plasticityTab.
            here, each beam soft beam has it's own plasticity parameter.
        ============================================================"""

    def usePlasticityWithTable(self):
        positions = self.mechanicalObjectMO.position.value
        with self.mechanicalObjectMO.rest_position.writeable() as restPos:
            for idx, pos in enumerate(positions):
                if abs(pos[self.insertionIndex]) > plasticityTab[idx]:
                    restPos[idx][self.insertionIndex] = pos[self.insertionIndex] - \
                        plasticityTab[idx]


# #######################
# Function and class
# #######################
def BuildCosseratGeometry(config):
    # Define: the number of section, the total length and the length of each beam.

    [x, y, z] = config['init_pos']
    totalLength = config['tot_length']
    nbSection = config['nbSectionS']
    nbFrames = config['nbFramesF']

    lengthS = totalLength / nbSection
    # Define: the length of each beam in a list, the positions of each beam
    # (flexion, torsion), the abs of each section
    positionS = []
    longeurS = []
    temp = x
    curv_abs_inputS = [x]
    for i in range(nbSection):
        positionS.append([0, 0, 0])
        longeurS.append((((i + 1) * lengthS) - i * lengthS))
        temp += longeurS[i]
        curv_abs_inputS.append(temp)
    curv_abs_inputS[nbSection] = totalLength + x

    # ###############################################
    # Define: the number of frame and the length between each frame.
    # Define: the abs of each frame and the position of each frame.
    lengthF = totalLength / nbFrames
    framesF = []
    curv_abs_outputF = []
    cable_positionF = []
    for i in range(nbFrames):
        sol = i * lengthF
        framesF.append([sol + x, y, z, 0, 0, 0, 1])
        cable_positionF.append([sol + x, y, z])
        curv_abs_outputF.append(sol + x)

    framesF.append([totalLength + x, y, z, 0, 0, 0, 1])
    cable_positionF.append([totalLength + x, y, z])
    curv_abs_outputF.append(totalLength + x)

    return [positionS, curv_abs_inputS, longeurS, framesF, curv_abs_outputF, cable_positionF]


class MoveTargetProcess(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.targetMO = args[0]
        self.targetIndex = 0
        self.targetList = [[85.0, 0, 0.35857], [85.0, 0.5, 0.35857], [
            86.0, 1.5, 0.35857], [87.0, 2.5, 0.35857]]
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

            self.inputConstraintPointsMO.position.value = np.array(
                [[40, 0, 0], [45, 0, 0]])

            with self.inputConstraintPointsMO.position.writeable() as posA:
                print("1. List The position :", posA)
                # tip = self.needleCollisionMO.position[self.tipIndex]
                # print("Tip pose :", tip)
                # posA = [[40.0, 0., 0.], tip]
                # print("1. List The position :", posA)
                # posA.append(self.needleCollisionMO.position[self.tipIndex])

            indice = self.diffMapping.findData('indices')
            print("====> indices :",
                  self.diffMapping.addPointProcess([40.0, 0., 0.]))

        if key == "-":  # -
            with self.inputConstraintPointsMO.rest_position.writeable() as posA:
                sz = len(posA)
                print("Tip pose :", posA[sz - 1])


def Cube(
        name="Cube",
        surfaceMeshFileName="mesh/cube.obj",
        translation=None,
        rotation=None,
        uniformScale=1.,
        totalMass=300.,
        volume=20.,
        inertiaMatrix=None,
        color=None,
        isAStaticObject=False, parent=None):
    if color is None:
        color = [1., 1., 0.]
    if inertiaMatrix is None:
        inertiaMatrix = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    if rotation is None:
        rotation = [0., 0., 0.]
    if translation is None:
        translation = [0., 0., 0.]
    cubeNode = parent.addChild(name)
    cubeNode.addObject('MechanicalObject', name="mstate", template="Rigid3",
                       translation=translation, rotation=rotation)
    cubeNode.addObject('UniformMass', name="mass", vertexMass=[totalMass, volume, inertiaMatrix[:]],
                       showAxisSizeFactor='0.')
    # cubeNode.addObject('UncoupledConstraintCorrection')

    objectCollis = cubeNode.addChild('collision')
    objectCollis.addObject('MeshObjLoader', name="loader",
                           filename=surfaceMeshFileName, triangulate="true")
    objectCollis.addObject('MeshTopology', src="@loader")
    objectCollis.addObject('MechanicalObject')
    objectCollis.addObject('TriangleCollisionModel')
    objectCollis.addObject('LineCollisionModel')
    objectCollis.addObject('PointCollisionModel')
    objectCollis.addObject('RigidMapping')
    return

    # visu = cubeNode.addChild('visualisation')
    # path = surfaceMeshFileName
    # if path.endswith('.stl'):
    #     visu.addObject('MeshSTLLoader', name='loader', filename=path)
    # elif path.endswith('.obj'):
    #     visu.addObject('MeshObjLoader', name='loader', filename=path)
    # else:
    #     print("Extension not handled in STLIB/python/stlib/visuals for file: " + str(path))
    #
    # visu.addObject('MeshTopology', src='@loader', name='topo')
    # visu.addObject('OglModel', name="OglModel", src="@loader",
    #                rotation=rotation,
    #                translation=translation,
    #                scale3d=[1, 1, 1],
    #                color=color,
    #                updateNormals=False)
