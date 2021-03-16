# -*- coding: utf-8 -*-
"""Basic scene using Cosserat in SofaPython3.

Based on the work done with SofaPython. See POEMapping.py
"""

__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com, yinoussa.adagolodjo@inria.fr"
__version__ = "1.0.0"
__copyright__ = "(c) 2020,Inria"
__date__ = "March 8 2020"

import Sofa.Core
import os

from stlib3.scene import MainHeader
from splib3.numerics import Quat
# from stlib3.physics.collision import CollisionMesh

path = os.path.dirname(os.path.abspath(__file__))+'/../mesh/'


class Animation(Sofa.Core.Controller):

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.rigidBaseMO = args[0]
        self.rateAngularDeformMO = args[1]

        self.rate = 0.2
        self.angularRate = 0.02
        return

    def onKeypressedEvent(self, event):
        key = event['key']
        # ######## Rate angular #########
        if key == "I":  # +
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[5][1] += self.angularRate

        if key == "K":  # -
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[5][1] -= self.angularRate

        # ######## Reste rigid position #########
        if key == "+":  # up
            with self.rigidBaseMO.rest_position.writeable() as posA:
                qOld = Quat()
                for i in range(0, 4):
                    qOld[i] = posA[0][i+3]
                qNew = Quat.createFromEuler([0., self.angularRate, 0.], 'ryxz')
                qNew.normalize()
                qNew.rotateFromQuat(qOld)
                for i in range(0, 4):
                    posA[0][i+3] = qNew[i]

        if key == "-":  # down
            with self.rigidBaseMO.rest_position.writeable() as posA:
                qOld = Quat()
                for i in range(0, 4):
                    qOld[i] = posA[0][i+3]

                qNew = Quat.createFromEuler([0., -self.angularRate,  0.], 'ryxz')
                qNew.normalize()
                qNew.rotateFromQuat(qOld)
                for i in range(0, 4):
                    posA[0][i+3] = qNew[i]

        if ord(key) == 18:  # left
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][0] -= self.rate

        if ord(key) == 20:  # right
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][0] += self.rate

        if ord(key) == 21:  # down
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][1] -= self.rate

        if ord(key) == 19:  # up
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][1] += self.rate

def createFemCube(parentNode):
    FemNode = parentNode.addChild("FemNode")
    gelVolume = FemNode.addChild("gelVolume")
    gelVolume.addObject("RegularGridTopology", name="HexaTop", n="6 6 6", min="40 -6 -10", max="100 14 10")
    gelVolume.addObject("TetrahedronSetTopologyContainer", name="Container", position="@HexaTop.position")
    gelVolume.addObject("TetrahedronSetTopologyModifier", name="Modifier")
    gelVolume.addObject("Hexa2TetraTopologicalMapping", input="@HexaTop", output="@Container", swapping="false")

    GelSurface = FemNode.addChild("GelSurface")
    GelSurface.addObject("TriangleSetTopologyContainer", name="Container", position="@../GelVolume/HexaTop.position")
    # GelSurface.addObject("TriangleSetTopologyModifier", input="@../GelVolume/Container", output="@Container",
    #                      flipNormals="false")

    gelNode = FemNode.addChild("gelNode")
    # gelNode.addObject('VisualStyle', displayFlags='showVisualModels hideBehaviorModels showCollisionModels '
    #                                                  'hideMappings hideForceFields showWireframe '
    #                                                  'showInteractionForceFields hideForceFields')
    gelNode.addObject("EulerImplicitSolver", rayleighMass="0.1", rayleighStiffness="0.1")
    gelNode.addObject('SparseLDLSolver', name='preconditioner')
    gelNode.addObject('TetrahedronSetTopologyContainer', src="@../gelVolume/Container", name='container')
    # gelNode.addObject('TetrahedronSetTopologyModifier')
    gelNode.addObject('MechanicalObject', name='tetras', template='Vec3d')
    gelNode.addObject('TetrahedronFEMForceField', template='Vec3d', name='FEM', method='large',
                      poissonRatio='0.45', youngModulus='2000')
    # gelNode.addObject('UniformMass', totalMass='5')
    gelNode.addObject('BoxROI', name='ROI1', box='40 -6 -10 100 -4 10', drawBoxes='true')
    gelNode.addObject('RestShapeSpringsForceField', points='@ROI1.indices', stiffness='1e12')

    surfaceNode = gelNode.addChild("surfaceNode")
    surfaceNode.addObject('TriangleSetTopologyContainer', name="surfContainer", src="@../../GelSurface/Container")
    surfaceNode.addObject('MechanicalObject', name='msSurface')
    surfaceNode.addObject('TriangleCollisionModel', name='surface')
    surfaceNode.addObject('BarycentricMapping')

    gelNode.addObject('LinearSolverConstraintCorrection')

    return FemNode


def createScene(rootNode):

    MainHeader(rootNode, plugins=["SoftRobots", "SofaSparseSolver", "SofaPreconditioner",
                                  "SofaOpenglVisual", "CosseratPlugin", "BeamAdapter", "SofaDeformable",
                                  "SofaImplicitOdeSolver", 'SofaEngine', 'SofaMeshCollision', 'SofaSimpleFem',
                                  'SofaConstraint', 'SofaTopologyMapping'],
               repositoryPaths=[os.getcwd()])

    # rootNode.addObject('VisualStyle', displayFlags='showVisualModels showInteractionForceFields')
    rootNode.addObject('VisualStyle', displayFlags='showBehaviorModels hideCollisionModels hideBoundingCollisionModels '
                                                   'showForceFields hideInteractionForceFields showWireframe')
    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('GenericConstraintSolver', tolerance="1e-20", maxIterations="500", printLog="0")

    gravity = [0, 0, 0]
    rootNode.gravity.value = gravity
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.addObject('OglSceneFrame', style="Arrows", alignment="TopRight")
    # ###############
    # New adds to use the sliding Actuator
    ###############
    cableNode = rootNode.addChild('cableNode')
    cableNode.addObject('EulerImplicitSolver', firstOrder="0", rayleighStiffness="0.1", rayleighMass='0.1')
    cableNode.addObject('SparseLUSolver', name='solver')
    cableNode.addObject('GenericConstraintCorrection')

    # ###############
    # RigidBase
    ###############
    rigidBaseNode = cableNode.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode.addObject('MechanicalObject', template='Rigid3d', name="RigidBaseMO",
                                          position="0 0 0  0 0 0 1", showObject='1', showObjectScale='5.')
    rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="50000",
                            angularStiffness="50000", external_points="0", mstate="@RigidBaseMO", points="0",
                            template="Rigid3d")

    ###############
    # Rate of angular Deformation  (2 sections)
    ###############
    totalLen = 80.
    nbSectionS = 6
    lengthS = totalLen/nbSectionS
    positionS = []
    longeurS = []
    sumS = 0.
    curv_abs_inputS = [0.0]
    for i in range(nbSectionS):
        positionS.append([0, 0, 0])
        longeurS.append((((i + 1) * lengthS) - i * lengthS))
        sumS += longeurS[i]
        curv_abs_inputS.append(sumS)
    longeurS[nbSectionS - 1] = longeurS[nbSectionS - 1] + 1.
    curv_abs_inputS[nbSectionS] = 81.

    longeur = '15 15 15 15 6 15'  # beams size
    rateAngularDeformNode = cableNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject(
        'MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=positionS, showIndices="1")
    rateAngularDeformNode.addObject(
        'BeamHookeLawForceField', crossSectionShape='circular', length=longeur, radius='2.0', youngModulus='1.e12')
    ################################
    # Animation (to move the dofs) #
    ################################
    rootNode.addObject(Animation(RigidBaseMO, rateAngularDeformMO))

    ##############
    #   Frames   #
    ##############
    nbFramesF = 14
    lengthF = totalLen / nbFramesF
    framesF = []
    curv_abs_outputF = []
    cable_positionF = []
    for i in range(nbFramesF):
        sol = i * lengthF
        framesF.append([sol, 0, 0, 0, 0, 0, 1])
        cable_positionF.append([sol, 0, 0])
        curv_abs_outputF.append(sol)
    framesF.append([totalLen, 0, 0, 0, 0, 0, 1])
    print("=============> framesF : ", framesF)
    cable_positionF.append([totalLen, 0, 0])
    curv_abs_outputF.append(totalLen)

    # the node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.addChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)
    framesMO = mappedFrameNode.addObject(
        'MechanicalObject', template='Rigid3d', name="FramesMO", position=framesF, showObject='1', showObjectScale='1')

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=curv_abs_inputS,
                              curv_abs_output=curv_abs_outputF, input1=inputMO, input2=inputMO_rigid,
                              output=outputMO, debug='0', max=6.e-2, deformationAxis=1, nonColored="0", radius=5)

    slidingPoint = mappedFrameNode.addChild('slidingPoint')
    slidingPointMO = slidingPoint.addObject('MechanicalObject', name="cablePos", position=cable_positionF,
                                            showObject="1", showIndices="0")
    slidingPoint.addObject('IdentityMapping')

    # Create FEM Node
    femPos = [" 41.0 0 0 45 0 0 50 0 0 55 0 0 60 0 0 60 0 0  70 0 0 "]
    cubeNode = createFemCube(rootNode)
    gelNode = cubeNode.getChild('gelNode')
    femPoints = gelNode.addChild('femPoints')
    inputFEMCable = femPoints.addObject('MechanicalObject', name="pointsInFEM", position=femPos, showIndices="1")
    femPoints.addObject('BarycentricMapping')

    mappedPointsNode = slidingPoint.addChild('MappedPoints')
    femPoints.addChild(mappedPointsNode)
    mappedPoints = mappedPointsNode.addObject('MechanicalObject', template='Vec3d', position=femPos, name="FramesMO")

    inputCableMO = slidingPointMO.getLinkPath()
    inputFEMCableMO = inputFEMCable.getLinkPath()
    outputPointMO = mappedPoints.getLinkPath()

    mappedPointsNode.addObject('CosseratNeedleSlidingConstraint', name="QPConstraint")
    mappedPointsNode.addObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO, lastPointIsFixed=0,
                               input2=inputCableMO, output=outputPointMO, direction="@../../FramesMO.position")

    return rootNode

