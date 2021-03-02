# -*- coding: utf-8 -*-

import os
import Sofa
from stlib.scene import MainHeader, ContactHeader
from splib.numerics import Vec3, Quat
from stlib.physics.collision import CollisionMesh

path = os.path.dirname(os.path.abspath(__file__))+'/../mesh/'


class Animation(Sofa.PythonScriptController):

    def __init__(self, rigidBaseNode, rateAngularDeformNode):
        self.rigidBaseNode = rigidBaseNode
        self.rateAngularDeformNode = rateAngularDeformNode

        self.rate = 0.2
        self.angularRate = 0.02
        return

    def initGraph(self, nodeRigid):
        self.rigidBaseMO = self.rigidBaseNode.getObject('RigidBaseMO')
        self.rateAngularDeformMO = self.rateAngularDeformNode.getObject(
            'rateAngularDeformMO')

    def onKeyPressed(self, c):

        ######### Rate angular #########
        if c == "I":  # +
            posA = self.rateAngularDeformMO.findData('rest_position').value
            posA[5][1] += self.angularRate
            self.rateAngularDeformMO.findData('rest_position').value = posA

        if c == "K":  # -
            posA = self.rateAngularDeformMO.findData('rest_position').value
            posA[5][1] -= self.angularRate
            self.rateAngularDeformMO.findData('rest_position').value = posA

        ######### Reste rigid position #########
        if c == "+":  # up
            posA = self.rigidBaseMO.findData('rest_position').value
            qOld = Quat()
            for i in range(0, 4):
                qOld[i] = posA[0][i+3]

            qNew = Quat.createFromEuler([ 0., self.angularRate, 0.], 'ryxz')
            qNew.normalize()
            qNew.rotateFromQuat(qOld)
            for i in range(0, 4):
                posA[0][i+3] = qNew[i]

            self.rigidBaseMO.findData('rest_position').value = posA

        if c == "-":  # down
            posA = self.rigidBaseMO.findData('rest_position').value
            qOld = Quat()
            for i in range(0, 4):
                qOld[i] = posA[0][i+3]

            qNew = Quat.createFromEuler([0., -self.angularRate,  0.], 'ryxz')
            qNew.normalize()
            qNew.rotateFromQuat(qOld)
            for i in range(0, 4):
                posA[0][i+3] = qNew[i]

            self.rigidBaseMO.findData('rest_position').value = posA

        if ord(c) == 18:  # left
            pos = self.rigidBaseMO.findData('rest_position').value
            pos[0][0] -= self.rate
            self.rigidBaseMO.findData('rest_position').value = pos

            #posA = self.rateAngularDeformMO.findData('position').value
            #for i in range(len(posA)):
                #posA[i][2] -= self.angularRate
            #self.rateAngularDeformMO.findData('position').value = posA

        if ord(c) == 20:  # right
            pos = self.rigidBaseMO.findData('rest_position').value
            pos[0][0] += self.rate
            self.rigidBaseMO.findData('rest_position').value = pos

        if ord(c) == 21:  # down
            pos = self.rigidBaseMO.findData('rest_position').value
            pos[0][1] -= self.rate
            self.rigidBaseMO.findData('rest_position').value = pos

        if ord(c) == 19:  # up
            pos = self.rigidBaseMO.findData('rest_position').value
            pos[0][1] += self.rate
            self.rigidBaseMO.findData('rest_position').value = pos

            #for i in range(len(posA)):
                #posA[i][2] += self.angularRate
            #self.rateAngularDeformMO.findData('position').value = posA


def createFemCube(parentNode):
    FemNode = parentNode.createChild("FemNode")
    gelVolume = FemNode.createChild("gelVolume")
    gelVolume.createObject("RegularGridTopology", name="HexaTop", n="6 6 6", min="40 -6 -10", max="100 14 10")
    gelVolume.createObject("TetrahedronSetTopologyContainer", name="Container", position="@HexaTop.position")
    gelVolume.createObject("TetrahedronSetTopologyModifier", name="Modifier")
    gelVolume.createObject("Hexa2TetraTopologicalMapping", input="@HexaTop", output="@Container", swapping="false")

    GelSurface = FemNode.createChild("GelSurface")
    GelSurface.createObject("TriangleSetTopologyContainer", name="Container", position="@../GelVolume/HexaTop.position")
    GelSurface.createObject("TriangleSetTopologyModifier", nput="@../GelVolume/Container", output="@Container",
                            flipNormals="false")

    gelNode = FemNode.createChild("gelNode")
    # gelNode.createObject('VisualStyle', displayFlags='showVisualModels hideBehaviorModels showCollisionModels '
    #                                                  'hideMappings hideForceFields showWireframe '
    #                                                  'showInteractionForceFields hideForceFields')
    gelNode.createObject("EulerImplicitSolver", rayleighMass="0.1", rayleighStiffness="0.1")
    gelNode.createObject('SparseLDLSolver', name='preconditioner')
    gelNode.createObject('TetrahedronSetTopologyContainer', src="@../gelVolume/Container", name='container')
    # gelNode.createObject('TetrahedronSetTopologyModifier')
    gelNode.createObject('MechanicalObject', name='tetras', template='Vec3d')
    gelNode.createObject('TetrahedronFEMForceField', template='Vec3d', name='FEM', method='large', poissonRatio='0.45',
                         youngModulus='2000')
    # gelNode.createObject('UniformMass', totalMass='5')
    gelNode.createObject('BoxROI', name='ROI1', box='40 -6 -10 100 -4 10', drawBoxes='true')
    gelNode.createObject('RestShapeSpringsForceField', points='@ROI1.indices', stiffness='1e12')

    surfaceNode = gelNode.createChild("surfaceNode")
    surfaceNode.createObject('TriangleSetTopologyContainer', name="surfContainer", src="@../../GelSurface/Container")
    surfaceNode.createObject('MechanicalObject', name='msSurface')
    surfaceNode.createObject('TriangleCollisionModel', name='surface')
    surfaceNode.createObject('BarycentricMapping')

    gelNode.createObject('LinearSolverConstraintCorrection')

    return FemNode


def createScene(rootNode):

    MainHeader(rootNode, plugins=["SoftRobots", "SofaPython", "SofaSparseSolver", "SofaPreconditioner",
                                  "SofaOpenglVisual", "CosseratPlugin", "BeamAdapter", "SofaDeformable",
                                  "SofaImplicitOdeSolver", 'SofaEngine', 'SofaMeshCollision', 'SofaSimpleFem',
                                  'SofaConstraint', 'SofaTopologyMapping'],
               repositoryPaths=[os.getcwd()])

    # rootNode.createObject('VisualStyle', displayFlags='showVisualModels showInteractionForceFields')
    rootNode.createObject('VisualStyle', displayFlags='showBehaviorModels hideCollisionModels '
                                                      'hideBoundingCollisionModels showForceFields '
                                                      'hideInteractionForceFields showWireframe')
    rootNode.createObject('FreeMotionAnimationLoop')
    rootNode.createObject('GenericConstraintSolver', tolerance="1e-20", maxIterations="500", printLog="0")

    rootNode.gravity = "0 0 0"
    rootNode.createObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.createObject('OglSceneFrame', style="Arrows", alignment="TopRight")
    # ###############
    # New adds to use the sliding Actuator
    ###############
    cableNode = rootNode.createChild('cableNode')
    cableNode.createObject('EulerImplicitSolver', firstOrder="0", rayleighStiffness="0.1", rayleighMass='0.1')
    cableNode.createObject('SparseLUSolver', name='solver')
    cableNode.createObject('GenericConstraintCorrection')

    # ###############
    # RigidBase
    ###############
    rigidBaseNode = cableNode.createChild('rigidBase')
    RigidBaseMO = rigidBaseNode.createObject('MechanicalObject', template='Rigid3d',
                                             name="RigidBaseMO", position="0 0 0  0 0 0 1", showObject='1',
                                             showObjectScale='5.')
    rigidBaseNode.createObject('RestShapeSpringsForceField', name='spring', stiffness="50000",
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
    sumS = 0.; curv_abs_inputS = []
    curv_abs_inputS.append(0.0)
    for i in range(nbSectionS):
        positionS.append([0, 0, 0])
        longeurS.append((((i + 1) * lengthS) - i * lengthS))
        sumS += longeurS[i]
        curv_abs_inputS.append(sumS)
    longeurS[nbSectionS - 1] = longeurS[nbSectionS - 1] + 1.
    curv_abs_inputS[nbSectionS] = 81.

    longeur = '15 15 15 15 6 15'  # beams size
    rateAngularDeformNode = cableNode.createChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.createObject(
        'MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=positionS, showIndices="1")
    BeamHookeLawForce = rateAngularDeformNode.createObject(
        'BeamHookeLawForceField', crossSectionShape='circular', length=longeur, radius='2.0', youngModulus='1.e12')
    ################################
    # Animation (to move the dofs) #
    ################################
    anim = Animation(rigidBaseNode, rateAngularDeformNode)

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
    mappedFrameNode = rigidBaseNode.createChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)
    framesMO = mappedFrameNode.createObject(
        'MechanicalObject', template='Rigid3d', name="FramesMO", position=framesF, showObject='1', showObjectScale='1')

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.createObject('DiscretCosseratMapping', curv_abs_input=curv_abs_inputS,
                                 curv_abs_output=curv_abs_outputF, input1=inputMO, input2=inputMO_rigid,
                                 output=outputMO, debug='0', max=6.e-2, deformationAxis=1, nonColored="0", radius=5)

    slidingPoint = mappedFrameNode.createChild('slidingPoint')
    slidingPointMO = slidingPoint.createObject('MechanicalObject', name="cablePos", position=cable_positionF,
                                               showObject="1", showIndices="0")
    slidingPoint.createObject('IdentityMapping')

    # Create FEM Node
    femPos = [" 41.0 0 0 45 0 0 50 0 0 55 0 0 60 0 0 60 0 0  70 0 0 "]
    cubeNode = createFemCube(rootNode)
    gelNode = cubeNode.getChild('gelNode')
    femPoints = gelNode.createChild('femPoints')
    inputFEMCable = femPoints.createObject('MechanicalObject', name="pointsInFEM",
                                           position=femPos, showIndices="1")
    femPoints.createObject('BarycentricMapping')

    mappedPointsNode = slidingPoint.createChild('MappedPoints')
    femPoints.addChild(mappedPointsNode)
    mappedPoints = mappedPointsNode.createObject('MechanicalObject', template='Vec3d', position=femPos,
                                                 name="FramesMO")

    inputCableMO = slidingPointMO.getLinkPath()
    inputFEMCableMO = inputFEMCable.getLinkPath()
    outputPointMO = mappedPoints.getLinkPath()

    mappedPointsNode.createObject('CosseratNeedleSlidingConstraint', name="QPConstraint")
    mappedPointsNode.createObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO,
                                  lastPointIsFixed=0,
                                  input2=inputCableMO, output=outputPointMO, direction="@../../FramesMO.position")

    return rootNode

    ##########################################
    # FEM Model                              #
    ##########################################
    # fingerTrans = [-12.5, -12.5, 7.5]
    # fingerRot = [0, 180, 0]
    # finger = rootNode.createChild('finger')
    # finger.createObject('EulerImplicitSolver', name='odesolver', firstOrder='0', rayleighMass=0.1,
    #                     rayleighStiffness=0.1)
    # finger.createObject('SparseLDLSolver', name='preconditioner')
    #
    # # Add a componant to load a VTK tetrahedral mesh and expose the resulting topology in the scene .
    # finger.createObject('MeshVTKLoader', name='loader', filename=path + 'finger.vtk', translation=fingerTrans,
    #                     rotation=fingerRot)
    # finger.createObject('TetrahedronSetTopologyContainer', src='@loader', name='container')
    # finger.createObject('TetrahedronSetTopologyModifier')
    #
    # # Create a mechanicaobject component to stores the DoFs of the model
    # finger.createObject('MechanicalObject', name='tetras', template='Vec3d', showIndices='false',
    #                     position="@loader.position", showIndicesScale='4e-5', rx='0', dz='0')
    #
    # finger.createObject('TetrahedronFEMForceField', template='Vec3d',
    #                     name='FEM', method='large', poissonRatio='0.45', youngModulus='500')
    # # Gives a mass to the model
    # finger.createObject('UniformMass', totalMass='0.075')
    #
    # # Add a TetrahedronFEMForceField componant which implement an elastic material model
    # # solved using the Finite Element Method on
    # # tetrahedrons.
    # finger.createObject('BoxROI', name='ROI1', box='-18 -15 -8 2 -3 8', drawBoxes='true')
    # finger.createObject('RestShapeSpringsForceField', points='@ROI1.indices', stiffness='1e12')
    #
    # ##########################################
    # # Cable points                           #
    # ##########################################
    # # Mappe points inside the meca, this points will be use for the bilateral mapping
    # FEMpos = [" 0.0 0 0 15 0 0 30 0 0 45 0 0 60 0 0 66 0 0 "]
    #
    # femPoints = finger.createChild('femPoints')
    # inputFEMCable = femPoints.createObject('MechanicalObject', name="pointsInFEM",
    #                                        position=FEMpos, showIndices="1")
    # femPoints.createObject('BarycentricMapping')
    #
    # ##########################################
    # # Visualization                          #
    # ##########################################
    # fingerVisu = finger.createChild('visu')
    # fingerVisu.createObject(
    #     'MeshSTLLoader', filename=path + "finger.stl", name="loader", translation=fingerTrans,
    #     rotation=fingerRot)
    # fingerVisu.createObject('OglModel', src="@loader",
    #                         template='ExtVec3f', color="0.0 0.7 0.7")
    # fingerVisu.createObject('BarycentricMapping')
    #
    # finger.createObject('LinearSolverConstraintCorrection')

    # mappedPointsNode = slidingPoint.createChild('MappedPoints')

    # femPoints.addChild(mappedPointsNode)
    # mappedPoints = mappedPointsNode.createObject('MechanicalObject', template='Vec3d', position=FEMpos,
    #                                              name="FramesMO", showObject='1', showObjectScale='1')
    #
    #
    #
    # inputCableMO = slidingPointMO.getLinkPath()
    # inputFEMCableMO = inputFEMCable.getLinkPath()
    # outputPointMO = mappedPoints.getLinkPath()
    #
    # mappedPointsNode.createObject('QPSlidingConstraint', name="QPConstraint")
    # mappedPointsNode.createObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO,
    #                               input2=inputCableMO, output=outputPointMO, direction="@../../FramesMO.position")


