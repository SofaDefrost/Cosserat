# -*- coding: utf-8 -*-

import os
import Sofa
from stlib.scene import MainHeader, ContactHeader
from splib.numerics import Vec3, Quat
from stlib.physics.collision import CollisionMesh

path = os.path.dirname(os.path.abspath(__file__)) + '/../mesh/'


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

    def onBeginAnimationStep(self, dt):
        pos = self.rigidBaseMO.findData('rest_position').value
        if pos[0][0] >= -17.0654:
            pos[0][0] -= self.rate
            self.rigidBaseMO.findData('rest_position').value = pos

    def onKeyPressed(self, c):

        ######### Rate angular #########
        if c == "+":  # +
            posA = self.rateAngularDeformMO.findData('rest_position').value
            posA[5][1] += self.angularRate
            self.rateAngularDeformMO.findData('rest_position').value = posA

        if c == "-":  # -
            posA = self.rateAngularDeformMO.findData('rest_position').value
            posA[5][1] -= self.angularRate
            self.rateAngularDeformMO.findData('rest_position').value = posA

        ######### Reste rigid position #########
        if ord(c) == 19:  # up
            posA = self.rigidBaseMO.findData('rest_position').value
            qOld = Quat()
            for i in range(0, 4):
                qOld[i] = posA[0][i + 3]

            qNew = Quat.createFromEuler([0., self.angularRate, 0.], 'ryxz')
            qNew.normalize()
            qNew.rotateFromQuat(qOld)
            for i in range(0, 4):
                posA[0][i + 3] = qNew[i]

            self.rigidBaseMO.findData('rest_position').value = posA

        if ord(c) == 21:  # down
            posA = self.rigidBaseMO.findData('rest_position').value
            qOld = Quat()
            for i in range(0, 4):
                qOld[i] = posA[0][i + 3]

            qNew = Quat.createFromEuler([0., -self.angularRate, 0.], 'ryxz')
            qNew.normalize()
            qNew.rotateFromQuat(qOld)
            for i in range(0, 4):
                posA[0][i + 3] = qNew[i]

            self.rigidBaseMO.findData('rest_position').value = posA

        if ord(c) == 18:  # left
            pos = self.rigidBaseMO.findData('rest_position').value
            pos[0][0] -= self.rate
            self.rigidBaseMO.findData('rest_position').value = pos

            # posA = self.rateAngularDeformMO.findData('position').value
            # for i in range(len(posA)):
            # posA[i][2] -= self.angularRate
            # self.rateAngularDeformMO.findData('position').value = posA

        if ord(c) == 20:  # right
            pos = self.rigidBaseMO.findData('rest_position').value
            pos[0][0] += self.rate
            self.rigidBaseMO.findData('rest_position').value = pos

            # for i in range(len(posA)):
            # posA[i][2] += self.angularRate
            # self.rateAngularDeformMO.findData('position').value = posA


def createScene(rootNode):
    MainHeader(rootNode, plugins=["SoftRobots", "SoftRobots.Inverse", "SofaPython", "SofaSparseSolver", "SofaConstraint",
                                  "SofaPreconditioner", "SofaOpenglVisual", "CosseratPlugin", "BeamAdapter"],
               repositoryPaths=[os.getcwd()])

    # rootNode.createObject('VisualStyle', displayFlags='showVisualModels hideBehaviorModels showCollisionModels '
    #                                                   'hideBoundingCollisionModels hideForceFields '
    #                                                   'showInteractionForceFields showWireframe')
    rootNode.createObject('VisualStyle', displayFlags='showVisualModels showInteractionForceFields showWireframe')

    rootNode.createObject('FreeMotionAnimationLoop')
    # rootNode.createObject('QPInverseProblemSolver', printLog='0')
    rootNode.createObject('GenericConstraintSolver', tolerance="1e-20", maxIterations="500", printLog="0")

    ContactHeader(rootNode, alarmDistance=2.5, contactDistance=2, frictionCoef=0.08)

    rootNode.gravity = "0 981 0"
    rootNode.createObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.createObject('OglSceneFrame', style="Arrows",
                          alignment="TopRight")

    showObject = "0"
    showIndices = "0"

    ##########################################
    # FEM Model                              #
    ##########################################
    finger = rootNode.createChild('finger')
    finger.createObject('EulerImplicitSolver', name='odesolver', firstOrder='0', rayleighMass=0.1,
                        rayleighStiffness=0.1)
    finger.createObject('SparseLDLSolver', name='preconditioner')

    # Add a componant to load a VTK tetrahedral mesh and expose the resulting topology in the scene .
    finger.createObject('MeshVTKLoader', name='loader', filename=path + 'finger.vtk', translation="-17.5 -12.5 7.5",
                        rotation="0 180 0")

    # finger.createObject('MeshExporter', name='loader', filename=path +'transFinger.vtk', exportAtEnd="true")

    finger.createObject('TetrahedronSetTopologyContainer', src='@loader', name='container')
    finger.createObject('TetrahedronSetTopologyModifier')
    # finger.createObject('TetrahedronSetTopologyAlgorithms', template='Vec3d')
    # finger.createObject('TetrahedronSetGeometryAlgorithms', template='Vec3d')

    # Create a mechanicaobject component to stores the DoFs of the model
    finger.createObject('MechanicalObject', name='tetras', template='Vec3d', showIndices='false',
                        showIndicesScale='4e-5', rx='0', dz='0')

    # Gives a mass to the model
    finger.createObject('UniformMass', totalMass='0.075')

    # Add a TetrahedronFEMForceField componant which implement an elastic material model
    # solved using the Finite Element Method on
    # tetrahedrons.
    finger.createObject('TetrahedronFEMForceField', template='Vec3d',
                        name='FEM', method='large', poissonRatio='0.45', youngModulus='1200')

    finger.createObject('BoxROI', name='ROI1',
                        box='-18 -15 -8 2 -3 8', drawBoxes='true')
    finger.createObject('RestShapeSpringsForceField',
                        points='@ROI1.indices', stiffness='1e12')

    ##########################################
    # Cable points                           #
    ##########################################
    # Mappe points inside the meca, this points will be use for the bilateral mapping
    FEMpos = [" 0.0 0 0 15 0 0 30 0 0 45 0 0 60 0 0 66 0 0 81 0.0 0.0"]
    # FEMpos = [" 81 0.0 0.0"]

    femPoints = finger.createChild('femPoints')
    inputFEMCable = femPoints.createObject('MechanicalObject', name="pointsInFEM", position=FEMpos, showObject="0",
                                           showIndices="0")
    femPoints.createObject('BarycentricMapping')

    # spheres = ["30. 0 0 48. 0 0 66 0 0 81. 0.0 0.0"]
    # spheresPos = finger.createChild('spheresPos')
    # spheresPosMec = spheresPos.createObject('MechanicalObject', name="spheresGoal", position=spheres, showObject="0",
    #                                        showIndices="1")
    # spheresPos.createObject('BarycentricMapping')
    finger.createObject('LinearSolverConstraintCorrection')

    ##########################################
    # Effector goal for interactive control  #
    ##########################################
    GoalPos = [[20.6307, 5.57305, -0.494896], [32.5759, 17.6405, -1.11956], [35.3802, 28.458, -1.52895],
               [36.3606, 42.5902, -1.85686]]
    goal = rootNode.createChild('goal')
    goal.createObject('EulerImplicitSolver', firstOrder='1')
    goal.createObject('CGLinearSolver', iterations='100', tolerance="1e-5", threshold="1e-5")
    # goal.createObject('MechanicalObject', name='goalMO', position='76.591 3.13521 0.35857')
    goal.createObject('MechanicalObject', name='goalMO', position=GoalPos)

    goal.createObject('SphereCollisionModel', radius='1.5')
    goal.createObject('UncoupledConstraintCorrection')
    goal.createObject('VisualStyle', displayFlags='showVisualModels showInteractionForceFields showCollisionModels')

    ##########################################
    # Visualization                          #
    ##########################################
    fingerVisu = finger.createChild('visu')
    fingerVisu.createObject(
        'MeshSTLLoader', filename=path + "finger.stl", name="loader", translation="-17.5 -12.5 7.5",
        rotation="0 180 0")
    # fingerVisu.createObject('STLExporter', filename=path+"transFinger", exportAtEnd="true")
    fingerVisu.createObject('OglModel', src="@loader", template='ExtVec3f', color="0.0 0.7 0.7")
    fingerVisu.createObject('BarycentricMapping')

    # ###############
    # New adds to use the sliding Actuator
    ###############
    cableNode = rootNode.createChild('cableNode')
    cableNode.createObject('EulerImplicitSolver', firstOrder="0", rayleighStiffness="1.0", rayleighMass='0.1')
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
    nbSectionS = 14
    lengthS = 80.0 / nbSectionS
    positionS = []
    longeurS = []
    sum = 0.
    # curv_abs_input = '0 15 30 45 60 66 81'
    curv_abs_inputS = []
    curv_abs_inputS.append(0.0)
    for i in range(nbSectionS):
        positionS.append([0, 0, 0])
        longeurS.append((((i+1)*lengthS) - i*lengthS))
        sum += longeurS[i]
        curv_abs_inputS.append(sum)
    longeurS[nbSectionS-1] = longeurS[nbSectionS-1] + 1.
    curv_abs_inputS[nbSectionS] = 81.

    print("=============> positionS : ", positionS)
    print("=============> longeurS : ", longeurS)
    print("=============> curv_abs_inputS : ", curv_abs_inputS)

    # position = ["0 0 0 " + "0 0 0 " + "0 0 0 " +
    #             "0 0 0 " + "0 0 0 " + "0 0 0 "]
    # longeur = '15 15 15 15 6 15'  # beams size
    rateAngularDeformNode = cableNode.createChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.createObject(
        'MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=positionS, showIndices="0")
    BeamHookeLawForce = rateAngularDeformNode.createObject(
        'BeamHookeLawForceField', crossSectionShape='circular', length=longeurS, radius='0.50', youngModulus='5e6')

    ################################
    # Animation (to move the dofs) #
    ################################
    anim = Animation(rigidBaseNode, rateAngularDeformNode)

    ##############
    #   Frames   #
    ##############
    nbFramesF = 28
    lengthF = 80.0/nbFramesF
    framesF = []
    curv_abs_outputF = []
    cable_positionF = []
    for i in range(nbFramesF):
        sol = i * lengthF
        framesF.append([sol, 0, 0,  0, 0, 0, 1])
        cable_positionF.append([sol, 0, 0])
        curv_abs_outputF.append(sol)
    framesF.append([81., 0, 0, 0, 0, 0, 1])
    cable_positionF.append([81., 0, 0])
    curv_abs_outputF.append(81.)

    # the node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.createChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)
    framesMO = mappedFrameNode.createObject('MechanicalObject', template='Rigid3d', name="FramesMO", position=framesF,
                                            showObject='0', showObjectScale='1')

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.createObject('DiscretCosseratMapping', curv_abs_input=curv_abs_inputS,
                                 curv_abs_output=curv_abs_outputF, input1=inputMO, input2=inputMO_rigid,
                                 output=outputMO, debug='0', max=6.e-2, deformationAxis=2, nonColored="0", radius=5)

    #  This create a new node in the scene. This node is appended to the finger's node.
    slidingPoint = mappedFrameNode.createChild('slidingPoint')

    # This create a MechanicalObject, a componant holding the degree of freedom of our
    # mechanical modelling. In the case of a cable it is a set of positions specifying
    # the points where the cable is passing by.
    slidingPointMO = slidingPoint.createObject('MechanicalObject', name="cablePos", position=cable_positionF,
                                               showObject="0", showIndices="0")
    slidingPoint.createObject('IdentityMapping')

    mappedPointsNode = slidingPoint.createChild('MappedPoints')
    femPoints.addChild(mappedPointsNode)
    mappedPoints = mappedPointsNode.createObject('MechanicalObject', template='Vec3d', position=FEMpos,
                                                 name="FramesMO", showObject='0', showObjectScale='0')

    inputCableMO = slidingPointMO.getLinkPath()
    inputFEMCableMO = inputFEMCable.getLinkPath()
    outputPointMO = mappedPoints.getLinkPath()

    mappedPointsNode.createObject('QPSlidingConstraint', name="QPConstraint")

    mappedPointsNode.createObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO,
                                  input2=inputCableMO, output=outputPointMO, direction="@../../FramesMO.position")
    ## Get the tree mstate links for the mapping

    return rootNode
