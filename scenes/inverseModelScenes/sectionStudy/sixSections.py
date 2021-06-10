# -*- coding: utf-8 -*-

import os
import Sofa
#from stlib3.scene import MainHeader, ContactHeader
#from splib.numerics import Vec3, Quat
#from stlib.physics.collision import CollisionMesh
from splib3.numerics import Vec3, Quat

path = os.path.dirname(os.path.abspath(__file__)) + '/../mesh/'


class Animation(Sofa.Core.Controller):
    
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
    #def __init__(self, rigidBaseNode, rateAngularDeformNode):
        self.rigidBaseNode = args[0]
        self.rateAngularDeformNode = args[1]

        self.rate = 0.2
        self.angularRate = 0.02
        
        self.initGraph()
        return

    def initGraph(self):
        self.rigidBaseMO = self.rigidBaseNode.getObject('RigidBaseMO')
        self.rateAngularDeformMO = self.rateAngularDeformNode.getObject(
            'rateAngularDeformMO')

    def onAnimateBeginEvent(self, dt):
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
            print("============> The rate angular :", posA)

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
        

    pluginNameList = 'SofaConstraint SofaDeformable SofaImplicitOdeSolver SofaMeshCollision SofaPreconditioner' \
                 ' SofaGeneralTopology SofaOpenglVisual SofaGeneralRigid SoftRobots SofaSparseSolver' \
                 ' CosseratPlugin BeamAdapter SofaBoundaryCondition'
    
    rootNode.addObject('RequiredPlugin', pluginName=pluginNameList, printLog='0')
    
    #MainHeader(rootNode, plugins=["SoftRobots", "SoftRobots.Inverse", "SofaPython", "SofaSparseSolver", "SofaConstraint",
                                  #"SofaPreconditioner", "SofaOpenglVisual", "CosseratPlugin", "BeamAdapter"],
               #repositoryPaths=[os.getcwd()])

    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showInteractionForceFields showWireframe')

    rootNode.findData('gravity').value = [9810., 0., 0.]
    rootNode.findData('dt').value = 0.01

    rootNode.addObject('FreeMotionAnimationLoop')
    #rootNode.addObject('DefaultPipeline', verbose="0")
    rootNode.addObject('CollisionPipeline', verbose="0")
    rootNode.addObject('BruteForceDetection', name="N2")
    #rootNode.addObject('DefaultContactManager', response="FrictionContact", responseParams="mu=0.08")
    rootNode.addObject('CollisionResponse', response="FrictionContact", responseParams="mu=0.08")
    rootNode.addObject('LocalMinDistance', name="Proximity", alarmDistance=0.3, contactDistance=0.26)
    # rootNode.addObject('QPInverseProblemSolver', printLog='0')
    rootNode.addObject('GenericConstraintSolver', tolerance="1e-20", maxIterations="500", printLog="0")
    rootNode.addObject('Camera', position="-35 0 280", lookAt="0 0 0")

    rootNode.gravity = "0 -0 0"
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.addObject('OglSceneFrame', style="Arrows",
                          alignment="TopRight")

    showObject = "0"
    showIndices = "0"
    return 

    """
        FEM Model
    """
    finger = rootNode.addChild('finger')
    finger.addObject('EulerImplicitSolver', name='odesolver', firstOrder='0', rayleighMass=0.1,
                        rayleighStiffness=0.1)
    finger.addObject('SparseLDLSolver', name='preconditioner')

    # Add a componant to load a VTK tetrahedral mesh and expose the resulting topology in the scene .
    finger.addObject('MeshVTKLoader', name='loader', filename=path + 'finger.vtk', translation="-17.5 -12.5 7.5",
                        rotation="0 180 0")

    # finger.addObject('MeshExporter', name='loader', filename=path +'transFinger.vtk', exportAtEnd="true")

    finger.addObject('TetrahedronSetTopologyContainer', src='@loader', name='container')
    finger.addObject('TetrahedronSetTopologyModifier')
    # finger.addObject('TetrahedronSetTopologyAlgorithms', template='Vec3d')
    # finger.addObject('TetrahedronSetGeometryAlgorithms', template='Vec3d')

    # Create a mechanicaobject component to stores the DoFs of the model
    finger.addObject('MechanicalObject', name='tetras', template='Vec3d', showIndices='false',
                        showIndicesScale='4e-5', rx='0', dz='0')

    # Gives a mass to the model
    finger.addObject('UniformMass', totalMass='0.005')

    # Add a TetrahedronFEMForceField componant which implement an elastic material model
    # solved using the Finite Element Method on
    # tetrahedrons.
    finger.addObject('TetrahedronFEMForceField', template='Vec3d',
                        name='FEM', method='large', poissonRatio='0.45', youngModulus='500')

    finger.addObject('BoxROI', name='ROI1',
                        box='-18 -15 -8 2 -3 8', drawBoxes='true')
    finger.addObject('RestShapeSpringsForceField',
                        points='@ROI1.indices', stiffness='1e12')

    ##########################################
    # Cable points                           #
    ##########################################
    # Mappe points inside the meca, this points will be use for the bilateral mapping
    FEMpos = [" 0.0 0 0 15 0 0 30 0 0 45 0 0 60 0 0 66 0 0 81 0.0 0.0"]
    # FEMpos = [" 81 0.0 0.0"]

    femPoints = finger.addChild('femPoints')
    inputFEMCable = femPoints.addObject('MechanicalObject', name="pointsInFEM", position=FEMpos, showObject="0",
                                           showIndices="0")
    femPoints.addObject('BarycentricMapping')

    # spheres = ["30. 0 0 48. 0 0 66 0 0 81. 0.0 0.0"]
    # spheresPos = finger.addChild('spheresPos')
    # spheresPosMec = spheresPos.addObject('MechanicalObject', name="spheresGoal", position=spheres, showObject="0",
    #                                        showIndices="1")
    # spheresPos.addObject('BarycentricMapping')
    finger.addObject('LinearSolverConstraintCorrection')

    ##########################################
    # Effector goal for interactive control  #
    ##########################################
    GoalPos = [[20.6307, 5.57305, -0.494896], [32.5759, 17.6405, -1.11956], [35.3802, 28.458, -1.52895],
               [36.3606, 42.5902, -1.85686]]
    goal = rootNode.addChild('goal')
    goal.addObject('EulerImplicitSolver', firstOrder='1')
    goal.addObject('CGLinearSolver', iterations='100', tolerance="1e-5", threshold="1e-5")
    # goal.addObject('MechanicalObject', name='goalMO', position='76.591 3.13521 0.35857')
    goal.addObject('MechanicalObject', name='goalMO', position=GoalPos)

    goal.addObject('SphereCollisionModel', radius='1.5')
    goal.addObject('UncoupledConstraintCorrection')
    goal.addObject('VisualStyle', displayFlags='showVisualModels showInteractionForceFields showCollisionModels')

    ##########################################
    # Visualization                          #
    ##########################################
    fingerVisu = finger.addChild('visu')
    fingerVisu.addObject(
        'MeshSTLLoader', filename=path + "finger.stl", name="loader", translation="-17.5 -12.5 7.5",
        rotation="0 180 0")
    # fingerVisu.addObject('STLExporter', filename=path+"transFinger", exportAtEnd="true")
    fingerVisu.addObject('OglModel', src="@loader", template='ExtVec3f', color="0.0 0.7 0.7")
    fingerVisu.addObject('BarycentricMapping')

    # ###############
    # New adds to use the sliding Actuator
    ###############
    cableNode = rootNode.addChild('cableNode')
    cableNode.addObject('EulerImplicitSolver', firstOrder="0", rayleighStiffness="1.0", rayleighMass='0.1')
    cableNode.addObject('SparseLUSolver', name='solver')
    cableNode.addObject('GenericConstraintCorrection')

    # ###############
    # RigidBase
    ###############
    rigidBaseNode = cableNode.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode.addObject('MechanicalObject', template='Rigid3d',
                                             name="RigidBaseMO", position="0 0 0  0 0 0 1", showObject='1',
                                             showObjectScale='5.')
    rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="50000",
                               angularStiffness="50000", external_points="0", mstate="@RigidBaseMO", points="0",
                               template="Rigid3d")

    ###############
    # Rate of angular Deformation  (2 sections)
    ###############
    nbSectionS = 3
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
    rateAngularDeformNode = cableNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject(
        'MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=positionS, showIndices="0")
    BeamHookeLawForce = rateAngularDeformNode.addObject(
        'BeamHookeLawForceField', crossSectionShape='circular', length=longeurS, radius='0.50', youngModulus='5e6')

    ################################
    # Animation (to move the dofs) #
    ################################
    anim = Animation(rigidBaseNode, rateAngularDeformNode)

    ##############
    #   Frames   #
    ##############
    nbFramesF = 6
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
    print("=============> framesF : ", framesF)
    cable_positionF.append([81., 0, 0])
    curv_abs_outputF.append(81.)

    # the node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.addChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)
    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d', name="FramesMO", position=framesF,
                                            showObject='0', showObjectScale='1')

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscretCosseratMapping', curv_abs_input=curv_abs_inputS,
                                 curv_abs_output=curv_abs_outputF, input1=inputMO, input2=inputMO_rigid,
                                 output=outputMO, debug='0', max=6.e-2, deformationAxis=2, nonColored="0", radius=5)

    #  This create a new node in the scene. This node is appended to the finger's node.
    slidingPoint = mappedFrameNode.addChild('slidingPoint')

    # This create a MechanicalObject, a componant holding the degree of freedom of our
    # mechanical modelling. In the case of a cable it is a set of positions specifying
    # the points where the cable is passing by.
    slidingPointMO = slidingPoint.addObject('MechanicalObject', name="cablePos", position=cable_positionF,
                                               showObject="0", showIndices="0")
    slidingPoint.addObject('IdentityMapping')

    mappedPointsNode = slidingPoint.addChild('MappedPoints')
    femPoints.addChild(mappedPointsNode)
    mappedPoints = mappedPointsNode.addObject('MechanicalObject', template='Vec3d', position=FEMpos,
                                                 name="FramesMO", showObject='0', showObjectScale='1')

    inputCableMO = slidingPointMO.getLinkPath()
    inputFEMCableMO = inputFEMCable.getLinkPath()
    outputPointMO = mappedPoints.getLinkPath()

    mappedPointsNode.addObject('QPSlidingConstraint', name="QPConstraint")
    mappedPointsNode.addObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO,
                                  input2=inputCableMO, output=outputPointMO, direction="@../../FramesMO.position")
    ## Get the tree mstate links for the mapping

    return rootNode
