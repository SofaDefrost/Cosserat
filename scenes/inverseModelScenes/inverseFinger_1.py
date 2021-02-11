# -*- coding: utf-8 -*-

import os
import Sofa
from stlib.scene import MainHeader, ContactHeader
from math import sqrt
from math import pi, sin, cos
from stlib.physics.collision import CollisionMesh

path = os.path.dirname(os.path.abspath(__file__))+'/mesh/'


class Animation(Sofa.PythonScriptController):

    def __init__(self, goal):
        self.goal = goal
        self.rateX = 0.2 
        self.rate = 0.02
        self.time = 0.01
        self.iter = 0
        self.goalPos = [85.591, 1.13521, 0.35857]
        
        #self.positions = [ [85.591  , 1.13521 , 0.35857],[83.8771269207  , 12.6336038452 , 0.35857], [79.4561670934  , 23.363504223 , 0.35857],
                         #[72.4872145094  , 30.7610327247 , 0.35857],[53.7702080239  , 43.9900691183 , 0.35857],[42.9324103894  , 48.8920846724 , 0.35857],
                         #[36.9324103894  , 49.8920846724 , 0.35857],[33.9324103894  , 51.8920846724 , 0.35857],
                         #[29.5022454485  , 53.5927137725 , -1.35857],[20.5022454485  , 51.5927137725 , -2.35857],[7.90460823185  , 47.6481704134 , -5.40]]
        
        self.positions = [ [85.591  , 1.13521 , 0.35857],[83.8771269207  , 12.6336038452 , 0.35857], [79.4561670934  , 23.363504223 , 0.35857],
                         [72.4872145094  , 30.7610327247 , 0.35857]]
        
        self.radius = 85.
        self.theta = 0.01
        self.dt = 0.01
        self.step = 1
        self.maxTheta = 1.2
        self.minTheta = 0.01
        self.coeff = 0.9
        self.seuil = 0.3
        
        return
    
    def initGraph(self, nodeRigid):
        self.goalMO = self.goal.getObject('goalMO')     
        
    
    def moveForward(self, axes=[0],dt=0.01):
        pos = [0.,0.,0.]
        if (self.theta < self.maxTheta):
            pos[0] = self.radius * cos(self.theta)
            for axe in axes:
                pos[axe] = self.coeff * self.radius * sin(self.theta)* (0.9 - 0.21*self.theta)
            self.theta += dt 
        else :
            self.step +=1 ;
            pos = self.goalMO.findData('position').value
        return pos
    
    def moveBack(self, axes=[0],dt=0.01):
        pos = [0.,0.,0.]
        if (self.theta > self.minTheta):
            pos[0] = self.radius * cos(self.theta)
            for axe in axes:
                pos[axe] = self.coeff * self.radius * sin(self.theta)* (0.9 - 0.21*self.theta)                
            self.theta -= dt 
            
        else :
            self.step +=1 ;
            pos = self.goalMO.findData('position').value

        return pos
    
    def onBeginAnimationStep(self, dt):
        
        pos = self.goalMO.findData('position').value
                
        if self.step == 1:
            axes = [1]
            pos = self.moveForward(axes,0.005)
            self.goalMO.findData('position').value = pos
        elif self.step == 2:
            axes = [1]
            pos = self.moveBack(axes,0.005)
            self.goalMO.findData('position').value = pos
    
#    def onBeginAnimationStep(self, dt):
#        
#        position = self.positions[self.iter]
#        pos = self.goalMO.findData('position').value
#        y1 = position[1]
#        y2 = pos[0][1] + self.rate
#        if y1 > y2 :
#            pos[0][0] -= self.rate
#            if self.rate < 3:
#                pos[0][1] += 2*self.rate
#            else:
#                pos[0][1] += 2*self.rate
#                
#            self.goalMO.findData('position').value = pos
#        else:
#            position = self.positions[self.iter]
#            if self.iter < len(self.positions)-1 :
#                self.iter +=1
#            self.goalMO.findData('position').value = position
#        
#        position = self.positions[self.iter]
    
    def onKeyPressed(self, c):

        if ord(c) == 19:  # up
            pos = self.goalMO.findData('position').value
            pos[0][0] -= self.rateX
            pos[0][1] += self.rateY
            self.goalMO.findData('position').value = pos
            #x = pos[0][0]            
            #y = 3*( sqrt((x*x)/49.) - 1)                        
            #pos[0][1] = y
            
            print("=======> Position :", pos)
            
            

        if c == 'I':  # down
            position = self.positions[self.iter]
            if self.iter < len(self.positions)-1 :
                self.iter +=1
            self.goalMO.findData('position').value = position
            
            #print(" [" +str(pos[0][0])+ "  , "+str(pos[0][1])+ " , "+str(pos[0][2])+ "],")


        #if ord(c) == 18:  # left
            #pos = self.rigidBaseMO.findData('position').value
            #pos[0][2] -= self.rate
            #self.rigidBaseMO.findData('position').value = pos
            #print("=======> Position :", pos)

            #posA = self.rateAngularDeformMO.findData('position').value
            #for i in range(len(posA)):
                #posA[i][2] -= self.angularRate
            #self.rateAngularDeformMO.findData('position').value = posA

        #if ord(c) == 20:  # right
            #pos = self.rigidBaseMO.findData('position').value
            #pos[0][2] += self.rate
            #self.rigidBaseMO.findData('position').value = pos
            #print("=======> Position :", pos)

            #posA = self.rateAngularDeformMO.findData('position').value
            #for i in range(len(posA)):
                #posA[i][2] += self.angularRate
            #self.rateAngularDeformMO.findData('position').value = posA


def createScene(rootNode):

    MainHeader(rootNode, plugins=["SoftRobots", "SoftRobots.Inverse", "SofaPython",
                                  "SofaSparseSolver", "SofaPreconditioner", "SofaOpenglVisual", "CosseratPlugin",
                                  "BeamAdapter"], repositoryPaths=[os.getcwd()])

    #rootNode.createObject('VisualStyle', displayFlags='showVisualModels hideBehaviorModels showCollisionModels
    # hideBoundingCollisionModels hideForceFields showInteractionForceFields showWireframe')
    rootNode.createObject('VisualStyle', displayFlags='showVisualModels showInteractionForceFields')

    rootNode.createObject('FreeMotionAnimationLoop')
    rootNode.createObject('QPInverseProblemSolver', printLog='0')
    # rootNode.createObject('GenericConstraintSolver', tolerance="1e-20", maxIterations="5 00", printLog="0")

    rootNode.gravity = "0 0 0"
    rootNode.createObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.createObject('OglSceneFrame', style="Arrows",
                          alignment="TopRight")    
    #ContactHeader(rootNode, alarmDistance=4,contactDistance=3, frictionCoef=0.08)

    ##########################################
    # Effector goal for interactive control  #
    ##########################################
    goal = rootNode.createChild('goal')
    goal.createObject('EulerImplicitSolver', firstOrder='1')
    goal.createObject('CGLinearSolver', iterations='100', tolerance="1e-5", threshold="1e-5")
    # goal.createObject('MechanicalObject', name='goalMO', position='76.591 3.13521 0.35857')
    goal.createObject('MechanicalObject', name='goalMO', position='85.591 1.13521 0.35857')
    
    goal.createObject('SphereCollisionModel', radius='5')
    goal.createObject('UncoupledConstraintCorrection')
    goal.createObject('VisualStyle', displayFlags='showVisualModels showInteractionForceFields showCollisionModels')
    #9.81   48.6217  -4.43186
    
    ################################
    # Animation (to move the dofs) #
    ################################
    anim = Animation(goal)

    

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
    # rigidBaseNode.createObject('EulerImplicitSolver', firstOrder="0", rayleighStiffness="1.0", rayleighMass='0.1')
    # rigidBaseNode.createObject('SparseLUSolver', name='solver')
    # rigidBaseNode.createObject('GenericConstraintCorrection')
    RigidBaseMO = rigidBaseNode.createObject('MechanicalObject', template='Rigid3d', name="RigidBaseMO",
                                             position="0 0 0  0 0 0 1", showObject='0', showObjectScale='0.1')
    rigidBaseNode.createObject('PartialFixedConstraint', fixedDirections="0 1 1 1 1 1", indices="0")
#    rigidBaseNode.createObject('RestShapeSpringsForceField', name='spring', stiffness="100",angularStiffness="100", external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")
    rigidBaseNode.createObject('UniformMass', totalMass="0.001", template="Rigid3d" )
    rigidBaseNode.createObject('SlidingActuator', name="SlidingActuator0", template='Rigid3d', direction='1 0 0 0 0 0', indices=0,  maxForce='10000', minForce='-2000') 
    # rigidBaseNode.createObject('SlidingActuator', name="SlidingActuator1", template='Rigid3d',direction='0 0 0 1 0 0', indices=0,  maxForce='10000', minForce='-3000')
    # rigidBaseNode.createObject('SlidingActuator', name="SlidingActuator2", template='Rigid3d',direction='0 0 1 0 0 0', indices=0,  maxForce='10000', minForce='-3000')                                                                                                               
    # rigidBaseNode.createObject('SlidingActuator', name="SlidingActuator3", template='Rigid3d',direction='0 1 0 0 0 0', indices=0,  maxForce='10000', minForce='-3000')                                                                                                               

    ###############
    # Rate of angular Deformation  (2 sections)
    ###############
    position = ["0 0 0 " + "0 0 0 " + "0 0 0 " +
                "0 0 0 " + "0 0 0 " + "0 0 0 "]
    longeur = '15 15 15 15 6 15'  # beams size
    rateAngularDeformNode = cableNode.createChild('rateAngularDeform')
    # rateAngularDeformNode.createObject('EulerImplicitSolver', firstOrder="0", rayleighStiffness="1.0", rayleighMass='0.1')
    # rateAngularDeformNode.createObject('SparseLUSolver', name='solver')
    # rateAngularDeformNode.createObject('GenericConstraintCorrection')
    rateAngularDeformMO = rateAngularDeformNode.createObject('MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=position)
    
    BeamHookeLawForce = rateAngularDeformNode.createObject('BeamHookeLawForceField', crossSectionShape='circular', length=longeur, radius='0.5', youngModulus='5e6')
    
    # rateAngularDeformNode.createObject('SlidingActuator', name="SlidingActuatorRate0", template='Vec3d', direction='1 0 0 ', indices=0, maxForce='100000', minForce='-30000')
    # rateAngularDeformNode.createObject('SlidingActuator', name="SlidingActuator0", template='Vec3d',direction='0 1 0 ', indices=3, maxForce='100000', minForce='-30000')
    # rateAngularDeformNode.createObject('SlidingActuator', name="SlidingActuator0", template='Vec3d',direction='0 1 0 ', indices=5, maxForce='100000', minForce='-30000')

    

    ##############
    #   Frames   #
    ##############
    frames = ["0.0 0 0  0 0 0 1   5 0 0  0 0 0 1  10.0 0 0  0 0 0 1    15.0 0 0  0 0 0 1   20.0 0 0  0 0 0 1" +
              " 30.0 0 0  0 0 0 1  35.0 0 0  0 0 0 1   40.0 0 0  0 0 0 1   45.0 0 0  0 0 0 1 55.0 0 0  0 0 0 1 60.0 0 0  0 0 0 1" +
              " 66.0 0 0  0 0 0 1   71.0 0 0  0 0 0 1   76.0 0 0  0 0 0 1  81.0 0 0  0 0 0 1"]
    # the node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.createChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)
    framesMO = mappedFrameNode.createObject(
        'MechanicalObject', template='Rigid3d', name="FramesMO", position=frames, showObject='0', showObjectScale='1')

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    curv_abs_input = '0 15 30 45 60 66 81'
    curv_abs_output = '0.0 5 10 15 20 30 35 40 45 55 60 66 71 76 81'
    mappedFrameNode.createObject('DiscretCosseratMapping', curv_abs_input=curv_abs_input,
                                 curv_abs_output=curv_abs_output, input1=inputMO, input2=inputMO_rigid, output=outputMO, debug='0')

    # actuators = mappedFrameNode.createChild('actuators')
    # actuator0 = actuators.createObject('SlidingActuator', name="SlidingActuator0", template='Rigid3d',
    #                                    direction='0 0 0 1 0 0', indices=1, maxForce='100000', minForce='-30000')
    cable_position = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0], [15.0, 0.0, 0.0], [20.0, 0.0, 0.0], [30.0, 0.0, 0.0], [35.0, 0.0, 0.0], [40.0, 0.0, 0.0], [45.0, 0.0, 0.0],
                      [55.0, 0.0, 0.0], [60.0, 0.0, 0.0], [66.0, 0.0, 0.0], [71.0, 0.0, 0.0], [76.0, 0.0, 0.0], [81.0, 0.0, 0.0]]
    #  This create a new node in the scene. This node is appended to the finger's node.
    

    ##################################################################################################################################################################################################################

    ##########################################
    # FEM Model                              #
    ##########################################
    finger = rootNode.createChild('finger')
    finger.createObject('EulerImplicitSolver', name='odesolver',firstOrder='0', rayleighMass=0.1, rayleighStiffness=0.1)
    finger.createObject('SparseLDLSolver', name='preconditioner')

    # Add a componant to load a VTK tetrahedral mesh and expose the resulting topology in the scene .
    finger.createObject('MeshVTKLoader', name='loader', filename=path + 'finger.vtk', translation="-17.5 -12.5 7.5", rotation="0 180 0")
    finger.createObject('TetrahedronSetTopologyContainer', src='@loader', name='container', poissonRatio='0.3',  youngModulus='2000')
    finger.createObject('TetrahedronSetTopologyModifier')
    finger.createObject('TetrahedronSetTopologyAlgorithms', template='Vec3d')
    finger.createObject('TetrahedronSetGeometryAlgorithms', template='Vec3d')

    # Create a mechanicaobject component to stores the DoFs of the model
    finger.createObject('MechanicalObject', name='tetras', template='Vec3d',showIndices='false', showIndicesScale='4e-5', rx='0', dz='0')
    # Gives a mass to the model
    finger.createObject('UniformMass', totalMass='0.075')

    # Add a TetrahedronFEMForceField componant which implement an elastic material model solved using the Finite Element Method on
    # tetrahedrons.
    finger.createObject('TetrahedronFEMForceField', template='Vec3d',name='FEM', method='large', poissonRatio='0.45',  youngModulus='2000')
    finger.createObject('BoxROI', name='ROI1',box='-18 -15 -8 2 -3 8', drawBoxes='true')
    finger.createObject('RestShapeSpringsForceField',points='@ROI1.indices', stiffness='1e12')

    ##########################################
    # Cable points                           #
    ##########################################
    # Mappe points inside the meca, this points will be use for the bilateral mapping
    FEMpos = [" 0.0 0 0 15 0 0 30 0 0 45 0 0 60 0 0 66 0 0 81 0.0 0.0"]
    #FEMpos = ["66. 0. 0.  81.0 0.0 0.0"]
    femPoints = finger.createChild('femPoints')
    inputFEMCable = femPoints.createObject('MechanicalObject', name="pointsInFEM", position=FEMpos, showObject="0", showIndices="0")
    femPoints.createObject('BarycentricMapping')

    ##########################################
    # Visualization                          #
    ##########################################
    fingerVisu = finger.createChild('visu')
    fingerVisu.createObject('MeshSTLLoader', filename=path+"finger.stl", name="loader", translation="-17.5 -12.5 7.5", rotation="0 180 0",)
    fingerVisu.createObject('OglModel', src="@loader", template='ExtVec3f', color="0.0 0.7 0.7")
    fingerVisu.createObject('BarycentricMapping')
    
     ##########################################
    # Visualization                          #
    ##########################################
    # In Sofa, visualization is handled by adding a rendering model.
    # Create an empty child node to store this rendering model.
    CollisionMesh(finger, surfaceMeshFileName="mesh/finger.stl", name="part0", translation="-17.5 -12.5 7.5",
                  rotation="0 180 0", collisionGroup=[1, 2])

    ##########################################
    #  Finger auto-Collision            #
    ##########################################
    CollisionMesh(finger,
                  surfaceMeshFileName="mesh/fingerCollision_part1.stl",
                  name="CollisionMeshAuto1", translation="-17.5 -12.5 7.5", rotation="0 180 0", collisionGroup=[1])

    #CollisionMesh(finger,
                  #surfaceMeshFileName="mesh/fingerCollision_part2.stl",
                  #name="CollisionMeshAuto2", translation="-17.5 -12.5 7.5", rotation="0 180 0", collisionGroup=[2])

    
    finger.createObject('LinearSolverConstraintCorrection')

    ##########################################
    # Effector                               #
    ##########################################
    effector = finger.createChild('fingertip')
    effector.createObject('MechanicalObject',position=("85.591 -5.13521 0.35857"))
    effector.createObject('PositionEffector', template='Vec3d',indices="0", effectorGoal="@../../goal/goalMO.position")
    effector.createObject('BarycentricMapping', mapForces="false", mapMasses="false")


    ##########################################
    # Multi mapped mstats                    #
    ##########################################
    slidingPoint = mappedFrameNode.createChild('slidingPoint')
    # This create a MechanicalObject, a componant holding the degree of freedom of our
    # mechanical modelling. In the case of a cable it is a set of positions specifying
    # the points where the cable is passing by.
    slidingPointMO = slidingPoint.createObject('MechanicalObject', name="cablePos", position=cable_position, showObject="0", showIndices="0")
    slidingPoint.createObject('IdentityMapping')
    mappedPointsNode = slidingPoint.createChild('MappedPoints')
    mappedPoints = mappedPointsNode.createObject('MechanicalObject', template='Vec3d', position=FEMpos, name="FramesMO", showObject='0', showObjectScale='1')
    mappedPointsNode.createObject('CosseratEquality', name="QPConstraint", eqDisp='0.0')

    ## Get the tree mstate links for the mapping
    inputCableMO = slidingPointMO.getLinkPath()
    inputFEMCableMO = inputFEMCable.getLinkPath()
    outputPointMO = mappedPoints.getLinkPath()


    femPoints.addChild(mappedPointsNode)
    # mappedPointsNode.createObject('CosseratActuator', name="QPConstraint", eqDisp='1.50', maxPositiveDisp='1.0', maxNegativeDisp='-1.', maxDispVariation="2")
    mappedPointsNode.createObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO, input2=inputCableMO, 
        output=outputPointMO, direction="@../../FramesMO.position")

    return rootNode
