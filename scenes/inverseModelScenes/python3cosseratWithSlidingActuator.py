# -*- coding: utf-8 -*-

import os
import Sofa
# from stlib.scene import MainHeader, ContactHeader

from math import sin, cos, exp
import os
import numpy as np

path = os.path.dirname(os.path.abspath(__file__))+'/mesh/'

constA = 4.9
constB = 0.125

class CostController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        # These are needed (and the normal way to override from a python class)
        Sofa.Core.Controller.__init__(self, *args, **kwargs)        

        self.solver   = kwargs.get("solver", None)
        self.goal_pos = kwargs.get("goal_pos",None)
        # self.goalPos  = kwargs.get("goalPos",None)
        self.iter     = 0
        self.effMO = kwargs.get("effMO",None)

        
        self.cost = np.inf
        self.best_cost = 1e-5
        self.qp_error = False
                       
        # print ("************* *****************************************")        
        # print(">>>>>>>>>>>>>>>>>>>>>>> CostController  <<<<<<<<<<<<<<<<<" )
        # print ("************* *****************************************")
        
        thetaList = np.linspace(5, 15, 15)
        zL        = np.linspace(90, 150, 15)
        # zL        = np.linspace(-35, -5, 15)    
        position = []
        i = 0
        for theta in thetaList:
            x , y= computeRtheta(theta)
            position.append([zL[i], x, y])  
            #position +=[[x+170.0, y+8, zL[i]]]
            i+=1

        
        # self.goalPos.position  = position        
        # print("The goal Before : ", self.goal_pos.position.value)
        
        goal = position[self.iter]
        var = 10.0 # float(goal[0])
        print (" Goal 0 : ", type(goal[0]))
        print (" Goal 1 : ", type(0.0))
        print (" Goal 2 : ", type(var))
        print ("The goal is : ", goal[0])

        # self.goal_pos.position = [[0.0,1.0,2.0]]
        # goal = str(90.0 1.0 2.0)
        # self.goal_pos.position = "0.0 1.0 2.0" #[[90.,1.0,2.0]]
        # self.goal_pos.position = [[0.0,1.0,2.0]]
        
        self.goal_pos.position = [[goal[0],goal[1],goal[2]]]

        # print("The goal After XXXXXXXXXXXXXXXXX: ", self.goal_pos.position.value)


    # def onAnimateEndEvent(self, dt):

    #     qp_error_messages = self.solver.objective.value

    #     if qp_error_messages and qp_error_messages != "":
    #         self.qp_error = True
    #         if qp_error_messages < self.best_cost :
    #             self.iter +=1
    #             goal = self.goalPos.value[self.iter]
    #             self.goal_pos.position = [[goal[0],goal[1],goal[2]]]
            # print ("************* *****************************************")
            # print("1 >>>>>>>>>>>>>>>>>>>>>>> Caught QP Error ", qp_error_messages)
            # print ("************* *****************************************")

        # if not self.qp_error:
        #     self.cost = get_distance_to_goal(self.effMO.position.value[0], self.goal_pos.position.value[0])            
        #     if self.cost < self.best_cost:
        #         self.best_cost = self.cost

    def onKeypressedEvent(self, c):
        key = c['key']
        if key == "0":
            print("You switch to X axis")
            self.iter -=1
            # goal = self.goalPos.position[self.iter]
            # self.goal_pos.position = [[goal[0],goal[1],goal[2]]]
        if key == "1":
            print("You switch to Y axis")
            self.iter +=1
            # goal = self.goalPos.position[self.iter]
            # self.goal_pos.position = [[goal[0],goal[1],goal[2]]]       


def get_distance_to_goal(kt_endpoint, goal_pos):
    print ("1 - ************* *****************************************")
    print("inside the function get_distance_to_goal")
    print("goal_pos : " ,goal_pos[:3])
    print("endpoint : " , kt_endpoint[:3])

    result = np.linalg.norm(goal_pos[:3]-kt_endpoint[:3])
    print ("The result is ", result)
    print ("2 - ************* *****************************************")
    return result

def computeRtheta(theta):   
    x=constA*np.exp(constB*theta)*np.cos(theta)
    y=constA*np.exp(constB*theta)*np.sin(-theta)    
    #z = np.linspace(0,2, 10000)
    #rTeta = constA * exp(-constB * teta)
    return x , y

def createScene(rootNode):

    # MainHeader(rootNode, plugins=["SoftRobots", "SoftRobots.Inverse", "SofaPython", "SofaSparseSolver", "SofaPreconditioner", "SofaOpenglVisual", "CosseratPlugin", "BeamAdapter", "SofaShells"],
    #            repositoryPaths=[os.getcwd()])
    rootNode.addObject(
        'VisualStyle', displayFlags='showVisualModels hideBehaviorModels showCollisionModels hideBoundingCollisionModels hideForceFields showInteractionForceFields showWireframe')

    rootNode.addObject('RequiredPlugin', pluginName=["SoftRobots", "SoftRobots.Inverse", "SofaSparseSolver", "SofaPreconditioner", "SofaOpenglVisual", "CosseratPlugin", "BeamAdapter"])
    # rootNode.addObject('RequiredPlugin', pluginName='SoftRobots')
    # rootNode.addObject('RequiredPlugin', pluginName='SoftRobots.Inverse')
    # rootNode.addObject('RequiredPlugin', pluginName='SofaSparseSolver')

    # rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings '
    #                                            'showForceFields')


    

    rootNode.addObject('FreeMotionAnimationLoop')    
    # qp_solver = rootNode.addObject('QPInverseProblemSolver', epsilon=0.0, printLog=False, displayTime=0, tolerance=1e-10,
    #                            maxIterations=10000) 
    qp_solver = rootNode.createObject('QPInverseProblemSolver', printLog='0')                                   
    # rootNode.addObject('GenericConstraintSolver', tolerance="1e-20", maxIterations="5 00", printLog="0")
    rootNode.addObject('CollisionPipeline', depth="6", verbose="0", draw="1")
    rootNode.addObject('BruteForceDetection', name="N2")
    rootNode.addObject('CollisionResponse', response="FrictionContact", responseParams="mu=0.65")
    rootNode.addObject('LocalMinDistance', name="Proximity", alarmDistance="0.6", contactDistance="0.44", angleCone="0.01")

    rootNode.gravity = [0.0, 0.0, 0.0]
    rootNode.dt = 0.01
    
    # rootNode.gravity = "0 -9180 0"
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.addObject('OglSceneFrame', style="Arrows",
                          alignment="TopRight")

    ##########################################
    # Effector goal for interactive control  #
    ##########################################
    goal = rootNode.addChild('goal')
    goal.addObject('EulerImplicitSolver', firstOrder='1')
    goal.addObject('CGLinearSolver', iterations='100', tolerance="1e-5", threshold="1e-5")
    goal_pos = goal.addObject('MechanicalObject', name='goalMO', position='90.0 3.0  0.35857')
    goal.addObject('SphereCollisionModel', radius='2', color="red")
    goal.addObject('UncoupledConstraintCorrection')
        

    # #######################################
    # New adds to use the sliding Actuator  #
    #########################################
    cableNode = rootNode.addChild('cableNode')
    cableNode.addObject('EulerImplicitSolver', firstOrder="0", rayleighStiffness="1.0", rayleighMass='0.1')
    cableNode.addObject('SparseLUSolver', name='solver')
    cableNode.addObject('GenericConstraintCorrection')

    # ###############
    # RigidBase
    ###############
    rigidBaseNode = cableNode.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode.addObject('MechanicalObject', template='Rigid3d', name="RigidBaseMO", 
        position="0 0 0  0 0 0 1", showObject='1', showObjectScale='0.1')
    rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="5",
                               angularStiffness="5", external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")
    
    for j in range(0,3):
        direction = [0,0,0,0,0,0]      
        direction[j] = 1  
        rigidBaseNode.addObject('SlidingActuator', name="SlidingActuator"+str(j), template='Rigid3d',
                                       direction=direction, indices=0, maxForce='100000', minForce='-30000')
    
    ###############
    # Rate of angular Deformation  (2 sections)
    ###############
    position = ["0 0 0 " + "0 0 0 " + "0 0 0 " +
                "0 0 0 " + "0 0 0 " + "0 0 0 "]
    longeur = '15 15 15 15 6 15'  # beams size
    rateAngularDeformNode = cableNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject(
        'MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=position)
    BeamHookeLawForce = rateAngularDeformNode.addObject(
        'BeamHookeLawForceField', crossSectionShape='circular', length=longeur, radius='0.5', youngModulus='5e6')
    
    for i in range(1,6):
        rateAngularDeformNode.addObject('SlidingActuator', name="SlidingActuatory"+str(i), template='Vec3d',
                                       direction='0 1 0 ', indices=i, maxForce='100000', minForce='-30000')
    for i in range(1,6):
        rateAngularDeformNode.addObject('SlidingActuator', name="SlidingActuatorz"+str(i), template='Vec3d',
                                       direction='0 0 1', indices=i, maxForce='100000', minForce='-30000')

    ################################
    # Animation (to move the dofs) #
    ################################
    # anim = Animation(rigidBaseNode, rateAngularDeformNode)

    ##############
    #   Frames   #
    ##############
    frames = ["0.0 0 0  0 0 0 1   5 0 0  0 0 0 1  10.0 0 0  0 0 0 1    15.0 0 0  0 0 0 1   20.0 0 0  0 0 0 1" +
              " 30.0 0 0  0 0 0 1  35.0 0 0  0 0 0 1   40.0 0 0  0 0 0 1   45.0 0 0  0 0 0 1 55.0 0 0  0 0 0 1 60.0 0 0  0 0 0 1" +
              " 66.0 0 0  0 0 0 1   71.0 0 0  0 0 0 1   76.0 0 0  0 0 0 1  81.0 0 0  0 0 0 1"]
    # the node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.addChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)
    framesMO = mappedFrameNode.addObject(
        'MechanicalObject', template='Rigid3d', name="FramesMO", position=frames, showObject='1', showObjectScale='1')

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    curv_abs_input = '0 15 30 45 60 66 81'
    curv_abs_output = '0.0 5 10 15 20 30 35 40 45 55 60 66 71 76 81'
    mappedFrameNode.addObject('DiscretCosseratMapping', curv_abs_input=curv_abs_input,
                                 curv_abs_output=curv_abs_output, input1=inputMO, input2=inputMO_rigid, output=outputMO, debug='0')

    cable_position = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0], [15.0, 0.0, 0.0], [20.0, 0.0, 0.0], [30.0, 0.0, 0.0], [35.0, 0.0, 0.0], [40.0, 0.0, 0.0], [45.0, 0.0, 0.0],
                      [55.0, 0.0, 0.0], [60.0, 0.0, 0.0], [66.0, 0.0, 0.0], [71.0, 0.0, 0.0], [76.0, 0.0, 0.0], [81.0, 0.0, 0.0]]
    #  This create a new node in the scene. This node is appended to the finger's node.
  
    ##########################################
    # Effector                               #
    ##########################################
    effector = mappedFrameNode.addChild('fingertip')
    effMO  = effector.addObject('MechanicalObject',position=("89 3  0.35857"))
    effector.addObject('PositionEffector', template='Vec3d',indices="0", effectorGoal="@../../../../goal/goalMO.position")
    # effector.addObject('BarycentricMapping', mapForces="false", mapMasses="false")
    effector.addObject('SkinningMapping', nbRef='1',  mapForces='false', mapMasses='false')
    rootNode.addObject(CostController(name="CostController", goal_pos=goal_pos, effMO=effMO, solver=qp_solver))

    return rootNode
