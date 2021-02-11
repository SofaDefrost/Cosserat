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
        #Sofa.PythonScriptController.__init__(self, *args, **kwargs)

        self.solver   = kwargs["solver"]
        self.posEffector = kwargs["posEffector"]
        self.goalPos  = kwargs["goalPos"]
        self.effMO    = None
        self.iter     = 0        
        self.effMO = kwargs["effMO"]
        
        self.cost = np.inf
        self.best_cost = 1e-5
        self.qp_error = False
                       
        print ("************* *****************************************")        
        print(">>>>>>>>>>>>>>>>>>>>>>> CostController  <<<<<<<<<<<<<<<<<" )
        print ("************* *****************************************")

        # # th=np.linspace(475, 500, 10)
        # thetaList = np.linspace(5, 15, 20)
        # # zL        = np.linspace(90, 150, 15)
        # zL        = np.linspace(-35, -5, 20)    
        # position = []
        # i = 0
        # for theta in thetaList:
        #     x , y= computeRtheta(theta)
        #     # position.append([zL[i], x, y])  
        #     position +=[[x+170.0, y+8, zL[i]]]
        #     i+=1

        # position += [[130.72644140849795, -2.778005970892028, -3.0],[118.72644140849795, 6.778005970892028, -2.0]]

        # self.goalPos.position  = position
        # print("_________________________________________________________")
        # print ( position)
        # print("_________________________________________________________")
        
        # self.iter = len(position) -1
        # print("============W   Size  : ",len(position)) 

        # print("The goal Before : ", self.goalos.position.value)
        
        # goal = position[self.iter]
        
        # print(" The goal is : ", goal)
        # var = 90.0 # float(goal[0])
        # print (" Goal 0 : ", type(goal[0]))
        # print (" Goal 1 : ", type(0.0))
        # print (" Goal 2 : ", type(var))

        # self.goal_pos.position = [[0.0,1.0,2.0]]
        # goal = str(90.0 1.0 2.0)
        # print ("The goal is : ", goal[0])
        # self.goal_pos.position = "0.0 1.0 2.0" #[[90.,1.0,2.0]]
        # self.goal_pos.position = [[0.0,1.0,2.0]]

        self.iter = len(self.goalPos.position) -1
        goal = self.goalPos.position[self.iter]

        print( "iter :  ", self.iter)
        print( "Goal :  ", goal)
        
        self.posEffector.effectorGoal = [[goal[0],goal[1],goal[2]]]
                
    def onAnimateEndEvent(self, dt):
        
        # qp_error_messages = self.solver.getLoggedMessagesAsString()
        qp_error_messages = self.solver.getData('objective').value 
        # info = self.solver.getData('info').value.toList()
        
        if qp_error_messages and qp_error_messages != "":
            self.qp_error = True
            if qp_error_messages < self.best_cost :
                self.iter -=1
                goal = self.goalPos.position[self.iter]
                self.posEffector.effectorGoal = [[goal[0],goal[1],goal[2]]]
            print ("************* *****************************************")
            print("1 >>>>>>>>>>>>>>>>>>>>>>> Caught QP Error ", qp_error_messages)
            print ("************* *****************************************")


    def onKeypressedEvent(self, c):
        key = c['key']
        if key == "0":
            self.iter -=1
            print("You switch to : ", self.iter)
            goal = self.goalPos.position[self.iter]
            print("You switch to : ", self.iter)
            self.posEffector.effectorGoal = [[goal[0],goal[1],goal[2]]]
        if key == "1":
            print("You switch to Y axis")
            self.iter +=1
            goal = self.goalPos.position[self.iter]
            self.posEffector.effectorGoal = [[goal[0],goal[1],goal[2]]]       


def computeRtheta(theta):   
    x=constA*np.exp(constB*theta)*np.cos(theta)
    y=constA*np.exp(constB*theta)*np.sin(-theta)    
    #z = np.linspace(0,2, 10000)
    #rTeta = constA * exp(-constB * teta)
    return x , y

def createScene(rootNode):

    rootNode.addObject(
        'VisualStyle', displayFlags='showVisualModels hideBehaviorModels showCollisionModels hideBoundingCollisionModels hideForceFields showInteractionForceFields showWireframe')

    rootNode.addObject('RequiredPlugin', pluginName=["SoftRobots", "SoftRobots.Inverse", "SofaSparseSolver", "SofaPreconditioner", "SofaOpenglVisual", "CosseratPlugin", "BeamAdapter", "SofaShells"])
      

    rootNode.addObject('FreeMotionAnimationLoop')    
    # rootNode.addObject('GenericConstraintSolver', tolerance="1e-20", maxIterations="5 00", printLog="0")
    rootNode.addObject('CollisionPipeline', depth="6", verbose="0", draw="1")
    rootNode.addObject('BruteForceDetection', name="N2")
    rootNode.addObject('CollisionResponse', response="FrictionContact", responseParams="mu=0.0")
    rootNode.addObject('LocalMinDistance', name="Proximity", alarmDistance="0.6", contactDistance=0.1)

    qp_solver = rootNode.addObject('QPInverseProblemSolver', printLog='0')                                   

    rootNode.gravity = [0.0, 0.0, 0.0]
    rootNode.dt = 0.01
    
    # rootNode.gravity = "0 -9180 0"
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.addObject('OglSceneFrame', style="Arrows",
                          alignment="TopRight")

    listGoal = [[131.43185217,  48.79604369, -35.        ]
                , [136.21484585,  50.11369713, -33.42105263]
                , [141.33869767,  48.76436724, -31.84210526]
                , [145.34645072,  44.76930399, -30.26315789]
                , [146.90406087,  38.92980569, -28.68421053]
                , [145.20949465,  32.7015719,  -27.10526316]
                , [140.30310446,  27.85918702, -25.52631579]
                , [133.17424001,  26.01976351, -23.94736842]
                , [125.60412878,  28.14584607, -22.36842105]
                , [119.75408416,  34.17065867, -20.78947368]
                , [117.58417176,  42.87296571, -19.21052632]
                , [120.24923099,  52.07340331, -17.63157895]
                , [127.64645792,  59.13982565, -16.05263158]
                , [138.26873923,  61.69670049, -14.47368421]
                , [149.44981092,  58.35889736, -12.89473684]
                , [157.98433929,  49.27778283, -11.31578947]
                , [160.99352905,  36.31283813,  -9.73684211]
                , [156.81655482,  22.72573244,  -8.15789474]
                , [145.66963129,  12.41954738,  -6.57894737]
                , [129.84647457,   8.88263047,  -5.        ]
                , [111.92793183,   6.90126076,  -3.        ]]   
    ###############################################
    # Effector List goal for interactive control  #
    ###############################################
    goalList = rootNode.addChild('goalList')    
    goalListMecha = goalList.addObject('MechanicalObject', name='goalMO', position=listGoal, drawMode=1, showObject=1, showIndices=1)    
    # goalList.addObject('SphereModel', radius='1', color="blue")
    # goalList.addObject('BarycentricMapping' )
    # goalList.addObject('SkinningMapping', nbRef='1',  mapForces='false', mapMasses='false')

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
    # rateAngularDeformNode.addObject('SlidingActuator', name="SlidingActuator0", template='Vec3d',
    #                                    direction='0 1 0 ', indices=2, maxForce='100000', minForce='-30000')    
    # for i in range(1,6):
    #     rateAngularDeformNode.addObject('SlidingActuator', name="SlidingActuatorY"+str(i), template='Vec3d',
    #                                    direction='0 1 0 ', indices=i, maxForce='100000', minForce='-30000')
    for i in range(1,6):
        rateAngularDeformNode.addObject('SlidingActuator', name="SlidingActuatorZ"+str(i), template='Vec3d',
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

    # actuators = mappedFrameNode.addChild('actuators')
    # actuator0 = actuators.addObject('SlidingActuator', name="SlidingActuator0", template='Rigid3d',
    #                                    direction='0 0 0 1 0 0', indices=1, maxForce='100000', minForce='-30000')
    cable_position = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0], [15.0, 0.0, 0.0], [20.0, 0.0, 0.0], [30.0, 0.0, 0.0], [35.0, 0.0, 0.0], [40.0, 0.0, 0.0], [45.0, 0.0, 0.0],
                      [55.0, 0.0, 0.0], [60.0, 0.0, 0.0], [66.0, 0.0, 0.0], [71.0, 0.0, 0.0], [76.0, 0.0, 0.0], [81.0, 0.0, 0.0]]
    #  This create a new node in the scene. This node is appended to the finger's node.
  
    ##########################################
    # Effector                               #
    ##########################################
    effector = mappedFrameNode.addChild('fingertip') #81.0, 0.0, 0.0
    effMO  = effector.addObject('MechanicalObject',position=("83 0  0.0"))
    posEffector  = effector.addObject('PositionEffector', template='Vec3d',indices="0", effectorGoal="90.0 3.0  0.35857")
    # effector.addObject('BarycentricMapping', mapForces="false", mapMasses="false")
    effector.addObject('SkinningMapping', nbRef='1',  mapForces='false', mapMasses='false')
    ##########################################
    #  Implant collision                   #
    ##########################################
    CollisInstrumentCombined = mappedFrameNode.addChild('CollisInstrumentCombined')
    CollisInstrumentCombined.addObject('EdgeSetTopologyContainer', name="collisEdgeSet", position=cable_position, edges="0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14")
    CollisInstrumentCombined.addObject('EdgeSetTopologyModifier', name="colliseEdgeModifier")
    CollisInstrumentCombined.addObject('MechanicalObject', name="CollisionDOFs", position=cable_position)
    CollisInstrumentCombined.addObject('LineCollisionModel', bothSide="1",group='2' )
    CollisInstrumentCombined.addObject('PointCollisionModel', bothSide="1", group='2')
    CollisInstrumentCombined.addObject('SkinningMapping', name="multimapp")

    ##########################################
    #         Cochlea                        #
    ##########################################

    cochleaNode = rootNode.addChild('cochleaNode')
    cochleaNode.addObject('MeshObjLoader', name='loader', filename='mesh/cochleeCompleteTroueeSimpleOrientationMegaTrouTransformed.obj', translation="-15 5 -3", flipNormals="false") # , scale3d="10 10 10", translation="130 0 0", flipNormals="false", rotation="0 0 48"
    cochleaNode.addObject('MeshTopology',src = '@loader')
    cochleaNode.addObject('MechanicalObject', name='dofs', template='Vec3d', showIndices='false', showIndicesScale='4e-5', rx='0',printLog="0")
    cochleaNode.addObject('TriangleCollisionModel', group='1')
    cochleaNode.addObject('LineCollisionModel', group='1')
    cochleaNode.addObject('PointCollisionModel', group='1')
        
    
    visuCochleaNode = cochleaNode.addChild('visuCochleaNode')
    visuCochleaNode.addObject('OglModel', name="VisualModel", color="3.0 0.5 0.0 0.9")
    
    membraneNode = rootNode.addChild('membraneNode')
    membraneNode.addObject('EulerImplicitSolver', rayleighStiffness='0.0', rayleighMass='0.0')
    membraneNode.addObject('SparseLDLSolver')
    membraneNode.addObject('MeshObjLoader', name='loader', filename='mesh/membraneBasilaireBetterFitTransformed.obj', flipNormals="false", translation="-15 5 -3")
    membraneNode.addObject('Vertex2Frame', name="frames", template="Rigid3d.", position="@loader.position", normals="@loader.normals", invertNormals="false")
    #membraneNode.addObject('Mesh',src ='@loader')
    membraneNode.addObject('MechanicalObject', name='MO', template="Rigid3d", position="@frames.frames", showIndices="false", showIndicesScale="0.000005")
    membraneNode.addObject('TriangleSetTopologyContainer', name="coarseTopo", src="@loader")

    membraneNode.addObject('UniformMass', showAxisSizeFactor="0.1", totalMass="0.001")
    membraneNode.addObject('TriangularShellForceField', name="FEM", youngModulus="1e4", poissonRatio="0.33", rayleighStiffness="0", thickness="4.0e-1", measure="Strain (norm)")
    membraneNode.addObject('RestShapeSpringsForceField', points="7 8 10 11 12 13 14 19 20 22 23 35 36 51 52 53 56 57 58 60 61 62 64 67 80 81 82 83 89 90 91 99 100 101 114 115 116 121 122 124 126 130 131 134 152 153 159 160 164 172 173 188 189 197 198 199 201 212 213 214 215 216 217 218 220 224 225 226 234 235 236 239 243 244 245 246 248 250 254 258 259 267 271 277 284 285 287 294 296 297 298 299 300 301 302 303 310 311 317 329 330 331 332 334 336 338 339 340 341 343 347 348 350 357 359 361 364 366 367 369 370 371 372 374 376 381 383 385 386 390 396 399 403 404 405 406 407 433 435 450 453 454 456 460 464 465 470 471 472 474 478 479 480 488 491 505 506 510 512 513 514 516 520 523 525 526 529 534 535 536 554 557 558 561 562 563 564 565 567 568 569 577 586 590 592 593 598 599 600 605 610 612 616 627 630 631 632 639 647 649 650 651 652 653 654 655 662 663 668 675 687 692 693 694 695 696 697 698", stiffness="1000000", angularStiffness="100000")
    membraneNode.addObject('LinearSolverConstraintCorrection')

    
    #collMembraneNode = membraneNode.addChild('collMembraneNode')
    #collMembraneNode.addObject('MechanicalObject', name="Coll_MO", template="Vec3d", src="@../loader")
    #collMembraneNode.addObject('TriangleCollisionModel', bothSide="1", group='1')
    #collMembraneNode.addObject('LineCollisionModel', bothSide="1", group='1')
    #collMembraneNode.addObject('PointCollisionModel', bothSide="1", group='1')
    #collMembraneNode.addObject('IdentityMapping', input="@../MO", output="@Coll_MO")

    #visuMembraneNode = membraneNode.addChild('visuMembraneNode')
    #visuMembraneNode.addObject('OglModel', name="Visual", color="red")

    #controle = CostController(rootNode, name="CostController", goal_pos=goal_pos, effMO=effMO, solver=qp_solver)
    rootNode.addObject(CostController(name="CostController", goalPos=goalListMecha, posEffector=posEffector, effMO=effMO, solver=qp_solver))

    return rootNode
