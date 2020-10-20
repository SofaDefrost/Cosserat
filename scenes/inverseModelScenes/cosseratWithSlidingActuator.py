# -*- coding: utf-8 -*-

import os
import Sofa
from stlib.scene import MainHeader, ContactHeader

path = os.path.dirname(os.path.abspath(__file__))+'/mesh/'

def createScene(rootNode):

    MainHeader(rootNode, plugins=["SoftRobots", "SoftRobots.Inverse", "SofaPython", "SofaSparseSolver", "SofaPreconditioner", "SofaOpenglVisual", "CosseratPlugin", "BeamAdapter", "SofaShells"],
               repositoryPaths=[os.getcwd()])

    rootNode.createObject(
        'VisualStyle', displayFlags='showVisualModels hideBehaviorModels showCollisionModels hideBoundingCollisionModels hideForceFields showInteractionForceFields showWireframe')

    rootNode.createObject('FreeMotionAnimationLoop')
    rootNode.createObject('QPInverseProblemSolver', printLog='0')    
    # rootNode.createObject('GenericConstraintSolver', tolerance="1e-20", maxIterations="5 00", printLog="0")
    rootNode.createObject('CollisionPipeline', depth="6", verbose="0", draw="1")
    rootNode.createObject('BruteForceDetection', name="N2")
    rootNode.createObject('CollisionResponse', response="FrictionContact", responseParams="mu=0.65")
    rootNode.createObject('LocalMinDistance', name="Proximity", alarmDistance="0.6", contactDistance="0.44", angleCone="0.01")
    
    rootNode.gravity = "0 -9180 0"
    rootNode.createObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.createObject('OglSceneFrame', style="Arrows",
                          alignment="TopRight")

    ##########################################
    # Effector goal for interactive control  #
    ##########################################
    goal = rootNode.createChild('goal')
    goal.createObject('EulerImplicitSolver', firstOrder='1')
    goal.createObject('CGLinearSolver', iterations='100', tolerance="1e-5", threshold="1e-5")
    goal.createObject('MechanicalObject', name='goalMO', position='90 3  0.35857')
    goal.createObject('SphereCollisionModel', radius='5')
    goal.createObject('UncoupledConstraintCorrection')

    

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
    RigidBaseMO = rigidBaseNode.createObject('MechanicalObject', template='Rigid3d', name="RigidBaseMO", 
        position="0 0 0  0 0 0 1", showObject='1', showObjectScale='0.1')
    rigidBaseNode.createObject('RestShapeSpringsForceField', name='spring', stiffness="5",
                               angularStiffness="5", external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")
    
    rigidBaseNode.createObject('SlidingActuator', name="SlidingActuator0", template='Rigid3d',
                                       direction='1 0 0 0 0 0', indices=0, maxForce='100000', minForce='-30000')
    # rigidBaseNode.createObject('SlidingActuator', name="SlidingActuator0", template='Rigid3d',
    #                                    direction='0 1 0 0 0 0', indices=0, maxForce='100000', minForce='-30000')

    ###############
    # Rate of angular Deformation  (2 sections)
    ###############
    position = ["0 0 0 " + "0 0 0 " + "0 0 0 " +
                "0 0 0 " + "0 0 0 " + "0 0 0 "]
    longeur = '15 15 15 15 6 15'  # beams size
    rateAngularDeformNode = cableNode.createChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.createObject(
        'MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=position)
    BeamHookeLawForce = rateAngularDeformNode.createObject(
        'BeamHookeLawForceField', crossSectionShape='circular', length=longeur, radius='0.5', youngModulus='5e6')
    rateAngularDeformNode.createObject('SlidingActuator', name="SlidingActuator0", template='Vec3d',
                                       direction='0 1 0 ', indices=2, maxForce='100000', minForce='-30000')
    # rateAngularDeformNode.createObject('SlidingActuator', name="SlidingActuator0", template='Vec3d',
    #                                    direction='0 1 0 ', indices=3, maxForce='100000', minForce='-30000')
    # rateAngularDeformNode.createObject('SlidingActuator', name="SlidingActuator0", template='Vec3d',
    #                                    direction='0 1 0 ', indices=5, maxForce='100000', minForce='-30000')
    for i in range(1,6):
        rateAngularDeformNode.createObject('SlidingActuator', name="SlidingActuator"+str(i), template='Vec3d',
                                       direction='0 1 0 ', indices=i, maxForce='100000', minForce='-30000')
    for i in range(1,6):
        rateAngularDeformNode.createObject('SlidingActuator', name="SlidingActuator"+str(i), template='Vec3d',
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
    mappedFrameNode = rigidBaseNode.createChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)
    framesMO = mappedFrameNode.createObject(
        'MechanicalObject', template='Rigid3d', name="FramesMO", position=frames, showObject='1', showObjectScale='1')

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
  
    ##########################################
    # Effector                               #
    ##########################################
    effector = mappedFrameNode.createChild('fingertip')
    effector.createObject('MechanicalObject',position=("89 3  0.35857"))
    effector.createObject('PositionEffector', template='Vec3d',indices="0", effectorGoal="@../../../../goal/goalMO.position")
    # effector.createObject('BarycentricMapping', mapForces="false", mapMasses="false")
    effector.createObject('SkinningMapping', nbRef='1',  mapForces='false', mapMasses='false')


    # slidingPoint = mappedFrameNode.createChild('slidingPoint')

    # # This create a MechanicalObject, a componant holding the degree of freedom of our
    # # mechanical modelling. In the case of a cable it is a set of positions specifying
    # # the points where the cable is passing by.
    # slidingPointMO = slidingPoint.createObject('MechanicalObject', name="cablePos",
    #                           position=cable_position, showObject="1", showIndices="1")
    # slidingPoint.createObject('IdentityMapping')


    # mappedPointsNode = slidingPoint.createChild('MappedPoints')
    # femPoints.addChild(mappedPointsNode)
    # mappedPoints = mappedPointsNode.createObject('MechanicalObject', template='Vec3d', position=FEMpos, name="FramesMO", showObject='1', showObjectScale='1')

    # ## Get the tree mstate links for the mapping
    # inputCableMO = slidingPointMO.getLinkPath()
    # inputFEMCableMO = inputFEMCable.getLinkPath()
    # outputPointMO = mappedPoints.getLinkPath()

    # mappedPointsNode.createObject('CosseratActuator', nodeame="QPConstraint")

    # mappedPointsNode.createObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO, input2=inputCableMO, output=outputPointMO)

    ##########################################
    #         Cochlea                        #
    ##########################################

    cochleaNode = rootNode.createChild('cochleaNode')
    cochleaNode.createObject('MeshObjLoader', name='loader', filename='mesh/cochleeCompleteTroueeSimpleOrientationMegaTrou.obj', flipNormals="false", scale3d="10 10 10", translation="130 0 0")
    cochleaNode.createObject('MeshTopology',src = '@loader')
    cochleaNode.createObject('MechanicalObject', name='dofs', template='Vec3d', showIndices='false', showIndicesScale='4e-5', rx='0',printLog="0")
    cochleaNode.createObject('TriangleCollisionModel', group='1')
    cochleaNode.createObject('LineCollisionModel', group='1')
    cochleaNode.createObject('PointCollisionModel', group='1')
    
    visuCochleaNode = cochleaNode.createChild('visuCochleaNode')
    visuCochleaNode.createObject('OglModel', name="VisualModel", color="3.0 0.5 0.0 0.9")
    
    membraneNode = rootNode.createChild('membraneNode')
    membraneNode.createObject('EulerImplicitSolver', rayleighStiffness='0.0', rayleighMass='0.0')
    membraneNode.createObject('SparseLDLSolver')
    membraneNode.createObject('MeshObjLoader', name='loader', filename='mesh/membraneBasilaireBetterFit.obj', flipNormals="false", scale3d="10 10 10", translation="130 0 0")
    membraneNode.createObject('Vertex2Frame', name="frames", template="Rigid3d.", position="@loader.position", normals="@loader.normals", invertNormals="false")
    #membraneNode.createObject('Mesh',src ='@loader')
    membraneNode.createObject('MechanicalObject', name='MO', template="Rigid3d", position="@frames.frames", showIndices="false", showIndicesScale="0.000005")
    membraneNode.createObject('TriangleSetTopologyContainer', name="coarseTopo", src="@loader")

    membraneNode.createObject('UniformMass', showAxisSizeFactor="0.1", totalMass="0.001")
    membraneNode.createObject('TriangularShellForceField', name="FEM", youngModulus="1e4", poissonRatio="0.33", rayleighStiffness="0", thickness="4.0e-1", measure="Strain (norm)")
    membraneNode.createObject('RestShapeSpringsForceField', points="7 8 10 11 12 13 14 19 20 22 23 35 36 51 52 53 56 57 58 60 61 62 64 67 80 81 82 83 89 90 91 99 100 101 114 115 116 121 122 124 126 130 131 134 152 153 159 160 164 172 173 188 189 197 198 199 201 212 213 214 215 216 217 218 220 224 225 226 234 235 236 239 243 244 245 246 248 250 254 258 259 267 271 277 284 285 287 294 296 297 298 299 300 301 302 303 310 311 317 329 330 331 332 334 336 338 339 340 341 343 347 348 350 357 359 361 364 366 367 369 370 371 372 374 376 381 383 385 386 390 396 399 403 404 405 406 407 433 435 450 453 454 456 460 464 465 470 471 472 474 478 479 480 488 491 505 506 510 512 513 514 516 520 523 525 526 529 534 535 536 554 557 558 561 562 563 564 565 567 568 569 577 586 590 592 593 598 599 600 605 610 612 616 627 630 631 632 639 647 649 650 651 652 653 654 655 662 663 668 675 687 692 693 694 695 696 697 698", stiffness="1000000", angularStiffness="100000")
    membraneNode.createObject('LinearSolverConstraintCorrection')

    
    collMembraneNode = membraneNode.createChild('collMembraneNode')
    collMembraneNode.createObject('MechanicalObject', name="Coll_MO", template="Vec3d", src="@../loader")
    collMembraneNode.createObject('TriangleCollisionModel', bothSide="1", group='1')
    collMembraneNode.createObject('LineCollisionModel', bothSide="1", group='1')
    collMembraneNode.createObject('PointCollisionModel', bothSide="1", group='1')
    collMembraneNode.createObject('IdentityMapping', input="@../MO", output="@Coll_MO")

    visuMembraneNode = membraneNode.createChild('visuMembraneNode')
    visuMembraneNode.createObject('OglModel', name="Visual", color="red")

    return rootNode
