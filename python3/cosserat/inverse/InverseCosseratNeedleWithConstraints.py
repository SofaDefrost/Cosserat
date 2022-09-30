# -*- coding: utf-8 -*-
"""Basic scene using Cosserat in SofaPython3.

Based on the work done with SofaPython. See POEMapping.py.
"""

# from usefulFunctions import MoveTargetProcess
from cosserat.usefulFunctions import buildEdges, pluginList, BuildCosseratGeometry, MoveTargetProcess
from createFemRegularGrid import createFemCube
__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com, yinoussa.adagolodjo@inria.fr"
__version__ = "1.0.0"
__copyright__ = "(c) 2020,Inria"
__date__ = "March 16 2021"

import os
import sys

sys.path.append('../')


def effectorTarget(parentNode, position=[80., 0., 0.35857]):
    target = parentNode.createChild("Target")
    target.addObject("EulerImplicitSolver", firstOrder=True)
    target.addObject(
        'VisualStyle', displayFlags='showVisualModels showInteractionForceFields showCollisionModels')
    target.addObject("CGLinearSolver")
    targetMO = target.addObject("MechanicalObject", name="dofs", position=position,
                                showObject=0, showObjectScale=0, drawMode=0, showColor=[1., 1., 1., 1.])
    target.addObject('SphereCollisionModel', radius='2')
    target.addObject("UncoupledConstraintCorrection")
    return [target, targetMO]


def createScene(rootNode):
    from stlib3.scene import MainHeader
    MainHeader(rootNode, plugins=pluginList +
               ["SoftRobots.Inverse", "Sofa.Component.LinearSolver.Iterative", "Sofa.Component.Collision.Response.Contact",
                "Sofa.Component.Collision.Geometry", "Sofa.Component.Collision.Detection.Intersection", "Sofa.Component.Collision.Detection.Algorithm"],
               repositoryPaths=[os.getcwd()])
    rootNode.addObject('VisualStyle',
                       displayFlags='showVisualModels hideBehaviorModels showCollisionModels '
                                    'hideBoundingCollisionModels hideForceFields showInteractionForceFields '
                                    'showWireframe')

    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.createObject('QPInverseProblemSolver',
                          printLog='0', epsilon=0.0002)
    rootNode.addObject('CollisionPipeline', depth="6", verbose="0", draw="1")
    rootNode.addObject('BruteForceDetection', name="N2")
    rootNode.addObject('DefaultContactManager',
                       response="FrictionContact", responseParams="mu=0.65")
    rootNode.addObject('LocalMinDistance', name="Proximity", alarmDistance="0.6", contactDistance="0.44",
                       angleCone="0.01")

    rootNode.gravity = [0.0, 0.0, 0.0]
    rootNode.dt = 0.01

    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.addObject('OglSceneFrame', style="Arrows", alignment="TopRight")

    # #######################################
    # New adds to use the sliding Actuator  #
    #########################################
    cableNode = rootNode.addChild('cableNode')
    cableNode.addObject('EulerImplicitSolver', firstOrder="0",
                        rayleighStiffness="1.0", rayleighMass='0.1')
    cableNode.addObject('SparseLUSolver', name='solver')
    cableNode.addObject('GenericConstraintCorrection')

    #################################
    ##           RigidBase         ##
    #################################
    rigidBaseNode = cableNode.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode.addObject('MechanicalObject', template='Rigid3d', name="RigidBaseMO",
                                          position="0 0 0  0 0 0 1", showObject='1', showObjectScale='0.1')
    rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="5",
                            angularStiffness="5", external_points="0", mstate="@RigidBaseMO", points="0",
                            template="Rigid3d")

    #################################################################
    ##  Sliding actuator to guide the base of the needle, only the  ##
    # translation are tacking into account here.
    #################################################################
    for j in range(6):
        direction = [0, 0, 0, 0, 0, 0]
        direction[j] = 1
        rigidBaseNode.addObject('SlidingActuator', name=f"SlidingActuator{str(j)}", template='Rigid3d',
                                direction=direction, indices=0, maxForce='1e6', minForce='-30000')

    cosserat_config = {'init_pos': [0., 0., 0.], 'tot_length': 80, 'nbSectionS': 8,
                       'nbFramesF': 16, 'beamMass': 0.22}

    [positionS, curv_abs_inputS, longeurS, framesF, curv_abs_outputF, cable_positionF] = \
        BuildCosseratGeometry(cosserat_config)

    #############################################
    # Rate of angular Deformation  (2 sections)
    #############################################
    rateAngularDeformNode = cableNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject(
        'MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=positionS)
    rateAngularDeformNode.addObject(
        'BeamHookeLawForceField', crossSectionShape='circular', length=longeurS, radius='0.5', youngModulus='1.0e12')

    # for i in range(0, nbSectionS):
    #     rateAngularDeformNode.addObject('SlidingActuator', name="SlidingActuatorY" + str(i), template='Vec3d',
    #                                     direction='0 1 0 ', indices=i, maxForce='100000', minForce='-30000')
    # for i in range(0, nbSectionS):
    #     rateAngularDeformNode.addObject('SlidingActuator', name="SlidingActuatorZ" + str(i), template='Vec3d',
    #                                     direction='0 0 1', indices=i, maxForce='100000', minForce='-30000')

    ##############
    #   Frames   #
    ##############
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

    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=curv_abs_inputS, output=outputMO,
                              curv_abs_output=curv_abs_outputF, input1=inputMO, input2=inputMO_rigid,
                              debug='0', max=6.e-2, deformationAxis=1, nonColored="0", radius=5)

    # #Mapped mechanical points
    #  This create a new node in the scene. This node is appended to the finger's node.
    framePoints = mappedFrameNode.addChild('framePoints')
    framePointsMO = framePoints.addObject('MechanicalObject', name="cablePos", position=cable_positionF,
                                          showObject="1", showIndices="0")
    framePoints.addObject('IdentityMapping')

    ##########################################
    # ## create FEM grid
    ##########################################
    cubeNode = createFemCube(rootNode)
    ####################################################################################
    # Attached the target to FEM model
    # targetNode = cubeNode.gelNode.addChild('target')
    # targetNode.addObject('VisualStyle', displayFlags='showCollisionModels ')
    # targetMO = targetNode.addObject('MechanicalObject', name="targetMO", position='85.0 1.5  0.35857', showObject="1",
    #                                 showIndices="1", template="Vec3d")
    # targetNode.addObject('SphereCollisionModel', radius='2', color="red")
    # targetNode.addObject('BarycentricMapping')
    [targetNode, targetMO] = effectorTarget(rootNode)
    ####################################################################################
    # Create constraint Points inside the volume of the deformable object
    # These points are created
    femPos = [" 41.0 0 0  45 0 0 50 0 0 55 0 0 "]
    gelNode = cubeNode.getChild('gelNode')
    femPoints = gelNode.addChild('femPoints')
    inputFEMCable = femPoints.addObject(
        'MechanicalObject', name="pointsInFEM", position=femPos, showIndices="1")
    femPoints.addObject('BarycentricMapping')

    ####################################################################################
    # Effector constraint,
    # first define the tip on then needle then attach the effector constraint which
    # leads the tip of the needle to the defined target #
    ####################################################################################
    effector = mappedFrameNode.addChild('fingertip')
    effMO = effector.addObject('MechanicalObject', position=[80., 0., 0.35857])
    posEffector = effector.addObject('PositionEffector', template='Vec3d', indices="0",
                                     effectorGoal=targetMO.getLinkPath()+'.position')  # "@../../../../FemNode/gelNode/target/targetMO.position"
    targetNode.addObject(MoveTargetProcess(posEffector))
    effector.addObject('SkinningMapping', nbRef='1',
                       mapForces='false', mapMasses='false')

    mappedPointsNode = framePoints.addChild('MappedPoints')
    mappedPoints = mappedPointsNode.addObject('MechanicalObject', template='Vec3d', position=femPos, name="FramesMO",
                                              showObject='1', showIndices='1', showObjectScale='1')
    mappedPointsNode.addObject(
        'CosseratEquality', name="QPConstraint", eqDisp='0.0')

    # Get the tree mstate links for the mapping
    inputCableMO = framePointsMO.getLinkPath()
    inputFEMCableMO = inputFEMCable.getLinkPath()
    outputPointMO = mappedPoints.getLinkPath()

    femPoints.addChild(mappedPointsNode)

    mappedPointsNode.addObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO, lastPointIsFixed=0,
                               input2=inputCableMO, output=outputPointMO, direction="@../../FramesMO.position")

    # rootNode.addObject(CostController(name="CostController", goal_pos=goalPos, effMO=effMO, solver=qp_solver))

    # return rootNode
