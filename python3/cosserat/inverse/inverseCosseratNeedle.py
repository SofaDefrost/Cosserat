# -*- coding: utf-8 -*-
"""Basic scene using Cosserat in SofaPython3.

Based on the work done with SofaPython. See POEMapping.pyscn.
"""

__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com, yinoussa.adagolodjo@inria.fr"
__version__ = "1.0.0"
__copyright__ = "(c) 2020,Inria"
__date__ = "March 16 2021"

import os
import sys
from cosserat.usefulFunctions import buildEdges, pluginList, BuildCosseratGeometry
sys.path.append('../')
path = f'{os.path.dirname(os.path.abspath(__file__))}/../mesh/'


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
    qp_solver = rootNode.addObject('QPInverseProblemSolver', printLog='0')
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

    ##########################################
    # Effector goal for interactive control  #
    ##########################################
    goal = rootNode.addChild('goal')
    goal.addObject('EulerImplicitSolver', firstOrder='1')
    goal.addObject('CGLinearSolver', iterations='100',
                   tolerance="1e-5", threshold="1e-5")
    goalPos = goal.addObject(
        'MechanicalObject', name='goalMO', position='90.0 3.0  0.35857')
    goal.addObject('SphereCollisionModel', radius='2', color="red")
    goal.addObject('UncoupledConstraintCorrection')

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
    for j in range(3):
        direction = [0, 0, 0, 0, 0, 0]
        direction[j] = 1
        rigidBaseNode.addObject('SlidingActuator', name=f'SlidingActuator{str(j)}', template='Rigid3d',
                                direction=direction, indices=0, maxForce='100000', minForce='-30000')

    #############################################
    # Rate of angular Deformation  (2 sections)
    #############################################

    # Define: the number of section, the total lenght and the lenght of each beam.
    nbSectionS = 6
    totalLength = 80.0
    lengthS = totalLength / nbSectionS

    # Define: the length of each beam in a list, the positions of eahc beam
    # (flexion, torsion), the abs of each section
    positionS = []
    longeurS = []
    temp = 0.
    curv_abs_inputS = [0.0]
    for i in range(nbSectionS):
        positionS.append([0, 0, 0])
        longeurS.append((((i + 1) * lengthS) - i * lengthS))
        temp += longeurS[i]
        curv_abs_inputS.append(temp)
    curv_abs_inputS[nbSectionS] = totalLength

    rateAngularDeformNode = cableNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject(
        'MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=positionS)
    BeamHookeLawForce = rateAngularDeformNode.addObject(
        'BeamHookeLawForceField', crossSectionShape='circular', length=longeurS, radius='0.5', youngModulus='5e6')

    for i in range(1, 6):
        rateAngularDeformNode.addObject('SlidingActuator', name=f"SlidingActuatory{str(i)}", template='Vec3d',
                                        direction='0 1 0 ', indices=i, maxForce='100000', minForce='-30000')
    for i in range(1, 6):
        rateAngularDeformNode.addObject('SlidingActuator', name=f"SlidingActuatorz{str(i)}", template='Vec3d',
                                        direction='0 0 1', indices=i, maxForce='100000', minForce='-30000')

    ##############
    #   Frames   #
    ##############
    # Define: the number of frame and the length between each frame.
    nbFramesF = 16
    lengthF = totalLength / nbFramesF

    # Define: the abs of each frame and the position of each frame.
    framesF = []
    curv_abs_outputF = []
    cable_positionF = []
    for i in range(nbFramesF):
        sol = i * lengthF
        framesF.append([sol, 0, 0, 0, 0, 0, 1])
        cable_positionF.append([sol, 0, 0])
        curv_abs_outputF.append(sol)

    framesF.append([totalLength, 0, 0, 0, 0, 0, 1])
    cable_positionF.append([totalLength, 0, 0])
    curv_abs_outputF.append(totalLength)

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
    slidingPoint = mappedFrameNode.addChild('slidingPoint')
    slidingPoint.addObject('MechanicalObject', name="cablePos", position=cable_positionF, showObject="1",
                           showIndices="0")
    slidingPoint.addObject('IdentityMapping')

    ##########################################
    # Effector                               #
    ##########################################
    effector = mappedFrameNode.addChild('fingertip')
    effMO = effector.addObject('MechanicalObject', position=[81., 0., 0.35857])
    effector.addObject('PositionEffector', template='Vec3d', indices="0",
                       effectorGoal="@../../../../goal/goalMO.position")
    effector.addObject('SkinningMapping', nbRef='1',
                       mapForces='false', mapMasses='false')

    return rootNode
