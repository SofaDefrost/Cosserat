# -*- coding: utf-8 -*-
"""Basic scene using Cosserat in SofaPython3.

Based on the work done with SofaPython. See POEMapping.pyscn.
"""

__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com, yinoussa.adagolodjo@inria.fr"
__version__ = "1.0.0"
__copyright__ = "(c) 2020,Inria"
__date__ = "March 8 2020"

import SofaRuntime
import Sofa
import os
from splib3.numerics import Quat

# from stlib3.scene import Node
path = os.path.dirname(os.path.abspath(__file__))+'/../mesh/'


class Animation(Sofa.Core.Controller):
    """
        Implements the AnimationManager as a PythonScriptController
    """
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.rigidBaseMO = args[0]
        self.rateAngularDeformMO = args[1]

        self.rate = 0.2
        self.angularRate = 0.02
        return

    def onKeypressedEvent(self, event):
        key = event['key']
        if ord(key) == 19:  # up
            with self.rigidBaseMO.rest_position.writeable() as posA:
                qOld = Quat()
                for i in range(0, 4):
                    qOld[i] = posA[0][i+3]

                qNew = Quat.createFromEuler([0., self.angularRate, 0.], 'ryxz')
                qNew.normalize()
                qNew.rotateFromQuat(qOld)
                for i in range(0, 4):
                    posA[0][i+3] = qNew[i]

        if ord(key) == 21:  # down
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


def createScene(rootNode):

    from stlib3.scene import MainHeader
    # from stlib3.physics.collision import CollisionMesh

    MainHeader(rootNode, plugins=["SoftRobots", "SofaSparseSolver", 'SofaDeformable', 'SofaEngine',
                                  'SofaImplicitOdeSolver', 'SofaLoader', 'SofaSimpleFem', "SofaPreconditioner",
                                  "SofaOpenglVisual", "CosseratPlugin", "BeamAdapter", "SofaConstraint"],
               repositoryPaths=[os.getcwd()])
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showInteractionForceFields showWireframe')

    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('GenericConstraintSolver', tolerance="1e-20", maxIterations="500", printLog="0")
    # ContactHeader(rootNode, alarmDistance=2.5, contactDistance=2, frictionCoef=0.08)

    gravity = [0, 0, 0]
    rootNode.gravity.value = gravity
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.addObject('OglSceneFrame', style="Arrows", alignment="TopRight")

    ##########################################
    # FEM Model                              #
    ##########################################
    finger = rootNode.addChild('finger')
    finger.addObject('EulerImplicitSolver', name='odesolver', firstOrder='0', rayleighMass=0.1, rayleighStiffness=0.1)
    finger.addObject('SparseLDLSolver', name='preconditioner')

    # Add a componant to load a VTK tetrahedral mesh and expose the resulting topology in the scene .
    finger.addObject('MeshVTKLoader', name='loader', filename=path + 'finger.vtk', translation="-17.5 -12.5 7.5",
                     rotation="0 180 0")
    finger.addObject('TetrahedronSetTopologyContainer', src='@loader', name='container')
    finger.addObject('TetrahedronSetTopologyModifier')
    # Create a mechanicaobject component to stores the DoFs of the model
    finger.addObject('MechanicalObject', name='tetras', template='Vec3d', showIndices='false', showIndicesScale='4e-5',
                     rx='0', dz='0')
    # Gives a mass to the model
    finger.addObject('UniformMass', totalMass='0.075')
    # Add a TetrahedronFEMForceField component which implement an elastic material model
    # solved using the Finite Element Method on
    # tetrahedrons.
    finger.addObject('TetrahedronFEMForceField', template='Vec3d',
                        name='FEM', method='large', poissonRatio='0.45',  youngModulus='500')

    finger.addObject('BoxROI', name='ROI1', box='-18 -15 -8 2 -3 8', drawBoxes='true')
    finger.addObject('RestShapeSpringsForceField', points='@ROI1.indices', stiffness='1e12')

    ##########################################
    # Cable points                           #
    ##########################################
    # Mappe points inside the meca, this points will be use for the bilateral mapping
    FEMpos = [" 0.0 0 0 15 0 0 30 0 0 45 0 0 60 0 0 66 0 0 81 0.0 0.0"]

    femPoints = finger.addChild('femPoints')
    inputFEMCable = femPoints.addObject('MechanicalObject', name="pointsInFEM", position=FEMpos, showObject="1",
                                        showIndices="1")
    femPoints.addObject('BarycentricMapping')

    ##########################################
    #  Finger auto-Collision            #
    ##########################################
    # CollisionMesh(finger,
    #               surfaceMeshFileName="mesh/fingerCollision_part1.stl",
    #               name="CollisionMeshAuto1", translation="-17.5 -12.5 7.5", rotation="0 180 0", collisionGroup=[1])
    #
    # CollisionMesh(finger,
    #               surfaceMeshFileName="mesh/fingerCollision_part2.stl",
    #               name="CollisionMeshAuto2", translation="-17.5 -12.5 7.5", rotation="0 180 0", collisionGroup=[2])

    finger.addObject('LinearSolverConstraintCorrection')

    ##########################################
    # Visualization                          #
    ##########################################
    fingerVisu = finger.addChild('visu')
    fingerVisu.addObject(
        'MeshSTLLoader', filename=path+"finger.stl", name="loader", translation="-17.5 -12.5 7.5",
        rotation="0 180 0")
    fingerVisu.addObject('OglModel', src="@loader", template='ExtVec3f', color="0.0 0.7 0.7")
    fingerVisu.addObject('BarycentricMapping')

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
    position = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.], [0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]
    longeur = [15, 15, 15, 15, 6, 15]  # beams size
    rateAngularDeformNode = cableNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject('MechanicalObject', template='Vec3d',
                                                          name='rateAngularDeformMO', position=position,
                                                          showIndices="1")
    BeamHookeLawForce = rateAngularDeformNode.addObject('BeamHookeLawForceField', crossSectionShape='circular',
                                                        length=longeur, radius='0.50', youngModulus='5e6')

    ################################
    # Animation (to move the dofs) #
    ################################
    animate = Animation(RigidBaseMO, rateAngularDeformMO)
    rootNode.addObject(animate)

    ##############
    #   Frames   #
    ##############
    frames = [
        "0.0 0 0  0 0 0 1   5 0 0  0 0 0 1  10.0 0 0  0 0 0 1    15.0 0 0  0 0 0 1   20.0 0 0  0 0 0 1  "
        "30.0 0 0  0 0 0 1  35.0 0 0  0 0 0 1   40.0 0 0  0 0 0 1   45.0 0 0  0 0 0 1 55.0 0 0  0 0 0 1 "
        "60.0 0 0  0 0 0 1  66.0 0 0  0 0 0 1   71.0 0 0  0 0 0 1   76.0 0 0  0 0 0 1  81.0 0 0  0 0 0 1"]
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
    # mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=curv_abs_input,
    #                              curv_abs_output=curv_abs_output, input1=inputMO, input2=inputMO_rigid,
    #                              output=outputMO, debug='0', max=2.e-3, deformationAxis=1)
    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=curv_abs_input,
                                 curv_abs_output=curv_abs_output, input1=inputMO, input2=inputMO_rigid,
                                 output=outputMO, debug='0', max=6.e-2, deformationAxis=1, nonColored="0", radius=5)

    # actuators = mappedFrameNode.addObject('actuators')
    # actuator0 = actuators.addObject('SlidingActuator', name="SlidingActuator0", template='Rigid3d',
    #                                    direction='0 0 0 1 0 0', indices=1, maxForce='100000', minForce='-30000')
    cable_position = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0], [15.0, 0.0, 0.0], [20.0, 0.0, 0.0],
                      [30.0, 0.0, 0.0], [35.0, 0.0, 0.0], [40.0, 0.0, 0.0], [45.0, 0.0, 0.0],
                      [55.0, 0.0, 0.0], [60.0, 0.0, 0.0], [66.0, 0.0, 0.0], [71.0, 0.0, 0.0], [76.0, 0.0, 0.0],
                      [81.0, 0.0, 0.0]]
    #  This create a new node in the scene. This node is appended to the finger's node.
    slidingPoint = mappedFrameNode.addChild('slidingPoint')

    # This create a MechanicalObject, a componant holding the degree of freedom of our
    # mechanical modelling. In the case of a cable it is a set of positions specifying
    # the points where the cable is passing by.
    slidingPointMO = slidingPoint.addObject('MechanicalObject', name="cablePos",
                              position=cable_position, showObject="1", showIndices="0")
    slidingPoint.addObject('IdentityMapping')


    mappedPointsNode = slidingPoint.addChild('MappedPoints')
    femPoints.addChild(mappedPointsNode)
    mappedPoints = mappedPointsNode.addObject('MechanicalObject', template='Vec3d', position=FEMpos,
                                                 name="FramesMO", showObject='1', showObjectScale='1')

    inputCableMO = slidingPointMO.getLinkPath()
    inputFEMCableMO = inputFEMCable.getLinkPath()
    outputPointMO = mappedPoints.getLinkPath()

    mappedPointsNode.addObject('QPSlidingConstraint', name="QPConstraint")
    mappedPointsNode.addObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO,
                                  input2=inputCableMO, output=outputPointMO, direction="@../../FramesMO.position")
    ## Get the tree mstate links for the mapping


    return rootNode
