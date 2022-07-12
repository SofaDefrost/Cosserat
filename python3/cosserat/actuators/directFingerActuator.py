# -*- coding: utf-8 -*-

"""_summary_ Basic scene using Cosserat in SofaPython3.
    The Finger is modeled usind FEM modèle while de cable is modeled using cosserat theory.
    The link between these two meachanical models is constraint based using Lagrangian Multiplier

Returns:
    _type_: _description_
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
from cosserat.usefulFunctions import buildEdges, pluginList, BuildCosseratGeometry
from cosserat.cosseratObject import Cosserat

# from stlib3.scene import Node
path = f'{os.path.dirname(os.path.abspath(__file__))}/../mesh/'


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

    def _extracted_from_onKeypressedEvent_10(self, qOld, posA, angularRate):
        qNew = Quat.createFromEuler([0., angularRate, 0.], 'ryxz')
        qNew.normalize()
        qNew.rotateFromQuat(qOld)
        for i in range(4):
            posA[0][i+3] = qNew[i]

    def onKeypressedEvent(self, event):
        key = event['key']
        if ord(key) == 19:  # up
            with self.rigidBaseMO.rest_position.writeable() as posA:
                qOld = Quat()
                for i in range(4):
                    qOld[i] = posA[0][i+3]

                self._extracted_from_onKeypressedEvent_10(
                    self, qOld, posA, self.angularRate)

        if ord(key) == 21:  # down
            with self.rigidBaseMO.rest_position.writeable() as posA:
                qOld = Quat()
                for i in range(4):
                    qOld[i] = posA[0][i+3]

                self._extracted_from_onKeypressedEvent_10(
                    qOld, posA, -self.angularRate)

        if ord(key) == 18:  # left
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][0] -= self.rate
        if ord(key) == 20:  # right
            with self.rigidBaseMO.rest_position.writeable() as posA:
                posA[0][0] += self.rate


def createScene(rootNode):

    from stlib3.scene import MainHeader

    MainHeader(rootNode, plugins=pluginList,
               repositoryPaths=[os.getcwd()])
    rootNode.addObject(
        'VisualStyle', displayFlags='showVisualModels showInteractionForceFields showWireframe showMechanicalMappings')

    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('GenericConstraintSolver',
                       tolerance="1e-20", maxIterations="500", printLog="0")

    gravity = [0, 0, 0]
    rootNode.gravity.value = gravity
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.addObject('OglSceneFrame', style="Arrows", alignment="TopRight")

    ##########################################
    # FEM Model                              #
    ##########################################
    finger = rootNode.addChild('finger')
    finger.addObject('EulerImplicitSolver', name='odesolver',
                     firstOrder='0', rayleighMass=0.1, rayleighStiffness=0.1)
    finger.addObject('SparseLDLSolver', name='preconditioner')

    # Add a componant to load a VTK tetrahedral mesh and expose the resulting topology in the scene .
    finger.addObject('MeshVTKLoader', name='loader', filename=path +
                     'finger.vtk', translation="-17.5 -12.5 7.5", rotation="0 180 0")
    finger.addObject('TetrahedronSetTopologyContainer',
                     src='@loader', name='container')
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

    finger.addObject('BoxROI', name='ROI1',
                     box='-18 -15 -8 2 -3 8', drawBoxes='true')
    finger.addObject('RestShapeSpringsForceField',
                     points='@ROI1.indices', stiffness='1e12')

    ##########################################
    # Cable points                           #
    ##########################################
    # Mappe points inside the meca, this points will be use for the bilateral mapping
    FEMpos = [" 0.0 0 0 15 0 0 30 0 0 45 0 0 60 0 0 66 0 0 81 0.0 0.0"]

    femPoints = finger.addChild('femPoints')
    inputFEMCable = femPoints.addObject('MechanicalObject', name="pointsInFEM", position=FEMpos, showObject="1",
                                        showIndices="1")
    femPoints.addObject('BarycentricMapping')

    finger.addObject('LinearSolverConstraintCorrection')

    ##########################################
    # Visualization                          #
    ##########################################
    fingerVisu = finger.addChild('visu')
    fingerVisu.addObject('MeshSTLLoader', filename=f'{path}finger.stl',
                         name="loader", translation="-17.5 -12.5 7.5", rotation="0 180 0")
    fingerVisu.addObject('OglModel', src="@loader",
                         template='ExtVec3f', color="0.0 0.7 0.7")
    fingerVisu.addObject('BarycentricMapping')

    # ###############
    # New adds to use the sliding Actuator
    ###############
    cableNode = rootNode.addChild('cableNode')
    cableNode.addObject('EulerImplicitSolver', firstOrder="0",
                        rayleighStiffness="0.1", rayleighMass='0.1')
    cableNode.addObject('SparseLUSolver', name='solver')
    cableNode.addObject('GenericConstraintCorrection')

    beamGeometrie = {'init_pos': [0., 0., 0.], 'tot_length': 81, 'nbSectionS': 12,
                     'nbFramesF': 30, 'buildCollisionModel': 0, 'beamMass': 0.}
    cosserat = cableNode.addChild(Cosserat(parent=cableNode, cosseratGeometry=beamGeometrie, radius=0.5,
                                           useCollisionModel=True, name="cosserat", youngModulus=5e6, poissonRatio=0.4))

    cableNode.addObject(Animation(cosserat.rigidBaseNode.RigidBaseMO,
                                  cosserat.cosseratCoordinateNode.cosseratCoordinateMO))
    mappedFrameNode = cosserat.cosseratFrame

    #  This create a new node in the scene. This node is appended to the finger's node.
    slidingPoint = mappedFrameNode.addChild('slidingPoint')

    # This create a MechanicalObject, a componant holding the degree of freedom of our
    # mechanical modelling. In the case of a cable it is a set of positions specifying
    # the points where the cable is passing by.
    slidingPointMO = slidingPoint.addObject('MechanicalObject', name="cablePos",
                                            position=cosserat.frames3D, showObject="1", showIndices="0")
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
    # Get the tree mstate links for the mapping

    return rootNode
