# -*- coding: utf-8 -*-

"""_summary_ Basic scene using Cosserat in SofaPython3.
    The Finger is modeled using FEM model while de cable is modeled using cosserat theory.
    The link between these two mechanical models is performed by constraint based on Lagrangian Multiplier

Returns:
    _type_: _description_
"""

__authors__ = "Younes"
__contact__ = "adagolodjo@protonmail.com, yinoussa.adagolodjo@inria.fr"
__version__ = "1.0.0"
__copyright__ = "(c) 2020,Inria"
__date__ = "july 2023"

import Sofa
import os
from cosserat.cosseratObject import Cosserat
from useful.header import addHeader, addVisual, addSolverNode, Finger
from controler import FingerController

path = f'{os.path.dirname(os.path.abspath(__file__))}/mesh/'

femPos = [" 0.0 0 0 15 0 0 30 0 0 45 0 0 60 0 0 66 0 0 81 0.0 0.0"]


def createScene(rootNode):
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    addVisual(rootNode)

    addHeader(rootNode, isConstrained=False)
    rootNode.addObject('FreeMotionAnimationLoop', parallelCollisionDetectionAndFreeMotion=False,
                       parallelODESolving=False)
    rootNode.addObject('GenericConstraintSolver', name='ConstraintSolver', tolerance=1e-20, maxIterations=1000,
                       multithreading=False)
    rootNode.findData('gravity').value = [0., 0., 0.]

    femFingerNode = rootNode.addChild('femFingerNode')

    """ Add FEM finger to the scene"""
    fingerNode, femPointsNode = Finger(femFingerNode, name="Finger", rotation=[0.0, 180.0, 0.0],
                                       translation=[-17.5, -12.5, 7.5], path=path)
    """ 
        Add Cosserat cable to the scene
    """
    beamGeometrie = {'init_pos': [0., 0., 0.], 'tot_length': 81, 'nbSectionS': 12,
                     'nbFramesF': 30, 'buildCollisionModel': 0, 'beamMass': 0.}

    cableNode = addSolverNode(rootNode, name="cableSolver", isConstrained=True)
    cosseratCable = cableNode.addChild(
        Cosserat(parent=cableNode, cosseratGeometry=beamGeometrie, name="cosserat", radius=0.5,
                 youngModulus=5e6, poissonRatio=0.4))

    mappedFrameNode = cosseratCable.cosseratFrame
    #  This creates a new node in the scene. This node is appended to the finger's node.
    cable3DPosNode = mappedFrameNode.addChild('cable3DPosNode')

    # This creates a MechanicalObject, a componant holding the degree of freedom of our
    # mechanical modelling. In the case of a cable it is a set of positions specifying
    # the points where the cable is passing by.
    cable3DPosNodeMO = cable3DPosNode.addObject('MechanicalObject', name="cablePos", position=cosseratCable.frames3D,
                                                showObject="1", showIndices="0")
    cable3DPosNode.addObject('IdentityMapping')

    """ These positions are in fact the distance between fem points and the cable points"""
    distancePointsNode = cable3DPosNode.addChild('distancePoints')
    femPointsNode.addChild(distancePointsNode)
    mappedPoints = distancePointsNode.addObject('MechanicalObject', template='Vec3d', position=femPos,
                                                name="distancePointsMO", showObject='1', showObjectScale='1')
    # cableBaseMo = cosseratCable.rigidBaseNode.getObject('MechanicalObject')
    # cosseratCoordinateMo = cosseratCable.cosseratCoordinateNode.getObject('MechanicalObject')

    """The controller of the cable is added to the scene"""
    cableNode.addObject(FingerController(cosseratCable.rigidBaseNode.RigidBaseMO,
                                         cosseratCable.cosseratCoordinateNode.cosseratCoordinateMO))

    inputCableMO = cable3DPosNodeMO.getLinkPath()
    inputFEMCableMO = femPointsNode.getLinkPath()
    outputPointMO = mappedPoints.getLinkPath()
    """ This constraint is used to compute the distance between the cable and the fem points"""
    distancePointsNode.addObject('QPSlidingConstraint', name="QPConstraint")
    distancePointsNode.addObject('DifferenceMultiMapping', name="pointsMulti", input1=inputFEMCableMO, indices="5",
                                 input2=inputCableMO, output=outputPointMO, direction="@../../FramesMO.position")

    return

    # ###############
    # New adds to use the sliding Actuator
    ###############
    cableNode = rootNode.addChild('cableNode')
    cableNode.addObject('EulerImplicitSolver', firstOrder="0",
                        rayleighStiffness="0.1", rayleighMass='0.1')
    cableNode.addObject('SparseLUSolver', name='solver')
    cableNode.addObject('GenericConstraintCorrection')

    cosserat = cableNode.addChild(Cosserat(parent=cableNode, cosseratGeometry=beamGeometrie, radius=0.5,
                                           useCollisionModel=True, name="cosserat", youngModulus=5e6, poissonRatio=0.4))

    cableNode.addObject(Animation(cosserat.rigidBaseNode.RigidBaseMO,
                                  cosserat.cosseratCoordinateNode.cosseratCoordinateMO))

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
