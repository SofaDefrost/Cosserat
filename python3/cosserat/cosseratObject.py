# -*- coding: utf-8 -*-
"""
    Cosserat class in SofaPython3.
"""

__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "October, 26 2021"

from dataclasses import dataclass
import Sofa
from usefulFunctions import buildEdges, pluginList, BuildCosseratGeometry
import numpy as np

cosserat_config = {'init_pos': [0., 0., 0.], 'tot_length': 6, 'nbSectionS': 6,
                   'nbFramesF': 12, 'buildCollisionModel': 1}


# @dataclass
class Cosserat(Sofa.Prefab):
    """ActuatedArm is a reusable sofa model of a S90 servo motor and the tripod actuation arm.
           Parameters:
               -parent:        node where the ServoArm will be attached
                - translation the position in space of the structure
                - eulerRotation the orientation of the structure
                - attachingTo (MechanicalObject)    a rest shape force field will constraint the object
                                                 to follow arm position
           Structure:
           Node : {
                name : 'Cosserat'
                Node0 MechanicalObject :     // Rigid position of the base of the beam
                Node1 MechanicalObject :    // Vec3d, The rate angular composed of the twist and the bending along y and z
                Node1 ForceField          //
                    MechanicalObject     //  The child of the two precedent nodes, Rigid positions
                    Cosserat Mapping  //  it allow the transfer from the local to the global frame
            }
    """
    properties = [
        {'name': 'name', 'type': 'string', 'help': 'Node name', 'default': 'Cosserat'},
        {'name': 'translation', 'type': 'Vec3d', 'help': 'Cosserat base position', 'default': [0., 0., 0.]},
        {'name': 'rotation', 'type': 'Vec3d', 'help': 'Cosserat base Rotation', 'default': [0., 0., 0.]},
        {'name': 'youngModulus', 'type': 'double', 'help': 'Beam Young modulus', 'default': 1.e6},
        {'name': 'poissonRatio', 'type': 'double', 'help': 'Beam poisson ratio', 'default': 0.4},
        {'name': 'shape', 'type': 'string', 'help': 'beam section', 'default': "circular"},
        {'name': 'radius', 'type': 'double', 'help': 'the radius in case of circular section', 'default': 1.0},
        {'name': 'length_Y', 'type': 'double', 'help': 'the radius in case of circular section', 'default': 1.0},
        {'name': 'length_Z', 'type': 'double', 'help': 'the radius in case of circular section', 'default': 1.0},
        {'name': 'attachingToLink', 'type': 'string', 'help': 'a rest shape force field will constraint the object '
                                                              'to follow arm position', 'default': '1'}]

    def __init__(self, *args, **kwargs):
        Sofa.Prefab.__init__(self, *args, **kwargs)
        self.cosserat_geometry = args[0]

        self.solverNode = self.addSolverNode()
        self.rigidBaseNode = self.addRigidBaseNode()
        [positionS, curv_abs_inputS, longeurS, framesF, curv_abs_outputF, frames3D] = \
            BuildCosseratGeometry(self.cosserat_geometry)
        self.cosseratCoordinateNode = self.addCosseratCoordinate(positionS, longeurS)
        self.addCosseratFrame(framesF, curv_abs_inputS, curv_abs_outputF, frames3D)

    def init(self):
        pass

    def addSolverNode(self):
        solverNode = self.addChild('solverNode')
        solverNode.addObject('EulerImplicitSolver', rayleighStiffness="0.2", rayleighMass='0.1')
        solverNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")
        # solverNode.addObject('GenericConstraintCorrection')
        return solverNode

    def addRigidBaseNode(self):
        if self.solverNode is not None:
            rigidBaseNode = self.solverNode.addChild('rigidBase')
        else:
            rigidBaseNode = self.addChild('rigidBase')
        rigidBaseNode.addObject('MechanicalObject', template='Rigid3d', name="RigidBaseMO",
                                translation=self.translation.value, showObjectScale=0.2,
                                rotation=self.rotation.value, showObject=1)
        if int(self.attachingToLink.value):
            rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring',
                                    stiffness=1e8, angularStiffness=1.e8, external_points=0,
                                    mstate="@RigidBaseMO", points=0, template="Rigid3d")
        return rigidBaseNode

    def addCosseratCoordinate(self, positionS, longeurS):
        if self.solverNode is not None:
            cosseratCoordinateNode = self.solverNode.addChild('cosseratCoordinate')
        else:
            cosseratCoordinateNode = self.addChild('cosseratCoordinate')
        cosseratCoordinateNode.addObject('MechanicalObject',
                                         template='Vec3d', name='cosseratCoordinateMO',
                                         position=positionS,
                                         showIndices=0)
        cosseratCoordinateNode.addObject('BeamHookeLawForceField', crossSectionShape=self.shape.value,
                                         length=longeurS, youngModulus=self.youngModulus.value,
                                         poissonRatio=self.poissonRatio.value,
                                         radius=self.radius.value,
                                         lengthY=self.length_Y.value, lengthZ=self.length_Z.value)
        return cosseratCoordinateNode

    def addCosseratFrame(self, framesF, curv_abs_inputS, curv_abs_outputF, frames3D):

        cosseratInSofaFrameNode = self.rigidBaseNode.addChild('cosseratInSofaFrameNode')
        self.cosseratCoordinateNode.addChild(cosseratInSofaFrameNode)
        framesMO = cosseratInSofaFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                                     name="FramesMO", position=framesF,
                                                     showObject=1, showObjectScale=0.1)
        cosseratInSofaFrameNode.addObject('UniformMass', totalMass="0.00022", showAxisSizeFactor='0')
        print(" Rigid Pos :== > ", self.cosseratCoordinateNode.cosseratCoordinateMO.getLinkPath())
        print("  Rate Pos :== > ", self.rigidBaseNode.RigidBaseMO.getLinkPath())
        cosseratInSofaFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=curv_abs_inputS,
                                          curv_abs_output=curv_abs_outputF, name='cosseratMapping',
                                          input1=self.cosseratCoordinateNode.cosseratCoordinateMO.getLinkPath(),
                                          input2=self.rigidBaseNode.RigidBaseMO.getLinkPath(),
                                          output=framesMO.getLinkPath(), debug=0, radius=0)

        self.addObject('MechanicalMatrixMapper', template='Vec3,Rigid3',
                       object1=self.cosseratCoordinateNode.cosseratCoordinateMO.getLinkPath(),
                       object2=self.rigidBaseNode.RigidBaseMO.getLinkPath(),
                       nodeToParse=cosseratInSofaFrameNode.getLinkPath())

        # @todo add the node collision condition
        return cosseratInSofaFrameNode

    def addEdgeCollision(self, position3D, edges):
        collisInstrumentCombined = self.addChild('collisInstrumentCombined')
        collisInstrumentCombined.addObject('EdgeSetTopologyContainer', name="collisEdgeSet", position=position3D,
                                           edges=edges)
        collisInstrumentCombined.addObject('EdgeSetTopologyModifier', name="collisEdgeModifier")
        collisInstrumentCombined.addObject('MechanicalObject', name="CollisionDOFs")
        collisInstrumentCombined.addObject('LineCollisionModel', bothSide="1", group='2')
        collisInstrumentCombined.addObject('PointCollisionModel', bothSide="1", group='2')
        collisInstrumentCombined.addObject('IdentityMapping', name="mapping")
        return collisInstrumentCombined


def createScene(rootNode):
    from stlib3.scene import Scene
    scene = Scene(rootNode, gravity=[0.0, -9810, 0.0], plugins=pluginList + ["SoftRobots.Inverse", "SofaConstraint"],
                  iterative=True)
    rootNode.dt = 0.003
    rootNode.gravity = [0., -9.810, 0.]
    scene.addMainHeader()

    cosserat = scene.Modelling.addChild(Cosserat(cosserat_config, name="cosseroot", radius=0.2))

    return rootNode
