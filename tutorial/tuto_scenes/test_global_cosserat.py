
from re import X
from splib3.numerics import Quat, Vec3, vsub
from math import pi
import numpy as np
from scipy.linalg import logm, inv
from scipy.spatial.transform import Rotation
import Sofa
import Sofa.Core
from dataclasses import dataclass, field
from math import pi


class BaseObject:

    node: Sofa.Core.Node
    base: Sofa.Core.Node
    deformable: Sofa.Core.Node
    rod: Sofa.Core.Node

    def __init__(self, modelling, simulation, params, positions, length, name='BaseBeam', collisionGroup=0):
        self.modelling = modelling
        self.simulation = simulation
        self.params = params
        self.positions = positions
        self.length = length
        self.name = name
        self.collisionGroup = collisionGroup

    def addCylinderTopology(self):
        positions = self.rod.MechanicalObject.position.value

        topology = self.modelling.Topology.addChild(f'{self.name}CylinderTopo')
        edgetopo = topology.addChild('Edge')
        edgetopo.addObject('EdgeSetTopologyContainer', edges=[[k, k + 1] for k in range(len(positions) - 1)])
        edgetopo.addObject('EdgeSetTopologyModifier')
        edgetopo.addObject('MechanicalObject', template='Rigid3',
                           position=self.rod.MechanicalObject.position.getLinkPath())

        quadtopo = edgetopo.addChild('Quad')
        quadtopo.addObject('QuadSetTopologyContainer')
        quadtopo.addObject('QuadSetTopologyModifier')
        quadtopo.addObject('MechanicalObject')
        quadtopo.addObject('Edge2QuadTopologicalMapping',
                           input=edgetopo.EdgeSetTopologyContainer.getLinkPath(),
                           output=quadtopo.QuadSetTopologyContainer.getLinkPath(),
                           flipNormals=True, nbPointsOnEachCircle=10, radius=self.params.radius)

    def addVisualModel(self):
        quadtopo = self.modelling.Topology.getChild(f'{self.name}CylinderTopo').Edge.Quad

        visual = self.rod.addChild('VisualModel')
        visual.addObject('MeshTopology', position=quadtopo.MechanicalObject.position.getLinkPath(),
                         quads=quadtopo.QuadSetTopologyContainer.quads.getLinkPath())
        visual.addObject('OglModel', color=[0.1, 0.1, 0.1, 1.0])
        visual.addObject('SkinningMapping')

    def addCollisionModel(self):
        quadtopo = self.modelling.Topology.getChild(f'{self.name}CylinderTopo').Edge.Quad

        collision = self.rod.addChild('CollisionModel')
        collision.addObject('MeshTopology', position=quadtopo.MechanicalObject.position.getLinkPath(),
                            quads=quadtopo.QuadSetTopologyContainer.quads.getLinkPath())
        collision.addObject('MechanicalObject')
        collision.addObject('TriangleCollisionModel', group=self.collisionGroup)
        collision.addObject('PointCollisionModel', group=self.collisionGroup)
        collision.addObject('LineCollisionModel', group=self.collisionGroup)
        collision.addObject('SkinningMapping')
    
    def addLocal : 

class BaseCosserat(BaseObject):
    """
    Base rod, based on Cosserat. With a visual and a collision model.
    """

    deformabletemplate = 'Vec3'

    def __init__(self, modelling, simulation, params, positions, length, name='BaseCosserat', collisionGroup=0):
        super().__init__(modelling, simulation, params, positions, length, name, collisionGroup)
        self.__addRod()
        #self.addCylinderTopology()
        #self.addVisualModel()

    def __addRod(self):
        self.node = self.modelling.addChild(self.name)
        self.simulation.addChild(self.node)
        self.node.addObject('RequiredPlugin', pluginName=['Cosserat'])

        nbSections = CableParameters.nbSections
        nbPoints = nbSections + 1

        self.base = self.node.addChild('RigidBase')
        self.base.addObject('MechanicalObject', template='Rigid3', name="rigid_base", position=self.positions[0])

        rod = self.base.addChild('Rod')
        self.rod = rod
        self.deformable = self.node.addChild('Deformable')
        self.deformable.addChild(rod)

        rod.addObject('EdgeSetTopologyContainer',
                      position=[pos[:3] for pos in self.positions],
                      edges=[[i, i+1] for i in range(nbSections)])
        rod.addObject('MechanicalObject', template='Rigid3', name="frame_mstate",
                      position=self.positions, showIndices=False, showIndicesScale=0.005, showObject=False,
                      showObjectScale=0.05)
        rod.addObject("BeamInterpolation")

