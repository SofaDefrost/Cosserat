
from splib3.numerics import Quat, Vec3, vsub
from math import pi
import numpy as np
from scipy.linalg import logm, inv
from scipy.spatial.transform import Rotation
import Sofa
import Sofa.Core
from dataclasses import dataclass, field
from math import pi


# Units are m, kg, s


@dataclass
class CableParameters:
    youngModulus: float = 1.205e11
    poissonRatio: float = 0.499
    radius: float = 0.004
    density: float = 7.850e3
    nbSections: int = 5 
    length: float = 10


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

def getStrainFromQuat(frame, curvAbs, gXp):
    """

    Args:
        frame:
        curvAbs: abscissa curve
        gXp: previous matrix

    Returns:

    """

    gX = np.zeros((4, 4), dtype=float)
    gX[0:3, 0:3] = Rotation.from_quat(frame[3:7]).as_matrix()
    gX[:3, 3] = frame[:3]
    gX[3, 3] = 1

    if curvAbs <= 0.:
        xi = [0., 0., 0.]
    else:
        gXpInv = inv(gXp)
        xiHat = 1 / curvAbs * logm(np.dot(gXpInv, gX))
        xi = [xiHat[2][1], xiHat[0][2], xiHat[1][0], xiHat[0][3], xiHat[1][3], xiHat[2][3]]

    return xi[:3], gX


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

        # Convert Rigid3 orientation description to Cosserat bending description
        # [[torsion strain, y_bending strain, z_bending strain]]
        gX = np.zeros((4, 4), dtype=float)
        gX[:3, :3] = Rotation.from_quat(Quat(self.positions[0][3:7])).as_matrix()
        gX[:3, 3] = self.positions[0][:3]
        gX[3, 3] = 1

        lengths = []
        strain = []
        for i in range(len(self.positions) - 1):
            length = Vec3(vsub(self.positions[i][:3], self.positions[i + 1][:3])).getNorm()
            lengths.append(length)
            xi, gX = getStrainFromQuat(self.positions[i + 1], length, gX)
            print(f'the strain xi : {xi}')
            strain.append([float(xi[0]),float(xi[1]),float(xi[2])])

        self.deformable.addObject('MechanicalObject', position=strain, name="strain_state", rest_position=[0, 0, 0] * nbSections)
        self.deformable.addObject('BeamHookeLawForceField',
                                  youngModulus=self.params.youngModulus,
                                  poissonRatio=self.params.poissonRatio,
                                  radius=self.params.radius,
                                  length=lengths)
        self.node.addData(name="indexExtremity", type='int', value=nbPoints - 1)

        totalMass = self.params.density * self.length * self.params.radius * self.params.radius * pi
        rod.addObject('UniformMass', totalMass=totalMass, showAxisSizeFactor=0.01)
        l = 0
        curv_abs = [l]
        for length in lengths:
            l += length
            curv_abs.append(l)
        rod.addObject('DiscreteCosseratMapping',
                      curv_abs_input=curv_abs,
                      curv_abs_output=curv_abs,
                      input1=self.deformable.strain_state.getLinkPath(),
                      input2=self.base.rigid_base.getLinkPath(),
                      output=rod.getLinkPath(),
                      debug=False, baseIndex=0)

        return rod


# Test scene
def createScene(rootnode):
    from header import addHeader, addSolvers, modelling, simulation = addHeader(rootnode)
    rootnode.VisualStyle.displayFlags = ['hideBehavior']
    addSolvers(simulation, iterative=False)

    nbSections = CableParameters.nbSections
    length = CableParameters.length
    dx = length / nbSections
    ##positions = [[float(dx * i), float(0.), float(0.), float(0.), float(0.), float(0.), 1] for i in range(nbSections + 1)]
    
    # input1 to test the two torsion
    #positions = [[ -0.00,-0.00,-0.00,0.00,0.00,0.00,1.00 ], [ 1.97,0.20,-0.20,0.00,0.10,0.10,0.99 ], [ 3.79,0.78,-0.78,-0.00,0.20,0.20,0.96 ], 
    #             [ 5.31,1.69,-1.69,-0.00,0.29,0.29,0.91 ], [ 6.40,2.87,-2.87,-0.00,0.38,0.38,0.84 ], [ 6.98,4.22,-4.22,-0.00,0.46,0.46,0.76 ]]

    # input2 for to test the tosion 
    positions=[[-0.00,-0.00,-0.00,0.00,0.00,0.00,1.00 ], [ 2.00,-0.00,-0.00,0.10,0.00,0.00,1.00 ], [ 4.00,-0.00,-0.00,0.20,0.00,0.00,0.98 ], 
               [ 6.00,-0.00,-0.00,0.30,0.00,0.00,0.96 ], [ 8.00,-0.00,-0.00,0.39,0.00,0.00,0.92 ], [ 10.00,-0.00,-0.00,0.48,0.00,0.00,0.88 ]]
    beam = BaseCosserat(modelling, simulation, CableParameters, positions, length)

    beam.node.RigidBase.addObject('FixedProjectiveConstraint', indices=0)
