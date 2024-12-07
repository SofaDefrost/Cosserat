from scripts.utils.baseobject import BaseObject
from splib3.numerics import Quat, Vec3, vsub
from math import pi
import numpy as np
from scipy.linalg import logm, inv
from scipy.spatial.transform import Rotation
import Sofa


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
    gX[0:3, 3] = frame[0:3]
    gX[3, 3] = 1

    if curvAbs <= 0.:
        xi = [0., 0., 0.]
    else:
        gXpInv = inv(gXp)
        xiHat = 1 / curvAbs * logm(np.dot(gXpInv, gX))
        xi = [xiHat[2][1], xiHat[0][2], xiHat[1][0], xiHat[0][3], xiHat[1][3], xiHat[2][3]]

    return xi[0:3], gX


class BaseCosserat(BaseObject):
    """
    Base rod, based on Cosserat. With a visual and a collision model.
    """

    deformabletemplate = 'Vec3'

    def __init__(self, modelling, simulation, params, positions, length, name='BaseCosserat', collisionGroup=0):
        super().__init__(modelling, simulation, params, positions, length, name, collisionGroup)
        self.__addRod()
        self.addCylinderTopology()
        self.addVisualModel()

    def __addRod(self):
        self.node = self.modelling.addChild(self.name)
        self.simulation.addChild(self.node)
        self.node.addObject('RequiredPlugin', pluginName=['Cosserat'])

        nbSections = self.params.nbSections
        nbPoints = nbSections + 1

        self.base = self.node.addChild('RigidBase')
        self.base.addObject('MechanicalObject', template='Rigid3', position=self.positions[0])

        rod = self.base.addChild('Rod')
        self.rod = rod
        self.deformable = self.node.addChild('Deformable')
        self.deformable.addChild(rod)

        rod.addObject('EdgeSetTopologyContainer',
                      position=[pos[0:3] for pos in self.positions],
                      edges=[[i, i+1] for i in range(nbSections)])
        rod.addObject('MechanicalObject', template='Rigid3',
                      position=self.positions, showIndices=False, showIndicesScale=0.005, showObject=False,
                      showObjectScale=0.05)
        rod.addObject("BeamInterpolation")

        # Convert Rigid3 orientation description to Cosserat bending description
        # [[torsion strain, y_bending strain, z_bending strain]]
        gX = np.zeros((4, 4), dtype=float)
        gX[0:3, 0:3] = Rotation.from_quat(Quat(self.positions[0][3:7])).as_matrix()
        gX[0:3, 3] = self.positions[0][0:3]
        gX[3, 3] = 1

        lengths = []
        strain = []
        for i in range(len(self.positions) - 1):
            length = Vec3(vsub(self.positions[i][0:3], self.positions[i + 1][0:3])).getNorm()
            lengths.append(length)
            xi, gX = getStrainFromQuat(self.positions[i + 1], length, gX)
            strain.append(xi)

        self.deformable.addObject('MechanicalObject', position=strain, rest_position=[0, 0, 0] * nbSections)
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
                      input1=self.deformable.MechanicalObject.getLinkPath(),
                      input2=self.base.MechanicalObject.getLinkPath(),
                      output=rod.getLinkPath(),
                      debug=False, baseIndex=0)

        return rod


# Test scene
def createScene(rootnode):
    from scripts.utils.header import addHeader, addSolvers
    import params

    settings, modelling, simulation = addHeader(rootnode)
    rootnode.VisualStyle.displayFlags = ['hideBehavior']
    addSolvers(simulation, iterative=False)

    nbSections = params.CableParameters.nbSections
    length = 5
    dx = length / nbSections
    positions = [[dx * i, 0, 0, 0, 0, 0, 1] for i in range(nbSections + 1)]
    beam = BaseCosserat(modelling, simulation, params.CableParameters, positions, length)
    beam.node.RigidBase.addObject('FixedConstraint', indices=0)
