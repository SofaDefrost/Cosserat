# -*- coding: utf-8 -*-
"""
    Cosserat class in SofaPython3.
"""

__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "October, 26 2021"

# from dataclasses import dataclass
import Sofa
from cosserat.usefulFunctions import buildEdges, pluginList, BuildCosseratGeometry
from splib3.numerics import Quat
from cosserat.utils import addEdgeCollision, addPointsCollision

cosserat_config = {'init_pos': [0., 0., 0.], 'tot_length': 6, 'nbSectionS': 6,
                   'nbFramesF': 12, 'buildCollisionModel': 1, 'beamMass': 0.22}


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
    prefabParameters = [
        {'name': 'name', 'type': 'string',
            'help': 'Node name', 'default': 'Cosserat'},
        {'name': 'position', 'type': 'Rigid3d::VecCoord', 'help': 'Cosserat base position',
         'default': [[0., 0., 0., 0, 0, 0, 1.]]},
        {'name': 'translation', 'type': 'Vec3d',
            'help': 'Cosserat base Rotation', 'default': [0., 0., 0.]},
        {'name': 'rotation', 'type': 'Vec3d',
            'help': 'Cosserat base Rotation', 'default': [0., 0., 0.]},
        {'name': 'youngModulus', 'type': 'double',
            'help': 'Beam Young modulus', 'default': 1.e6},
        {'name': 'poissonRatio', 'type': 'double',
            'help': 'Beam poisson ratio', 'default': 0.4},
        {'name': 'shape', 'type': 'string',
            'help': 'beam section', 'default': "circular"},
        {'name': 'radius', 'type': 'double',
            'help': 'the radius in case of circular section', 'default': 0.02},
        {'name': 'length_Y', 'type': 'double',
            'help': 'the radius in case of circular section', 'default': 1.0},
        {'name': 'length_Z', 'type': 'double',
            'help': 'the radius in case of circular section', 'default': 1.0},
        {'name': 'rayleighStiffness', 'type': 'double', 'help': 'Rayleigh damping - stiffness matrix coefficient',
         'default': 0.0},
        {'name': 'attachingToLink', 'type': 'string', 'help': 'a rest shape force field will constraint the object '
                                                              'to follow arm position', 'default': '1'},
        {'name': 'showObject', 'type': 'string', 'help': ' Draw object arrow ', 'default': '0'}]

    def __init__(self, *args, **kwargs):
        Sofa.Prefab.__init__(self, *args, **kwargs)
        self.cosseratGeometry = kwargs['cosseratGeometry']
        self.beamMass = self.cosseratGeometry['beamMass']
        self.parent = kwargs.get('parent')
        self.useInertiaParams = False
        self.radius = kwargs.get('radius', )

        if self.parent.hasObject("EulerImplicitSolver") is False:
            print('The code does not have parent EulerImplicite')
            self.solverNode = self.addSolverNode()
        else:
            self.solverNode = self.parent

        if 'inertialParams' in kwargs:
            self.useInertiaParams = True
            self.inertialParams = kwargs['inertialParams']

        self.rigidBaseNode = self.addRigidBaseNode()
        [positionS, curv_abs_inputS, longeurS, framesF, curv_abs_outputF, self.frames3D] = \
            BuildCosseratGeometry(self.cosseratGeometry)

        self.cosseratCoordinateNode = self.addCosseratCoordinate(
            positionS, longeurS)
        self.cosseratFrame = self.addCosseratFrame(
            framesF, curv_abs_inputS, curv_abs_outputF)
        # print(f'=== > {curv_abs_inputS}')

    def init(self):
        pass

    def addCollisionModel(self):
        tab_edges = buildEdges(self.frames3D)
        return addEdgeCollision(self.cosseratFrame, self.frames3D, tab_edges)

    def addPointCollisionModel(self, nodeName='CollisionPoints'):
        tab_edges = buildEdges(self.frames3D)
        return addPointsCollision(self.cosseratFrame, self.frames3D, tab_edges, nodeName)

    def addSlidingPoints(self):
        slidingPoint = self.cosseratFrame.addChild('slidingPoint')
        slidingPoint.addObject('MechanicalObject', name="slidingPointMO", position=self.frames3D,
                               showObject="0", showIndices="0")
        slidingPoint.addObject('IdentityMapping')
        return slidingPoint

    def addSlidingPointsWithContainer(self):
        slidingPoint = self.cosseratFrame.addChild('slidingPoint')
        container = slidingPoint.addObject("PointSetTopologyContainer")
        modifier = slidingPoint.addObject("PointSetTopologyModifier")
        slidingPoint.addObject('MechanicalObject', name="slidingPointMO", position=self.frames3D,
                               showObject="1", showIndices="0")
        slidingPoint.addObject('IdentityMapping')
        return slidingPoint

    def addSolverNode(self):
        solverNode = self.addChild('solverNode')
        solverNode.addObject('EulerImplicitSolver',
                             rayleighStiffness="0.2", rayleighMass='0.1')
        solverNode.addObject('SparseLDLSolver', name='solver',
                             template="CompressedRowSparseMatrixd")
        solverNode.addObject('GenericConstraintCorrection')
        return solverNode

    def addRigidBaseNode(self):
        rigidBaseNode = self.addChild('rigidBase')

        # trans = [t for t in self.translation.value]
        trans = list(self.translation.value)
        rot = list(self.rotation.value)
        # @todo converter
        positions = [list(pos) for pos in self.position.value]

        rigidBaseNode.addObject('MechanicalObject', template='Rigid3d', name="RigidBaseMO", showObjectScale=0.2,
                                translation=trans, position=positions, rotation=rot,
                                showObject=int(self.showObject.value))

        # one can choose to set this to false and directly attach the beam base
        # to a control object in order to be able to drive it.
        if int(self.attachingToLink.value):
            print("Adding the rest shape to the base")
            rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness=1e8, angularStiffness=1.e8,
                                    external_points=0, mstate="@RigidBaseMO", points=0, template="Rigid3d")
        return rigidBaseNode

    def addCosseratCoordinate(self, bendingStates, listOfSectionsLength):
        cosseratCoordinateNode = self.addChild('cosseratCoordinate')
        cosseratCoordinateNode.addObject('MechanicalObject',
                                         template='Vec3d', name='cosseratCoordinateMO',
                                         position=bendingStates,
                                         showIndices=0)

        if self.useInertiaParams is False:
            cosseratCoordinateNode.addObject('BeamHookeLawForceField', crossSectionShape=self.shape.value,
                                             length=listOfSectionsLength, radius=self.radius.value,
                                             youngModulus=self.youngModulus.value, poissonRatio=self.poissonRatio.value,
                                             rayleighStiffness=self.rayleighStiffness.value,
                                             lengthY=self.length_Y.value, lengthZ=self.length_Z.value)
        else:
            GA = self.inertialParams['GA']
            GI = self.inertialParams['GI']
            EA = self.inertialParams['EA']
            EI = self.inertialParams['EI']
            print(f'{GA}')
            cosseratCoordinateNode.addObject('BeamHookeLawForceField', crossSectionShape=self.shape.value,
                                             length=listOfSectionsLength, radius=self.radius.value, useInertiaParams=True,
                                             GI=GI, GA=GA, EI=EI, EA=EA, rayleighStiffness=self.rayleighStiffness.value,
                                             lengthY=self.length_Y.value, lengthZ=self.length_Z.value)
        return cosseratCoordinateNode

    def addCosseratFrame(self, framesF, curv_abs_inputS, curv_abs_outputF):

        cosseratInSofaFrameNode = self.rigidBaseNode.addChild(
            'cosseratInSofaFrameNode')
        self.cosseratCoordinateNode.addChild(cosseratInSofaFrameNode)
        framesMO = cosseratInSofaFrameNode.addObject(
            'MechanicalObject', template='Rigid3d', name="FramesMO", position=framesF, showObject=0, showObjectScale=0.1)
        if self.beamMass != 0.:
            cosseratInSofaFrameNode.addObject(
                'UniformMass', totalMass=self.beamMass, showAxisSizeFactor='0')

        cosseratInSofaFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=curv_abs_inputS,
                                          curv_abs_output=curv_abs_outputF, name='cosseratMapping',
                                          input1=self.cosseratCoordinateNode.cosseratCoordinateMO.getLinkPath(),
                                          input2=self.rigidBaseNode.RigidBaseMO.getLinkPath(),
                                          output=framesMO.getLinkPath(), debug=0, radius=self.radius.value)
        return cosseratInSofaFrameNode


def createScene(rootNode):
    rootNode.addObject('RequiredPlugin', name='plugins', pluginName=[pluginList,
                                                                     ['SofaEngine', 'SofaLoader', 'SofaSimpleFem',
                                                                      'SofaExporter']])
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideCollisionModels '
                                                   'hideBoundingCollisionModels hideForceFields '
                                                   'hideInteractionForceFields hideWireframe showMechanicalMappings')
    rootNode.findData('dt').value = 0.01
    rootNode.findData('gravity').value = [0., -9.81, 0.]
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('GenericConstraintSolver',
                       tolerance=1e-5, maxIterations=5e2)
    rootNode.addObject('Camera', position="-35 0 280", lookAt="0 0 0")

    solverNode = rootNode.addChild('solverNode')
    solverNode.addObject('EulerImplicitSolver',
                         rayleighStiffness="0.2", rayleighMass='0.1')
    solverNode.addObject('SparseLDLSolver', name='solver',
                         template="CompressedRowSparseMatrixd")
    solverNode.addObject('GenericConstraintCorrection')

    cosserat = solverNode.addChild(
        Cosserat(parent=solverNode, cosseratGeometry=cosserat_config, name="cosserat", radius=0.15))

    # use this to add the collision if the beam will interact with another object
    collisionModel = cosserat.addCollisionModel()

    # Attach a force at the beam tip,
    # we can attach this force to non mechanical node thanks to the MechanicalMatrixMapper component
    beamFrame = cosserat.cosseratFrame
    beamFrame.addObject('ConstantForceField', name='constForce', showArrowSize=1.e-2, indices=12,
                       forces=[0., -100., 0., 0., 0., 0.])

    return rootNode
