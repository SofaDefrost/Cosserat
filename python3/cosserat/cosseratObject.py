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
from cosserat.usefulFunctions import buildEdges, pluginList, BuildCosseratGeometry

cosserat_config = {'init_pos': [0., 0., 0.], 'tot_length': 6, 'nbSectionS': 6,
                   'nbFramesF': 12, 'buildCollisionModel': 1, 'beamMass': 0.22}


# @dataclass
def addEdgeCollision(parentNode, position3D, edges):
    collisInstrumentCombined = parentNode.addChild('collisInstrumentCombined')
    collisInstrumentCombined.addObject('EdgeSetTopologyContainer', name="collisEdgeSet", position=position3D,
                                       edges=edges)
    collisInstrumentCombined.addObject('EdgeSetTopologyModifier', name="collisEdgeModifier")
    collisInstrumentCombined.addObject('MechanicalObject', name="CollisionDOFs")
    collisInstrumentCombined.addObject('LineCollisionModel', bothSide="1", group='2')
    collisInstrumentCombined.addObject('PointCollisionModel', bothSide="1", group='2')
    collisInstrumentCombined.addObject('IdentityMapping', name="mapping")
    return collisInstrumentCombined


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
        {'name': 'position', 'type': 'Rigid3d::VecCoord', 'help': 'Cosserat base position',
         'default': [[0., 0., 0., 0, 0, 0, 1.]]},
        {'name': 'translation', 'type': 'Vec3d', 'help': 'Cosserat base Rotation', 'default': [0., 0., 0.]},
        {'name': 'rotation', 'type': 'Vec3d', 'help': 'Cosserat base Rotation', 'default': [0., 0., 0.]},
        {'name': 'youngModulus', 'type': 'double', 'help': 'Beam Young modulus', 'default': 1.e6},
        {'name': 'poissonRatio', 'type': 'double', 'help': 'Beam poisson ratio', 'default': 0.4},
        {'name': 'shape', 'type': 'string', 'help': 'beam section', 'default': "circular"},
        {'name': 'radius', 'type': 'double', 'help': 'the radius in case of circular section', 'default': 1.0},
        {'name': 'length_Y', 'type': 'double', 'help': 'the radius in case of circular section', 'default': 1.0},
        {'name': 'length_Z', 'type': 'double', 'help': 'the radius in case of circular section', 'default': 1.0},
        {'name': 'attachingToLink', 'type': 'string', 'help': 'a rest shape force field will constraint the object '
                                                              'to follow arm position', 'default': '1'},
        {'name': 'showObject', 'type': 'string', 'help': ' Draw object arrow ', 'default': '0'}]


    def __init__(self, *args, **kwargs):
        Sofa.Prefab.__init__(self, *args, **kwargs)
        self.cosseratGeometry = kwargs['cosseratGeometry']
        self.beamMass = self.cosseratGeometry['beamMass']
        self.parent = kwargs.get('parent', None)

        if self.parent.hasObject("EulerImplicitSolver") is False:
            self.solverNode = self.addSolverNode()
        else:
            self.solverNode = self.parent

        self.rigidBaseNode = self.addRigidBaseNode()
        [positionS, curv_abs_inputS, longeurS, framesF, curv_abs_outputF, self.frames3D] = \
            BuildCosseratGeometry(self.cosseratGeometry)
        self.cosseratCoordinateNode = self.addCosseratCoordinate(positionS, longeurS)
        self.cosseratFrame = self.addCosseratFrame(framesF, curv_abs_inputS, curv_abs_outputF)
        # print(f'=== > {curv_abs_inputS}')

    def init(self):
        pass

    def addCollisionModel(self):
        tab_edges = buildEdges(self.frames3D)
        return addEdgeCollision(self.cosseratFrame, self.frames3D, tab_edges)

    def addSolverNode(self):
        solverNode = self.addChild('solverNode')
        solverNode.addObject('EulerImplicitSolver', rayleighStiffness="0.2", rayleighMass='0.1')
        solverNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")
        solverNode.addObject('GenericConstraintCorrection')
        return solverNode

    def addRigidBaseNode(self):
        rigidBaseNode = self.solverNode.addChild('rigidBase')

        trans = [t for t in self.translation.value]
        rot = [r for r in self.rotation.value]
        positions = []
        for pos in self.position.value:
            _pos = [p for p in pos]
            positions.append(_pos)
        rigidBaseNode.addObject('MechanicalObject', template='Rigid3d', name="RigidBaseMO",
                                showObjectScale=0.2, translation=trans,
                                position=positions, rotation=rot, showObject=int(self.showObject.value))
        # one can choose to set this to false and directly attach the beam base
        # to a control object in order to be able to drive it.
        if int(self.attachingToLink.value):
            rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring',
                                    stiffness=1e8, angularStiffness=1.e8, external_points=0,
                                    mstate="@RigidBaseMO", points=0, template="Rigid3d")
        return rigidBaseNode

    def addCosseratCoordinate(self, positionS, longeurS):
        cosseratCoordinateNode = self.solverNode.addChild('cosseratCoordinate')
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

    def addCosseratFrame(self, framesF, curv_abs_inputS, curv_abs_outputF):

        cosseratInSofaFrameNode = self.rigidBaseNode.addChild('cosseratInSofaFrameNode')
        self.cosseratCoordinateNode.addChild(cosseratInSofaFrameNode)
        framesMO = cosseratInSofaFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                                     name="FramesMO", position=framesF,
                                                     showObject=int(self.showObject.value), showObjectScale=0.1)
        if self.beamMass != 0.:
            cosseratInSofaFrameNode.addObject('UniformMass', totalMass=self.beamMass, showAxisSizeFactor='0')
        cosseratInSofaFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=curv_abs_inputS,
                                          curv_abs_output=curv_abs_outputF, name='cosseratMapping',
                                          input1=self.cosseratCoordinateNode.cosseratCoordinateMO.getLinkPath(),
                                          input2=self.rigidBaseNode.RigidBaseMO.getLinkPath(),
                                          output=framesMO.getLinkPath(), debug=0, radius=0)

        if self.beamMass != 0.:
            self.solverNode.addObject('MechanicalMatrixMapper', template='Vec3,Rigid3',
                                      object1=self.cosseratCoordinateNode.cosseratCoordinateMO.getLinkPath(),
                                      object2=self.rigidBaseNode.RigidBaseMO.getLinkPath(),
                                      nodeToParse=cosseratInSofaFrameNode.getLinkPath())
        return cosseratInSofaFrameNode


def createScene(rootNode):
    rootNode.addObject('RequiredPlugin', name='plugins', pluginName=[pluginList,
                                                                     ['SofaEngine', 'SofaLoader', 'SofaSimpleFem',
                                                                      'SofaExporter']])
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideCollisionModels '
                                                   'hideBoundingCollisionModels hireForceFields '
                                                   'hideInteractionForceFields hideWireframe')
    rootNode.findData('dt').value = 0.01
    rootNode.findData('gravity').value = [0., -9.81, 0.]
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('GenericConstraintSolver', tolerance=1e-5, maxIterations=5e2)
    rootNode.addObject('Camera', position="-35 0 280", lookAt="0 0 0")

    solverNode = rootNode.addChild('solverNode')
    solverNode.addObject('EulerImplicitSolver', rayleighStiffness="0.2", rayleighMass='0.1')
    solverNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")
    solverNode.addObject('GenericConstraintCorrection')

    cosserat = solverNode.addChild(
        Cosserat(parent=solverNode, cosseratGeometry=cosserat_config, name="cosserat", radius=0.2))

    # use this to add the collision if the beam will interact with another object
    collisionModel = cosserat.addCollisionModel()

    # attach force at the beam tip,
    # we can attach this force to non mechanical node thanks to the MechanicalMatrixMapper component
    beamFrame = cosserat.cosseratFrame
    beamFrame.addObject('ConstantForceField', name='constForce', showArrowSize=1.e-2, indices=12,
                        force=[0., -100., 0., 0., 0., 0.])

    return rootNode
