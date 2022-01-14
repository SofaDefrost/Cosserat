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

from cosserat.cosseratObject import Cosserat

linearConfig = {'init_pos': [0., 0., 0.], 'tot_length': 1, 'nbSectionS': 15,
                'nbFramesF': 30, 'buildCollisionModel': 1, 'beamMass': 0.22}

nonLinearConfig = {'init_pos': [0., 0., 0.], 'tot_length': 1, 'nbSectionS': 8,
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


class NonLinearCosserat(Sofa.Prefab):
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
                                                              'to follow arm position', 'default': '1'}]

    def __init__(self, *args, **kwargs):
        Sofa.Prefab.__init__(self, *args, **kwargs)
        self.cosseratGeometry = kwargs['cosseratGeometry']
        self.needCollisionModel = kwargs['useCollisionModel']
        self.beamMass = self.cosseratGeometry['beamMass']
        self.parent = kwargs['parent']
        self.legendreControlPos = kwargs['legendreControlPoints']

        if self.parent.hasObject("EulerImplicitSolver") is False:
            # print("===> The EulerImplicit is not in the node Yet ")
            self.solverNode = self.addSolverNode()
        else:
            self.solverNode = self.parent
            # print("===> The EulerImplicit is in the node Yet ")
        # self.solverNode = self.parent
        self.rigidBaseNode = self.addRigidBaseNode()
        [positionS, curv_abs_inputS, sectionLength, framesF, curv_abs_outputF, frames3D] = \
            BuildCosseratGeometry(self.cosseratGeometry)
        self.legendreControlPointsNode = self.addLegendrePolynomialsNode()
        self.cosseratCoordinateNode = self.addCosseratCoordinate(positionS, sectionLength, curv_abs_inputS)
        self.cosseratFrame = self.addCosseratFrame(framesF, curv_abs_inputS, curv_abs_outputF)
        if self.needCollisionModel:
            tab_edges = buildEdges(frames3D)
            self.cosseratFrameCollision = addEdgeCollision(self.cosseratFrame, frames3D, tab_edges)

    def init(self):
        pass

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
                                position=positions, rotation=rot, showObject=1)
        # one can choose to set this to false and directly attach the beam base
        # to a control object in order to be able to drive it.
        if int(self.attachingToLink.value):
            rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring',
                                    stiffness=1e8, angularStiffness=1.e8, external_points=0,
                                    mstate="@RigidBaseMO", points=0, template="Rigid3d")
        return rigidBaseNode

    def addLegendrePolynomialsNode(self):
        legendreControlPointsNode = self.solverNode.addChild('legendreControlPointsNode')
        # legendreControlPointsNode.addObject('EulerImplicitSolver', rayleighStiffness="0.2", rayleighMass='0.1')
        # legendreControlPointsNode.addObject('CGLinearSolver', tolerance=1.e-12, iterations=1, threshold=1.e-18)
        legendreControlPointsNode.addObject('MechanicalObject',
                                            template='Vec3d', name='legendreControlPointsMO',
                                            position=self.legendreControlPos, rest_position=self.legendreControlPos,
                                            showIndices=0)
        return legendreControlPointsNode

    def addCosseratCoordinate(self, positionS, longeurS, curv_abs_inputS):
        cosseratCoordinateNode = self.legendreControlPointsNode.addChild('cosseratCoordinate')
        positionXi = [[0., 0., 0.] for _ in range(len(curv_abs_inputS) - 1)]
        cosseratCoordinateNode.addObject('MechanicalObject',
                                         template='Vec3d', name='cosseratCoordinateMO', position=positionXi,
                                         showIndices=0)
        cosseratCoordinateNode.addObject('BeamHookeLawForceField', crossSectionShape=self.shape.value,
                                         length=longeurS,
                                         youngModulus=self.youngModulus.value,
                                         poissonRatio=self.poissonRatio.value,
                                         radius=self.radius.value,
                                         lengthY=self.length_Y.value, lengthZ=self.length_Z.value)
        # print(f'the curv_abs_inputS is : {curv_abs_inputS}')
        # print(f'the length is : {longeurS}')
        localCurv = curv_abs_inputS
        # localCurv.pop(0)
        controlPointsAbs = [0.3333333333333333, 0.6666666666666666, 1.0]
        cosseratCoordinateNode.addObject('LegendrePolynomialsMapping', curvAbscissa=localCurv, order=3,
                                         controlPointsAbs=controlPointsAbs, applyRestPosition=True)
        return cosseratCoordinateNode

    def addCosseratFrame(self, framesF, curv_abs_inputS, curv_abs_outputF):

        cosseratInSofaFrameNode = self.rigidBaseNode.addChild('cosseratInSofaFrameNode')
        self.cosseratCoordinateNode.addChild(cosseratInSofaFrameNode)
        framesMO = cosseratInSofaFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                                     name="FramesMO", position=framesF,
                                                     showObject=1, showObjectScale=0.05)
        # print(f'curvAbs inside frame :{curv_abs_inputS}')
        cosseratInSofaFrameNode.addObject('UniformMass', totalMass=self.beamMass, showAxisSizeFactor='0')
        cosseratInSofaFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=curv_abs_inputS,
                                          curv_abs_output=curv_abs_outputF, name='cosseratMapping',
                                          input1=self.cosseratCoordinateNode.cosseratCoordinateMO.getLinkPath(),
                                          input2=self.rigidBaseNode.RigidBaseMO.getLinkPath(),
                                          output=framesMO.getLinkPath(), debug=0, radius=0)

        self.solverNode.addObject('MechanicalMatrixMapper', template='Vec3,Rigid3',
                                  object1=self.cosseratCoordinateNode.cosseratCoordinateMO.getLinkPath(),
                                  object2=self.rigidBaseNode.RigidBaseMO.getLinkPath(),
                                  nodeToParse=cosseratInSofaFrameNode.getLinkPath())
        return cosseratInSofaFrameNode


initialStrain1 = [[0., 0., 0], [0., 0., 0], [0., 0., 0]]
initialStrain2 = [[0., 0., -0.52475341], [0., 0., -0.3098944], [0., 0., -0.10211416]]
initialStrain3 = [[0., 0., -0.96779204], [0., 0., -0.55894208], [0., 0., -0.18167142]]
# initialStrain4 = [[0., 0., -0.96770587], [0., 0., -0.55875284], [0., 0., -0.18155108]]
initialStrain4 = [[0., 0., 0], [0., 0., 0], [0., 0., 0]]


def createScene(rootNode):
    rootNode.addObject('RequiredPlugin', name='plugins', pluginName=[pluginList,
                                                                     ['SofaEngine', 'SofaLoader', 'SofaSimpleFem',
                                                                      'SofaExporter']])
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideCollisionModels '
                                                   'hideBoundingCollisionModels hireForceFields '
                                                   'hideInteractionForceFields hideWireframe')
    rootNode.findData('dt').value = 0.01
    # rootNode.findData('gravity').value = [0., -9.81, 0.]
    rootNode.findData('gravity').value = [0., 0., 0.]
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    # rootNode.addObject('FreeMotionAnimationLoop')
    # rootNode.addObject('GenericConstraintSolver', tolerance=1e-5, maxIterations=5e2)
    rootNode.addObject('Camera', position="-35 0 280", lookAt="0 0 0")

    solverNode = rootNode.addChild('solverNode')
    solverNode.addObject('EulerImplicitSolver', rayleighStiffness="0.2", rayleighMass='0.1')
    solverNode.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")
    # solverNode.addObject('SparseLUSolver', name='solver', template="CompressedRowSparseMatrixd")
    # solverNode.addObject('CGLinearSolver', tolerance=1.e-12, iterations=1000, threshold=1.e-18)

    needCollisionModel = 0  # use this if the collision model if the beam will interact with another object
    nonLinearCosserat = solverNode.addChild(
        NonLinearCosserat(parent=solverNode, cosseratGeometry=nonLinearConfig, useCollisionModel=needCollisionModel,
                          name="cosserat", radius=0.1, legendreControlPoints=initialStrain4))
    cosseratNode = nonLinearCosserat.legendreControlPointsNode
    cosseratNode.addObject('MechanicalMatrixMapper', template='Vec3,Vec3',
                           object1=cosseratNode.getLinkPath(),
                           object2=cosseratNode.getLinkPath(),
                           name='cosseratCoordinateNodeMapper',
                           nodeToParse=nonLinearCosserat.cosseratCoordinateNode.getLinkPath())

    beamFrame = nonLinearCosserat.cosseratFrame
    beamFrame.addObject('ConstantForceField', name='constForce', showArrowSize=1.e-8, indices=12,
                        force=[0., 0., 0., 0., 0., 450.])
    return rootNode
