# -*- coding: utf-8 -*-
"""
    Cosserat class in SofaPython3.
"""

__authors__ = "adagolodjo"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "October, 26 2021"

import Sofa
# from cosserat.usefulFunctions import buildEdges
from cosserat.utils import addEdgeCollision, addPointsCollision
from useful.header import addHeader, addVisual, addSolverNode
# from useful.params import Parameters
from useful.params import BeamPhysicsParameters, BeamGeometryParameters, SimulationParameters, ContactParameters

cosserat_config = {'init_pos': [0., 0., 0.], 'tot_length': 6, 'nbSectionS': 6,
                   'nbFramesF': 12, 'buildCollisionModel': 1, 'beamMass': 0.22}


class CosseratBase(Sofa.Prefab):

    """
    CosseratBase model prefab class. It is a prefab class that allow to create a cosserat beam/rod in Sofa.
           Structure:
           Node : {
                name : 'CosseratBase'
                Node0 MechanicalObject :     // Rigid position of the base of the beam
                Node1 MechanicalObject :    // Vec3d, cosserat local parameters composed of the twist and the bending along y and z
                Node1 ForceField          // Base on Hook's law, it compute the force applied on the beam
                (Node0-Node1)-child MechanicalObject     //  The child of the two precedent nodes, Rigid positions
                Allow to compute the cosserat frame in the world frame (Sofa frame)
                    Cosserat Mapping  //  it allow the transfer from the local to the word frame
            }
            params

    """
    prefabParameters = [
        {'name': 'name', 'type': 'string', 'help': 'Node name', 'default': 'Cosserat'},
        {'name': 'position', 'type': 'Rigid3d::VecCoord', 'help': 'Cosserat base position',
         'default': [[0., 0., 0., 0, 0, 0, 1.]]},
        {'name': 'translation', 'type': 'Vec3d', 'help': 'Cosserat base Rotation', 'default': [0., 0., 0.]},
        {'name': 'rotation', 'type': 'Vec3d',
            'help': 'Cosserat base Rotation', 'default': [0., 0., 0.]},
        {'name': 'youngModulus', 'type': 'double',
            'help': 'Beam Young modulus', 'default': 1.e6},
        {'name': 'poissonRatio', 'type': 'double',
            'help': 'Beam poisson ratio', 'default': 0.4},
        {'name': 'shape', 'type': 'string', 'help': 'beam section', 'default': "circular"},
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
        beamPhysicsParams = BeamPhysicsParameters()
        beamGeometryParams = BeamGeometryParameters(init_pos=[0., 0., 0.])
        # self.bGeometryParams = kwargs.get('beamGeometryParams')
        # params = Parameters()

        self.cosseratGeometry = kwargs['cosseratGeometry']
        self.beamMass = beamPhysicsParams.beamMass # self.cosseratGeometry['beamMass']
        self.parent = kwargs.get('parent')
        self.useInertiaParams = beamPhysicsParams.useInertia # False
        self.radius = beamPhysicsParams.beamRadius # kwargs.get('radius')

        if self.parent.hasObject("EulerImplicitSolver") is False:
            print('The code does not have parent EulerImplicit')
            self.solverNode = addSolverNode(self.parent)
        else:
            self.solverNode = self.parent

        if 'inertialParams' in kwargs:
            self.useInertiaParams = True
            self.inertialParams = kwargs['inertialParams']

        self.rigidBaseNode = self.addRigidBaseNode()

        cosserat_geometry = CosseratGeometry(beamGeometryParams)
        self.frames3D = cosserat_geometry.cable_positionF

        # [positionS, curv_abs_inputS, sectionsLength, framesF, curv_abs_outputF, self.frames3D] = __buildCosseratGeometry(params.beamGeoParams)

        self.cosseratCoordinateNode = self.addCosseratCoordinate(cosserat_geometry.bendingState,
                                                                 cosserat_geometry.sectionsLengthList)
        self.cosseratFrame = self.addCosseratFrame(cosserat_geometry.framesF, cosserat_geometry.curv_abs_inputS,
                                                   cosserat_geometry.curv_abs_outputF)

    def init(self):
        pass

    def addCollisionModel(self):
        tab_edges = generate_edge_list(self.frames3D)
        return addEdgeCollision(self.cosseratFrame, self.frames3D, tab_edges)

    def addPointCollisionModel(self, nodeName='CollisionPoints'):
        tab_edges = generate_edge_list(self.frames3D)
        return addPointsCollision(self.cosseratFrame, self.frames3D, tab_edges, nodeName)

    def addSlidingPoints(self):
        slidingPoint = self.cosseratFrame.addChild('slidingPoint')
        slidingPoint.addObject('MechanicalObject', name="slidingPointMO", position=self.frames3D,
                               showObject="0", showIndices="0")
        slidingPoint.addObject('IdentityMapping')
        return slidingPoint

    def addSlidingPointsWithContainer(self):
        slidingPoint = self.cosseratFrame.addChild('slidingPoint')
        slidingPoint.addObject("PointSetTopologyContainer")
        slidingPoint.addObject("PointSetTopologyModifier")
        slidingPoint.addObject('MechanicalObject', name="slidingPointMO", position=self.frames3D,
                               showObject="1", showIndices="0")
        slidingPoint.addObject('IdentityMapping')
        return slidingPoint

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
                                             length=listOfSectionsLength, radius=self.radius,
                                             youngModulus=self.youngModulus.value, poissonRatio=self.poissonRatio.value,
                                             rayleighStiffness=self.rayleighStiffness.value,
                                             lengthY=self.length_Y.value, lengthZ=self.length_Z.value)
        else:
            self._extracted_from_addCosseratCoordinate_15(
                cosseratCoordinateNode, listOfSectionsLength
            )
        return cosseratCoordinateNode

    # TODO Rename this here and in `addCosseratCoordinate`
    def _extracted_from_addCosseratCoordinate_15(self, cosseratCoordinateNode, listOfSectionsLength):
        GA = self.inertialParams['GA']
        GI = self.inertialParams['GI']
        EA = self.inertialParams['EA']
        EI = self.inertialParams['EI']
        print(f'{GA}')
        cosseratCoordinateNode.addObject('BeamHookeLawForceField', crossSectionShape=self.shape.value,
                                         length=listOfSectionsLength, radius=self.radius, useInertiaParams=True,
                                         GI=GI, GA=GA, EI=EI, EA=EA, rayleighStiffness=self.rayleighStiffness.value,
                                         lengthY=self.length_Y.value, lengthZ=self.length_Z.value)

    def addCosseratFrame(self, framesF, curv_abs_inputS, curv_abs_outputF):

        cosseratInSofaFrameNode = self.rigidBaseNode.addChild(
            'cosseratInSofaFrameNode')
        self.cosseratCoordinateNode.addChild(cosseratInSofaFrameNode)
        framesMO = cosseratInSofaFrameNode.addObject(
            'MechanicalObject', template='Rigid3d', name="FramesMO", position=framesF, showObject=1, showObjectScale=0.001)
        if self.beamMass != 0.:
            cosseratInSofaFrameNode.addObject(
                'UniformMass', totalMass=self.beamMass, showAxisSizeFactor='0')

        cosseratInSofaFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=curv_abs_inputS,
                                          curv_abs_output=curv_abs_outputF, name='cosseratMapping',
                                          input1=self.cosseratCoordinateNode.cosseratCoordinateMO.getLinkPath(),
                                          input2=self.rigidBaseNode.RigidBaseMO.getLinkPath(),
                                          output=framesMO.getLinkPath(), debug=0, radius=self.radius)
        return cosseratInSofaFrameNode


class CosseratGeometry:
    def __init__(self, beamGeoParams):
        # Data validation checks for beamGeoParams
        if not isinstance(beamGeoParams, BeamGeometryParameters):
            raise ValueError("beamGeoParams must be an instance of BeamGeoParams.")

        self.bendingState, self.curv_abs_inputS, self.sectionsLengthList = self.__calculate_beam_parameters(beamGeoParams)
        self.framesF, self.curv_abs_outputF, self.cable_positionF = self.__calculate_frame_parameters(beamGeoParams)

    def __calculate_beam_parameters(self, beamGeoParams):
        # Data validation checks for beamGeoParams attributes
        if not all(hasattr(beamGeoParams, attr) for attr in ['init_pos', 'beamLength', 'nbSection']):
            raise ValueError("beamGeoParams must have 'init_pos', 'beamLength', and 'nbSection' attributes.")

        total_length = beamGeoParams.beamLength
        nb_sections = beamGeoParams.nbSection
        x, y, z = beamGeoParams.init_pos

        if not all(isinstance(val, (int, float)) for val in [x, y, z, total_length]):
            raise ValueError("init_pos and beamLength in beamGeoParams must be numeric values.")

        if not isinstance(nb_sections, int) or nb_sections <= 0:
            raise ValueError("nbSection in beamGeoParams must be a positive integer.")

        length_s = total_length / nb_sections
        bendingState = []
        listOfSectionsLength = []
        temp = x
        curv_abs_input_s = [x]

        for i in range(nb_sections):
            bendingState.append([0, 0, 0])
            listOfSectionsLength.append((((i + 1) * length_s) - i * length_s))
            temp += listOfSectionsLength[i]
            curv_abs_input_s.append(temp)
        curv_abs_input_s[nb_sections] = total_length + x

        return bendingState, curv_abs_input_s, listOfSectionsLength

    def __calculate_frame_parameters(self, beamGeoParams):
        # Data validation checks for beamGeoParams attributes
        # if not all(hasattr(beamGeoParams, attr) for attr in ['init_pos', 'beamLength', 'nbFrames']):
        #     raise ValueError("beamGeoParams must have 'init_pos', 'beamLength', and 'nbFrames' attributes.")

        x, y, z = beamGeoParams.init_pos
        total_length = beamGeoParams.beamLength
        nb_frames = beamGeoParams.nbFrames

        if not all(isinstance(val, (int, float)) for val in [x, y, z, total_length]):
            raise ValueError("init_pos and beamLength in beamGeoParams must be numeric values.")

        if not isinstance(nb_frames, int) or nb_frames <= 0:
            raise ValueError("nbFrames in beamGeoParams must be a positive integer.")

        length_f = total_length / nb_frames
        frames_f = []
        curv_abs_output_f = []
        cable_position_f = []

        for i in range(nb_frames):
            sol = i * length_f
            frames_f.append([sol + x, y, z, 0, 0, 0, 1])
            cable_position_f.append([sol + x, y, z])
            curv_abs_output_f.append(sol + x)

        frames_f.append([total_length + x, y, z, 0, 0, 0, 1])
        cable_position_f.append([total_length + x, y, z])
        curv_abs_output_f.append(total_length + x)

        return frames_f, curv_abs_output_f, cable_position_f

from typing import List

def generate_edge_list(cable3DPos: List[List[float]]) -> List[int]:
    """
    Generate an edge list required in the EdgeSetTopologyContainer component.

    Parameters:
        cable3DPos (List[List[float]]): A list of 3D points representing the cable positions.

    Returns:
        List[int]: A list of indices forming edges in the EdgeSetTopologyContainer.
    """
    number_of_points = len(cable3DPos)
    edge_list = [i for i in range(number_of_points - 1) for _ in range(2)]
    return edge_list

def createScene(rootNode):
    addHeader(rootNode)
    addVisual(rootNode)

    rootNode.findData('dt').value = 0.01
    rootNode.findData('gravity').value = [0., -9.81, 0.]
    rootNode.addObject('BackgroundSetting', color='0 0.168627 0.211765')
    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('GenericConstraintSolver', tolerance=1e-5, maxIterations=5e2)
    rootNode.addObject('Camera', position="-35 0 280", lookAt="0 0 0")

    solverNode = rootNode.addChild('solverNode')
    solverNode.addObject('EulerImplicitSolver',
                         rayleighStiffness="0.2", rayleighMass='0.1')
    solverNode.addObject('SparseLDLSolver', name='solver',
                         template="CompressedRowSparseMatrixd")
    solverNode.addObject('GenericConstraintCorrection')

    cosserat = solverNode.addChild(
        CosseratBase(parent=solverNode, cosseratGeometry=cosserat_config, name="cosserat", radius=0.15))

    # use this to add the collision if the beam will interact with another object
    cosserat.addCollisionModel()

    # Attach a force at the beam tip,
    # we can attach this force to a non-mechanical node to control the beam in order to be able to drive it.
    beamFrame = cosserat.cosseratFrame
    beamFrame.addObject('ConstantForceField', name='constForce', showArrowSize=1.e-2, indices=12,
                        force=[0., -100., 0., 0., 0., 0.])

    return rootNode
