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
from useful.utils import addEdgeCollision, addPointsCollision, _create_rigid_node
from useful.header import addHeader, addVisual, addSolverNode
from useful.params import Parameters, BeamGeometryParameters
from useful.geometry import CosseratGeometry, generate_edge_list
from numpy import array
from typing import List



class CosseratBase(Sofa.Prefab):
    """
    CosseratBase model prefab class. It is a prefab class that allow to create a cosserat beam/rod in Sofa.
           Structure:
           Node : {
                name : 'CosseratBase'
                Node0 MechanicalObject :     // Rigid position of the base of the beam
                Node1 MechanicalObject :    // Vec3d, cosserat local parameters composed of the twist and the bending along y and z
                Node1 ForceField          // Base on Hook's law, it computed the force applied on the beam
                (Node0-Node1)-child MechanicalObject     //  The child of the two precedent nodes, Rigid positions
                Allow to compute the cosserat frame in the world frame (Sofa frame)
                    Cosserat Mapping  //  it allow the transfer from the locial to the word frame
            }
            params

    """

    prefabParameters = [
        {"name": "name", "type": "string", "help": "Node name", "default": "Cosserat"},
        {
            "name": "translation",
            "type": "Vec3d",
            "help": "Cosserat base Rotation",
            "default": array([0.0, 0.0, 0.0]),
        },
        {
            "name": "rotation",
            "type": "Vec3d",
            "help": "Cosserat base Rotation",
            "default": array([0.0, 0.0, 0.0]),
        }
    ]

    def __init__(self, *args, **kwargs):
        Sofa.Prefab.__init__(self, *args, **kwargs)
        self.params = kwargs.get(
            "beam_params", Parameters()
        )  # Use the Parameters class with default values

        beamPhysicsParams = self.params.beam_physics_params
        self.beam_mass = beamPhysicsParams.beam_mass  # self.cosseratGeometry['beamMass']
        self.use_inertia_params = beamPhysicsParams.useInertia  # False
        self.radius = beamPhysicsParams.beam_radius  # kwargs.get('radius')

        self.solverNode = kwargs.get("parent")

        if "inertialParams" in kwargs:
            self.use_inertia_params = True
            self.inertialParams = kwargs["inertialParams"]

        self.rigidBaseNode = self._addRigidBaseNode()

        cosserat_geometry = CosseratGeometry(self.params.beam_geo_params)
        self.frames3D = cosserat_geometry.cable_positionF

        self.cosseratCoordinateNode = self._add_cosserat_coordinate(
            cosserat_geometry.bendingState, cosserat_geometry.sectionsLengthList
        )

        self.cosseratFrame = self._addCosseratFrame(
            cosserat_geometry.framesF,
            cosserat_geometry.curv_abs_inputS,
            cosserat_geometry.curv_abs_outputF,
        )

    def init(self):
        pass

    def addCollisionModel(self):
        tab_edges = generate_edge_list(self.frames3D)
        return addEdgeCollision(self.cosseratFrame, self.frames3D, tab_edges)

    def _addPointCollisionModel(self, nodeName="CollisionPoints"):
        tab_edges = generate_edge_list(self.frames3D)
        return addPointsCollision(
            self.cosseratFrame, self.frames3D, tab_edges, nodeName
        )

    def _addSlidingPoints(self):
        slidingPoint = self.cosseratFrame.addChild("slidingPoint")
        slidingPoint.addObject("MechanicalObject", name="slidingPointMO", position=self.frames3D)
        slidingPoint.addObject("IdentityMapping")
        return slidingPoint

    def _addSlidingPointsWithContainer(self):
        slidingPoint = self._addSlidingPoints()
        slidingPoint.addObject("PointSetTopologyContainer")
        slidingPoint.addObject("PointSetTopologyModifier")
        return slidingPoint

    def _addRigidBaseNode(self):
        rigidBaseNode = _create_rigid_node(self, "RigidBase",
                           self.translation, self.rotation)
        return rigidBaseNode

    def _add_cosserat_coordinate(self, initial_curvature: List[float], section_lengths: List[float]) -> None:
        """
        Adds a cosserat coordinate node with a BeamHookeLawForceField object to the graph.

        Args:
            initial_curvature: Initial curvature of the cosserat coordinate.
            section_lengths: Length of each section in the cosserat coordinate.

        Returns:
            The cosserat coordinate node added to the model.
        """
        cosserat_coordinate_node = self.addChild("cosseratCoordinate")
        cosserat_coordinate_node.addObject(
            "MechanicalObject",
            template="Vec3d",
            name="cosseratCoordinateMO",
            position=initial_curvature
        )

        if not self.use_inertia_params:
            self._add_beam_hooke_law_without_inertia(cosserat_coordinate_node, section_lengths)
        else:
            self._add_beam_hooke_law_with_inertia(cosserat_coordinate_node, section_lengths)

        return cosserat_coordinate_node

    def _add_beam_hooke_law_without_inertia(self, cosserat_coordinate_node: None,
                                            section_lengths: list[float]) -> None:
        """
        Adds a BeamHookeLawForceField object to the cosserat coordinate node without inertia parameters.

        Args:
            cosserat_coordinate_node: The cosserat coordinate node to add the object to.
            section_lengths: Length of each section in the cosserat coordinate.
        """
        cosserat_coordinate_node.addObject(
            "BeamHookeLawForceField",
            crossSectionShape=self.params.beam_physics_params.beam_shape,
            length=section_lengths,
            radius=self.params.beam_physics_params.beam_radius,
            youngModulus=self.params.beam_physics_params.young_modulus,
            poissonRatio=self.params.beam_physics_params.poisson_ratio,
            rayleighStiffness=self.params.simu_params.rayleigh_stiffness,
            lengthY=self.params.beam_physics_params.length_Y,
            lengthZ=self.params.beam_physics_params.length_Z,
        )

    def _add_beam_hooke_law_with_inertia(self, cosserat_coordinate_node: None, section_lengths: List[float]) -> None:
        """
        Adds a BeamHookeLawForceField object to the cosserat coordinate node with inertia parameters.

        Args:
            cosserat_coordinate_node: The cosserat coordinate node to add the object to.
            section_lengths: Length of each section in the cosserat coordinate.
        """
        GA = self.params.beam_physics_params.GA
        GI = self.params.beam_physics_params.GI
        EA = self.params.beam_physics_params.EA
        EI = self.params.beam_physics_params.EI
        cosseratCoordinateNode.addObject(
            "BeamHookeLawForceField",
            crossSectionShape=self.params.beam_physics_params.beam_shape,
            length=section_lengths,
            radius=self.params.beam_physics_params.beam_radius,
            useInertiaParams=True,
            GI=GI,
            GA=GA,
            EI=EI,
            EA=EA,
            rayleighStiffness=self.params.simu_params.rayleigh_stiffness,
            lengthY=self.params.beam_physics_params.length_Y,
            lengthZ=self.params.beam_physics_params.length_Z,
        )

    # TODO Rename this here and in `addCosseratCoordinate`


    def _addCosseratFrame(self, framesF, curv_abs_inputS, curv_abs_outputF):
        cosseratInSofaFrameNode = self.rigidBaseNode.addChild("cosseratInSofaFrameNode")
        self.cosseratCoordinateNode.addChild(cosseratInSofaFrameNode)
        framesMO = cosseratInSofaFrameNode.addObject(
            "MechanicalObject",
            template="Rigid3d",
            name="FramesMO",
            position=framesF
        )

        cosseratInSofaFrameNode.addObject(
            "UniformMass", totalMass=self.beam_mass, showAxisSizeFactor="0"
        )

        cosseratInSofaFrameNode.addObject(
            "DiscreteCosseratMapping",
            curv_abs_input=curv_abs_inputS,
            curv_abs_output=curv_abs_outputF,
            name="cosseratMapping",
            input1=self.cosseratCoordinateNode.cosseratCoordinateMO.getLinkPath(),
            input2=self.rigidBaseNode.RigidBaseMO.getLinkPath(),
            output=framesMO.getLinkPath(),
            debug=0,
            radius=self.radius,
        )
        return cosseratInSofaFrameNode


Params = Parameters(beam_geo_params=BeamGeometryParameters())


def createScene(rootNode):
    addHeader(rootNode)
    addVisual(rootNode)

    rootNode.findData("dt").value = 0.01
    rootNode.findData("gravity").value = [0.0, -9.81, 0.0]
    rootNode.addObject("BackgroundSetting", color="0 0.168627 0.211765")
    rootNode.addObject("FreeMotionAnimationLoop")
    rootNode.addObject("GenericConstraintSolver", tolerance=1e-5, maxIterations=5e2)
    rootNode.addObject("Camera", position="-35 0 280", lookAt="0 0 0")

    solverNode = rootNode.addChild("solverNode")
    solverNode.addObject(
        "EulerImplicitSolver", rayleighStiffness="0.2", rayleighMass="0.1"
    )
    solverNode.addObject(
        "SparseLDLSolver", name="solver", template="CompressedRowSparseMatrixd"
    )
    solverNode.addObject("GenericConstraintCorrection")

    # Create a
    cosserat = solverNode.addChild(CosseratBase(parent=solverNode, beam_params=Params))
    cosserat.rigidBaseNode.addObject(
            "RestShapeSpringsForceField",
            name="spring",
            stiffness=1e8,
            angularStiffness=1.0e8,
            external_points=0,
            # mstate="@RigidBaseMO",
            points=0,
            template="Rigid3d"
        )


    # use this to add the collision if the beam will interact with another object
    cosserat.addCollisionModel()

    # Attach a force at the beam tip,
    # we can attach this force to a non-mechanical node to control the beam in order to be able to drive it.
    cosserat.cosseratFrame

    return rootNode
