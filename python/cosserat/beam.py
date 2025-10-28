# -*- coding: utf-8 -*-
"""
Cosserat class in SofaPython3.

This module provides a prefab class to create and manipulate Cosserat beam/rod models in SOFA.
The CosseratBase class encapsulates the physics and geometry of a beam, handling
the creation of frames, coordinates, and physical properties needed for simulation.
"""

__authors__ = "adagolodjo"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "October, 26 2021"

import logging
from typing import List

import Sofa
from numpy import array

from .geometry import CosseratGeometry, generate_edge_list
from .header import addHeader, addVisual
from .params import (BeamGeometryParameters, BeamPhysicsBaseParameters,
                     Parameters)
from .utils import addEdgeCollision, addPointsCollision, create_rigid_node


class CosseratBase(Sofa.Prefab):
    """
    CosseratBase model prefab class that creates a cosserat beam/rod in SOFA.

    This class creates a complete beam model with the following structure:
        Node : {
            name : 'CosseratBase'
            Node0 MechanicalObject :     // Rigid position of the base of the beam
            Node1 MechanicalObject :     // Vec3d, cosserat local parameters composed of the twist and the bending along y and z
            Node1 ForceField :           // Based on Hooke's law, computes the forces applied on the beam
            (Node0-Node1)-child MechanicalObject :  // Child of the two precedent nodes, Rigid positions
                                                    // Allows computing the cosserat frame in the world frame (SOFA frame)
            Cosserat Mapping :           // Allows the transfer from the local to the world frame
        }

    Parameters:
        name (str): Node name for the CosseratBase prefab
        translation (numpy.ndarray): 3D vector representing the initial position of the beam base
        rotation (numpy.ndarray): 3D vector representing the initial orientation of the beam base
        params (Parameters): Physics and geometry parameters for the beam
        parent (Sofa.Node): Parent node in the SOFA scene graph
        inertialParams (Dict[str, Any], optional): Custom inertia parameters if needed
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
        },
    ]

    def __init__(self, *args, **kwargs):
        """
        Initialize the CosseratBase prefab with the given parameters.

        Args:
            *args: Variable length argument list passed to the parent class
            **kwargs: Arbitrary keyword arguments including:
                - params (Parameters): Beam physics and geometry parameters
                - parent (Sofa.Node): Parent node in the SOFA scene graph
                - inertialParams (Dict[str, Any], optional): Custom inertia parameters
        """
        Sofa.Prefab.__init__(self, *args, **kwargs)

        # Parameter validation
        if "params" not in kwargs:
            raise ValueError("The 'params' parameter is required for CosseratBase")
        if "parent" not in kwargs:
            raise ValueError("The 'parent' parameter is required for CosseratBase")

        self.params = kwargs.get("params")
        self.solverNode = kwargs.get("parent")

        # Initialize translation and rotation from prefab parameters or defaults
        # Try to get from Prefab parameters first, then from kwargs, then default
        translation_param = getattr(self, 'translation', kwargs.get("translation", [0.0, 0.0, 0.0]))
        rotation_param = getattr(self, 'rotation', kwargs.get("rotation", [0.0, 0.0, 0.0]))

        # Handle SOFA DataContainer objects and extract actual values
        if hasattr(translation_param, 'value'):
            self._translation_value = translation_param.value
        else:
            self._translation_value = translation_param

        if hasattr(rotation_param, 'value'):
            self._rotation_value = rotation_param.value
        else:
            self._rotation_value = rotation_param

        # Ensure they are lists (convert numpy arrays if needed)
        if hasattr(self._translation_value, 'tolist'):
            self._translation_value = self._translation_value.tolist()
        if hasattr(self._rotation_value, 'tolist'):
            self._rotation_value = self._rotation_value.tolist()

        # Extract physics parameters
        self.beam_physics_params = self.params.beam_physics_params
        self.beam_mass = self.beam_physics_params.beam_mass
        self.use_inertia_params = self.beam_physics_params.useInertia
        self.radius = self.beam_physics_params.beam_radius

        # Log parameters instead of print
        logging.info(f"The beam mass is: {self.beam_mass}")
        logging.info(f"The beam radius is: {self.radius}")

        # Override inertia params if provided
        if "inertialParams" in kwargs:
            self.use_inertia_params = True
            self.inertial_params = kwargs["inertialParams"]

        # Create the beam structure
        self.rigid_base_node = self._add_rigid_base_node()

        cosserat_geometry = CosseratGeometry(self.params.beam_geo_params)
        self.frames3D = cosserat_geometry.cable_positionF

        self.cosserat_coordinate_node = self._add_cosserat_coordinate(
            cosserat_geometry.bendingState, cosserat_geometry.sectionsLengthList
        )

        self.cosserat_frame = self._add_cosserat_frame(
            cosserat_geometry.framesF,
            cosserat_geometry.curv_abs_inputS,
            cosserat_geometry.curv_abs_outputF,
        )

    def __repr__(self) -> str:
        """
        Return a string representation of the CosseratBase object.

        Returns:
            str: A string representation including key properties
        """
        return (
            f"CosseratBase(name='{self.name}', "
                f"mass={self.beam_mass}, "
                f"radius={self.radius}, "
                f"use_inertia={self.use_inertia_params})"
        )

    def add_collision_model(self) -> "Sofa.Node":
        """
        Add an edge-based collision model to the cosserat beam.

        Returns:
            Sofa.Node: The created collision node
        """
        tab_edges = generate_edge_list(self.frames3D)
        return addEdgeCollision(self.cosserat_frame, self.frames3D, tab_edges)

    def _add_point_collision_model(
        self, node_name: str = "CollisionPoints"
    ) -> "Sofa.Core.Node":
        """
        Add a point-based collision model to the cosserat beam.

        Args:
            node_name: Name of the collision node

        Returns:
            Sofa.Node: The created collision node
        """
        tab_edges = generate_edge_list(self.frames3D)
        return addPointsCollision(
            self.cosserat_frame, self.frames3D, tab_edges, node_name
        )

    def _add_sliding_points(self) -> "Sofa.Core.Node":
        """
        Add sliding points to the cosserat frame.

        These points can be used for interaction or visualization.

        Returns:
            Sofa.Node: The created sliding point node
        """
        sliding_point = self.cosserat_frame.addChild("slidingPoint")
        sliding_point.addObject(
            "MechanicalObject", name="slidingPointMO", position=self.frames3D
        )
        sliding_point.addObject("IdentityMapping")
        return sliding_point

    def _add_sliding_points_with_container(self) -> "Sofa.Core.Node":
        """
        Add sliding points with topology container and modifier.

        This extends the basic sliding points with topology objects that
        allow modifying the point set during simulation.

        Returns:
            Sofa.Node: The created sliding point node with topology container
        """
        sliding_point = self._add_sliding_points()
        sliding_point.addObject("PointSetTopologyContainer")
        sliding_point.addObject("PointSetTopologyModifier")
        return sliding_point

    def _add_rigid_base_node(self) -> Sofa.Core.Node:
        """
        Create a rigid node at the base of the beam.

        This node defines the global position and orientation of the beam's base.

        Returns:
            Sofa.Node: The created rigid base node
        """
        rigid_base_node = create_rigid_node(
            self, "RigidBase", self._translation_value, self._rotation_value
        )
        return rigid_base_node

    def _add_cosserat_coordinate(
        self, initial_curvature: List[float], section_lengths: List[float]
    ):
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
            position=initial_curvature,
        )

        if not self.use_inertia_params:
            self._add_beam_hooke_law_without_inertia(
                cosserat_coordinate_node, section_lengths
            )
        else:
            self._add_beam_hooke_law_with_inertia(
                cosserat_coordinate_node, section_lengths
            )

        return cosserat_coordinate_node

    def _add_beam_hooke_law_without_inertia(
        self, cosserat_coordinate_node: "Sofa.Node", section_lengths: List[float]
    ) -> None:
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

    def _add_beam_hooke_law_with_inertia(
        self, cosserat_coordinate_node: "Sofa.Node", section_lengths: List[float]
    ) -> None:
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
        cosserat_coordinate_node.addObject(
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

    def _add_cosserat_frame(
        self,
        frames_f: List,
        curv_abs_input_s: List[float],
        curv_abs_output_f: List[float],
    ) -> "Sofa.Node":
        """
        Create the node that represents the cosserat frames in the SOFA world frame.

        This method creates the mapping between local cosserat coordinates and
        world frame rigid positions.

        Args:
            frames_f: List of frame positions
            curv_abs_input_s: Curvilinear abscissa input values
            curv_abs_output_f: Curvilinear abscissa output values

        Returns:
            Sofa.Node: The node containing the cosserat frames in SOFA world frame
        """
        cosserat_in_sofa_frame_node = self.rigid_base_node.addChild(
            "cosseratInSofaFrameNode"
        )
        self.cosserat_coordinate_node.addChild(cosserat_in_sofa_frame_node)
        frames_mo = cosserat_in_sofa_frame_node.addObject(
            "MechanicalObject",
            template="Rigid3d",
            name="FramesMO",
            position=frames_f,
        )

        cosserat_in_sofa_frame_node.addObject(
            "UniformMass",
            totalMass=self.beam_mass,
            showAxisSizeFactor="0",
        )

        cosserat_in_sofa_frame_node.addObject(
            "DiscreteCosseratMapping",
            curv_abs_input=curv_abs_input_s,
            curv_abs_output=curv_abs_output_f,
            name="cosseratMapping",
            input1=self.cosserat_coordinate_node.cosseratCoordinateMO.getLinkPath(),
            input2=self.rigid_base_node.RigidBaseMO.getLinkPath(),
            output=frames_mo.getLinkPath(),
            debug=0,
            radius=self.radius,
        )
        return cosserat_in_sofa_frame_node


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

    # Create a Cosserat beam
    cosserat = solverNode.addChild(CosseratBase(parent=solverNode, params=Params))
    cosserat.rigid_base_node.addObject(
        "RestShapeSpringsForceField",
        name="spring",
        stiffness=1e8,
        angularStiffness=1.0e8,
        external_points=0,
        # mstate="@RigidBaseMO",
        points=0,
        template="Rigid3d",
    )

    # use this to add the collision if the beam will interact with another object
    cosserat.add_collision_model()

    # Attach a force at the beam tip,
    # we can attach this force to a non-mechanical node to control the beam in order to be able to drive it.
    cosserat.cosserat_frame

    return rootNode
