"""
Cosserat Beam Geometry Module
=============================

This module defines parameter classes for configuring Cosserat beam simulations.
Cosserat beam theory is an extension of classical beam theory that accounts for
micro-rotations and is particularly useful for modeling slender structures with
complex behaviors such as medical instruments, cables, and soft robotics components.

The parameters are organized into several dataclasses:
- BeamGeometryParameters: Defines the beam's physical dimensions and discretization
- CosseratGeometry: Manages the geometric representation of the beam

Mathematical Concepts:
- Curvilinear abscissa: The parametric coordinate along the beam's centerline
- Bending state: Local curvature vectors [kx, ky, kz] at beam sections
- Frames: Position and orientation (quaternion) along the beam

This module provides functions to:
1. Calculate beam section parameters (calculate_beam_parameters)
2. Calculate frame parameters for visualization (calculate_frame_parameters)
3. Generate edge lists for topological representation (generate_edge_list)
"""

from typing import Dict, List, Optional, Tuple, Union, cast

import numpy as np
from numpy.typing import NDArray

from useful.params import BeamGeometryParameters


def calculate_beam_parameters(beamGeoParams: BeamGeometryParameters) -> Tuple[List[List[float]], List[float], List[float]]:
    """
    Calculate beam section parameters based on geometry parameters.

    This function discretizes the beam into sections and calculates:
    1. The initial bending state (zero curvature initially)
    2. The curvilinear abscissa for each section node
    3. The length of each section

    Parameters:
        beamGeoParams: Geometry parameters defining beam dimensions and discretization.

    Returns:
        Tuple containing:
        - bendingState: List of [kx, ky, kz] curvature values for each section (initially zeros)
        - curv_abs_input_s: List of curvilinear abscissa values at section nodes
        - listOfSectionsLength: List of section lengths

    Raises:
        ValueError: If beam geometry parameters are invalid.
    """
    # Data validation checks for beamGeoParams attributes
    if not all(hasattr(beamGeoParams, attr) for attr in ['beam_length', 'nb_section']):
        raise ValueError("beamGeoParams must have 'beam_length' and 'nb_section' attributes.")

    total_length = beamGeoParams.beam_length
    nb_sections = beamGeoParams.nb_section

    if not isinstance(total_length, (int, float)) or total_length <= 0:
        raise ValueError(f"beam_length must be a positive number, got {total_length}")

    if not isinstance(nb_sections, int) or nb_sections <= 0:
        raise ValueError(f"nb_section must be a positive integer, got {nb_sections}")

    # Calculate section length
    length_s = total_length / nb_sections

    # Initialize lists
    bendingState: List[List[float]] = []
    listOfSectionsLength: List[float] = []
    temp = 0.0
    curv_abs_input_s: List[float] = [0.0]

    # Calculate for each section
    for i in range(nb_sections):
        # Initial bending state is zero curvature in all directions
        bendingState.append([0.0, 0.0, 0.0])

        # All sections have equal length
        section_length = length_s
        listOfSectionsLength.append(section_length)

        # Calculate curvilinear abscissa
        temp += section_length
        curv_abs_input_s.append(temp)

    # Ensure the final abscissa matches the total length exactly
    curv_abs_input_s[nb_sections] = total_length

    return bendingState, curv_abs_input_s, listOfSectionsLength

def calculate_frame_parameters(beamGeoParams: BeamGeometryParameters) -> Tuple[List[List[float]], List[float], List[List[float]]]:
    """
    Calculate frame parameters for visualization and computation.

    This function generates frames along the beam and calculates:
    1. The frame positions and orientations (as position + quaternion)
    2. The curvilinear abscissa for each frame
    3. The cable positions (x,y,z) for each frame

    Each frame consists of [x, y, z, qx, qy, qz, qw] where:
    - (x,y,z) is the position
    - (qx,qy,qz,qw) is the quaternion representing orientation

    Parameters:
        beamGeoParams: Geometry parameters defining beam dimensions and discretization.

    Returns:
        Tuple containing:
        - frames_f: List of [x, y, z, qx, qy, qz, qw] for each frame
        - curv_abs_output_f: List of curvilinear abscissa values at frame positions
        - cable_position_f: List of [x, y, z] positions for each frame

    Raises:
        ValueError: If beam geometry parameters are invalid.
    """
    # Data validation checks for beamGeoParams attributes
    if not all(hasattr(beamGeoParams, attr) for attr in ['beam_length', 'nb_frames']):
        raise ValueError("beamGeoParams must have 'beam_length' and 'nb_frames' attributes.")

    total_length = beamGeoParams.beam_length
    nb_frames = beamGeoParams.nb_frames

    if not isinstance(total_length, (int, float)) or total_length <= 0:
        raise ValueError(f"beam_length must be a positive number, got {total_length}")

    if not isinstance(nb_frames, int) or nb_frames <= 0:
        raise ValueError(f"nb_frames must be a positive integer, got {nb_frames}")

    # Calculate frame spacing
    length_f = total_length / nb_frames

    # Initialize frame data structures
    frames_f: List[List[float]] = []
    curv_abs_output_f: List[float] = []
    cable_position_f: List[List[float]] = []

    # Generate frames along the beam
    for i in range(nb_frames):
        # Calculate curvilinear abscissa for this frame
        sol = i * length_f

        # Create frame with position [sol,0,0] and identity quaternion [0,0,0,1]
        frames_f.append([sol, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0])

        # Create cable position (initially straight along X-axis)
        cable_position_f.append([sol, 0.0, 0.0])

        # Store curvilinear abscissa
        curv_abs_output_f.append(sol)

    # Add the final frame at the end of the beam
    frames_f.append([total_length, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0])
    cable_position_f.append([total_length, 0.0, 0.0])
    curv_abs_output_f.append(total_length)

    return frames_f, curv_abs_output_f, cable_position_f

def generate_edge_list(cable3DPos: List[List[float]]) -> List[List[int]]:
    """
    Generate an edge list required in the EdgeSetTopologyContainer component.

    This function creates connectivity information between adjacent points,
    allowing for visualization and simulation of the beam as a series of
    connected segments.

    Parameters:
        cable3DPos: A list of 3D points [x,y,z] representing the cable positions.

    Returns:
        A list of index pairs [start_idx, end_idx] defining edges between adjacent points.
        For example, [[0,1], [1,2], ...] connects point 0 to 1, point 1 to 2, etc.
    """
    if not cable3DPos:
        return []

    number_of_points = len(cable3DPos)
    edges: List[List[int]] = []

    for i in range(number_of_points - 1):
        edges.append([i, i + 1])

    return edges

class CosseratGeometry:
    """
    A class that encapsulates the geometric aspects of a Cosserat beam.

    This class handles:
    - Section discretization for physics modeling
    - Frame generation for visualization and interaction
    - Maintaining curvilinear abscissa values
    - Access to geometric properties and transformations

    Attributes:
        bendingState: List of [kx, ky, kz] curvature vectors at each section
        curv_abs_sections: List of curvilinear abscissa values at section nodes
        section_lengths: List of section lengths
        frames: List of [x, y, z, qx, qy, qz, qw] for each frame (position + quaternion)
        curv_abs_frames: List of curvilinear abscissa values at frame positions
        cable_positions: List of [x, y, z] positions for each frame
        edge_list: List of index pairs [start_idx, end_idx] defining topology
        params: The beam geometry parameters used to initialize this object
    """

    def __init__(self, params: BeamGeometryParameters):
        """
        Initialize the CosseratGeometry with beam geometry parameters.

        Parameters:
            params: Geometry parameters defining beam dimensions and discretization.
        """
        self.params = params
        
        # Calculate and store beam section parameters
        self.bendingState, self.curv_abs_sections, self.section_lengths = calculate_beam_parameters(params)
        
        # Calculate and store frame parameters
        self.frames, self.curv_abs_frames, self.cable_positions = calculate_frame_parameters(params)
        
        # Generate and store edge list for topology
        self.edge_list = generate_edge_list(self.cable_positions)
    
    def get_beam_length(self) -> float:
        """
        Get the total length of the beam.

        Returns:
            The beam length in model units.
        """
        return self.params.beam_length
    
    def get_number_of_sections(self) -> int:
        """
        Get the number of sections used in beam discretization.

        Returns:
            The number of sections.
        """
        return self.params.nb_section
    
    def get_number_of_frames(self) -> int:
        """
        Get the number of frames used for visualization.

        Returns:
            The number of frames.
        """
        return self.params.nb_frames
    
    def get_bending_state(self) -> List[List[float]]:
        """
        Get the current bending state of the beam.

        Returns:
            List of [kx, ky, kz] curvature vectors for each section.
        """
        return self.bendingState
    
    def get_curvilinear_abscissa_sections(self) -> List[float]:
        """
        Get the curvilinear abscissa values at section nodes.

        Returns:
            List of curvilinear abscissa values.
        """
        return self.curv_abs_sections
    
    def get_section_lengths(self) -> List[float]:
        """
        Get the lengths of each section.

        Returns:
            List of section lengths.
        """
        return self.section_lengths
    
    def get_frames(self) -> List[List[float]]:
        """
        Get the frames (position + orientation) along the beam.

        Returns:
            List of [x, y, z, qx, qy, qz, qw] for each frame.
        """
        return self.frames
    
    def get_curvilinear_abscissa_frames(self) -> List[float]:
        """
        Get the curvilinear abscissa values at frame positions.

        Returns:
            List of curvilinear abscissa values.
        """
        return self.curv_abs_frames
    
    def get_cable_positions(self) -> List[List[float]]:
        """
        Get the cable positions for each frame.

        Returns:
            List of [x, y, z] positions.
        """
        return self.cable_positions
    
    def get_edge_list(self) -> List[List[int]]:
        """
        Get the edge list defining beam topology.

        Returns:
            List of index pairs [start_idx, end_idx] defining connectivity.
        """
        return self.edge_list
    
    def update_bending_state(self, new_bending_state: List[List[float]]) -> None:
        """
        Update the bending state of the beam.

        Parameters:
            new_bending_state: New list of [kx, ky, kz] curvature vectors.
            
        Raises:
            ValueError: If the length of new_bending_state doesn't match the number of sections.
        """
        if len(new_bending_state) != self.params.nb_section:
            raise ValueError(f"Expected {self.params.nb_section} bending state vectors, got {len(new_bending_state)}")
        self.bendingState = new_bending_state
    
    def update_frames(self, new_frames: List[List[float]]) -> None:
        """
        Update the frames along the beam.

        Parameters:
            new_frames: New list of [x, y, z, qx, qy, qz, qw] for each frame.
            
        Raises:
            ValueError: If the length of new_frames doesn't match the number of frames + 1.
        """
        expected_length = self.params.nb_frames + 1  # Account for the frame at the end
        if len(new_frames) != expected_length:
            raise ValueError(f"Expected {expected_length} frames, got {len(new_frames)}")
        self.frames = new_frames
        
        # Also update cable positions based on the new frame positions
        self.cable_positions = [[frame[0], frame[1], frame[2]] for frame in new_frames]
        
        # Re-generate edge list if cable positions changed
        self.edge_list = generate_edge_list(self.cable_positions)
    
    def to_dict(self) -> Dict:
        """
        Convert the geometry data to a dictionary.

        Useful for serialization or inspection.

        Returns:
            Dictionary containing all geometry data.
        """
        return {
            "beam_length": self.params.beam_length,
            "nb_section": self.params.nb_section,
            "nb_frames": self.params.nb_frames,
            "bendingState": self.bendingState,
            "curv_abs_sections": self.curv_abs_sections,
            "section_lengths": self.section_lengths,
            "frames": self.frames,
            "curv_abs_frames": self.curv_abs_frames,
            "cable_positions": self.cable_positions,
            "edge_list": self.edge_list
        }
