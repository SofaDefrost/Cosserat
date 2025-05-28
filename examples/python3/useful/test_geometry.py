"""
Tests for the Cosserat Beam Geometry module.

This module contains tests for the CosseratGeometry class, ensuring that:
1. The class initializes properly with valid parameters
2. All getter methods return the expected values
3. Update methods properly modify the internal state
4. Edge cases and error conditions are handled correctly
5. Calculated values match expected results

Run with: pytest -xvs test_geometry.py
"""

import pytest
import numpy as np
from typing import List

from useful.params import BeamGeometryParameters
from useful.geometry import CosseratGeometry, calculate_beam_parameters, calculate_frame_parameters, generate_edge_list


@pytest.fixture
def default_params() -> BeamGeometryParameters:
    """
    Fixture that provides default beam geometry parameters.
    
    Returns:
        BeamGeometryParameters with default values
    """
    return BeamGeometryParameters()


@pytest.fixture
def custom_params() -> BeamGeometryParameters:
    """
    Fixture that provides custom beam geometry parameters.
    
    Returns:
        BeamGeometryParameters with custom values
    """
    return BeamGeometryParameters(
        beam_length=2.0,
        nb_section=10,
        nb_frames=20,
        build_collision_model=1
    )


@pytest.fixture
def default_geometry(default_params) -> CosseratGeometry:
    """
    Fixture that provides a CosseratGeometry instance with default parameters.
    
    Returns:
        CosseratGeometry instance with default parameters
    """
    return CosseratGeometry(default_params)


@pytest.fixture
def custom_geometry(custom_params) -> CosseratGeometry:
    """
    Fixture that provides a CosseratGeometry instance with custom parameters.
    
    Returns:
        CosseratGeometry instance with custom parameters
    """
    return CosseratGeometry(custom_params)


class TestCosseratGeometry:
    """Test suite for the CosseratGeometry class."""

    def test_initialization(self, default_params, default_geometry):
        """Test that the CosseratGeometry class initializes correctly with default parameters."""
        # Verify the parameters were stored
        assert default_geometry.params == default_params

        # Verify internal state was calculated
        assert hasattr(default_geometry, 'bendingState')
        assert hasattr(default_geometry, 'curv_abs_sections')
        assert hasattr(default_geometry, 'section_lengths')
        assert hasattr(default_geometry, 'frames')
        assert hasattr(default_geometry, 'curv_abs_frames')
        assert hasattr(default_geometry, 'cable_positions')
        assert hasattr(default_geometry, 'edge_list')

        # Verify lengths of calculated arrays
        assert len(default_geometry.bendingState) == default_params.nb_section
        assert len(default_geometry.curv_abs_sections) == default_params.nb_section + 1
        assert len(default_geometry.section_lengths) == default_params.nb_section
        assert len(default_geometry.frames) == default_params.nb_frames + 1
        assert len(default_geometry.curv_abs_frames) == default_params.nb_frames + 1
        assert len(default_geometry.cable_positions) == default_params.nb_frames + 1
        assert len(default_geometry.edge_list) == default_params.nb_frames

    def test_custom_initialization(self, custom_params, custom_geometry):
        """Test that the CosseratGeometry class initializes correctly with custom parameters."""
        # Verify the parameters were stored
        assert custom_geometry.params == custom_params

        # Verify internal state was calculated with correct sizes
        assert len(custom_geometry.bendingState) == custom_params.nb_section
        assert len(custom_geometry.curv_abs_sections) == custom_params.nb_section + 1
        assert len(custom_geometry.section_lengths) == custom_params.nb_section
        assert len(custom_geometry.frames) == custom_params.nb_frames + 1
        assert len(custom_geometry.curv_abs_frames) == custom_params.nb_frames + 1
        assert len(custom_geometry.cable_positions) == custom_params.nb_frames + 1
        assert len(custom_geometry.edge_list) == custom_params.nb_frames

    def test_get_beam_length(self, default_geometry, custom_geometry):
        """Test get_beam_length returns the correct beam length."""
        assert default_geometry.get_beam_length() == 1.0
        assert custom_geometry.get_beam_length() == 2.0

    def test_get_number_of_sections(self, default_geometry, custom_geometry):
        """Test get_number_of_sections returns the correct number of sections."""
        assert default_geometry.get_number_of_sections() == 5
        assert custom_geometry.get_number_of_sections() == 10

    def test_get_number_of_frames(self, default_geometry, custom_geometry):
        """Test get_number_of_frames returns the correct number of frames."""
        assert default_geometry.get_number_of_frames() == 30
        assert custom_geometry.get_number_of_frames() == 20

    def test_get_bending_state(self, default_geometry):
        """Test get_bending_state returns the correct bending state."""
        bending_state = default_geometry.get_bending_state()
        assert len(bending_state) == default_geometry.get_number_of_sections()
        # Initially all sections should have zero curvature
        for section in bending_state:
            assert section == [0.0, 0.0, 0.0]

    def test_get_curvilinear_abscissa_sections(self, default_geometry):
        """Test get_curvilinear_abscissa_sections returns the correct abscissa values."""
        abscissa = default_geometry.get_curvilinear_abscissa_sections()
        assert len(abscissa) == default_geometry.get_number_of_sections() + 1
        # First value should be 0
        assert abscissa[0] == 0.0
        # Last value should be beam length
        assert abscissa[-1] == default_geometry.get_beam_length()
        # Values should be evenly spaced
        section_length = default_geometry.get_beam_length() / default_geometry.get_number_of_sections()
        for i in range(len(abscissa) - 1):
            assert pytest.approx(abscissa[i + 1] - abscissa[i]) == section_length

    def test_get_section_lengths(self, default_geometry):
        """Test get_section_lengths returns the correct section lengths."""
        section_lengths = default_geometry.get_section_lengths()
        assert len(section_lengths) == default_geometry.get_number_of_sections()
        # All sections should have equal length
        expected_length = default_geometry.get_beam_length() / default_geometry.get_number_of_sections()
        for length in section_lengths:
            assert pytest.approx(length) == expected_length

    def test_get_frames(self, default_geometry):
        """Test get_frames returns the correct frame data."""
        frames = default_geometry.get_frames()
        assert len(frames) == default_geometry.get_number_of_frames() + 1
        # Each frame should have 7 components: [x, y, z, qx, qy, qz, qw]
        for frame in frames:
            assert len(frame) == 7
        # First frame should be at [0,0,0] with identity quaternion [0,0,0,1]
        assert frames[0] == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
        # Last frame should be at [beam_length,0,0] with identity quaternion
        assert frames[-1] == [default_geometry.get_beam_length(), 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]

    def test_get_curvilinear_abscissa_frames(self, default_geometry):
        """Test get_curvilinear_abscissa_frames returns the correct abscissa values."""
        abscissa = default_geometry.get_curvilinear_abscissa_frames()
        assert len(abscissa) == default_geometry.get_number_of_frames() + 1
        # First value should be 0
        assert abscissa[0] == 0.0
        # Last value should be beam length
        assert abscissa[-1] == default_geometry.get_beam_length()
        # Values should be evenly spaced
        frame_spacing = default_geometry.get_beam_length() / default_geometry.get_number_of_frames()
        for i in range(len(abscissa) - 1):
            assert pytest.approx(abscissa[i + 1] - abscissa[i]) == frame_spacing

    def test_get_cable_positions(self, default_geometry):
        """Test get_cable_positions returns the correct cable positions."""
        cable_positions = default_geometry.get_cable_positions()
        assert len(cable_positions) == default_geometry.get_number_of_frames() + 1
        # Each position should have 3 components: [x, y, z]
        for position in cable_positions:
            assert len(position) == 3
        # First position should be at [0,0,0]
        assert cable_positions[0] == [0.0, 0.0, 0.0]
        # Last position should be at [beam_length,0,0]
        assert cable_positions[-1] == [default_geometry.get_beam_length(), 0.0, 0.0]

    def test_get_edge_list(self, default_geometry):
        """Test get_edge_list returns the correct edge list."""
        edge_list = default_geometry.get_edge_list()
        assert len(edge_list) == default_geometry.get_number_of_frames()
        # Check edge connections
        for i, edge in enumerate(edge_list):
            assert edge == [i, i + 1]

    def test_update_bending_state(self, default_geometry):
        """Test update_bending_state correctly updates the bending state."""
        # Create a new bending state with non-zero curvature
        new_bending_state = []
        for i in range(default_geometry.get_number_of_sections()):
            new_bending_state.append([0.1 * i, 0.2 * i, 0.3 * i])
        
        # Update the bending state
        default_geometry.update_bending_state(new_bending_state)
        
        # Verify the update
        updated_state = default_geometry.get_bending_state()
        assert len(updated_state) == len(new_bending_state)
        for i in range(len(updated_state)):
            assert updated_state[i] == new_bending_state[i]

    def test_update_frames(self, default_geometry):
        """Test update_frames correctly updates frames and cable positions."""
        # Create new frames with modified positions and orientations
        new_frames = []
        for i in range(default_geometry.get_number_of_frames() + 1):
            s = i / default_geometry.get_number_of_frames() * default_geometry.get_beam_length()
            # Bend the beam into a quarter circle in the XY plane
            theta = np.pi/2 * (s / default_geometry.get_beam_length())
            x = np.sin(theta) * default_geometry.get_beam_length() / np.pi * 2
            y = (1 - np.cos(theta)) * default_geometry.get_beam_length() / np.pi * 2
            # Simple quaternion for rotation around Z axis
            qx, qy, qz = 0.0, 0.0, np.sin(theta/2)
            qw = np.cos(theta/2)
            new_frames.append([x, y, 0.0, qx, qy, qz, qw])
        
        # Update the frames
        default_geometry.update_frames(new_frames)
        
        # Verify the frames update
        updated_frames = default_geometry.get_frames()
        assert len(updated_frames) == len(new_frames)
        for i in range(len(updated_frames)):
            for j in range(7):  # 7 components per frame
                assert pytest.approx(updated_frames[i][j]) == new_frames[i][j]
        
        # Verify the cable positions update
        cable_positions = default_geometry.get_cable_positions()
        assert len(cable_positions) == len(new_frames)
        for i in range(len(cable_positions)):
            for j in range(3):  # 3 components per position
                assert pytest.approx(cable_positions[i][j]) == new_frames[i][j]
        
        # Verify the edge list update
        edge_list = default_geometry.get_edge_list()
        assert len(edge_list) == default_geometry.get_number_of_frames()
        for i, edge in enumerate(edge_list):
            assert edge == [i, i + 1]

    def test_to_dict(self, default_geometry):
        """Test to_dict returns a complete dictionary representation."""
        geo_dict = default_geometry.to_dict()
        
        # Check that all expected keys are present
        expected_keys = [
            "beam_length", "nb_section", "nb_frames", 
            "bendingState", "curv_abs_sections", "section_lengths",
            "frames", "curv_abs_frames", "cable_positions", "edge_list"
        ]
        for key in expected_keys:
            assert key in geo_dict
        
        # Check that values match the geometry properties
        assert geo_dict["beam_length"] == default_geometry.get_beam_length()
        assert geo_dict["nb_section"] == default_geometry.get_number_of_sections()
        assert geo_dict["nb_frames"] == default_geometry.get_number_of_frames()
        assert geo_dict["bendingState"] == default_geometry.get_bending_state()
        assert geo_dict["curv_abs_sections"] == default_geometry.get_curvilinear_abscissa_sections()
        assert geo_dict["section_lengths"] == default_geometry.get_section_lengths()
        assert geo_dict["frames"] == default_geometry.get_frames()
        assert geo_dict["curv_abs_frames"] == default_geometry.get_curvilinear_abscissa_frames()
        assert geo_dict["cable_positions"] == default_geometry.get_cable_positions()
        assert geo_dict["edge_list"] == default_geometry.get_edge_list()

    def test_invalid_bending_state_update(self, default_geometry):
        """Test that update_bending_state raises ValueError for invalid input."""
        # Create a bending state with wrong number of sections
        invalid_bending_state = [[0.0, 0.0, 0.0]] * (default_geometry.get_number_of_sections() + 1)
        
        # Update should raise ValueError
        with pytest.raises(ValueError) as excinfo:
            default_geometry.update_bending_state(invalid_bending_state)
        assert "Expected" in str(excinfo.value)
        assert "bending state vectors" in str(excinfo.value)

    def test_invalid_frames_update(self, default_geometry):
        """Test that update_frames raises ValueError for invalid input."""
        # Create frames with wrong number of frames
        invalid_frames = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]] * default_geometry.get_number_of_frames()
        
        # Update should raise ValueError
        with pytest.raises(ValueError) as excinfo:
            default_geometry.update_frames(invalid_frames)
        assert "Expected" in str(excinfo.value)
        assert "frames" in str(excinfo.value)


def test_calculate_beam_parameters():
    """Test the calculate_beam_parameters function directly."""
    params = BeamGeometryParameters(beam_length=1.0, nb_section=5)
    bending_state, curv_abs, section_lengths = calculate_beam_parameters(params)
    
    # Check dimensions
    assert len(bending_state) == params.nb_section
    assert len(curv_abs) == params.nb_section + 1
    assert len(section_lengths) == params.nb_section
    
    # Check values
    for section in bending_state:
        assert section == [0.0, 0.0, 0.0]
    
    assert curv_abs[0] == 0.0
    assert curv_abs[-1] == params.beam_length
    
    expected_section_length = params.beam_length / params.nb_section
    for length in section_lengths:
        assert pytest.approx(length) == expected_section_length


def test_calculate_frame_parameters():
    """Test the calculate_frame_parameters function directly."""
    params = BeamGeometryParameters(beam_length=1.0, nb_frames=10)
    frames, curv_abs, cable_positions = calculate_frame_parameters(params)
    
    # Check dimensions
    assert len(frames) == params.nb_frames + 1
    assert len(curv_abs) == params.nb_frames + 1
    assert len(cable_positions) == params.nb_frames + 1
    
    # Check values
    for i, frame in enumerate(frames):
        s = i / params.nb_frames * params.beam_length
        assert pytest.approx(frame[0]) == s  # x position
        assert frame[1:3] == [0.0, 0.0]  # y, z positions
        assert frame[3:6] == [0.0, 0.0, 0.0]  # qx, qy, qz
        assert frame[6] == 1.0  # qw
    
    assert curv_abs[0] == 0.0
    assert curv_abs[-1] == params.beam_length
    
    for i, pos in enumerate(cable_positions):
        s = i / params.nb_frames * params.beam_length
        assert pytest.approx(pos[0]) == s  # x position
        assert pos[1:3] == [0.0, 0.0]  # y, z positions


def test_generate_edge_list():
    """Test the generate_edge_list function directly."""
    # Create a list of positions
    positions = [[i, 0.0, 0.0] for i in range(5)]
    edge_list = generate_edge_list(positions)
    
    # Check result
    assert len(edge_list) == len(positions) - 1
    for i, edge in enumerate(edge_list):
        assert edge == [i, i + 1]
    
    # Test empty list case
    assert generate_edge_list([]) == []

