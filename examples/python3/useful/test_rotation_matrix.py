import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal, assert_allclose
from compute_rotation_matrix import (
    rotation_matrix_x,
    rotation_matrix_y,
    rotation_matrix_z,
    compute_rotation_matrix,
    euler_angles_from_rotation_matrix,
    _validate_angle
)


# Helper functions for tests
def is_rotation_matrix(R, tolerance=1e-6):
    """Check if a matrix is a valid rotation matrix."""
    # Check if determinant is 1
    det_check = abs(np.linalg.det(R) - 1.0) < tolerance
    
    # Check if matrix is orthogonal (R.T @ R = I)
    ortho_check = np.allclose(R.T @ R, np.eye(3), atol=tolerance)
    
    return det_check and ortho_check


class TestBasicRotationMatrices:
    """Tests for basic rotation matrices around X, Y, and Z axes."""
    
    def test_rotation_matrix_x_90_degrees(self):
        """Test rotation matrix around X axis for 90 degrees."""
        angle = np.pi / 2
        R = rotation_matrix_x(angle)
        
        # Known result for 90-degree rotation around X
        expected = np.array([
            [1, 0, 0],
            [0, 0, -1],
            [0, 1, 0]
        ])
        
        assert_array_almost_equal(R, expected)
        
    def test_rotation_matrix_y_90_degrees(self):
        """Test rotation matrix around Y axis for 90 degrees."""
        angle = np.pi / 2
        R = rotation_matrix_y(angle)
        
        # Known result for 90-degree rotation around Y
        expected = np.array([
            [0, 0, 1],
            [0, 1, 0],
            [-1, 0, 0]
        ])
        
        assert_array_almost_equal(R, expected)
        
    def test_rotation_matrix_z_90_degrees(self):
        """Test rotation matrix around Z axis for 90 degrees."""
        angle = np.pi / 2
        R = rotation_matrix_z(angle)
        
        # Known result for 90-degree rotation around Z
        expected = np.array([
            [0, -1, 0],
            [1, 0, 0],
            [0, 0, 1]
        ])
        
        assert_array_almost_equal(R, expected)
    
    def test_rotation_matrix_x_identity(self):
        """Test that 0 angle rotation is identity."""
        R = rotation_matrix_x(0)
        assert_array_equal(R, np.eye(3))

    def test_rotation_matrix_y_identity(self):
        """Test that 0 angle rotation is identity."""
        R = rotation_matrix_y(0)
        assert_array_equal(R, np.eye(3))

    def test_rotation_matrix_z_identity(self):
        """Test that 0 angle rotation is identity."""
        R = rotation_matrix_z(0)
        assert_array_equal(R, np.eye(3))
    
    def test_rotation_matrix_x_360_degrees(self):
        """Test that 360-degree rotation returns to identity."""
        R = rotation_matrix_x(2 * np.pi)
        assert_array_almost_equal(R, np.eye(3))

    def test_rotation_matrix_y_360_degrees(self):
        """Test that 360-degree rotation returns to identity."""
        R = rotation_matrix_y(2 * np.pi)
        assert_array_almost_equal(R, np.eye(3))

    def test_rotation_matrix_z_360_degrees(self):
        """Test that 360-degree rotation returns to identity."""
        R = rotation_matrix_z(2 * np.pi)
        assert_array_almost_equal(R, np.eye(3))
        

class TestCompositeRotationMatrices:
    """Tests for composite rotation matrices."""
    
    def test_composite_rotation_properties(self):
        """Test that composite rotations maintain rotation matrix properties."""
        angles = [(0.1, 0.2, 0.3), (np.pi/4, np.pi/3, np.pi/6), (0, 0, 0)]
        
        for x, y, z in angles:
            R = compute_rotation_matrix(x, y, z)
            
            # Check rotation matrix properties
            assert is_rotation_matrix(R)
    
    def test_composite_rotation_decomposition(self):
        """Test composite rotation decomposition into individual rotations."""
        x, y, z = 0.1, 0.2, 0.3
        
        # Compute combined rotation matrix
        R_combined = compute_rotation_matrix(x, y, z)
        
        # Compute individual rotation matrices and combine
        R_x = rotation_matrix_x(x)
        R_y = rotation_matrix_y(y)
        R_z = rotation_matrix_z(z)
        
        # ZYX convention means R_z * R_y * R_x
        R_manual = R_z @ R_y @ R_x
        
        assert_array_almost_equal(R_combined, R_manual)


class TestRoundTripConversion:
    """Tests for round-trip conversion between Euler angles and rotation matrices."""
    
    @pytest.mark.parametrize("angles", [
        (0.1, 0.2, 0.3),
        (np.pi/4, np.pi/3, np.pi/6),
        (0, 0, 0),
        (-0.1, -0.2, -0.3),
    ])
    def test_angles_to_matrix_to_angles(self, angles):
        """Test round-trip conversion from angles to matrix and back."""
        x, y, z = angles
        
        # Convert angles to rotation matrix
        R = compute_rotation_matrix(x, y, z)
        
        # Convert rotation matrix back to angles
        recovered_angles = euler_angles_from_rotation_matrix(R)
        
        # Check that the recovered angles produce the same rotation matrix
        R_recovered = compute_rotation_matrix(*recovered_angles)
        
        # The recovered angles might not be exactly the same as the original ones
        # due to multiple representations, but they should produce the same matrix
        assert_array_almost_equal(R, R_recovered)


class TestEdgeCases:
    """Tests for edge cases, including gimbal lock."""
    
    def test_gimbal_lock_positive_y_90_degrees(self):
        """Test gimbal lock case where y is +90 degrees."""
        # Create a rotation with y = +90 degrees (gimbal lock)
        R = compute_rotation_matrix(0.1, np.pi/2, 0.3)
        
        # In gimbal lock, we lose one degree of freedom
        # Extract Euler angles
        recovered_angles = euler_angles_from_rotation_matrix(R)
        
        # Check that the recovered angles produce the same rotation matrix
        R_recovered = compute_rotation_matrix(*recovered_angles)
        
        assert_array_almost_equal(R, R_recovered)
    
    def test_gimbal_lock_negative_y_90_degrees(self):
        """Test gimbal lock case where y is -90 degrees."""
        # Create a rotation with y = -90 degrees (gimbal lock)
        R = compute_rotation_matrix(0.1, -np.pi/2, 0.3)
        
        # Extract Euler angles
        recovered_angles = euler_angles_from_rotation_matrix(R)
        
        # Check that the recovered angles produce the same rotation matrix
        R_recovered = compute_rotation_matrix(*recovered_angles)
        
        assert_array_almost_equal(R, R_recovered)
    
    def test_large_angles(self):
        """Test rotation with very large angles."""
        # Angles larger than 2π should work correctly
        x, y, z = 10*np.pi, -5*np.pi, 7*np.pi
        
        R = compute_rotation_matrix(x, y, z)
        
        # Check rotation matrix properties
        assert is_rotation_matrix(R)


class TestInputValidation:
    """Tests for input validation."""
    
    def test_validate_angle_with_valid_inputs(self):
        """Test _validate_angle function with valid inputs."""
        valid_inputs = [0, 1.0, np.pi, np.float32(1.5), np.float64(2.0)]
        
        for angle in valid_inputs:
            _validate_angle(angle, "test_angle")  # Should not raise an exception
    
    def test_validate_angle_with_invalid_inputs(self):
        """Test _validate_angle function with invalid inputs."""
        invalid_inputs = ["string", [1, 2, 3], {'a': 1}, None]
        
        for angle in invalid_inputs:
            with pytest.raises(TypeError):
                _validate_angle(angle, "test_angle")
    
    def test_rotation_matrix_x_with_invalid_input(self):
        """Test rotation_matrix_x function with invalid input."""
        with pytest.raises(TypeError):
            rotation_matrix_x("not a number")
    
    def test_rotation_matrix_y_with_invalid_input(self):
        """Test rotation_matrix_y function with invalid input."""
        with pytest.raises(TypeError):
            rotation_matrix_y("not a number")
    
    def test_rotation_matrix_z_with_invalid_input(self):
        """Test rotation_matrix_z function with invalid input."""
        with pytest.raises(TypeError):
            rotation_matrix_z("not a number")
    
    def test_compute_rotation_matrix_with_invalid_inputs(self):
        """Test compute_rotation_matrix function with invalid inputs."""
        with pytest.raises(TypeError):
            compute_rotation_matrix("not a number", 0, 0)
        
        with pytest.raises(TypeError):
            compute_rotation_matrix(0, "not a number", 0)
        
        with pytest.raises(TypeError):
            compute_rotation_matrix(0, 0, "not a number")
    
    def test_euler_angles_from_rotation_matrix_with_invalid_inputs(self):
        """Test euler_angles_from_rotation_matrix function with invalid inputs."""
        with pytest.raises(TypeError):
            euler_angles_from_rotation_matrix("not a matrix")
        
        with pytest.raises(ValueError):
            euler_angles_from_rotation_matrix(np.array([1, 2, 3]))  # Wrong shape
        
        with pytest.raises(ValueError):
            euler_angles_from_rotation_matrix(np.zeros((2, 2)))  # Wrong shape


class TestMatrixProperties:
    """Tests for rotation matrix properties."""
    
    @pytest.mark.parametrize("angles", [
        (0.1, 0.2, 0.3),
        (np.pi/4, np.pi/3, np.pi/6),
        (0, 0, 0),
        (-0.1, -0.2, -0.3),
    ])
    def test_determinant_is_one(self, angles):
        """Test that the determinant of rotation matrices is 1."""
        x, y, z = angles
        
        # Individual rotations
        assert abs(np.linalg.det(rotation_matrix_x(x)) - 1.0) < 1e-10
        assert abs(np.linalg.det(rotation_matrix_y(y)) - 1.0) < 1e-10
        assert abs(np.linalg.det(rotation_matrix_z(z)) - 1.0) < 1e-10
        
        # Composite rotation
        R = compute_rotation_matrix(x, y, z)
        assert abs(np.linalg.det(R) - 1.0) < 1e-10
    
    @pytest.mark.parametrize("angles", [
        (0.1, 0.2, 0.3),
        (np.pi/4, np.pi/3, np.pi/6),
        (0, 0, 0),
        (-0.1, -0.2, -0.3),
    ])
    def test_orthogonality(self, angles):
        """Test that rotation matrices are orthogonal (R.T @ R = I)."""
        x, y, z = angles
        
        # Individual rotations
        for R in [rotation_matrix_x(x), rotation_matrix_y(y), rotation_matrix_z(z)]:
            assert_array_almost_equal(R.T @ R, np.eye(3))
            assert_array_almost_equal(R @ R.T, np.eye(3))
        
        # Composite rotation
        R = compute_rotation_matrix(x, y, z)
        assert_array_almost_equal(R.T @ R, np.eye(3))
        assert_array_almost_equal(R @ R.T, np.eye(3))
    
    def test_rotational_invariance(self):
        """Test rotational invariance of vectors with same length (preserves length)."""
        vectors = [
            np.array([1, 0, 0]),
            np.array([0, 1, 0]),
            np.array([0, 0, 1]),
            np.array([1, 1, 1]) / np.sqrt(3),
        ]
        
        angles = (0.7, 0.8, 0.9)
        R = compute_rotation_matrix(*angles)
        
        for v in vectors:
            original_length = np.linalg.norm(v)
            rotated_v = R @ v
            rotated_length = np.linalg.norm(rotated_v)
            
            assert abs(original_length - rotated_length) < 1e-10
    
    def test_successive_rotations(self):
        """Test that successive rotations combine correctly."""
        # First rotation
        angles1 = (0.1, 0.2, 0.3)
        R1 = compute_rotation_matrix(*angles1)
        
        # Second rotation
        angles2 = (0.4, 0.5, 0.6)
        R2 = compute_rotation_matrix(*angles2)
        
        # Combined rotation by matrix multiplication
        R_combined = R2 @ R1
        
        # Applying combined rotation to a vector
        v = np.array([1, 2, 3])
        result1 = R_combined @ v
        
        # Applying rotations sequentially
        result2 = R2 @ (R1 @ v)
        
        assert_array_almost_equal(result1, result2)

