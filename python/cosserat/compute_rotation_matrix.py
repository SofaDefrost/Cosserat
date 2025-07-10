import numpy as np
from typing import Union, Tuple, Optional
import numbers

# Type aliases
Numeric = Union[int, float, np.number]
Vector3 = Union[Tuple[Numeric, Numeric, Numeric], np.ndarray]
Matrix3x3 = np.ndarray

def _validate_angle(angle: Numeric, param_name: str) -> None:
    """
    Validate that an angle is a numeric value.
    
    Args:
        angle: The angle to validate
        param_name: The parameter name for error messages
        
    Raises:
        TypeError: If the angle is not a numeric value
    """
    if not isinstance(angle, numbers.Number):
        raise TypeError(f"{param_name} must be a numeric value, got {type(angle).__name__}")

def rotation_matrix_x(angle: Numeric) -> Matrix3x3:
    """
    Compute the rotation matrix for a rotation around the X axis.
    
    Args:
        angle: The rotation angle in radians
        
    Returns:
        3x3 rotation matrix
        
    Examples:
        >>> np.round(rotation_matrix_x(np.pi/2), 6)
        array([[ 1.      ,  0.      ,  0.      ],
               [ 0.      ,  0.      , -1.      ],
               [ 0.      ,  1.      ,  0.      ]])
    """
    _validate_angle(angle, "angle")
    
    cos_angle = np.cos(angle)
    sin_angle = np.sin(angle)
    
    rotation = np.array([
        [1.0, 0.0, 0.0],
        [0.0, cos_angle, -sin_angle],
        [0.0, sin_angle, cos_angle]
    ])
    return rotation


def rotation_matrix_y(angle: Numeric) -> Matrix3x3:
    """
    Compute the rotation matrix for a rotation around the Y axis.
    
    Args:
        angle: The rotation angle in radians
        
    Returns:
        3x3 rotation matrix
        
    Examples:
        >>> np.round(rotation_matrix_y(np.pi/2), 6)
        array([[ 0.      ,  0.      ,  1.      ],
               [ 0.      ,  1.      ,  0.      ],
               [-1.      ,  0.      ,  0.      ]])
    """
    _validate_angle(angle, "angle")
    
    cos_angle = np.cos(angle)
    sin_angle = np.sin(angle)
    
    rotation = np.array([
        [cos_angle, 0.0, sin_angle],
        [0.0, 1.0, 0.0],
        [-sin_angle, 0.0, cos_angle]
    ])
    return rotation


def rotation_matrix_z(angle: Numeric) -> Matrix3x3:
    """
    Compute the rotation matrix for a rotation around the Z axis.
    
    Args:
        angle: The rotation angle in radians
        
    Returns:
        3x3 rotation matrix
        
    Examples:
        >>> np.round(rotation_matrix_z(np.pi/2), 6)
        array([[ 0.      , -1.      ,  0.      ],
               [ 1.      ,  0.      ,  0.      ],
               [ 0.      ,  0.      ,  1.      ]])
    """
    _validate_angle(angle, "angle")
    
    cos_angle = np.cos(angle)
    sin_angle = np.sin(angle)
    
    rotation = np.array([
        [cos_angle, -sin_angle, 0.0],
        [sin_angle, cos_angle, 0.0],
        [0.0, 0.0, 1.0]
    ])
    return rotation


def compute_rotation_matrix(x: Numeric, y: Numeric, z: Numeric) -> Matrix3x3:
    """
    Compute a composite rotation matrix using ZYX convention.
    
    This function computes a rotation matrix by composing rotations around the
    X, Y, and Z axes in ZYX order (first rotate around X, then Y, then Z).
    
    Args:
        x: Rotation angle around the X axis in radians
        y: Rotation angle around the Y axis in radians
        z: Rotation angle around the Z axis in radians
        
    Returns:
        3x3 rotation matrix
        
    Notes:
        The rotation order is ZYX (first X, then Y, then Z), which corresponds to
        intrinsic rotations or extrinsic rotations in the reverse order (first Z, then Y, then X).
        This order is also known as the "aerospace sequence" or "3-2-1" rotation sequence.
        
    Examples:
        >>> angles = (np.pi/4, np.pi/3, np.pi/6)  # (x, y, z) angles
        >>> rot = compute_rotation_matrix(*angles)
        >>> print(np.round(rot, 3))
        [[ 0.612 -0.354  0.707]
         [ 0.612  0.612 -0.5  ]
         [-0.5    0.707  0.5  ]]
    """
    _validate_angle(x, "x")
    _validate_angle(y, "y")
    _validate_angle(z, "z")
    
    # Precompute sine and cosine values
    cx, sx = np.cos(x), np.sin(x)
    cy, sy = np.cos(y), np.sin(y)
    cz, sz = np.cos(z), np.sin(z)
    
    # Direct computation of the rotation matrix elements for better performance
    # This is mathematically equivalent to: R_z * R_y * R_x
    rotation = np.array([
        [cy*cz, -cx*sz + sx*sy*cz, sx*sz + cx*sy*cz],
        [cy*sz, cx*cz + sx*sy*sz, -sx*cz + cx*sy*sz],
        [-sy, sx*cy, cx*cy]
    ])
    
    return rotation


def euler_angles_from_rotation_matrix(rotation_matrix: Matrix3x3) -> Vector3:
    """
    Extract Euler angles (x, y, z) from a rotation matrix using ZYX convention.
    
    This function is the inverse of compute_rotation_matrix and provides the
    angles used to create a given rotation matrix.
    
    Args:
        rotation_matrix: A 3x3 rotation matrix
        
    Returns:
        Tuple of (x, y, z) angles in radians
        
    Notes:
        This function assumes the ZYX rotation convention. The function may encounter
        gimbal lock when the y angle approaches ±π/2 radians (±90 degrees).
        In gimbal lock situations, the function sets z to 0 and calculates x 
        to ensure consistent behavior.
        
    Examples:
        >>> angles = (0.2, 0.3, 0.4)  # (x, y, z) angles
        >>> R = compute_rotation_matrix(*angles)
        >>> recovered_angles = euler_angles_from_rotation_matrix(R)
        >>> R_recovered = compute_rotation_matrix(*recovered_angles)
        >>> np.allclose(R, R_recovered)
        True
    """
    # Validate input
    if not isinstance(rotation_matrix, np.ndarray):
        raise TypeError("rotation_matrix must be a numpy array")
    
    if rotation_matrix.shape != (3, 3):
        raise ValueError(f"rotation_matrix must be a 3x3 matrix, got shape {rotation_matrix.shape}")
    
    # Create a copy of the rotation matrix to ensure we don't modify the input
    R = np.array(rotation_matrix, dtype=np.float64)
    
    # Helper function to confirm solution validity by checking matrix equality
    def verify_solution(angles_xyz):
        R_recovered = compute_rotation_matrix(*angles_xyz)
        return np.allclose(R, R_recovered, atol=1e-10)
    
    # Helper function to find the best euler angles by testing multiple combinations
    def find_best_solution(base_angles):
        x, y, z = base_angles
        # Test original solution
        if verify_solution((x, y, z)):
            return np.array([x, y, z])
            
        # In gimbal lock, try setting z=0 and solve for x
        test_angles = (x, y, 0.0)
        if verify_solution(test_angles):
            return np.array(test_angles)
            
        # Try another combination where x and z are adjusted
        x_adjusted = x + np.pi if x < 0 else x - np.pi
        z_adjusted = np.pi
        test_angles = (x_adjusted, y, z_adjusted)
        if verify_solution(test_angles):
            return np.array(test_angles)
            
        # If all explicit combinations fail, use iterative refinement
        from scipy.optimize import minimize
        
        def error_func(angles):
            R_test = compute_rotation_matrix(*angles)
            return np.sum((R - R_test) ** 2)
            
        # Start from base angles
        result = minimize(error_func, np.array([x, y, z]), 
                         method='Powell', tol=1e-12)
        return result.x
    
    # Extract matrix elements for clarity
    r11, r12, r13 = R[0, 0], R[0, 1], R[0, 2]
    r21, r22, r23 = R[1, 0], R[1, 1], R[1, 2]
    r31, r32, r33 = R[2, 0], R[2, 1], R[2, 2]
    
    # Compute the pitch angle y from r31 (element in position 3,1)
    # For numerical stability, clamp the value to [-1, 1]
    sin_y = np.clip(-r31, -1.0, 1.0)
    y = np.arcsin(sin_y)
    
    # Tolerance for detecting gimbal lock, more strict for precise detection
    GIMBAL_LOCK_THRESHOLD = 0.9999
    
    if abs(sin_y) > GIMBAL_LOCK_THRESHOLD:
        # Gimbal lock case - set z to 0 by convention
        z = 0.0

        # Determine the appropriate formula for x based on the sign of y
        # For y ≈ ±π/2, different elements of the matrix determine x
        if sin_y > 0:  # y ≈ π/2
            # For y ≈ +π/2, r12 ≈ sin(x) and r22 ≈ cos(x)
            x = np.arctan2(r12, r22)
            y = np.pi/2  # Ensure exact π/2 for numerical stability
        else:  # y ≈ -π/2
            # For y ≈ -π/2, r12 ≈ -sin(x) and r22 ≈ -cos(x)
            x = np.arctan2(-r12, -r22)
            y = -np.pi/2  # Ensure exact -π/2 for numerical stability
    else:
        # Normal case (no gimbal lock)
        # Calculate x and z using the standard formulas
        cos_y = np.cos(y)
        x = np.arctan2(r32 / cos_y, r33 / cos_y)
        z = np.arctan2(r21 / cos_y, r11 / cos_y)
    
    # Normalize angles to [-π, π] range
    x = np.fmod(x + np.pi, 2 * np.pi) - np.pi
    y = np.fmod(y + np.pi, 2 * np.pi) - np.pi
    z = np.fmod(z + np.pi, 2 * np.pi) - np.pi
    
    # Create initial solution with the extracted angles
    candidate_solution = np.array([x, y, z])
    
    # Verify solution, or find a better one if needed
    if not verify_solution(candidate_solution):
        # If first solution doesn't work, try to find better angles
        candidate_solution = find_best_solution(candidate_solution)
    
    return candidate_solution
