#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example demonstrating the usage of Lie group bindings in Cosserat plugin.

This example shows how to use the SO2, SO3, SE2, and SE3 Lie groups for various
operations including rotations, transformations, and conversions.
"""

import Sofa
import numpy as np
import math
from Cosserat.LieGroups import SO2, SO3, SE2, SE3

def print_section(title):
    """Helper to print section titles with formatting"""
    print("\n" + "="*80)
    print(f" {title} ".center(80, "-"))
    print("="*80)

def example_SO2():
    """Examples with SO2 (2D rotations)"""
    print_section("SO2 - 2D Rotations")
    
    # Create identity rotation
    identity = SO2()
    print("Identity angle:", identity.angle)
    
    # Create rotation from angle (in radians)
    r1 = SO2(math.pi/4)  # 45 degrees
    print("r1 angle (45°):", r1.angle)
    print("r1 matrix:\n", r1.matrix())
    
    # Compose rotations
    r2 = SO2(math.pi/6)  # 30 degrees
    r3 = r1 * r2
    print("r3 = r1 * r2 angle:", r3.angle)  # Should be 75 degrees (π/4 + π/6)
    print("r3 matrix:\n", r3.matrix())
    
    # Inverse rotation
    r_inv = r1.inverse()
    print("r1 inverse angle:", r_inv.angle)  # Should be -45 degrees (-π/4)
    
    # Verify identity property
    r_identity = r1 * r1.inverse()
    print("r1 * r1.inverse() is identity?", r_identity.isApprox(identity))
    
    # Rotate a point
    point = np.array([1.0, 0.0])
    rotated_point = r1.act(point)
    print("Original point:", point)
    print("Rotated point (45°):", rotated_point)
    
    # Logarithm (angle in radians)
    log_r1 = r1.log()
    print("log(r1):", log_r1)
    
    # Exponential map (create rotation from angle)
    exp_angle = SO2.exp(np.array([math.pi/3]))  # 60 degrees
    print("exp(π/3) angle:", exp_angle.angle)

def example_SO3():
    """Examples with SO3 (3D rotations)"""
    print_section("SO3 - 3D Rotations")
    
    # Create identity rotation
    identity = SO3()
    print("Identity quaternion:", identity.quaternion.coeffs())
    
    # Create rotation from angle-axis
    angle = math.pi/4  # 45 degrees
    axis = np.array([0.0, 0.0, 1.0])  # z-axis
    r1 = SO3(angle, axis)
    print("r1 matrix (45° around z-axis):\n", r1.matrix())
    
    # Create rotation from quaternion (using Eigen's format: w, x, y, z)
    from math import sin, cos
    half_angle = angle/2
    quat_coeffs = np.array([cos(half_angle), 0, 0, sin(half_angle)])
    r2 = SO3(quat_coeffs)
    print("r2 is approximately r1?", r2.isApprox(r1))
    
    # Compose rotations - rotate around z, then around y
    r_z = SO3(math.pi/4, np.array([0.0, 0.0, 1.0]))  # 45° around z
    r_y = SO3(math.pi/6, np.array([0.0, 1.0, 0.0]))  # 30° around y
    r3 = r_z * r_y
    print("Combined rotation matrix:\n", r3.matrix())
    
    # Inverse rotation
    r_inv = r_z.inverse()
    print("r_z inverse matrix:\n", r_inv.matrix())
    
    # Verify identity property
    r_identity = r_z * r_z.inverse()
    print("r_z * r_z.inverse() is identity?", r_identity.isApprox(identity))
    
    # Rotate a point
    point = np.array([1.0, 0.0, 0.0])
    rotated_point = r_z.act(point)
    print("Original point:", point)
    print("Rotated point (45° around z):", rotated_point)
    
    # Logarithm (rotation vector)
    log_r_z = r_z.log()
    print("log(r_z) (rotation vector):", log_r_z)
    
    # Exponential map (create rotation from rotation vector)
    rot_vector = np.array([0.0, math.pi/3, 0.0])  # 60° around y-axis
    exp_rot = SO3.exp(rot_vector)
    print("exp(rot_vector) matrix:\n", exp_rot.matrix())
    
    # Converting between representations
    angle_axis = r3.angleAxis()
    print("r3 as angle-axis: angle =", angle_axis.angle(), "axis =", angle_axis.axis())
    
    # Skew-symmetric matrix (hat operator)
    omega = np.array([1.0, 2.0, 3.0])
    omega_hat = SO3.hat(omega)
    print("hat([1,2,3]) (skew-symmetric matrix):\n", omega_hat)
    
    # Vee operator (inverse of hat)
    omega_recovered = SO3.vee(omega_hat)
    print("vee(hat([1,2,3])) =", omega_recovered)

def example_SE2():
    """Examples with SE2 (2D rigid transformations)"""
    print_section("SE2 - 2D Rigid Transformations")
    
    # Create identity transformation
    identity = SE2()
    print("Identity transformation - rotation:", identity.rotation.angle)
    print("Identity transformation - translation:", identity.translation)
    
    # Create from rotation and translation
    r = SO2(math.pi/4)  # 45 degrees
    t = np.array([1.0, 2.0])
    g1 = SE2(r, t)
    print("g1 rotation angle:", g1.rotation.angle)
    print("g1 translation:", g1.translation)
    print("g1 matrix:\n", g1.matrix())
    
    # Create from angle and translation
    g2 = SE2(math.pi/6, np.array([3.0, 1.0]))
    print("g2 matrix:\n", g2.matrix())
    
    # Compose transformations
    g3 = g1 * g2
    print("g3 = g1 * g2 matrix:\n", g3.matrix())
    
    # Inverse transformation
    g_inv = g1.inverse()
    print("g1 inverse matrix:\n", g_inv.matrix())
    
    # Verify identity property
    g_identity = g1 * g1.inverse()
    print("g1 * g1.inverse() is identity?", g_identity.isApprox(identity))
    
    # Transform a point
    point = np.array([1.0, 0.0])
    transformed_point = g1.act(point)
    print("Original point:", point)
    print("Transformed point:", transformed_point)
    
    # Logarithm (twist coordinates: vx, vy, omega)
    log_g1 = g1.log()
    print("log(g1) (twist coordinates):", log_g1)
    
    # Exponential map (create transformation from twist)
    twist = np.array([0.5, 1.0, math.pi/3])  # translation + 60° rotation
    exp_g = SE2.exp(twist)
    print("exp(twist) matrix:\n", exp_g.matrix())

def example_SE3():
    """Examples with SE3 (3D rigid transformations)"""
    print_section("SE3 - 3D Rigid Transformations")
    
    # Create identity transformation
    identity = SE3()
    print("Identity transformation - rotation quaternion:", identity.rotation.quaternion.coeffs())
    print("Identity transformation - translation:", identity.translation)
    
    # Create from rotation and translation
    r = SO3(math.pi/4, np.array([0.0, 0.0, 1.0]))  # 45° around z
    t = np.array([1.0, 2.0, 3.0])
    g1 = SE3(r, t)
    print("g1 rotation matrix:\n", g1.rotation.matrix())
    print("g1 translation:", g1.translation)
    print("g1 matrix:\n", g1.matrix())
    
    # Create from homogeneous transformation matrix
    import numpy as np
    T = np.eye(4)
    T[0:3, 0:3] = r.matrix()
    T[0:3, 3] = t
    g2 = SE3(T)
    print("g2 is approximately g1?", g2.isApprox(g1))
    
    # Compose transformations
    r_y = SO3(math.pi/6, np.array([0.0, 1.0, 0.0]))  # 30° around y
    t2 = np.array([0.0, 0.0, 1.0])
    g3 = SE3(r_y, t2)
    g4 = g1 * g3
    print("g4 = g1 * g3 matrix:\n", g4.matrix())
    
    # Inverse transformation
    g_inv = g1.inverse()
    print("g1 inverse matrix:\n", g_inv.matrix())
    
    # Verify identity property
    g_identity = g1 * g1.inverse()
    print("g1 * g1.inverse() is identity?", g_identity.isApprox(identity))
    
    # Transform a point
    point = np.array([1.0, 0.0, 0.0])
    transformed_point = g1.act(point)
    print("Original point:", point)
    print("Transformed point:", transformed_point)
    
    # Logarithm (twist coordinates: vx, vy, vz, wx, wy, wz)
    log_g1 = g1.log()
    print("log(g1) (twist coordinates):", log_g1)
    
    # Exponential map (create transformation from twist)
    twist = np.array([0.1, 0.2, 0.3, 0.0, 0.0, math.pi/3])  # translation + 60° around z
    exp_g = SE3.exp(twist)
    print("exp(twist) matrix:\n", exp_g.matrix())
    
    # Baker-Campbell-Hausdorff formula
    twist1 = np.array([0.1, 0.0, 0.0, 0.0, 0.0, math.pi/6])
    twist2 = np.array([0.0, 0.1, 0.0, 0.0, math.pi/6, 0.0])
    bch = SE3.BCH(twist1, twist2)
    print("BCH(twist1, twist2):", bch)
    
    # Compare with direct composition
    g_bch = SE3.exp(bch)
    g_direct = SE3.exp(twist1) * SE3.exp(twist2)
    print("exp(BCH) is approximately exp(t1)*exp(t2)?", g_bch.isApprox(g_direct))

def practical_example():
    """Practical example: robot kinematics"""
    print_section("Practical Example: Robot Kinematics")
    
    # Define a simple 2-link planar robot
    link1_length = 1.0
    link2_length = 0.8
    
    # Joint angles
    theta1 = math.pi/4  # 45 degrees
    theta2 = math.pi/3  # 60 degrees
    
    # Forward kinematics using SE2
    T1 = SE2(theta1, np.array([link1_length, 0.0]))  # First joint
    T2 = SE2(theta2, np.array([link2_length, 0.0]))  # Second joint
    
    # End effector position
    T_ee = T1 * T2
    print("End effector transformation:\n", T_ee.matrix())
    print("End effector position:", T_ee.translation)
    print("End effector orientation (degrees):", T_ee.rotation.angle * 180 / math.pi)
    
    # Demonstration with 3D robot (simplified)
    print("\n3D Robot Example:")
    
    # Create joint transformations
    j1 = SE3(SO3(theta1, np.array([0.0, 0.0, 1.0])), np.array([0.0, 0.0, 0.0]))
    j2 = SE3(SO3(0.0, np.array([0.0, 1.0, 0.0])), np.array([link1_length, 0.0, 0.0]))
    j3 = SE3(SO3(theta2, np.array([0.0, 1.0, 0.0])), np.array([link2_length, 0.0, 0.0]))
    
    # End effector transformation
    T3D_ee = j1 * j2 * j3
    print("3D End effector transformation:\n", T3D_ee.matrix())
    print("3D End effector position:", T3D_ee.translation)

def main():
    """Main function to run all examples"""
    print("Cosserat Lie Groups Examples")
    
    example_SO2()
    example_SO3()
    example_SE2()
    example_SE3()
    practical_example()
    
    print("\nAll examples completed successfully!")

if __name__ == "__main__":
    main()

