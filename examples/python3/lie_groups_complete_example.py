#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comprehensive example demonstrating the Lie group bindings in the Cosserat plugin.

This example covers:
1. Basic operations with SO2 and SO3 (rotations)
2. Rigid transformations with SE2 and SE3
3. Spacetime transformations with SGal3
4. Bundle usage for combined transformations
5. Utility functions for angles and vectors
6. Practical examples like robot kinematics and trajectory interpolation
"""

import Sofa
import numpy as np
import math
from Cosserat.LieGroups import SO2, SO3, SE2, SE3, SGal3, PoseVel, RotTrans
import Cosserat.LieGroups as lg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def print_section(title):
    """Helper to print section titles with formatting"""
    print("\n" + "="*80)
    print(f" {title} ".center(80, "-"))
    print("="*80)

def example_SO2():
    """Demonstration of SO2 (2D rotations)"""
    print_section("SO2 - 2D Rotations")
    
    # Create rotations
    identity = SO2()
    r1 = SO2(math.pi/4)  # 45 degrees
    r2 = SO2(math.pi/3)  # 60 degrees
    
    print("r1 angle:", r1.angle, "radians =", r1.angle * 180/math.pi, "degrees")
    print("r2 angle:", r2.angle, "radians =", r2.angle * 180/math.pi, "degrees")
    
    # Composition
    r3 = r1 * r2
    print("r3 = r1 * r2 angle:", r3.angle, "radians =", r3.angle * 180/math.pi, "degrees")
    
    # Inverse
    r1_inv = r1.inverse()
    print("r1 inverse angle:", r1_inv.angle, "radians =", r1_inv.angle * 180/math.pi, "degrees")
    
    # Group properties
    r_identity = r1 * r1_inv
    print("r1 * r1^(-1) ≈ identity?", r_identity.isApprox(identity))
    
    # Rotate points
    points = np.array([[1, 0], [0, 1], [1, 1], [-1, 0]])
    rotated_points = np.array([r1.act(p) for p in points])
    
    # Plot original and rotated points
    plt.figure(figsize=(8, 8))
    plt.scatter(points[:, 0], points[:, 1], label='Original Points')
    plt.scatter(rotated_points[:, 0], rotated_points[:, 1], label='Rotated Points')
    plt.grid(True)
    plt.axis('equal')
    plt.legend()
    plt.title('SO2 Rotation (45°)')
    plt.savefig('so2_rotation.png')
    
    # Lie algebra and exponential map
    v = np.array([math.pi/6])  # 30 degrees in Lie algebra
    r_exp = SO2.exp(v)
    print("exp(π/6) angle:", r_exp.angle, "radians =", r_exp.angle * 180/math.pi, "degrees")
    
    # Logarithm
    log_r1 = r1.log()
    print("log(r1):", log_r1, "should be ≈", r1.angle)
    
    # Interpolation between rotations using utility functions
    angles = [lg.slerp_angle(0, math.pi/2, t) for t in np.linspace(0, 1, 5)]
    print("Interpolation from 0 to 90°:", [a * 180/math.pi for a in angles])

def example_SO3():
    """Demonstration of SO3 (3D rotations)"""
    print_section("SO3 - 3D Rotations")
    
    # Create rotations
    identity = SO3()
    
    # Rotation around Z-axis by 45 degrees
    r1 = SO3(math.pi/4, np.array([0, 0, 1]))
    print("r1 (45° around Z) matrix:\n", r1.matrix())
    
    # Rotation around Y-axis by 30 degrees
    r2 = SO3(math.pi/6, np.array([0, 1, 0]))
    print("r2 (30° around Y) matrix:\n", r2.matrix())
    
    # Composition of rotations
    r3 = r1 * r2
    print("r3 = r1 * r2 matrix:\n", r3.matrix())
    
    # Convert to angle-axis representation
    aa = r3.angleAxis()
    print("r3 as angle-axis: angle =", aa.angle(), "radians, axis =", aa.axis())
    
    # Inverse rotation
    r1_inv = r1.inverse()
    print("r1 inverse matrix:\n", r1_inv.matrix())
    
    # Group properties
    r_identity = r1 * r1_inv
    print("r1 * r1^(-1) ≈ identity?", r_identity.isApprox(identity))
    
    # Rotate a 3D point
    point = np.array([1, 0, 0])
    rotated_point = r1.act(point)
    print(f"Rotating point {point} with r1 gives {rotated_point}")
    
    # Generate a set of points on a circle in the XY plane
    theta = np.linspace(0, 2*np.pi, 20)
    circle_points = np.column_stack((np.cos(theta), np.sin(theta), np.zeros_like(theta)))
    
    # Rotate all points
    rotated_circle = np.array([r1.act(p) for p in circle_points])
    
    # Plot original and rotated points in 3D
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(circle_points[:, 0], circle_points[:, 1], circle_points[:, 2], 
               label='Original Circle', color='blue')
    ax.scatter(rotated_circle[:, 0], rotated_circle[:, 1], rotated_circle[:, 2], 
               label='Rotated Circle', color='red')
    
    # Add coordinate axes
    ax.quiver(0, 0, 0, 1, 0, 0, color='r', label='X')
    ax.quiver(0, 0, 0, 0, 1, 0, color='g', label='Y')
    ax.quiver(0, 0, 0, 0, 0, 1, color='b', label='Z')
    
    ax.set_xlim([-1.5, 1.5])
    ax.set_ylim([-1.5, 1.5])
    ax.set_zlim([-1.5, 1.5])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    ax.set_title('SO3 Rotation (45° around Z-axis)')
    plt.savefig('so3_rotation.png')
    
    # Lie algebra operations
    # Create a rotation vector (element of the Lie algebra so(3))
    omega = np.array([0.1, 0.2, 0.3])
    
    # Convert to skew-symmetric matrix
    omega_hat = SO3.hat(omega)
    print("Skew-symmetric matrix (hat operator):\n", omega_hat)
    
    # Convert back to vector (vee operator)
    omega_vee = SO3.vee(omega_hat)
    print("Vector recovered (vee operator):", omega_vee)
    
    # Exponential map
    r_exp = SO3.exp(omega)
    print("exp(omega) quaternion:", r_exp.quaternion.coeffs())
    
    # Logarithm
    log_r1 = r1.log()
    print("log(r1):", log_r1, "- this is the rotation vector corresponding to r1")

def example_SE2():
    """Demonstration of SE2 (2D rigid transformations)"""
    print_section("SE2 - 2D Rigid Transformations")
    
    # Create transformations
    identity = SE2()
    
    # Create from rotation and translation
    r = SO2(math.pi/4)  # 45 degrees
    t = np.array([2.0, 1.0])
    g1 = SE2(r, t)
    
    print("g1 rotation angle:", g1.rotation.angle, "radians =", g1.rotation.angle * 180/math.pi, "degrees")
    print("g1 translation:", g1.translation)
    print("g1 matrix:\n", g1.matrix())
    
    # Create directly from angle and translation
    g2 = SE2(math.pi/6, np.array([1.0, 3.0]))  # 30 degrees
    print("g2 matrix:\n", g2.matrix())
    
    # Composition of transformations
    g3 = g1 * g2
    print("g3 = g1 * g2 matrix:\n", g3.matrix())
    print("g3 rotation angle:", g3.rotation.angle, "radians =", g3.rotation.angle * 180/math.pi, "degrees")
    print("g3 translation:", g3.translation)
    
    # Inverse transformation
    g1_inv = g1.inverse()
    print("g1 inverse matrix:\n", g1_inv.matrix())
    
    # Group properties
    g_identity = g1 * g1_inv
    print("g1 * g1^(-1) ≈ identity?", g_identity.isApprox(identity))
    
    # Transform points
    points = np.array([[0, 0], [1, 0], [0, 1], [1, 1]])
    transformed_points = np.array([g1.act(p) for p in points])
    
    # Plot original and transformed points
    plt.figure(figsize=(8, 8))
    plt.scatter(points[:, 0], points[:, 1], label='Original Points')
    plt.scatter(transformed_points[:, 0], transformed_points[:, 1], label='Transformed Points')
    
    # Draw coordinate frames
    origin = np.array([0, 0])
    x_axis = np.array([1, 0])
    y_axis = np.array([0, 1])
    
    # Original frame
    plt.arrow(origin[0], origin[1], x_axis[0], x_axis[1], color='red', width=0.02, head_width=0.1)
    plt.arrow(origin[0], origin[1], y_axis[0], y_axis[1], color='green', width=0.02, head_width=0.1)
    
    # Transformed frame
    new_origin = g1.translation
    new_x = g1.act(x_axis) - new_origin
    new_y = g1.act(y_axis) - new_origin
    plt.arrow(new_origin[0], new_origin[1], new_x[0], new_x[1], color='darkred', width=0.02, head_width=0.1)
    plt.arrow(new_origin[0], new_origin[1], new_y[0], new_y[1], color='darkgreen', width=0.02, head_width=0.1)
    
    plt.grid(True)
    plt.axis('equal')
    plt.xlim(-1, 4)
    plt.ylim(-1, 4)
    plt.legend()
    plt.title('SE2 Rigid Transformation')
    plt.savefig('se2_transform.png')
    
    # Lie algebra and exponential map
    # SE(2) Lie algebra: [vx, vy, omega]
    twist = np.array([0.5, 1.0, math.pi/4])
    g_exp = SE2.exp(twist)
    print("exp(twist) matrix:\n", g_exp.matrix())
    
    # Logarithm
    log_g1 = g1.log()
    print("log(g1):", log_g1, "- this is the twist coordinates for g1")

def example_SE3():
    """Demonstration of SE3 (3D rigid transformations)"""
    print_section("SE3 - 3D Rigid Transformations")
    
    # Create transformations
    identity = SE3()
    
    # Create from rotation and translation
    r = SO3(math.pi/4, np.array([0, 0, 1]))  # 45 degrees around Z
    t = np.array([2.0, 1.0, 0.5])
    g1 = SE3(r, t)
    
    print("g1 rotation matrix:\n", g1.rotation.matrix())
    print("g1 translation:", g1.translation)
    print("g1 homogeneous matrix:\n", g1.matrix())
    
    # Create from homogeneous transformation matrix
    T = np.eye(4)
    T[0:3, 0:3] = r.matrix()
    T[0:3, 3] = t
    g2 = SE3(T)
    print("g2 is approximately g1?", g2.isApprox(g1))
    
    # Composition of transformations
    r2 = SO3(math.pi/6, np.array([1, 0, 0]))  # 30 degrees around X
    t2 = np.array([0.0, 0.0, 1.0])
    g3 = SE3(r2, t2)
    g4 = g1 * g3
    print("g4 = g1 * g3 matrix:\n", g4.matrix())
    
    # Inverse transformation
    g1_inv = g1.inverse()
    print("g1 inverse matrix:\n", g1_inv.matrix())
    
    # Group properties
    g_identity = g1 * g1_inv
    print("g1 * g1^(-1) ≈ identity?", g_identity.isApprox(identity))
    
    # Transform a 3D point
    point = np.array([1, 0, 0])
    transformed_point = g1.act(point)
    print(f"Transforming point {point} with g1 gives {transformed_point}")
    
    # Generate a set of points in 3D
    x = np.linspace(-1, 1, 5)
    y = np.linspace(-1, 1, 5)
    z = np.zeros((5, 5))
    X, Y = np.meshgrid(x, y)
    points = np.column_stack((X.flatten(), Y.flatten(), z.flatten()))
    
    # Transform all points
    transformed_points = np.array([g1.act(p) for p in points])
    
    # Plot original and transformed points in 3D
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], 
               label='Original Points', color='blue')
    ax.scatter(transformed_points[:, 0], transformed_points[:, 1], transformed_points[:, 2], 
               label='Transformed Points', color='red')
    
    # Draw coordinate frames
    origin = np.zeros(3)
    x_axis = np.array([1, 0, 0])
    y_axis = np.array([0, 1, 0])
    z_axis = np.array([0, 0, 1])
    
    # Original frame
    ax.quiver(origin[0], origin[1], origin[2], 
              x_axis[0], x_axis[1], x_axis[2], color='r', label='X')
    ax.quiver(origin[0], origin[1], origin[2], 
              y_axis[0], y_axis[1], y_axis[2], color='g', label='Y')
    ax.quiver(origin[0], origin[1], origin[2], 
              z_axis[0], z_axis[1], z_axis[2], color='b', label='Z')
    
    # Transformed frame
    new_origin = g1.translation
    new_x = g1.rotation.act(x_axis)
    new_y = g1.rotation.act(y_axis)
    new_z = g1.rotation.act(z_axis)
    
    ax.quiver(new_origin[0], new_origin[1], new_origin[2], 
              new_x[0], new_x[1], new_x[2], color='darkred')
    ax.quiver(new_origin[0], new_origin[1], new_origin[2], 
              new_y[0], new_y[1], new_y[2], color='darkgreen')
    ax.quiver(new_origin[0], new_origin[1], new_origin[2], 
              new_z[0], new_z[1], new_z[2], color='darkblue')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim([-2, 3])
    ax.set_ylim([-2, 3])
    ax.set_zlim([-1, 2])
    ax.legend()
    ax.set_title('SE3 Rigid Transformation')
    plt.savefig('se3_transform.png')
    
    # Lie algebra and exponential map
    # SE(3) Lie algebra: [vx, vy, vz, wx, wy, wz]
    twist = np.array([0.1, 0.2, 0.3, 0.0, 0.0, math.pi/3])
    g_exp = SE3.exp(twist)
    print("exp(twist) matrix:\n", g_exp.matrix())
    
    # Logarithm
    log_g1 = g1.log()
    print("log(g1):", log_g1, "- this is the twist coordinates for g1")
    
    # Baker-Campbell-Hausdorff formula demonstration
    twist1 = np.array([0.1, 0.0, 0.0, 0.0, 0.0, math.pi/6])
    twist2 = np.array([0.0, 0.1, 0.0, 0.0, math.pi/6, 0.0])
    
    # Compose transformations
    g_t1 = SE3.exp(twist1)
    g_t2 = SE3.exp(twist2)
    g_composed = g_t1 * g_t2
    
    # Calculate BCH approximation
    bch = SE3.BCH(twist1, twist2)
    g_bch = SE3.exp(bch)
    
    print("Direct composition vs BCH approximation:")
    print("g_composed.matrix():\n", g_composed.matrix())
    print("g_bch.matrix():\n", g_bch.matrix())
    print("Difference:", np.linalg.norm(g_composed.matrix() - g_bch.matrix()))

def example_SGal3():
    """Demonstration of SGal3 (Special Galilean group - spacetime transformations)"""
    print_section("SGal3 - Special Galilean Group")
    
    # Create transformations
    identity = SGal3()
    
    # Create from SE3 pose, velocity, and time
    r = SO3(math.pi/4, np.array([0, 0, 1]))  # 45 degrees around Z
    t = np.array([1.0, 2.0, 0.0])
    pose = SE3(r, t)
    velocity = np.array([0.5, 0.2, 0.1])  # Linear velocity
    time = 0.0
    
    g1 = SGal3(pose, velocity, time)
    
    print("g1 pose:\n", g1.pose.matrix())
    print("g1 velocity:", g1.velocity)
    print("g1 time:", g1.time)
    
    # Create using factory functions
    g2 = lg.fromComponents(t, r, velocity, time=1.0)
    print("g2 created with fromComponents:")
    print("- pose:\n", g2.pose.matrix())
    print("- velocity:", g2.velocity)
    print("- time:", g2.time)
    
    # Create using Euler angles
    g3 = lg.fromPositionEulerVelocityTime(
        position=np.array([2.0, 0.0, 1.0]),
        roll=0.0,
        pitch=math.pi/6,  # 30 degrees
        yaw=math.pi/4,    # 45 degrees
        velocity=np.array([0.1, 0.2, 0.3]),
        time=2.0
    )
    print("g3 created with fromPositionEulerVelocityTime:")
    print("- pose:\n", g3.pose.matrix())
    print("- velocity:", g3.velocity)
    print("- time:", g3.time)
    
    # Convert to position, Euler angles, velocity, and time
    params = lg.toPositionEulerVelocityTime(g3)
    print("g3 converted to parameters:")
    print("- position:", params[0:3])
    print("- euler angles (rad):", params[3:6])
    print("- euler angles (deg):", params[3:6] * 180/math.pi)
    print("- velocity:", params[6:9])
    print("- time:", params[9])
    
    # Composition of transformations
    g4 = g1 * g2
    print("\nComposed transformation g4 = g1 * g2:")
    print("- pose:\n", g4.pose.matrix())
    print("- velocity:", g4.velocity)
    print("- time:", g4.time)
    
    # Inverse transformation
    g1_inv = g1.inverse()
    print("\ng1 inverse:")
    print("- pose:\n", g1_inv.pose.matrix())
    print("- velocity:", g1_inv.velocity)
    print("- time:", g1_inv.time)
    
    # Group properties
    g_identity = g1 * g1_inv
    print("\ng1 * g1^(-1) ≈ identity?", g_identity.isApprox(identity))
    
    # Transform a point-velocity-time tuple
    point_vel_time = np.array([1.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0])
    transformed = g1.transform(point_vel_time)
    print("\nTransforming point-velocity-time tuple:")
    print("- original:", point_vel_time)
    print("- transformed:", transformed)
    
    # Inverse transformation using division operator
    back_transformed = transformed / g1
    print("- back-transformed:", back_transformed)
    print("- original ≈ back-transformed?", np.allclose(point_vel_time, back_transformed))
    
    # Interpolation between two Galilean transformations
    print("\nInterpolation between g1 and g3:")
    
    for t in np.linspace(0, 1, 5):
        g_interp = lg.interpolate(g1, g3, t)
        print(f"Interpolation at t={t}:")
        print(f"- position: {g_interp.pose.translation}")
        print(f"- velocity: {g_interp.velocity}")
        print(f"- time: {g_interp.time}")

def example_Bundle():
    """Demonstration of Bundle (product manifolds)"""
    print_section("Bundle - Product Manifolds")
    
    # PoseVel - Bundle of SE3 pose and R3 velocity
    identity_pose_vel = PoseVel()
    
    # Create a pose and velocity
    r = SO3(math.pi/4, np.array([0, 1, 0]))  # 45 degrees around Y
    t = np.array([1.0, 0.0, 2.0])
    pose = SE3(r, t)
    velocity = np.array([0.5, 0.1, 0.2])
    
    # Create bundle
    pose_vel = PoseVel(pose, velocity)
    
    # Access components
    print("PoseVel bundle:")
    print("- pose matrix:\n", pose_vel.get_pose().matrix())
    print("- velocity:", pose_vel.get_velocity())
    
    # Group operations
    pose_vel2 = PoseVel(
        SE3(SO3(math.pi/6, np.array([1, 0, 0])), np.array([0.0, 1.0, 0.0])),
        np.array([0.1, 0.3, 0.0])
    )
    
    # Composition
    composed = pose_vel * pose_vel2
    print("\nComposed PoseVel:")
    print("- pose matrix:\n", composed.get_pose().matrix())
    print("- velocity:", composed.get_velocity())
    
    # Inverse
    inverse = pose_vel.inverse()
    print("\nInverse PoseVel:")
    print("- pose matrix:\n", inverse.get_pose().matrix())
    print("- velocity:", inverse.get_velocity())
    
    # RotTrans - Bundle of SO3 rotation and R3 translation
    identity_rot_trans = RotTrans()
    
    # Create a rotation and translation
    rot = SO3(math.pi/3, np.array([0, 0, 1]))  # 60 degrees around Z
    trans = np.array([2.0, 1.0, 0.5])
    
    # Create bundle
    rot_trans = RotTrans(rot, trans)
    
    # Access components
    print("\nRotTrans bundle:")
    print("- rotation matrix:\n", rot_trans.get_rotation().matrix())
    print("- translation:", rot_trans.get_translation())
    
    # Action on a point (using the product structure)
    point = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # [position, velocity]
    transformed = pose_vel.act(point)
    print("\nAction of PoseVel on a point-velocity:")
    print("- original:", point)
    print("- transformed:", transformed)

def example_utils():
    """Demonstration of utility functions"""
    print_section("Utility Functions")
    
    # Angle utilities
    angle1 = 3.5  # > π
    angle2 = -4.0  # < -π
    
    print("Original angles:", angle1, angle2)
    print("Normalized angles:", lg.normalize_angle(angle1), lg.normalize_angle(angle2))
    
    angle_a = math.pi/4  # 45 degrees
    angle_b = -math.pi/6  # -30 degrees
    
    print(f"Difference between {angle_a} and {angle_b}:", lg.angle_difference(angle_a, angle_b))
    print(f"Distance between {angle_a} and {angle_b}:", lg.angle_distance(angle_a, angle_b))
    
    small_angle = 1e-12
    print(f"Is {small_angle} near zero?", lg.is_angle_near_zero(small_angle))
    
    angle_c = angle_a + 1e-12
    print(f"Are {angle_a} and {angle_c} nearly equal?", lg.are_angles_nearly_equal(angle_a, angle_c))
    
    # Interpolation utilities
    scalar1 = 10
    scalar2 = 20
    print(f"Linear interpolation between {scalar1} and {scalar2}:")
    for t in [0, 0.25, 0.5, 0.75, 1.0]:
        print(f"  t={t}: {lg.lerp(scalar1, scalar2, t)}")
    
    angle_start = 0
    angle_end = math.pi  # 180 degrees
    print(f"Spherical linear interpolation from {angle_start} to {angle_end}:")
    for t in [0, 0.25, 0.5, 0.75, 1.0]:
        interp = lg.slerp_angle(angle_start, angle_end, t)
        print(f"  t={t}: {interp} rad = {interp * 180/math.pi} deg")
    
    # Numerical utilities
    x_values = [0.001, 0.01, 0.1, 1.0]
    print("\nNumerical utilities for small angles:")
    for x in x_values:
        print(f"  x={x}:")
        print(f"    sinc(x): {lg.sinc(x)} vs {np.sin(x)/x if x != 0 else 1}")
        print(f"    1-cos(x): {lg.one_minus_cos(x)} vs {1-np.cos(x)}")
    
    # Vector utilities
    v1 = np.array([1.0, 2.0, 3.0])
    v2 = np.array([0.0, 0.0, 1.0])
    
    print("\nVector utilities:")
    print(f"  v1 = {v1}, v2 = {v2}")
    print(f"  Normalized v1: {lg.safe_normalize(v1)}")
    print(f"  Projection of v1 onto v2: {lg.project_vector(v1, v2)}")
    
    # SE(2) path interpolation
    start_config = np.array([0, 0, 0])  # [angle, x, y]
    end_config = np.array([math.pi/2, 1, 1])  # [angle, x, y]
    
    print("\nSE(2) path interpolation from [0, 0, 0] to [π/2, 1, 1]:")
    for t in [0, 0.25, 0.5, 0.75, 1.0]:
        interp = lg.interpolate_se2_path(start_config, end_config, t)
        print(f"  t={t}: angle={interp[0]} rad = {interp[0] * 180/math.pi} deg, position=({interp[1]}, {interp[2]})")

def robot_kinematics_example():
    """Practical example: robot kinematics"""
    print_section("Practical Example: Robot Kinematics")
    
    # Define a 3-link planar robot
    link1_length = 1.0
    link2_length = 0.8
    link3_length = 0.6
    
    # Joint angles (in radians)
    theta1 = math.pi/4   # 45 degrees
    theta2 = math.pi/6   # 30 degrees
    theta3 = -math.pi/3  # -60 degrees
    
    # Forward kinematics using SE2
    print("2D Robot Kinematics with SE2:")
    
    # Base frame
    T0 = SE2()  # Identity transformation
    
    # Transformations for each joint
    T1 = SE2(theta1, np.array([link1_length, 0.0]))  # First joint
    T2 = SE2(theta2, np.array([link2_length, 0.0]))  # Second joint
    T3 = SE2(theta3, np.array([link3_length, 0.0]))  # Third joint
    
    # Compute end-effector pose
    T_ee = T0 * T1 * T2 * T3
    
    print("End-effector pose:")
    print("- matrix:\n", T_ee.matrix())
    print("- position:", T_ee.translation)
    print("- orientation:", T_ee.rotation.angle * 180/math.pi, "degrees")
    
    # Visualize the robot
    plt.figure(figsize=(10, 8))
    
    # Plot the links
    p0 = np.array([0, 0])
    p1 = T1.translation
    p2 = (T1 * T2).translation
    p3 = T_ee.translation
    
    plt.plot([p0[0], p1[0]], [p0[1], p1[1]], 'b-', linewidth=3, label='Link 1')
    plt.plot([p1[0], p2[0]], [p1[1], p2[1]], 'g-', linewidth=3, label='Link 2')
    plt.plot([p2[0], p3[0]], [p2[1], p3[1]], 'r-', linewidth=3, label='Link 3')
    
    # Plot the joints
    plt.scatter([p0[0], p1[0], p2[0]], [p0[1], p1[1], p2[1]], 
                color='black', s=80, label='Joints')
    plt.scatter([p3[0]], [p3[1]], color='orange', s=100, label='End-effector')
    
    plt.grid(True)
    plt.axis('equal')
    plt.xlim(-0.5, 2.5)
    plt.ylim(-0.5, 2.5)
    plt.legend()
    plt.title('2D Robot Kinematics with SE2')
    plt.savefig('robot_kinematics_2d.png')
    
    # 3D Robot Kinematics with SE3
    print("\n3D Robot Kinematics with SE3:")
    
    # Base frame
    T0_3d = SE3()  # Identity transformation
    
    # Transformations for each joint (now in 3D)
    # First joint: rotate around Z, then extend along X
    T1_3d = SE3(
        SO3(theta1, np.array([0, 0, 1])),
        np.array([0, 0, 0])
    ) * SE3(
        SO3(),  # Identity rotation
        np.array([link1_length, 0, 0])
    )
    
    # Second joint: rotate around Y, then extend along X
    T2_3d = SE3(
        SO3(theta2, np.array([0, 1, 0])),
        np.array([0, 0, 0])
    ) * SE3(
        SO3(),  # Identity rotation
        np.array([link2_length, 0, 0])
    )
    
    # Third joint: rotate around Y, then extend along X
    T3_3d = SE3(
        SO3(theta3, np.array([0, 1, 0])),
        np.array([0, 0, 0])
    ) * SE3(
        SO3(),  # Identity rotation
        np.array([link3_length, 0, 0])
    )
    
    # Compute transformations to each joint
    T_to_1 = T0_3d * T1_3d
    T_to_2 = T_to_1 * T2_3d
    T_to_3 = T_to_2 * T3_3d  # End-effector
    
    print("End-effector pose:")
    print("- matrix:\n", T_to_3.matrix())
    print("- position:", T_to_3.translation)
    
    # Get positions for visualization
    p0_3d = np.array([0, 0, 0])
    p1_3d = T_to_1.translation
    p2_3d = T_to_2.translation
    p3_3d = T_to_3.translation
    
    # Visualize the 3D robot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the links
    ax.plot([p0_3d[0], p1_3d[0]], [p0_3d[1], p1_3d[1]], [p0_3d[2], p1_3d[2]], 
             'b-', linewidth=3, label='Link 1')
    ax.plot([p1_3d[0], p2_3d[0]], [p1_3d[1], p2_3d[1]], [p1_3d[2], p2_3d[2]], 
             'g-', linewidth=3, label='Link 2')
    ax.plot([p2_3d[0], p3_3d[0]], [p2_3d[1], p3_3d[1]], [p2_3d[2], p3_3d[2]], 
             'r-', linewidth=3, label='Link 3')
    
    # Plot the joints
    ax.scatter([p0_3d[0], p1_3d[0], p2_3d[0]], 
               [p0_3d[1], p1_3d[1], p2_3d[1]], 
               [p0_3d[2], p1_3d[2], p2_3d[2]], 
               color='black', s=80, label='Joints')
    ax.scatter([p3_3d[0]], [p3_3d[1]], [p3_3d[2]], 
               color='orange', s=100, label='End-effector')
    
    # Draw coordinate frames at each joint
    for i, T in enumerate([T0_3d, T_to_1, T_to_2, T_to_3]):
        origin = T.translation
        x_dir = T.rotation.act(np.array([0.2, 0, 0]))
        y_dir = T.rotation.act(np.array([0, 0.2, 0]))
        z_dir = T.rotation.act(np.array([0, 0, 0.2]))
        
        ax.quiver(origin[0], origin[1], origin[2], 
                  x_dir[0], x_dir[1], x_dir[2], color='r')
        ax.quiver(origin[0], origin[1], origin[2], 
                  y_dir[0], y_dir[1], y_dir[2], color='g')
        ax.quiver(origin[0], origin[1], origin[2], 
                  z_dir[0], z_dir[1], z_dir[2], color='b')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim([0, 3])
    ax.set_ylim([-1.5, 1.5])
    ax.set_zlim([-1.5, 1.5])
    ax.legend()
    ax.set_title('3D Robot Kinematics with SE3')
    plt.savefig('robot_kinematics_3d.png')

def trajectory_interpolation_example():
    """Practical example: trajectory interpolation"""
    print_section("Practical Example: Trajectory Interpolation")
    
    # Define waypoints for a 3D trajectory using SE3
    waypoints = [
        SE3(SO3(), np.array([0, 0, 0])),  # Starting point (identity)
        SE3(SO3(math.pi/4, np.array([0, 0, 1])), np.array([1, 1, 0])),  # Waypoint 1
        SE3(SO3(math.pi/2, np.array([0, 1, 0])), np.array([2, 0, 1])),  # Waypoint 2
        SE3(SO3(math.pi/2, np.array([1, 0, 0])), np.array([3, 1, 0]))   # End point
    ]
    
    # Generate a smooth trajectory
    num_points = 100
    trajectory = []
    
    # Interpolate between waypoints
    segments = len(waypoints) - 1
    points_per_segment = num_points // segments
    
    for i in range(segments):
        start = waypoints[i]
        end = waypoints[i+1]
        
        for t in np.linspace(0, 1, points_per_segment):
            # Interpolate between waypoints
            # For a real implementation, you might want a more sophisticated spline interpolation
            start_log = start.log()
            end_log = end.log()
            interp_log = start_log + t * (end_log - start_log)
            interp_pose = SE3.exp(interp_log)
            trajectory.append(interp_pose)
    
    # Extract positions for visualization
    positions = np.array([pose.translation for pose in trajectory])
    
    # Visualize the trajectory
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the trajectory
    ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], 'b-', linewidth=2, label='Trajectory')
    
    # Plot waypoints
    waypoint_positions = np.array([pose.translation for pose in waypoints])
    ax.scatter(waypoint_positions[:, 0], waypoint_positions[:, 1], waypoint_positions[:, 2], 
               color='red', s=100, label='Waypoints')
    
    # Draw coordinate frames at each waypoint
    for i, pose in enumerate(waypoints):
        origin = pose.translation
        x_dir = pose.rotation.act(np.array([0.3, 0, 0]))
        y_dir = pose.rotation.act(np.array([0, 0.3, 0]))
        z_dir = pose.rotation.act(np.array([0, 0, 0.3]))
        
        ax.quiver(origin[0], origin[1], origin[2], 
                  x_dir[0], x_dir[1], x_dir[2], color='r')
        ax.quiver(origin[0], origin[1], origin[2], 
                  y_dir[0], y_dir[1], y_dir[2], color='g')
        ax.quiver(origin[0], origin[1], origin[2], 
                  z_dir[0], z_dir[1], z_dir[2], color='b')
    
    # Also draw some frames along the trajectory
    step = len(trajectory) // 10
    for i in range(0, len(trajectory), step):
        pose = trajectory[i]
        origin = pose.translation
        x_dir = pose.rotation.act(np.array([0.2, 0, 0]))
        y_dir = pose.rotation.act(np.array([0, 0.2, 0]))
        z_dir = pose.rotation.act(np.array([0, 0, 0.2]))
        
        ax.quiver(origin[0], origin[1], origin[2], 
                  x_dir[0], x_dir[1], x_dir[2], color='darkred', alpha=0.5)
        ax.quiver(origin[0], origin[1], origin[2], 
                  y_dir[0], y_dir[1], y_dir[2], color='darkgreen', alpha=0.5)
        ax.quiver(origin[0], origin[1], origin[2], 
                  z_dir[0], z_dir[1], z_dir[2], color='darkblue', alpha=0.5)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    ax.set_title('SE3 Trajectory Interpolation')
    plt.savefig('trajectory_interpolation.png')
    
    print("Generated a smooth trajectory with", len(trajectory), "poses")
    print("Trajectory starts at", trajectory[0].translation, "and ends at", trajectory[-1].translation)
    
    # Add velocity to create a SGal3 trajectory (spacetime)
    print("\nExtending to a spacetime trajectory with SGal3:")
    
    # Define velocities at waypoints
    velocities = [
        np.array([0.0, 0.0, 0.0]),  # Start with zero velocity
        np.array([0.5, 0.5, 0.0]),  # Waypoint 1
        np.array([0.0, -0.5, 0.5]), # Waypoint 2
        np.array([0.0, 0.0, 0.0])   # End with zero velocity
    ]
    
    # Define times for waypoints
    times = [0.0, 2.0, 4.0, 6.0]
    
    # Create SGal3 waypoints
    sgal_waypoints = [SGal3(pose, vel, t) for pose, vel, t in zip(waypoints, velocities, times)]
    
    # Generate a smooth spacetime trajectory
    sgal_trajectory = []
    
    for i in range(segments):
        start = sgal_waypoints[i]
        end = sgal_waypoints[i+1]
        
        for t in np.linspace(0, 1, points_per_segment):
            # Interpolate between waypoints
            interp_sgal = lg.interpolate(start, end, t)
            sgal_trajectory.append(interp_sgal)
    
    # Extract information for visualization
    sgal_positions = np.array([g.pose.translation for g in sgal_trajectory])
    sgal_velocities = np.array([g.velocity for g in sgal_trajectory])
    sgal_times = np.array([g.time for g in sgal_trajectory])
    
    # Print some information
    print("Generated a spacetime trajectory with", len(sgal_trajectory), "poses")
    print("Time spans from", sgal_times[0], "to", sgal_times[-1])
    
    # Plot velocity magnitude along the trajectory
    plt.figure(figsize=(10, 6))
    velocity_magnitudes = np.linalg.norm(sgal_velocities, axis=1)
    plt.plot(sgal_times, velocity_magnitudes, 'b-', linewidth=2)
    plt.xlabel('Time')
    plt.ylabel('Velocity Magnitude')
    plt.grid(True)
    plt.title('Velocity Profile of the Trajectory')
    plt.savefig('velocity_profile.png')

def main():
    """Main function to run all examples"""
    print_section("Cosserat Lie Groups - Comprehensive Examples")
    
    # Run the examples
    example_SO2()
    example_SO3()
    example_SE2()
    example_SE3()
    example_SGal3()
    example_Bundle()
    example_utils()
    robot_kinematics_example()
    trajectory_interpolation_example()
    
    print("\nAll examples completed successfully!")
    print("Generated the following visualization files:")
    print("- so2_rotation.png - Demonstration of SO2 rotation")
    print("- so3_rotation.png - Demonstration of SO3 rotation")
    print("- se2_transform.png - Demonstration of SE2 transformation")
    print("- se3_transform.png - Demonstration of SE3 transformation")
    print("- robot_kinematics_2d.png - 2D robot kinematics with SE2")
    print("- robot_kinematics_3d.png - 3D robot kinematics with SE3")
    print("- trajectory_interpolation.png - Trajectory interpolation with SE3")
    print("- velocity_profile.png - Velocity profile of a SGal3 trajectory")

if __name__ == "__main__":
    main()

