#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example demonstrating the SE23 Lie group in the Cosserat plugin.

SE23 represents the Special Euclidean group in (2+1) dimensions, 
which is used for rigid body transformations with velocity in 3D space.

This example demonstrates:
1. Basic SE23 operations
2. Transformation of points with velocities
3. Interpolation between configurations
4. Integration with robot kinematics
5. Visualization of trajectories with velocity vectors
"""

import Sofa
import numpy as np
import math
from Cosserat.LieGroups import SO3, SE3, SE23
import Cosserat.LieGroups as lg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def print_section(title):
    """Helper to print section titles with formatting"""
    print("\n" + "="*80)
    print(f" {title} ".center(80, "-"))
    print("="*80)

def basic_se23_operations():
    """Demonstration of basic SE23 operations"""
    print_section("Basic SE23 Operations")
    
    # Create SE23 elements
    identity = SE23()  # Identity transformation with zero velocity
    
    # Create from rotation, position, and velocity
    r = SO3(math.pi/4, np.array([0, 0, 1]))  # 45 degrees around Z
    t = np.array([1.0, 2.0, 0.0])
    vel = np.array([0.5, 0.2, 0.1])  # Linear velocity
    
    # Method 1: Create from individual components
    g1 = SE23(r, t, vel)
    
    # Method 2: Create from SE3 and velocity
    pose = SE3(r, t)
    g2 = SE23(pose, vel)
    
    # Method 3: Create using utility functions
    g3 = lg.fromComponents(t, r, vel)
    
    # Method 4: Create using Euler angles and velocity
    g4 = lg.fromPositionEulerVelocity(
        position=np.array([1.0, 2.0, 0.0]),
        roll=0.0,
        pitch=0.0,
        yaw=math.pi/4,  # 45 degrees around Z
        velocity=vel
    )
    
    # Verify they're all equivalent
    print("g1 ≈ g2?", g1.isApprox(g2))
    print("g1 ≈ g3?", g1.isApprox(g3))
    print("g1 ≈ g4?", g1.isApprox(g4))
    
    # Access components
    print("\nComponents of g1:")
    print("- Pose matrix:\n", g1.pose.matrix())
    print("- Rotation matrix:\n", g1.pose.rotation.matrix())
    print("- Position:", g1.pose.translation)
    print("- Velocity:", g1.velocity)
    print("- Extended matrix:\n", g1.matrix())
    
    # Convert to parameters
    params = lg.toPositionEulerVelocity(g1)
    print("\nParameters of g1:")
    print("- position:", params[0:3])
    print("- euler angles (rad):", params[3:6])
    print("- euler angles (deg):", params[3:6] * 180/math.pi)
    print("- velocity:", params[6:9])
    
    # Group operations
    # Create another SE23 element
    r2 = SO3(math.pi/6, np.array([1, 0, 0]))  # 30 degrees around X
    t2 = np.array([0.0, 0.0, 1.0])
    vel2 = np.array([0.1, 0.3, 0.2])
    h = SE23(r2, t2, vel2)
    
    # Composition
    g_composed = g1 * h
    print("\nComposition (g1 * h):")
    print("- Position:", g_composed.pose.translation)
    print("- Velocity:", g_composed.velocity)
    
    # Inverse
    g_inv = g1.inverse()
    print("\nInverse of g1:")
    print("- Position:", g_inv.pose.translation)
    print("- Velocity:", g_inv.velocity)
    
    # Verify g1 * g1^(-1) = identity
    g_identity = g1 * g_inv
    print("\ng1 * g1^(-1) ≈ identity?", g_identity.isApprox(identity))
    
    # Lie algebra and exponential map
    # SE23 Lie algebra: [vx, vy, vz, wx, wy, wz, ax, ay, az]
    # where (vx, vy, vz) is linear velocity, (wx, wy, wz) is angular velocity,
    # and (ax, ay, az) is linear acceleration
    
    twist = np.array([0.1, 0.2, 0.3, 0.0, 0.0, math.pi/3, 0.5, 0.1, 0.2])
    g_exp = SE23.exp(twist)
    print("\nExponential map from twist:")
    print("- Position:", g_exp.pose.translation)
    print("- Velocity:", g_exp.velocity)
    
    # Logarithmic map
    log_g1 = g1.log()
    print("\nLogarithm of g1 (twist coordinates):", log_g1)
    
    # Baker-Campbell-Hausdorff formula demonstration
    twist1 = np.array([0.1, 0.0, 0.0, 0.0, 0.0, math.pi/6, 0.1, 0.0, 0.0])
    twist2 = np.array([0.0, 0.1, 0.0, 0.0, math.pi/6, 0.0, 0.0, 0.1, 0.0])
    
    # Compose transformations
    g_t1 = SE23().exp(twist1)
    g_t2 = SE23().exp(twist2)
    g_composed = g_t1 * g_t2
    
    # Calculate BCH approximation
    bch = SE23.BCH(twist1, twist2)
    g_bch = SE23().exp(bch)
    
    print("\nDirect composition vs BCH approximation:")
    print("g_composed.matrix():\n", g_composed.matrix())
    print("g_bch.matrix():\n", g_bch.matrix())
    print("Difference:", np.linalg.norm(g_composed.matrix() - g_bch.matrix()))

def transformation_of_points_with_velocities():
    """Demonstration of transforming points with velocities using SE23"""
    print_section("Transformation of Points with Velocities")
    
    # Create an SE23 transformation
    r = SO3(math.pi/4, np.array([0, 0, 1]))  # 45 degrees around Z
    t = np.array([1.0, 2.0, 0.0])
    vel = np.array([0.5, 0.2, 0.1])  # Linear velocity
    g = SE23(r, t, vel)
    
    # Create a point with velocity
    point = np.array([1.0, 0.0, 0.0])  # Position
    point_vel = np.array([0.1, 0.0, 0.0])  # Velocity
    
    # Combine into a single vector [position, velocity, ...]
    point_vel_full = np.zeros(9)
    point_vel_full[0:3] = point
    point_vel_full[3:6] = point_vel
    
    # Transform the point with velocity
    transformed = g.transform(point_vel_full)
    
    # Extract components
    transformed_point = transformed[0:3]
    transformed_vel = transformed[3:6]
    
    print("Original point:", point)
    print("Original velocity:", point_vel)
    print("Transformed point:", transformed_point)
    print("Transformed velocity:", transformed_vel)
    
    # Verify how the transformation affects velocity
    # Velocity transformation includes both:
    # 1. Rotation of the original velocity
    # 2. Addition of the SE23 velocity component
    expected_rotated_vel = r.act(point_vel)
    expected_vel = expected_rotated_vel + vel
    
    print("\nVerification:")
    print("Expected velocity after rotation:", expected_rotated_vel)
    print("Expected velocity after full transform:", expected_vel)
    print("Actual transformed velocity:", transformed_vel)
    print("Matches:", np.allclose(expected_vel, transformed_vel))
    
    # Visualize the transformation
    plt.figure(figsize=(10, 8))
    
    # Plot original point and velocity
    plt.scatter(point[0], point[1], color='blue', label='Original Point')
    plt.arrow(point[0], point[1], point_vel[0], point_vel[1], 
              color='blue', width=0.02, head_width=0.1, label='Original Velocity')
    
    # Plot transformed point and velocity
    plt.scatter(transformed_point[0], transformed_point[1], color='red', label='Transformed Point')
    plt.arrow(transformed_point[0], transformed_point[1], transformed_vel[0], transformed_vel[1], 
              color='red', width=0.02, head_width=0.1, label='Transformed Velocity')
    
    # Draw coordinate frames
    origin = np.array([0, 0])
    x_axis = np.array([0.5, 0])
    y_axis = np.array([0, 0.5])
    
    # Original frame
    plt.arrow(origin[0], origin[1], x_axis[0], x_axis[1], color='darkblue', width=0.02, head_width=0.1)
    plt.arrow(origin[0], origin[1], y_axis[0], y_axis[1], color='darkgreen', width=0.02, head_width=0.1)
    
    # Transformed frame
    new_origin = g.pose.translation[0:2]  # Just the x,y components
    new_x = g.pose.rotation.act(np.array([0.5, 0, 0]))[0:2]
    new_y = g.pose.rotation.act(np.array([0, 0.5, 0]))[0:2]
    
    plt.arrow(new_origin[0], new_origin[1], new_x[0], new_x[1], color='darkred', width=0.02, head_width=0.1)
    plt.arrow(new_origin[0], new_origin[1], new_y[0], new_y[1], color='darkgreen', width=0.02, head_width=0.1)
    
    # Velocity of the frame itself
    plt.arrow(new_origin[0], new_origin[1], vel[0], vel[1], 
              color='purple', width=0.02, head_width=0.1, label='Frame Velocity')
    
    plt.grid(True)
    plt.axis('equal')
    plt.xlim(-1, 4)
    plt.ylim(-1, 4)
    plt.legend()
    plt.title('SE23 Transformation of Point with Velocity')
    plt.savefig('se23_point_transformation.png')
    
    # Transform multiple points with velocities
    print("\nTransforming multiple points with velocities:")
    
    # Create a grid of points with varying velocities
    x, y = np.meshgrid(np.linspace(-1, 1, 5), np.linspace(-1, 1, 5))
    points = np.column_stack((x.flatten(), y.flatten(), np.zeros_like(x.flatten())))
    
    # Assign velocities based on position (e.g., radial velocities)
    velocities = np.zeros_like(points)
    for i in range(len(points)):
        # Normalize to get direction
        direction = points[i] / (np.linalg.norm(points[i]) + 1e-10)
        # Scale magnitude with distance from origin
        magnitude = np.linalg.norm(points[i]) * 0.2
        velocities[i] = direction * magnitude
    
    # Transform all points with velocities
    transformed_points = np.zeros_like(points)
    transformed_velocities = np.zeros_like(velocities)
    
    for i in range(len(points)):
        point_vel_full = np.zeros(9)
        point_vel_full[0:3] = points[i]
        point_vel_full[3:6] = velocities[i]
        
        transformed = g.transform(point_vel_full)
        transformed_points[i] = transformed[0:3]
        transformed_velocities[i] = transformed[3:6]
    
    # Visualize in 3D
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot original points and velocities
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], color='blue', label='Original Points')
    
    # Plot transformed points and velocities
    ax.scatter(transformed_points[:, 0], transformed_points[:, 1], transformed_points[:, 2], 
               color='red', label='Transformed Points')
    
    # Draw some velocity vectors (only a subset for clarity)
    step = 4
    for i in range(0, len(points), step):
        ax.quiver(points[i, 0], points[i, 1], points[i, 2],
                  velocities[i, 0], velocities[i, 1], velocities[i, 2],
                  color='blue', label='Original Velocity' if i == 0 else None)
        
        ax.quiver(transformed_points[i, 0], transformed_points[i, 1], transformed_points[i, 2],
                  transformed_velocities[i, 0], transformed_velocities[i, 1], transformed_velocities[i, 2],
                  color='red', label='Transformed Velocity' if i == 0 else None)
    
    # Plot the frame velocity
    ax.quiver(g.pose.translation[0], g.pose.translation[1], g.pose.translation[2],
              vel[0], vel[1], vel[2], color='purple', linewidth=2, label='Frame Velocity')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    ax.set_title('SE23 Transformation of Points with Velocities')
    plt.savefig('se23_points_transformation_3d.png')

def interpolation_between_configurations():
    """Demonstration of interpolation between SE23 configurations"""
    print_section("Interpolation Between SE23 Configurations")
    
    # Define two SE23 configurations
    # Start: At origin with 0 velocity
    start = SE23()
    
    # End: Rotated, translated, and with velocity
    r_end = SO3(math.pi/2, np.array([0, 0, 1]))  # 90 degrees around Z
    t_end = np.array([2.0, 1.0, 0.5])
    vel_end = np.array([0.5, 0.2, 0.1])
    end = SE23(r_end, t_end, vel_end)
    
    # Generate a trajectory by interpolating between configurations
    num_steps = 50
    trajectory = []
    times = np.linspace(0, 1, num_steps)
    
    for t in times:
        # Interpolate
        g_interp = lg.interpolate(start, end, t)
        trajectory.append(g_interp)
    
    # Extract data for visualization
    positions = np.array([g.pose.translation for g in trajectory])
    velocities = np.array([g.velocity for g in trajectory])
    
    # Visualization in 3D
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the trajectory
    ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], 'b-', linewidth=2, label='Trajectory')
    
    # Add velocity vectors at regular intervals
    step = num_steps // 10
    for i in range(0, num_steps, step):
        g = trajectory[i]
        ax.quiver(g.pose.translation[0], g.pose.translation[1], g.pose.translation[2],
                  g.velocity[0], g.velocity[1], g.velocity[2],
                  color='red', length=0.2, label='Velocity' if i == 0 else None)
    
    # Also visualize the orientation changes along the trajectory
    for i in range(0, num_steps, step):
        g = trajectory[i]
        pos = g.pose.translation
        
        # Get the orientation axes
        x_dir = g.pose.rotation.act(np.array([0.2, 0, 0]))
        y_dir = g.pose.rotation.act(np.array([0, 0.2, 0]))
        z_dir = g.pose.rotation.act(np.array([0, 0, 0.2]))
        
        # Draw coordinate axes
        ax.quiver(pos[0], pos[1], pos[2], x_dir[0], x_dir[1], x_dir[2], color='darkred')
        ax.quiver(pos[0], pos[1], pos[2], y_dir[0], y_dir[1], y_dir[2], color='darkgreen')
        ax.quiver(pos[0], pos[1], pos[2], z_dir[0], z_dir[1], z_dir[2], color='darkblue')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    ax.set_title('Interpolation Between SE23 Configurations')
    plt.savefig('se23_interpolation.png')
    
    # Plot velocity magnitude over the trajectory
    plt.figure(figsize=(10, 6))
    velocity_magnitudes = np.linalg.norm(velocities, axis=1)
    plt.plot(times, velocity_magnitudes, 'b-', linewidth=2)
    plt.xlabel('Interpolation Parameter (t)')
    plt.ylabel('Velocity Magnitude')
    plt.grid(True)
    plt.title('Velocity Profile During Interpolation')
    plt.savefig('se23_velocity_profile.png')
    
    # Analyze the twist coordinates during interpolation
    twists = np.array([g.log() for g in trajectory])
    
    plt.figure(figsize=(12, 8))
    
    # Plot linear velocity components (first 3)
    plt.subplot(3, 1, 1)
    plt.plot(times, twists[:, 0], 'r-', label='vx')
    plt.plot(times, twists[:, 1], 'g-', label='vy')
    plt.plot(times, twists[:, 2], 'b-', label='vz')
    plt.ylabel('Linear Velocity')
    plt.grid(True)
    plt.legend()
    
    # Plot angular velocity components (middle 3)
    plt.subplot(3, 1, 2)
    plt.plot(times, twists[:, 3], 'r-', label='wx')
    plt.plot(times, twists[:, 4], 'g-', label='wy')
    plt.plot(times, twists[:, 5], 'b-', label='wz')
    plt.ylabel('Angular Velocity')
    plt.grid(True)
    plt.legend()
    
    # Plot acceleration components (last 3)
    plt.subplot(3, 1, 3)
    plt.plot(times, twists[:, 6], 'r-', label='ax')
    plt.plot(times, twists[:, 7], 'g-', label='ay')
    plt.plot(times, twists[:, 8], 'b-', label='az')
    plt.xlabel('Interpolation Parameter (t)')
    plt.ylabel('Acceleration')
    plt.grid(True)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('se23_twist_components.png')

def robot_kinematics_with_velocity():
    """Demonstration of robot kinematics with velocity using SE23"""
    print_section("Robot Kinematics with Velocity Using SE23")
    
    # Define a 2-link robot arm
    link1_length = 1.0
    link2_length = 0.8
    
    # Joint angles and velocities
    theta1 = math.pi/4  # 45 degrees
    theta2 = math.pi/3  # 60 degrees
    
    omega1 = 0.5  # Angular velocity of joint 1 (rad/s)
    omega2 = 0.3  # Angular velocity of joint 2 (rad/s)
    
    # Base frame (identity with zero velocity)
    T0 = SE23()
    
    # First joint: Rotation around Z and extension along X
    # Velocity is induced by the angular velocity
    R1 = SO3(theta1, np.array([0, 0, 1]))
    p1 = np.array([link1_length, 0, 0])
    
    # Calculate velocity for first link
    # v = ω × r where r is the position vector from rotation axis
    # For rotation around Z, the cross product is [-omega*y, omega*x, 0]
    v1 = np.array([-omega1 * 0, omega1 * link1_length/2, 0])
    
    # Create SE23 for first joint
    T1 = SE23(R1, p1, v1)
    
    # Second joint: Rotation around Z and extension along X
    # This is relative to the first joint
    R2 = SO3(theta2, np.array([0, 0, 1]))
    p2 = np.array([link2_length, 0, 0])
    
    # Calculate velocity for second link
    # This includes the velocity from the first joint plus its own rotation
    v2 = np.array([-omega2 * 0, omega2 * link2_length/2, 0])
    
    # Create SE23 for second joint (relative to first joint)
    T2 = SE23(R2, p2, v2)
    
    # Calculate absolute SE23 for each link
    T_to_1 = T0 * T1
    
    # For the second link, we need to compose SE23 carefully
    # The velocity will automatically combine properly
    T_to_2 = T_to_1 * T2
    
    # Calculate positions and velocities at each joint
    p0 = T0.pose.translation
    v0 = T0.velocity
    
    p1_absolute = T_to_1.pose.translation
    v1_absolute = T_to_1.velocity
    
    p2_absolute = T_to_2.pose.translation
    v2_absolute = T_to_2.velocity
    
    print("Base position:", p0)
    print("Base velocity:", v0)
    print("Joint 1 position:", p1_absolute)
    print("Joint 1 velocity:", v1_absolute)
    print("Joint 2 (end-effector) position:", p2_absolute)
    print("Joint 2 (end-effector) velocity:", v2_absolute)
    
    # Visualize the robot with velocities
    plt.figure(figsize=(10, 8))
    
    # Plot the links
    plt.plot([p0[0], p1_absolute[0]], [p0[1], p1_absolute[1]], 'b-', linewidth=3, label='Link 1')
    plt.plot([p1_absolute[0], p2_absolute[0]], [p1_absolute[1], p2_absolute[1]], 
             'g-', linewidth=3, label='Link 2')
    
    # Plot velocity vectors
    # Use a scale factor to make velocities visible
    scale = 0.5
    
    # Plot velocity vectors at the links' centers
    link1_center = (p0 + p1_absolute) / 2
    link2_center = (p1_absolute + p2_absolute) / 2
    
    plt.arrow(link1_center[0], link1_center[1], 
              v1_absolute[0] * scale, v1_absolute[1] * scale, 
              color='blue', width=0.02, head_width=0.1, label='Link 1 Velocity')
    
    plt.arrow(link2_center[0], link2_center[1], 
              v2_absolute[0] * scale, v2_absolute[1] * scale, 
              color='green', width=0.02, head_width=0.1, label='Link 2 Velocity')
    
    # Plot the end-effector velocity
    plt.arrow(p2_absolute[0], p2_absolute[1], 
              v2_absolute[0] * scale, v2_absolute[1] * scale, 
              color='red', width=0.02, head_width=0.1, label='End-Effector Velocity')
    
    # Plot the joints
    plt.scatter([p0[0], p1_absolute[0]], [p0[1], p1_absolute[1]], 
                color='black', s=80, label='Joints')
    plt.scatter([p2_absolute[0]], [p2_absolute[1]], 
                color='orange', s=100, label='End-Effector')
    
    plt.grid(True)
    plt.axis('equal')
    plt.xlim(-0.5, 2.5)
    plt.ylim(-0.5, 2.5)
    plt.legend()
    plt.title('Robot Kinematics with Velocity Using SE23')
    plt.savefig('se23_robot_kinematics.png')
    
    # Generate a trajectory of the robot over time
    print("\nGenerating robot trajectory...")
    
    # Simulate for 10 seconds
    simulation_time = 5.0
    dt = 0.1
    steps = int(simulation_time / dt)
    
    # Storage for trajectory
    joint_angles = []
    end_effector_positions = []
    end_effector_velocities = []
    
    # Initial joint angles
    current_theta1 = theta1
    current_theta2 = theta2
    
    for i in range(steps):
        time = i * dt
        
        # Update joint angles based on angular velocities
        current_theta1 += omega1 * dt
        current_theta2 += omega2 * dt
        
        joint_angles.append([current_theta1, current_theta2])
        
        # Calculate new SE23 for the robot
        R1 = SO3(current_theta1, np.array([0, 0, 1]))
        R2 = SO3(current_theta2, np.array([0, 0, 1]))
        
        # Calculate velocities
        v1 = np.array([-omega1 * 0, omega1 * link1_length/2, 0])
        v2 = np.array([-omega2 * 0, omega2 * link2_length/2, 0])
        
        T1 = SE23(R1, p1, v1)
        T2 = SE23(R2, p2, v2)
        
        T_to_1 = T0 * T1
        T_to_2 = T_to_1 * T2
        
        end_effector_positions.append(T_to_2.pose.translation)
        end_effector_velocities.append(T_to_2.velocity)
    
    # Convert to numpy arrays
    end_effector_positions = np.array(end_effector_positions)
    end_effector_velocities = np.array(end_effector_velocities)
    
    # Visualize the trajectory
    plt.figure(figsize=(10, 8))
    
    # Plot the trajectory
    plt.plot(end_effector_positions[:, 0], end_effector_positions[:, 1], 
             'b-', linewidth=1, label='End-Effector Trajectory')
    
    # Plot velocity vectors at intervals
    step = steps // 10
    for i in range(0, steps, step):
        plt.arrow(end_effector_positions[i, 0], end_effector_positions[i, 1],
                  end_effector_velocities[i, 0] * scale, end_effector_velocities[i, 1] * scale,
                  color='red', width=0.01, head_width=0.05)
    
    # Plot the initial robot position
    plt.plot([p0[0], p1_absolute[0]], [p0[1], p1_absolute[1]], 'k-', linewidth=1, alpha=0.5)
    plt.plot([p1_absolute[0], p2_absolute[0]], [p1_absolute[1], p2_absolute[1]], 
             'k-', linewidth=1, alpha=0.5)
    
    # Plot some intermediate robot positions
    for i in range(0, steps, steps // 5):
        theta1_i = joint_angles[i][0]
        theta2_i = joint_angles[i][1]
        
        # Calculate positions
        R1_i = SO3(theta1_i, np.array([0, 0, 1]))
        p1_i = R1_i.act(p1)
        
        R2_i = SO3(theta2_i, np.array([0, 0, 1]))
        p2_local_i = R2_i.act(p2)
        
        p1_abs_i = p1_i
        p2_abs_i = p1_abs_i + R1_i.act(p2_local_i)
        
        # Plot links
        plt.plot([p0[0], p1_abs_i[0]], [p0[1], p1_abs_i[1]], 'k-', linewidth=1, alpha=0.3)
        plt.plot([p1_abs_i[0], p2_abs_i[0]], [p1_abs_i[1], p2_abs_i[1]], 
                 'k-', linewidth=1, alpha=0.3)
    
    plt.scatter(p0[0], p0[1], color='black', s=80, label='Base Joint')
    
    plt.grid(True)
    plt.axis('equal')
    plt.xlim(-2.5, 2.5)
    plt.ylim(-2.5, 2.5)
    plt.legend()
    plt.title('Robot End-Effector Trajectory with Velocity')
    plt.savefig('se23_robot_trajectory.png')
    
    # Plot velocity magnitude over time
    plt.figure(figsize=(10, 6))
    time_steps = np.arange(steps) * dt
    velocity_magnitudes = np.linalg.norm(end_effector_velocities, axis=1)
    
    plt.plot(time_steps, velocity_magnitudes, 'b-', linewidth=2)
    plt.xlabel('Time (s)')
    plt.ylabel('End-Effector Velocity Magnitude')
    plt.grid(True)
    plt.title('End-Effector Velocity Magnitude Over Time')
    plt.savefig('se23_robot_velocity_profile.png')

def trajectory_with_velocity_vectors():
    """Demonstration of trajectory visualization with velocity vectors"""
    print_section("Trajectory Visualization with Velocity Vectors")
    
    # Define waypoints with velocities using SE23
    waypoints = [
        SE23(
            SE3(SO3(), np.array([0, 0, 0])),  # Identity pose
            np.array([0, 0, 0])  # Zero velocity
        ),
        SE23(
            SE3(SO3(math.pi/4, np.array([0, 0, 1])), np.array([1, 1, 0])),
            np.array([0.5, 0.5, 0])  # Velocity in x,y
        ),
        SE23(
            SE3(SO3(math.pi/2, np.array([0, 1, 0])), np.array([2, 0, 1])),
            np.array([0, -0.5, 0.5])  # Velocity in -y,z
        ),
        SE23(
            SE3(SO3(math.pi/2, np.array([1, 0, 0])), np.array([3, 1, 0])),
            np.array([0, 0, 0])  # Zero velocity
        )
    ]
    
    # Generate a smooth trajectory
    num_points = 200
    trajectory = []
    
    # Interpolate between waypoints
    segments = len(waypoints) - 1
    points_per_segment = num_points // segments
    
    for i in range(segments):
        start = waypoints[i]
        end = waypoints[i+1]
        
        for t in np.linspace(0, 1, points_per_segment):
            # Interpolate between waypoints using SE23 interpolation
            interp_pose = lg.interpolate(start, end, t)
            trajectory.append(interp_pose)
    
    # Extract positions and velocities for visualization
    positions = np.array([pose.pose.translation for pose in trajectory])
    velocities = np.array([pose.velocity for pose in trajectory])
    
    # Create 3D visualization
    fig = plt.figure(figsize=(15, 12))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the trajectory
    ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], 
            'b-', linewidth=2, label='Trajectory')
    
    # Plot waypoints
    waypoint_positions = np.array([pose.pose.translation for pose in waypoints])
    ax.scatter(waypoint_positions[:, 0], waypoint_positions[:, 1], waypoint_positions[:, 2], 
               color='red', s=100, label='Waypoints')
    
    # Plot velocity vectors at intervals
    step = num_points // 20
    for i in range(0, num_points, step):
        pos = positions[i]
        vel = velocities[i]
        
        # Scale velocity for visibility
        scale = 0.3
        ax.quiver(pos[0], pos[1], pos[2], 
                  vel[0] * scale, vel[1] * scale, vel[2] * scale, 
                  color='green', label='Velocity' if i == 0 else None)
    
    # Draw coordinate frames at waypoints
    for i, pose in enumerate(waypoints):
        origin = pose.pose.translation
        r = pose.pose.rotation
        
        # Create coordinate axes
        x_dir = r.act(np.array([0.3, 0, 0]))
        y_dir = r.act(np.array([0, 0.3, 0]))
        z_dir = r.act(np.array([0, 0, 0.3]))
        
        # Draw the axes
        ax.quiver(origin[0], origin[1], origin[2], 
                  x_dir[0], x_dir[1], x_dir[2], color='red')
        ax.quiver(origin[0], origin[1], origin[2], 
                  y_dir[0], y_dir[1], y_dir[2], color='green')
        ax.quiver(origin[0], origin[1], origin[2], 
                  z_dir[0], z_dir[1], z_dir[2], color='blue')
        
        # Draw the velocity vector
        vel = pose.velocity
        ax.quiver(origin[0], origin[1], origin[2], 
                  vel[0], vel[1], vel[2], color='purple', linewidth=2,
                  label='Waypoint Velocity' if i == 0 else None)
    
    # Also draw some coordinate frames along the trajectory
    step = num_points // 10
    for i in range(0, num_points, step):
        pose = trajectory[i]
        origin = pose.pose.translation
        r = pose.pose.rotation
        
        # Create coordinate axes
        x_dir = r.act(np.array([0.2, 0, 0]))
        y_dir = r.act(np.array([0, 0.2, 0]))
        z_dir = r.act(np.array([0, 0, 0.2]))
        
        # Draw the axes
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
    ax.set_title('SE23 Trajectory with Velocity Vectors')
    plt.savefig('se23_trajectory_velocity.png')
    
    # Analyze how velocity changes along the trajectory
    
    # Plot velocity components
    plt.figure(figsize=(12, 8))
    t = np.linspace(0, 1, num_points)
    
    plt.subplot(3, 1, 1)
    plt.plot(t, velocities[:, 0], 'r-', label='X Velocity')
    plt.grid(True)
    plt.ylabel('X Velocity')
    plt.legend()
    
    plt.subplot(3, 1, 2)
    plt.plot(t, velocities[:, 1], 'g-', label='Y Velocity')
    plt.grid(True)
    plt.ylabel('Y Velocity')
    plt.legend()
    
    plt.subplot(3, 1, 3)
    plt.plot(t, velocities[:, 2], 'b-', label='Z Velocity')
    plt.grid(True)
    plt.xlabel('Trajectory Parameter (t)')
    plt.ylabel('Z Velocity')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('se23_velocity_components.png')
    
    # Plot velocity magnitude
    plt.figure(figsize=(10, 6))
    velocity_magnitude = np.linalg.norm(velocities, axis=1)
    plt.plot(t, velocity_magnitude, 'k-', linewidth=2)
    plt.grid(True)
    plt.xlabel('Trajectory Parameter (t)')
    plt.ylabel('Velocity Magnitude')
    plt.title('Velocity Magnitude Along Trajectory')
    plt.savefig('se23_velocity_magnitude.png')
    
    print("Trajectory visualization complete. Files saved:")
    print("- se23_trajectory_velocity.png")
    print("- se23_velocity_components.png")
    print("- se23_velocity_magnitude.png")

def main():
    """Main function to run all examples"""
    print_section("SE23 Examples - Special Euclidean Group in (2+1)D")
    
    # Run the examples
    basic_se23_operations()
    transformation_of_points_with_velocities()
    interpolation_between_configurations()
    robot_kinematics_with_velocity()
    trajectory_with_velocity_vectors()
    
    print("\nAll examples completed successfully!")
    print("Visualization files saved to the current directory")

if __name__ == "__main__":
    main()

