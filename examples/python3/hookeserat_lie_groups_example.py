#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HookeSerat Mapping with Lie Groups Example

This example demonstrates the integration of Lie groups with the enhanced
HookeSeratDiscretMapping for advanced beam control and state estimation.

Key features demonstrated:
1. HookeSeratDiscretMapping with Lie groups support
2. SE3 pose representation for beam states
3. Uncertainty propagation using Lie algebra
4. Advanced beam control using Lie group operations
"""

import Sofa
import numpy as np
import math
from Cosserat.LieGroups import SO3, SE3, PoseVel
from cosserat import (BeamGeometryParameters, BeamPhysicsBaseParameters,
                      CosseratGeometry)
from introduction_and_setup import (_add_cosserat_frame, _add_cosserat_state,
                                    _add_rigid_base, add_mini_header)

def print_section(title):
    """Helper to print section titles with formatting"""
    print("\n" + "="*80)
    print(f" {title} ".center(80, "-"))
    print("="*80)

def create_enhanced_beam_scene(root_node, use_lie_groups=True):
    """Create a beam scene using HookeSeratDiscretMapping with Lie groups support"""

    print_section("Creating Enhanced Beam Scene")

    # Configure scene
    add_mini_header(root_node)
    root_node.gravity = [0, -9.81, 0]  # Add gravity

    # Add time integration solver
    solver_node = root_node.addChild("solver_node")
    solver_node.addObject("EulerImplicitSolver", firstOrder="0")
    solver_node.addObject("SparseLDLSolver", name="solver")

    # Define beam geometry with more sections for better control
    beam_geometry_params = BeamGeometryParameters(
        beam_length=20.0,
        nb_section=5,    # 5 sections for physics
        nb_frames=15     # 15 frames for smooth visualization
    )

    beam_geometry = CosseratGeometry(beam_geometry_params)

    print("Beam configuration:")
    print(f"  - Length: {beam_geometry.get_beam_length()}")
    print(f"  - Sections: {beam_geometry.get_number_of_sections()}")
    print(f"  - Frames: {beam_geometry.get_number_of_frames()}")
    print(f"  - Section lengths: {beam_geometry.section_lengths}")

    # Create rigid base
    base_node = _add_rigid_base(solver_node)

    # Define initial bending states with some curvature
    # Using Lie groups concepts: each section has material curvature
    custom_bending_states = []
    for i in range(beam_geometry.get_number_of_sections()):
        # Create a gentle S-curve: alternating positive/negative curvature
        curvature_z = 0.05 * (-1)**i  # Alternating curvature
        custom_bending_states.append([0.0, 0.0, curvature_z])

    print(f"Initial bending states: {custom_bending_states}")

    # Create cosserat state
    bending_node = _add_cosserat_state(solver_node, beam_geometry,
                                       custom_bending_states=custom_bending_states)

    # Create cosserat frame with HookeSeratDiscretMapping
    frame_node = _add_cosserat_frame(base_node, bending_node, beam_geometry,
                                     beam_mass=2.0)

    if use_lie_groups:
        # Add Lie groups-based controller
        add_lie_groups_controller(frame_node, beam_geometry)

    return solver_node, beam_geometry

def add_lie_groups_controller(frame_node, beam_geometry):
    """Add a controller that uses Lie groups for beam control"""

    print_section("Adding Lie Groups Controller")

    # Create controller node
    controller_node = frame_node.addChild("lie_groups_controller")

    # Add a Python controller that uses Lie groups
    controller_node.addObject("PythonScriptController",
                              filename="hookeserat_controller.py",
                              classname="LieGroupsController")

    # Store beam geometry for the controller
    controller_node.addData(name="beam_geometry", value=beam_geometry,
                           help="Beam geometry object for Lie groups operations")

    print("Lie groups controller added to beam")

def demonstrate_lie_groups_operations():
    """Demonstrate Lie groups operations for beam control"""

    print_section("Lie Groups Operations for Beam Control")

    # Example: Represent beam tip pose using SE3
    print("1. Beam Pose Representation using SE3:")

    # Initial beam pose (straight along x-axis)
    initial_pose = SE3()
    print(f"Initial pose (identity): translation={initial_pose.translation}")

    # Apply some deformation (bending)
    # Simulate a beam bent with constant curvature
    curvature = 0.1  # rad/m
    beam_length = 20.0
    total_angle = curvature * beam_length

    # Create rotation around y-axis (bending in x-z plane)
    rotation = SO3(total_angle, np.array([0.0, 1.0, 0.0]))
    translation = np.array([beam_length, 0.0, 0.0])

    bent_pose = SE3(rotation, translation)
    print(f"Bent pose: rotation_angle={total_angle:.3f} rad, translation={bent_pose.translation}")

    # 2. Uncertainty Propagation
    print("\n2. Uncertainty Propagation using Lie Algebra:")

    # Represent pose uncertainty as element of se(3)
    pose_uncertainty = np.array([0.01, 0.01, 0.01,  # position uncertainty
                                0.001, 0.001, 0.001])  # orientation uncertainty

    print(f"Pose uncertainty vector: {pose_uncertainty}")

    # Propagate uncertainty through exponential map
    uncertain_pose = SE3.exp(pose_uncertainty) * bent_pose
    print(f"Uncertain pose translation: {uncertain_pose.translation}")

    # 3. Control using Lie Algebra
    print("\n3. Control using Lie Algebra Operations:")

    # Desired pose (target configuration)
    target_pose = SE3(SO3(math.pi/6, np.array([0.0, 1.0, 0.0])),
                     np.array([18.0, 2.0, 0.0]))

    # Compute pose error in Lie algebra
    pose_error = (target_pose * bent_pose.inverse()).log()
    print(f"Pose error (twist coordinates): {pose_error}")

    # Control law: apply correction proportional to error
    gain = 0.1
    control_input = -gain * pose_error
    print(f"Control input: {control_input}")

    # 4. Velocity Estimation
    print("\n4. Velocity Estimation using PoseVel Bundle:")

    # Create pose-velocity bundle
    velocity = np.array([0.1, 0.05, 0.02])  # Linear velocity
    pose_vel = PoseVel(bent_pose, velocity)

    print(f"Pose-velocity: position={pose_vel.get_pose().translation}, velocity={pose_vel.get_velocity()}")

    # Transform a point using the bundle
    point = np.array([1.0, 0.0, 0.0, 0.1, 0.0, 0.0])  # [position, velocity]
    transformed = pose_vel.act(point)
    print(f"Transformed point-velocity: {transformed}")

def create_comparison_scene():
    """Create a scene comparing traditional vs Lie groups enhanced control"""

    print_section("Creating Comparison Scene")

    root_node = Sofa.Core.Node("comparison_scene")

    # Create two beams side by side
    # Left beam: Traditional approach
    traditional_node = root_node.addChild("traditional_beam")
    traditional_node.addObject("MechanicalObject", template="Rigid3d",
                              position=[-5, 0, 0, 0, 0, 0, 1])  # Offset to the left

    # Right beam: Lie groups enhanced
    enhanced_node = root_node.addChild("enhanced_beam")
    enhanced_node.addObject("MechanicalObject", template="Rigid3d",
                           position=[5, 0, 0, 0, 0, 0, 1])   # Offset to the right

    print("Comparison scene created with traditional and enhanced beams")

    return root_node

def main():
    """Main function demonstrating HookeSerat with Lie groups"""

    print("HookeSerat Mapping with Lie Groups - Advanced Example")
    print("="*80)

    # Demonstrate Lie groups operations
    demonstrate_lie_groups_operations()

    # Create enhanced beam scene
    root_node = Sofa.Core.Node("hookeserat_lie_groups_demo")
    solver_node, beam_geometry = create_enhanced_beam_scene(root_node, use_lie_groups=True)

    # Create comparison scene
    comparison_root = create_comparison_scene()

    print("\n" + "="*80)
    print("Demo setup complete!")
    print("This example demonstrates:")
    print("- HookeSeratDiscretMapping with Lie groups integration")
    print("- SE3 pose representation for beam states")
    print("- Uncertainty propagation using Lie algebra")
    print("- Advanced control using Lie group operations")
    print("- Comparison between traditional and enhanced approaches")
    print("="*80)

    return root_node

if __name__ == "__main__":
    # Run the demonstration
    scene = main()

    # In a real SOFA environment, you would run the simulation:
    # Sofa.Simulation.init(scene)
    # Sofa.Simulation.animate(scene)</content>
</xai:function_call">Now I need to create the controller script that the example references. Let me create the hookeserat_controller.py file. 

<xai:function_call name="write">
<parameter name="filePath">examples/python3/hookeserat_controller.py