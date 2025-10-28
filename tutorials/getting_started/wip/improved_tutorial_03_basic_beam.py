# -*- coding: utf-8 -*-
"""
Tutorial 03: Applied Forces and Interactions
==========================================

This tutorial builds on the previous tutorials by adding:
- Different types of forces applied to the beam
- A force controller for dynamic force adjustment
- Keyboard interaction to modify forces in real-time
- Visualization of applied forces

Key concepts:
- Using ConstantForceField for simple forces
- Creating a custom controller for complex forces
- Applying forces based on beam orientation
- Interactive force control
"""

import os
import sys
from math import sqrt

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

import Sofa
from splib3.numerics import Quat

from cosserat import BeamGeometryParameters, CosseratGeometry

# Import helper functions from the previous tutorials
from tutorials.getting_started.wip.improved_tutorial_00_basic_beam import (add_required_plugins,
                                                                           create_rigid_base,
                                                                           create_cosserat_state,
                                                                           create_cosserat_frame)
from tutorials.getting_started.wip.improved_tutorial_02_basic_beam import create_solver_node

# Damping parameter for dynamics
v_damping_param: float = 8.e-1

class ForceController(Sofa.Core.Controller):
    """
    Controller to apply and adjust forces to the beam.
    
    This controller can:
    1. Apply different types of forces
    2. Adjust force magnitude over time
    3. Respond to keyboard input
    """
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        
        # Store references to scene objects
        self.forceNode = kwargs["forceNode"]
        self.frames = kwargs["frame_node"].FramesMO
        self.force_type = kwargs["force_type"]
        self.tip_controller = kwargs.get("tip_controller", None)
        self.geoParams = kwargs["geoParams"]
        
        # Configure force parameters
        self.nb_frames = self.geoParams.nb_frames - 1  # Last frame index
        self.applyForce = True                         # Whether to apply force
        self.forceCoeff = 0.0                          # Initial force magnitude
        self.theta = 0.1                               # Rotation parameter
        self.incremental = 0.01                        # Force increment per step
        
        print(f"🎮 Force controller initialized with force type {self.force_type}")
        print(f"   Press '+' to increase force, '-' to decrease force")

    def onAnimateEndEvent(self, event):
        """Called at the end of each animation step to update forces."""
        # Update force coefficient based on whether we're applying force
        if self.applyForce:
            self.forceCoeff += self.incremental
        else:
            self.forceCoeff = max(0.0, self.forceCoeff - self.incremental)
        
        # Apply the appropriate force type
        if self.force_type == 1:
            # Simple angular force
            self.incremental = 0.1
            self.compute_force()
            
        elif self.force_type == 2:
            # Force orthogonal to beam orientation
            self.incremental = 0.4
            self.compute_orthogonal_force()
            
        elif self.force_type == 3:
            # Rotate the tip controller
            self.rotate_force()

    def compute_force(self):
        """Apply a simple angular force to the beam tip."""
        with self.forceNode.forces.writeable() as force:
            # Apply force with equal components in Y and Z angular directions
            vec = [
                0.0,  # No X force
                0.0,  # No Y force
                0.0,  # No Z force
                0.0,  # No X angular force
                self.forceCoeff / sqrt(2),  # Y angular force
                self.forceCoeff / sqrt(2),  # Z angular force
            ]
            # Update the force vector
            for i, v in enumerate(vec):
                force[0][i] = v
                
        if self.forceCoeff > 0:
            print(f"💪 Applied Type 1 force: Angular force with magnitude {self.forceCoeff:.2f}")

    def compute_orthogonal_force(self):
        """
        Apply a force that remains orthogonal to the beam's current axis.
        This uses the orientation of the last frame to determine force direction.
        """
        # Get the last frame's position and orientation
        position = self.frames.position[self.nb_frames]
        orientation = Quat(
            position[3], position[4], position[5], position[6]
        )
        
        # Calculate force vector in frame's local coordinates, then rotate to global
        with self.forceNode.forces.writeable() as force:
            # Apply force in the local Y direction (orthogonal to beam axis)
            vec = orientation.rotate([0.0, self.forceCoeff * 5.0e-2, 0.0])
            
            # Update the linear force components
            for count in range(3):
                force[0][count] = vec[count]
                
        if self.forceCoeff > 0:
            print(f"💪 Applied Type 2 force: Orthogonal force with magnitude {self.forceCoeff:.2f}")

    def rotate_force(self):
        """
        Rotate the tip controller around the X axis.
        This is used for more complex beam manipulation.
        """
        if self.tip_controller and self.forceCoeff <= 100.0 * 3.14159:
            with self.tip_controller.position.writeable() as position:
                # Get orientation of the last frame
                last_frame = self.frames.position[self.nb_frames]
                vec = Quat(
                    last_frame[3], last_frame[4], last_frame[5], last_frame[6]
                )
                
                # Apply rotation around X axis
                vec.rotateFromEuler([self.theta, 0.0, 0.0])
                vec.normalize()
                
                # Update the controller orientation
                for i, v in enumerate(vec):
                    position[0][i + 3] = v
                    
            print(f"🔄 Rotated tip controller, angle: {self.forceCoeff:.2f}")

    def onKeypressedEvent(self, event):
        """Handle keyboard input to control force application."""
        key = event["key"]
        if key == "+":
            # Increase force
            self.applyForce = True
            print("⬆️ Increasing force")
        elif key == "-":
            # Decrease force
            self.applyForce = False
            print("⬇️ Decreasing force")


def createScene(root_node):
    """
    Create a scene with a beam and various force interactions.
    """
    # Configure scene with all required plugins
    add_required_plugins(root_node)

    # Add gravity
    root_node.gravity = [0, -9.81, 0]

    # Create solver node for time integration
    solver_node = create_solver_node(root_node)
    
    # Create a beam with moderate discretization
    beam_geometry_params = BeamGeometryParameters(
        beam_length=30.0,  # Total beam length
        nb_section=10,     # Number of sections for physics
        nb_frames=20,      # Number of frames for visualization
    )

    # Create geometry object
    beam_geometry = CosseratGeometry(beam_geometry_params)

    print(f"🚀 Created beam with:")
    print(f"   - Length: {beam_geometry.get_beam_length()}")
    print(f"   - Sections: {beam_geometry.get_number_of_sections()}")
    print(f"   - Frames: {beam_geometry.get_number_of_frames()}")
    print(f"   - Mass will be distributed across frames")

    # Create rigid base
    base_node = create_rigid_base(solver_node)

    # Initialize with straight beam
    custom_bending_states = []
    for i in range(beam_geometry.get_number_of_sections()):
        custom_bending_states.append([0, 0.0, 0.0])

    # Create cosserat state
    bending_node = create_cosserat_state(
        solver_node, 
        beam_geometry, 
        node_name="cosserat_states",
        custom_bending_states=custom_bending_states
    )

    # Create cosserat frame with mass
    frame_node = create_cosserat_frame(
        base_node, 
        bending_node, 
        beam_geometry, 
        beam_mass=5.0
    )

    # === FORCE APPLICATION ===
    # Set up initial force values (no force)
    force_null = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
    # Add a constant force field to the tip
    tip_frame_index = beam_geometry.get_number_of_frames() - 1
    const_force_node = frame_node.addObject(
        'ConstantForceField', 
        name='constForce', 
        showArrowSize=1.0,  # Make the force visible
        indices=[tip_frame_index],  # Apply to the last frame
        forces=[force_null]  # Initial force is zero
    )
    
    # Create a rigid body to visualize the tip controller (for force type 3)
    tip_controller = root_node.addChild('tip_controller')
    controller_state = tip_controller.addObject(
        'MechanicalObject', 
        template='Rigid3d', 
        name="controlEndEffector",
        showObjectScale=0.5,  # Make it visible
        position=[beam_geometry.get_beam_length(), 0, 0, 0, 0, 0, 1],
        showObject=True
    )
    
    # Add the force controller to dynamically update forces
    # Try different force types: 1, 2, or 3
    force_type = 1  # Change this to experiment with different force types
    
    controller = ForceController(
        name="forceController",
        forceNode=const_force_node,
        frame_node=frame_node,
        force_type=force_type,
        tip_controller=tip_controller,
        geoParams=beam_geometry_params
    )
    root_node.addObject(controller)
    
    # Instructions for the user
    print("\n⚡ Force Application Tutorial")
    print(f"  - Force Type: {force_type}")
    print("  - Use '+' key to increase force")
    print("  - Use '-' key to decrease force")
    print("  - Watch how the beam responds to different force types!\n")

    return root_node
