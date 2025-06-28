# -*- coding: utf-8 -*-
"""
Tutorial 04: Advanced Applications
=================================

This final tutorial demonstrates more advanced uses of the Cosserat plugin:
- Multiple beams working together
- More complex force interactions
- Higher-level beam creation functions
- Realistic visualization

Key concepts:
- Creating reusable beam creation functions
- Working with multiple beams
- Different force control approaches
- Extending the plugin to real applications
"""

import os
import sys
from math import pi

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

def create_cosserat_beam(parent_node, position, beam_length, nb_section, nb_frames, 
                         beam_mass=5.0, young_modulus=1.0e3, poisson_ratio=0.4, 
                         beam_radius=1.0, fixed_base=True):
    """
    Create a complete Cosserat beam with all components.
    
    This function encapsulates all the steps to create a beam in one call,
    making it easier to create multiple beams.
    
    Parameters:
        parent_node: The SOFA node to attach this to
        position: Base position and orientation [x,y,z,qx,qy,qz,qw]
        beam_length: Length of the beam
        nb_section: Number of sections for physics
        nb_frames: Number of frames for visualization
        beam_mass: Mass of the beam
        young_modulus: Material stiffness parameter
        poisson_ratio: Material property
        beam_radius: Radius of the beam visualization
        fixed_base: Whether to fix the base with springs
        
    Returns:
        Dictionary with all created nodes
    """
    # Create a container node for the beam
    beam_node = parent_node.addChild(f"beam_{position[0]}_{position[1]}_{position[2]}")
    
    # Create beam geometry
    beam_geometry_params = BeamGeometryParameters(
        beam_length=beam_length,
        nb_section=nb_section,
        nb_frames=nb_frames
    )
    beam_geometry = CosseratGeometry(beam_geometry_params)
    
    # Create the rigid base
    if fixed_base:
        # Use stiffness springs to fix the base
        base_node = create_rigid_base(
            beam_node, 
            node_name="rigid_base",
            positions=position
        )
    else:
        # Create a movable base without springs
        base_node = beam_node.addChild("rigid_base")
        base_node.addObject(
            "MechanicalObject",
            template="Rigid3d",
            name="cosserat_base_mo",
            position=position,
            showObject=True,
            showObjectScale="0.1",
        )
    
    # Initialize straight beam
    custom_bending_states = []
    for i in range(beam_geometry.get_number_of_sections()):
        custom_bending_states.append([0, 0.0, 0.0])
    
    # Create Cosserat state
    bending_node = create_cosserat_state(
        beam_node, 
        beam_geometry, 
        node_name="cosserat_states",
        custom_bending_states=custom_bending_states
    )
    
    # Customize the force field with provided material properties
    bending_node.BeamHookeLawForceField.youngModulus = young_modulus
    bending_node.BeamHookeLawForceField.poissonRatio = poisson_ratio
    bending_node.BeamHookeLawForceField.radius = beam_radius
    
    # Create Cosserat frames
    frame_node = create_cosserat_frame(
        base_node, 
        bending_node, 
        beam_geometry, 
        node_name="cosserat_frames",
        beam_mass=beam_mass
    )
    
    # Return all components for later reference
    return {
        "node": beam_node,
        "base": base_node,
        "state": bending_node,
        "frames": frame_node,
        "geometry": beam_geometry
    }

class MultiForceController(Sofa.Core.Controller):
    """
    Advanced controller that can manipulate multiple beams.
    
    This controller demonstrates how to:
    1. Control multiple beams at once
    2. Apply different force patterns
    3. Create more complex animations
    """
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        
        # Store a list of beams to control
        self.beams = kwargs["beams"]
        self.force_nodes = kwargs["force_nodes"]
        self.animation_type = kwargs.get("animation_type", 1)
        
        # Animation parameters
        self.time = 0.0
        self.frequency = 0.5  # Hz
        self.amplitude = 5.0  # Force magnitude
        
        print(f"🎮 Multi-force controller initialized with animation type {self.animation_type}")
        print(f"   Controlling {len(self.beams)} beams")

    def onAnimateBeginEvent(self, event):
        """Update forces at the beginning of each animation step."""
        # Increment time
        dt = self.getContext().getDt()
        self.time += dt
        
        # Apply different animation patterns based on type
        if self.animation_type == 1:
            # Sinusoidal pattern - each beam gets a phase-shifted force
            self._apply_wave_pattern()
        elif self.animation_type == 2:
            # Circular pattern - beams move in a circular pattern
            self._apply_circular_pattern()
        elif self.animation_type == 3:
            # Alternating pattern - beams take turns being active
            self._apply_alternating_pattern()

    def _apply_wave_pattern(self):
        """Apply a sinusoidal wave pattern to the beams."""
        for i, (beam, force_node) in enumerate(zip(self.beams, self.force_nodes)):
            # Calculate phase offset for this beam
            phase = (2 * pi * i) / len(self.beams)
            
            # Calculate force based on sine wave
            force_value = self.amplitude * sin(2 * pi * self.frequency * self.time + phase)
            
            # Get the last frame's orientation
            frames = beam["frames"].FramesMO
            last_frame_idx = beam["geometry"].get_number_of_frames() - 1
            position = frames.position[last_frame_idx]
            orientation = Quat(position[3], position[4], position[5], position[6])
            
            # Apply force in local Y direction
            local_force = [0.0, force_value, 0.0]
            global_force = orientation.rotate(local_force)
            
            # Update the force field
            with force_node.forces.writeable() as force:
                for j in range(3):
                    force[0][j] = global_force[j]

    def _apply_circular_pattern(self):
        """Apply forces to make the beams move in a circular pattern."""
        for i, (beam, force_node) in enumerate(zip(self.beams, self.force_nodes)):
            # Calculate angle around the circle
            angle = 2 * pi * self.frequency * self.time
            
            # Calculate force components to create circular motion
            force_x = self.amplitude * cos(angle)
            force_y = self.amplitude * sin(angle)
            
            # Apply force
            with force_node.forces.writeable() as force:
                force[0][0] = force_x
                force[0][1] = force_y
                force[0][2] = 0.0

    def _apply_alternating_pattern(self):
        """Each beam takes turns being active."""
        # Determine which beam should be active
        period = 2.0  # seconds per beam
        active_idx = int((self.time / period) % len(self.beams))
        
        # Apply forces only to the active beam
        for i, (beam, force_node) in enumerate(zip(self.beams, self.force_nodes)):
            with force_node.forces.writeable() as force:
                if i == active_idx:
                    # Active beam gets Y force
                    force[0][1] = self.amplitude
                else:
                    # Inactive beams get no force
                    force[0][0] = 0.0
                    force[0][1] = 0.0
                    force[0][2] = 0.0


def createScene(root_node):
    """
    Create a scene with multiple beams demonstrating advanced concepts.
    """
    # Configure scene with all required plugins
    add_required_plugins(root_node)

    # Add gravity
    root_node.gravity = [0, -9.81, 0]

    # Create solver node for time integration
    solver_node = create_solver_node(root_node)
    
    # === CREATE MULTIPLE BEAMS ===
    # We'll create three beams in different positions
    
    # Beam 1: Base configuration
    beam1 = create_cosserat_beam(
        parent_node=solver_node,
        position=[-20, 0, 0, 0, 0, 0, 1],  # Left side
        beam_length=30.0,
        nb_section=15,
        nb_frames=30,
        beam_mass=5.0,
        young_modulus=1.0e3,
        poisson_ratio=0.4
    )
    
    # Beam 2: Stiffer material
    beam2 = create_cosserat_beam(
        parent_node=solver_node,
        position=[0, 0, 0, 0, 0, 0, 1],  # Center
        beam_length=30.0,
        nb_section=15,
        nb_frames=30,
        beam_mass=5.0,
        young_modulus=5.0e3,  # 5x stiffer
        poisson_ratio=0.4
    )
    
    # Beam 3: Thicker beam
    beam3 = create_cosserat_beam(
        parent_node=solver_node,
        position=[20, 0, 0, 0, 0, 0, 1],  # Right side
        beam_length=30.0,
        nb_section=15,
        nb_frames=30,
        beam_mass=10.0,  # Heavier
        young_modulus=1.0e3,
        poisson_ratio=0.4,
        beam_radius=2.0  # Thicker
    )
    
    # Store all beams in a list for easy access
    beams = [beam1, beam2, beam3]
    
    # === ADD FORCE FIELDS TO EACH BEAM ===
    force_nodes = []
    
    for beam in beams:
        # Get the tip frame index
        tip_frame_index = beam["geometry"].get_number_of_frames() - 1
        
        # Add a constant force field to the tip
        force_node = beam["frames"].addObject(
            'ConstantForceField', 
            name='constForce', 
            showArrowSize=1.0,
            indices=[tip_frame_index],
            forces=[[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]  # Initial force is zero
        )
        
        force_nodes.append(force_node)
    
    # Add the multi-force controller to orchestrate all beams
    animation_type = 1  # Try different animation types (1, 2, or 3)
    
    controller = MultiForceController(
        name="multiForceController",
        beams=beams,
        force_nodes=force_nodes,
        animation_type=animation_type
    )
    root_node.addObject(controller)
    
    # === VISUAL ENHANCEMENTS ===
    # Add a floor for reference
    floor = root_node.addChild("floor")
    floor.addObject("MeshOBJLoader", name="loader", filename="mesh/floor.obj", scale=100)
    floor.addObject("OglModel", src="@loader", color=[0.5, 0.5, 0.5, 1.0])
    
    # === INSTRUCTIONS ===
    print("\n🚀 Advanced Cosserat Tutorial")
    print(f"  - Created 3 beams with different properties:")
    print(f"    - Beam 1: Standard configuration")
    print(f"    - Beam 2: Stiffer material (5x)")
    print(f"    - Beam 3: Thicker and heavier")
    print(f"  - Animation Type: {animation_type}")
    print(f"    - Type 1: Sinusoidal wave pattern")
    print(f"    - Type 2: Circular motion pattern")
    print(f"    - Type 3: Alternating active beam")
    print("  - Press 'Animate' to start the simulation\n")

    return root_node
    
# Import missing functions
def sin(angle):
    """Sine function for the MultiForceController."""
    import math
    return math.sin(angle)
    
def cos(angle):
    """Cosine function for the MultiForceController."""
    import math
    return math.cos(angle)
