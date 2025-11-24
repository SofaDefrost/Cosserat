#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Multi-Section Beam Support Example

This example demonstrates enhanced multi-section beam support with section-specific
properties in the HookeSerat mapping. It shows how to:

1. Define beams with different properties per section
2. Handle material discontinuities
3. Optimize computation for heterogeneous beams
4. Validate section continuity and compatibility
"""

import numpy as np
import math
from typing import List, Dict, Any

def print_section(title):
    """Helper to print section titles with formatting"""
    print("\n" + "="*80)
    print(f" {title} ".center(80, "-"))
    print("="*80)

class SectionProperties:
    """Properties for a single beam section"""

    def __init__(self, length: float, young_modulus: float, radius: float,
                 poisson_ratio: float = 0.4, density: float = 1000.0):
        self.length = length
        self.young_modulus = young_modulus
        self.radius = radius
        self.poisson_ratio = poisson_ratio
        self.density = density

        # Derived properties
        self.area = math.pi * radius**2
        self.second_moment = math.pi * radius**4 / 4  # For circular cross-section
        self.polar_moment = math.pi * radius**4 / 2

    def get_stiffness_matrix(self):
        """Compute local stiffness matrix for this section"""
        E = self.young_modulus
        G = E / (2 * (1 + self.poisson_ratio))  # Shear modulus
        I = self.second_moment
        J = self.polar_moment
        A = self.area

        # Local stiffness matrix for Cosserat beam (simplified)
        # [k_xx, k_xy, k_xz, k_θx, k_θy, k_θz]
        k_local = np.array([
            [E * A / self.length, 0, 0, 0, 0, 0],  # Extension
            [0, 12 * E * I / self.length**3, 0, 0, 0, 6 * E * I / self.length**2],  # Bending y
            [0, 0, 12 * E * I / self.length**3, 0, -6 * E * I / self.length**2, 0],  # Bending z
            [0, 0, 0, G * J / self.length, 0, 0],  # Torsion
            [0, 0, -6 * E * I / self.length**2, 0, 4 * E * I / self.length, 0],  # Bending y-z coupling
            [0, 6 * E * I / self.length**2, 0, 0, 0, 4 * E * I / self.length]   # Bending z-y coupling
        ])

        return k_local

class MultiSectionBeam:
    """Beam composed of multiple sections with different properties"""

    def __init__(self, sections: List[SectionProperties]):
        self.sections = sections
        self.validate_sections()
        self.compute_global_properties()

    def validate_sections(self):
        """Validate section compatibility and continuity"""
        if not self.sections:
            raise ValueError("Beam must have at least one section")

        # Check for reasonable property ranges
        for i, section in enumerate(self.sections):
            if section.length <= 0:
                raise ValueError(f"Section {i}: length must be positive")
            if section.young_modulus <= 0:
                raise ValueError(f"Section {i}: Young modulus must be positive")
            if section.radius <= 0:
                raise ValueError(f"Section {i}: radius must be positive")

        print(f"Validated {len(self.sections)} sections")

    def compute_global_properties(self):
        """Compute global beam properties"""
        self.total_length = sum(s.length for s in self.sections)
        self.section_positions = []  # Start position of each section

        current_pos = 0.0
        for section in self.sections:
            self.section_positions.append(current_pos)
            current_pos += section.length

        # Compute center of mass
        total_mass = sum(s.length * s.density * s.area for s in self.sections)
        com_pos = 0.0
        for i, section in enumerate(self.sections):
            section_mass = section.length * section.density * section.area
            com_pos += section_mass * (self.section_positions[i] + section.length / 2)
        self.center_of_mass = com_pos / total_mass if total_mass > 0 else 0.0

        print(f"Total length: {self.total_length:.3f}")
        print(f"Center of mass at: {self.center_of_mass:.3f}")

    def get_section_at_position(self, position: float) -> tuple:
        """Get section index and local position for global coordinate"""
        if position < 0 or position > self.total_length:
            raise ValueError(f"Position {position} outside beam range [0, {self.total_length}]")

        cumulative_length = 0.0
        for i, section in enumerate(self.sections):
            if position <= cumulative_length + section.length:
                local_pos = position - cumulative_length
                return i, local_pos
            cumulative_length += section.length

        # Should not reach here
        raise RuntimeError("Position mapping failed")

    def get_global_stiffness_matrix(self):
        """Assemble global stiffness matrix from section contributions"""
        n_sections = len(self.sections)
        # For a beam with n sections, we have n+1 nodes (including endpoints)
        n_dofs = 6 * (n_sections + 1)  # 6 DOFs per node (3 translation + 3 rotation)

        K_global = np.zeros((n_dofs, n_dofs))

        # Assemble stiffness matrices for each section
        for i, section in enumerate(self.sections):
            K_local = section.get_stiffness_matrix()

            # Map local DOFs to global DOFs
            # Node i and node i+1 for section i
            global_dofs = [6*i + j for j in range(6)] + [6*(i+1) + j for j in range(6)]

            for local_i in range(12):
                for local_j in range(12):
                    global_i = global_dofs[local_i]
                    global_j = global_dofs[local_j]
                    K_global[global_i, global_j] += K_local[local_i % 6, local_j % 6]

        return K_global

    def apply_boundary_conditions(self, stiffness_matrix, fixed_nodes=None):
        """Apply boundary conditions to stiffness matrix"""
        if fixed_nodes is None:
            fixed_nodes = [0]  # Fix first node by default

        K_modified = stiffness_matrix.copy()
        n_dofs = K_modified.shape[0]

        # Remove rows and columns corresponding to fixed DOFs
        free_dofs = []
        for i in range(n_dofs // 6):  # For each node
            if i not in fixed_nodes:
                free_dofs.extend([6*i + j for j in range(6)])

        K_free = K_modified[np.ix_(free_dofs, free_dofs)]
        return K_free, free_dofs

class BeamOptimizer:
    """Optimizer for multi-section beam performance"""

    def __init__(self, beam: MultiSectionBeam):
        self.beam = beam

    def optimize_section_distribution(self, target_stiffness_profile):
        """Optimize section properties for desired stiffness profile"""
        print("Optimizing section distribution...")

        # Simple optimization: adjust Young modulus to match target stiffness
        for i, section in enumerate(self.beam.sections):
            target_k = target_stiffness_profile[i] if i < len(target_stiffness_profile) else target_stiffness_profile[-1]

            # Scale Young modulus to achieve target stiffness
            current_k = section.young_modulus * section.area / section.length
            scale_factor = target_k / current_k if current_k > 0 else 1.0

            section.young_modulus *= scale_factor
            print(f"Section {i}: E scaled by {scale_factor:.3f}")

    def check_material_continuity(self):
        """Check for material discontinuities that might cause stress concentrations"""
        print("Checking material continuity...")

        discontinuities = []
        for i in range(len(self.beam.sections) - 1):
            curr_section = self.beam.sections[i]
            next_section = self.beam.sections[i+1]

            # Check Young modulus discontinuity
            e_ratio = next_section.young_modulus / curr_section.young_modulus
            if abs(e_ratio - 1.0) > 0.5:  # More than 50% change
                discontinuities.append({
                    'type': 'Young_modulus',
                    'location': i,
                    'ratio': e_ratio
                })

            # Check radius discontinuity
            r_ratio = next_section.radius / curr_section.radius
            if abs(r_ratio - 1.0) > 0.2:  # More than 20% change
                discontinuities.append({
                    'type': 'radius',
                    'location': i,
                    'ratio': r_ratio
                })

        if discontinuities:
            print(f"Found {len(discontinuities)} material discontinuities:")
            for disc in discontinuities:
                print(f"  {disc['type']} change at section {disc['location']}: {disc['ratio']:.3f}")
        else:
            print("No significant material discontinuities found")

        return discontinuities

def create_heterogeneous_beam():
    """Create a beam with different properties in each section"""

    print_section("Creating Heterogeneous Multi-Section Beam")

    # Define sections with different properties
    sections = [
        SectionProperties(length=5.0, young_modulus=200e9, radius=0.01, density=7800),  # Steel
        SectionProperties(length=8.0, young_modulus=70e9, radius=0.015, density=2700),  # Aluminum
        SectionProperties(length=3.0, young_modulus=400e9, radius=0.008, density=1600), # Carbon fiber
        SectionProperties(length=4.0, young_modulus=100e9, radius=0.012, density=1200), # Composite
    ]

    # Create multi-section beam
    beam = MultiSectionBeam(sections)

    # Display section information
    print("\nSection Properties:")
    for i, section in enumerate(beam.sections):
        print(f"Section {i}: L={section.length:.1f}m, E={section.young_modulus/1e9:.0f}GPa, "
              f"r={section.radius*1000:.0f}mm, ρ={section.density:.0f}kg/m³")

    return beam

def demonstrate_beam_analysis():
    """Demonstrate analysis of multi-section beam"""

    print_section("Multi-Section Beam Analysis")

    beam = create_heterogeneous_beam()

    # Test position mapping
    test_positions = [2.0, 7.0, 12.0, 18.0]
    print("\nPosition Mapping Test:")
    for pos in test_positions:
        try:
            section_idx, local_pos = beam.get_section_at_position(pos)
            print(f"Global pos {pos:.1f}m -> Section {section_idx}, local pos {local_pos:.1f}m")
        except ValueError as e:
            print(f"Global pos {pos:.1f}m -> {e}")

    # Compute global stiffness matrix
    print("\nComputing Global Stiffness Matrix...")
    K_global = beam.get_global_stiffness_matrix()
    print(f"Global stiffness matrix shape: {K_global.shape}")

    # Apply boundary conditions (fix first node)
    K_free, free_dofs = beam.apply_boundary_conditions(K_global, fixed_nodes=[0])
    print(f"Free DOFs after boundary conditions: {len(free_dofs)}")

    # Check condition number (measure of numerical stability)
    cond_num = np.linalg.cond(K_free) if K_free.size > 0 else float('inf')
    print(f"Condition number of free stiffness matrix: {cond_num:.2e}")

    return beam

def demonstrate_beam_optimization():
    """Demonstrate beam optimization"""

    print_section("Beam Optimization")

    beam = create_heterogeneous_beam()
    optimizer = BeamOptimizer(beam)

    # Check for material discontinuities
    optimizer.check_material_continuity()

    # Define target stiffness profile (increasing stiffness)
    target_stiffness = [1.0, 1.5, 2.0, 2.5]  # Relative stiffness values

    print(f"\nTarget stiffness profile: {target_stiffness}")

    # Optimize section distribution
    optimizer.optimize_section_distribution(target_stiffness)

    # Re-check discontinuities after optimization
    print("\nAfter optimization:")
    optimizer.check_material_continuity()

def simulate_beam_loading():
    """Simulate loading on multi-section beam"""

    print_section("Beam Loading Simulation")

    beam = create_heterogeneous_beam()

    # Define load case: point load at center
    load_position = beam.total_length / 2
    load_magnitude = 1000.0  # 1000 N downward

    print(f"Applying {load_magnitude}N load at position {load_position:.1f}m")

    # Find which section contains the load
    section_idx, local_pos = beam.get_section_at_position(load_position)
    print(f"Load applied to section {section_idx} at local position {local_pos:.1f}m")

    # Simple beam deflection calculation (Euler-Bernoulli)
    # This is a simplified analysis - real implementation would solve the system
    deflections = []
    for i, section in enumerate(beam.sections):
        # Simplified: assume uniform load distribution for demo
        EI = section.young_modulus * section.second_moment
        L = section.length

        # Maximum deflection for cantilever beam with end load
        # w_max = (P * L^3) / (3 * EI) for cantilever
        deflection = (load_magnitude * L**3) / (3 * EI) if i == section_idx else 0.0
        deflections.append(deflection)
        print(f"Section {i} max deflection: {deflection*1000:.3f}mm")

    return beam

def main():
    """Main function demonstrating multi-section beam support"""

    print("Multi-Section Beam Support Example")
    print("="*80)

    # Create and analyze heterogeneous beam
    beam = demonstrate_beam_analysis()

    # Demonstrate optimization
    demonstrate_beam_optimization()

    # Simulate loading
    simulate_beam_loading()

    print("\n" + "="*80)
    print("Multi-section beam demonstration completed!")
    print("Features demonstrated:")
    print("- Heterogeneous beam construction")
    print("- Section property validation")
    print("- Global stiffness matrix assembly")
    print("- Material discontinuity analysis")
    print("- Beam optimization")
    print("- Load case simulation")
    print("="*80)

if __name__ == "__main__":
    main()</content>
</xai:function_call">Now let me create the performance optimization example. This will show techniques for optimizing the HookeSerat mapping for real-time performance. 

<xai:function_call name="write">
<parameter name="filePath">examples/python3/performance_optimization_example.py