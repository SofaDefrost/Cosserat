# Cosserat Plugin

[![Gitter](https://img.shields.io/badge/chat-on_Gitter-ff69b4.svg)](https://app.gitter.im/#/room/#sofa-framework_cosserat-needle-insertion:gitter.im)
[![Support](https://img.shields.io/badge/support-on_GitHub_Discussions-blue.svg)](https://github.com/sofa-framework/sofa/discussions/categories/cosserat)
![download](https://img.shields.io/github/downloads/SofaDefrost/Cosserat/total.svg)
![forks](https://img.shields.io/github/forks/SofaDefrost/Cosserat.svg)
![stars](https://img.shields.io/github/stars/SofaDefrost/Cosserat.svg)

## Overview

An open-source plugin compatible with the Sofa framework for simulating 1D objects. It enables modeling of rigid and flexible 1D entities such as rods, wires, or needles using Cosserat beam theory. You can create scenes in Python or XML, or develop new C++ components. Contributions are welcome!

## Theory and Background

The Cosserat model simulates deformation of robot bodies with rod-like geometry and mechanical properties. It replicates nonlinear deformations including bending, torsion, extension, and shearing.

Key feature: Each material point is treated as a rigid body with 6 degrees of freedom (3 translations + 3 rotations), unlike traditional models that treat points as particles with 3 translations.

This creates a framework similar to articulated solids, where relative positions are defined by strain states. Ideal for modeling concentric tube robots, cable-actuated continuum robots, or pneumatic soft robots with constant cross-sections.

For theoretical details: [Theory](examples/python3/tutorial/Writerside/topics/Theory.md)

## Features

### Piece-wise Constant Strain (PCS)
Based on:
- **Discrete cosserat approach for soft robot dynamics: A new piece-wise constant strain model with torsion and shears** [Link](https://ieeexplore.ieee.org/document/7759808)
- **Coupling numerical deformable models in global and reduced coordinates for the simulation of the direct and inverse kinematics of Soft Robots** [Link](https://hal.archives-ouvertes.fr/hal-03192168/document)

<div align="center">
  <a href="https://www.youtube.com/watch?v=qwzKAgw31pU"><img src="https://img.youtube.com/vi/qwzKAgw31pU/0.jpg" alt="Video demonstration"></a>
</div>

### Piece-wise Non-constant Strain
Allows strain to vary continuously within each beam section, enabling more accurate modeling of complex deformations compared to the constant strain assumption. Useful for beams with varying material properties or under non-uniform loads.

### DCM with Plastic Model
Extends the Discrete Cosserat Model to include plastic behavior, simulating permanent deformations in materials that exceed elastic limits. This is essential for modeling materials that undergo irreversible changes under stress.

## Use Cases and Applications

### Cochlear Implant Modeling
Using Discrete Cosserat Model (DCM):

| View 1 | View 2 | View 3 |
|--------|--------|--------|
| ![DCM Implant View 1](tutorial/images/multiSectionWithColorMap1.png) | ![DCM Implant View 2](tutorial/images/multiSectionWithColorMap2.png) | ![DCM Implant View 3](tutorial/images/multiSectionWithColorMap3.png) |

### Cable Modeling for Deformable Robot Control

| Soft Gripper Simulation | Model Convergence Study |
|-------------------------|-------------------------|
| ![Soft Gripper](tutorial/images/cosseratgripper_2.png) | ![Convergence Study](tutorial/images/tenCossseratSections.png) |

### Actuation Examples

| Cable Actuation 1 | Cable Actuation 2 | Cable Actuation 3 |
|-------------------|-------------------|-------------------|
| ![Actuation 1](tutorial/images/actuationConstraint_2.png) | ![Actuation 2](tutorial/images/circleActuationConstraint.png) | ![Actuation 3](tutorial/images/actuationConstraint_1.png) |

### Animation Examples

| Example 1 | Example 2 | Example 3 |
|-----------|-----------|-----------|
| ![PCS Example 1](tutorial/images/example1.gif) | ![PCS Example 2](tutorial/images/example2.gif) | ![PCS Example 3](tutorial/images/example3.gif) |

### Tripod with Bending Sensors
![Tripod Sensor](tutorial/images/tripod_sensor.png)

## Getting Started

Follow the tutorial: [Cosserat Tutorial](tutorial/text/cosserat_tutorial.md)

## Publications

1. **Discrete cosserat approach for soft robot dynamics: A new piece-wise constant strain model with torsion and shears**  
   [IEEE Link](https://ieeexplore.ieee.org/document/7759808)

2. **Coupling numerical deformable models in global and reduced coordinates for the simulation of the direct and the inverse kinematics of Soft Robots**  
   [IEEE Link](https://ieeexplore.ieee.org/abstract/document/9362217)  
   [PDF Download](https://hal.archives-ouvertes.fr/hal-03192168/document)

## Contributing

We welcome contributions! Please see our contribution guidelines and feel free to open issues or pull requests.
