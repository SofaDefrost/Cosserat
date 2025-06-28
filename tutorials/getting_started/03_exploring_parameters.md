---
title: Basic Slide 03
date: 2025-06-17
---

# Comparison regarding the number of sections 

## First configuration
```python
beam_geometry_params = BeamGeometryParameters(
        beam_length=30.0,  # Same beam length
        nb_section=30,  # 30 sections for good physics resolution
        nb_frames=30,  # 30 frames for smooth visualization
    )
```

## second configuration

```python
beam_geometry_params = BeamGeometryParameters(
        beam_length=30.0,  # Same beam length
        nb_section=3,  # 3 sections for good physics resolution
        nb_frames=3,  # 3 frames for smooth visualization
    )
```
---

## The two beams are falling under the influence of gravity


--- 

## Let's show that this is not only a matter of visualisation

The number of sections is not just about making the beam look smoother. It has a direct impact on the accuracy of the physics simulation. Each section has its own set of physical properties, and the interactions between these sections determine the overall behavior of the beam. A higher number of sections allows for a more detailed and accurate representation of the beam's deformation, but it also increases the computational cost. It's a trade-off between accuracy and performance.

```python
beam_geometry_params = BeamGeometryParameters(
        beam_length=30.0,  # Same beam length
        nb_section=3,  # 3 sections for good physics resolution
        nb_frames=15,  # 3 frames for smooth visualization
    )
```
--- 

## Let's show that this is not only a matter of visualisation

```python
beam_geometry_params = BeamGeometryParameters(
        beam_length=30.0,  # Same beam length
        nb_section=3,  # 3 sections for good physics resolution
        nb_frames=30,  # 3 frames for smooth visualization
    )
```