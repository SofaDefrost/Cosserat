# Rigid Implementation

## Overview
This directory contains the rigid body variant of the Cosserat beam force field, specialized for simulations involving rigid body mechanics.

## Files
- `BeamHookeLawForceFieldRigid.h`: Header file for the rigid implementation
- `BeamHookeLawForceFieldRigid.cpp`: Implementation file
- `BeamHookeLawForceFieldRigid.inl`: Template implementation details

## Features
- Specialized implementation for rigid body mechanics
- Optimized for rigid body constraints
- Built on base implementation with rigid-specific optimizations

## Usage
Use this implementation when dealing with rigid body mechanics or when the beam segments should be treated as rigid bodies.

## Limitations
- Only suitable for rigid body simulations
- May not be appropriate for highly flexible beam simulations
