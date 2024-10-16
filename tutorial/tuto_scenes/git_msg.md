Refactor event key handling and add new cable-driven Cosserat beam scene

- Updated key event handling in geo_cable_driven_cosserat_beam.py to use string-based key comparison for + and - keys.
- Fixed minor formatting issue by adding missing comma in RestShapeSpringsForceField.
- Added new scene file geo_cosserat_cable_driven_cosserat_beam.py for simulating cable-driven Cosserat beams, with detailed parameter definitions for beam geometry and physics.
- Introduced PullingCosseratCable controller in pulling_cosserat_cable.py to handle cable pulling mechanics.
- Replaced BeamPhysicsParameters with BeamPhysicsParametersNoInertia in tuto_5.py for more specific physical simulations.
