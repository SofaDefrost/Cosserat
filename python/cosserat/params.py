"""
Cosserat Beam Parameters Module
===============================

This module defines parameter classes for configuring Cosserat beam simulations.
Cosserat beam theory is an extension of classical beam theory that accounts for
micro-rotations and is particularly useful for modeling slender structures with
complex behaviors such as medical instruments, cables, and soft robotics components.

The parameters are organized into several dataclasses:
- BeamGeometryParameters: Defines the beam's physical dimensions and discretization
- BeamPhysicsBaseParameters: Defines material properties and cross-section
- BeamPhysicsParametersNoInertia: Simplified physics model without inertia
- BeamPhysicsParametersWithInertia: Advanced physics model with inertia terms
- SimulationParameters: Controls the simulation solver behavior
- VisualParameters: Controls the visual representation of the beam
- ContactParameters: Controls collision detection and response
- Parameters: Aggregates all parameter types into a single configuration

Default Values & Their Justification
-----------------------------------

BeamGeometryParameters:
- beam_length = 1.0m: Standard unit length for ease of scaling
- nb_section = 5: Balances computational cost with accuracy for most applications
- nb_frames = 30: Provides smoother visualization than the minimum required frames
- build_collision_model = 0: Collision detection disabled by default for performance

BeamPhysicsBaseParameters:
- young_modulus = 1.205e8 Pa: Approximates the elasticity of a moderately stiff plastic
- poisson_ratio = 0.4: Typical value for many common materials
- beam_mass = 1.0 kg: Unit mass for simple scaling
- beam_radius = 0.01m: 1cm radius, appropriate for visualization at meter scale
- beam_shape = 'circular': Most common and computationally efficient cross-section

BeamPhysicsParametersWithInertia (additional parameters):
- GI = 1.5708 N·m²: Torsional rigidity for a typical 1cm radius beam
- GA = 3.1416e4 N: Shear stiffness derived from default material properties
- EI = 0.7854 N·m²: Bending stiffness derived from default material properties
- EA = 3.1416e4 N: Axial stiffness derived from default material properties

SimulationParameters:
- rayleigh_stiffness = 0.2: Moderate stiffness damping for numerical stability
- rayleigh_mass = 0.1: Light mass damping to preserve dynamic behaviors
- firstOrder = False: Second-order integration for better accuracy

ContactParameters:
- alarmDistance = 0.05m: Detection range of 5cm for efficient collision detection
- contactDistance = 0.01m: 1cm contact response threshold
- tolerance = 1.0e-8: Precision level appropriate for meter-scale simulations
- maxIterations = 100: Balance between accuracy and performance for contact resolution

Usage Examples
-------------

Example 1: Creating default parameters
```python
from params import Parameters

# Create default parameters
default_params = Parameters()

# Validate all parameters
default_params.validate()
```

Example 2: Customizing beam properties
```python
from params import Parameters, BeamPhysicsParametersNoInertia, BeamGeometryParameters

# Create custom physics parameters for a softer material
physics_params = BeamPhysicsParametersNoInertia(
    young_modulus=5.0e6,  # 5 MPa, like a soft rubber
    poisson_ratio=0.45,   # Higher value for more rubber-like behavior
    beam_mass=0.5,        # Lighter 0.5kg beam
    beam_radius=0.005     # 5mm radius
)

# Create custom geometry with higher resolution
geo_params = BeamGeometryParameters(
    beam_length=2.0,      # 2m long beam
    nb_section=10,        # More sections for higher accuracy
    nb_frames=50          # More frames for smoother visualization
)

# Create parameters with custom physics and geometry
custom_params = Parameters(
    beam_physics_params=physics_params,
    beam_geo_params=geo_params
)
```

Example 3: Working with inertia
```python
from params import Parameters, BeamPhysicsParametersWithInertia

# Create parameters with inertia for more dynamic simulations
inertia_params = Parameters(
    beam_physics_params=BeamPhysicsParametersWithInertia(
        # Custom inertia parameters if needed
        GI=2.0,
        EI=1.0
    )
)
```

Example 4: Configuring contact parameters
```python
from params import Parameters, ContactParameters

# Setup parameters with friction contact
contact_params = ContactParameters(
    responseParams="mu=0.3",  # Friction coefficient of 0.3
    alarmDistance=0.1,        # 10cm alarm distance
    contactDistance=0.02      # 2cm contact distance
)

# Create parameters with custom contact handling
params_with_contact = Parameters(
    contact_params=contact_params
)
```
"""

# @todo use this dataclass to create the cosserat object

from dataclasses import dataclass, field
from typing import List, Literal, Optional, Type, Union, cast


#
@dataclass(frozen=True)
class BeamGeometryParameters:
    """
    Cosserat Beam Geometry parameters.

    Parameters:
        beam_length: Length of the beam in meters.
        nb_section: Number of sections along the beam length.
        nb_frames: Number of frames along the beam.
        build_collision_model: Flag to determine if collision model should be built (0: no, 1: yes).
    """

    beam_length: float = 1.0  # beam length in m
    nb_section: int = 5  # number of sections, sections along the beam length
    nb_frames: int = 30  # number of frames along the beam
    build_collision_model: int = 0

    def validate(self) -> None:
        """Validate the beam geometry parameters."""
        if self.beam_length <= 0:
            raise ValueError(f"Beam length must be positive, got {self.beam_length}")
        if self.nb_section <= 0:
            raise ValueError(f"Number of sections must be positive, got {self.nb_section}")
        if self.nb_frames <= 0:
            raise ValueError(f"Number of frames must be positive, got {self.nb_frames}")
        if self.nb_frames < self.nb_section:
            raise ValueError(f"Number of frames ({self.nb_frames}) must be greater than or equal to number of sections ({self.nb_section})")

    def __str__(self) -> str:
        """Return a string representation of the beam geometry parameters."""
        return (f"BeamGeometryParameters(length={self.beam_length}m, "
                f"sections={self.nb_section}, frames={self.nb_frames})")

@dataclass(frozen=True)
class BeamPhysicsBaseParameters:
    """
    Base class for Cosserat Beam Physics parameters.

    Parameters:
        young_modulus: Young's modulus of the beam material (Pa).
        poisson_ratio: Poisson's ratio of the beam material (dimensionless).
        beam_mass: Mass of the beam (kg).
        beam_radius: Radius of the beam for circular cross-section (m).
        beam_length: Length of the beam along the X axis (m).
        beam_shape: Shape of the beam cross-section ('circular' or 'rectangular').
        length_Y: Length in Y direction for rectangular beam (m).
        length_Z: Length in Z direction for rectangular beam (m).
        useInertia: Flag to determine if inertia should be considered.
    """
    young_modulus: float = 1.205e8  # ~120 MPa: comparable to some plastics like PMMA
    poisson_ratio: float = 0.3      # Common value for many materials (0.3-0.4 for plastics)
    beam_mass: float = 1.0          # 1kg for simplicity and easy scaling
    beam_radius: float = 0.01       # 1cm radius, suitable for visualization at meter scale
    beam_length: float = 1.0        # 1m length along the X axis for standard unit reference
    beam_shape: Literal['circular', 'rectangular'] = 'circular'
    length_Y: float = 0.1  # length in Y direction for rectangular beam
    length_Z: float = 0.1  # length in Z direction for rectangular beam
    useInertia: bool = False

    def validate(self) -> None:
        """Validate the beam physics parameters."""
        if self.young_modulus <= 0:
            raise ValueError(f"Young's modulus must be positive, got {self.young_modulus}")
        if not (0 < self.poisson_ratio < 0.5):
            raise ValueError(f"Poisson's ratio must be between 0 and 0.5, got {self.poisson_ratio}")
        if self.beam_mass <= 0:
            raise ValueError(f"Beam mass must be positive, got {self.beam_mass}")
        if self.beam_radius <= 0:
            raise ValueError(f"Beam radius must be positive, got {self.beam_radius}")
        if self.beam_length <= 0:
            raise ValueError(f"Beam length must be positive, got {self.beam_length}")
        if self.beam_shape not in ['circular', 'rectangular']:
            raise ValueError(f"Beam shape must be either 'circular' or 'rectangular', got '{self.beam_shape}'")
        if self.beam_shape == 'rectangular' and (self.length_Y <= 0 or self.length_Z <= 0):
            raise ValueError(f"For rectangular beam, length_Y and length_Z must be positive, got {self.length_Y} and {self.length_Z}")

    @property
    def cross_sectional_area(self) -> float:
        """Calculate the cross-sectional area of the beam."""
        if self.beam_shape == 'circular':
            import math
            return math.pi * self.beam_radius ** 2
        else:  # rectangular
            return self.length_Y * self.length_Z

    @property
    def moment_of_inertia(self) -> float:
        """Calculate the moment of inertia of the beam cross-section."""
        if self.beam_shape == 'circular':
            import math
            return (math.pi * self.beam_radius ** 4) / 4
        else:  # rectangular
            return (self.length_Y * self.length_Z ** 3) / 12

    @property
    def shear_modulus(self) -> float:
        """Calculate the shear modulus based on Young's modulus and Poisson's ratio."""
        return self.young_modulus / (2 * (1 + self.poisson_ratio))

    def __str__(self) -> str:
        """Return a string representation of the beam physics parameters."""
        return (f"BeamPhysicsParameters(E={self.young_modulus:.2e}Pa, v={self.poisson_ratio}, "
                f"m={self.beam_mass}kg, shape={self.beam_shape}, length={self.beam_length}m)")

@dataclass(frozen=True)
class BeamPhysicsParametersNoInertia(BeamPhysicsBaseParameters):
    """
    Parameters for a Cosserat Beam without inertia.

    This class inherits all parameters from BeamPhysicsBaseParameters
    and sets useInertia to False by default.
    """
    pass

@dataclass(frozen=True)
class BeamPhysicsParametersWithInertia(BeamPhysicsBaseParameters):
    """
    Parameters for a Cosserat Beam with inertia.

    Parameters:
        GI: Torsional rigidity (N·m²).
        GA: Shear stiffness (N).
        EI: Bending stiffness (N·m²).
        EA: Axial stiffness (N).

    This class inherits all parameters from BeamPhysicsBaseParameters
    and additionally includes inertia-related parameters.
    """
    GI: float = 1.5708    # π/2: Torsional rigidity for 1cm radius beam (G*J where J=πr⁴/2)
    GA: float = 3.1416e4  # πr²G: Shear stiffness (G=E/2(1+ν) for default material)
    EI: float = 0.7854    # π/4: Bending stiffness for 1cm radius beam (E*I where I=πr⁴/4)
    EA: float = 3.1416e4  # πr²E: Axial stiffness for default material and radius
    useInertia: bool = True

    def validate(self) -> None:
        """Validate the beam physics parameters including inertia parameters."""
        super().validate()
        if self.GI <= 0:
            raise ValueError(f"GI (torsional rigidity) must be positive, got {self.GI}")
        if self.GA <= 0:
            raise ValueError(f"GA (shear stiffness) must be positive, got {self.GA}")
        if self.EI <= 0:
            raise ValueError(f"EI (bending stiffness) must be positive, got {self.EI}")
        if self.EA <= 0:
            raise ValueError(f"EA (axial stiffness) must be positive, got {self.EA}")

    def __str__(self) -> str:
        """Return a string representation of the beam physics parameters with inertia."""
        base_str = super().__str__()
        return base_str[:-1] + f", GI={self.GI:.2e}, GA={self.GA:.2e}, EI={self.EI:.2e}, EA={self.EA:.2e})"

@dataclass(frozen=True)
class SimulationParameters:
    """
    Simulation parameters for the Cosserat Beam simulation.

    Parameters:
        rayleigh_stiffness: Rayleigh damping coefficient for stiffness.
        rayleigh_mass: Rayleigh damping coefficient for mass.
        firstOrder: Flag to use first-order integration scheme instead of second-order.
    """
    rayleigh_stiffness: float = 0.2  # Moderate stiffness-proportional damping for stability
    rayleigh_mass: float = 0.1       # Light mass-proportional damping to control low-frequency motion
    firstOrder: bool = False         # False = second-order integration (better accuracy)

    def validate(self) -> None:
        """Validate the simulation parameters."""
        if self.rayleigh_stiffness < 0:
            raise ValueError(f"Rayleigh stiffness must be non-negative, got {self.rayleigh_stiffness}")
        if self.rayleigh_mass < 0:
            raise ValueError(f"Rayleigh mass must be non-negative, got {self.rayleigh_mass}")

    def __str__(self) -> str:
        """Return a string representation of the simulation parameters."""
        return (f"SimulationParameters(rayleigh_stiffness={self.rayleigh_stiffness}, "
                f"rayleigh_mass={self.rayleigh_mass}, firstOrder={self.firstOrder})")

@dataclass(frozen=True)
class VisualParameters:
    """
    Visual parameters for the Cosserat Beam visualization.

    Parameters:
        showObject: Flag to determine if object should be shown (0: no, 1: yes).
        show_object_scale: Scale factor for visualization.
        show_object_color: RGBA color for visualization (values between 0.0 and 1.0).
    """

    showObject: int = 1
    show_object_scale: float = 1.0
    show_object_color: List[float] = field(default_factory=lambda: [1.0, 0.0, 0.0, 1.0])
    def validate(self) -> None:
        """Validate the visual parameters."""
        if len(self.show_object_color) != 4:
            raise ValueError(f"Color must have four components (RGBA), got {len(self.show_object_color)}")
        if not all(0.0 <= x <= 1.0 for x in self.show_object_color):
            raise ValueError(f"Color components must be in range [0, 1], got {self.show_object_color}")
        if self.show_object_scale <= 0:
            raise ValueError(f"Show object scale must be positive, got {self.show_object_scale}")

    def __str__(self) -> str:
        """Return a string representation of the visual parameters."""
        return (f"VisualParameters(showObject={self.showObject}, scale={self.show_object_scale}, "
                f"color={self.show_object_color})")


@dataclass(frozen=True)
class ContactParameters:
    """
    Contact parameters for the Cosserat Beam simulation.

    Parameters:
        responseParams: Parameters for the contact response (e.g., "mu=0.0" for friction coefficient).
        response: Type of contact constraint to use.
        alarmDistance: Distance at which collision detection is triggered.
        contactDistance: Distance at which contact response is activated.
        isMultithreading: Flag to use multithreading for collision detection.
        tolerance: Tolerance value for contact resolution.
        maxIterations: Maximum number of iterations for contact resolution.
        epsilon: Epsilon value for numerical stability in contact calculations.
    """

    responseParams: str = "mu=0.0"
    response: str = "FrictionContactConstraint"
    alarmDistance: float = 0.05     # 5cm detection range for efficient broad-phase collision detection
    contactDistance: float = 0.01   # 1cm contact response threshold (objects closer than this will generate contact forces)
    isMultithreading: bool = False  # Single-threaded collision detection by default for deterministic behavior
    tolerance: float = 1.0e-8       # Convergence criterion for constraint solvers
    maxIterations: int = 100        # Maximum solver iterations for contact resolution
    epsilon: float = 1.0e-6         # Small value for numerical stability in calculations
    def validate(self) -> None:
        """Validate the contact parameters."""
        if self.alarmDistance <= 0:
            raise ValueError(f"Alarm distance must be positive, got {self.alarmDistance}")
        if self.contactDistance <= 0:
            raise ValueError(f"Contact distance must be positive, got {self.contactDistance}")
        if self.alarmDistance < self.contactDistance:
            raise ValueError(f"Alarm distance ({self.alarmDistance}) must be greater than or equal to contact distance ({self.contactDistance})")
        if self.tolerance <= 0:
            raise ValueError(f"Tolerance must be positive, got {self.tolerance}")
        if self.maxIterations <= 0:
            raise ValueError(f"Maximum iterations must be positive, got {self.maxIterations}")
        if self.epsilon <= 0:
            raise ValueError(f"Epsilon must be positive, got {self.epsilon}")

    def __str__(self) -> str:
        """Return a string representation of the contact parameters."""
        return (f"ContactParameters(response={self.response}, responseParams={self.responseParams}, "
                f"alarmDist={self.alarmDistance}, contactDist={self.contactDistance}, "
                f"multithreading={self.isMultithreading})")

@dataclass(frozen=True)
class Parameters:
    """
    Comprehensive parameters for the Cosserat Beam simulation.

    This class aggregates all parameter sets needed for a complete Cosserat Beam
    simulation, providing a single entry point for configuration and validation.

    Parameters:
        beam_physics_params: Physics parameters for the beam material properties,
            including Young's modulus, Poisson's ratio, mass, and cross-section properties.
        simu_params: Simulation parameters controlling the numerical solver,
            including Rayleigh damping coefficients and integration order.
        contact_params: Contact parameters for collision handling,
            including alarm distance, contact distance, and solver settings.
        beam_geo_params: Geometry parameters defining beam dimensions and discretization,
            including length, number of sections, and number of frames.
        visual_params: Visual parameters for rendering the beam,
            including visibility flags, scale, and color.
    """

    beam_physics_params: BeamPhysicsBaseParameters = field(default_factory=BeamPhysicsParametersNoInertia)
    simu_params: SimulationParameters = field(
        default_factory=SimulationParameters
    )
    contact_params: ContactParameters = field(
        default_factory=ContactParameters
    )
    beam_geo_params: BeamGeometryParameters = field(
        default_factory=BeamGeometryParameters
    )
    visual_params: VisualParameters = field(default_factory=VisualParameters)

    def validate(self) -> None:
        """
        Validate all parameter sets in this parameters object.

        This method calls the validate method on each constituent parameter object.
        If any validation fails, a ValueError will be raised with appropriate message.

        Returns:
            None

        Raises:
            ValueError: If any parameter validation fails.
        """
        self.beam_physics_params.validate()
        self.simu_params.validate()
        self.contact_params.validate()
        self.beam_geo_params.validate()
        self.visual_params.validate()

    def __str__(self) -> str:
        """Return a comprehensive string representation of all parameters."""
        return (f"Parameters(\n"
                f"  Physics: {self.beam_physics_params}\n"
                f"  Geometry: {self.beam_geo_params}\n"
                f"  Simulation: {self.simu_params}\n"
                f"  Contact: {self.contact_params}\n"
                f"  Visual: {self.visual_params}\n"
                f")")
