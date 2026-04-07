# Cosserat Plugin Tutorial Guide

This guide explains the progression through the Cosserat plugin tutorials and demonstrates the different API levels available.

## Tutorial Progression

### 🌟 Tutorial 01: Basic Cosserat Beam
**File:** `getting_started/tutorial_01_basic_beam.py`

**What you'll learn:**
- How to use `BeamGeometryParameters` to define beam dimensions
- How to create a `CosseratGeometry` object that automatically calculates all geometric properties
- Basic beam creation with manual geometry setup functions
- Clean, modular code structure

**API Level:** Medium-level (manual functions + CosseratGeometry)

**Key Code:**
```python
# Define beam parameters
beam_geometry_params = BeamGeometryParameters(
    beam_length=30.0,
    nb_section=3,
    nb_frames=4
)

# Automatic geometry calculation
beam_geometry = CosseratGeometry(beam_geometry_params)

# Use geometry in beam creation
_add_cosserat_state(root_node, beam_geometry, custom_bending_states)
```

### 🚀 Tutorial 02: Cosserat Beam with Forces
**File:** `getting_started/tutorial_02_with_forces.py`

**What you'll learn:**
- Adding dynamic simulation with gravity and applied forces
- Configuring solvers for dynamic behavior
- Mass distribution across beam frames
- Force application at specific beam locations

**API Level:** Medium-level (builds on Tutorial 01)

**Key Code:**
```python
# Dynamic scene setup
root_node.gravity = [0, -9.81, 0]
root_node.addObject("EulerImplicitSolver")
root_node.addObject("SparseLDLSolver")

# Beam with mass
frame_node = _add_cosserat_frame(base_node, bending_node, beam_geometry, beam_mass=5.0)

# Apply forces
frame_node.addObject('ConstantForceField', indices=[tip_index], forces=[force_vector])
```

### 🎮 Tutorial 03: Interactive Cosserat Beam
**File:** `getting_started/tutorial_03_interaction.py`

**What you'll learn:**
- Using the highest-level API: `CosseratBase` prefab
- Complete beam setup with a single class
- Spring attachments and interactive forces
- Physics parameter configuration

**API Level:** High-level (CosseratBase prefab)

**Key Code:**
```python
# Complete beam setup in one line!
beam = solver_node.addChild(CosseratBase(parent=solver_node, params=Params))

# Everything is automatically created:
# - Rigid base
# - Cosserat coordinates  
# - Frame mappings
# - Geometry calculations
```

## API Comparison

### Three Levels of API

| Level | Complexity | Control | Use Case |
|-------|------------|---------|----------|
| **High-level**<br>`CosseratBase` | Lowest | Least | Quick prototyping, standard beams |
| **Medium-level**<br>`CosseratGeometry` + functions | Medium | Medium | Custom beam setups, learning |
| **Low-level**<br>Manual calculations | Highest | Most | Advanced customization, research |

### When to Use Each API

#### Use CosseratBase (Tutorial 03) when:
- ✅ You want a complete beam quickly
- ✅ Standard physics parameters work for you
- ✅ You're prototyping or learning
- ✅ You need collision detection (built-in)

#### Use CosseratGeometry + functions (Tutorials 01-02) when:
- ✅ You want to understand the beam construction process
- ✅ You need custom force fields or solvers
- ✅ You're building educational content
- ✅ You want modular, reusable code

#### Use manual calculations when:
- ✅ You need complete control over every parameter
- ✅ You're doing research with non-standard setups
- ✅ You're extending the plugin with new features

## Code Evolution Showcase

### Old Manual Approach (Before Reorganization)
```python
# Manual geometry calculations - error-prone and verbose
nb_sections = 6
beam_length = 30
length_s = beam_length / float(nb_sections)
bending_states = []
list_sections_length = []
temp = 0.0
section_curv_abs = [0.0]

for i in range(nb_sections):
    bending_states.append([0, 0.0, 0.2])
    list_sections_length.append((((i + 1) * length_s) - i * length_s))
    temp += list_sections_length[i]
    section_curv_abs.append(temp)

# ... more manual calculations for frames ...
```

### New CosseratGeometry Approach (Tutorials 01-02)
```python
# Clean, automatic geometry - no manual calculations!
beam_geometry_params = BeamGeometryParameters(
    beam_length=30.0,
    nb_section=6,
    nb_frames=32
)
beam_geometry = CosseratGeometry(beam_geometry_params)

# All geometry data automatically available:
# - beam_geometry.section_lengths
# - beam_geometry.frames
# - beam_geometry.curv_abs_sections
# - beam_geometry.curv_abs_frames
```

### Highest Level CosseratBase Approach (Tutorial 03)
```python
# Ultimate simplicity - everything in one prefab!
beam = solver_node.addChild(CosseratBase(parent=solver_node, params=Params))
# Done! Complete beam with physics, visualization, and interaction ready.
```

## Migration Benefits

### Before Reorganization:
- ❌ Manual geometry calculations prone to errors
- ❌ Code scattered across multiple directories
- ❌ Inconsistent import patterns
- ❌ Hard to find and reuse functionality
- ❌ No clear learning progression

### After Reorganization:
- ✅ Automatic geometry calculations with `CosseratGeometry`
- ✅ Clean Python package structure
- ✅ Consistent, simple imports: `from cosserat import ...`
- ✅ Clear API progression from simple to complex
- ✅ Comprehensive documentation and examples
- ✅ Backward compatibility maintained

## Next Steps

1. **Start with Tutorial 01** to understand the basic concepts
2. **Progress to Tutorial 02** to learn about forces and dynamics
3. **Try Tutorial 03** to see the power of the high-level API
4. **Explore the examples/** directory for more complex scenarios
5. **Read the API documentation** in `tutorials/documentation/`

## Getting Help

- Check the `tutorials/README.md` for structure overview
- Look at `examples/` for more complex use cases
- Refer to `tutorials/documentation/` for detailed API docs
- All tutorials include extensive comments explaining the concepts

Happy coding with Cosserat beams! 🎉

