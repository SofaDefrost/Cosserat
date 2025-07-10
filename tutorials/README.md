# Cosserat Plugin Tutorials

This directory contains the reorganized tutorial structure for the Cosserat plugin.

## New Directory Structure

```
plugin.Cosserat/
├── python/                        # Main Python package
│   ├── cosserat/                  # Core Cosserat Python package
│   │   ├── __init__.py            # Main module exports
│   │   ├── beam.py                # CosseratBase class (main beam class)
│   │   ├── geometry.py            # CosseratGeometry and helper functions
│   │   ├── params.py              # Parameter classes
│   │   ├── utils.py               # Utility functions
│   │   ├── header.py              # Scene setup helpers
│   │   └── *.py                   # Other core modules
│   └── tests/                     # Unit tests
├── tutorials/                     # Tutorial content
│   ├── getting_started/           # Step-by-step beginner tutorials
│   │   ├── tutorial_01_basic_beam.py
│   │   ├── tutorial_02_with_forces.py
│   │   └── tutorial_03_interaction.py
│   ├── documentation/             # Tutorial documentation
│   │   ├── cosserat_tutorial.md
│   │   └── api_reference.md
│   └── assets/                    # Meshes, textures, data files
├── examples/                      # Clean examples directory
│   ├── basic/                     # Simple demonstration examples
│   ├── advanced/                  # Complex application examples
│   └── benchmarks/                # Performance and validation examples
└── docs/                          # Project documentation
```

## Migration from Old Structure

### What Changed

1. **Core utilities moved**: Files from `examples/python3/useful/` are now in `python/cosserat/`
2. **CosseratBase renamed**: `CosseratBase.py` is now `beam.py` in the `python/cosserat/` package
3. **Tutorials reorganized**: Tutorial files moved from `tutorial/tuto_scenes/` to `tutorials/getting_started/`
4. **Examples categorized**: Examples now organized by complexity level
5. **Proper Python package**: The codebase now follows Python packaging standards

### Import Changes

**Old imports:**
```python
from useful.geometry import CosseratGeometry
from useful.params import Parameters
from examples.python3.cosserat.CosseratBase import CosseratBase
```

**New imports:**
```python
from cosserat import CosseratGeometry, Parameters, CosseratBase
# or more specific:
from cosserat.geometry import CosseratGeometry
from cosserat.params import Parameters
from cosserat.beam import CosseratBase
```

## Getting Started

### For New Users

1. Start with `tutorials/getting_started/tutorial_01_basic_beam.py`
2. Progress through the numbered tutorials in order
3. Refer to `tutorials/documentation/` for detailed explanations

### For Existing Users

1. Update your imports as shown above
2. The `CosseratGeometry` class now has improved property names but maintains backward compatibility
3. All functionality is preserved, just better organized

## Tutorial Progression

### Getting Started Series

1. **tutorial_01_basic_beam.py**: Create a simple Cosserat beam
2. **tutorial_02_with_forces.py**: Add forces and constraints
3. **tutorial_03_interaction.py**: Interactive manipulation

### Examples by Category

- **basic/**: Simple, standalone examples for learning concepts
- **advanced/**: Complex scenarios with multiple components
- **benchmarks/**: Performance testing and validation scenarios

## Development

### Adding New Tutorials

1. Place beginner tutorials in `tutorials/getting_started/`
2. Use descriptive, numbered filenames
3. Add documentation in `tutorials/documentation/`
4. Include any required assets in `tutorials/assets/`

### Testing

- Unit tests are located in `python/tests/`
- Run tests to ensure the reorganization didn't break functionality

## Backward Compatibility

- The `CosseratGeometry` class includes compatibility properties for old property names
- Most existing code should work with minimal import changes
- Examples include path setup code for finding the new Python package

## Questions?

Refer to the documentation in `tutorials/documentation/` or check the examples for usage patterns.

