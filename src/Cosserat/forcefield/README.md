# Cosserat Force Field Components

This directory contains the implementation of force fields for the Cosserat beam model, organized for maintainability and clarity.

## Directory Structure

```
.
├── base/               # Base implementation
├── standard/          # Standard implementation
├── rigid/             # Rigid body variant
├── experimental/      # Components under development
├── archive/           # Deprecated components
├── docs/              # Documentation
│   ├── api/           # API documentation
│   ├── design/        # Design documents
│   └── implementation/# Implementation details
└── maintain.sh        # Maintenance script
```

## Components

### Core Components
- `BaseBeamHookeLawForceField`: Base implementation with Lie Group integration
- `BeamHookeLawForceField`: Standard implementation for general use
- `BeamHookeLawForceFieldRigid`: Specialized implementation for rigid bodies

### Experimental
- `CosseratInternalActuation`: Internal actuation component (under development)

### Archived
- Previous implementations and deprecated components

## Development

### Maintenance
Use the maintain.sh script for common maintenance tasks:
```bash
./maintain.sh check-docs  # Check documentation completeness
./maintain.sh check-code  # Check code structure
./maintain.sh clean       # Clean temporary files
```

### Contributing
1. Follow the established directory structure
2. Update documentation appropriately
3. Add tests for new features
4. Use experimental/ for work in progress

### Documentation
- See docs/README.md for documentation structure
- Keep API documentation up to date
- Document design decisions
- Maintain implementation details

## License
See individual component files for license information.
