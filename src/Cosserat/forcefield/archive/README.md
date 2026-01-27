# Archived Components

This directory contains components that have been deprecated or replaced by better alternatives.

## MyUniformVelocityDampingForceField

A modified version of SOFA's UniformVelocityDampingForceField.

### Status
- Archived
- Replaced by standard SOFA implementation
- Kept for reference

### Original Features
- Velocity-based damping
- Support for multiple DOF types
- Implicit/explicit integration options
- Selective DOF damping via indices

### Replacement
Use SOFA's standard UniformVelocityDampingForceField instead, which provides:
- Better maintained codebase
- More extensive testing
- Better integration with SOFA framework
- Regular updates and bug fixes

## Improvement Recommendations

For future force field implementations:

1. Code Organization
   - Maintain clear separation between base, standard, and specialized implementations
   - Use consistent naming conventions
   - Keep experimental features in dedicated branches until ready

2. Documentation
   - Maintain up-to-date API documentation
   - Include usage examples
   - Document performance characteristics
   - Clear status indicators for experimental features

3. Testing
   - Implement comprehensive unit tests
   - Add performance benchmarks
   - Include validation tests
   - Document test coverage

4. Code Quality
   - Follow SOFA coding standards
   - Regular code review process
   - Clear deprecation process
   - Performance optimization guidelines

