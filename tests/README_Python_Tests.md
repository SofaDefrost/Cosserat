# Cosserat Python Binding Tests

This directory contains Python unit tests for the Cosserat plugin's Python bindings. The tests verify that the C++ components are properly exposed to Python via pybind11.

## Test Structure

The Python tests are organized into several test classes:

### TestPointsManager
Tests the `PointsManager` class functionality:
- Object creation and initialization
- Adding points to state (`addNewPointToState()`)
- Removing points from state (`removeLastPointfromState()`)
- Multiple point operations

### TestLieGroups
Tests various Lie group classes:
- `SO2` (2D rotations): identity, composition, inverse
- `SO3` (3D rotations): identity, matrix operations
- `SE3` (rigid body transformations): exponential/logarithm maps
- Group inverse operations

### TestBundleOperations
Tests Bundle (product manifold) functionality:
- Bundle identity elements
- Group composition operations
- Mixed Lie group bundles

### TestCosseratIntegration
Integration tests:
- Module import verification
- Available binding attributes
- Scene creation with Cosserat components

## Running the Tests

### Method 1: Using the Test Runner Script

The easiest way to run the tests is using the provided test runner:

```bash
# Run all tests
./run_python_tests.py

# Run with verbose output
./run_python_tests.py --verbose

# Run only PointsManager tests
./run_python_tests.py --pattern PointsManager

# Check dependencies only
./run_python_tests.py --check-deps

# Run a specific test file
./run_python_tests.py --test-file unit/test_cosserat_bindings.py
```

### Method 2: Direct Python Execution

```bash
# Navigate to the Tests directory
cd Tests

# Run the test file directly
python3 unit/test_cosserat_bindings.py
```

### Method 3: Using CMake/CTest

If the project is built with CMake and tests are enabled:

```bash
# From the build directory
ctest -R CosseratPythonBindings

# Or run all tests
ctest
```

## Requirements

### Required Dependencies
- Python 3.6+
- NumPy
- unittest (part of Python standard library)

### Optional Dependencies (for full testing)
- SOFA Framework
- Cosserat plugin properly built and installed
- Python bindings compiled and available in PYTHONPATH

## Environment Setup

The tests require the Cosserat Python bindings to be available. This typically means:

1. **Build the Cosserat plugin** with Python bindings enabled
2. **Set PYTHONPATH** to include the location of the compiled bindings:
   ```bash
   export PYTHONPATH=/path/to/build/lib/python3/site-packages:$PYTHONPATH
   ```
3. **Ensure SOFA is properly installed** and accessible to Python

## Test Behavior

The tests are designed to be robust and handle missing dependencies gracefully:

- **Missing SOFA/Cosserat modules**: Tests will be skipped with appropriate messages
- **Missing optional components**: Individual tests skip themselves if components aren't available
- **Environment issues**: Clear error messages help diagnose setup problems

## Understanding Test Results

### Test Outcomes
- **PASS**: Test executed successfully
- **SKIP**: Test was skipped due to missing dependencies or components
- **FAIL**: Test found an actual problem with the bindings
- **ERROR**: Test couldn't run due to environment or setup issues

### Interpreting Skipped Tests
Skipped tests are normal and expected when:
- SOFA is not installed
- Cosserat bindings haven't been compiled
- Specific Lie group classes aren't available
- Optional components are missing

### Common Issues

1. **Import Errors**: Usually indicate PYTHONPATH isn't set correctly
2. **Missing Module Attributes**: May indicate incomplete binding compilation
3. **Scene Creation Failures**: Often related to missing SOFA plugins or components

## Extending the Tests

To add new tests:

1. **Add test methods** to existing test classes or create new test classes
2. **Follow naming convention**: Test methods should start with `test_`
3. **Handle missing components gracefully** using `skipTest()` for missing dependencies
4. **Document expected behavior** in test docstrings
5. **Update this README** if adding new test categories

### Example Test Method

```python
def test_new_functionality(self):
    """Test description of what this verifies."""
    try:
        # Test implementation
        result = self.component.new_method()
        self.assertIsNotNone(result)
    except AttributeError:
        self.skipTest("new_method not available")
    except Exception as e:
        self.fail(f"new_method failed: {e}")
```

## Integration with CI/CD

These tests can be integrated into continuous integration pipelines:

1. **CMake Integration**: Tests are automatically added if Python is found
2. **Dependency Checking**: The test runner can verify environment setup
3. **Return Codes**: Proper exit codes for automation systems
4. **Detailed Output**: JSON/XML output can be added for CI systems

## Troubleshooting

### "Could not import Sofa.Cosserat"
This usually means:
- SOFA is not installed
- Cosserat plugin wasn't built with Python bindings
- PYTHONPATH doesn't include the binding location

### "No module named 'numpy'"
Install NumPy:
```bash
pip install numpy
```

### "PointsManager not available"
This indicates the PointsManager binding wasn't compiled or isn't exposed properly.

### Tests Pass but Skip Everything
This usually means the bindings exist but the required components (like specific Lie groups) aren't available. Check the binding compilation logs.

## Contributing

When contributing to these tests:

1. **Test both positive and negative cases**
2. **Ensure tests are independent** and can run in any order
3. **Add appropriate cleanup** in tearDown methods
4. **Document any special setup requirements**
5. **Test the tests** on systems both with and without SOFA/Cosserat

