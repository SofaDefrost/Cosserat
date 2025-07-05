#!/usr/bin/env python3
"""
Simple test to verify the codebase reorganization works correctly.

This script tests that:
1. The new Python package can be imported
2. Key classes and functions are accessible
3. Backward compatibility is maintained
"""

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'python'))

def test_basic_imports():
    """Test that basic imports work."""
    print("Testing basic imports...")

    try:
        from cosserat import (BeamGeometryParameters, BeamPhysicsParameters,
                              CosseratBase, CosseratGeometry, Parameters,
                              addHeader, addSolverNode, addVisual)
        print("✓ Basic imports successful")
    except ImportError as e:
        print(f"✗ Basic import failed: {e}")
        return False

    return True

def test_specific_imports():
    """Test that specific module imports work."""
    print("Testing specific imports...")

    try:
        from cosserat.beam import CosseratBase
        from cosserat.geometry import (CosseratGeometry,
                                       calculate_beam_parameters)
        from cosserat.params import Parameters
        from cosserat.utils import addEdgeCollision
        print("✓ Specific imports successful")
    except ImportError as e:
        print(f"✗ Specific import failed: {e}")
        return False

    return True

def test_geometry_class():
    """Test that CosseratGeometry class works and has backward compatibility."""
    print("Testing CosseratGeometry class...")

    try:
        from cosserat.geometry import CosseratGeometry
        from cosserat.params import BeamGeometryParameters

        # Create test parameters
        params = BeamGeometryParameters(
            beam_length=30.0,
            nb_section=6,
            nb_frames=12
        )

        # Create geometry
        geometry = CosseratGeometry(params)

        # Test new property names
        assert hasattr(geometry, 'cable_positions'), "Missing cable_positions property"
        assert hasattr(geometry, 'section_lengths'), "Missing section_lengths property"
        assert hasattr(geometry, 'frames'), "Missing frames property"

        # Test backward compatibility properties
        assert hasattr(geometry, 'cable_positionF'), "Missing backward compatibility cable_positionF"
        assert hasattr(geometry, 'sectionsLengthList'), "Missing backward compatibility sectionsLengthList"
        assert hasattr(geometry, 'framesF'), "Missing backward compatibility framesF"

        # Test that they return the same data
        assert geometry.cable_positions == geometry.cable_positionF, "Compatibility property mismatch"
        assert geometry.section_lengths == geometry.sectionsLengthList, "Compatibility property mismatch"

        print("✓ CosseratGeometry class working correctly")
    except Exception as e:
        print(f"✗ CosseratGeometry test failed: {e}")
        return False

    return True

def main():
    """Run all tests."""
    print("=== Testing Cosserat Plugin Reorganization ===")
    print()

    tests = [
        test_basic_imports,
        test_specific_imports,
        test_geometry_class
    ]

    passed = 0
    total = len(tests)

    for test in tests:
        if test():
            passed += 1
        print()

    print(f"=== Results: {passed}/{total} tests passed ===")

    if passed == total:
        print("🎉 All tests passed! Reorganization successful.")
        return 0
    else:
        print("❌ Some tests failed. Check the reorganization.")
        return 1

if __name__ == "__main__":
    sys.exit(main())

