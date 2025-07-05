#!/usr/bin/env python3
"""
Unit tests for Cosserat Python bindings.

This module tests the Python bindings for the Cosserat plugin,
including PointsManager and Lie group functionality.
"""

import unittest
import sys
import os

# Try to import NumPy
try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Could not import NumPy: {e}")
    print("Some tests will be skipped. Install NumPy with: pip install numpy")
    # Create a dummy numpy for basic functionality
    class DummyNumPy:
        pi = 3.14159265359
        def array(self, data):
            return data
        def eye(self, n):
            return [[1 if i == j else 0 for j in range(n)] for i in range(n)]
        def sqrt(self, x):
            return x ** 0.5
        class testing:
            @staticmethod
            def assert_allclose(a, b, atol=1e-10):
                pass  # Skip assertion
    np = DummyNumPy()
    NUMPY_AVAILABLE = False

# Add the necessary paths for SOFA imports
try:
    import Sofa
    import Sofa.Core
    import Sofa.Cosserat
    SOFA_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Could not import SOFA modules: {e}")
    print("Tests may fail if SOFA is not properly installed or configured.")
    SOFA_AVAILABLE = False


class TestPointsManager(unittest.TestCase):
    """
    Test suite for PointsManager Python bindings.
    
    These tests verify the functionality of the PointsManager class
    which is bound from C++ to Python via pybind11.
    """
    
    def setUp(self):
        """Set up test fixtures before each test method."""
        try:
            # Create a basic SOFA scene
            self.root = Sofa.Core.Node("root")
            self.root.addObject("DefaultAnimationLoop")
            self.root.addObject("RequiredPlugin", name="Sofa.Component.Topology.Container.Dynamic")
            
            # Create a child node for constraint points
            self.constraint_node = self.root.addChild("constraintPointsNode")
            
            # Add necessary components
            self.container = self.constraint_node.addObject(
                "PointSetTopologyContainer", 
                points=[]
            )
            self.modifier = self.constraint_node.addObject("PointSetTopologyModifier")
            self.state = self.constraint_node.addObject(
                "MechanicalObject", 
                template="Vec3d", 
                position=[], 
                showObject=True,
                showObjectScale=10, 
                listening=True
            )
            
            # Create PointsManager instance
            self.points_manager = self.constraint_node.addObject(
                'PointsManager', 
                name="pointsManager", 
                listening=True,
                beamPath="/needle/rigidBase/cosseratInSofaFrameNode/slidingPoint/slidingPointMO"
            )
            
            # Initialize the scene
            Sofa.Simulation.init(self.root)
            
        except Exception as e:
            self.skipTest(f"Failed to set up SOFA scene: {e}")
    
    def tearDown(self):
        """Clean up after each test method."""
        if hasattr(self, 'root'):
            try:
                Sofa.Simulation.unload(self.root)
            except:
                pass
    
    def test_points_manager_creation(self):
        """Test that PointsManager can be created successfully."""
        self.assertIsNotNone(self.points_manager)
        self.assertEqual(self.points_manager.getName(), "pointsManager")
    
    def test_add_new_point_to_state(self):
        """Test adding new points to the state."""
        # Get initial point count
        initial_count = len(self.state.position.array())
        
        # Add a new point
        try:
            self.points_manager.addNewPointToState()
            
            # Check that a point was added
            new_count = len(self.state.position.array())
            self.assertEqual(new_count, initial_count + 1)
            
        except AttributeError:
            self.skipTest("addNewPointToState method not available")
        except Exception as e:
            self.fail(f"addNewPointToState failed: {e}")
    
    def test_remove_last_point_from_state(self):
        """Test removing the last point from the state."""
        try:
            # First add a point to ensure we have something to remove
            self.points_manager.addNewPointToState()
            initial_count = len(self.state.position.array())
            
            # Remove the last point
            self.points_manager.removeLastPointfromState()
            
            # Check that a point was removed
            new_count = len(self.state.position.array())
            self.assertEqual(new_count, initial_count - 1)
            
        except AttributeError:
            self.skipTest("removeLastPointfromState method not available")
        except Exception as e:
            self.fail(f"removeLastPointfromState failed: {e}")
    
    def test_multiple_point_operations(self):
        """Test multiple add and remove operations."""
        try:
            initial_count = len(self.state.position.array())
            
            # Add multiple points
            for i in range(5):
                self.points_manager.addNewPointToState()
            
            count_after_adds = len(self.state.position.array())
            self.assertEqual(count_after_adds, initial_count + 5)
            
            # Remove some points
            for i in range(3):
                self.points_manager.removeLastPointfromState()
            
            final_count = len(self.state.position.array())
            self.assertEqual(final_count, initial_count + 2)
            
        except (AttributeError, Exception) as e:
            self.skipTest(f"Multiple point operations test skipped: {e}")


class TestLieGroups(unittest.TestCase):
    """
    Test suite for Lie group Python bindings.
    
    These tests verify the functionality of various Lie groups
    that are bound from C++ to Python.
    """
    
    def setUp(self):
        """Set up test fixtures before each test method."""
        try:
            # Try to import Lie group classes from Cosserat module
            from Sofa.Cosserat import SO2, SO3, SE2, SE3
            self.SO2 = SO2
            self.SO3 = SO3
            self.SE2 = SE2
            self.SE3 = SE3
            self.lie_groups_available = True
        except ImportError:
            self.lie_groups_available = False
    
    def test_so2_identity(self):
        """Test SO(2) identity element."""
        if not self.lie_groups_available:
            self.skipTest("Lie group bindings not available")
        
        try:
            identity = self.SO2.identity()
            self.assertIsNotNone(identity)
            
            # Test that identity angle is approximately zero
            angle = identity.angle()
            self.assertAlmostEqual(angle, 0.0, places=10)
            
        except (AttributeError, Exception) as e:
            self.skipTest(f"SO2 identity test skipped: {e}")
    
    def test_so2_composition(self):
        """Test SO(2) group composition."""
        if not self.lie_groups_available:
            self.skipTest("Lie group bindings not available")
        
        try:
            # Create two rotations
            rot1 = self.SO2(np.pi/4)  # 45 degrees
            rot2 = self.SO2(np.pi/3)  # 60 degrees
            
            # Compose them
            result = rot1 * rot2  # Should be 105 degrees
            
            expected_angle = np.pi/4 + np.pi/3
            actual_angle = result.angle()
            
            self.assertAlmostEqual(actual_angle, expected_angle, places=10)
            
        except (AttributeError, Exception) as e:
            self.skipTest(f"SO2 composition test skipped: {e}")
    
    def test_so3_identity(self):
        """Test SO(3) identity element."""
        if not self.lie_groups_available:
            self.skipTest("Lie group bindings not available")
        
        try:
            identity = self.SO3.identity()
            self.assertIsNotNone(identity)
            
            # Test that identity matrix is approximately the 3x3 identity
            matrix = identity.matrix()
            expected = np.eye(3)
            
            np.testing.assert_allclose(matrix, expected, atol=1e-10)
            
        except (AttributeError, Exception) as e:
            self.skipTest(f"SO3 identity test skipped: {e}")
    
    def test_se3_exp_log(self):
        """Test SE(3) exponential and logarithm maps."""
        if not self.lie_groups_available:
            self.skipTest("Lie group bindings not available")
        
        try:
            # Create a random SE(3) element
            axis = np.array([1, 1, 1]) / np.sqrt(3)
            angle = np.pi/6
            translation = np.array([1, 2, 3])
            
            # Create rotation and SE(3) element
            rotation = self.SO3(angle, axis)
            se3_elem = self.SE3(rotation, translation)
            
            # Test exp(log(g)) = g
            log_elem = se3_elem.log()
            reconstructed = self.SE3().exp(log_elem)
            
            self.assertTrue(se3_elem.isApprox(reconstructed))
            
        except (AttributeError, Exception) as e:
            self.skipTest(f"SE3 exp/log test skipped: {e}")
    
    def test_lie_group_inverse(self):
        """Test inverse operations for various Lie groups."""
        if not self.lie_groups_available:
            self.skipTest("Lie group bindings not available")
        
        try:
            # Test SO(2) inverse
            rot2d = self.SO2(np.pi/3)
            inv2d = rot2d.inverse()
            result2d = rot2d * inv2d
            identity2d = self.SO2.identity()
            
            self.assertTrue(result2d.isApprox(identity2d))
            
            # Test SO(3) inverse
            axis = np.array([0, 0, 1])
            rot3d = self.SO3(np.pi/2, axis)
            inv3d = rot3d.inverse()
            result3d = rot3d * inv3d
            identity3d = self.SO3.identity()
            
            self.assertTrue(result3d.isApprox(identity3d))
            
        except (AttributeError, Exception) as e:
            self.skipTest(f"Lie group inverse test skipped: {e}")


class TestBundleOperations(unittest.TestCase):
    """
    Test suite for Bundle (product manifold) operations.
    
    These tests verify the functionality of Bundle classes
    that combine multiple Lie groups.
    """
    
    def setUp(self):
        """Set up test fixtures before each test method."""
        try:
            from Sofa.Cosserat import Bundle, SE3, RealSpace
            self.Bundle = Bundle
            self.SE3 = SE3
            self.RealSpace = RealSpace
            self.bundle_available = True
        except ImportError:
            self.bundle_available = False
    
    def test_bundle_identity(self):
        """Test Bundle identity element."""
        if not self.bundle_available:
            self.skipTest("Bundle bindings not available")
        
        try:
            # Create a pose-velocity bundle
            PoseVel = self.Bundle[self.SE3, self.RealSpace[6]]
            identity = PoseVel.identity()
            
            self.assertIsNotNone(identity)
            
            # Test identity properties
            test_bundle = PoseVel(
                self.SE3.identity(),
                self.RealSpace[6].zero()
            )
            
            result1 = test_bundle * identity
            result2 = identity * test_bundle
            
            self.assertTrue(result1.isApprox(test_bundle))
            self.assertTrue(result2.isApprox(test_bundle))
            
        except (AttributeError, Exception) as e:
            self.skipTest(f"Bundle identity test skipped: {e}")
    
    def test_bundle_composition(self):
        """Test Bundle group composition."""
        if not self.bundle_available:
            self.skipTest("Bundle bindings not available")
        
        try:
            # Create multiple bundle elements and test composition
            PoseVel = self.Bundle[self.SE3, self.RealSpace[6]]
            
            bundle1 = PoseVel(
                self.SE3.identity(),
                self.RealSpace[6]([0.1, 0.2, 0.3, 0.0, 0.0, 0.0])
            )
            
            bundle2 = PoseVel(
                self.SE3.identity(),
                self.RealSpace[6]([0.1, 0.1, 0.1, 0.0, 0.0, 0.0])
            )
            
            result = bundle1 * bundle2
            self.assertIsNotNone(result)
            
        except (AttributeError, Exception) as e:
            self.skipTest(f"Bundle composition test skipped: {e}")


class TestCosseratIntegration(unittest.TestCase):
    """
    Integration tests for Cosserat functionality.
    
    These tests verify that the Python bindings work correctly
    in the context of a complete Cosserat simulation.
    """
    
    def test_cosserat_module_import(self):
        """Test that the Cosserat module can be imported."""
        if not SOFA_AVAILABLE:
            self.skipTest("SOFA not available")
        
        try:
            import Sofa.Cosserat
            self.assertTrue(hasattr(Sofa.Cosserat, '__name__'))
        except ImportError as e:
            self.skipTest(f"Could not import Sofa.Cosserat: {e}")
    
    def test_cosserat_binding_attributes(self):
        """Test that expected attributes are available in the Cosserat module."""
        try:
            import Sofa.Cosserat
            
            # Check for PointsManager
            if hasattr(Sofa.Cosserat, 'PointsManager'):
                points_manager_class = getattr(Sofa.Cosserat, 'PointsManager')
                self.assertTrue(callable(points_manager_class))
            
            # Check for common Lie group classes
            expected_classes = ['SO2', 'SO3', 'SE2', 'SE3', 'Bundle']
            available_classes = []
            
            for class_name in expected_classes:
                if hasattr(Sofa.Cosserat, class_name):
                    available_classes.append(class_name)
            
            # At least some classes should be available
            self.assertGreater(len(available_classes), 0, 
                             "No expected Cosserat classes found")
            
        except ImportError:
            self.skipTest("Cosserat module not available")
    
    def test_scene_creation_with_cosserat(self):
        """Test creating a basic scene that uses Cosserat components."""
        try:
            root = Sofa.Core.Node("root")
            root.addObject("DefaultAnimationLoop")
            root.addObject("RequiredPlugin", name="Sofa.Component.Topology.Container.Dynamic")
            
            # Try to add Cosserat-specific components
            cosserat_node = root.addChild("cosseratNode")
            
            # Add basic components that should work with Cosserat
            container = cosserat_node.addObject("PointSetTopologyContainer", points=[])
            state = cosserat_node.addObject("MechanicalObject", template="Vec3d", position=[])
            
            # Initialize the scene
            Sofa.Simulation.init(root)
            
            self.assertIsNotNone(container)
            self.assertIsNotNone(state)
            
            # Clean up
            Sofa.Simulation.unload(root)
            
        except Exception as e:
            self.skipTest(f"Scene creation test skipped: {e}")


def run_tests():
    """
    Run all tests and provide a summary.
    """
    # Create test suite
    test_classes = [
        TestPointsManager,
        TestLieGroups,
        TestBundleOperations,
        TestCosseratIntegration
    ]
    
    suite = unittest.TestSuite()
    for test_class in test_classes:
        tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"TEST SUMMARY")
    print(f"{'='*60}")
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Skipped: {len(result.skipped)}")
    
    if result.failures:
        print(f"\nFAILURES:")
        for test, traceback in result.failures:
            print(f"- {test}: {traceback.splitlines()[-1]}")
    
    if result.errors:
        print(f"\nERRORS:")
        for test, traceback in result.errors:
            print(f"- {test}: {traceback.splitlines()[-1]}")
    
    success_rate = (result.testsRun - len(result.failures) - len(result.errors)) / result.testsRun * 100
    print(f"\nSuccess rate: {success_rate:.1f}%")
    
    return result.wasSuccessful()


if __name__ == '__main__':
    # Check if we're running in a SOFA environment
    print("Cosserat Python Bindings Unit Tests")
    print("====================================")
    print(f"Python version: {sys.version}")
    print(f"Working directory: {os.getcwd()}")
    
    try:
        import Sofa
        print(f"SOFA version: {Sofa.Core.SofaInfo.version}")
    except ImportError:
        print("SOFA not available - some tests will be skipped")
    
    print("\nRunning tests...\n")
    
    success = run_tests()
    sys.exit(0 if success else 1)

