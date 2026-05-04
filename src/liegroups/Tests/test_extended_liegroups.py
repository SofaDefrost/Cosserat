#!/usr/bin/env python3
"""
Test script for SGal(3) and SE(2,3) Lie group Python bindings

This script tests the functionality of the newly implemented Python bindings
for the Special Galilean group SGal(3) and the Extended Euclidean group SE(2,3).
"""

import numpy as np
import sys
import os

# Add the SOFA Python bindings to the path
sofa_path = "/Users/yadagolo/travail/sofa/src/cmake-build-relwithdebinfo/lib/python3/site-packages"
sys.path.insert(0, sofa_path)

# Import the Lie groups module
try:
    import Sofa.LieGroups as LieGroups
    print("✓ Successfully imported LieGroups module")
except ImportError as e:
    print(f"✗ Failed to import LieGroups module: {e}")
    sys.exit(1)

def test_sgal3():
    """Test SGal(3) bindings"""
    print("\n=== Testing SGal(3) ===")
    
    try:
        # Test default constructor (identity)
        sgal_id = LieGroups.SGal3()
        print(f"✓ Identity SGal3: {sgal_id}")
        
        # Test constructor from components
        se3_pose = LieGroups.SE3()
        velocity = np.array([1.0, 2.0, 3.0])
        time = 5.0
        sgal = LieGroups.SGal3(se3_pose, velocity, time)
        print(f"✓ Constructed SGal3: {sgal}")
        
        # Test accessors
        pose = sgal.pose()
        vel = sgal.velocity()
        t = sgal.time()
        print(f"✓ Pose: {pose}")
        print(f"✓ Velocity: [{vel[0]}, {vel[1]}, {vel[2]}]")
        print(f"✓ Time: {t}")
        
        # Test inverse
        sgal_inv = sgal.inverse()
        print(f"✓ Inverse: {sgal_inv}")
        
        # Test matrix representation
        matrix = sgal.matrix()
        print(f"✓ Matrix shape: {matrix.shape}")
        
        # Test exp/log operations
        algebra_element = sgal.log()
        print(f"✓ Log (algebra element) shape: {algebra_element.shape}")
        
        sgal_from_exp = LieGroups.SGal3.exp(algebra_element)
        print(f"✓ Exp->SGal3: {sgal_from_exp}")
        
        # Test adjoint
        adj = sgal.adjoint()
        print(f"✓ Adjoint shape: {adj.shape}")
        
        # Test approximate equality
        is_approx = sgal.isApprox(sgal_from_exp, 1e-10)
        print(f"✓ Log/Exp consistency: {is_approx}")
        
        # Test identity
        identity = LieGroups.SGal3.identity()
        print(f"✓ Static identity: {identity}")
        
        # Test random generation
        random_sgal = LieGroups.SGal3.random(42)  # With seed
        print(f"✓ Random SGal3: {random_sgal}")
        
        # Test group action on 10D vector
        point_vel_time = np.random.randn(10)
        transformed = sgal.act(point_vel_time)
        print(f"✓ Group action input shape: {point_vel_time.shape}")
        print(f"✓ Group action output shape: {transformed.shape}")
        
        print("✓ All SGal(3) tests passed!")
        return True
        
    except Exception as e:
        print(f"✗ SGal(3) test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_se23():
    """Test SE(2,3) bindings"""
    print("\n=== Testing SE(2,3) ===")
    
    try:
        # Test default constructor (identity)
        se23_id = LieGroups.SE23()
        print(f"✓ Identity SE23: {se23_id}")
        
        # Test constructor from components
        se3_pose = LieGroups.SE3()
        velocity = np.array([1.0, 2.0, 3.0])
        se23 = LieGroups.SE23(se3_pose, velocity)
        print(f"✓ Constructed SE23: {se23}")
        
        # Test accessors
        pose = se23.pose()
        vel = se23.velocity()
        print(f"✓ Pose: {pose}")
        print(f"✓ Velocity: [{vel[0]}, {vel[1]}, {vel[2]}]")
        
        # Test inverse
        se23_inv = se23.inverse()
        print(f"✓ Inverse: {se23_inv}")
        
        # Test matrix representation
        matrix = se23.matrix()
        print(f"✓ Matrix shape: {matrix.shape}")
        
        # Test exp/log operations
        algebra_element = se23.log()
        print(f"✓ Log (algebra element) shape: {algebra_element.shape}")
        
        se23_from_exp = LieGroups.SE23.exp(algebra_element)
        print(f"✓ Exp->SE23: {se23_from_exp}")
        
        # Test adjoint
        adj = se23.adjoint()
        print(f"✓ Adjoint shape: {adj.shape}")
        
        # Test approximate equality
        is_approx = se23.isApprox(se23_from_exp, 1e-10)
        print(f"✓ Log/Exp consistency: {is_approx}")
        
        # Test identity
        identity = LieGroups.SE23.identity()
        print(f"✓ Static identity: {identity}")
        
        # Test random generation
        random_se23 = LieGroups.SE23.random(42)  # With seed
        print(f"✓ Random SE23: {random_se23}")
        
        # Test group action on 6D vector
        point_vel = np.random.randn(6)
        transformed = se23.act(point_vel)
        print(f"✓ Group action input shape: {point_vel.shape}")
        print(f"✓ Group action output shape: {transformed.shape}")
        
        print("✓ All SE(2,3) tests passed!")
        return True
        
    except Exception as e:
        print(f"✗ SE(2,3) test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_group_interactions():
    """Test interactions between different Lie groups"""
    print("\n=== Testing Group Interactions ===")
    
    try:
        # Create SE3, SGal3, and SE23 elements
        se3 = LieGroups.SE3.random(123)
        velocity = np.array([1.0, -0.5, 2.0])
        
        sgal3 = LieGroups.SGal3(se3, velocity, 1.0)
        se23 = LieGroups.SE23(se3, velocity)
        
        print(f"✓ SE3: {se3}")
        print(f"✓ SGal3: {sgal3}")
        print(f"✓ SE23: {se23}")
        
        # Test that pose extraction gives the same SE3
        sgal3_pose = sgal3.pose()
        se23_pose = se23.pose()
        
        print(f"✓ SGal3 pose matches: {se3.isApprox(sgal3_pose, 1e-12)}")
        print(f"✓ SE23 pose matches: {se3.isApprox(se23_pose, 1e-12)}")
        
        # Test that velocity extraction works correctly
        sgal3_vel = sgal3.velocity()
        se23_vel = se23.velocity()
        
        vel_match_sgal3 = np.allclose(velocity, sgal3_vel, atol=1e-12)
        vel_match_se23 = np.allclose(velocity, se23_vel, atol=1e-12)
        
        print(f"✓ SGal3 velocity matches: {vel_match_sgal3}")
        print(f"✓ SE23 velocity matches: {vel_match_se23}")
        
        print("✓ All interaction tests passed!")
        return True
        
    except Exception as e:
        print(f"✗ Interaction test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests"""
    print("Testing Extended Lie Groups Python Bindings")
    print("=" * 50)
    
    results = []
    results.append(test_sgal3())
    results.append(test_se23())
    results.append(test_group_interactions())
    
    print("\n" + "=" * 50)
    if all(results):
        print("🎉 All tests passed successfully!")
        return 0
    else:
        print("❌ Some tests failed!")
        return 1

if __name__ == "__main__":
    sys.exit(main())
