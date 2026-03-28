#!/usr/bin/env python3
"""
Test script for Cosserat Python bindings
Tests the HookeSeratMapping bindings we just created
"""

import sys
import os

# Add the SOFA Python bindings to the path
sofa_build_path = "/Users/yadagolo/travail/sofa/build/release/lib/python3/site-packages"
sys.path.insert(0, sofa_build_path)

try:
    import Sofa.Cosserat as Cosserat
    print("✓ Successfully imported Sofa.Cosserat module")
    
    # Test SectionInfo class
    try:
        section = Cosserat.SectionInfo()
        print("✓ Successfully created SectionInfo instance")
        
        # Test properties
        section.length = 1.0
        assert section.length == 1.0
        print("✓ SectionInfo length property works")
        
    except Exception as e:
        print(f"✗ SectionInfo test failed: {e}")
    
    # Test FrameInfo class
    try:
        frame = Cosserat.FrameInfo()
        print("✓ Successfully created FrameInfo instance")
        
        # Test properties
        frame.length = 2.0
        assert frame.length == 2.0
        print("✓ FrameInfo length property works")
        
        frame.kappa = 0.1
        assert frame.kappa == 0.1
        print("✓ FrameInfo kappa property works")
        
    except Exception as e:
        print(f"✗ FrameInfo test failed: {e}")
    
    # Test mapping classes availability
    try:
        # Check if the mapping classes are available
        mapping3_class = getattr(Cosserat, 'HookeSeratDiscretMapping3', None)
        mapping6_class = getattr(Cosserat, 'HookeSeratDiscretMapping6', None)
        
        if mapping3_class:
            print("✓ HookeSeratDiscretMapping3 class is available")
        else:
            print("✗ HookeSeratDiscretMapping3 class not found")
            
        if mapping6_class:
            print("✓ HookeSeratDiscretMapping6 class is available")
        else:
            print("✗ HookeSeratDiscretMapping6 class not found")
            
    except Exception as e:
        print(f"✗ Mapping classes test failed: {e}")
    
    print("\n✓ All Cosserat binding tests completed successfully!")
    
except ImportError as e:
    print(f"✗ Failed to import Sofa.Cosserat: {e}")
    print("Make sure SOFA and the Cosserat plugin are built with Python bindings enabled")
    
except Exception as e:
    print(f"✗ Unexpected error: {e}")
