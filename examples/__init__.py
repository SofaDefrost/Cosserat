"""Examples compatibility module

This module provides backward compatibility imports for examples.
"""

# Legacy imports for backward compatibility
import sys
import os

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))

# Import the main cosserat module
try:
    from cosserat import *
except ImportError:
    # Fallback for when running from examples directory
    pass

