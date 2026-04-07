#!/usr/bin/env python3
"""
Python test runner for Cosserat bindings.

This script provides a convenient way to run Python binding tests
with proper environment setup.
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path


def setup_environment():
    """Set up the environment for running tests."""
    # Get the current script directory
    script_dir = Path(__file__).parent.absolute()
    project_root = script_dir.parent
    
    # Add potential Python paths
    python_paths = [
        script_dir / "unit",
        project_root / "build" / "lib" / "python3" / "site-packages",
        project_root / "lib" / "python3" / "site-packages",
    ]
    
    # Update PYTHONPATH
    current_pythonpath = os.environ.get("PYTHONPATH", "")
    new_paths = [str(p) for p in python_paths if p.exists()]
    
    if current_pythonpath:
        new_paths.append(current_pythonpath)
    
    os.environ["PYTHONPATH"] = os.pathsep.join(new_paths)
    
    print(f"Updated PYTHONPATH: {os.environ['PYTHONPATH']}")
    
    # Set working directory to Tests directory
    os.chdir(script_dir)
    print(f"Working directory: {os.getcwd()}")


def find_python_executable():
    """Find the appropriate Python executable."""
    # Try different Python executables
    candidates = ["python3", "python"]
    
    for candidate in candidates:
        try:
            result = subprocess.run(
                [candidate, "--version"], 
                capture_output=True, 
                text=True, 
                check=True
            )
            print(f"Using Python: {candidate} - {result.stdout.strip()}")
            return candidate
        except (subprocess.CalledProcessError, FileNotFoundError):
            continue
    
    raise RuntimeError("No suitable Python executable found")


def run_tests(test_file=None, verbose=False, pattern=None):
    """Run the Python tests."""
    setup_environment()
    python_exe = find_python_executable()
    
    # Default test file
    if test_file is None:
        test_file = "unit/test_cosserat_bindings.py"
    
    test_path = Path(test_file)
    if not test_path.exists():
        raise FileNotFoundError(f"Test file not found: {test_path}")
    
    # Build command
    cmd = [python_exe, str(test_path)]
    
    if verbose:
        cmd.extend(["-v"])
    
    if pattern:
        cmd.extend(["-k", pattern])
    
    print(f"Running command: {' '.join(cmd)}")
    print("=" * 60)
    
    # Run the tests
    try:
        result = subprocess.run(cmd, check=False)
        return result.returncode
    except KeyboardInterrupt:
        print("\nTests interrupted by user")
        return 1
    except Exception as e:
        print(f"Error running tests: {e}")
        return 1


def check_dependencies():
    """Check if required dependencies are available."""
    python_exe = find_python_executable()
    
    required_modules = ["unittest"]
    recommended_modules = ["numpy"]
    optional_modules = ["Sofa", "Sofa.Core", "Sofa.Cosserat"]
    
    print("Checking dependencies...")
    print("-" * 40)
    
    # Check required modules
    for module in required_modules:
        try:
            subprocess.run(
                [python_exe, "-c", f"import {module}"], 
                check=True, 
                capture_output=True
            )
            print(f"✓ {module} - available")
        except subprocess.CalledProcessError:
            print(f"✗ {module} - MISSING (required)")
            return False
    
    # Check optional modules
    for module in optional_modules:
        try:
            subprocess.run(
                [python_exe, "-c", f"import {module}"], 
                check=True, 
                capture_output=True
            )
            print(f"✓ {module} - available")
        except subprocess.CalledProcessError:
            print(f"✗ {module} - missing (optional - tests will be skipped)")
    
    return True


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Run Cosserat Python binding tests",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                           # Run all tests
  %(prog)s --verbose                 # Run with verbose output
  %(prog)s --pattern PointsManager   # Run only PointsManager tests
  %(prog)s --check-deps              # Check dependencies only
  %(prog)s --test-file my_test.py    # Run specific test file
        """
    )
    
    parser.add_argument(
        "--test-file", 
        help="Specific test file to run (default: unit/test_cosserat_bindings.py)"
    )
    parser.add_argument(
        "--verbose", "-v", 
        action="store_true", 
        help="Verbose output"
    )
    parser.add_argument(
        "--pattern", "-k", 
        help="Run only tests matching pattern"
    )
    parser.add_argument(
        "--check-deps", 
        action="store_true", 
        help="Check dependencies and exit"
    )
    
    args = parser.parse_args()
    
    print("Cosserat Python Binding Test Runner")
    print("====================================")
    
    # Check dependencies
    if not check_dependencies():
        print("\nError: Required dependencies are missing")
        return 1
    
    if args.check_deps:
        print("\nDependency check completed")
        return 0
    
    print("\nRunning tests...")
    return run_tests(
        test_file=args.test_file,
        verbose=args.verbose,
        pattern=args.pattern
    )


if __name__ == "__main__":
    sys.exit(main())

