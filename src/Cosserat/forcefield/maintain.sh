#!/bin/bash

# Script to help maintain the forcefield codebase organization

# Function to check documentation completeness
check_docs() {
    local missing=0
    
    # Check required documentation files
    for dir in docs/api docs/implementation docs/design; do
        if [ ! -d "$dir" ]; then
            echo "Missing directory: $dir"
            missing=1
        fi
    done

    # Check README files
    for dir in . docs experimental archive base standard rigid; do
        if [ ! -f "$dir/README.md" ]; then
            echo "Missing README.md in $dir"
            missing=1
        fi
    done

    if [ $missing -eq 0 ]; then
        echo "Documentation structure is complete."
    fi
}

# Function to check code organization
check_code_structure() {
    local missing=0
    
    # Check required code directories
    for dir in base standard rigid experimental archive; do
        if [ ! -d "$dir" ]; then
            echo "Missing directory: $dir"
            missing=1
        fi
    done

    # Check core implementation files
    for file in base/BaseBeamHookeLawForceField.{h,inl} standard/BeamHookeLawForceField.{cpp,h,inl} rigid/BeamHookeLawForceFieldRigid.{cpp,h,inl}; do
        if [ ! -f "$file" ]; then
            echo "Missing file: $file"
            missing=1
        fi
    done

    if [ $missing -eq 0 ]; then
        echo "Code structure is complete."
    fi
}

# Function to clean temporary files
clean_temp_files() {
    find . -name "*.swp" -delete
    find . -name "*~" -delete
    find . -name "*.bak" -delete
    find . -name ".DS_Store" -delete
    echo "Temporary files cleaned."
}

# Main menu
case "$1" in
    "check-docs")
        check_docs
        ;;
    "check-code")
        check_code_structure
        ;;
    "clean")
        clean_temp_files
        ;;
    *)
        echo "Usage: $0 {check-docs|check-code|clean}"
        echo "  check-docs  : Check documentation completeness"
        echo "  check-code  : Check code structure"
        echo "  clean       : Clean temporary files"
        ;;
esac
