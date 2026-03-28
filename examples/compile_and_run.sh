#!/bin/bash

# Script pour compiler et exécuter l'exemple d'optimisation de trajectoire

echo "=== Compilation de simple_trajectory_optimization ==="

# Chemins
EIGEN_PATH="/opt/homebrew/Cellar/eigen@3/3.4.1/include"
EIGEN3_PATH="/opt/homebrew/Cellar/eigen@3/3.4.1/include/eigen3"
SRC_PATH="../src"

# Compiler
clang++ -std=c++20 -O2 \
    -I${SRC_PATH} \
    -I${EIGEN_PATH} \
    -I${EIGEN3_PATH} \
    simple_trajectory_optimization.cpp \
    -o simple_trajectory_optimization

if [ $? -eq 0 ]; then
    echo "✓ Compilation réussie!"
    echo ""
    echo "=== Exécution de l'exemple ==="
    echo ""
    ./simple_trajectory_optimization
else
    echo "✗ Erreur de compilation"
    exit 1
fi
