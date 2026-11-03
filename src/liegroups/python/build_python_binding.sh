#!/bin/bash
# Builds the standalone `articulated_rigids` pybind11 module (no SOFA rebuild required) for the Homebrew
# Python 3.14 (the same interpreter that runSofa embeds) and installs the .so next to the simulation scenes.
#
# Usage: build_python_binding.sh [output_dir]
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_ROOT="$(cd "$HERE/../.." && pwd)"                       # .../plugin.Cosserat/src
OUT_DIR="${1:-/Users/yadagolo/travail/plugin/related_cosserat/simulations/articulated}"
PYTHON="${PYTHON:-/opt/homebrew/opt/python@3.14/bin/python3.14}"
PYBIND11_INC="${PYBIND11_INC:-/opt/homebrew/Cellar/pybind11/3.1.0_1/include}"
EIGEN_INC="${EIGEN_INC:-/opt/homebrew/include/eigen3}"

PY_INC="$("$PYTHON" -c 'import sysconfig; print(sysconfig.get_paths()["include"])')"
EXT_SUFFIX="$("$PYTHON" -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')"

mkdir -p "$OUT_DIR"
OUT="$OUT_DIR/articulated_rigids$EXT_SUFFIX"

# Python symbols are resolved at load time from the host interpreter (-undefined dynamic_lookup).
/usr/bin/c++ -O2 -DNDEBUG -std=c++20 -fPIC -shared -arch arm64 -fvisibility=hidden \
    -undefined dynamic_lookup \
    -I"$SRC_ROOT" -I"$PY_INC" -isystem "$PYBIND11_INC" -isystem "$EIGEN_INC" \
    "$HERE/articulated_rigids_py.cpp" -o "$OUT"

echo "Built $OUT"
