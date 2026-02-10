#!/bin/bash
# Test script to validate conda environment setup for OpenGV on Windows (Git Bash)
# This mirrors the CI workflow environment validation
# Usage: bash scripts/test_conda_env_windows.sh

set -e  # Exit on error

# Change to project root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

echo "Testing OpenGV Conda Environment (Windows)"
echo "==========================================="

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda is not installed or not in PATH"
    exit 1
fi

# Check if environment exists
if ! conda env list | grep -q "^opengv "; then
    echo "ERROR: 'opengv' conda environment does not exist"
    echo "Please create it first: conda env create -f environment.yml"
    exit 1
fi

echo "[OK] Conda environment 'opengv' exists"

# Activate environment
eval "$(conda shell.bash hook)"
conda activate opengv

echo "[OK] Activated conda environment"

# Check Python
echo -n "Checking Python... "
PYTHON_VERSION=$(python --version 2>&1)
echo "$PYTHON_VERSION"
if ! python -c "import sys; assert sys.version_info >= (3, 7)"; then
    echo "ERROR: Python 3.7+ required"
    exit 1
fi
echo "[OK] Python version OK"

# Check required packages
echo "Checking required packages..."
MISSING=0

for pkg in numpy pybind11; do
    if python -c "import $pkg" 2>/dev/null; then
        echo "[OK] $pkg"
    else
        echo "[FAIL] $pkg MISSING"
        MISSING=1
    fi
done

if [ $MISSING -eq 1 ]; then
    echo "ERROR: Missing packages"
    exit 1
fi

# Check CMake
echo -n "Checking CMake... "
if command -v cmake &> /dev/null; then
    CMAKE_VERSION=$(cmake --version | head -n1)
    echo "$CMAKE_VERSION"
    echo "[OK] CMake found"
else
    echo "ERROR: CMake not found"
    exit 1
fi

# Check Eigen (via conda)
echo -n "Checking Eigen... "
if conda list eigen 2>/dev/null | grep -q eigen; then
    echo "[OK] Eigen found via conda"
elif [ -d "$CONDA_PREFIX/include/eigen3" ] || [ -d "$CONDA_PREFIX/include/eigen" ]; then
    echo "[OK] Eigen headers found"
elif [ -d "$CONDA_PREFIX/Library/include/eigen3" ]; then
    echo "[OK] Eigen headers found (Windows conda path)"
else
    echo "[WARN] Eigen headers not found in expected location (may still work)"
fi

# Check for python311.lib (critical for Windows builds)
echo -n "Checking Python lib file... "
pythonDir=$(python -c "import sys; print(sys.prefix)")
if [ -f "$pythonDir/libs/python311.lib" ]; then
    echo "[OK] python311.lib found"
elif ls "$pythonDir/libs/"python*.lib 1>/dev/null 2>&1; then
    echo "[OK] Python lib files found"
else
    echo "[WARN] Python lib files not found - may cause link errors"
fi

# Test basic build configuration
echo ""
echo "Testing CMake configuration..."
rm -rf build_test
mkdir -p build_test
cd build_test

if cmake .. -DBUILD_PYTHON=ON -DBUILD_TESTS=OFF -A x64 2>&1 | tee cmake_output.log; then
    echo "[OK] CMake configuration successful"

    # Check if Python bindings target was created
    if grep -q "pyopengv" cmake_output.log || [ -f CMakeCache.txt ]; then
        echo "[OK] Python bindings target configured"
    fi

    cd ..
    rm -rf build_test
else
    echo "ERROR: CMake configuration failed"
    cd ..
    rm -rf build_test
    exit 1
fi

echo ""
echo "==========================================="
echo "All tests passed! [OK]"
echo "The conda environment is ready for building OpenGV."
echo ""
echo "To build OpenGV:"
echo "  bash scripts/build_windows_conda.sh"
echo ""
echo "Or manually:"
echo "  1. mkdir build && cd build"
echo "  2. cmake .. -DBUILD_PYTHON=ON -A x64"
echo "  3. cmake --build . --config Release"
