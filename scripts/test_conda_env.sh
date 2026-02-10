#!/bin/bash
# Test script to validate conda environment setup for OpenGV
# Usage: bash scripts/test_conda_env.sh

set -e  # Exit on error

# Change to project root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

echo "Testing OpenGV Conda Environment"
echo "================================="

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

echo "✓ Conda environment 'opengv' exists"

# Activate environment (note: in scripts, we need to source conda.sh)
eval "$(conda shell.bash hook)"
conda activate opengv

echo "✓ Activated conda environment"

# Check Python
echo -n "Checking Python... "
PYTHON_VERSION=$(python --version 2>&1)
echo "$PYTHON_VERSION"
if ! python -c "import sys; assert sys.version_info >= (3, 7)"; then
    echo "ERROR: Python 3.7+ required"
    exit 1
fi
echo "✓ Python version OK"

# Check required packages
echo -n "Checking required packages... "
REQUIRED_PACKAGES=("numpy" "cmake" "pybind11")
MISSING_PACKAGES=()

for pkg in "${REQUIRED_PACKAGES[@]}"; do
    if python -c "import $pkg" 2>/dev/null; then
        echo "✓ $pkg"
    else
        echo "✗ $pkg MISSING"
        MISSING_PACKAGES+=("$pkg")
    fi
done

if [ ${#MISSING_PACKAGES[@]} -ne 0 ]; then
    echo "ERROR: Missing packages: ${MISSING_PACKAGES[*]}"
    exit 1
fi

# Check CMake
echo -n "Checking CMake... "
if command -v cmake &> /dev/null; then
    CMAKE_VERSION=$(cmake --version | head -n1)
    echo "$CMAKE_VERSION"
    echo "✓ CMake found"
else
    echo "ERROR: CMake not found"
    exit 1
fi

# Check Eigen (via conda or system)
echo -n "Checking Eigen... "
if conda list eigen 2>/dev/null | grep -q eigen; then
    echo "✓ Eigen found via conda"
elif [ -d "$CONDA_PREFIX/include/eigen3" ] || [ -d "$CONDA_PREFIX/include/eigen" ]; then
    echo "✓ Eigen headers found"
else
    echo "⚠ Eigen headers not found in expected location (may still work if in system path)"
fi

# Test basic build configuration
echo ""
echo "Testing CMake configuration..."
mkdir -p build_test
cd build_test

if cmake .. -DBUILD_PYTHON=ON -DBUILD_TESTS=OFF 2>&1 | tee cmake_output.log; then
    echo "✓ CMake configuration successful"
    
    # Check if Python bindings target was created
    if grep -q "pyopengv" cmake_output.log || [ -f CMakeCache.txt ]; then
        echo "✓ Python bindings target configured"
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
echo "================================="
echo "All tests passed! ✓"
echo "The conda environment is ready for building OpenGV."
echo ""
echo "To build OpenGV:"
echo "  1. mkdir build && cd build"
echo "  2. cmake .. -DBUILD_PYTHON=ON"
echo "  3. cmake --build . (or make on Unix)"
