#!/bin/bash
# Build script for Ubuntu 22.04 - mirrors the GitHub Actions CI workflow
# Usage: ./scripts/build_ubuntu22.04.sh
#
# This script performs the exact same steps as the build-wheel-ubuntu workflow

set -e  # Exit on any error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

echo "========================================"
echo "OpenGV Ubuntu 22.04 Build"
echo "========================================"
echo ""

# Step 1: Install dependencies
echo "Step 1: Installing system dependencies..."
echo "----------------------------------------"
sudo apt-get update
sudo apt-get install -y \
    cmake \
    build-essential \
    libeigen3-dev \
    python3-dev \
    python3-pip \
    python3-venv \
    pybind11-dev
echo "[OK] System dependencies installed"
echo ""

# Step 2: Create virtual environment and install Python dependencies
echo "Step 2: Setting up Python environment..."
echo "-----------------------------------------"
rm -rf .venv
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install numpy scipy pytest wheel build
echo "[OK] Python dependencies installed"
echo ""

# Step 3: Build wheel
echo "Step 3: Building wheel..."
echo "-------------------------"
rm -rf dist
pip wheel . --no-deps --no-build-isolation -w dist/
echo "[OK] Wheel built successfully"
ls -la dist/
echo ""

# Step 4: Install wheel and test
echo "Step 4: Installing wheel and running tests..."
echo "----------------------------------------------"
pip install dist/pyopengv-*.whl
pytest python/ -v
echo "[OK] Wheel installed and Python tests completed"
echo ""

# Step 5: Build with CMake for C++ tests
echo "Step 5: Building C++ tests (CMake)..."
echo "--------------------------------------"
rm -rf build
mkdir -p build
cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTS=ON \
    -DBUILD_PYTHON=OFF
cmake --build . --parallel
echo "[OK] C++ build completed"
echo ""

# Step 6: Run C++ tests
echo "Step 6: Running C++ tests..."
echo "----------------------------"
ctest --output-on-failure
echo "[OK] C++ tests completed"
echo ""

echo "========================================"
echo "Ubuntu 22.04 build completed successfully!"
echo "========================================"
echo ""
echo "Wheel location: dist/pyopengv-*.whl"
echo ""
