#!/bin/bash
# Build script for macOS - mirrors the GitHub Actions CI workflow
# Usage: ./scripts/build_macos.sh
#
# This script performs the exact same steps as the build-and-test-macos workflow
# Compatible with both Intel (x86_64) and Apple Silicon (ARM64) Macs

set -e  # Exit on any error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

echo "========================================"
echo "OpenGV macOS Build"
echo "========================================"
echo ""
echo "Architecture: $(uname -m)"
echo "macOS version: $(sw_vers -productVersion)"
echo ""

# Step 1: Install dependencies
echo "Step 1: Installing dependencies via Homebrew..."
echo "-----------------------------------------------"
if ! command -v brew &> /dev/null; then
    echo "ERROR: Homebrew not found. Please install from https://brew.sh"
    exit 1
fi

brew install cmake eigen pybind11
echo "[OK] Dependencies installed"
echo ""

# Step 2: Create virtual environment and install Python dependencies
echo "Step 2: Setting up Python environment..."
echo "-----------------------------------------"
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install numpy scipy pytest wheel build
echo "[OK] Python dependencies installed"
echo ""

# Step 3: Configure CMake (with Python and tests)
echo "Step 3: Configuring CMake..."
echo "----------------------------"
rm -rf build
mkdir -p build
cd build
cmake .. \
    -DCMAKE_CXX_STANDARD=14 \
    -DBUILD_TESTS=ON \
    -DBUILD_PYTHON=ON \
    -DPYTHON_EXECUTABLE=$(which python3)
echo "[OK] CMake configured"
echo ""

# Step 4: Build everything (C++ and Python)
echo "Step 4: Building (C++ and Python)..."
echo "-------------------------------------"
make -j$(sysctl -n hw.ncpu)
echo "[OK] Build completed"
echo ""

# Step 5: Run C++ tests
echo "Step 5: Running C++ tests..."
echo "----------------------------"
ctest --output-on-failure
echo "[OK] C++ tests completed"
echo ""

# Step 6: Test Python module (direct build)
echo "Step 6: Testing Python module (direct build)..."
echo "------------------------------------------------"
cd ..
PYTHONPATH=build/lib pytest python/ -v
echo "[OK] Python tests completed"
echo ""

# Step 7: Build wheel
echo "Step 7: Building wheel..."
echo "-------------------------"
rm -rf dist
pip wheel . --no-deps -w dist/

if ls dist/pyopengv-*.whl 1> /dev/null 2>&1; then
    echo "[OK] Wheel built successfully"
    ls -lh dist/
else
    echo "[ERROR] Wheel build failed"
    exit 1
fi
echo ""

# Step 8: Install and test wheel
echo "Step 8: Testing wheel installation..."
echo "--------------------------------------"
# Unset PYTHONPATH to test clean wheel install
unset PYTHONPATH
pip install dist/pyopengv-*.whl --force-reinstall
pytest python/ -v
echo "[OK] Wheel installation and tests passed"
echo ""

echo "========================================"
echo "macOS build completed successfully!"
echo "========================================"
echo ""
echo "Outputs:"
echo "  C++ library: build/lib/libopengv.a"
echo "  Python module: build/lib/pyopengv.*.so"
echo "  Wheel: dist/pyopengv-*.whl"
echo ""
echo "To test Python module:"
echo "  python3 -c 'import pyopengv; print(\"✓ Success\")'"
