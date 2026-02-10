#!/bin/bash
# Build and test script for OpenGV on Windows with Clang and Ninja
# Usage (in Git Bash): bash scripts/build_windows_clang_ninja.sh

set -e  # Exit on error

echo "========================================"
echo "OpenGV Windows Build (Clang + Ninja)"
echo "========================================"
echo ""

# Check prerequisites
echo "Checking prerequisites..."
echo "-------------------------"

if ! command -v clang &> /dev/null; then
    echo "ERROR: clang not found. Install with: choco install llvm -y"
    exit 1
fi

if ! command -v ninja &> /dev/null; then
    echo "ERROR: ninja not found. Install with: choco install ninja -y"
    exit 1
fi

if ! command -v python &> /dev/null; then
    echo "ERROR: python not found"
    exit 1
fi

if [ ! -d "/c/vcpkg" ]; then
    echo "ERROR: vcpkg not found at C:\\vcpkg"
    echo "Please install: git clone https://github.com/Microsoft/vcpkg.git C:\\vcpkg"
    exit 1
fi

echo "[OK] clang: $(clang --version | head -n1)"
echo "[OK] ninja: $(ninja --version)"
echo "[OK] python: $(python --version)"
echo "[OK] vcpkg: found at C:\\vcpkg"
echo ""

# Step 1: Install vcpkg dependencies
echo "Step 1: Installing dependencies (vcpkg)..."
echo "------------------------------------------"
if [ ! -f "/c/vcpkg/installed/x64-windows/include/Eigen/Core" ]; then
    echo "Installing eigen3..."
    /c/vcpkg/vcpkg install eigen3:x64-windows
fi

if [ ! -f "/c/vcpkg/installed/x64-windows/include/pybind11/pybind11.h" ]; then
    echo "Installing pybind11..."
    /c/vcpkg/vcpkg install pybind11:x64-windows
fi
echo "[OK] vcpkg dependencies installed"
echo ""

# Step 2: Install Python dependencies (from pyproject.toml deps)
echo "Step 2: Installing Python dependencies..."
echo "-----------------------------------------"
python -m pip install --upgrade pip --quiet
pip install pybind11 numpy scipy pytest wheel build --quiet
echo "[OK] Python dependencies installed"
echo ""

# Step 3: Configure CMake with Clang and Ninja
echo "Step 3: Configuring CMake..."
echo "----------------------------"
export CC=clang
export CXX=clang++

rm -rf build_clang
mkdir -p build_clang
cd build_clang

cmake .. \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTS=ON \
  -DBUILD_PYTHON=ON \
  -DPYTHON_EXECUTABLE=$(which python) \
  -DCMAKE_TOOLCHAIN_FILE=/c/vcpkg/scripts/buildsystems/vcpkg.cmake

echo "[OK] CMake configured"
echo ""

# Step 4: Build with Ninja
echo "Step 4: Building with Ninja..."
echo "-------------------------------"
ninja -j$(nproc)
echo "[OK] Build completed"
echo ""

# Step 5: Run C++ tests
echo "Step 5: Running C++ tests..."
echo "----------------------------"
ctest --output-on-failure
if [ $? -eq 0 ]; then
    echo "[OK] All C++ tests passed"
else
    echo "[WARN] Some C++ tests failed"
fi
echo ""

# Step 6: Test Python module (direct build)
echo "Step 6: Testing Python module (direct build)..."
echo "------------------------------------------------"
cd ..
export PYTHONPATH="$PWD/build_clang/lib"
pytest python/ -v

if [ $? -eq 0 ]; then
    echo "[OK] All Python tests passed"
else
    echo "[WARN] Some Python tests failed"
fi
echo ""

echo "========================================"
echo "Build and test completed!"
echo "========================================"
echo ""
echo "Outputs:"
echo "  C++ library: build_clang/lib/opengv.lib"
echo "  Python module: build_clang/lib/pyopengv.pyd"
echo ""
echo "To test Python module:"
echo "  python -c 'import pyopengv; print(\"✓ Success\")'"
