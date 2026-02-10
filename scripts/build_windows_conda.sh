#!/bin/bash
# Build and test script for OpenGV on Windows using Conda (for Git Bash / CI)
# This script mirrors the GitHub Actions workflow for local testing
# Usage (in Git Bash with conda): bash scripts/build_windows_conda.sh

set -e  # Exit on error

echo "========================================"
echo "OpenGV Windows Build (Conda + MSVC)"
echo "========================================"
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda is not installed or not in PATH"
    exit 1
fi

# Check if environment exists and activate
if ! conda env list | grep -q "^opengv "; then
    echo "ERROR: 'opengv' conda environment does not exist"
    echo "Please create it first: conda env create -f environment.yml"
    exit 1
fi

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate opengv
echo "[OK] Conda environment 'opengv' activated"

# Step 1: Verify environment
echo ""
echo "Step 1: Verifying environment..."
echo "---------------------------------"
python --version
cmake --version
pythonDir=$(python -c "import sys; print(sys.prefix)")
pythonLibDir="$pythonDir/libs"
echo "Python directory: $pythonDir"
echo "Python lib directory: $pythonLibDir"
echo "Python lib files:"
ls -la "$pythonLibDir" | grep -E "python[0-9]+\.lib" || echo "(no python*.lib found)"
echo "[OK] Environment verified"
echo ""

# Step 2: Install Python dependencies
echo "Step 2: Installing Python dependencies..."
echo "-----------------------------------------"
pip install numpy scipy pytest wheel build --quiet
echo "[OK] Python dependencies installed"
echo ""

# Step 3: Configure CMake
echo "Step 3: Configuring CMake..."
echo "----------------------------"
rm -rf build
mkdir -p build
cd build

cmake .. \
  -DBUILD_TESTS=ON \
  -DBUILD_PYTHON=ON \
  -DPython3_ROOT_DIR="$pythonDir" \
  -A x64

echo "[OK] CMake configured"
echo ""

# Step 4: Build
echo "Step 4: Building (Release, parallel)..."
echo "---------------------------------------"
cmake --build . --config Release --parallel
echo "[OK] Build completed"
echo ""

# Step 5: Run C++ tests
echo "Step 5: Running C++ tests..."
echo "----------------------------"
ctest -C Release --output-on-failure
if [ $? -eq 0 ]; then
    echo "[OK] All C++ tests passed"
else
    echo "[WARN] Some C++ tests failed"
fi
echo ""

# Step 6: Check outputs
echo "Step 6: Checking outputs..."
echo "---------------------------"
if [ -f "lib/Release/opengv.lib" ]; then
    echo "[OK] Found: lib/Release/opengv.lib"
else
    echo "[WARN] opengv.lib not found in lib/Release"
fi

if [ -f "lib/Release/pyopengv.pyd" ]; then
    echo "[OK] Found: lib/Release/pyopengv.pyd"
else
    echo "[WARN] pyopengv.pyd not found in lib/Release"
fi
echo ""

# Step 7: Test Python module (direct build)
echo "Step 7: Testing Python module (direct build)..."
echo "------------------------------------------------"
cd ..
export PYTHONPATH="$PWD/build/lib/Release"
pytest python/ -v
if [ $? -eq 0 ]; then
    echo "[OK] All Python tests passed"
else
    echo "[WARN] Some Python tests failed"
fi
echo ""

# Step 8: Build and test wheel
echo "Step 8: Building and testing wheel..."
echo "--------------------------------------"
pip wheel . --no-deps -w dist/
pip install dist/pyopengv-*.whl --force-reinstall
unset PYTHONPATH
pytest python/ -v
if [ $? -eq 0 ]; then
    echo "[OK] Wheel tests passed"
else
    echo "[WARN] Some wheel tests failed"
fi
echo ""

echo "========================================"
echo "Windows Conda build completed!"
echo "========================================"
echo ""
echo "Outputs:"
echo "  C++ library: build/lib/Release/opengv.lib"
echo "  Python module: build/lib/Release/pyopengv.pyd"
echo "  Wheel: dist/pyopengv-*.whl"
