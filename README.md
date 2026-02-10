# OpenGV

OpenGV is a collection of computer vision methods for solving geometric vision problems. It contains absolute-pose, relative-pose, triangulation, and point-cloud alignment methods for the calibrated case.

**Author:** Laurent Kneip, ShanghaiTech, Mobile Perception Lab  
**Homepage:** http://laurentkneip.github.io/opengv

## Features

- **Absolute Pose Solvers**: P3P (Kneip, Gao), EPnP, UPnP, GPnP, **SQPnP**
- **Relative Pose Solvers**: 2-point, 5-point (Nister, Kneip), 7-point, 8-point, Eigensolver
- **RANSAC/LMedS Integration**: Robust estimation with outlier rejection
- **Omnidirectional Camera Support**: Full 360° bearing vector handling
- **Triangulation**: Multiple methods for 3D point reconstruction
- **Python Bindings**: Full API access via `pyopengv`

## Installation

### Installation on macOS

#### Prerequisites

Install the required dependencies using Homebrew:

```bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install dependencies
brew install cmake eigen pybind11

# For Python bindings, install Python packages
pip3 install numpy scipy
```

**Alternative: Using Conda** (optional, not required):

If you prefer using conda instead of Homebrew:

```bash
conda env create -f environment.yml
conda activate opengv
bash scripts/test_conda_env.sh
```

### Installation on Linux

#### Prerequisites

Install the required dependencies using your package manager:

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y cmake build-essential libeigen3-dev python3-dev python3-pip python3-venv pybind11-dev

# Fedora/RHEL
sudo dnf install cmake eigen3-devel python3-devel python3-pip pybind11-devel

# For Python bindings
pip3 install numpy scipy
```

**Alternative: Using Conda** (optional, not required):

If you prefer using conda:

```bash
conda env create -f environment.yml
conda activate opengv
bash scripts/test_conda_env.sh
```

### Installation on Windows

#### Recommended: Using Conda

Windows is easiest with Conda, which handles all dependencies:

```powershell
# Install Miniconda or Anaconda if you haven't already
# Download from: https://docs.conda.io/en/latest/miniconda.html

# Create a conda environment from the provided environment.yml file
conda env create -f environment.yml

# Activate the environment
conda activate opengv
# If needed: conda init powershell (then restart PowerShell)

# Verify the environment is set up correctly
python --version
cmake --version

# Test the conda environment (optional but recommended)
scripts\test_conda_env.bat
```

### Building the C++ Library

First, clone the repository and initialize submodules (all platforms):

```bash
git clone https://github.com/laurentkneip/opengv.git
cd opengv
git submodule update --init --recursive
```

#### macOS

```bash
mkdir build && cd build
cmake ..
make -j$(sysctl -n hw.ncpu)
make test
```

#### Linux

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
make test
```

#### Windows (PowerShell)

```powershell
mkdir build
cd build
cmake .. -A x64
cmake --build . --config Release --parallel
ctest -C Release --output-on-failure --parallel
```

#### Windows (Git Bash with Clang + Ninja) - Faster Builds

For significantly faster build times on Windows, use Clang and Ninja in Git Bash. This method is much faster than MSBuild and provides better compile-time diagnostics.

**Prerequisites:**
```bash
# Install LLVM (includes Clang) and Ninja via Chocolatey
choco install llvm ninja

# Install vcpkg for Eigen and pybind11
git clone https://github.com/Microsoft/vcpkg.git C:\vcpkg
C:\vcpkg\bootstrap-vcpkg.bat
C:\vcpkg\vcpkg install eigen3:x64-windows pybind11:x64-windows
```

**Build:**
```bash
mkdir build_clang && cd build_clang
cmake .. -G Ninja \
  -DCMAKE_CXX_COMPILER=clang++ \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTS=ON \
  -DBUILD_PYTHON=ON \
  -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake

# Build with 8 parallel jobs
ninja -j 8

# Run tests
ctest --output-on-failure
```

### Building with Python Bindings

If using conda, activate the environment first: `conda activate opengv`

#### macOS

```bash
mkdir build && cd build
cmake .. -DBUILD_PYTHON=ON -DPYTHON_EXECUTABLE=$(which python3)
make -j$(sysctl -n hw.ncpu)
# Output: build/lib/pyopengv.cpython-*.so
```

#### Linux

```bash
mkdir build && cd build
cmake .. -DBUILD_PYTHON=ON -DPYTHON_EXECUTABLE=$(which python3)
make -j$(nproc)
# Output: build/lib/pyopengv.cpython-*.so
```

#### Windows (PowerShell)

```powershell
conda activate opengv
mkdir build
cd build
cmake .. -DBUILD_PYTHON=ON -A x64 -DPYTHON_EXECUTABLE=$env:CONDA_PREFIX\python.exe
cmake --build . --config Release --parallel
# Output: build/lib/Release/pyopengv.pyd
```

**Note**: If using conda, CMake will auto-detect conda Python. You can also explicitly specify it with:
```powershell
cmake .. -DBUILD_PYTHON=ON -A x64 -DPYTHON_EXECUTABLE=$env:CONDA_PREFIX\python.exe
```

### Installing the Python Wheel

The easiest way to install the Python bindings:

#### macOS / Linux

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install numpy scipy
pip install .
```

#### Windows (PowerShell)

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install numpy scipy
pip install .
```

Or build a distributable wheel (all platforms):

```bash
pip install wheel build
pip wheel . --no-deps -w dist/
pip install dist/pyopengv-*.whl
```

### Verify Installation

#### macOS / Linux

```bash
cd build
./bin/test_absolute_pose
python3 -c "import pyopengv; print('pyopengv installed successfully!')"
pip install pytest
PYTHONPATH=build/lib pytest python/ -v
```

#### Windows (PowerShell)

```powershell
cd build
.\bin\Release\test_absolute_pose.exe
python -c "import pyopengv; print('pyopengv installed successfully!')"
pip install pytest
$env:PYTHONPATH="$PWD\lib\Release"; pytest ..\python\ -v
```

## Quick Start (Python)

```python
import numpy as np
import pyopengv

# Generate sample data
def normalized(x):
    return x / np.linalg.norm(x)

# Random camera pose
position = np.random.uniform(-2, 2, 3)
rotation = np.eye(3)  # Identity rotation

# Generate 3D points and their bearing vectors
num_points = 20
points = np.random.uniform(-5, 5, (num_points, 3))
bearings = np.array([normalized(rotation.T @ (p - position)) for p in points])

# Solve absolute pose with different methods
result_epnp = pyopengv.absolute_pose_epnp(bearings, points)
result_sqpnp = pyopengv.absolute_pose_sqpnp(bearings, points)

print("EPnP result:\n", result_epnp)
print("SQPnP result:\n", result_sqpnp)

# RANSAC with outliers
result_ransac = pyopengv.absolute_pose_ransac(
    bearings, points, "KNEIP", threshold=0.01, iterations=1000)
print("RANSAC result:\n", result_ransac)
```

## Available Python Functions

### Absolute Pose
- `absolute_pose_p3p_kneip(bearings, points)` - P3P Kneip solver
- `absolute_pose_p3p_gao(bearings, points)` - P3P Gao solver
- `absolute_pose_epnp(bearings, points)` - EPnP solver
- `absolute_pose_sqpnp(bearings, points)` - SQPnP solver (omnidirectional support)
- `absolute_pose_ransac(bearings, points, algo, threshold, iterations)` - RANSAC
- `absolute_pose_lmeds(bearings, points, algo, threshold, iterations)` - LMedS

**RANSAC/LMedS algorithms**: `"KNEIP"`, `"GAO"`, `"EPNP"`, `"SQPNP"`

### Relative Pose
- `relative_pose_twopt(b1, b2, rotation)` - 2-point with known rotation
- `relative_pose_fivept_nister(b1, b2)` - 5-point Nister
- `relative_pose_fivept_kneip(b1, b2)` - 5-point Kneip
- `relative_pose_sevenpt(b1, b2)` - 7-point algorithm
- `relative_pose_eightpt(b1, b2)` - 8-point algorithm
- `relative_pose_ransac(b1, b2, algo, threshold, iterations)` - RANSAC

### Triangulation
- `triangulation_triangulate(b1, b2, position, rotation)` - Triangulate points
- `triangulation_triangulate2(b1, b2, position, rotation)` - Alternative method

## Troubleshooting

### Eigen not found

**On macOS (Homebrew):**
```bash
# Install Eigen
brew install eigen

# If CMake still can't find it, specify the path:
cmake .. -DEIGEN_INCLUDE_DIR=/opt/homebrew/include/eigen3
```

**On Windows (Conda):**
```bash
# Eigen should be available through conda, but if not found:
conda install eigen
# Or specify the path manually:
cmake .. -DEIGEN_INCLUDE_DIR=$CONDA_PREFIX/include/eigen3
```

**On Linux:**
```bash
sudo apt-get install libeigen3-dev
# Or on Fedora:
sudo dnf install eigen3-devel
```

### C++14 required error
The library requires C++14. This should be handled automatically, but if you see errors:
```bash
cmake .. -DCMAKE_CXX_STANDARD=14
```

### Python module not found
Ensure you're using the correct Python:
```bash
# On macOS/Linux:
cmake .. -DBUILD_PYTHON=ON -DPYTHON_EXECUTABLE=$(which python3)

# On Windows (with conda):
cmake .. -DBUILD_PYTHON=ON -DPYTHON_EXECUTABLE=$CONDA_PREFIX\python.exe
```

### PyBind11 not found (Windows)
If using conda, PyBind11 should be available automatically. If not:
```bash
conda install pybind11
# Or install via pip (less recommended):
pip install pybind11
```

### Windows-specific issues

1. **MSVC compiler required**: On Windows, you need Visual Studio with C++ build tools.

   **Option A: Install via Visual Studio (Recommended)**
   1. Download [Visual Studio 2022 Community](https://visualstudio.microsoft.com/downloads/) (free)
   2. Run the installer
   3. Select **"Desktop development with C++"** workload
   4. Click Install (requires ~6-8 GB)
   5. Restart your terminal after installation

   **Option B: Install Build Tools only (smaller download)**
   1. Download [Build Tools for Visual Studio 2022](https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2022)
   2. Run the installer
   3. Select **"C++ build tools"** workload
   4. Make sure these are checked:
      - MSVC v143 - VS 2022 C++ x64/x86 build tools
      - Windows 10/11 SDK
   5. Click Install

   **Option C: Via Conda (if using conda environment)**
   ```powershell
   conda activate opengv
   conda install -c conda-forge compilers
   ```

2. **Path issues**: If CMake can't find Python on Windows:
   ```bash
   # Make sure Python is in PATH or use full path:
   cmake .. -DPYTHON_EXECUTABLE=C:\Users\YourName\miniconda3\envs\opengv\python.exe
   ```

3. **Library naming**: On Windows, the Python extension is `.pyd` instead of `.so`. The output location may also differ:
   - Debug: `build/lib/Debug/pyopengv.pyd`
   - Release: `build/lib/Release/pyopengv.pyd`

## Citation

If you use OpenGV in your research, please cite:

```bibtex
@inproceedings{kneip2014opengv,
  title={OpenGV: A unified and generalized approach to real-time calibrated geometric vision},
  author={Kneip, Laurent and Furgale, Paul},
  booktitle={IEEE International Conference on Robotics and Automation (ICRA)},
  year={2014}
}
```

For the SQPnP solver:

```bibtex
@inproceedings{terzakis2020sqpnp,
  title={A Consistently Fast and Globally Optimal Solution to the Perspective-n-Point Problem},
  author={Terzakis, George and Lourakis, Manolis},
  booktitle={European Conference on Computer Vision (ECCV)},
  year={2020}
}
```

## License

See [License.txt](License.txt) for details.

