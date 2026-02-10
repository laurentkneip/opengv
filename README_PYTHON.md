# PyOpenGV - Python Bindings for OpenGV

Python bindings for OpenGV, a library for solving calibrated central and non-central geometric vision problems.

## Features

- **Absolute Pose**: P3P, EPnP, UPnP, GPnP, SQPnP solvers
- **Relative Pose**: 5-point, 7-point, 8-point, eigensolver methods
- **Omnidirectional Support**: Angular consistency validation for panoramic cameras
- **RANSAC Integration**: Robust estimation with sample consensus
- **Triangulation**: Multiple methods for 3D point reconstruction

## Installation

### From PyPI (when published)

```bash
pip install pyopengv
```

### From Source

1. Install dependencies:
```bash
./build_wheel.sh install-deps
```

2. Build and install:
```bash
./build_wheel.sh build
./build_wheel.sh install-wheel
```

Or for development:
```bash
pip install -e .
```

## Quick Start

```python
import numpy as np
import pyopengv

# Example: Absolute pose estimation with P3P
bearings = np.random.randn(3, 3)  # 3D bearing vectors
points = np.random.randn(3, 3)     # 3D world points

# Normalize bearings
bearings = bearings / np.linalg.norm(bearings, axis=1, keepdims=True)

# Compute pose
transformation = pyopengv.absolute_pose.p3p_kneip(bearings, points)
```

## Algorithm Selection Guide

### Absolute Pose (PnP) Algorithms

**For Standard Pinhole Cameras (FOV < 100°):**
- **EPnP** (`absolute_pose_epnp`): Fast and accurate for standard cameras
- **P3P** (`absolute_pose_p3p_kneip`, `absolute_pose_p3p_gao`): Minimal solver (3 points)
- **SQPnP Hybrid** (`absolute_pose_sqpnp_hybrid`): Automatically switches between SQPnP and EPnP

**For Panoramic/Wide-Angle/360° Cameras (FOV > 120°):**
- **UPnP** (`absolute_pose_upnp`) - **RECOMMENDED**: 20-30x more accurate than EPnP for panoramic views
  - Best accuracy for wide-angle scenarios
  - Returns multiple solutions (pick best)
  - ~3-4x slower than EPnP but worth it for accuracy
- **UPnP + Nonlinear Refinement** - **ULTIMATE ACCURACY**: Use UPnP result as initial guess for optimization
  ```python
  transformations = pyopengv.absolute_pose_upnp(bearings, points)
  refined = pyopengv.absolute_pose_optimize_nonlinear(bearings, points, transformations[0])
  ```

**Algorithms to AVOID for Panoramic Views:**
- ❌ **SQPnP** (`absolute_pose_sqpnp`): Not designed for wide-angle/backward-facing vectors
- Use `upnp` or `sqpnp_hybrid` instead

### Example: Panoramic Camera Pose Estimation

```python
import numpy as np
import pyopengv

# 360° panoramic bearing vectors (can point in any direction)
bearings = np.array([
    [0.707, 0.707, 0.0],    # Forward-right
    [-0.707, 0.707, 0.0],   # Forward-left
    [0.0, -1.0, 0.0],       # Backward
    [0.0, 0.0, 1.0],        # Up
])
points = np.random.randn(4, 3) * 5.0  # 3D world points

# Use UPnP for panoramic cameras
transformations = pyopengv.absolute_pose.upnp(bearings, points)

# UPnP returns multiple solutions - pick the one with best reprojection error
best_transformation = transformations[0]  # Or evaluate each solution
```

### RANSAC for Outlier Rejection

```python
# For panoramic cameras with outliers
transformation = pyopengv.absolute_pose_ransac(
    bearings, 
    points,
    "UPNP",  # Use UPnP for panoramic views
    threshold=0.01,  # rad
    iterations=1000
)
```

## Building Wheels

The repository includes a helper script for building and distributing wheels:

```bash
# Build wheel locally
./build_wheel.sh build

# Install from built wheel
./build_wheel.sh install-wheel

# Upload to TestPyPI (for testing)
./build_wheel.sh upload-test

# Upload to PyPI (production)
./build_wheel.sh upload
```

## Requirements

- Python >= 3.7
- NumPy >= 1.19.0
- CMake >= 3.15 (for building from source)
- C++11 compatible compiler
- Eigen3 library

## License

BSD 3-Clause License. See LICENSE.txt for details.

## Citation

If you use OpenGV in your research, please cite:

```
@article{kneip2014opengv,
  title={OpenGV: A unified and generalized approach to real-time calibrated geometric vision},
  author={Kneip, Laurent and Furgale, Paul},
  journal={IEEE International Conference on Robotics and Automation (ICRA)},
  year={2014}
}
```

## Links

- Homepage: http://laurentkneip.github.io/opengv
- Documentation: http://laurentkneip.github.io/opengv
- Repository: https://github.com/laurentkneip/opengv
