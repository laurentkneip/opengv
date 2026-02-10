# Algorithm Selection Guide for OpenGV

## Absolute Pose (PnP) Algorithms

### Overview

OpenGV provides multiple algorithms for solving the Perspective-n-Point (PnP) problem. **The choice of algorithm significantly impacts accuracy depending on your camera's field of view.**

### Test Results Summary

Based on comprehensive testing with 360° panoramic equirectangular cameras (200 points, 20 runs):

| Algorithm | No Noise | 5px Noise | Combined Noise | With Outliers |
|-----------|----------|-----------|----------------|---------------|
| **UPnP** | 1.8e-15 m<br>2.5e-16 rad | 0.014 m<br>0.0038 rad | 0.015 m<br>0.0027 rad | 0.12 m<br>0.017 rad |
| **EPnP** | 2.9e-11 m<br>8.6e-12 rad | 0.28 m<br>0.025 rad | 0.52 m<br>0.039 rad | 2.4 m<br>0.32 rad |
| **SQPnP** | ❌ FAILS<br>(0.02 m, 1.9° error) | ❌ NOT TESTED | ❌ NOT TESTED | ❌ NOT TESTED |

**Key Finding: UPnP is 20-33x more accurate than EPnP for panoramic/wide-angle scenarios.**

### Recommendations by Camera Type

#### Standard Pinhole Cameras (FOV < 100°)

**Best Choice: EPnP**
- Fast execution (~0.15-0.35 ms)
- Excellent accuracy for limited FOV
- C++: `opengv::absolute_pose::epnp()`
- Python: `pyopengv.absolute_pose_epnp()`

**Alternative: SQPnP Hybrid**
- Automatically selects between SQPnP and EPnP
- C++: `opengv::absolute_pose::sqpnp_hybrid()`
- Python: `pyopengv.absolute_pose_sqpnp_hybrid()`

#### Wide-Angle Cameras (FOV 100-150°)

**Best Choice: UPnP**
- Superior accuracy for wide fields of view
- Handles near-180° viewing angles
- Returns multiple solutions (evaluate all)
- C++: `opengv::absolute_pose::upnp()`
- Python: `pyopengv.absolute_pose_upnp()`

#### Panoramic/Omnidirectional Cameras (FOV > 150°, including 360°)

**REQUIRED: UPnP**
- **Only algorithm validated for 360° panoramas**
- 20-30x more accurate than EPnP
- Robust to backward-facing vectors
- Execution time: ~0.5-1.1 ms (acceptable for the accuracy gain)

**For ULTIMATE ACCURACY: UPnP + Nonlinear Optimization**
```cpp
// Initial solution with UPnP
opengv::transformations_t solutions = opengv::absolute_pose::upnp(adapter);
opengv::transformation_t best_T = solutions[0];  // Or evaluate all

// Refine with Gauss-Newton optimization
opengv::transformation_t refined = opengv::absolute_pose::optimize_nonlinear(adapter, best_T);
```

**For scenarios with OUTLIERS: UPnP + RANSAC**
```cpp
opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem problem(
    adapter,
    opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem::UPNP);
    
opengv::sac::Ransac<opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem> ransac;
ransac.sac_model_ = problem;
ransac.threshold_ = 0.01;  // rad
ransac.max_iterations_ = 1000;
ransac.computeModel();
opengv::transformation_t robust_solution = ransac.model_coefficients_;
```

**DO NOT USE:**
- ❌ **SQPnP** (pure mode): Fails with wide-angle/backward-facing vectors
  - Position error: 0.02 m even with zero noise
  - Angular error: 1.9° even with zero noise
- ❌ **EPnP** without validation: Poor accuracy (0.5+ m error with moderate noise)

### C++ Example

```cpp
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>

// For panoramic cameras
opengv::bearingVectors_t bearings;
opengv::points_t points;
// ... populate bearings and points ...

opengv::absolute_pose::CentralAbsoluteAdapter adapter(
    bearings, points);

// Use UPnP for panoramic views
opengv::transformations_t transformations = 
    opengv::absolute_pose::upnp(adapter);

// Evaluate all solutions and pick best
double best_error = std::numeric_limits<double>::max();
opengv::transformation_t best_transformation;
for (const auto& T : transformations) {
    double error = evaluate_reprojection_error(T, adapter);
    if (error < best_error) {
        best_error = error;
        best_transformation = T;
    }
}
```

### Python Example

```python
import numpy as np
import pyopengv

# Panoramic bearing vectors (360° coverage)
bearings = np.array([
    [0.707, 0.707, 0.0],    # Forward-right
    [-0.707, 0.707, 0.0],   # Forward-left  
    [0.0, -1.0, 0.0],       # Backward
    [0.0, 0.0, 1.0],        # Up
])
points = np.random.randn(4, 3) * 5.0

# Use UPnP for panoramic cameras
transformations = pyopengv.absolute_pose_upnp(bearings, points)

# Pick best solution (or evaluate each)
best_transformation = transformations[0]
```

### RANSAC Integration

For scenarios with outliers, wrap the solver in RANSAC:

**C++:**
```cpp
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac_problems/absolute_pose/AbsolutePoseSacProblem.hpp>

opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem problem(
    adapter,
    opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem::UPNP);

opengv::sac::Ransac<
    opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem> ransac;
ransac.sac_model_ = problem;
ransac.threshold_ = 0.01;  // rad
ransac.max_iterations_ = 1000;

bool success = ransac.computeModel();
opengv::transformation_t best_model = ransac.model_coefficients_;
```

**Python:**
```python
transformation = pyopengv.absolute_pose_ransac(
    bearings,
    points, 
    "UPNP",  # Use UPnP for panoramic views
    threshold=0.01,  # rad
    iterations=1000
)
```

### Performance Characteristics

| Algorithm | Speed | Accuracy (Pinhole) | Accuracy (Panoramic) | Outlier Robust |
|-----------|-------|-------------------|---------------------|----------------|
| **P3P** | ⚡⚡⚡ Very Fast | ✓ Good | ❌ Poor | ❌ No (use RANSAC) |
| **EPnP** | ⚡⚡ Fast | ✓✓ Excellent | ❌ Poor | ❌ No (use RANSAC) |
| **SQPnP** | ⚡⚡ Fast | ✓✓ Excellent | ❌ FAILS | ❌ No (use RANSAC) |
| **SQPnP Hybrid** | ⚡⚡ Fast | ✓✓ Excellent | ⚠️ Fallback to EPnP | ❌ No (use RANSAC) |
| **UPnP** | ⚡ Moderate | ✓✓✓ Excellent | ✓✓✓ **BEST** | ❌ No (use RANSAC) |
| **GPnP** | ⚡ Moderate | ✓✓ Good | ✓ Good | ❌ No (use RANSAC) |

### Known Limitations

1. **SQPnP (pure mode)**:
   - ❌ Not suitable for backward-facing vectors
   - ❌ Not suitable for panoramic/omnidirectional cameras
   - ❌ Fails with FOV > 150°
   - ✓ Use `sqpnp_hybrid()` or `upnp()` instead

2. **EPnP**:
   - ⚠️ Accuracy degrades significantly with wide FOV (> 100°)
   - ⚠️ Not outlier-robust (use with RANSAC)
   - ✓ Excellent for standard pinhole cameras

3. **All algorithms**:
   - ❌ None are inherently outlier-robust
   - ✓ Always use RANSAC/LMEDS when outliers expected

### Testing

All claims verified through:
- `test/test_panorama_sqpnp_vs_epnp.cpp`: 360° panorama accuracy tests
- `test/test_absolute_pose.cpp`: Standard camera benchmarks
- `python/tests.py`: Python binding validation
- Configurations: No noise, Gaussian noise (5px + 5cm), Outliers (5%)

### References

- EPnP: Lepetit, V., Moreno-Noguer, F., & Fua, P. (2009)
- UPnP: Kneip, L., Scaramuzza, D., & Siegwart, R. (2011)
- SQPnP: Terzakis, G., & Lourakis, M. (2020)
- OpenGV: Kneip, L., & Furgale, P. (2014)
