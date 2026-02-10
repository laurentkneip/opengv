# OpenGV Python Tests

## Test Structure

The Python tests have been refactored into modular files for better organization and maintainability.

### Test Files

1. **`test_utils.py`** - Shared utility functions
   - `normalized()` - Vector normalization
   - `generateRandomPoint()` - Random 3D point generation
   - `generateRandomTranslation()` - Random translation generation
   - `generateRandomRotation()` - Random rotation matrix generation
   - `calculate_ransac_threshold()` - RANSAC threshold computation

2. **`test_relative_pose.py`** - Relative pose estimation (camera-to-camera)
   - `test_relative_pose()` - Basic 8-point algorithm
   - `test_relative_pose_ransac()` - RANSAC-based estimation
   - `test_relative_pose_ransac_rotation_only()` - Rotation-only RANSAC

3. **`test_triangulation.py`** - 3D point reconstruction
   - `test_triangulation()` - Triangulation from two views

4. **`test_absolute_pose_basic.py`** - Absolute pose estimation (camera pose from 3D-2D)
   - `test_absolute_pose()` - Tests P3P (Kneip, Gao), EPnP, UPnP, SQPnP
   - `test_absolute_pose_ransac()` - RANSAC with KNEIP and EPNP
   - `test_absolute_pose_lmeds()` - LMedS robust estimation

5. **`test_panorama.py`** - Panoramic 360° algorithm comparison
   - `test_panorama_360()` - Comprehensive test with varied depth ranges (0.5-50m)
   - Tests UPnP, EPnP, SQPnP on full sphere with backward-facing vectors
   - Validates that UPnP is superior for wide-angle/panoramic cameras

6. **`run_all_tests.py`** - Master test runner
   - Runs all test modules in sequence
   - Provides summary of pass/fail status

7. **`tests.py`** - Original monolithic test file (**DEPRECATED**)
   - Kept for backward compatibility until all tests migrated
   - **Use the modular test files instead**

## Running Tests

### Run All Tests
```bash
cd python
python run_all_tests.py
```

### Run Individual Test Modules
```bash
python test_relative_pose.py
python test_triangulation.py
python test_absolute_pose_basic.py
python test_panorama.py
```

### Run Original Monolithic Tests (Deprecated)
```bash
python tests.py
```

## Key Test Insights

### Panorama Test Results
The panorama test validates algorithm performance on 360° panoramic cameras:

- **UPnP**: ~1e-15 m error (perfect) - only truly universal global solver
- **EPnP**: ~1e-11 m error (good) - degrades with extreme depth variation
- **SQPnP**: ~0.02 m error (FAILS) - not designed for backward-facing vectors

**Conclusion**: Use UPnP for panoramic/omnidirectional cameras (FOV > 150°)


