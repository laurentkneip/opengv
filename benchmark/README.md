# OpenGV Benchmarks

This folder contains benchmark programs for evaluating the performance of various pose estimation algorithms in OpenGV.

## Benchmark Programs

### benchmark_panorama_algorithms

Benchmarks the performance of different absolute pose algorithms (SQPnP, EPnP, UPnP) on 360-degree equirectangular panorama data.

**Features:**
- Tests algorithms under various noise conditions
- Measures both accuracy (position and rotation errors) and timing performance
- Evaluates robustness to outliers
- Uses realistic panorama projection models

**Benchmark Scenarios:**
1. No Noise - Pure algorithmic performance baseline
2. 5 pixel noise - Bearing vector noise from pixel measurements
3. 5cm 3D point noise - Point localization uncertainty
4. Combined noise - Realistic scenario with both sources
5. 5% outliers - Moderate outlier contamination
6. 10% outliers - High outlier stress test
7. High noise - 10 pixel + 10cm noise stress test

**Usage:**
```bash
./benchmark_panorama_algorithms
```

**Output:**
- Averaged position error (meters)
- Averaged angular error (degrees)
- Averaged execution time (microseconds and milliseconds)
- Comparison table across all algorithms

**Building:**
The benchmark is built automatically with the project when `BUILD_TESTS` is enabled:
```bash
cmake -DBUILD_TESTS=ON ..
make benchmark_panorama_algorithms
```

## Notes

- Benchmarks are **not** run as part of `make test` / `ctest`
- They are intended for performance analysis, not pass/fail validation
- Each benchmark runs 100 iterations by default for statistical robustness
- Results may vary based on hardware and random seed
- **Important**: SQPnP is not suitable for mixed hemisphere bearing vectors (e.g., forward + backward facing views). UPnP is recommended for full 360° panoramic scenarios.
