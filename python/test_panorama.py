"""
Tests for panoramic 360° algorithm comparison.

This test validates UPnP, EPnP, and SQPnP performance on full sphere panoramic views
with varied depth ranges (near, medium, far) including backward-facing vectors.
"""

import time

import numpy as np
from test_utils import generateRandomTranslation, normalized

import pyopengv


def test_panorama_360():
    """Test algorithms on full 360° panoramic views with varied depth ranges.

    Uses the same approach as the C++ test_panorama_algorithms:
    - Generates points uniformly distributed across full 360° sphere
    - Uses spherical coordinates (azimuth, elevation) for uniform coverage
    - Depth ranges: near (0.5-2m), medium (2-10m), far (10-50m, tapered)
    - Tests UPnP, EPnP, and SQPnP accuracy
    - Includes backward-facing vectors (behind camera)

    Expected results (based on C++ test):
    - UPnP: ~1e-15m (perfect), ~0.01m with noise
    - EPnP: ~1e-11m to 1e-5m (good, varies by architecture)
    - SQPnP: ~0.02m (FAILS even without noise on panoramic)
    """
    print("\n" + "=" * 80)
    print("Testing Panorama 360° Algorithm Comparison")
    print("=" * 80)
    print("\nGenerating points uniformly distributed across full 360° sphere")
    print("Camera frame: ENU (East-North-Up), Identity rotation")
    print("Depth ranges: near (0.5-2m, 25%), medium (2-10m, 50%), far (10-50m tapered, 25%)")
    print(f"NumPy version: {np.__version__}, default dtype: {np.empty(1).dtype}")

    num_points = 200
    num_runs = 20

    # Statistics accumulators
    sqpnp_position_errors = []
    sqpnp_rotation_errors = []
    epnp_position_errors = []
    epnp_rotation_errors = []
    upnp_position_errors = []
    upnp_rotation_errors = []

    sqpnp_times = []
    epnp_times = []
    upnp_times = []

    for run in range(num_runs):
        # Random camera position
        position = generateRandomTranslation(2.0)
        rotation = np.eye(3, dtype=np.float64)  # Identity - camera frame = ENU frame

        # Generate random 3D points uniformly distributed across full 360° sphere
        points = np.empty((num_points, 3), dtype=np.float64)
        bearing_vectors = np.empty((num_points, 3), dtype=np.float64)

        for i in range(num_points):
            # Generate uniform distribution in spherical coordinates
            azimuth = np.random.uniform(0, 2.0 * np.pi)
            elevation = np.random.uniform(-np.pi / 2.0, np.pi / 2.0)

            # Depth distribution: near (0.5-2m), medium (2-10m), far (10-50m with taper)
            rand = np.random.random()
            if rand < 0.25:
                depth = np.random.uniform(0.5, 2.0)  # Near points
            elif rand < 0.75:
                depth = np.random.uniform(2.0, 10.0)  # Medium distance
            else:
                # Far points: quadratic taper from 10m to 50m
                u = np.random.uniform(0, 1)
                depth = 10.0 + (50.0 - 10.0) * (1.0 - u * u)

            # Convert spherical to Cartesian in ENU frame
            cos_elev = np.cos(elevation)
            point_camera = np.array(
                [
                    depth * cos_elev * np.sin(azimuth),  # East
                    depth * cos_elev * np.cos(azimuth),  # North
                    depth * np.sin(elevation),  # Up
                ]
            )

            # Transform from camera frame (ENU) to world frame
            world_point = rotation @ point_camera + position
            points[i] = world_point

            # Bearing vector in camera frame (ENU)
            bearing_vectors[i] = normalized(point_camera)

        # Test SQPnP
        t_start = time.perf_counter()
        sqpnp_result = pyopengv.absolute_pose_sqpnp(bearing_vectors, points)
        sqpnp_time = (time.perf_counter() - t_start) * 1000.0  # ms

        # Test EPnP
        t_start = time.perf_counter()
        epnp_result = pyopengv.absolute_pose_epnp(bearing_vectors, points)
        epnp_time = (time.perf_counter() - t_start) * 1000.0  # ms

        # Test UPnP (returns multiple solutions, pick best)
        t_start = time.perf_counter()
        upnp_results = pyopengv.absolute_pose_upnp(bearing_vectors, points)
        upnp_time = (time.perf_counter() - t_start) * 1000.0  # ms

        # Find best UPnP solution
        best_upnp_error = float("inf")
        upnp_result = None
        for T in upnp_results:
            err = np.linalg.norm(T[:, 3] - position)
            if err < best_upnp_error:
                best_upnp_error = err
                upnp_result = T

        # Calculate errors for SQPnP
        sqpnp_t = sqpnp_result[:, 3]
        sqpnp_R = sqpnp_result[:, :3]
        sqpnp_t_err = np.linalg.norm(sqpnp_t - position)
        sqpnp_R_rel = sqpnp_R.T @ rotation
        trace = np.clip(np.trace(sqpnp_R_rel), -1.0, 3.0)
        sqpnp_r_err = np.arccos((trace - 1.0) / 2.0)

        # Calculate errors for EPnP
        epnp_t = epnp_result[:, 3]
        epnp_R = epnp_result[:, :3]
        epnp_t_err = np.linalg.norm(epnp_t - position)
        epnp_R_rel = epnp_R.T @ rotation
        trace = np.clip(np.trace(epnp_R_rel), -1.0, 3.0)
        epnp_r_err = np.arccos((trace - 1.0) / 2.0)

        # Calculate errors for UPnP
        upnp_t = upnp_result[:, 3]
        upnp_R = upnp_result[:, :3]
        upnp_t_err = np.linalg.norm(upnp_t - position)
        upnp_R_rel = upnp_R.T @ rotation
        trace = np.clip(np.trace(upnp_R_rel), -1.0, 3.0)
        upnp_r_err = np.arccos((trace - 1.0) / 2.0)

        # Accumulate statistics
        sqpnp_position_errors.append(sqpnp_t_err)
        sqpnp_rotation_errors.append(sqpnp_r_err)
        epnp_position_errors.append(epnp_t_err)
        epnp_rotation_errors.append(epnp_r_err)
        upnp_position_errors.append(upnp_t_err)
        upnp_rotation_errors.append(upnp_r_err)

        sqpnp_times.append(sqpnp_time)
        epnp_times.append(epnp_time)
        upnp_times.append(upnp_time)

    # Print results
    print("\n" + "-" * 80)
    print("Results (averaged over {} runs):".format(num_runs))
    print("-" * 80)

    print("\nSQPnP:")
    print(f"  Position error: {np.mean(sqpnp_position_errors):.6e} m")
    print(
        f"  Rotation error: {np.mean(sqpnp_rotation_errors):.6e} rad "
        f"({np.degrees(np.mean(sqpnp_rotation_errors)):.4f}°)"
    )
    print(f"  Time: {np.mean(sqpnp_times):.3f} ms")

    print("\nEPnP:")
    print(f"  Position error: {np.mean(epnp_position_errors):.6e} m")
    print(
        f"  Rotation error: {np.mean(epnp_rotation_errors):.6e} rad "
        f"({np.degrees(np.mean(epnp_rotation_errors)):.4f}°)"
    )
    print(f"  Time: {np.mean(epnp_times):.3f} ms")

    print("\nUPnP:")
    print(f"  Position error: {np.mean(upnp_position_errors):.6e} m")
    print(
        f"  Rotation error: {np.mean(upnp_rotation_errors):.6e} rad "
        f"({np.degrees(np.mean(upnp_rotation_errors)):.4f}°)"
    )
    print(f"  Time: {np.mean(upnp_times):.3f} ms")

    # Validation
    print("\n" + "-" * 80)
    print("Validation:")
    print("-" * 80)

    # UPnP should be nearly perfect (1e-14 level)
    assert (
        np.mean(upnp_position_errors) < 1e-6
    ), f"UPnP position error {np.mean(upnp_position_errors):.6e} too large (expected ~1e-15)"
    assert (
        np.mean(upnp_rotation_errors) < 1e-6
    ), f"UPnP rotation error {np.mean(upnp_rotation_errors):.6e} too large (expected ~1e-16)"
    print("  [PASS] UPnP: Perfect accuracy (< 1e-6)")

    # EPnP should be good but not perfect (~1e-11 level typically, but can be 1e-5 on
    # some architectures)
    # Relax threshold to 1e-4 to account for numerical precision variations across platforms
    assert np.mean(epnp_position_errors) < 1e-4, (
        f"EPnP position error {np.mean(epnp_position_errors):.6e} too large "
        "  (expected ~1e-5 to 1e-11)"
    )
    assert np.mean(epnp_rotation_errors) < 1e-5, (
        f"EPnP rotation error {np.mean(epnp_rotation_errors):.6e} too large "
        "(expected ~1e-6 to 1e-12)"
    )
    print("  [PASS] EPnP: Good accuracy (< 1e-4 m, < 1e-5 rad)")

    # SQPnP typically fails on panoramic views (errors > 0.01m even without noise)
    if np.mean(sqpnp_position_errors) > 0.01:
        print(
            f"  [EXPECTED] SQPnP: Poor accuracy on panoramic ({np.mean(sqpnp_position_errors):.2f}m)"
        )
        print("             Note: SQPnP not designed for 360° views with backward vectors")

    # Calculate relative performance
    if np.mean(epnp_position_errors) > 0:
        upnp_advantage = np.mean(epnp_position_errors) / np.mean(upnp_position_errors)
        print(f"\n  UPnP is {upnp_advantage:.1f}x more accurate than EPnP for panoramic views")

    print("\n" + "=" * 80)
    print("Done testing Panorama 360°")
    print("=" * 80)


if __name__ == "__main__":
    test_panorama_360()
    print("\nPanorama test passed!")
