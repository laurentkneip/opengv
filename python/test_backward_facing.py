#!/usr/bin/env python3
"""
Test backward-facing points for 360° cameras.

This test validates the hypothesis that EPnP struggles with backward-facing
points (negative Z in camera frame) while UPnP handles them correctly.

Test setup:
- All 3D points are behind the camera (negative Z / backward hemisphere)
- Compare RANSAC with GAO/KNEIP minimal solvers
- Compare EPnP vs UPnP refinement on the inliers
"""

import numpy as np
from scipy.spatial.transform import Rotation

import pyopengv


def rotation_matrix_from_euler(roll: float, pitch: float, yaw: float) -> np.ndarray:
    """Create rotation matrix from Euler angles (radians)."""
    return Rotation.from_euler("xyz", [roll, pitch, yaw]).as_matrix()


def rotation_error(R1: np.ndarray, R2: np.ndarray) -> float:
    """Compute rotation error in radians."""
    R_diff = R1.T @ R2
    trace = np.trace(R_diff)
    trace = np.clip(trace, -1.0, 3.0)
    return np.arccos((trace - 1.0) / 2.0)


def generate_backward_points(
    num_points: int, min_depth: float = 2.0, max_depth: float = 20.0
) -> np.ndarray:
    """
    Generate 3D points in the backward hemisphere (negative Y in ENU frame).

    In ENU (East-North-Up) convention used by OpenGV:
    - X = East (right)
    - Y = North (forward)
    - Z = Up

    Backward = negative Y hemisphere
    """
    points = []
    for _ in range(num_points):
        # Generate random direction in backward hemisphere
        # Azimuth: 90° to 270° (backward facing)
        azimuth = np.random.uniform(np.pi / 2, 3 * np.pi / 2)
        # Elevation: -60° to +60° (avoid poles)
        elevation = np.random.uniform(-np.pi / 3, np.pi / 3)

        depth = np.random.uniform(min_depth, max_depth)

        # Convert to ENU Cartesian
        cos_elev = np.cos(elevation)
        x = depth * cos_elev * np.sin(azimuth)  # East
        y = depth * cos_elev * np.cos(azimuth)  # North (negative for backward)
        z = depth * np.sin(elevation)  # Up

        points.append([x, y, z])

    return np.array(points)


def generate_forward_points(
    num_points: int, min_depth: float = 2.0, max_depth: float = 20.0
) -> np.ndarray:
    """
    Generate 3D points in the forward hemisphere (positive Y in ENU frame).
    """
    points = []
    for _ in range(num_points):
        # Azimuth: -90° to +90° (forward facing)
        azimuth = np.random.uniform(-np.pi / 2, np.pi / 2)
        elevation = np.random.uniform(-np.pi / 3, np.pi / 3)

        depth = np.random.uniform(min_depth, max_depth)

        cos_elev = np.cos(elevation)
        x = depth * cos_elev * np.sin(azimuth)
        y = depth * cos_elev * np.cos(azimuth)  # Positive for forward
        z = depth * np.sin(elevation)

        points.append([x, y, z])

    return np.array(points)


def add_pixel_noise(
    bearings: np.ndarray, noise_std_pixels: float, focal_length: float = 800.0
) -> np.ndarray:
    """Add pixel-level noise to bearing vectors."""
    noisy_bearings = []
    for b in bearings:
        # Convert to angles
        theta = np.arctan2(np.sqrt(b[0] ** 2 + b[1] ** 2), b[2])
        phi = np.arctan2(b[1], b[0])

        # Add angular noise (pixel noise / focal length)
        noise_rad = noise_std_pixels / focal_length
        theta += np.random.normal(0, noise_rad)
        phi += np.random.normal(0, noise_rad)

        # Convert back to unit vector
        noisy_b = np.array(
            [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
        )
        noisy_bearings.append(noisy_b / np.linalg.norm(noisy_b))

    return np.array(noisy_bearings)


def run_test(
    points_world: np.ndarray,
    R_true: np.ndarray,
    t_true: np.ndarray,
    noise_pixels: float = 0.0,
    test_name: str = "Test",
) -> dict:
    """
    Run pose estimation test comparing EPnP vs UPnP refinement.

    Returns dict with position and rotation errors for each method.
    """
    # Transform points to camera frame
    points_camera = (R_true.T @ (points_world.T - t_true.reshape(3, 1))).T

    # Generate bearing vectors (normalized directions)
    bearings = points_camera / np.linalg.norm(points_camera, axis=1, keepdims=True)

    # Add noise if specified
    if noise_pixels > 0:
        bearings = add_pixel_noise(bearings, noise_pixels)

    results = {}

    # Test EPnP directly
    try:
        T = pyopengv.absolute_pose_epnp(bearings, points_world)

        R_est = T[:3, :3]
        t_est = T[:3, 3]

        pos_err = np.linalg.norm(t_est - t_true)
        rot_err = rotation_error(R_est, R_true)

        results["EPnP"] = {
            "position_error": pos_err,
            "rotation_error": rot_err,
            "rotation_error_deg": np.degrees(rot_err),
            "success": True,
        }
    except Exception as e:
        results["EPnP"] = {
            "position_error": float("inf"),
            "rotation_error": float("inf"),
            "rotation_error_deg": float("inf"),
            "success": False,
            "error": str(e),
        }

    # Test UPnP - returns list of solutions, pick best one
    try:
        solutions = pyopengv.absolute_pose_upnp(bearings, points_world)

        best_err = float("inf")
        best_T = None
        for T in solutions:
            T = np.array(T)
            R_est = T[:3, :3]
            t_est = T[:3, 3]
            err = np.linalg.norm(t_est - t_true)
            if err < best_err:
                best_err = err
                best_T = T

        if best_T is not None:
            R_est = best_T[:3, :3]
            t_est = best_T[:3, 3]
            pos_err = np.linalg.norm(t_est - t_true)
            rot_err = rotation_error(R_est, R_true)

            results["UPnP"] = {
                "position_error": pos_err,
                "rotation_error": rot_err,
                "rotation_error_deg": np.degrees(rot_err),
                "success": True,
            }
        else:
            results["UPnP"] = {
                "position_error": float("inf"),
                "rotation_error": float("inf"),
                "rotation_error_deg": float("inf"),
                "success": False,
                "error": "No solutions",
            }
    except Exception as e:
        results["UPnP"] = {
            "position_error": float("inf"),
            "rotation_error": float("inf"),
            "rotation_error_deg": float("inf"),
            "success": False,
            "error": str(e),
        }

    # Test RANSAC with different minimal solvers + refinement algorithms
    # This is the proper comparison: same RANSAC, different refinement
    ransac_configs = [
        ("KNEIP", "EPNP", "KNEIP→EPnP"),
        ("KNEIP", "UPNP", "KNEIP→UPnP"),
        ("GAO", "EPNP", "GAO→EPnP"),
        ("GAO", "UPNP", "GAO→UPnP"),
    ]

    for minimal, refine, name in ransac_configs:
        try:
            T = pyopengv.absolute_pose_ransac_with_refine(
                bearings,
                points_world,
                minimal,
                refine,
                threshold=1e-5,
                iterations=1000,
                probability=0.99,
            )

            R_est = T[:3, :3]
            t_est = T[:3, 3]

            pos_err = np.linalg.norm(t_est - t_true)
            rot_err = rotation_error(R_est, R_true)

            results[name] = {
                "position_error": pos_err,
                "rotation_error": rot_err,
                "rotation_error_deg": np.degrees(rot_err),
                "success": True,
            }
        except Exception as e:
            results[name] = {
                "position_error": float("inf"),
                "rotation_error": float("inf"),
                "rotation_error_deg": float("inf"),
                "success": False,
                "error": str(e),
            }

    return results


def print_results(results: dict, test_name: str):
    """Pretty print test results."""
    print(f"\n{'='*60}")
    print(f"  {test_name}")
    print(f"{'='*60}")
    print(f"{'Method':<20} {'Pos Error (m)':<15} {'Rot Error (deg)':<15}")
    print(f"{'-'*50}")

    for method, data in results.items():
        if data["success"]:
            print(
                f"{method:<20} {data['position_error']:<15.6f} {data['rotation_error_deg']:<15.4f}"
            )
        else:
            print(f"{method:<20} {'FAILED':<15} {'FAILED':<15}")


def main():
    """Run backward-facing vs forward-facing comparison tests."""
    np.random.seed(42)

    # Ground truth pose
    R_true = rotation_matrix_from_euler(0.1, 0.2, 0.3)
    t_true = np.array([1.0, 2.0, 0.5])

    num_points = 50
    num_trials = 10

    print("\n" + "=" * 70)
    print("  360° BACKWARD-FACING POINTS TEST")
    print("  Comparing EPnP vs UPnP refinement after RANSAC")
    print("=" * 70)

    # Aggregate results
    agg_results = {"backward": {}, "forward": {}, "backward_noisy": {}, "forward_noisy": {}}

    for trial in range(num_trials):
        # Test 1: Backward-facing points, no noise
        points_back = generate_backward_points(num_points)
        results_back = run_test(
            points_back, R_true, t_true, noise_pixels=0.0, test_name="Backward (no noise)"
        )

        # Test 2: Forward-facing points, no noise
        points_fwd = generate_forward_points(num_points)
        results_fwd = run_test(
            points_fwd, R_true, t_true, noise_pixels=0.0, test_name="Forward (no noise)"
        )

        # Test 3: Backward-facing with noise
        results_back_noisy = run_test(
            points_back, R_true, t_true, noise_pixels=2.0, test_name="Backward (2px noise)"
        )

        # Test 4: Forward-facing with noise
        results_fwd_noisy = run_test(
            points_fwd, R_true, t_true, noise_pixels=2.0, test_name="Forward (2px noise)"
        )

        # Accumulate
        for method in results_back:
            if method not in agg_results["backward"]:
                agg_results["backward"][method] = {"pos": [], "rot": []}
                agg_results["forward"][method] = {"pos": [], "rot": []}
                agg_results["backward_noisy"][method] = {"pos": [], "rot": []}
                agg_results["forward_noisy"][method] = {"pos": [], "rot": []}

            if results_back[method]["success"]:
                agg_results["backward"][method]["pos"].append(
                    results_back[method]["position_error"]
                )
                agg_results["backward"][method]["rot"].append(
                    results_back[method]["rotation_error_deg"]
                )

            if results_fwd[method]["success"]:
                agg_results["forward"][method]["pos"].append(results_fwd[method]["position_error"])
                agg_results["forward"][method]["rot"].append(
                    results_fwd[method]["rotation_error_deg"]
                )

            if results_back_noisy[method]["success"]:
                agg_results["backward_noisy"][method]["pos"].append(
                    results_back_noisy[method]["position_error"]
                )
                agg_results["backward_noisy"][method]["rot"].append(
                    results_back_noisy[method]["rotation_error_deg"]
                )

            if results_fwd_noisy[method]["success"]:
                agg_results["forward_noisy"][method]["pos"].append(
                    results_fwd_noisy[method]["position_error"]
                )
                agg_results["forward_noisy"][method]["rot"].append(
                    results_fwd_noisy[method]["rotation_error_deg"]
                )

    # Print summary statistics
    print("\n" + "=" * 70)
    print("  SUMMARY (mean ± std over {} trials)".format(num_trials))
    print("=" * 70)

    for scenario, label in [
        ("backward", "BACKWARD-FACING (no noise)"),
        ("forward", "FORWARD-FACING (no noise)"),
        ("backward_noisy", "BACKWARD-FACING (2px noise)"),
        ("forward_noisy", "FORWARD-FACING (2px noise)"),
    ]:
        print(f"\n{label}:")
        print(f"{'Method':<20} {'Pos Error (m)':<20} {'Rot Error (deg)':<20}")
        print("-" * 60)

        for method in agg_results[scenario]:
            pos_vals = agg_results[scenario][method]["pos"]
            rot_vals = agg_results[scenario][method]["rot"]

            if pos_vals:
                pos_mean, pos_std = np.mean(pos_vals), np.std(pos_vals)
                rot_mean, rot_std = np.mean(rot_vals), np.std(rot_vals)
                print(
                    f"{method:<20} {pos_mean:.6f} ± {pos_std:.6f}    {rot_mean:.4f} ± {rot_std:.4f}"
                )
            else:
                print(f"{method:<20} {'ALL FAILED':<20} {'ALL FAILED':<20}")

    # Conclusion
    print("\n" + "=" * 70)
    print("  ANALYSIS")
    print("=" * 70)

    # Compare EPnP vs UPnP refinement on backward-facing data (same RANSAC, different refinement)
    back_kneip_epnp_pos = (
        np.mean(agg_results["backward_noisy"]["KNEIP→EPnP"]["pos"])
        if agg_results["backward_noisy"].get("KNEIP→EPnP", {}).get("pos")
        else float("inf")
    )
    back_kneip_upnp_pos = (
        np.mean(agg_results["backward_noisy"]["KNEIP→UPnP"]["pos"])
        if agg_results["backward_noisy"].get("KNEIP→UPnP", {}).get("pos")
        else float("inf")
    )
    fwd_kneip_epnp_pos = (
        np.mean(agg_results["forward_noisy"]["KNEIP→EPnP"]["pos"])
        if agg_results["forward_noisy"].get("KNEIP→EPnP", {}).get("pos")
        else float("inf")
    )
    fwd_kneip_upnp_pos = (
        np.mean(agg_results["forward_noisy"]["KNEIP→UPnP"]["pos"])
        if agg_results["forward_noisy"].get("KNEIP→UPnP", {}).get("pos")
        else float("inf")
    )

    back_kneip_epnp_rot = (
        np.mean(agg_results["backward_noisy"]["KNEIP→EPnP"]["rot"])
        if agg_results["backward_noisy"].get("KNEIP→EPnP", {}).get("rot")
        else float("inf")
    )
    back_kneip_upnp_rot = (
        np.mean(agg_results["backward_noisy"]["KNEIP→UPnP"]["rot"])
        if agg_results["backward_noisy"].get("KNEIP→UPnP", {}).get("rot")
        else float("inf")
    )
    fwd_kneip_epnp_rot = (
        np.mean(agg_results["forward_noisy"]["KNEIP→EPnP"]["rot"])
        if agg_results["forward_noisy"].get("KNEIP→EPnP", {}).get("rot")
        else float("inf")
    )
    fwd_kneip_upnp_rot = (
        np.mean(agg_results["forward_noisy"]["KNEIP→UPnP"]["rot"])
        if agg_results["forward_noisy"].get("KNEIP→UPnP", {}).get("rot")
        else float("inf")
    )

    # Direct algorithm results
    back_epnp_pos = (
        np.mean(agg_results["backward_noisy"]["EPnP"]["pos"])
        if agg_results["backward_noisy"].get("EPnP", {}).get("pos")
        else float("inf")
    )
    back_upnp_pos = (
        np.mean(agg_results["backward_noisy"]["UPnP"]["pos"])
        if agg_results["backward_noisy"].get("UPnP", {}).get("pos")
        else float("inf")
    )
    fwd_epnp_pos = (
        np.mean(agg_results["forward_noisy"]["EPnP"]["pos"])
        if agg_results["forward_noisy"].get("EPnP", {}).get("pos")
        else float("inf")
    )
    fwd_upnp_pos = (
        np.mean(agg_results["forward_noisy"]["UPnP"]["pos"])
        if agg_results["forward_noisy"].get("UPnP", {}).get("pos")
        else float("inf")
    )

    back_epnp_rot = (
        np.mean(agg_results["backward_noisy"]["EPnP"]["rot"])
        if agg_results["backward_noisy"].get("EPnP", {}).get("rot")
        else float("inf")
    )
    back_upnp_rot = (
        np.mean(agg_results["backward_noisy"]["UPnP"]["rot"])
        if agg_results["backward_noisy"].get("UPnP", {}).get("rot")
        else float("inf")
    )
    fwd_epnp_rot = (
        np.mean(agg_results["forward_noisy"]["EPnP"]["rot"])
        if agg_results["forward_noisy"].get("EPnP", {}).get("rot")
        else float("inf")
    )
    fwd_upnp_rot = (
        np.mean(agg_results["forward_noisy"]["UPnP"]["rot"])
        if agg_results["forward_noisy"].get("UPnP", {}).get("rot")
        else float("inf")
    )

    print("\n" + "-" * 70)
    print("  DIRECT ALGORITHMS (no RANSAC, all points)")
    print("-" * 70)

    print("\nBackward-facing with noise:")
    print(f"  EPnP:  pos={back_epnp_pos:.6f} m,  rot={back_epnp_rot:.4f}°")
    print(f"  UPnP:  pos={back_upnp_pos:.6f} m,  rot={back_upnp_rot:.4f}°")
    if back_upnp_pos > 0 and back_upnp_rot > 0:
        print(
            f"  Ratio (EPnP/UPnP):  pos={back_epnp_pos/back_upnp_pos:.2f}x,  "
            f"rot={back_epnp_rot/back_upnp_rot:.2f}x"
        )

    print("\nForward-facing with noise:")
    print(f"  EPnP:  pos={fwd_epnp_pos:.6f} m,  rot={fwd_epnp_rot:.4f}°")
    print(f"  UPnP:  pos={fwd_upnp_pos:.6f} m,  rot={fwd_upnp_rot:.4f}°")
    if fwd_upnp_pos > 0 and fwd_upnp_rot > 0:
        print(
            f"  Ratio (EPnP/UPnP):  pos={fwd_epnp_pos/fwd_upnp_pos:.2f}x,  "
            f"rot={fwd_epnp_rot/fwd_upnp_rot:.2f}x"
        )

    print("\n" + "-" * 70)
    print("  RANSAC + REFINEMENT (KNEIP P3P → refinement)")
    print("-" * 70)

    print(f"\nBackward-facing with noise:")
    print(f"  KNEIP→EPnP:  pos={back_kneip_epnp_pos:.6f} m,  rot={back_kneip_epnp_rot:.4f}°")
    print(f"  KNEIP→UPnP:  pos={back_kneip_upnp_pos:.6f} m,  rot={back_kneip_upnp_rot:.4f}°")
    if back_kneip_upnp_pos > 0 and back_kneip_upnp_rot > 0:
        print(
            f"  Ratio (EPnP/UPnP):  pos={back_kneip_epnp_pos/back_kneip_upnp_pos:.2f}x,  rot={back_kneip_epnp_rot/back_kneip_upnp_rot:.2f}x"
        )

    print(f"\nForward-facing with noise:")
    print(f"  KNEIP→EPnP:  pos={fwd_kneip_epnp_pos:.6f} m,  rot={fwd_kneip_epnp_rot:.4f}°")
    print(f"  KNEIP→UPnP:  pos={fwd_kneip_upnp_pos:.6f} m,  rot={fwd_kneip_upnp_rot:.4f}°")
    if fwd_kneip_upnp_pos > 0 and fwd_kneip_upnp_rot > 0:
        print(
            f"  Ratio (EPnP/UPnP):  pos={fwd_kneip_epnp_pos/fwd_kneip_upnp_pos:.2f}x,  rot={fwd_kneip_epnp_rot/fwd_kneip_upnp_rot:.2f}x"
        )

    # Conclusions
    print("\n" + "=" * 70)
    print("  CONCLUSIONS")
    print("=" * 70)

    if back_epnp_pos > 2 * back_upnp_pos or back_epnp_rot > 2 * back_upnp_rot:
        print("\n⚠️  Direct EPnP is significantly worse than UPnP on all points!")
        print(f"   Position: EPnP {back_epnp_pos/back_upnp_pos:.1f}x worse")
        print(f"   Rotation: EPnP {back_epnp_rot/back_upnp_rot:.1f}x worse")

    if back_epnp_pos > 1.5 * fwd_epnp_pos:
        print("\n⚠️  EPnP is worse on backward-facing than forward-facing!")
        print(f"   Backward: {back_epnp_pos:.4f} m, {back_epnp_rot:.4f}°")
        print(f"   Forward:  {fwd_epnp_pos:.4f} m, {fwd_epnp_rot:.4f}°")

    if abs(back_kneip_epnp_pos - back_kneip_upnp_pos) / max(back_kneip_epnp_pos, 1e-9) < 0.3:
        print("\n✓  With RANSAC, EPnP and UPnP refinement perform similarly.")
        print("   RANSAC filters out problematic point configurations.")


if __name__ == "__main__":
    main()
