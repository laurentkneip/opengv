#!/usr/bin/env python
"""
Master test runner for OpenGV Python bindings.

This runs all modular test files in sequence.
"""
import sys

import test_absolute_pose_basic
import test_panorama

# Import all test modules
import test_relative_pose
import test_triangulation


def run_all_tests():
    """Run all test suites."""
    print("=" * 80)
    print("Running All OpenGV Python Tests")
    print("=" * 80)

    failures = []

    # Test 1: Relative Pose
    print("\n" + "=" * 80)
    print("TEST SUITE 1: Relative Pose")
    print("=" * 80)
    try:
        test_relative_pose.test_relative_pose()
        test_relative_pose.test_relative_pose_ransac()
        test_relative_pose.test_relative_pose_ransac_rotation_only()
        print("[PASS] Relative pose tests")
    except Exception as e:
        print(f"[FAIL] Relative pose tests - {e}")
        failures.append(("Relative Pose", str(e)))

    # Test 2: Triangulation
    print("\n" + "=" * 80)
    print("TEST SUITE 2: Triangulation")
    print("=" * 80)
    try:
        test_triangulation.test_triangulation()
        print("[PASS] Triangulation tests")
    except Exception as e:
        print(f"[FAIL] Triangulation tests - {e}")
        failures.append(("Triangulation", str(e)))

    # Test 3: Absolute Pose (Basic)
    print("\n" + "=" * 80)
    print("TEST SUITE 3: Absolute Pose (Basic)")
    print("=" * 80)
    try:
        test_absolute_pose_basic.test_absolute_pose()
        test_absolute_pose_basic.test_absolute_pose_ransac()
        test_absolute_pose_basic.test_absolute_pose_lmeds()
        print("[PASS] Absolute pose basic tests")
    except Exception as e:
        print(f"[FAIL] Absolute pose basic tests - {e}")
        failures.append(("Absolute Pose Basic", str(e)))

    # Test 4: Panorama 360°
    print("\n" + "=" * 80)
    print("TEST SUITE 4: Panorama 360°")
    print("=" * 80)
    try:
        test_panorama.test_panorama_360()
        print("[PASS] Panorama tests")
    except Exception as e:
        print(f"[FAIL] Panorama tests - {e}")
        failures.append(("Panorama", str(e)))

    # Summary
    print("\n" + "=" * 80)
    print("TEST SUMMARY")
    print("=" * 80)

    if not failures:
        print("[PASS] ALL TESTS PASSED!")
        return 0
    else:
        print(f"[FAIL] {len(failures)} test suite(s) failed:")
        for name, error in failures:
            print(f"  - {name}: {error}")
        return 1


if __name__ == "__main__":
    sys.exit(run_all_tests())
