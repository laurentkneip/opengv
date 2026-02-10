"""
Tests for relative pose estimation (camera-to-camera pose).
"""

import numpy as np
from test_utils import (
    generateRandomPoint,
    generateRandomRotation,
    generateRandomTranslation,
    normalized,
)

import pyopengv


def test_relative_pose():
    """Test basic relative pose estimation."""
    print("Testing relative pose")

    num_points = 100

    position1 = np.zeros(3)
    rotation1 = np.eye(3)

    position2 = generateRandomTranslation(2.0)
    rotation2 = generateRandomRotation(0.5)

    points = np.array([generateRandomPoint(10.0, 4.0) for _ in range(num_points)])

    bearing_vectors1 = np.array(
        [normalized(rotation1.T @ (points[i] - position1)) for i in range(num_points)]
    )
    bearing_vectors2 = np.array(
        [normalized(rotation2.T @ (points[i] - position2)) for i in range(num_points)]
    )

    result = pyopengv.relative_pose_eightpt(bearing_vectors1, bearing_vectors2)

    # Verify result shape (eightpt returns essential matrix)
    assert result.shape == (3, 3), f"Expected (3,3) essential matrix, got {result.shape}"

    print("Done testing relative pose")


def test_relative_pose_ransac():
    """Test RANSAC-based relative pose estimation."""
    print("Testing relative pose ransac")

    result = pyopengv.relative_pose_ransac(
        np.random.randn(100, 3), np.random.randn(100, 3), "STEWENIUS", 0.01, 1000
    )

    assert result.shape == (3, 4), f"Expected (3,4), got {result.shape}"

    print("Done testing relative pose ransac")


def test_relative_pose_ransac_rotation_only():
    """Test RANSAC-based rotation-only relative pose estimation."""
    print("Testing relative pose ransac rotation only")

    result = pyopengv.relative_pose_ransac_rotation_only(
        np.random.randn(100, 3), np.random.randn(100, 3), 0.01, 1000
    )

    assert result.shape == (3, 3), f"Expected (3,3), got {result.shape}"

    print("Done testing relative pose ransac rotation only")


if __name__ == "__main__":
    test_relative_pose()
    test_relative_pose_ransac()
    test_relative_pose_ransac_rotation_only()
    print("\nAll relative pose tests passed!")
