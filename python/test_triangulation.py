"""
Tests for triangulation (3D point reconstruction from multiple views).
"""

import numpy as np
from test_utils import (
    generateRandomPoint,
    generateRandomRotation,
    generateRandomTranslation,
    normalized,
)

import pyopengv


def test_triangulation():
    """Test triangulation of 3D points from two views."""
    print("Testing triangulation")

    num_points = 100

    # First camera at origin
    position1 = np.zeros(3)
    rotation1 = np.eye(3)

    # Second camera with random pose
    position2 = generateRandomTranslation(2.0)
    rotation2 = generateRandomRotation(0.5)

    # Generate random 3D points
    points = np.array([generateRandomPoint(10.0, 4.0) for _ in range(num_points)])

    # Compute bearing vectors from both cameras
    bearing_vectors1 = np.array(
        [normalized(rotation1.T @ (points[i] - position1)) for i in range(num_points)]
    )
    bearing_vectors2 = np.array(
        [normalized(rotation2.T @ (points[i] - position2)) for i in range(num_points)]
    )

    # Triangulate points
    triangulated = pyopengv.triangulation_triangulate(
        bearing_vectors1, bearing_vectors2, position2, rotation2
    )

    # Verify shape
    assert triangulated.shape == (
        num_points,
        3,
    ), f"Expected ({num_points}, 3), got {triangulated.shape}"

    # Verify triangulation accuracy
    errors = np.linalg.norm(triangulated - points, axis=1)
    mean_error = np.mean(errors)

    assert mean_error < 1e-6, f"Triangulation error {mean_error:.2e} exceeds threshold 1e-6"

    print(f"  Mean triangulation error: {mean_error:.2e}")
    print("Done testing triangulation")


if __name__ == "__main__":
    test_triangulation()
    print("\nTriangulation test passed!")
