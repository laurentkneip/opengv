"""
Basic tests for absolute pose estimation (camera pose from 3D-2D correspondences).
"""

import numpy as np
from test_utils import (
    generateRandomPoint,
    generateRandomRotation,
    generateRandomTranslation,
    normalized,
)

import pyopengv


def test_absolute_pose():
    """Test basic absolute pose estimation with various algorithms."""
    print("Testing absolute pose")

    num_points = 50
    position = generateRandomTranslation(2.0)
    rotation = generateRandomRotation(0.5)

    points = np.array([generateRandomPoint(10.0, 4.0) for _ in range(num_points)])
    bearing_vectors = np.array(
        [normalized(rotation.T @ (points[i] - position)) for i in range(num_points)]
    )

    # Test P3P Kneip
    result_kneip = pyopengv.absolute_pose_p3p_kneip(bearing_vectors[:3], points[:3])
    assert len(result_kneip) > 0, "P3P Kneip returned no solutions"
    assert result_kneip[0].shape == (3, 4), f"Expected (3,4), got {result_kneip[0].shape}"

    # Test P3P Gao
    result_gao = pyopengv.absolute_pose_p3p_gao(bearing_vectors[:3], points[:3])
    assert len(result_gao) > 0, "P3P Gao returned no solutions"

    # Test EPnP
    result_epnp = pyopengv.absolute_pose_epnp(bearing_vectors, points)
    assert result_epnp.shape == (3, 4), f"Expected (3,4), got {result_epnp.shape}"

    # Test UPnP
    result_upnp = pyopengv.absolute_pose_upnp(bearing_vectors, points)
    assert len(result_upnp) > 0, "UPnP returned no solutions"
    assert result_upnp[0].shape == (3, 4), f"Expected (3,4), got {result_upnp[0].shape}"

    # Test SQPnP
    result_sqpnp = pyopengv.absolute_pose_sqpnp(bearing_vectors, points)
    assert result_sqpnp.shape == (3, 4), f"Expected (3,4), got {result_sqpnp.shape}"

    print("Done testing absolute pose")


def test_absolute_pose_ransac():
    """Test RANSAC-based absolute pose estimation."""
    print("Testing absolute pose RANSAC")

    num_points = 100
    position = generateRandomTranslation(2.0)
    rotation = generateRandomRotation(0.5)

    points = np.array([generateRandomPoint(10.0, 4.0) for _ in range(num_points)])
    bearing_vectors = np.array(
        [normalized(rotation.T @ (points[i] - position)) for i in range(num_points)]
    )

    # Test RANSAC with KNEIP
    result_kneip = pyopengv.absolute_pose_ransac(bearing_vectors, points, "KNEIP", 0.01, 1000)
    assert result_kneip.shape == (3, 4), f"Expected (3,4), got {result_kneip.shape}"
    print("  RANSAC KNEIP: OK")

    # Test RANSAC with EPNP
    result_epnp = pyopengv.absolute_pose_ransac(bearing_vectors, points, "EPNP", 0.01, 1000)
    assert result_epnp.shape == (3, 4), f"Expected (3,4), got {result_epnp.shape}"
    print("  RANSAC EPNP: OK")

    print("Done testing absolute pose RANSAC")


def test_absolute_pose_lmeds():
    """Test LMedS-based absolute pose estimation."""
    print("Testing absolute pose LMedS")

    num_points = 100
    position = generateRandomTranslation(2.0)
    rotation = generateRandomRotation(0.5)

    points = np.array([generateRandomPoint(10.0, 4.0) for _ in range(num_points)])
    bearing_vectors = np.array(
        [normalized(rotation.T @ (points[i] - position)) for i in range(num_points)]
    )

    result = pyopengv.absolute_pose_lmeds(bearing_vectors, points, "KNEIP", 1000)
    assert result.shape == (3, 4), f"Expected (3,4), got {result.shape}"
    print("  LMedS KNEIP: OK")

    print("Done testing absolute pose LMedS")


if __name__ == "__main__":
    test_absolute_pose()
    test_absolute_pose_ransac()
    test_absolute_pose_lmeds()
    print("\nAll basic absolute pose tests passed!")
