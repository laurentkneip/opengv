"""
Utility functions for OpenGV Python tests.
"""

import numpy as np


def calculate_ransac_threshold(pixel_error_std, focal_length):
    """
    Calculate RANSAC threshold for absolute pose estimation.

    Args:
        pixel_error_std: Standard deviation of pixel error (e.g., 0.5 pixels)
        focal_length: Focal length in pixels (e.g., 800 for standard camera)

    Returns:
        threshold: Angular threshold for RANSAC (1 - cos(angle))

    The formula computes the angular error corresponding to pixel error:
    - sqrt(2.0) * pixel_error accounts for 2D pixel error in both x and y
    - arctan(error / focal_length) converts to angular error
    - 1 - cos(angle) is the RANSAC threshold metric
    """
    angular_error = np.arctan(np.sqrt(2.0) * pixel_error_std / focal_length)
    threshold = 1.0 - np.cos(angular_error)
    return threshold


def normalized(x):
    """Normalize a vector to unit length."""
    return x / np.linalg.norm(x)


def generateRandomPoint(maximumDepth, minimumDepth):
    """Generate a random 3D point within depth range."""
    cleanPoint = np.random.uniform(-1.0, 1.0, 3)
    direction = normalized(cleanPoint)
    return (maximumDepth - minimumDepth) * cleanPoint + minimumDepth * direction


def generateRandomTranslation(maximumParallax):
    """Generate a random translation vector."""
    return np.random.uniform(-maximumParallax, maximumParallax, 3)


def generateRandomRotation(maxAngle):
    """Generate a random rotation matrix from roll-pitch-yaw angles."""
    rpy = np.random.uniform(-maxAngle, maxAngle, 3)

    R1 = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, np.cos(rpy[0]), -np.sin(rpy[0])],
            [0.0, np.sin(rpy[0]), np.cos(rpy[0])],
        ]
    )

    R2 = np.array(
        [
            [np.cos(rpy[1]), 0.0, np.sin(rpy[1])],
            [0.0, 1.0, 0.0],
            [-np.sin(rpy[1]), 0.0, np.cos(rpy[1])],
        ]
    )

    R3 = np.array(
        [
            [np.cos(rpy[2]), -np.sin(rpy[2]), 0.0],
            [np.sin(rpy[2]), np.cos(rpy[2]), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )

    return R3 @ R2 @ R1


def addNoise(std_dev, cleanValue):
    """Add Gaussian noise to a value."""
    return cleanValue + np.random.normal(0.0, std_dev)
