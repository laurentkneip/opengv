"""
Lambda Twist P3P Implementation
Based on: "Lambda Twist: An Accurate Fast Robust Perspective Three Point (P3P) Solver"
by Mikael Persson and Klas Nordberg, ECCV 2018

This implementation works directly with bearing vectors (unit rays), making it suitable
for any calibrated camera model including:
- Perspective cameras
- Fisheye/omnidirectional cameras
- Equirectangular (spherical) cameras
"""

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np


@dataclass
class P3PSolution:
    """Holds a single P3P solution"""

    R: np.ndarray  # 3x3 rotation matrix
    t: np.ndarray  # 3x1 translation vector
    lambdas: np.ndarray  # depths λ1, λ2, λ3


def solve_cubic_single_real_root(coeffs: np.ndarray) -> float:
    """
    Find a single real root of cubic: c3*x^3 + c2*x^2 + c1*x + c0 = 0
    Uses Newton-Raphson with a robust initialization heuristic.

    Args:
        coeffs: [c0, c1, c2, c3] polynomial coefficients

    Returns:
        A single real root
    """
    c0, c1, c2, c3 = coeffs

    # Handle degenerate cases
    if abs(c3) < 1e-14:
        # Quadratic case
        if abs(c2) < 1e-14:
            # Linear case
            if abs(c1) < 1e-14:
                return 0.0
            return -c0 / c1
        disc = c1 * c1 - 4 * c2 * c0
        if disc >= 0:
            return (-c1 + np.sqrt(disc)) / (2 * c2)
        return -c1 / (2 * c2)

    # Normalize to monic polynomial: x^3 + a*x^2 + b*x + c = 0
    a = c2 / c3
    b = c1 / c3
    c = c0 / c3

    # Depressed cubic substitution: x = t - a/3
    # t^3 + p*t + q = 0
    p = b - a * a / 3.0
    q = c - a * b / 3.0 + 2.0 * a * a * a / 27.0

    # Discriminant
    disc = q * q / 4.0 + p * p * p / 27.0

    if disc >= 0:
        # One real root (Cardano's formula)
        sqrt_disc = np.sqrt(disc)
        u = np.cbrt(-q / 2.0 + sqrt_disc)
        v = np.cbrt(-q / 2.0 - sqrt_disc)
        t = u + v
    else:
        # Three real roots - use trigonometric method
        # Pick the one with largest magnitude for numerical stability
        m = 2.0 * np.sqrt(-p / 3.0)
        theta = np.arccos(3.0 * q / (p * m)) / 3.0
        t = m * np.cos(theta)

    root = t - a / 3.0

    # Polish with Newton-Raphson
    for _ in range(5):
        f = c0 + root * (c1 + root * (c2 + root * c3))
        df = c1 + root * (2 * c2 + 3 * root * c3)
        if abs(df) < 1e-14:
            break
        delta = f / df
        root -= delta
        if abs(delta) < 1e-12:
            break

    return root


def eig3x3_known_zero(M: np.ndarray) -> Tuple[np.ndarray, float, float]:
    """
    Eigendecomposition of 3x3 symmetric matrix known to have one zero eigenvalue.
    Returns eigenvectors as columns of E, and the two non-zero eigenvalues.

    Implements Algorithm 2 from the paper.

    Args:
        M: 3x3 symmetric matrix with det(M) = 0

    Returns:
        E: 3x3 orthogonal matrix with eigenvectors as columns
        sigma1, sigma2: the two non-zero eigenvalues (sigma1 >= 0, sigma2 <= 0 or both same sign)
    """
    # Find the null space vector (eigenvector for eigenvalue 0)
    # Use cross product of columns
    b3 = np.cross(M[:, 1], M[:, 2])
    norm_b3 = np.linalg.norm(b3)
    if norm_b3 < 1e-10:
        b3 = np.cross(M[:, 0], M[:, 1])
        norm_b3 = np.linalg.norm(b3)
    if norm_b3 < 1e-10:
        b3 = np.cross(M[:, 0], M[:, 2])
        norm_b3 = np.linalg.norm(b3)
    if norm_b3 > 1e-10:
        b3 = b3 / norm_b3
    else:
        b3 = np.array([0, 0, 1.0])

    # Characteristic polynomial for 2 remaining eigenvalues
    # det(M - σI) = 0, but we know one root is 0
    # Reduced to quadratic: σ^2 + p1*σ + p0 = 0
    m = M.flatten()  # row-major
    m1, m2, m3 = m[0], m[1], m[2]
    m4, m5, m6 = m[3], m[4], m[5]
    m7, m8, m9 = m[6], m[7], m[8]

    # Trace and sum of 2x2 minors
    p1 = -(m1 + m5 + m9)
    p0 = m1 * m5 + m1 * m9 + m5 * m9 - m2 * m4 - m3 * m7 - m6 * m8

    # Solve quadratic
    disc = p1 * p1 - 4 * p0
    if disc < 0:
        disc = 0
    sqrt_disc = np.sqrt(disc)
    sigma1 = (-p1 + sqrt_disc) / 2.0
    sigma2 = (-p1 - sqrt_disc) / 2.0

    def get_eigenvector(M, sigma):
        """Get eigenvector for eigenvalue sigma"""
        A = M - sigma * np.eye(3)
        # Find null vector of A
        # Use the row with largest norm for numerical stability
        norms = [np.linalg.norm(A[i]) for i in range(3)]
        idx = np.argsort(norms)[::-1]

        # Cross product of two largest rows
        v = np.cross(A[idx[0]], A[idx[1]])
        norm_v = np.linalg.norm(v)
        if norm_v > 1e-10:
            return v / norm_v

        # Fallback: use SVD
        _, _, Vh = np.linalg.svd(A)
        return Vh[-1]

    b1 = get_eigenvector(M, sigma1)
    b2 = get_eigenvector(M, sigma2)

    # Ensure orthogonality via Gram-Schmidt
    b2 = b2 - np.dot(b2, b1) * b1
    norm_b2 = np.linalg.norm(b2)
    if norm_b2 > 1e-10:
        b2 = b2 / norm_b2
    else:
        # b1 and b2 are parallel, generate orthogonal
        b2 = np.cross(b1, b3)
        b2 = b2 / np.linalg.norm(b2)

    b3 = np.cross(b1, b2)  # Ensure right-handed

    E = np.column_stack([b1, b2, b3])

    # Reorder so |sigma1| >= |sigma2| as per paper
    if abs(sigma1) < abs(sigma2):
        E = np.column_stack([b2, b1, b3])
        sigma1, sigma2 = sigma2, sigma1

    return E, sigma1, sigma2


def refine_lambdas(
    lambdas: np.ndarray,
    M12: np.ndarray,
    M13: np.ndarray,
    M23: np.ndarray,
    a12: float,
    a13: float,
    a23: float,
    max_iters: int = 2,
) -> np.ndarray:
    """
    Gauss-Newton refinement of lambda values.

    Minimizes: (Λ^T M12 Λ - a12)^2 + (Λ^T M13 Λ - a13)^2 + (Λ^T M23 Λ - a23)^2
    """
    L = lambdas.copy()

    for _ in range(max_iters):
        # Residuals
        r12 = L @ M12 @ L - a12
        r13 = L @ M13 @ L - a13
        r23 = L @ M23 @ L - a23

        # Jacobian (derivatives of quadratic forms)
        J12 = 2 * M12 @ L
        J13 = 2 * M13 @ L
        J23 = 2 * M23 @ L

        J = np.vstack([J12, J13, J23])
        r = np.array([r12, r13, r23])

        # Gauss-Newton step
        try:
            delta = np.linalg.lstsq(J, -r, rcond=None)[0]
            L = L + delta
            if np.linalg.norm(delta) < 1e-12:
                break
        except:
            break

    return L


def lambda_twist_p3p(y: np.ndarray, x: np.ndarray, refine: bool = True) -> List[P3PSolution]:
    """
    Lambda Twist P3P Solver

    Solves for camera pose given 3 bearing vectors and their corresponding 3D points.

    Args:
        y: 3x3 array where y[i] is the i-th bearing vector (will be normalized)
        x: 3x3 array where x[i] is the i-th 3D world point
        refine: Whether to apply Gauss-Newton refinement (default True)

    Returns:
        List of P3PSolution objects, each containing R, t, and lambda values
    """
    solutions = []

    # Normalize bearing vectors
    y = np.array(y, dtype=np.float64)
    x = np.array(x, dtype=np.float64)

    y_norm = np.zeros_like(y)
    for i in range(3):
        norm = np.linalg.norm(y[i])
        if norm < 1e-10:
            return solutions  # Degenerate
        y_norm[i] = y[i] / norm

    # Compute a_ij (squared distances between 3D points)
    a12 = np.sum((x[0] - x[1]) ** 2)
    a13 = np.sum((x[0] - x[2]) ** 2)
    a23 = np.sum((x[1] - x[2]) ** 2)

    # Check for collinear 3D points
    if a12 < 1e-10 or a13 < 1e-10 or a23 < 1e-10:
        return solutions  # Degenerate

    # Compute b_ij (cosines of angles between bearing vectors)
    b12 = np.dot(y_norm[0], y_norm[1])
    b13 = np.dot(y_norm[0], y_norm[2])
    b23 = np.dot(y_norm[1], y_norm[2])

    # Construct M matrices (Eq. 4)
    M12 = np.array([[1, -b12, 0], [-b12, 1, 0], [0, 0, 0]], dtype=np.float64)

    M13 = np.array([[1, 0, -b13], [0, 0, 0], [-b13, 0, 1]], dtype=np.float64)

    M23 = np.array([[0, 0, 0], [0, 1, -b23], [0, -b23, 1]], dtype=np.float64)

    # Construct D1 and D2 (Eq. 5-6)
    D1 = M12 * a23 - M23 * a12
    D2 = M13 * a23 - M23 * a13

    # Compute cubic coefficients (Eq. 10)
    # det(D1 + γ*D2) = c3*γ^3 + c2*γ^2 + c1*γ + c0 = 0
    # Using multilinearity of determinant
    d10, d11, d12 = D1[:, 0], D1[:, 1], D1[:, 2]  # columns of D1
    d20, d21, d22 = D2[:, 0], D2[:, 1], D2[:, 2]  # columns of D2

    def det3(a, b, c):
        return np.dot(a, np.cross(b, c))

    # c0: all columns from D1
    c0 = det3(d10, d11, d12)

    # c1: one column from D2, two from D1
    c1 = det3(d20, d11, d12) + det3(d10, d21, d12) + det3(d10, d11, d22)

    # c2: two columns from D2, one from D1
    c2 = det3(d20, d21, d12) + det3(d20, d11, d22) + det3(d10, d21, d22)

    # c3: all columns from D2
    c3 = det3(d20, d21, d22)

    # Handle special cases where D1 or D2 is already rank 2
    if abs(c0) < 1e-10:
        gamma = 0.0
    elif abs(c3) < 1e-10:
        # Use inverse polynomial
        gamma = 1.0 / solve_cubic_single_real_root(np.array([c3, c2, c1, c0]))
    else:
        gamma = solve_cubic_single_real_root(np.array([c0, c1, c2, c3]))

    # D0 = D1 + γ*D2 (Eq. 7)
    D0 = D1 + gamma * D2

    # Eigendecomposition (Eq. 11-12)
    E, sigma1, sigma2 = eig3x3_known_zero(D0)

    # Handle case where both eigenvalues have same sign (shouldn't happen for valid P3P)
    if sigma1 * sigma2 > 0:
        # Try to find another root
        return solutions

    # Compute s values (Eq. 9)
    if abs(sigma1) < 1e-14:
        return solutions

    s_squared = -sigma2 / sigma1
    if s_squared < 0:
        return solutions

    s_values = [np.sqrt(s_squared), -np.sqrt(s_squared)]

    e = E.flatten()  # Column-major to match paper notation
    e0, e3, e6 = E[0, 0], E[1, 0], E[2, 0]  # First eigenvector
    e1, e4, e7 = E[0, 1], E[1, 1], E[2, 1]  # Second eigenvector

    tau_solutions = []

    for s in s_values:
        # Compute w0, w1 (Eq. 13)
        denom = s * e1 - e0
        if abs(denom) < 1e-14:
            continue

        w0 = (e3 - s * e4) / denom
        w1 = (e6 - s * e7) / denom

        # Quadratic coefficients for τ (Eq. 15)
        # Actually, let's compute directly from D1 constraint
        # λ₁ = (w₀ + w₁*τ)*λ₂, λ₃ = τ*λ₂
        # Substitute into Λᵀ*D₁*Λ = 0

        L0 = np.array([w0, 1.0, 0.0])
        L1 = np.array([w1, 0.0, 1.0])

        a_coef = L1 @ D1 @ L1
        b_coef = 2 * L0 @ D1 @ L1
        c_coef = L0 @ D1 @ L0

        # Solve quadratic
        if abs(a_coef) < 1e-14:
            if abs(b_coef) < 1e-14:
                continue
            taus = [-c_coef / b_coef]
        else:
            disc = b_coef * b_coef - 4 * a_coef * c_coef
            if disc < 0:
                continue
            sqrt_disc = np.sqrt(disc)
            taus = [(-b_coef + sqrt_disc) / (2 * a_coef), (-b_coef - sqrt_disc) / (2 * a_coef)]

        for tau in taus:
            if tau <= 0:  # τ must be positive (geometric feasibility)
                continue
            tau_solutions.append((tau, w0, w1))

    # For each valid τ, compute full Λ
    for tau, w0, w1 in tau_solutions:
        # λ₂ from Eq. 16: λ₂² * [1,τ]ᵀ M₂₃[2x2] [1,τ] = a₂₃
        # M₂₃[2x2] = [[1, -b₂₃], [-b₂₃, 1]]
        # Expanding: 1 - 2*τ*b₂₃ + τ²
        denom = 1.0 + tau * tau - 2.0 * tau * b23
        if denom <= 1e-10:
            continue

        lambda2_sq = a23 / denom
        if lambda2_sq <= 0:
            continue

        lambda2 = np.sqrt(lambda2_sq)
        lambda3 = tau * lambda2
        lambda1 = (w0 + w1 * tau) * lambda2

        # Check geometric feasibility
        if lambda1 <= 0:
            continue

        lambdas = np.array([lambda1, lambda2, lambda3])

        # Optional refinement
        if refine:
            lambdas = refine_lambdas(lambdas, M12, M13, M23, a12, a13, a23)
            if np.any(lambdas <= 0):
                continue

        # Recover R and t (Eq. 17-22)
        z1 = lambdas[0] * y_norm[0] - lambdas[1] * y_norm[1]
        z2 = lambdas[1] * y_norm[1] - lambdas[2] * y_norm[2]

        Y = np.column_stack([z1, z2, np.cross(z1, z2)])

        dx12 = x[0] - x[1]
        dx23 = x[1] - x[2]
        X = np.column_stack([dx12, dx23, np.cross(dx12, dx23)])

        # R = Y * X^(-1)
        try:
            X_inv = np.linalg.inv(X)
        except np.linalg.LinAlgError:
            continue

        R = Y @ X_inv

        # Ensure R is a proper rotation matrix
        U, _, Vt = np.linalg.svd(R)
        R = U @ Vt
        if np.linalg.det(R) < 0:
            R = U @ np.diag([1, 1, -1]) @ Vt

        # t = λ₁*y₁ - R*x₁
        t = lambdas[0] * y_norm[0] - R @ x[0]

        solutions.append(P3PSolution(R=R, t=t, lambdas=lambdas))

    return solutions
