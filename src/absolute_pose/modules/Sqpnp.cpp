/******************************************************************************
 * Author:   Laurent Kneip                                                    *
 * Contact:  kneip.laurent@gmail.com                                          *
 * License:  Copyright (c) 2013 Laurent Kneip, ANU. All rights reserved.      *
 *                                                                            *
 * Redistribution and use in source and binary forms, with or without         *
 * modification, are permitted provided that the following conditions         *
 * are met:                                                                   *
 * * Redistributions of source code must retain the above copyright           *
 *   notice, this list of conditions and the following disclaimer.            *
 * * Redistributions in binary form must reproduce the above copyright        *
 *   notice, this list of conditions and the following disclaimer in the      *
 *   documentation and/or other materials provided with the distribution.     *
 * * Neither the name of ANU nor the names of its contributors may be         *
 *   used to endorse or promote products derived from this software without   *
 *   specific prior written permission.                                       *
 *                                                                            *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"*
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE *
 * ARE DISCLAIMED. IN NO EVENT SHALL ANU OR THE CONTRIBUTORS BE LIABLE        *
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL *
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR *
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT         *
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY  *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF     *
 * SUCH DAMAGE.                                                               *
 ******************************************************************************/

// SQPnP implementation adapted for bearing vectors and omnidirectional cameras
// Based on EPnP but modified to use angular error instead of reprojection error
// and to support backward-facing bearing vectors


#include <iostream>
#include <cmath>
using namespace std;

#include <opengv/absolute_pose/modules/Sqpnp.hpp>


opengv::absolute_pose::modules::Sqpnp::Sqpnp(void)
{
  maximum_number_of_correspondences = 0;
  number_of_correspondences = 0;

  pws = 0;
  us = 0;
  alphas = 0;
  pcs = 0;

  this->uc = 0.0;
  this->vc = 0.0;
  this->fu = 1.0;
  this->fv = 1.0;
  
  // Default: hybrid mode OFF (pure SQPnP only)
  this->use_hybrid_mode = false;
}

opengv::absolute_pose::modules::Sqpnp::~Sqpnp()
{
  delete [] pws;
  delete [] us;
  delete [] alphas;
  delete [] pcs;
}

void
opengv::absolute_pose::modules::Sqpnp::
    set_maximum_number_of_correspondences(int n)
{
  if (maximum_number_of_correspondences < n)
  {
    if (pws != 0) delete [] pws;
    if (us != 0) delete [] us;
    if (alphas != 0) delete [] alphas;
    if (pcs != 0) delete [] pcs;

    maximum_number_of_correspondences = n;
    pws = new double[3 * maximum_number_of_correspondences];
    us = new double[3 * maximum_number_of_correspondences];
    alphas = new double[4 * maximum_number_of_correspondences];
    pcs = new double[3 * maximum_number_of_correspondences];
  }
}

void
opengv::absolute_pose::modules::Sqpnp::reset_correspondences(void)
{
  number_of_correspondences = 0;
}

void
opengv::absolute_pose::modules::Sqpnp::add_correspondence(
    double X,
    double Y,
    double Z,
    double x,
    double y,
    double z)
{
  // Store full normalized bearing vector [x, y, z] instead of [x/z, y/z]
  // Normalize the bearing vector to ensure it's a unit vector
  double norm = sqrt(x*x + y*y + z*z);
  if(norm < EPSILON_ZERO_NORM)
  {
    // Degenerate case: zero bearing vector, skip this correspondence entirely
    return;
  }

  // Only store data after validation passes
  pws[3 * number_of_correspondences    ] = X;
  pws[3 * number_of_correspondences + 1] = Y;
  pws[3 * number_of_correspondences + 2] = Z;

  us[3 * number_of_correspondences    ] = x / norm;
  us[3 * number_of_correspondences + 1] = y / norm;
  us[3 * number_of_correspondences + 2] = z / norm;

  number_of_correspondences++;
}

void
opengv::absolute_pose::modules::Sqpnp::choose_control_points(void)
{
  // Take C0 as the reference points centroid:
  cws[0][0] = cws[0][1] = cws[0][2] = 0;
  for(int i = 0; i < number_of_correspondences; i++)
    for(int j = 0; j < 3; j++)
      cws[0][j] += pws[3 * i + j];

  for(int j = 0; j < 3; j++)
    cws[0][j] /= number_of_correspondences;


  // Take C1, C2, and C3 from PCA on the reference points:
  Eigen::MatrixXd PW0(number_of_correspondences,3);

  for(int i = 0; i < number_of_correspondences; i++)
    for(int j = 0; j < 3; j++)
      PW0(i,j) = pws[3 * i + j] - cws[0][j];

  Eigen::MatrixXd PW0tPW0 = PW0.transpose() * PW0;
  Eigen::JacobiSVD< Eigen::MatrixXd > SVD(
      PW0tPW0,
      Eigen::ComputeFullV | Eigen::ComputeFullU );
  Eigen::MatrixXd D = SVD.singularValues();
  Eigen::MatrixXd Ut = SVD.matrixU().transpose();

  for(int i = 1; i < 4; i++)
  {
    double k = sqrt(D(i - 1,0) / number_of_correspondences);
    for(int j = 0; j < 3; j++)
      cws[i][j] = cws[0][j] + k * Ut((i - 1),j);
  }
}

void
opengv::absolute_pose::modules::Sqpnp::
    compute_barycentric_coordinates(void)
{
  Eigen::Matrix3d CC;

  for(int i = 0; i < 3; i++)
    for(int j = 1; j < 4; j++)
      CC(i,j-1) = cws[j][i] - cws[0][i];

  Eigen::Matrix3d CC_inv = CC.inverse();

  for(int i = 0; i < number_of_correspondences; i++)
  {
    double * pi = pws + 3 * i;
    double * a = alphas + 4 * i;

    for(int j = 0; j < 3; j++)
      a[1 + j] =
        CC_inv(j,0) * (pi[0] - cws[0][0]) +
        CC_inv(j,1) * (pi[1] - cws[0][1]) +
        CC_inv(j,2) * (pi[2] - cws[0][2]);
    a[0] = 1.0f - a[1] - a[2] - a[3];
  }
}

void
opengv::absolute_pose::modules::Sqpnp::fill_M(
    Eigen::MatrixXd & M,
    const int row,
    const double * as,
    const double x,
    const double y,
    const double z)
{
  // For normalized bearing vectors, we use the constraint:
  // The bearing vector [x, y, z] should be parallel to the point in camera frame
  // This gives: [x, y, z] × (R*P + t) = 0
  // Expressing P in barycentric coordinates: P = Σ α_i * C_i
  // And P_c = Σ α_i * c_i (control points in camera frame)
  // So: [x, y, z] × (Σ α_i * c_i) = 0
  // This gives 2 independent constraints per correspondence
  
  // Cross product constraint: [x, y, z] × [cx, cy, cz] = [y*cz - z*cy, z*cx - x*cz, x*cy - y*cx]
  // We use the first two components (they are independent)
  
  for(int i = 0; i < 4; i++)
  {
    // First constraint: y*cz - z*cy = 0
    // For control point c_i = [cx_i, cy_i, cz_i], this becomes:
    // y * cz_i - z * cy_i = 0
    // Expressing in barycentric coordinates: Σ α_j * (y * cz_j - z * cy_j) = 0
    M(row, 3*i)   = 0.0;
    M(row, 3*i+1) = -as[i] * z;  // -z * cy component
    M(row, 3*i+2) = as[i] * y;   // y * cz component
    
    // Second constraint: z*cx - x*cz = 0
    // z * cx_i - x * cz_i = 0
    M(row+1, 3*i)   = as[i] * z;   // z * cx component
    M(row+1, 3*i+1) = 0.0;
    M(row+1, 3*i+2) = -as[i] * x;  // -x * cz component
  }
}

void
opengv::absolute_pose::modules::Sqpnp::compute_ccs(
    const double * betas,
    const Eigen::MatrixXd & ut)
{
  for(int i = 0; i < 4; i++)
    ccs[i][0] = ccs[i][1] = ccs[i][2] = 0.0f;

  for(int i = 0; i < 4; i++)
  {
    for(int j = 0; j < 4; j++)
      for(int k = 0; k < 3; k++)
        ccs[j][k] += betas[i] * ut(11-i,3 * j + k);
  }
}

void
opengv::absolute_pose::modules::Sqpnp::compute_pcs(void)
{
  for(int i = 0; i < number_of_correspondences; i++)
  {
    double * a = alphas + 4 * i;
    double * pc = pcs + 3 * i;

    for(int j = 0; j < 3; j++)
      pc[j] =
          a[0] * ccs[0][j] +
          a[1] * ccs[1][j] +
          a[2] * ccs[2][j] +
          a[3] * ccs[3][j];
  }
}

double
opengv::absolute_pose::modules::Sqpnp::compute_pose(
    double R[3][3],
    double t[3])
{
  // =================================================================
  // TRUE SQPnP ALGORITHM
  // =================================================================
  // Uses Omega matrix formulation that directly supports 360° bearing vectors
  // Based on Terzakis & Lourakis ECCV 2020 paper
  // This is the primary solver for omnidirectional/panorama cameras
  
  Eigen::MatrixXd Omega;
  compute_omega_matrix(Omega);
  
  Eigen::Matrix3d R_sqp;
  Eigen::Vector3d t_sqp;
  double sqp_error = sqp_solve(Omega, R_sqp, t_sqp);
  
  // =================================================================
  // EPnP-BASED FALLBACK (only if hybrid mode is enabled)
  // =================================================================
  if (!use_hybrid_mode)
  {
    // Pure SQPnP mode: use SQPnP result directly
    for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
        R[i][j] = R_sqp(i, j);
      t[i] = t_sqp(i);
    }
    return sqp_error;
  }
  
  // Hybrid mode: Also run the EPnP-based approach as a fallback
  
  choose_control_points();
  compute_barycentric_coordinates();

  Eigen::MatrixXd M(2*number_of_correspondences,12);

  for(int i = 0; i < number_of_correspondences; i++)
    fill_M(M, 2 * i, alphas + 4 * i, us[3 * i], us[3 * i + 1], us[3 * i + 2]);

  Eigen::MatrixXd MtM = M.transpose() * M;
  Eigen::JacobiSVD< Eigen::MatrixXd > SVD(
      MtM,
      Eigen::ComputeFullV | Eigen::ComputeFullU );
  Eigen::MatrixXd Ut = SVD.matrixU().transpose();
  
  // Apply Gram-Schmidt orthogonalization to null space vectors
  Eigen::MatrixXd Ut_orthogonal = Ut;
  gram_schmidt_orthogonalize(Ut_orthogonal, 6, 6);
  Ut = Ut_orthogonal;

  Eigen::Matrix<double,6,10> L_6x10;
  Eigen::Matrix<double,6,1> Rho;

  compute_L_6x10(Ut,L_6x10);
  compute_rho(Rho);

  double Betas[4][4], rep_errors[4];
  double Rs[4][3][3], ts[4][3];

  find_betas_approx_1(L_6x10, Rho, Betas[1]);
  gauss_newton(L_6x10, Rho, Betas[1]);
  rep_errors[1] = compute_R_and_t(Ut, Betas[1], Rs[1], ts[1]);

  find_betas_approx_2(L_6x10, Rho, Betas[2]);
  gauss_newton(L_6x10, Rho, Betas[2]);
  rep_errors[2] = compute_R_and_t(Ut, Betas[2], Rs[2], ts[2]);

  find_betas_approx_3(L_6x10, Rho, Betas[3]);
  gauss_newton(L_6x10, Rho, Betas[3]);
  rep_errors[3] = compute_R_and_t(Ut, Betas[3], Rs[3], ts[3]);
  
  // Additional initialization using average of approximations
  for(int i = 0; i < 4; i++)
  {
    Betas[0][i] = (Betas[1][i] + Betas[2][i] + Betas[3][i]) / 3.0;
  }
  gauss_newton(L_6x10, Rho, Betas[0]);
  rep_errors[0] = compute_R_and_t(Ut, Betas[0], Rs[0], ts[0]);

  // =================================================================
  // SELECT BEST SOLUTION
  // =================================================================
  // Compare true SQPnP result with EPnP-based results
  
  // Find best EPnP-based solution
  int N_epnp = 0;
  for(int i = 1; i < 4; i++)
  {
    if (rep_errors[i] < rep_errors[N_epnp]) N_epnp = i;
  }
  double epnp_best_error = rep_errors[N_epnp];
  
  // Choose between SQPnP and EPnP-based solution
  if(sqp_error < epnp_best_error)
  {
    // Use true SQPnP solution
    for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
        R[i][j] = R_sqp(i, j);
      t[i] = t_sqp(i);
    }
    return sqp_error;
  }
  else
  {
    // Use EPnP-based solution
    copy_R_and_t(Rs[N_epnp], ts[N_epnp], R, t);
    return epnp_best_error;
  }
}

void
opengv::absolute_pose::modules::Sqpnp::copy_R_and_t(
    const double R_src[3][3],
    const double t_src[3],
    double R_dst[3][3],
    double t_dst[3])
{
  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 3; j++)
      R_dst[i][j] = R_src[i][j];
    t_dst[i] = t_src[i];
  }
}

double
opengv::absolute_pose::modules::Sqpnp::dist2(
    const double * p1,
    const double * p2)
{
  return
    (p1[0] - p2[0]) * (p1[0] - p2[0]) +
    (p1[1] - p2[1]) * (p1[1] - p2[1]) +
    (p1[2] - p2[2]) * (p1[2] - p2[2]);
}

double
opengv::absolute_pose::modules::Sqpnp::dot(
    const double * v1,
    const double * v2)
{
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

double
opengv::absolute_pose::modules::Sqpnp::angular_error(
    const double R[3][3],
    const double t[3])
{
  double sum_angular_error = 0.0;

  for(int i = 0; i < number_of_correspondences; i++)
  {
    // Get world point
    double * pw = pws + 3 * i;
    
    // Transform world point to camera frame
    double Xc = dot(R[0], pw) + t[0];
    double Yc = dot(R[1], pw) + t[1];
    double Zc = dot(R[2], pw) + t[2];
    
    // Compute norm squared to check for degenerate case
    double norm_sq = Xc * Xc + Yc * Yc + Zc * Zc;
    if (norm_sq < EPSILON_ZERO_NORM * EPSILON_ZERO_NORM)
    {
      // Degenerate case: point at camera center, skip this correspondence
      continue;
    }
    
    // Normalize to get estimated bearing vector
    double norm = sqrt(norm_sq);
    double b_est_x = Xc / norm;
    double b_est_y = Yc / norm;
    double b_est_z = Zc / norm;
    
    // Get observed bearing vector (already normalized and stored)
    double b_obs_x = us[3 * i];
    double b_obs_y = us[3 * i + 1];
    double b_obs_z = us[3 * i + 2];
    
    // Compute dot product (cosine of angle)
    double cos_angle = b_est_x * b_obs_x + b_est_y * b_obs_y + b_est_z * b_obs_z;
    
    // Clamp to [-1, 1] to avoid numerical issues with acos
    if (cos_angle > 1.0) cos_angle = 1.0;
    if (cos_angle < -1.0) cos_angle = -1.0;
    
    // Compute angular error
    double angle = acos(cos_angle);
    sum_angular_error += angle;
  }

  return sum_angular_error / number_of_correspondences;
}

void
opengv::absolute_pose::modules::Sqpnp::gram_schmidt_orthogonalize(
    Eigen::MatrixXd & V,
    int start_col,
    int num_cols)
{
  // Gram-Schmidt orthogonalization for columns [start_col, start_col+num_cols)
  // This ensures the null space vectors form an orthogonal basis, improving numerical stability
  // and helping achieve global solutions by maintaining proper subspace structure
  
  for(int j = start_col; j < start_col + num_cols && j < V.cols(); j++)
  {
    // Normalize the current column
    double norm = V.col(j).norm();
    if(norm > EPSILON_ZERO_NORM)
    {
      V.col(j) /= norm;
    }
    
    // Subtract projections onto previous columns (Gram-Schmidt process)
    for(int k = start_col; k < j; k++)
    {
      double dot_product = V.col(j).dot(V.col(k));
      V.col(j) -= dot_product * V.col(k);
    }
    
    // Renormalize after orthogonalization
    norm = V.col(j).norm();
    if(norm > EPSILON_ZERO_NORM)
    {
      V.col(j) /= norm;
    }
    else
    {
      // If column becomes zero (degenerate case), set to unit vector in a standard direction
      V.col(j).setZero();
      if(j < V.rows())
        V(j, j) = 1.0;
    }
  }
}

void
opengv::absolute_pose::modules::Sqpnp::project_to_nullspace(
    const Eigen::MatrixXd & M,
    Eigen::MatrixXd & Ut_orthogonal)
{
  // Project solution onto null space for improved numerical stability
  // This helps ensure solutions stay in the correct subspace and find global solutions
  
  Eigen::MatrixXd MtM = M.transpose() * M;
  Eigen::JacobiSVD< Eigen::MatrixXd > SVD(
      MtM,
      Eigen::ComputeFullV | Eigen::ComputeFullU );
  
  Ut_orthogonal = SVD.matrixU().transpose();
  
  // Apply Gram-Schmidt to ensure orthogonality of null space vectors
  gram_schmidt_orthogonalize(Ut_orthogonal, 6, 6);
}

void
opengv::absolute_pose::modules::Sqpnp::compute_omega_matrix(Eigen::MatrixXd & Omega)
{
  // Build the Omega matrix for true SQPnP formulation
  // This directly uses bearing vectors without perspective camera assumption
  // 
  // For each correspondence i with bearing vector u_i and world point M_i:
  // The constraint is: u_i × (R*M_i + t) = 0 (bearing parallel to transformed point)
  // This can be written as: [u_i]× * (R*M_i + t) = 0
  // where [u_i]× is the skew-symmetric matrix of u_i
  //
  // Expanding: [u_i]× * R * M_i + [u_i]× * t = 0
  // This gives us 3 linear constraints per correspondence (2 are independent)
  
  // The Omega matrix is 9x9, formed from sum of outer products
  Omega = Eigen::MatrixXd::Zero(9, 9);
  
  // Q matrix accumulator for translation solving
  Eigen::Matrix3d Q_sum = Eigen::Matrix3d::Zero();
  
  for(int i = 0; i < number_of_correspondences; i++)
  {
    // Get bearing vector (already normalized)
    double ux = us[3 * i];
    double uy = us[3 * i + 1];
    double uz = us[3 * i + 2];
    
    // Get world point
    double px = pws[3 * i];
    double py = pws[3 * i + 1];
    double pz = pws[3 * i + 2];
    
    // Projection matrix: I - u*u^T (projects onto plane perpendicular to u)
    Eigen::Matrix3d P;
    P(0,0) = 1.0 - ux*ux; P(0,1) = -ux*uy;      P(0,2) = -ux*uz;
    P(1,0) = -uy*ux;      P(1,1) = 1.0 - uy*uy; P(1,2) = -uy*uz;
    P(2,0) = -uz*ux;      P(2,1) = -uz*uy;      P(2,2) = 1.0 - uz*uz;
    
    // World point as vector
    Eigen::Vector3d M_i(px, py, pz);
    
    // Accumulate Q_sum for translation estimation
    Q_sum += P;
    
    // Build the 9x9 contribution to Omega from this correspondence
    // For rotation vector r = [r1 r2 r3 r4 r5 r6 r7 r8 r9]^T (column-major of R)
    // The constraint P * R * M_i = 0 can be written as A_i * r = 0
    // where A_i is a 3x9 matrix
    
    Eigen::Matrix<double, 3, 9> A_i = Eigen::Matrix<double, 3, 9>::Zero();
    
    // A_i = P * [M_i^T ⊗ I_3] where ⊗ is Kronecker product
    // This expands to:
    // Column 0-2: P * M_i(0) * I_3 = P * px
    // Column 3-5: P * M_i(1) * I_3 = P * py
    // Column 6-8: P * M_i(2) * I_3 = P * pz
    
    A_i.block<3,3>(0,0) = P * px;
    A_i.block<3,3>(0,3) = P * py;
    A_i.block<3,3>(0,6) = P * pz;
    
    // Accumulate Omega = sum of A_i^T * A_i
    Omega += A_i.transpose() * A_i;
  }
}

void
opengv::absolute_pose::modules::Sqpnp::nearest_rotation_matrix(
    const Eigen::Matrix3d & M_in,
    Eigen::Matrix3d & R_out)
{
  // Find the nearest rotation matrix to M_in using SVD
  // R = U * V^T where M = U * S * V^T
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(M_in, Eigen::ComputeFullU | Eigen::ComputeFullV);
  
  Eigen::Matrix3d U = svd.matrixU();
  Eigen::Matrix3d V = svd.matrixV();
  
  R_out = U * V.transpose();
  
  // Ensure proper rotation (det = 1, not -1)
  if(R_out.determinant() < 0)
  {
    // Flip sign of last column of U
    U.col(2) = -U.col(2);
    R_out = U * V.transpose();
  }
}

double
opengv::absolute_pose::modules::Sqpnp::sqp_solve(
    const Eigen::MatrixXd & Omega,
    Eigen::Matrix3d & R_out,
    Eigen::Vector3d & t_out)
{
  // Sequential Quadratic Programming solver for rotation
  // Minimizes r^T * Omega * r subject to R being a rotation matrix
  //
  // Uses iterative projection: 
  // 1. Solve for r that minimizes quadratic cost
  // 2. Project r onto SO(3) (nearest rotation matrix)
  // 3. Repeat until convergence
  
  const int max_sqp_iterations = 15;
  const double convergence_threshold = 1e-10;
  
  // Initialize with identity rotation
  Eigen::Matrix3d R_current = Eigen::Matrix3d::Identity();
  double prev_cost = std::numeric_limits<double>::max();
  
  for(int iter = 0; iter < max_sqp_iterations; iter++)
  {
    // Vectorize current rotation (column-major)
    Eigen::Map<Eigen::Matrix<double, 9, 1>> r_vec(R_current.data());
    
    // Compute current cost
    double cost = r_vec.transpose() * Omega * r_vec;
    
    // Check convergence
    if(std::abs(prev_cost - cost) < convergence_threshold)
    {
      break;
    }
    prev_cost = cost;
    
    // Compute gradient: g = 2 * Omega * r
    Eigen::Matrix<double, 9, 1> gradient = 2.0 * Omega * r_vec;
    
    // Compute step using Newton-like update with regularization
    // H = 2 * Omega + lambda * I
    double lambda = 1e-6;
    Eigen::MatrixXd H = 2.0 * Omega + lambda * Eigen::MatrixXd::Identity(9, 9);
    
    // Solve H * delta_r = -gradient
    Eigen::Matrix<double, 9, 1> delta_r = H.ldlt().solve(-gradient);
    
    // Update rotation vector
    Eigen::Matrix<double, 9, 1> r_new = r_vec + delta_r;
    
    // Reshape to 3x3 matrix
    Eigen::Matrix3d R_new = Eigen::Map<Eigen::Matrix3d>(r_new.data());
    
    // Project onto SO(3) - find nearest rotation matrix
    nearest_rotation_matrix(R_new, R_current);
  }
  
  R_out = R_current;
  
  // Compute translation given rotation
  // t = -Q_sum^{-1} * sum_i(P_i * R * M_i)
  // where Q_sum = sum_i(P_i) and P_i = I - u_i * u_i^T
  
  Eigen::Matrix3d Q_sum = Eigen::Matrix3d::Zero();
  Eigen::Vector3d q_sum = Eigen::Vector3d::Zero();
  
  for(int i = 0; i < number_of_correspondences; i++)
  {
    // Get bearing vector
    double ux = us[3 * i];
    double uy = us[3 * i + 1];
    double uz = us[3 * i + 2];
    
    // Get world point
    double px = pws[3 * i];
    double py = pws[3 * i + 1];
    double pz = pws[3 * i + 2];
    
    // Projection matrix
    Eigen::Matrix3d P;
    P(0,0) = 1.0 - ux*ux; P(0,1) = -ux*uy;      P(0,2) = -ux*uz;
    P(1,0) = -uy*ux;      P(1,1) = 1.0 - uy*uy; P(1,2) = -uy*uz;
    P(2,0) = -uz*ux;      P(2,1) = -uz*uy;      P(2,2) = 1.0 - uz*uz;
    
    Eigen::Vector3d M_i(px, py, pz);
    
    Q_sum += P;
    q_sum += P * R_out * M_i;
  }
  
  // Solve for translation: Q_sum * t = -q_sum
  t_out = -Q_sum.ldlt().solve(q_sum);
  
  // Compute final angular error
  double total_error = 0.0;
  for(int i = 0; i < number_of_correspondences; i++)
  {
    Eigen::Vector3d u_i(us[3*i], us[3*i+1], us[3*i+2]);
    Eigen::Vector3d M_i(pws[3*i], pws[3*i+1], pws[3*i+2]);

    Eigen::Vector3d p_cam = R_out * M_i + t_out;
    double norm = p_cam.norm();
    if(norm > EPSILON_ZERO_NORM)
    {
      Eigen::Vector3d p_dir = p_cam / norm;
      double cos_angle = u_i.dot(p_dir);
      if(cos_angle > 1.0) cos_angle = 1.0;
      if(cos_angle < -1.0) cos_angle = -1.0;
      total_error += acos(cos_angle);
    }
  }

  return total_error / number_of_correspondences;
}

void
opengv::absolute_pose::modules::Sqpnp::estimate_R_and_t(
    double R[3][3],
    double t[3])
{
  double pc0[3], pw0[3];

  pc0[0] = pc0[1] = pc0[2] = 0.0;
  pw0[0] = pw0[1] = pw0[2] = 0.0;

  for(int i = 0; i < number_of_correspondences; i++)
  {
    const double * pc = pcs + 3 * i;
    const double * pw = pws + 3 * i;

    for(int j = 0; j < 3; j++)
    {
      pc0[j] += pc[j];
      pw0[j] += pw[j];
    }
  }
  for(int j = 0; j < 3; j++)
  {
    pc0[j] /= number_of_correspondences;
    pw0[j] /= number_of_correspondences;
  }

  Eigen::MatrixXd Abt(3,3);

  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 3; j++)
      Abt(i,j) = 0.0;
  }

  for(int i = 0; i < number_of_correspondences; i++)
  {
    double * pc = pcs + 3 * i;
    double * pw = pws + 3 * i;

    for(int j = 0; j < 3; j++)
    {
      Abt(j,0) += (pc[j] - pc0[j]) * (pw[0] - pw0[0]);
      Abt(j,1) += (pc[j] - pc0[j]) * (pw[1] - pw0[1]);
      Abt(j,2) += (pc[j] - pc0[j]) * (pw[2] - pw0[2]);
    }
  }

  Eigen::JacobiSVD< Eigen::MatrixXd > SVD(
      Abt,
      Eigen::ComputeFullV | Eigen::ComputeFullU );
  Eigen::MatrixXd Abt_u = SVD.matrixU();
  Eigen::MatrixXd Abt_v = SVD.matrixV();

  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      R[i][j] = Abt_u.row(i) * Abt_v.row(j).transpose();

  const double det =
      R[0][0] * R[1][1] * R[2][2] +
      R[0][1] * R[1][2] * R[2][0] +
      R[0][2] * R[1][0] * R[2][1] -
      R[0][2] * R[1][1] * R[2][0] -
      R[0][1] * R[1][0] * R[2][2] -
      R[0][0] * R[1][2] * R[2][1];

  //change 1: negative determinant problem is solved by changing Abt_v, not R

  if (det < 0)
  {
    //R[2][0] = -R[2][0];
    //R[2][1] = -R[2][1];
    //R[2][2] = -R[2][2];
    Eigen::MatrixXd Abt_v_prime = Abt_v;
    Abt_v_prime.col(2) = -Abt_v.col(2);
    for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++)
        R[i][j] = Abt_u.row(i) * Abt_v_prime.row(j).transpose();
  }

  t[0] = pc0[0] - dot(R[0], pw0);
  t[1] = pc0[1] - dot(R[1], pw0);
  t[2] = pc0[2] - dot(R[2], pw0);
}

void
opengv::absolute_pose::modules::Sqpnp::print_pose(
    const double R[3][3],
    const double t[3])
{
  cout << R[0][0] << " " << R[0][1] << " " << R[0][2] << " " << t[0] << endl;
  cout << R[1][0] << " " << R[1][1] << " " << R[1][2] << " " << t[1] << endl;
  cout << R[2][0] << " " << R[2][1] << " " << R[2][2] << " " << t[2] << endl;
}

double
opengv::absolute_pose::modules::Sqpnp::compute_R_and_t(
    const Eigen::MatrixXd & Ut,
    const double * betas,
    double R[3][3],
    double t[3])
{
  compute_ccs(betas, Ut);
  compute_pcs();

  estimate_R_and_t(R, t);

  return angular_error(R, t);
}

// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
// betas_approx_1 = [B11 B12     B13         B14]

void
opengv::absolute_pose::modules::Sqpnp::find_betas_approx_1(
    const Eigen::Matrix<double,6,10> & L_6x10,
    const Eigen::Matrix<double,6,1> & Rho,
    double * betas)
{
  Eigen::MatrixXd L_6x4(6,4);

  for(int i = 0; i < 6; i++)
  {
    L_6x4(i,0) = L_6x10(i,0);
    L_6x4(i,1) = L_6x10(i,1);
    L_6x4(i,2) = L_6x10(i,3);
    L_6x4(i,3) = L_6x10(i,6);
  }

  Eigen::JacobiSVD<Eigen::MatrixXd> SVD(
      L_6x4,
      Eigen::ComputeFullV | Eigen::ComputeFullU);
  Eigen::VectorXd Rho_temp = Rho;
  Eigen::VectorXd b4 = SVD.solve(Rho_temp);

  if (b4[0] < 0)
  {
    betas[0] = sqrt(-b4[0]);
    betas[1] = -b4[1] / betas[0];
    betas[2] = -b4[2] / betas[0];
    betas[3] = -b4[3] / betas[0];
  }
  else
  {
    betas[0] = sqrt(b4[0]);
    betas[1] = b4[1] / betas[0];
    betas[2] = b4[2] / betas[0];
    betas[3] = b4[3] / betas[0];
  }
}

// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
// betas_approx_2 = [B11 B12 B22                            ]

void
opengv::absolute_pose::modules::Sqpnp::find_betas_approx_2(
    const Eigen::Matrix<double,6,10> & L_6x10,
    const Eigen::Matrix<double,6,1> & Rho,
    double * betas)
{
  Eigen::MatrixXd L_6x3(6,3);

  for(int i = 0; i < 6; i++)
  {
    L_6x3(i,0) = L_6x10(i,0);
    L_6x3(i,1) = L_6x10(i,1);
    L_6x3(i,2) = L_6x10(i,2);
  }

  Eigen::JacobiSVD<Eigen::MatrixXd> SVD(
      L_6x3,
      Eigen::ComputeFullV | Eigen::ComputeFullU);
  Eigen::VectorXd Rho_temp = Rho;
  Eigen::VectorXd b3 = SVD.solve(Rho_temp);

  if (b3[0] < 0)
  {
    betas[0] = sqrt(-b3[0]);
    betas[1] = (b3[2] < 0) ? sqrt(-b3[2]) : 0.0;
  }
  else
  {
    betas[0] = sqrt(b3[0]);
    betas[1] = (b3[2] > 0) ? sqrt(b3[2]) : 0.0;
  }

  if (b3[1] < 0) betas[0] = -betas[0];

  betas[2] = 0.0;
  betas[3] = 0.0;
}

// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
// betas_approx_3 = [B11 B12 B22 B13 B23                    ]

void
opengv::absolute_pose::modules::Sqpnp::find_betas_approx_3(
    const Eigen::Matrix<double,6,10> & L_6x10,
    const Eigen::Matrix<double,6,1> & Rho,
    double * betas)
{
  Eigen::MatrixXd L_6x5(6,5);

  for(int i = 0; i < 6; i++)
  {
    L_6x5(i,0) = L_6x10(i,0);
    L_6x5(i,1) = L_6x10(i,1);
    L_6x5(i,2) = L_6x10(i,2);
    L_6x5(i,3) = L_6x10(i,3);
    L_6x5(i,4) = L_6x10(i,4);
  }

  Eigen::JacobiSVD<Eigen::MatrixXd> SVD(
      L_6x5,
      Eigen::ComputeFullV | Eigen::ComputeFullU);
  Eigen::VectorXd Rho_temp = Rho;
  Eigen::VectorXd b5 = SVD.solve(Rho_temp);

  if (b5[0] < 0)
  {
    betas[0] = sqrt(-b5[0]);
    betas[1] = (b5[2] < 0) ? sqrt(-b5[2]) : 0.0;
  }
  else
  {
    betas[0] = sqrt(b5[0]);
    betas[1] = (b5[2] > 0) ? sqrt(b5[2]) : 0.0;
  }
  if (b5[1] < 0) betas[0] = -betas[0];
  betas[2] = b5[3] / betas[0];
  betas[3] = 0.0;
}

void
opengv::absolute_pose::modules::Sqpnp::compute_L_6x10(
    const Eigen::MatrixXd & Ut,
    Eigen::Matrix<double,6,10> & L_6x10 )
{
  double dv[4][6][3];

  for(int i = 0; i < 4; i++)
  {
    int a = 0, b = 1;
    for(int j = 0; j < 6; j++)
    {
      dv[i][j][0] = Ut(11-i,3 * a    ) - Ut(11-i,3 * b);
      dv[i][j][1] = Ut(11-i,3 * a + 1) - Ut(11-i,3 * b + 1);
      dv[i][j][2] = Ut(11-i,3 * a + 2) - Ut(11-i,3 * b + 2);

      b++;
      if (b > 3)
      {
        a++;
        b = a + 1;
      }
    }
  }

  for(int i = 0; i < 6; i++)
  {
    L_6x10(i,0) =        dot(dv[0][i], dv[0][i]);
    L_6x10(i,1) = 2.0f * dot(dv[0][i], dv[1][i]);
    L_6x10(i,2) =        dot(dv[1][i], dv[1][i]);
    L_6x10(i,3) = 2.0f * dot(dv[0][i], dv[2][i]);
    L_6x10(i,4) = 2.0f * dot(dv[1][i], dv[2][i]);
    L_6x10(i,5) =        dot(dv[2][i], dv[2][i]);
    L_6x10(i,6) = 2.0f * dot(dv[0][i], dv[3][i]);
    L_6x10(i,7) = 2.0f * dot(dv[1][i], dv[3][i]);
    L_6x10(i,8) = 2.0f * dot(dv[2][i], dv[3][i]);
    L_6x10(i,9) =        dot(dv[3][i], dv[3][i]);
  }
}

void
opengv::absolute_pose::modules::Sqpnp::compute_rho(
    Eigen::Matrix<double,6,1> & Rho)
{
  Rho[0] = dist2(cws[0], cws[1]);
  Rho[1] = dist2(cws[0], cws[2]);
  Rho[2] = dist2(cws[0], cws[3]);
  Rho[3] = dist2(cws[1], cws[2]);
  Rho[4] = dist2(cws[1], cws[3]);
  Rho[5] = dist2(cws[2], cws[3]);
}

void
opengv::absolute_pose::modules::Sqpnp::compute_A_and_b_gauss_newton(
    const Eigen::Matrix<double,6,10> & L_6x10,
    const Eigen::Matrix<double,6,1> & Rho,
    double betas[4],
    Eigen::Matrix<double,6,4> & A,
    Eigen::Matrix<double,6,1> & b)
{
  for(int i = 0; i < 6; i++)
  {
    A(i,0) = 2*L_6x10(i,0)*betas[0] +   L_6x10(i,1)*betas[1] +
               L_6x10(i,3)*betas[2] +   L_6x10(i,6)*betas[3];
    A(i,1) =   L_6x10(i,1)*betas[0] + 2*L_6x10(i,2)*betas[1] +
               L_6x10(i,4)*betas[2] +   L_6x10(i,7)*betas[3];
    A(i,2) =   L_6x10(i,3)*betas[0] +   L_6x10(i,4)*betas[1] +
             2*L_6x10(i,5)*betas[2] +   L_6x10(i,8)*betas[3];
    A(i,3) =   L_6x10(i,6)*betas[0] +   L_6x10(i,7)*betas[1] +
               L_6x10(i,8)*betas[2] + 2*L_6x10(i,9)*betas[3];

    b(i,0) = Rho[i] - (
               L_6x10(i,0) * betas[0] * betas[0] +
               L_6x10(i,1) * betas[0] * betas[1] +
               L_6x10(i,2) * betas[1] * betas[1] +
               L_6x10(i,3) * betas[0] * betas[2] +
               L_6x10(i,4) * betas[1] * betas[2] +
               L_6x10(i,5) * betas[2] * betas[2] +
               L_6x10(i,6) * betas[0] * betas[3] +
               L_6x10(i,7) * betas[1] * betas[3] +
               L_6x10(i,8) * betas[2] * betas[3] +
               L_6x10(i,9) * betas[3] * betas[3]);
  }
}

void
opengv::absolute_pose::modules::Sqpnp::gauss_newton(
    const Eigen::Matrix<double,6,10> & L_6x10,
    const Eigen::Matrix<double,6,1> & Rho,
    double betas[4])
{
  // Increased iterations for better convergence to global solution
  const int iterations_number = 10;  // Increased from 5 to 10 for better convergence

  Eigen::Matrix<double,6,4> A;
  Eigen::Matrix<double,6,1> B;
  Eigen::Matrix<double,4,1> X;

  for(int k = 0; k < iterations_number; k++)
  {
    compute_A_and_b_gauss_newton(L_6x10,Rho,betas,A,B);
    qr_solve(A,B,X);

    for(int i = 0; i < 4; i++)
      betas[i] += X[i];
  }
}

void
opengv::absolute_pose::modules::Sqpnp::qr_solve(
    Eigen::Matrix<double,6,4> & A_orig,
    Eigen::Matrix<double,6,1> & b,
    Eigen::Matrix<double,4,1> & X)
{
  Eigen::Matrix<double,4,6> A = A_orig.transpose();

  static int max_nr = 0;
  static double * A1, * A2;

  const int nr = A_orig.rows();
  const int nc = A_orig.cols();

  if (max_nr != 0 && max_nr < nr)
  {
    delete [] A1;
    delete [] A2;
  }
  if (max_nr < nr)
  {
    max_nr = nr;
    A1 = new double[nr];
    A2 = new double[nr];
  }

  double * pA = A.data(), * ppAkk = pA;
  for(int k = 0; k < nc; k++)
  {
    double * ppAik = ppAkk, eta = fabs(*ppAik);
    for(int i = k + 1; i < nr; i++)
    {
      double elt = fabs(*ppAik);
      if (eta < elt) eta = elt;
      ppAik += nc;
    }

    if (eta == 0)
    {
      A1[k] = A2[k] = 0.0;
      cerr << "Warning: Matrix A is singular. Numerical instability detected." << endl;
      return;
    }
    else
    {
      double * ppAik = ppAkk, sum = 0.0, inv_eta = 1. / eta;
      for(int i = k; i < nr; i++)
      {
        *ppAik *= inv_eta;
        sum += *ppAik * *ppAik;
        ppAik += nc;
      }
      double sigma = sqrt(sum);
      if (*ppAkk < 0)
        sigma = -sigma;
      *ppAkk += sigma;
      A1[k] = sigma * *ppAkk;
      A2[k] = -eta * sigma;
      for(int j = k + 1; j < nc; j++)
      {
        double * ppAik = ppAkk, sum = 0;
        for(int i = k; i < nr; i++)
        {
          sum += *ppAik * ppAik[j - k];
          ppAik += nc;
        }
        double tau = sum / A1[k];
        ppAik = ppAkk;
        for(int i = k; i < nr; i++)
        {
          ppAik[j - k] -= tau * *ppAik;
          ppAik += nc;
        }
      }
    }
    ppAkk += nc + 1;
  }

  // b <- Qt b
  double * ppAjj = pA, * pb = b.data();
  for(int j = 0; j < nc; j++)
  {
    double * ppAij = ppAjj, tau = 0;
    for(int i = j; i < nr; i++)
    {
      tau += *ppAij * pb[i];
      ppAij += nc;
    }
    tau /= A1[j];
    ppAij = ppAjj;
    for(int i = j; i < nr; i++)
    {
      pb[i] -= tau * *ppAij;
      ppAij += nc;
    }
    ppAjj += nc + 1;
  }

  // X = R-1 b
  double * pX = X.data();
  pX[nc - 1] = pb[nc - 1] / A2[nc - 1];
  for(int i = nc - 2; i >= 0; i--)
  {
    double * ppAij = pA + i * nc + (i + 1), sum = 0;

    for(int j = i + 1; j < nc; j++)
    {
      sum += *ppAij * pX[j];
      ppAij++;
    }
    pX[i] = (pb[i] - sum) / A2[i];
  }
}



void
opengv::absolute_pose::modules::Sqpnp::relative_error(
    double & rot_err,
    double & transl_err,
    const double Rtrue[3][3],
    const double ttrue[3],
    const double Rest[3][3],
    const double test[3])
{
  double qtrue[4], qest[4];

  mat_to_quat(Rtrue, qtrue);
  mat_to_quat(Rest, qest);

  double rot_err1 = sqrt((qtrue[0] - qest[0]) * (qtrue[0] - qest[0]) +
                         (qtrue[1] - qest[1]) * (qtrue[1] - qest[1]) +
                         (qtrue[2] - qest[2]) * (qtrue[2] - qest[2]) +
                         (qtrue[3] - qest[3]) * (qtrue[3] - qest[3])) /
                    sqrt(qtrue[0] * qtrue[0] +
                         qtrue[1] * qtrue[1] +
                         qtrue[2] * qtrue[2] +
                         qtrue[3] * qtrue[3]);

  double rot_err2 = sqrt((qtrue[0] + qest[0]) * (qtrue[0] + qest[0]) +
                         (qtrue[1] + qest[1]) * (qtrue[1] + qest[1]) +
                         (qtrue[2] + qest[2]) * (qtrue[2] + qest[2]) +
                         (qtrue[3] + qest[3]) * (qtrue[3] + qest[3]) ) /
                    sqrt(qtrue[0] * qtrue[0] +
                         qtrue[1] * qtrue[1] +
                         qtrue[2] * qtrue[2] +
                         qtrue[3] * qtrue[3]);

  rot_err = min(rot_err1,rot_err2);

  transl_err = sqrt((ttrue[0] - test[0]) * (ttrue[0] - test[0]) +
                    (ttrue[1] - test[1]) * (ttrue[1] - test[1]) +
                    (ttrue[2] - test[2]) * (ttrue[2] - test[2])) /
               sqrt(ttrue[0] * ttrue[0] +
                    ttrue[1] * ttrue[1] +
                    ttrue[2] * ttrue[2]);
}

void
opengv::absolute_pose::modules::Sqpnp::mat_to_quat(
    const double R[3][3],
    double q[4])
{
  double tr = R[0][0] + R[1][1] + R[2][2];
  double n4;

  if (tr > 0.0f)
  {
    q[0] = R[1][2] - R[2][1];
    q[1] = R[2][0] - R[0][2];
    q[2] = R[0][1] - R[1][0];
    q[3] = tr + 1.0f;
    n4 = q[3];
  }
  else if ( (R[0][0] > R[1][1]) && (R[0][0] > R[2][2]) )
  {
    q[0] = 1.0f + R[0][0] - R[1][1] - R[2][2];
    q[1] = R[1][0] + R[0][1];
    q[2] = R[2][0] + R[0][2];
    q[3] = R[1][2] - R[2][1];
    n4 = q[0];
  }
  else if (R[1][1] > R[2][2])
  {
    q[0] = R[1][0] + R[0][1];
    q[1] = 1.0f + R[1][1] - R[0][0] - R[2][2];
    q[2] = R[2][1] + R[1][2];
    q[3] = R[2][0] - R[0][2];
    n4 = q[1];
  }
  else
  {
    q[0] = R[2][0] + R[0][2];
    q[1] = R[2][1] + R[1][2];
    q[2] = 1.0f + R[2][2] - R[0][0] - R[1][1];
    q[3] = R[0][1] - R[1][0];
    n4 = q[2];
  }
  double scale = 0.5f / double(sqrt(n4));

  q[0] *= scale;
  q[1] *= scale;
  q[2] *= scale;
  q[3] *= scale;
}
