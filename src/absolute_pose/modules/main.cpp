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


#include <Eigen/NonLinearOptimization>
#include <Eigen/NumericalDiff>

#include <opengv/absolute_pose/modules/main.hpp>
#include <opengv/absolute_pose/modules/mlpnp.hpp>
#include <opengv/absolute_pose/modules/gp3p/modules.hpp>
#include <opengv/absolute_pose/modules/gpnp1/modules.hpp>
#include <opengv/absolute_pose/modules/gpnp2/modules.hpp>
#include <opengv/absolute_pose/modules/gpnp3/modules.hpp>
#include <opengv/absolute_pose/modules/gpnp4/modules.hpp>
#include <opengv/absolute_pose/modules/gpnp5/modules.hpp>
#include <opengv/absolute_pose/modules/upnp2.hpp>
#include <opengv/absolute_pose/modules/upnp4.hpp>
#include <opengv/OptimizationFunctor.hpp>
#include <opengv/math/roots.hpp>
#include <opengv/math/arun.hpp>
#include <opengv/math/cayley.hpp>


/////////////////////
// MLPnP
/////////////////////
void
opengv::absolute_pose::modules::mlpnp_main(
	const bearingVectors_t& f,
	const points_t& p,
	const cov3_mats_t& covMats,
	const std::vector<int>& indices,
transformation_t& result)
{
	size_t numberCorrespondences = indices.size();
	assert(numberCorrespondences > 5);

	bool planar = false;
	// compute the nullspace of all vectors
	std::vector<Eigen::MatrixXd> nullspaces(numberCorrespondences);
	Eigen::MatrixXd points3(3, numberCorrespondences);
	opengv::points_t points3v(numberCorrespondences);
	opengv::points4_t points4v(numberCorrespondences);
	for (int i = 0; i < numberCorrespondences; i++)
	{
		bearingVector_t f_current = f[indices[i]];
		points3.col(i) = p[indices[i]];
		// nullspace of right vector
		Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::HouseholderQRPreconditioner>
			svd_f(f_current.transpose(), Eigen::ComputeFullV);
		nullspaces[i] = svd_f.matrixV().block(0, 1, 3, 2);
		points3v[i] = p[indices[i]];
	}

	//////////////////////////////////////
	// 1. test if we have a planar scene
	//////////////////////////////////////

	Eigen::Matrix3d planarTest = points3*points3.transpose();
	Eigen::FullPivHouseholderQR<Matrix3d> rankTest(planarTest);
	//int r, c;
	//double minEigenVal = abs(eigen_solver.eigenvalues().real().minCoeff(&r, &c));
	Eigen::Matrix3d eigenRot;
	eigenRot.setIdentity();

	// if yes -> transform points to new eigen frame
	//if (minEigenVal < 1e-3 || minEigenVal == 0.0)
	rankTest.setThreshold(1e-10);
	if (rankTest.rank() == 2)
	{
		planar = true;
		// self adjoint is faster and more accurate than general eigen solvers
		// also has closed form solution for 3x3 self-adjoint matrices
		// in addition this solver sorts the eigenvalues in increasing order
		Eigen::SelfAdjointEigenSolver<Matrix3d> eigen_solver(planarTest);
		eigenRot = eigen_solver.eigenvectors().real();
		eigenRot.transposeInPlace();
		for (int i = 0; i < numberCorrespondences; i++)
			points3.col(i) = eigenRot * points3.col(i);
	}
	//////////////////////////////////////
	// 2. stochastic model
	//////////////////////////////////////
	Eigen::SparseMatrix<double> P(2 * numberCorrespondences,
		2 * numberCorrespondences);
	bool use_cov = false;
	P.setIdentity(); // standard

	// if we do have covariance information 
	// -> fill covariance matrix
	if (covMats.size() == numberCorrespondences)
	{
		use_cov = true;
		int l = 0;
		for (int i = 0; i < numberCorrespondences; ++i)
		{
			// invert matrix
			cov2_mat_t temp = (nullspaces[i].transpose() * covMats[i] * nullspaces[i]).inverse();
			P.coeffRef(l, l) = temp(0, 0);
			P.coeffRef(l, l + 1) = temp(0, 1);
			P.coeffRef(l + 1, l) = temp(1, 0);
			P.coeffRef(l + 1, l + 1) = temp(1, 1);
			l += 2;
		}
	}

	//////////////////////////////////////
	// 3. fill the design matrix A
	//////////////////////////////////////
	const int rowsA = 2 * numberCorrespondences;
	int colsA = 12;
	Eigen::MatrixXd A;
	if (planar)
	{
		colsA = 9;
		A = MatrixXd(rowsA, 9);
	}
	else
		A = MatrixXd(rowsA, 12);
	A.setZero();

	int col = 0;
	// fill design matrix
	if (planar)
	{
		for (int i = 0; i < numberCorrespondences; i++)
		{
			point_t pt3_current = points3.col(i);
			//point4_t points4v = opengv::point4_t(pt3_current[0],
			//	pt3_current[1], pt3_current[2], 1.0);
			//points4v = points4v / points4v.norm();
			// r12
			A(2 * i, 0) = nullspaces[i](0, 0) * pt3_current[1];
			A(2 * i + 1, 0) = nullspaces[i](0, 1) * pt3_current[1];
			// r13
			A(2 * i, 1) = nullspaces[i](0, 0) * pt3_current[2];
			A(2 * i + 1, 1) = nullspaces[i](0, 1) * pt3_current[2];
			// r22
			A(2 * i, 2) = nullspaces[i](1, 0) * pt3_current[1];
			A(2 * i + 1, 2) = nullspaces[i](1, 1)* pt3_current[1];
			// r23
			A(2 * i, 3) = nullspaces[i](1, 0) * pt3_current[2];
			A(2 * i + 1, 3) = nullspaces[i](1, 1) * pt3_current[2];
			// r32
			A(2 * i, 4) = nullspaces[i](2, 0) * pt3_current[1];
			A(2 * i + 1, 4) = nullspaces[i](2, 1) * pt3_current[1];
			// r33
			A(2 * i, 5) = nullspaces[i](2, 0) * pt3_current[2];
			A(2 * i + 1, 5) = nullspaces[i](2, 1) * pt3_current[2];
			// t1
			A(2 * i, 6) = nullspaces[i](0, 0);
			A(2 * i + 1, 6) = nullspaces[i](0, 1);
			// t2
			A(2 * i, 7) = nullspaces[i](1, 0);
			A(2 * i + 1, 7) = nullspaces[i](1, 1);
			// t3
			A(2 * i, 8) = nullspaces[i](2, 0);
			A(2 * i + 1, 8) = nullspaces[i](2, 1);
		}
	}
	else
	{
		for (int i = 0; i < numberCorrespondences; i++)
		{
			point_t pt3_current = points3.col(i);
			//point4_t points4v = opengv::point4_t(pt3_current[0],
			//	pt3_current[1], pt3_current[2], 1.0);
			//points4v = points4v / points4v.norm();
			// r11
			A(2 * i, 0) = nullspaces[i](0, 0) * pt3_current[0];
			A(2 * i + 1, 0) = nullspaces[i](0, 1) * pt3_current[0];
			// r12
			A(2 * i, 1) = nullspaces[i](0, 0) * pt3_current[1];
			A(2 * i + 1, 1) = nullspaces[i](0, 1) * pt3_current[1];
			// r13
			A(2 * i, 2) = nullspaces[i](0, 0) * pt3_current[2];
			A(2 * i + 1, 2) = nullspaces[i](0, 1) * pt3_current[2];
			// r21
			A(2 * i, 3) = nullspaces[i](1, 0) * pt3_current[0];
			A(2 * i + 1, 3) = nullspaces[i](1, 1) * pt3_current[0];
			// r22
			A(2 * i, 4) = nullspaces[i](1, 0) * pt3_current[1];
			A(2 * i + 1, 4) = nullspaces[i](1, 1)* pt3_current[1];
			// r23
			A(2 * i, 5) = nullspaces[i](1, 0) * pt3_current[2];
			A(2 * i + 1, 5) = nullspaces[i](1, 1) * pt3_current[2];
			// r31
			A(2 * i, 6) = nullspaces[i](2, 0) * pt3_current[0];
			A(2 * i + 1, 6) = nullspaces[i](2, 1) * pt3_current[0];
			// r32
			A(2 * i, 7) = nullspaces[i](2, 0) * pt3_current[1];
			A(2 * i + 1, 7) = nullspaces[i](2, 1) * pt3_current[1];
			// r33
			A(2 * i, 8) = nullspaces[i](2, 0) * pt3_current[2];
			A(2 * i + 1, 8) = nullspaces[i](2, 1) * pt3_current[2];
			// t1
			A(2 * i, 9) = nullspaces[i](0, 0);
			A(2 * i + 1, 9) = nullspaces[i](0, 1);
			// t2
			A(2 * i, 10) = nullspaces[i](1, 0);
			A(2 * i + 1, 10) = nullspaces[i](1, 1);
			// t3
			A(2 * i, 11) = nullspaces[i](2, 0);
			A(2 * i + 1, 11) = nullspaces[i](2, 1);
		}
	}

	//////////////////////////////////////
	// 4. solve least squares
	//////////////////////////////////////
	Eigen::MatrixXd AtPA;
	if (use_cov)
		AtPA = A.transpose() * P * A;
	else
		AtPA = A.transpose() * A;

	Eigen::JacobiSVD<Eigen::MatrixXd> svd_A(AtPA, Eigen::ComputeFullV);
	Eigen::MatrixXd result1 = svd_A.matrixV().col(colsA - 1);

	////////////////////////////////
	// now we treat the results differently,
	// depending on the scene (planar or not)
	////////////////////////////////
	//transformation_t T_final;
	rotation_t Rout;
	translation_t tout;
	if (planar) // planar case
	{
		rotation_t tmp;
		// until now, we only estimated 
		// row one and two of the transposed rotation matrix
		tmp << 0.0, result1(0, 0), result1(1, 0),
			0.0, result1(2, 0), result1(3, 0),
			0.0, result1(4, 0), result1(5, 0);
		double scale = 1 / sqrt(tmp.col(1).norm() * tmp.col(2).norm());
		// row 3
		tmp.col(0) = tmp.col(1).cross(tmp.col(2));
		tmp.transposeInPlace();
		// find best rotation matrix in frobenius sense
		Eigen::JacobiSVD<Eigen::MatrixXd> svd_R_frob(tmp, Eigen::ComputeFullU | Eigen::ComputeFullV);
		rotation_t Rout1 = svd_R_frob.matrixU() * svd_R_frob.matrixV().transpose();
		// test if we found a good rotation matrix
		if (Rout1.determinant() < 0)
			Rout1 *= -1.0;
		// rotate this matrix back using the eigen frame
		Rout1 = eigenRot.transpose() * Rout1;

		translation_t t = scale *
			translation_t(result1(6, 0), result1(7, 0), result1(8, 0));
		Rout1.transposeInPlace();
		Rout1 *= -1;
		if (Rout1.determinant() < 0.0)
			Rout1.col(2) *= -1;
		// now we have to find the best out of 4 combinations
		rotation_t R1, R2;
		R1.col(0) = Rout1.col(0); R1.col(1) = Rout1.col(1); R1.col(2) = Rout1.col(2);
		R2.col(0) = -Rout1.col(0); R2.col(1) = -Rout1.col(1); R2.col(2) = Rout1.col(2);

		vector<transformation_t> Ts(4);
		Ts[0].block<3, 3>(0, 0) = R1; Ts[0].block<3, 1>(0, 3) = t;
		Ts[1].block<3, 3>(0, 0) = R1; Ts[1].block<3, 1>(0, 3) = -t;
		Ts[2].block<3, 3>(0, 0) = R2; Ts[2].block<3, 1>(0, 3) = t;
		Ts[3].block<3, 3>(0, 0) = R2; Ts[3].block<3, 1>(0, 3) = -t;

		vector<double> normVal(4);
		for (int i = 0; i < 4; ++i)
		{
			point_t reproPt;
			double norms = 0.0;
			for (int p = 0; p < 6; ++p)
			{
				reproPt = Ts[i].block<3, 3>(0, 0)*points3v[p] + Ts[i].block<3, 1>(0, 3);
				reproPt = reproPt / reproPt.norm();
				norms += 1.0 - reproPt.transpose()*f[indices[p]];
			}
			normVal[i] = norms / 6.0;
		}
		std::vector<double>::iterator
			findMinRepro = std::min_element(std::begin(normVal), std::end(normVal));
		int idx = std::distance(std::begin(normVal), findMinRepro);
		Rout = Ts[idx].block<3, 3>(0, 0);
		tout = Ts[idx].block<3, 1>(0, 3);
		//Rout.transposeInPlace();
		//tout = -Rout*tout;
	}
	else // non-planar
	{
		rotation_t tmp;
		// until now, we only estimated 
		// row one and two of the transposed rotation matrix
		tmp << result1(0, 0), result1(3, 0), result1(6, 0),
			result1(1, 0), result1(4, 0), result1(7, 0),
			result1(2, 0), result1(5, 0), result1(8, 0);
		// get the scale

		double scale = 1.0 / std::pow(std::abs(tmp.col(0).norm() * tmp.col(1).norm()* tmp.col(2).norm()), 1.0 / 3.0);
		// find best rotation matrix in frobenius sense
		Eigen::JacobiSVD<Eigen::MatrixXd> svd_R_frob(tmp, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Rout = svd_R_frob.matrixU() * svd_R_frob.matrixV().transpose();
		// test if we found a good rotation matrix
		if (Rout.determinant() < 0)
			Rout *= -1.0;
		// scale translation
		tout = Rout * (scale *
			translation_t(result1(9, 0), result1(10, 0), result1(11, 0)));
		// find correct direction in terms of reprojection error, just take the first 6 correspondences
		Rout.transposeInPlace();
		double diff1 = 0.0;
		double diff2 = 0.0;
		for (int p = 0; p < 6; ++p)
		{
			bearingVector_t v1a = Rout*(points3v[p] - tout);
			bearingVector_t v1b = Rout*(points3v[p] + tout);
			v1a = v1a / v1a.norm();
			v1b = v1b / v1b.norm();
			diff1 += 1.0 - v1a.transpose()*f[indices[p]];
			diff2 += 1.0 - v1b.transpose()*f[indices[p]];
		}
		diff1 /= 6.0;
		diff2 /= 6.0;
		if (diff1 > diff2)
			tout = -tout;
		tout = -Rout*tout;
	}

	//////////////////////////////////////
	// 5. gauss newton
	//////////////////////////////////////
	cayley_t cay = math::rot2cayley(Rout.transpose());
	Eigen::VectorXd minx(6);
	minx[0] = cay[0];
	minx[1] = cay[1];
	minx[2] = cay[2];
	tout = -Rout.transpose()*tout;
	minx[3] = tout[0];
	minx[4] = tout[1];
	minx[5] = tout[2];

	modules::mlpnp::mlpnp_gn(minx,
		points3v, nullspaces, P, use_cov);

	Rout = math::cayley2rot(cayley_t(minx[0], minx[1], minx[2]));
	tout = translation_t(minx[3], minx[4], minx[5]);
	// result
	result.block<3, 3>(0, 0) = Rout;
	result.block<3, 1>(0, 3) = tout;
}


void
opengv::absolute_pose::modules::p3p_kneip_main(
    const bearingVectors_t & f,
    const points_t & p,
    transformations_t & solutions )
{
  point_t P1 = p[0];
  point_t P2 = p[1];
  point_t P3 = p[2];

  Eigen::Vector3d temp1 = P2 - P1;
  Eigen::Vector3d temp2 = P3 - P1;

  if( temp1.cross(temp2).norm() == 0)
    return;

  bearingVector_t f1 = f[0];
  bearingVector_t f2 = f[1];
  bearingVector_t f3 = f[2];

  Eigen::Vector3d e1 = f1;
  Eigen::Vector3d e3 = f1.cross(f2);
  e3 = e3/e3.norm();
  Eigen::Vector3d e2 = e3.cross(e1);

  rotation_t T;
  T.row(0) = e1.transpose();
  T.row(1) = e2.transpose();
  T.row(2) = e3.transpose();

  f3 = T*f3;

  if( f3(2,0) > 0)
  {
    f1 = f[1];
    f2 = f[0];
    f3 = f[2];

    e1 = f1;
    e3 = f1.cross(f2);
    e3 = e3/e3.norm();
    e2 = e3.cross(e1);

    T.row(0) = e1.transpose();
    T.row(1) = e2.transpose();
    T.row(2) = e3.transpose();

    f3 = T*f3;

    P1 = p[1];
    P2 = p[0];
    P3 = p[2];
  }

  Eigen::Vector3d n1 = P2-P1;
  n1 = n1/n1.norm();
  Eigen::Vector3d n3 = n1.cross(P3-P1);
  n3 = n3/n3.norm();
  Eigen::Vector3d n2 = n3.cross(n1);

  rotation_t N;
  N.row(0) = n1.transpose();
  N.row(1) = n2.transpose();
  N.row(2) = n3.transpose();

  P3 = N*(P3-P1);

  double d_12 = temp1.norm();
  double f_1 = f3(0,0)/f3(2,0);
  double f_2 = f3(1,0)/f3(2,0);
  double p_1 = P3(0,0);
  double p_2 = P3(1,0);

  double cos_beta = f1.dot(f2);
  double b = 1/( 1 - pow( cos_beta, 2 ) ) - 1;

  if( cos_beta < 0 )
    b = -sqrt(b);
  else
    b = sqrt(b);

  double f_1_pw2 = pow(f_1,2);
  double f_2_pw2 = pow(f_2,2);
  double p_1_pw2 = pow(p_1,2);
  double p_1_pw3 = p_1_pw2 * p_1;
  double p_1_pw4 = p_1_pw3 * p_1;
  double p_2_pw2 = pow(p_2,2);
  double p_2_pw3 = p_2_pw2 * p_2;
  double p_2_pw4 = p_2_pw3 * p_2;
  double d_12_pw2 = pow(d_12,2);
  double b_pw2 = pow(b,2);

  Eigen::Matrix<double,5,1> factors;

  factors(0,0) = -f_2_pw2*p_2_pw4
                 -p_2_pw4*f_1_pw2
                 -p_2_pw4;

  factors(1,0) = 2*p_2_pw3*d_12*b
                 +2*f_2_pw2*p_2_pw3*d_12*b
                 -2*f_2*p_2_pw3*f_1*d_12;

  factors(2,0) = -f_2_pw2*p_2_pw2*p_1_pw2
                 -f_2_pw2*p_2_pw2*d_12_pw2*b_pw2
                 -f_2_pw2*p_2_pw2*d_12_pw2
                 +f_2_pw2*p_2_pw4
                 +p_2_pw4*f_1_pw2
                 +2*p_1*p_2_pw2*d_12
                 +2*f_1*f_2*p_1*p_2_pw2*d_12*b
                 -p_2_pw2*p_1_pw2*f_1_pw2
                 +2*p_1*p_2_pw2*f_2_pw2*d_12
                 -p_2_pw2*d_12_pw2*b_pw2
                 -2*p_1_pw2*p_2_pw2;

  factors(3,0) = 2*p_1_pw2*p_2*d_12*b
                 +2*f_2*p_2_pw3*f_1*d_12
                 -2*f_2_pw2*p_2_pw3*d_12*b
                 -2*p_1*p_2*d_12_pw2*b;

  factors(4,0) = -2*f_2*p_2_pw2*f_1*p_1*d_12*b
                 +f_2_pw2*p_2_pw2*d_12_pw2
                 +2*p_1_pw3*d_12
                 -p_1_pw2*d_12_pw2
                 +f_2_pw2*p_2_pw2*p_1_pw2
                 -p_1_pw4
                 -2*f_2_pw2*p_2_pw2*p_1*d_12
                 +p_2_pw2*f_1_pw2*p_1_pw2
                 +f_2_pw2*p_2_pw2*d_12_pw2*b_pw2;

  std::vector<double> realRoots = math::o4_roots(factors);

  for( int i = 0; i < 4; i++ )
  {
    double cot_alpha =
        (-f_1*p_1/f_2-realRoots[i]*p_2+d_12*b)/
        (-f_1*realRoots[i]*p_2/f_2+p_1-d_12);

    double cos_theta = realRoots[i];
    double sin_theta = sqrt(1-pow(realRoots[i],2));
    double sin_alpha = sqrt(1/(pow(cot_alpha,2)+1));
    double cos_alpha = sqrt(1-pow(sin_alpha,2));

    if (cot_alpha < 0)
      cos_alpha = -cos_alpha;

    translation_t C;
    C(0,0) = d_12*cos_alpha*(sin_alpha*b+cos_alpha);
    C(1,0) = cos_theta*d_12*sin_alpha*(sin_alpha*b+cos_alpha);
    C(2,0) = sin_theta*d_12*sin_alpha*(sin_alpha*b+cos_alpha);

    C = P1 + N.transpose()*C;

    rotation_t R;
    R(0,0) = -cos_alpha;
    R(0,1) = -sin_alpha*cos_theta;
    R(0,2) = -sin_alpha*sin_theta;
    R(1,0) = sin_alpha;
    R(1,1) = -cos_alpha*cos_theta;
    R(1,2) = -cos_alpha*sin_theta;
    R(2,0) = 0.0;
    R(2,1) = -sin_theta;
    R(2,2) = cos_theta;

    R = N.transpose()*R.transpose()*T;

    transformation_t solution;
    solution.col(3) = C;
    solution.block<3,3>(0,0) = R;

    solutions.push_back(solution);
  }
}

void
opengv::absolute_pose::modules::p3p_gao_main(
    const bearingVectors_t & f,
    const points_t & points,
    transformations_t & solutions )
{
  point_t A = points[0];
  point_t B = points[1];
  point_t C = points[2];

  Eigen::Vector3d tempp;
  tempp = A-B;
  double AB = tempp.norm();
  tempp = B-C;
  double BC = tempp.norm();
  tempp = A-C;
  double AC = tempp.norm();

  bearingVector_t f1 = f[0];
  bearingVector_t f2 = f[1];
  bearingVector_t f3 = f[2];

  double cosalpha = f2.transpose()*f3;
  double cosbeta = f1.transpose()*f3;
  double cosgamma = f1.transpose()*f2;

  double a=pow((BC/AB),2);
  double b=pow((AC/AB),2);
  double p=2*cosalpha;
  double q=2*cosbeta;
  double r=2*cosgamma;

  double aSq = a * a;
  double bSq = b * b;
  double pSq = p*p;
  double qSq = q*q;
  double rSq = r*r;

  if ((pSq + qSq + rSq - p*q*r - 1) == 0)
    return;

  Eigen::Matrix<double,5,1> factors;

  factors[0] = -2*b + bSq + aSq + 1 - b*rSq*a + 2*b*a - 2*a;

  if (factors[0] == 0)
    return;

  factors[1] =
      -2*b*q*a - 2*aSq*q + b*rSq*q*a - 2*q + 2*b*q +
      4*a*q + p*b*r + b*r*p*a - bSq*r*p;
  factors[2] =
      qSq + bSq*rSq - b*pSq - q*p*b*r + bSq*pSq - b*rSq*a +
      2 - 2*bSq - a*b*r*p*q + 2*aSq - 4*a - 2*qSq*a + qSq*aSq;
  factors[3] =
      -bSq*r*p + b*r*p*a - 2*aSq*q + q*pSq*b +
      2*b*q*a + 4*a*q + p*b*r - 2*b*q - 2*q;
  factors[4] = 1 - 2*a + 2*b + bSq - b*pSq + aSq - 2*b*a;

  std::vector<double> x_temp = math::o4_roots(factors);
  Eigen::Matrix<double,4,1> x;
  for( size_t i = 0; i < 4; i++ ) x[i] = x_temp[i];

  double temp = (pSq*(a-1+b) + p*q*r - q*a*r*p + (a-1-b)*rSq);
  double b0 = b * temp * temp;

  double rCb = rSq*r;

  Eigen::Matrix<double,4,1> tempXP2;
  tempXP2[0] = x[0]*x[0];
  tempXP2[1] = x[1]*x[1];
  tempXP2[2] = x[2]*x[2];
  tempXP2[3] = x[3]*x[3];
  Eigen::Matrix<double,4,1> tempXP3;
  tempXP3[0] = tempXP2[0]*x[0];
  tempXP3[1] = tempXP2[1]*x[1];
  tempXP3[2] = tempXP2[2]*x[2];
  tempXP3[3] = tempXP2[3]*x[3];

  Eigen::Matrix<double,4,1> ones;
  for( size_t i = 0; i < 4; i++) ones[i] = 1.0;

  Eigen::Matrix<double,4,1> b1_part1 =
      (1-a-b)*tempXP2 + (q*a-q)*x + (1 - a + b)*ones;

  Eigen::Matrix<double,4,1> b1_part2 =
      (aSq*rCb + 2*b*rCb*a - b*rSq*rCb*a - 2*a*rCb + rCb + bSq*rCb
      - 2*rCb*b)*tempXP3
      +(p*rSq + p*aSq*rSq - 2*b*rCb*q*a + 2*rCb*b*q - 2*rCb*q - 2*p*(a+b)*rSq
      + rSq*rSq*p*b + 4*a*rCb*q + b*q*a*rCb*rSq - 2*rCb*aSq*q +2*rSq*p*b*a
      + bSq*rSq*p - rSq*rSq*p*bSq)*tempXP2
      +(rCb*qSq + rSq*rCb*bSq + r*pSq*bSq - 4*a*rCb - 2*a*rCb*qSq + rCb*qSq*aSq
      + 2*aSq*rCb - 2*bSq*rCb - 2*pSq*b*r + 4*p*a*rSq*q + 2*a*pSq*r*b
      - 2*a*rSq*q*b*p - 2*pSq*a*r + r*pSq - b*rSq*rCb*a + 2*p*rSq*b*q
      + r*pSq*aSq -2*p*q*rSq + 2*rCb - 2*rSq*p*aSq*q - rSq*rSq*q*b*p)*x
      +(4*a*rCb*q + p*rSq*qSq + 2*pSq*p*b*a - 4*p*a*rSq - 2*rCb*b*q - 2*pSq*q*r
      - 2*bSq*rSq*p + rSq*rSq*p*b + 2*p*aSq*rSq - 2*rCb*aSq*q - 2*pSq*p*a
      + pSq*p*aSq + 2*p*rSq + pSq*p + 2*b*rCb*q*a + 2*q*pSq*b*r + 4*q*a*r*pSq
      - 2*p*a*rSq*qSq - 2*pSq*aSq*r*q + p*aSq*rSq*qSq - 2*rCb*q - 2*pSq*p*b
      + pSq*p*bSq - 2*pSq*b*r*q*a)*ones;

  Eigen::Matrix<double,4,1> b1;
  b1[0] = b1_part1[0]*b1_part2[0];
  b1[1] = b1_part1[1]*b1_part2[1];
  b1[2] = b1_part1[2]*b1_part2[2];
  b1[3] = b1_part1[3]*b1_part2[3];

  Eigen::Matrix<double,4,1> y=b1/b0;
  Eigen::Matrix<double,4,1> tempYP2;
  tempYP2[0] = pow(y[0],2);
  tempYP2[1] = pow(y[1],2);
  tempYP2[2] = pow(y[2],2);
  tempYP2[3] = pow(y[3],2);

  Eigen::Matrix<double,4,1> tempXY;
  tempXY[0] = x[0]*y[0];
  tempXY[1] = x[1]*y[1];
  tempXY[2] = x[2]*y[2];
  tempXY[3] = x[3]*y[3];

  Eigen::Matrix<double,4,1> v= tempXP2 + tempYP2 - r*tempXY;

  Eigen::Matrix<double,4,1> Z;
  Z[0] = AB/sqrt(v[0]);
  Z[1] = AB/sqrt(v[1]);
  Z[2] = AB/sqrt(v[2]);
  Z[3] = AB/sqrt(v[3]);

  Eigen::Matrix<double,4,1> X;
  X[0] = x[0]*Z[0];
  X[1] = x[1]*Z[1];
  X[2] = x[2]*Z[2];
  X[3] = x[3]*Z[3];

  Eigen::Matrix<double,4,1> Y;
  Y[0] = y[0]*Z[0];
  Y[1] = y[1]*Z[1];
  Y[2] = y[2]*Z[2];
  Y[3] = y[3]*Z[3];

  for( int i = 0; i < 4; i++ )
  {
    //apply arun to find the transformation
    points_t p_cam;
    p_cam.push_back(X[i]*f1);
    p_cam.push_back(Y[i]*f2);
    p_cam.push_back(Z[i]*f3);

    transformation_t solution = math::arun_complete(points,p_cam);
    solutions.push_back(solution);
  }
}

void
opengv::absolute_pose::modules::gp3p_main(
    const Eigen::Matrix3d & f,
    const Eigen::Matrix3d & v,
    const Eigen::Matrix3d & p,
    transformations_t & solutions)
{
  Eigen::Matrix<double,48,85> groebnerMatrix =
      Eigen::Matrix<double,48,85>::Zero();
  gp3p::init(groebnerMatrix,f,v,p);
  gp3p::compute(groebnerMatrix);

  Eigen::Matrix<double,8,8> M = Eigen::Matrix<double,8,8>::Zero();
  M.block<6,8>(0,0) = -groebnerMatrix.block<6,8>(36,77);
  M(6,0) = 1.0;
  M(7,6) = 1.0;

  Eigen::EigenSolver< Eigen::Matrix<double,8,8> > Eig(M,true);
  Eigen::Matrix<std::complex<double>,8,1> D = Eig.eigenvalues();
  Eigen::Matrix<std::complex<double>,8,8> V = Eig.eigenvectors();

  for( int c = 0; c < V.cols(); c++ )
  {
    std::complex<double> eigValue = D[c];

    if( eigValue.imag() < 0.0001 )
    {
      cayley_t cayley;
      Eigen::Vector3d n;

      for(size_t i = 0; i < 3; i++)
      {
        std::complex<double> cay = V(i+4,c)/V(7,c);
        cayley[2-i] = cay.real();
        std::complex<double> depth = V(i+1,c)/V(7,c);
        n[2-i] = depth.real();
      }

      rotation_t rotation = math::cayley2rot(cayley);
      //the groebner problem was set up to find the transpose!
      rotation.transposeInPlace();

      point_t center_cam = Eigen::Vector3d::Zero();
      point_t center_world = Eigen::Vector3d::Zero();
      for( size_t i = 0; i < (size_t) f.cols(); i++ )
      {
        point_t temp = rotation*(n[i]*f.col(i)+v.col(i));
        center_cam = center_cam + temp;
        center_world = center_world + p.col(i);
      }

      center_cam = center_cam/f.cols();
      center_world = center_world/f.cols();
      translation_t translation = center_world - center_cam;

      transformation_t transformation;
      transformation.block<3,3>(0,0) = rotation;
      transformation.col(3) = translation;
      solutions.push_back(transformation);
    }
  }
}

void
opengv::absolute_pose::modules::gpnp_main(
    const Eigen::Matrix<double,12,1> & a,
    const Eigen::Matrix<double,12,12> & V,
    const points_t & c,
    transformation_t & transformation )
{
  //extracting the nullspace vectors
  Eigen::Matrix<double,12,1> vec_5 = V.col(7);
  Eigen::Matrix<double,12,1> vec_4 = V.col(8);
  Eigen::Matrix<double,12,1> vec_3 = V.col(9);
  Eigen::Matrix<double,12,1> vec_2 = V.col(10);
  Eigen::Matrix<double,12,1> vec_1 = V.col(11);

  point_t c0 = c[0];
  point_t c1 = c[0];
  point_t c2 = c[0];
  point_t c3 = c[0];

  Eigen::Matrix<double,12,1> solution;
  std::vector<double> errors;
  translation_t t;
  translations_t ts;
  rotation_t R;
  rotations_t Rs;
  std::vector<double> factors;

  solution = a;
  errors.push_back(gpnp_evaluate(solution,c,t,R));
  ts.push_back(t);
  Rs.push_back(R);

  //nice, now we just need to find the right combination
  //let's start with trying out the linear combination of the most right
  //null-space vector
  Eigen::Matrix<double,5,3> groebnerMatrix1 =
      Eigen::Matrix<double,5,3>::Zero();
  gpnp1::init(groebnerMatrix1,a,vec_1,c0,c1,c2,c3);
  gpnp1::compute(groebnerMatrix1);
  factors.push_back(-groebnerMatrix1(3,2)/groebnerMatrix1(3,1));
  gpnp_optimize( a, V, c, factors );
  solution = a;
  for(size_t i = 0; i < factors.size(); i++)
    solution += factors[i]*V.col(12-factors.size()+i);
  errors.push_back(gpnp_evaluate(solution,c,t,R));
  ts.push_back(t);
  Rs.push_back(R);

  //now let's compute the solution using two nullspace vectors
  Eigen::Matrix<double,10,6> groebnerMatrix2 =
      Eigen::Matrix<double,10,6>::Zero();
  gpnp2::init(groebnerMatrix2,a,vec_2,vec_1,c0,c1,c2,c3);
  gpnp2::compute(groebnerMatrix2);
  factors[0] = -groebnerMatrix2(8,5)/groebnerMatrix2(8,4);
  factors.push_back(
      -(groebnerMatrix2(7,4)*factors[0]+groebnerMatrix2(7,5))/
      groebnerMatrix2(7,3));
  gpnp_optimize( a, V, c, factors );
  solution = a;
  for(size_t i = 0; i < factors.size(); i++)
    solution += factors[i]*V.col(12-factors.size()+i);
  errors.push_back(gpnp_evaluate(solution,c,t,R));
  ts.push_back(t);
  Rs.push_back(R);

  //now let's compute the solution using three nullspace vectors
  Eigen::Matrix<double,15,18> groebnerMatrix3 =
      Eigen::Matrix<double,15,18>::Zero();
  gpnp3::init(groebnerMatrix3,a,vec_3,vec_2,vec_1,c0,c1,c2,c3);
  gpnp3::compute(groebnerMatrix3);
  factors[0] = -groebnerMatrix3(13,17)/groebnerMatrix3(13,16);
  factors[1] =
      -(groebnerMatrix3(12,16)*factors[0]+groebnerMatrix3(12,17))/
      groebnerMatrix3(12,15);
  factors.push_back(
      -(groebnerMatrix3(11,15)*factors[1]+groebnerMatrix3(11,16)*factors[0]+
      groebnerMatrix3(11,17))/groebnerMatrix3(11,14));
  gpnp_optimize( a, V, c, factors );
  solution = a;
  for(size_t i = 0; i < factors.size(); i++)
    solution += factors[i]*V.col(12-factors.size()+i);
  errors.push_back(gpnp_evaluate(solution,c,t,R));
  ts.push_back(t);
  Rs.push_back(R);

  //now let's compute the solution using four nullspace vectors
  Eigen::Matrix<double,25,37> groebnerMatrix4 =
      Eigen::Matrix<double,25,37>::Zero();
  gpnp4::init(groebnerMatrix4,a,vec_4,vec_3,vec_2,vec_1,c0,c1,c2,c3);
  gpnp4::compute(groebnerMatrix4);
  factors[0] = -groebnerMatrix4(23,36)/groebnerMatrix4(23,35);
  factors[1] =
      -(groebnerMatrix4(22,35)*factors[0]+groebnerMatrix4(22,36))/
      groebnerMatrix4(22,34);
  factors[2] =
      -(groebnerMatrix4(21,34)*factors[1]+groebnerMatrix4(21,35)*factors[0]+
      groebnerMatrix4(21,36))/groebnerMatrix4(21,33);
  factors.push_back(
      -(groebnerMatrix4(20,33)*factors[2]+groebnerMatrix4(20,34)*factors[1]+
      groebnerMatrix4(20,35)*factors[0]+groebnerMatrix4(20,36))/
      groebnerMatrix4(20,32));
  gpnp_optimize( a, V, c, factors );
  solution = a;
  for(size_t i = 0; i < factors.size(); i++)
    solution += factors[i]*V.col(12-factors.size()+i);
  errors.push_back(gpnp_evaluate(solution,c,t,R));
  ts.push_back(t);
  Rs.push_back(R);

  //now let's compute the solution using five nullspace vectors
  Eigen::Matrix<double,44,80> groebnerMatrix5 =
      Eigen::Matrix<double,44,80>::Zero();
  gpnp5::init(groebnerMatrix5,a,vec_5,vec_4,vec_3,vec_2,vec_1,c0,c1,c2,c3);
  gpnp5::compute(groebnerMatrix5);
  factors[0] = -groebnerMatrix5(42,79)/groebnerMatrix5(42,78);
  factors[1] =
      -(groebnerMatrix5(41,78)*factors[0]+groebnerMatrix5(41,79))/
      groebnerMatrix5(41,77);
  factors[2] =
      -(groebnerMatrix5(40,77)*factors[1]+groebnerMatrix5(40,78)*factors[0]+
      groebnerMatrix5(40,79))/groebnerMatrix5(40,76);
  factors[3] =
      -(groebnerMatrix5(39,76)*factors[2]+groebnerMatrix5(39,77)*factors[1]+
      groebnerMatrix5(39,78)*factors[0]+groebnerMatrix5(39,79))/
      groebnerMatrix5(39,75);
  factors.push_back(
      -(groebnerMatrix5(38,75)*factors[3]+groebnerMatrix5(38,76)*factors[1]+
      groebnerMatrix5(38,77)*factors[1]+groebnerMatrix5(38,78)*factors[0]+
      groebnerMatrix5(38,79))/groebnerMatrix5(38,74));
  gpnp_optimize( a, V, c, factors );
  solution = a;
  for(size_t i = 0; i < factors.size(); i++)
    solution += factors[i]*V.col(12-factors.size()+i);
  errors.push_back(gpnp_evaluate(solution,c,t,R));
  ts.push_back(t);
  Rs.push_back(R);

  //find best solution
  double smallestError = errors.at(0);
  int minimumIndex = 0;
  for( int i = 1; i < 6; i++ )
  {
    if( errors.at(i) < smallestError )
    {
      smallestError = errors.at(i);
      minimumIndex = i;
    }
  }

  transformation.col(3) = ts.at(minimumIndex);
  transformation.block<3,3>(0,0) = Rs.at(minimumIndex);
}

double
opengv::absolute_pose::modules::gpnp_evaluate(
    const Eigen::Matrix<double,12,1> & solution,
    const points_t & c,
    translation_t & t,
    rotation_t & R )
{
  points_t ccam;
  for(size_t i = 0; i<4; i++)
    ccam.push_back(solution.block<3,1>(i*3,0));

  transformation_t transformation = math::arun_complete(c,ccam);
  t = transformation.col(3);
  R = transformation.block<3,3>(0,0);

  //transform world points into camera frame and compute the error
  double error = 0.0;
  for(size_t i = 0; i<4; i++)
  {
    point_t ccam_reprojected = R.transpose() * (c[i] - t);
    error +=
        1.0 -
        (ccam_reprojected.dot(ccam[i])/(ccam[i].norm()*ccam_reprojected.norm()));
  }

  return error;
}

namespace opengv
{
namespace absolute_pose
{
namespace modules
{

struct GpnpOptimizationFunctor : OptimizationFunctor<double>
{
  const Eigen::Matrix<double,12,1> & _a;
  const Eigen::Matrix<double,12,12> & _V;
  const points_t & _c;
  size_t _dim;

  GpnpOptimizationFunctor(
      const Eigen::Matrix<double,12,1> & a,
      const Eigen::Matrix<double,12,12> & V,
      const points_t & c,
      size_t dim ) :
      OptimizationFunctor<double>(dim,6),
      _a(a),
      _V(V),
      _c(c),
      _dim(dim) {}

  int operator()(const VectorXd &x, VectorXd &fvec) const
  {
    assert( x.size() == _dim );
    assert( (unsigned int) fvec.size() == 6);

    Eigen::Matrix<double,12,1> solution = _a;
    for(size_t i = 0; i < _dim; i++)
      solution += x[i]*_V.col(12-_dim+i);

    points_t ccam;
    for(size_t i = 0; i<4; i++)
      ccam.push_back(solution.block<3,1>(i*3,0));

    Eigen::Vector3d diffw;
    Eigen::Vector3d diffc;
    size_t index = 0;

    for(size_t i = 0; i<3; i++)
    {
      for(size_t j = i+1; j < 4; j++)
      {
        diffw = _c[i]-_c[j];
        diffc = ccam[i]-ccam[j];
        fvec[index++] = diffw.dot(diffw)-diffc.dot(diffc);
      }
    }

    return 0;
  }
};

}
}
}

void
opengv::absolute_pose::modules::gpnp_optimize(
    const Eigen::Matrix<double,12,1> & a,
    const Eigen::Matrix<double,12,12> & V,
    const points_t & c,
    std::vector<double> & factors )
{
  const int n=factors.size();
  VectorXd x(n);

  for(size_t i = 0; i < factors.size(); i++)
    x[i] = factors[i];

  GpnpOptimizationFunctor functor( a, V, c, factors.size() );
  NumericalDiff<GpnpOptimizationFunctor> numDiff(functor);
  LevenbergMarquardt< NumericalDiff<GpnpOptimizationFunctor> > lm(numDiff);

  lm.resetParameters();
  lm.parameters.ftol = 1.E10*NumTraits<double>::epsilon();
  lm.parameters.xtol = 1.E10*NumTraits<double>::epsilon();
  lm.parameters.maxfev = 1000;
  lm.minimize(x);

  for(size_t i = 0; i < factors.size(); i++)
    factors[i] = x[i];
}

void
opengv::absolute_pose::modules::upnp_fill_s(
    const Eigen::Vector4d & quaternion,
    Eigen::Matrix<double,10,1> & s )
{
  s[0] = quaternion[0] * quaternion[0];
  s[1] = quaternion[1] * quaternion[1];
  s[2] = quaternion[2] * quaternion[2];
  s[3] = quaternion[3] * quaternion[3];
  s[4] = quaternion[0] * quaternion[1];
  s[5] = quaternion[0] * quaternion[2];
  s[6] = quaternion[0] * quaternion[3];
  s[7] = quaternion[1] * quaternion[2];
  s[8] = quaternion[1] * quaternion[3];
  s[9] = quaternion[2] * quaternion[3];
}

//we use this one if the number of correspondences is pretty low (more robust)
void
opengv::absolute_pose::modules::upnp_main(
    const Eigen::Matrix<double,10,10> & M,
    const Eigen::Matrix<double,1,10> & C,
    double gamma,
    std::vector<std::pair<double,Eigen::Vector4d>,Eigen::aligned_allocator< std::pair<double,Eigen::Vector4d> > > & quaternions )
{
  Eigen::Matrix<double,16,16> Action;
  upnp::setupAction_gj( M, C, gamma, Action );
  Eigen::EigenSolver< Eigen::Matrix<double,16,16> > Eig( Action, true );
  Eigen::Matrix<std::complex<double>,16,16> V = Eig.eigenvectors();
  
  //cut the double solutions
  double doubleSolThreshold = 0.00000001;
  
  for( int i = 0; i < 16; i++ )
  {
    //we decided to drop the test for imaginary part
    //I've noticed that when the number of points is really low, things get a little
    //weary with noise, and complex solutions might actually be pretty good
    
    Eigen::Vector4d quaternion;
    double norm = 0.0;
    for( int q = 0; q < 4; q++ )
    {
      quaternion[q] = V(11+q,i).real();
      norm += pow(quaternion[q],2.0);
    }
    norm = sqrt(norm);
    if(quaternion[0] < 0) // this here is maybe risky, what if quaternion[0] is very small
      norm *= -1.0;
    for( int q = 0; q < 4; q++ )
      quaternion[q] /= norm;
    
    bool alreadyThere = false;
    for( size_t s = 0; s < quaternions.size(); s++ )
    {
      Eigen::Vector4d diff = quaternion - quaternions[s].second;
      if( diff.norm() < doubleSolThreshold )
      {
        alreadyThere = true;
        break;
      }
    }
    
    if( !alreadyThere )
    {
      Eigen::Matrix<double,10,1> s;
      upnp_fill_s(quaternion,s);
      Eigen::Matrix<double,1,1> valueM = s.transpose() * M * s + 2.0 * C * s;
      double value = valueM[0] + gamma;
      
      std::vector<std::pair<double,Eigen::Vector4d>,Eigen::aligned_allocator< std::pair<double,Eigen::Vector4d> > >::iterator
          qidx = quaternions.begin();
      while( qidx != quaternions.end() && qidx->first < value )
        qidx++;
      
      quaternions.insert(qidx,std::pair<double,Eigen::Vector4d>(value,quaternion));
    }
  }
}

//this one is the really fast, symmetric version, that we use in the normal case
void
opengv::absolute_pose::modules::upnp_main_sym(
    const Eigen::Matrix<double,10,10> & M,
    const Eigen::Matrix<double,1,10> & C,
    double gamma,
    std::vector<std::pair<double,Eigen::Vector4d>,Eigen::aligned_allocator< std::pair<double,Eigen::Vector4d> > > & quaternions )
{
  Eigen::Matrix<double,8,8> Action;
  upnp::setupAction_sym_gj( M, C, gamma, Action );
  Eigen::EigenSolver< Eigen::Matrix<double,8,8> > Eig( Action, true );
  Eigen::Matrix<std::complex<double>,8,8> V = Eig.eigenvectors();
  
  //ok, let's cut the imaginary solutions (with a reasonable threshold!)
  // const double imagThreshold = 0.01;
  std::vector<std::pair<double,Eigen::Vector4d>,Eigen::aligned_allocator< std::pair<double,Eigen::Vector4d> > > bad_quaternions;
  
  for( int i = 0; i < 8; i++ )
  {
    Eigen::Vector4d quaternion;
    quaternion[3] = V(7,i).real();
    quaternion[2] = V(6,i).real();
    quaternion[1] = V(5,i).real();
    quaternion[0] = V(4,i).real();
    
    double norm = 0.0;
    for( int q = 0; q < 4; q++ )
      norm += pow(quaternion[q],2.0);
    norm = sqrt(norm);
    for( int q = 0; q < 4; q++ )
      quaternion[q] /= norm;
    
    Eigen::Matrix<double,10,1> s;
    upnp_fill_s(quaternion,s);
    Eigen::Matrix<double,1,1> valueM = s.transpose() * M * s + 2.0 * C * s;
    double value = valueM[0] + gamma;

    if( true )//fabs(D[i].imag()) < imagThreshold ) //use all results for the moment
    {
      std::vector<std::pair<double,Eigen::Vector4d>,Eigen::aligned_allocator< std::pair<double,Eigen::Vector4d> > >::iterator
          qidx = quaternions.begin();
      while( qidx != quaternions.end() && qidx->first < value )
        qidx++;
      
      quaternions.insert(qidx,std::pair<double,Eigen::Vector4d>(value,quaternion));
    }
    else
    {
      std::vector<std::pair<double,Eigen::Vector4d>,Eigen::aligned_allocator< std::pair<double,Eigen::Vector4d> > >::iterator
          qidx = bad_quaternions.begin();
      while( qidx != bad_quaternions.end() && qidx->first < value )
        qidx++;
      
      bad_quaternions.insert(qidx,std::pair<double,Eigen::Vector4d>(value,quaternion));
    }
  }
  if( quaternions.size() == 0 )
    quaternions = bad_quaternions;
}
