/******************************************************************************
* Author:   Steffen Urban                                              *
* Contact:  urbste@gmail.com                                          *
* License:  Copyright (c) 2016 Steffen Urban, ANU. All rights reserved.      *
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
 
/* if you use MLPnP consider citing our paper:
@INPROCEEDINGS {mlpnp2016,
	title={MLPNP - A REAL-TIME MAXIMUM LIKELIHOOD SOLUTION TO THE PERSPECTIVE-N-POINT PROBLEM},
	author={Urban, Steffen and Leitloff, Jens and Hinz, Stefan},
	booktitle={ISPRS Annals of Photogrammetry, Remote Sensing \& Spatial Information Sciences},
	pages={131-138},
	year={2016},
	volume={3}
}

29.06.2016 Steffen Urban
*/



#include <opengv/absolute_pose/modules/mlpnp.hpp>
#include <opengv/math/cayley.hpp>
#include <opengv/math/rodrigues.hpp>
#include <Eigen/Sparse>
#include <iostream>


void opengv::absolute_pose::modules::mlpnp::mlpnp_residuals_and_jacs(
	const Eigen::VectorXd& x, 
	const points_t& pts,
	const std::vector<Eigen::MatrixXd>& nullspaces,
	Eigen::VectorXd& r, 
	Eigen::MatrixXd& fjac,
	bool getJacs)
{
	rodrigues_t w(x[0], x[1], x[2]);
	translation_t T(x[3], x[4], x[5]);

	//rotation_t R = math::cayley2rot(c);
	rotation_t R = math::rodrigues2rot(w);
	int ii = 0;

	Eigen::MatrixXd jacs(2, 6);

	for (int i = 0; i < pts.size(); ++i)
	{
		Eigen::Vector3d ptCam = R*pts[i] + T;
		ptCam /= ptCam.norm();

		r[ii] = nullspaces[i].col(0).transpose()*ptCam;
		r[ii + 1] = nullspaces[i].col(1).transpose()*ptCam;
		if (getJacs)
		{
			// jacs
			modules::mlpnp::mlpnpJacs(pts[i],
				nullspaces[i].col(0), nullspaces[i].col(1),
				w, T,
				jacs);

			// r
			fjac(ii, 0) = jacs(0, 0);
			fjac(ii, 1) = jacs(0, 1);
			fjac(ii, 2) = jacs(0, 2);

			fjac(ii, 3) = jacs(0, 3);
			fjac(ii, 4) = jacs(0, 4);
			fjac(ii, 5) = jacs(0, 5);
			// s
			fjac(ii + 1, 0) = jacs(1, 0);
			fjac(ii + 1, 1) = jacs(1, 1);
			fjac(ii + 1, 2) = jacs(1, 2);

			fjac(ii + 1, 3) = jacs(1, 3);
			fjac(ii + 1, 4) = jacs(1, 4);
			fjac(ii + 1, 5) = jacs(1, 5);
			
		}
		ii += 2;
	}
}


void opengv::absolute_pose::modules::mlpnp::mlpnp_lm(
	Eigen::VectorXd& x,
	const points_t& pts,
	const std::vector<Eigen::MatrixXd>& nullspaces,
	const Eigen::SparseMatrix<double> Kll,
	bool use_cov)
{
	const int numObservations = pts.size();
	const int numUnknowns = 6;
	// check redundancy
	assert((2 * numObservations - 6) > 0);

	// =============
	// set all matrices up
	// =============

	Eigen::VectorXd r(2 * numObservations);
	Eigen::VectorXd rd(2 * numObservations);
	Eigen::MatrixXd Jac(2 * numObservations, numUnknowns);
	Eigen::VectorXd g(numUnknowns, 1);
	Eigen::VectorXd dx(numUnknowns, 1); // result vector
	Eigen::MatrixXd eyeMat(numUnknowns, numUnknowns);

	eyeMat.setIdentity();
	Jac.setZero();
	r.setZero();
	dx.setZero();
	g.setZero();

	int it_cnt = 0;
	bool stop = false;
	const int maxIt = 3;
	const double tau = 1e-3;
	double mu = 0.0, nu = 2.0;
	double epsP = 1e-6;
	Eigen::MatrixXd JacTSKll;
	Eigen::MatrixXd A;
	while (it_cnt < maxIt && !stop)
	{
		mlpnp_residuals_and_jacs(x, pts, nullspaces,
			r, Jac, true);
		
		if (use_cov)
			JacTSKll = Jac.transpose() * Kll;
		else
			JacTSKll = Jac.transpose();

		A = JacTSKll * Jac;

		// initalize mu
		if (it_cnt == 0 && !stop)
			mu = tau * sqrt(A.diagonal().maxCoeff());
		// get system matrix
		g = JacTSKll * r;


		A += mu*eyeMat;
		// solve
		Eigen::LDLT<Eigen::MatrixXd> chol(A);
		dx = chol.solve(g);
		//dx = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(g);
		// observation update
		Eigen::MatrixXd dl = Jac * dx;
		if (dl.array().abs().maxCoeff() < epsP || mu < 1e-10)
		{
			if (mu == 0.0)
			{
				it_cnt -= 2;
				stop = true;
				x = x - dx;
				break;
			}
			else
				mu = 0.0;
		}
		else
		{
			Eigen::VectorXd xnew = x - dx;

			rd.setZero();
			mlpnp_residuals_and_jacs(xnew, pts, nullspaces,
				rd, Jac, false);

			// levenberg marquardt
			double epsP = 0.0;
			double Sd = 0.0;
			double dS = 1e-15;
			if (use_cov)
			{
				epsP = r.transpose() * Kll * r;
				Sd = rd.transpose() * Kll * rd;
			}
			else
			{
				epsP = r.transpose() * r;
				Sd = rd.transpose() * rd; //Kll missing!!!
			}

			Eigen::VectorXd dS1 = dx;
			Eigen::VectorXd dS2 = (mu*dS1 + g);

			dS = (dS1.transpose()*dS2);

			double rho = (epsP - Sd) / dS;
			if (rho > 0.0)
			{
				double val1 = 1.0 - pow((2 * rho - 1), 3.0);
				double val2 = 1.0 / 3.0;

				if (val1 > val2)
					mu *= val1;
				else
					mu *= val2;

				nu = 2.0;
				x = xnew;
			}
			else
			{
				mu = mu*nu;
				nu = 2 * nu;
			}
		}
		++it_cnt;
	}//while
	// result
}

void opengv::absolute_pose::modules::mlpnp::mlpnp_gn(
	Eigen::VectorXd& x,
	const points_t& pts,
	const std::vector<Eigen::MatrixXd>& nullspaces,
	const Eigen::SparseMatrix<double> Kll,
	bool use_cov)
{
	const int numObservations = pts.size();
	const int numUnknowns = 6;
	// check redundancy
	assert((2 * numObservations - numUnknowns) > 0);

	// =============
	// set all matrices up
	// =============

	Eigen::VectorXd r(2 * numObservations);
	Eigen::VectorXd rd(2 * numObservations);
	Eigen::MatrixXd Jac(2 * numObservations, numUnknowns);
	Eigen::VectorXd g(numUnknowns, 1);
	Eigen::VectorXd dx(numUnknowns, 1); // result vector

	Jac.setZero();
	r.setZero();
	dx.setZero();
	g.setZero();

	int it_cnt = 0;
	bool stop = false;
	const int maxIt = 5;
	double epsP = 1e-6;

	Eigen::MatrixXd JacTSKll;
	Eigen::MatrixXd A;
	// solve simple gradient descent
	while (it_cnt < maxIt && !stop)
	{
		mlpnp_residuals_and_jacs(x, pts, 
			nullspaces,
			r, Jac, true);
		
		if (use_cov)
			JacTSKll = Jac.transpose() * Kll;
		else
			JacTSKll = Jac.transpose();

		A = JacTSKll * Jac;
		
		// get system matrix
		g = JacTSKll * r;

		// solve
		Eigen::LDLT<Eigen::MatrixXd> chol(A);
		dx = chol.solve(g);
		//dx = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(g);
		// this is to prevent the solution from falling into a wrong minimum
		// if the linear estimate is spurious
		if (dx.array().abs().maxCoeff() > 5.0 || dx.array().abs().minCoeff() > 1.0)
			break;
		// observation update
		Eigen::MatrixXd dl = Jac * dx;
		if (dl.array().abs().maxCoeff() < epsP)
		{
			stop = true;
			x = x - dx;
			break;
		}
		else
			x = x - dx;

		++it_cnt;
	}//while
	// result
}

void opengv::absolute_pose::modules::mlpnp::mlpnp_gn(
	Eigen::VectorXd& x,
	const points_t& pts,
	const std::vector<Eigen::MatrixXd>& nullspaces,
	const Eigen::SparseMatrix<double> Kll,
	Eigen::MatrixXd& Qldld,
	Eigen::MatrixXd& Qxx,
	bool use_cov)
{
	const int numObservations = pts.size();
	const int numUnknowns = 6;
	// check redundancy
	assert((2 * numObservations - 6) > 0);

	// =============
	// set all matrices up
	// =============
	
	Eigen::VectorXd r(2 * numObservations);
	Eigen::VectorXd rd(2 * numObservations);
	Eigen::MatrixXd Jac(2 * numObservations, numUnknowns);
	Eigen::VectorXd g(numUnknowns, 1);
	Eigen::VectorXd dx(numUnknowns, 1); // result vector

	Jac.setZero();
	r.setZero();
	dx.setZero();
	g.setZero();

	int it_cnt = 0;
	bool stop = false;
	const int maxIt = 5;
	double epsP = 1e-6;

	Eigen::MatrixXd JacTSKll;
	Eigen::MatrixXd A;
	// solve simple gradient descent
	while (it_cnt < maxIt && !stop)
	{
		mlpnp_residuals_and_jacs(x, pts,
			nullspaces,
			r, Jac, true);

		if (use_cov)
			JacTSKll = Jac.transpose() * Kll;
		else
			JacTSKll = Jac.transpose();

		A = JacTSKll * Jac;

		// get system matrix
		g = JacTSKll * r;

		// solve
		Eigen::LDLT<Eigen::MatrixXd> chol(A);
		dx = chol.solve(g);
		// this is to prevent the solution from falling into a wrong minimum
		// if the linear estimate is spurious
		if (dx.array().abs().maxCoeff() > 1.0 || dx.array().abs().minCoeff() > 1.0)
			break;
		// observation update
		Eigen::MatrixXd dl = Jac * dx;
		if (dl.array().abs().maxCoeff() < epsP)
		{
			stop = true;
			x = x - dx;
			break;
		}
		else
			x = x - dx;

		++it_cnt;
	}//while
	// statistics
	Qxx = A.inverse();
	Qldld = Jac * Qxx * Jac.transpose();
}