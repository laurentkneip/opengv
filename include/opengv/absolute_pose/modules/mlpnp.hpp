/******************************************************************************
* Author:   Steffen Urban													 *
* Contact:  urbste@gmail.com                                                 *
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


#ifndef OPENGV_ABSOLUTE_POSE_MODULES_MLPNP_HPP
#define OPENGV_ABSOLUTE_POSE_MODULES_MLPNP_HPP

#include <stdlib.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/src/Core/util/DisableStupidWarnings.h>
#include <opengv/types.hpp>

namespace opengv
{
namespace absolute_pose
{
namespace modules
{
namespace mlpnp
{
	void mlpnpJacs(
		const point_t& pt,
		const Eigen::Vector3d& nullspace_r,
		const Eigen::Vector3d& nullspace_s,
		const cayley_t& c,
		const translation_t& t,
		Eigen::MatrixXd& jacs);

	void do_scale(const point_t& pt,
	const rotation_t& rot,
	const translation_t& t,
	double& v1,
	Eigen::Vector2d& scales);

	void mlpnp_lm(Eigen::VectorXd& x,
		const points_t& pts,
		const std::vector<Eigen::MatrixXd>& nullspaces,
		const Eigen::SparseMatrix<double> Kll,
		bool use_cov);

	void mlpnp_gn(Eigen::VectorXd& x,
			const points_t& pts,
			const std::vector<Eigen::MatrixXd>& nullspaces,
			const Eigen::SparseMatrix<double> Kll,
			bool use_cov);

	void mlpnp_gn(Eigen::VectorXd& x,
		const points_t& pts,
		const std::vector<Eigen::MatrixXd>& nullspaces,
		const Eigen::SparseMatrix<double> Kll,
		Eigen::MatrixXd& Qldld,
		Eigen::MatrixXd& Qxx,
		bool use_cov);

	void mlpnp_residuals_and_jacs(
		const Eigen::VectorXd& x,
		const points_t& pts,
		const std::vector<Eigen::MatrixXd>& nullspaces,
		Eigen::VectorXd& r,
		Eigen::MatrixXd& fjac,
		bool getJacs);
}
}
}
}

#endif 
