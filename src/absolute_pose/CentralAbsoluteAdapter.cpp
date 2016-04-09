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


#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>


opengv::absolute_pose::CentralAbsoluteAdapter::CentralAbsoluteAdapter(
    const bearingVectors_t & bearingVectors,
    const points_t & points,
	const cov3_mats_t & covMats) :
    AbsoluteAdapterBase(),
    _bearingVectors(bearingVectors),
    _points(points),
	_cov_mats(covMats)
{}

opengv::absolute_pose::CentralAbsoluteAdapter::CentralAbsoluteAdapter(
    const bearingVectors_t & bearingVectors,
    const points_t & points,
	const rotation_t & R,
	const cov3_mats_t & covMats) :
    AbsoluteAdapterBase(R),
    _bearingVectors(bearingVectors),
	_points(points),
	_cov_mats(covMats)
{}

opengv::absolute_pose::CentralAbsoluteAdapter::CentralAbsoluteAdapter(
    const bearingVectors_t & bearingVectors,
    const points_t & points,
    const translation_t & t,
	const rotation_t & R,
	const cov3_mats_t & covMats) :
    AbsoluteAdapterBase(t,R),
    _bearingVectors(bearingVectors),
	_points(points),
	_cov_mats(covMats)
{}

opengv::absolute_pose::CentralAbsoluteAdapter::~CentralAbsoluteAdapter()
{}

opengv::bearingVector_t
opengv::absolute_pose::CentralAbsoluteAdapter::getBearingVector(
    size_t index ) const
{
  assert(index < _bearingVectors.size());
  return _bearingVectors[index];
}

opengv::bearingVectors_t
opengv::absolute_pose::CentralAbsoluteAdapter::getBearingVectors() const
{
	return _bearingVectors;
}

double
opengv::absolute_pose::CentralAbsoluteAdapter::
    getWeight( size_t index ) const
{
  return 1.0;
}

opengv::point_t
opengv::absolute_pose::CentralAbsoluteAdapter::getPoint(
    size_t index ) const
{
  assert(index < _bearingVectors.size());
  return _points[index];
}

opengv::points_t
opengv::absolute_pose::CentralAbsoluteAdapter::getPoints() const
{
	return _points;
}

opengv::translation_t
opengv::absolute_pose::CentralAbsoluteAdapter::getCamOffset(
    size_t index ) const
{
  //we could insert a check here that camIndex is 0, because this adapter is
  //for a single camera only
  return Eigen::Vector3d::Zero();
}

opengv::rotation_t
opengv::absolute_pose::CentralAbsoluteAdapter::getCamRotation(
    size_t index ) const
{
  //we could insert a check here that camIndex is 0, because this adapter is
  //for a single camera only
  return Eigen::Matrix3d::Identity();
}

size_t
opengv::absolute_pose::CentralAbsoluteAdapter::
    getNumberCorrespondences() const
{
  return _bearingVectors.size();
}

// covariances
opengv::cov3_mat_t
opengv::absolute_pose::CentralAbsoluteAdapter::getCovariance(
size_t index) const
{
	assert(index < _bearingVectors.size());
	return _cov_mats[index];
}

opengv::cov3_mats_t
opengv::absolute_pose::CentralAbsoluteAdapter::getCovariances() const
{
	return _cov_mats;
}