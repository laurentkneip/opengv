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


#include <opengv/absolute_pose/NoncentralAbsoluteAdapter.hpp>

opengv::absolute_pose::NoncentralAbsoluteAdapter::NoncentralAbsoluteAdapter(
    const bearingVectors_t & bearingVectors,
    const camCorrespondences_t & camCorrespondences,
    const points_t & points,
    const translations_t & camOffsets,
    const rotations_t & camRotations,
	const cov3_mats_t & covMats) :
    AbsoluteAdapterBase(),
    _bearingVectors(bearingVectors),
    _camCorrespondences(camCorrespondences),
    _points(points),
    _camOffsets(camOffsets),
    _camRotations(camRotations),
	_cov_mats(covMats)
{}

opengv::absolute_pose::NoncentralAbsoluteAdapter::NoncentralAbsoluteAdapter(
    const bearingVectors_t & bearingVectors,
    const camCorrespondences_t & camCorrespondences,
    const points_t & points,
    const translations_t & camOffsets,
    const rotations_t & camRotations,
    const rotation_t & R,
	const cov3_mats_t & covMats) :
    AbsoluteAdapterBase(R),
    _bearingVectors(bearingVectors),
    _camCorrespondences(camCorrespondences),
    _points(points),
    _camOffsets(camOffsets),
    _camRotations(camRotations),
	_cov_mats(covMats)
{}

opengv::absolute_pose::NoncentralAbsoluteAdapter::NoncentralAbsoluteAdapter(
    const bearingVectors_t & bearingVectors,
    const camCorrespondences_t & camCorrespondences,
    const points_t & points,
    const translations_t & camOffsets,
    const rotations_t & camRotations,
    const translation_t & t,
    const rotation_t & R,
	const cov3_mats_t & covMats) :
    AbsoluteAdapterBase(t,R),
    _bearingVectors(bearingVectors),
    _camCorrespondences(camCorrespondences),
    _points(points),
    _camOffsets(camOffsets),
    _camRotations(camRotations),
	_cov_mats(covMats)
{}

opengv::absolute_pose::NoncentralAbsoluteAdapter::~NoncentralAbsoluteAdapter()
{}

opengv::bearingVector_t
opengv::absolute_pose::NoncentralAbsoluteAdapter::
    getBearingVector( size_t index ) const
{
  assert(index < _bearingVectors.size());
  return _bearingVectors[index];
}
// edited by Steffen Urban / urbste@gmail.com
opengv::bearingVectors_t
opengv::absolute_pose::NoncentralAbsoluteAdapter::
getBearingVectors() const
{
	return _bearingVectors;
}
// edited by Steffen Urban / urbste@gmail.com
opengv::cov3_mat_t
opengv::absolute_pose::NoncentralAbsoluteAdapter::
getCovariance(size_t index) const
{
	assert(index < _bearingVectors.size());
	return _cov_mats[index];
}
// edited by Steffen Urban / urbste@gmail.com
opengv::cov3_mats_t
opengv::absolute_pose::NoncentralAbsoluteAdapter::
getCovariances() const
{
	return _cov_mats;
}

double
opengv::absolute_pose::NoncentralAbsoluteAdapter::
    getWeight( size_t index ) const
{
  return 1.0;
}

opengv::point_t
opengv::absolute_pose::NoncentralAbsoluteAdapter::
    getPoint( size_t index ) const
{
  assert(index < _bearingVectors.size());
  return _points[index];
}
// edited by Steffen Urban / urbste@gmail.com
opengv::points_t
opengv::absolute_pose::NoncentralAbsoluteAdapter::
getPoints() const
{
	return _points;
}

opengv::translation_t
opengv::absolute_pose::NoncentralAbsoluteAdapter::
    getCamOffset( size_t index ) const
{
  assert(_camCorrespondences[index] < _camOffsets.size());
  return _camOffsets[_camCorrespondences[index]];
}

opengv::rotation_t
opengv::absolute_pose::NoncentralAbsoluteAdapter::
    getCamRotation( size_t index ) const
{
  assert(_camCorrespondences[index] < _camRotations.size());
  return _camRotations[_camCorrespondences[index]];
}

size_t
opengv::absolute_pose::NoncentralAbsoluteAdapter::
    getNumberCorrespondences() const
{
  return _bearingVectors.size();
}
