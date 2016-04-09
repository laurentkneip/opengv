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

/**
 * \file MANoncentralAbsolute.hpp
 * \brief Adapter-class for passing bearing-vector-to-point correspondences to
 *        the non-central absolute-pose algorithms. It maps matlab
 *        types to opengv types.
 */

#ifndef OPENGV_ABSOLUTE_POSE_MANONCENTRALABSOLUTE_HPP_
#define OPENGV_ABSOLUTE_POSE_MANONCENTRALABSOLUTE_HPP_

#include <stdlib.h>
#include <vector>
#include <opengv/types.hpp>
#include <opengv/absolute_pose/AbsoluteAdapterBase.hpp>

/**
 * \brief The namespace of this library.
 */
namespace opengv
{
/**
 * \brief The namespace for the absolute pose methods.
 */
namespace absolute_pose
{

/**
 * Check the documentation of the parent-class to understand the meaning of
 * an AbsoluteAdapter. This child-class is for the non-central case and holds
 * data in form of pointers to matlab-data.
 */
class MANoncentralAbsolute : public AbsoluteAdapterBase
{
protected:
  using AbsoluteAdapterBase::_t;
  using AbsoluteAdapterBase::_R;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /**
   * \brief Constructor. See protected class-members to understand parameters
   */
  MANoncentralAbsolute(
      const double * points,
      const double * bearingVectors,
      int numberPoints,
      int numberBearingVectors );
  /**
  * \brief Constructor. See protected class-members to understand parameters
  */
  // edited by Steffen Urban / urbste@gmail.com
  MANoncentralAbsolute(
	  const double * points,
	  const double * bearingVectors,
	  const double * covMats,
	  int numberPoints,
	  int numberBearingVectors,
	  int numberCovMats);
  /**
   * Destructor
   */
  virtual ~MANoncentralAbsolute();

  //Access of correspondences
  
  /** See parent-class */
  virtual opengv::bearingVector_t getBearingVector(size_t index) const;
  /** See parent-class */
  virtual opengv::bearingVectors_t getBearingVectors() const;
  /** See parent-class */
  virtual opengv::cov3_mat_t getCovariance(size_t index) const;
  /** See parent-class */
  virtual opengv::cov3_mats_t getCovariances() const;
  /** See parent-class */
  virtual double getWeight( size_t index ) const;
  /** See parent-class */
  virtual opengv::translation_t getCamOffset( size_t index ) const;
  /** See parent-class */
  virtual opengv::rotation_t getCamRotation( size_t index ) const;
  /** See parent-class */
  virtual opengv::point_t getPoint( size_t index ) const;
  /** See parent-class */
  virtual opengv::points_t getPoints() const;
  /** See parent-class */
  virtual size_t getNumberCorrespondences() const;

protected:

  /** A pointer to the point data */
  const double * _points;
  /** A pointer to the bearing-vectors */
  const double * _bearingVectors;
  /** A pointer to the covariance-vectors */
  const double * _covMats;
  /** The number of points */
  int _numberPoints;
  /** The number of bearing vectors */
  int _numberBearingVectors;
  /** The number of covariance matrices */
  int _numberCovMats;
};

}
}

#endif /* OPENGV_ABSOLUTE_POSE_MANONCENTRALABSOLUTE_HPP_ */
