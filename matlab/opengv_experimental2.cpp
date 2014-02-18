static const char *copyright =
    " Copyright (c) 2013 Laurent Kneip, ANU. All rights reserved.";

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

// Matlab usage:
//
//    X = opengv_experimental2( data1, data2, algorithm )
//
// where
//    data1, data2 are matched points of dimension 6xn
//    algorithm is 0 for sixpt, 1 for ge, and 2 for seventeenpt
//    X is a 3x5 matrix returning the found transformation, plus the number of
//    Ransac-iterations and inliers
//

//matlab header

//standard headers
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "mex.h"

//include generic headers for opengv stuff
#include <opengv/types.hpp>

//include the matlab-adapters
#include <opengv/relative_pose/MANoncentralRelative.hpp>

//expose all ransac-facilities to matlab
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.hpp>

typedef opengv::sac_problems::relative_pose::NoncentralRelativePoseSacProblem nrelRansac;
typedef boost::shared_ptr<nrelRansac> nrelRansacPtr;

// The main mex-function
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{  
  // Characterize the type of the call
  const mxArray *data1 = prhs[0];
  const mxArray *data2 = prhs[1];
  
  const mxArray *temp1 = prhs[2];
  double *temp2 = (double*) mxGetData(temp1);
  int algorithm = floor(temp2[0]+0.01);
  
  const mwSize *data1dim = mxGetDimensions(data1);
  const mwSize *data2dim = mxGetDimensions(data2);
  
  //create three pointers to absolute, relative, and point_cloud adapters here
  opengv::relative_pose::RelativeAdapterBase* relativeAdapter =
      new opengv::relative_pose::MANoncentralRelative(
      (double*) mxGetData(data1),
      (double*) mxGetData(data2),
      data1dim[1],
      data2dim[1] );

  nrelRansacPtr problem;
  
  switch(algorithm)
  {
    case 0:
	  problem = nrelRansacPtr( new nrelRansac( *relativeAdapter, nrelRansac::SIXPT ) );
	  break;
	case 1:
	  problem = nrelRansacPtr( new nrelRansac( *relativeAdapter, nrelRansac::GE ) );
	  break;
	case 2:
	  problem = nrelRansacPtr( new nrelRansac( *relativeAdapter, nrelRansac::SEVENTEENPT ) );
	  break;
  }
  
  opengv::sac::Ransac<nrelRansac> ransac;
  ransac.sac_model_ = problem;
  ransac.threshold_ = 2.0*(1.0 - cos(atan(sqrt(2.0)*0.5/800.0)));
  ransac.max_iterations_ = 10000000;
  ransac.computeModel();
  
  Eigen::Matrix<double,3,5> result;
  result.block<3,4>(0,0) = ransac.model_coefficients_;
  result(0,4) = ransac.iterations_;
  result(1,4) = ransac.inliers_.size();
  
  int dims[2];
  dims[0] = 3;
  dims[1] = 5;
  plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  memcpy(mxGetData(plhs[0]), result.data(), 15*sizeof(double));
}
