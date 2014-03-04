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
//    X = opengv ( method, data1, data2 )
//    X = opengv ( method, indices, data1, data2 )
//    X = opengv ( method, indices, data1, data2, prior )
//
// where
//    method is a string that characterizes the algorithm to use
//    data1, data2 are matched points (each one of dimension 3xn or 6xn)
//    X is a 3xnxm matrix, where n is the second dimensionality of the solution space,
//    and m is the number of solutions
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
#include <opengv/absolute_pose/MACentralAbsolute.hpp>
#include <opengv/absolute_pose/MANoncentralAbsolute.hpp>
#include <opengv/relative_pose/MACentralRelative.hpp>
#include <opengv/relative_pose/MANoncentralRelative.hpp>
#include <opengv/point_cloud/MAPointCloud.hpp>

//expose all methods to matlab
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/point_cloud/methods.hpp>

//expose all ransac-facilities to matlab
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac_problems/absolute_pose/AbsolutePoseSacProblem.hpp>
#include <opengv/sac_problems/relative_pose/CentralRelativePoseSacProblem.hpp>
#include <opengv/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.hpp>
#include <opengv/sac_problems/relative_pose/EigensolverSacProblem.hpp>
#include <opengv/sac_problems/relative_pose/RotationOnlySacProblem.hpp>
#include <opengv/sac_problems/point_cloud/PointCloudSacProblem.hpp>

// The different methods that can be used within Matlab
static const char* methods[]=
{
  // absolute_pose methods
  "p2p",                      //  3
  "p3p_kneip",                //  9
  "p3p_gao",                  //  7
  "epnp",                     //  4
  "p3p_kneip_ransac",         // 16
  "p3p_gao_ransac",           // 14
  "epnp_ransac",              // 11
  "abs_nonlin_central",       // 18
  "gp3p",                     //  4
  "gp3p_ransac",              // 11
  "gpnp",                     //  4
  "abs_nonlin_noncentral",    // 21

  // relative_pose methods
  "twopt",                    //  5
  "twopt_rotationOnly",       // 18
  "rotationOnly",             // 12
  "fivept_stewenius",         // 16
  "fivept_nister",            // 13
  "fivept_kneip",             // 12
  "sevenpt",                  //  7
  "eightpt",                  //  7
  "eigensolver",              // 11
  "rotationOnly_ransac",      // 19
  "fivept_stewenius_ransac",  // 23
  "fivept_nister_ransac",     // 20
  "sevenpt_ransac",           // 14
  "eightpt_ransac",           // 14
  "eigensolver_ransac",       // 18
  "rel_nonlin_central",       // 18
  "sixpt",                    //  5
  "seventeenpt",              // 11
  "ge",                       //  2
  "sixpt_ransac",             // 12
  "seventeenpt_ransac",       // 18
  "ge_ransac",                //  9
  "rel_nonlin_noncentral",    // 21

  // point_cloud methods
  "threept_arun",             // 12
  "threept_arun_ransac"       // 19
};

// The length of the method strings (needed for comparison)
static const int methodsLengths[] =
    { 3,9,7,4,16,14,11,18,4,11,4,21,5,18,12,16,13,12,
	  7,7,11,19,23,20,14,14,18,18,5,11,2,12,18,9,21,12,19 };

static const int absCentralFirst    =  0;
static const int absCentralLast     =  7;
static const int absNoncentralFirst =  8;
static const int absNoncentralLast  = 11;
static const int relCentralFirst    = 12;
static const int relCentralLast     = 27;
static const int relNoncentralFirst = 28;
static const int relNoncentralLast  = 34;
static const int pointCloudFirst    = 35;
static const int pointCloudLast     = 36;

// The number of methods (needed for comparison)
static const int numberMethods = pointCloudLast + 1;

enum Method
{
  P2P,
  P3P_KNEIP,
  P3P_GAO,
  EPNP,
  P3P_KNEIP_RANSAC,
  P3P_GAO_RANSAC,
  EPNP_RANSAC,
  ABS_NONLIN_CENTRAL,
  GP3P,
  GP3P_RANSAC,
  GPNP,
  ABS_NONLIN_NONCENTRAL,
  TWOPT,
  TWOPT_ROTATIONONLY,
  ROTATIONONLY,
  FIVEPT_STEWENIUS,
  FIVEPT_NISTER,
  FIVEPT_KNEIP,
  SEVENPT,
  EIGHTPT,
  EIGENSOLVER,
  ROTATIONONLY_RANSAC,
  FIVEPT_STEWENIUS_RANSAC,
  FIVEPT_NISTER_RANSAC,
  SEVENPT_RANSAC,
  EIGHTPT_RANSAC,
  EIGENSOLVER_RANSAC,
  REL_NONLIN_CENTRAL,
  SIXPT,
  SEVENTEENPT,
  GE,
  SIXPT_RANSAC,
  SEVENTEENPT_RANSAC,
  GE_RANSAC,
  REL_NONLIN_NONCENTRAL,
  THREEPT_ARUN,
  THREEPT_ARUN_RANSAC
};

// Finds the method based on string comparison
int findCase( const char* input, int inputLength )
{
  int n = 0;
  while( n < numberMethods )
  {
    // First check the length of the string, it needs to be same
    if( inputLength == methodsLengths[n])
    {
      // Now check if all the elements are the same
      int allSame = 1;
      for( int i = 0; i < inputLength; i++ )
      {
        if( input[i] != methods[n][i] )
        {
          allSame = 0;
          break;
        }
      }
      
      // Break if method found
      if( allSame )
        return n;
      
      //Otherwise go on with the next one
    }
	
	n++;
  }
  
  // Return -1 if not found
  return -1;
}

// Print all possible cases
void printCases()
{
  mexPrintf("The known methods are:");
  
  for( int i = 0; i < numberMethods; i++ )
  {
    mexPrintf("\n");
    mexPrintf(methods[i]);
  }
  mexPrintf("\n");
}

typedef opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem absRansac;
typedef boost::shared_ptr<absRansac> absRansacPtr;

typedef opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem relRansac;
typedef boost::shared_ptr<relRansac> relRansacPtr;
typedef opengv::sac_problems::relative_pose::NoncentralRelativePoseSacProblem nrelRansac;
typedef boost::shared_ptr<nrelRansac> nrelRansacPtr;
typedef opengv::sac_problems::relative_pose::RotationOnlySacProblem rotRansac;
typedef boost::shared_ptr<rotRansac> rotRansacPtr;
typedef opengv::sac_problems::relative_pose::EigensolverSacProblem eigRansac;
typedef boost::shared_ptr<eigRansac> eigRansacPtr;

typedef opengv::sac_problems::point_cloud::PointCloudSacProblem ptRansac;
typedef boost::shared_ptr<ptRansac> ptRansacPtr;

// The main mex-function
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  // Check if right number of arguments
  if( nrhs < 3 || nrhs > 5 )
  {
    mexPrintf("opengv: Not an acceptable number of arguments\n");
    mexPrintf("Usage:  X = opengv( method, data1, data2 )\n");
    mexPrintf("Or:     X = opengv( method, indices, data1, data2 )\n");
    mexPrintf("Or:     X = opengv( method, indices, data1, data2, prior )\n");
    return;
  }
  
  // Get the method
  if( mxGetM(prhs[0]) != 1 )
  {
    mexPrintf("opengv: Bad input to mex function opengv\n");
    mexPrintf("Usage:  X = opengv( method, data1, data2 )\n");
    mexPrintf("Or:     X = opengv( method, indices, data1, data2 )\n");
    mexPrintf("Or:     X = opengv( method, indices, data1, data2, prior )\n");
    mexPrintf("Hint:   Method must be a string\n");
    return;
  }

  // Now get the string and find the caseNumber
  mwSize strlen = (mwSize) mxGetN(prhs[0]) + 1;  
  char * method = (char *) malloc(strlen);
  mxGetString(prhs[0], method, strlen);
  int caseNumber = findCase(method, (int) mxGetN(prhs[0]));
  
  // Return if method not found
  if( caseNumber < 0 )
  {
    mexPrintf("opengv: Unknown method\n");
    printCases();
    return;
  }
  
  // Characterize the type of the call
  int callCharacter = -1;
  const mxArray *data1;
  const mxArray *data2;
  const mwSize *data1dim;
  const mwSize *data2dim;
  
  if( nrhs == 3 ) // X = opengv( method, data1, data2 )
  {
    // Check the input
    data1 = prhs[1];
    data2 = prhs[2];

    // Check the dimensions of the arguments
    int ndimensions1 = mxGetNumberOfDimensions(data1);
    int ndimensions2 = mxGetNumberOfDimensions(data2);
    data1dim = mxGetDimensions(data1);
    data2dim = mxGetDimensions(data2);
    
    // Now check them
    if( ndimensions1 != 2 || ndimensions2 != 2 ||
        (data1dim[0] != 3 && data1dim[0] != 6) ||
        (data2dim[0] != 3 && data2dim[0] != 6) ||
        data1dim[1] != data2dim[1] ||
        data1dim[1] < 1 || data2dim[1] < 1 )
    {
      mexPrintf("opengv: Bad input to mex function\n");
      mexPrintf("Assuming signature: X = opengv( method, data1, data2 )\n");
      mexPrintf("Inputs data1 and data2 must have size (3,n) or (6,n),\n");
      mexPrintf("with an equal number of columns\n");
      return;
    }
    
    callCharacter = 0;
  }
  if( nrhs == 4 )
  {
    // X = opengv( method, indices, data1, data2 )

    // Check the input
    data1 = prhs[2];
    data2 = prhs[3];

    // Check the dimensions of the arguments
    int ndimensions1 = mxGetNumberOfDimensions(data1);
    int ndimensions2 = mxGetNumberOfDimensions(data2);
    int ndimensions3 = mxGetNumberOfDimensions(prhs[1]);
    data1dim = mxGetDimensions(data1);
    data2dim = mxGetDimensions(data2);
    const mwSize *indicesDim = mxGetDimensions(prhs[1]);
  
    // Now check them
    if( ndimensions1 != 2 || ndimensions2 != 2 || ndimensions3 != 2 ||
        (data1dim[0] != 3 && data1dim[0] != 6) ||
        (data2dim[0] != 3 && data2dim[0] != 6) ||
        indicesDim[0] != 1 ||
        data1dim[1] != data2dim[1] ||
        data1dim[1] < 1 || data2dim[1] < 1 ||
        data2dim[1] < indicesDim[1] )
    {
      mexPrintf("opengv: Bad input to mex function opengv\n");
      mexPrintf("Assuming signature: X = opengv( method, indices, data1, ");
      mexPrintf("data2 )\n");
      mexPrintf("Inputs data1 and data2 must have size (3,n) or (6,n),\n");
      mexPrintf("with an equal number of columns\n");
      mexPrintf("indices must be a 1xm vector, with m smaller or equal than n\n");
      return;
    }
    
    callCharacter = 1;
  }
  if(nrhs == 5)
  {
    // X = opengv( method, indices, data1, data2, prior )
    
    // Check the input
    data1 = prhs[2];
    data2 = prhs[3];

    // Check the dimensions of the arguments
    int ndimensions1 = mxGetNumberOfDimensions(data1);
    int ndimensions2 = mxGetNumberOfDimensions(data2);
    int ndimensions3 = mxGetNumberOfDimensions(prhs[1]);
    int ndimensions4 = mxGetNumberOfDimensions(prhs[4]);
    data1dim = mxGetDimensions(data1);
    data2dim = mxGetDimensions(data2);
    const mwSize *indicesDim = mxGetDimensions(prhs[1]);
    const mwSize *priorDim = mxGetDimensions(prhs[4]);
  
    // Now check them
    if( ndimensions1 != 2 || ndimensions2 != 2 || ndimensions3 != 2 || ndimensions4 != 2 ||
        (data1dim[0] != 3 && data1dim[0] != 6) ||
        (data2dim[0] != 3 && data2dim[0] != 6) ||
        indicesDim[0] != 1 ||
        priorDim[0] != 3 ||
        (priorDim[1] != 1 &&  priorDim[1] != 3 && priorDim[1] != 4) ||
        data1dim[1] != data2dim[1] ||
        data1dim[1] < 1 || data2dim[1] < 1 ||
        data2dim[1] < indicesDim[1] )
    {
      mexPrintf("opengv: Bad input to mex function opengv\n");
      mexPrintf("Assuming signature: X = opengv( method, indices, data1, ");
      mexPrintf("data2, prior )\n");
      mexPrintf("Inputs data1 and data2 must have size (3,n) or (6,n),\n");
      mexPrintf("with an equal number of columns\n");
      mexPrintf("indices must be a 1xm vector, with m smaller or equal than n\n");
      mexPrintf("prior must be a 3x1, 3x3, or 3x4 matrix\n");
      return;
    }
  
    callCharacter = 2;
  }
  
  //create three pointers to absolute, relative, and point_cloud adapters here
  opengv::absolute_pose::AbsoluteAdapterBase* absoluteAdapter;
  opengv::relative_pose::RelativeAdapterBase* relativeAdapter;
  opengv::point_cloud::PointCloudAdapterBase* pointCloudAdapter;
  
  int translationPrior = 0;
  int rotationPrior = 0;
  opengv::translation_t translation;
  opengv::rotation_t rotation;
  
  //set the prior if needed
  if( callCharacter == 2 )
  {
    const mxArray *prior;
    const mwSize *priorDim;
    
    prior = prhs[4];
    priorDim = mxGetDimensions(prhs[4]);
    
    if( priorDim[1] == 1 )
    {
      //set translation
      translationPrior = 1;
      double * ptr = (double*) mxGetData(prior);
      translation[0] = ptr[0];
      translation[1] = ptr[1];
      translation[2] = ptr[2];
    }
    if( priorDim[1] == 3 )
    {
      //set rotation
      rotationPrior = 1;
      double * ptr = (double*) mxGetData(prior);
      rotation(0,0) = ptr[0];
      rotation(1,0) = ptr[1];
      rotation(2,0) = ptr[2];
      rotation(0,1) = ptr[3];
      rotation(1,1) = ptr[4];
      rotation(2,1) = ptr[5];
      rotation(0,2) = ptr[6];
      rotation(1,2) = ptr[7];
      rotation(2,2) = ptr[8];
    }
    if( priorDim[1] == 4 )
    {
      translationPrior = 1;
      rotationPrior = 1;
      double * ptr = (double*) mxGetData(prior);
      rotation(0,0) = ptr[0];
      rotation(1,0) = ptr[1];
      rotation(2,0) = ptr[2];
      rotation(0,1) = ptr[3];
      rotation(1,1) = ptr[4];
      rotation(2,1) = ptr[5];
      rotation(0,2) = ptr[6];
      rotation(1,2) = ptr[7];
      rotation(2,2) = ptr[8];
      translation[0] = ptr[9];
      translation[1] = ptr[10];
      translation[2] = ptr[11];
    }
  }
  
  if( caseNumber >= absCentralFirst && caseNumber <= absCentralLast )
  {
    //central absolute case
    if( data1dim[0] != 3 || data2dim[0] != 3 )
    {
      mexPrintf("opengv: Bad input to mex function opengv\n");
      mexPrintf("Assuming method: ");
      mexPrintf(methods[caseNumber]);
      mexPrintf("\n");
      mexPrintf("Inputs data1 and data2 must have size (3,n) for a central ");
      mexPrintf("absolute method\n");
      return;
    }

    absoluteAdapter = new opengv::absolute_pose::MACentralAbsolute(
        (double*) mxGetData(data1),
        (double*) mxGetData(data2),
        data1dim[1],
        data2dim[1] );
    
    if( translationPrior == 1 )
      absoluteAdapter->sett(translation);
    if( rotationPrior == 1 )
      absoluteAdapter->setR(rotation);
  }
  
  if(caseNumber >= absNoncentralFirst && caseNumber <= absNoncentralLast )
  {
    //non-central absolute case    
    if( data1dim[0] != 3 || data2dim[0] != 6 )
    {
      mexPrintf("opengv: Bad input to mex function opengv\n");
      mexPrintf("Assuming method: ");
      mexPrintf(methods[caseNumber]);
      mexPrintf("\n");
      mexPrintf("Inputs data1 and data2 must have sizes (3,n) and (6,n) for ");
      mexPrintf("a noncentral absolute method\n");
      return;
    }
    
    absoluteAdapter = new opengv::absolute_pose::MANoncentralAbsolute(
        (double*) mxGetData(data1),
        (double*) mxGetData(data2),
        data1dim[1],
        data2dim[1] );
    
    if( translationPrior == 1 )
      absoluteAdapter->sett(translation);
    if( rotationPrior == 1 )
      absoluteAdapter->setR(rotation);
  }
  if(caseNumber >= relCentralFirst && caseNumber <= relCentralLast )
  {
    //central relative case    
    if( data1dim[0] != 3 || data2dim[0] != 3 )
    {
      mexPrintf("opengv: Bad input to mex function opengv\n");
      mexPrintf("Assuming method: ");
      mexPrintf(methods[caseNumber]);
      mexPrintf("\n");
      mexPrintf("Inputs data1 and data2 must have size (3,n) for a central ");
      mexPrintf("relative method\n");
      return;
    }
    
    relativeAdapter = new opengv::relative_pose::MACentralRelative(
        (double*) mxGetData(data1),
        (double*) mxGetData(data2),
        data1dim[1],
        data2dim[1] );
    
    if( translationPrior == 1 )
      relativeAdapter->sett12(translation);
    if( rotationPrior == 1 )
      relativeAdapter->setR12(rotation);
  }
  
  if(caseNumber >= relNoncentralFirst && caseNumber <= relNoncentralLast )
  {
    //noncentral relative case    
    if( data1dim[0] != 6 || data2dim[0] != 6 )
    {
      mexPrintf("opengv: Bad input to mex function opengv\n");
      mexPrintf("Assuming method: ");
      mexPrintf(methods[caseNumber]);
      mexPrintf("\n");
      mexPrintf("Inputs data1 and data2 must have size (6,n) for a ");
      mexPrintf("noncentral relative method\n");
      return;
    }
    
    relativeAdapter = new opengv::relative_pose::MANoncentralRelative(
        (double*) mxGetData(data1),
        (double*) mxGetData(data2),
        data1dim[1],
        data2dim[1] );
    
    if( translationPrior == 1 )
      relativeAdapter->sett12(translation);
    if( rotationPrior == 1 )
      relativeAdapter->setR12(rotation);
  }
    
  if(caseNumber >= pointCloudFirst && caseNumber <= pointCloudLast )
  {
    //point-cloud case    
    if( data1dim[0] != 3 || data2dim[0] != 3 )
    {
      mexPrintf("opengv: Bad input to mex function opengv\n");
      mexPrintf("Assuming method: ");
      mexPrintf(methods[caseNumber]);
      mexPrintf("\n");
      mexPrintf("Inputs data1 and data2 must have size (3,n) for a ");
      mexPrintf("point-cloud method\n");
      return;
    }
    
    pointCloudAdapter = new opengv::point_cloud::MAPointCloud(
        (double*) mxGetData(data1),
        (double*) mxGetData(data2),
        data1dim[1],
        data2dim[1] );
    
    if( translationPrior == 1 )
      pointCloudAdapter->sett12(translation);
    if( rotationPrior == 1 )
      pointCloudAdapter->setR12(rotation);
  }
  
  //check if a return argument is needed, otherwise we won't start computing
  if( nlhs != 1 )
  {
    if( nlhs > 1 )
      mexPrintf("opengv: Returns one parameter only\n");
    return;
  }

  //create the indices array (todo: check if there is a smarter way for doing this)
  std::vector<int> indices;
  int useIndices = 0;
  if( callCharacter > 0 )
  {
    useIndices = 1;
    const mwSize *indicesDim = mxGetDimensions(prhs[1]);
    int numberOfIndices = indicesDim[1];
    indices.reserve(numberOfIndices);
    double * mxIndices = (double*) mxGetData(prhs[1]);
    for( int i = 0; i < numberOfIndices; i++ )
      indices.push_back(floor(mxIndices[i]+0.01)-1);
  }
  
  Method methodEnum = static_cast<Method>(caseNumber);
  if( caseNumber != (int) methodEnum )
  {
    mexPrintf("opengv: This method is not yet implemented!\n");
    return;
  }
  
  // Finally, call the respective algorithm
  switch (methodEnum)
  {
    case P2P:
    {
      opengv::translation_t temp;
      if(useIndices)
        temp = opengv::absolute_pose::p2p(*absoluteAdapter,indices);
      else
        temp = opengv::absolute_pose::p2p(*absoluteAdapter);
      int dims[2];
      dims[0] = 3;
      dims[1] = 1;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), temp.data(), 3*sizeof(double));
      break;
    }
    case P3P_KNEIP:
    {
      opengv::transformations_t temp;
      if(useIndices)
        temp = opengv::absolute_pose::p3p_kneip(*absoluteAdapter,indices);
      else
        temp = opengv::absolute_pose::p3p_kneip(*absoluteAdapter);
      int dims[3];
      dims[0] = 3;
      dims[1] = 4;
      dims[2] = temp.size();
      plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      for( int i = 0; i < temp.size(); i++ )
      {
        void * targetAddress = ((char*) mxGetData(plhs[0])) + i*12*sizeof(double);
        memcpy(targetAddress, temp[i].data(), 12*sizeof(double));
      }
      break;
    }
    case P3P_GAO:
    {
      opengv::transformations_t temp;
      if(useIndices)
        temp = opengv::absolute_pose::p3p_gao(*absoluteAdapter,indices);
      else
        temp = opengv::absolute_pose::p3p_gao(*absoluteAdapter);
      int dims[3];
      dims[0] = 3;
      dims[1] = 4;
      dims[2] = temp.size();
      plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      for( int i = 0; i < temp.size(); i++ )
      {
        void * targetAddress = ((char*) mxGetData(plhs[0])) + i*12*sizeof(double);
        memcpy(targetAddress, temp[i].data(), 12*sizeof(double));
      }
      break;
    }
    case EPNP:
    {
      opengv::transformation_t temp;
      if(useIndices)
        temp = opengv::absolute_pose::epnp(*absoluteAdapter,indices);
      else
        temp = opengv::absolute_pose::epnp(*absoluteAdapter);
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), temp.data(), 12*sizeof(double));
      break;
    }
    case P3P_KNEIP_RANSAC:
    {
      absRansacPtr problem;
      if(useIndices)
        problem = absRansacPtr( new absRansac( *absoluteAdapter, absRansac::KNEIP, indices ) );
      else
        problem = absRansacPtr( new absRansac( *absoluteAdapter, absRansac::KNEIP ) );
      opengv::sac::Ransac<absRansac> ransac;
      ransac.sac_model_ = problem;
      ransac.threshold_ = 1.0 - cos(atan(sqrt(2.0)*0.5/800.0));
      ransac.max_iterations_ = 50;
      ransac.computeModel();
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), ransac.model_coefficients_.data(), 12*sizeof(double));
      break;
    }
    case P3P_GAO_RANSAC:
    {
      absRansacPtr problem;
      if(useIndices)
        problem = absRansacPtr( new absRansac( *absoluteAdapter, absRansac::GAO, indices ) );
      else
        problem = absRansacPtr( new absRansac( *absoluteAdapter, absRansac::GAO ) );
      opengv::sac::Ransac<absRansac> ransac;
      ransac.sac_model_ = problem;
      ransac.threshold_ = 1.0 - cos(atan(sqrt(2.0)*0.5/800.0));
      ransac.max_iterations_ = 50;
      ransac.computeModel();
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), ransac.model_coefficients_.data(), 12*sizeof(double));
      break;
    }
    case EPNP_RANSAC:
    {
      absRansacPtr problem;
      if(useIndices)
        problem = absRansacPtr( new absRansac( *absoluteAdapter, absRansac::EPNP, indices ) );
      else
        problem = absRansacPtr( new absRansac( *absoluteAdapter, absRansac::EPNP ) );
      opengv::sac::Ransac<absRansac> ransac;
      ransac.sac_model_ = problem;
      ransac.threshold_ = 1.0 - cos(atan(sqrt(2.0)*0.5/800.0));
      ransac.max_iterations_ = 50;
      ransac.computeModel();
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), ransac.model_coefficients_.data(), 12*sizeof(double));
      break;
    }
    case ABS_NONLIN_CENTRAL:
    {
      opengv::transformation_t temp;
      if(useIndices)
        temp = opengv::absolute_pose::optimize_nonlinear(*absoluteAdapter,indices);
      else
        temp = opengv::absolute_pose::optimize_nonlinear(*absoluteAdapter);
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), temp.data(), 12*sizeof(double));
      break;
    }
    case GP3P:
    {
      opengv::transformations_t temp;
      if(useIndices)
        temp = opengv::absolute_pose::gp3p(*absoluteAdapter,indices);
      else
        temp = opengv::absolute_pose::gp3p(*absoluteAdapter);
      int dims[3];
      dims[0] = 3;
      dims[1] = 4;
      dims[2] = temp.size();
      plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      for( int i = 0; i < temp.size(); i++ )
      {
        void * targetAddress = ((char*) mxGetData(plhs[0])) + i*12*sizeof(double);
        memcpy(targetAddress, temp[i].data(), 12*sizeof(double));
      }
      break;
    }
    case GP3P_RANSAC:
    {
      absRansacPtr problem;
      if(useIndices)
        problem = absRansacPtr( new absRansac( *absoluteAdapter, absRansac::GP3P, indices ) );
      else
        problem = absRansacPtr( new absRansac( *absoluteAdapter, absRansac::GP3P ) );
      opengv::sac::Ransac<absRansac> ransac;
      ransac.sac_model_ = problem;
      ransac.threshold_ = 1.0 - cos(atan(sqrt(2.0)*0.5/800.0));
      ransac.max_iterations_ = 50;
      ransac.computeModel();
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), ransac.model_coefficients_.data(), 12*sizeof(double));
      break;
    }
    case GPNP:
    {
      opengv::transformation_t temp;
      if(useIndices)
        temp = opengv::absolute_pose::gpnp(*absoluteAdapter,indices);
      else
        temp = opengv::absolute_pose::gpnp(*absoluteAdapter);
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), temp.data(), 12*sizeof(double));
      break;
    }
    case ABS_NONLIN_NONCENTRAL:
    {
      opengv::transformation_t temp;
      if(useIndices)
        temp = opengv::absolute_pose::optimize_nonlinear(*absoluteAdapter,indices);
      else
        temp = opengv::absolute_pose::optimize_nonlinear(*absoluteAdapter);
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), temp.data(), 12*sizeof(double));
      break;
    }
    case TWOPT:
    {
      opengv::translation_t temp;
      if(useIndices)
        temp = opengv::relative_pose::twopt(*relativeAdapter,false,indices);
      else
        temp = opengv::relative_pose::twopt(*relativeAdapter,false);
      int dims[2];
      dims[0] = 3;
      dims[1] = 1;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), temp.data(), 3*sizeof(double));
      break;
    }
    case TWOPT_ROTATIONONLY:
    {
      opengv::rotation_t temp;
      if(useIndices)
        temp = opengv::relative_pose::twopt_rotationOnly(*relativeAdapter,indices);
      else
        temp = opengv::relative_pose::twopt_rotationOnly(*relativeAdapter);
      int dims[2];
      dims[0] = 3;
      dims[1] = 3;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), temp.data(), 9*sizeof(double));
      break;
    }
    case ROTATIONONLY:
    {
      opengv::rotation_t temp;
      if(useIndices)
        temp = opengv::relative_pose::rotationOnly(*relativeAdapter,indices);
      else
        temp = opengv::relative_pose::rotationOnly(*relativeAdapter);
      int dims[2];
      dims[0] = 3;
      dims[1] = 3;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), temp.data(), 9*sizeof(double));
      break;
    }
    case FIVEPT_STEWENIUS:
    {
      opengv::complexEssentials_t temp2;
      if(useIndices)
        temp2 = opengv::relative_pose::fivept_stewenius(*relativeAdapter,indices);
      else
        temp2 = opengv::relative_pose::fivept_stewenius(*relativeAdapter);
      opengv::essentials_t temp;
      for(size_t i = 0; i < temp2.size(); i++)
      {
        opengv::essential_t essentialMatrix;
        for(size_t r = 0; r < 3; r++)
        {
          for(size_t c = 0; c < 3; c++)
            essentialMatrix(r,c) = temp2[i](r,c).real();
        }
        temp.push_back(essentialMatrix);
      }
      int dims[3];
      dims[0] = 3;
      dims[1] = 3;
      dims[2] = temp.size();
      plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      for( int i = 0; i < temp.size(); i++ )
      {
        void * targetAddress = ((char*) mxGetData(plhs[0])) + i*9*sizeof(double);
        memcpy(targetAddress, temp[i].data(), 9*sizeof(double));
      }
      break;
    }
    case FIVEPT_NISTER:
    {
      opengv::essentials_t temp;
      if(useIndices)
        temp = opengv::relative_pose::fivept_nister(*relativeAdapter,indices);
      else
        temp = opengv::relative_pose::fivept_nister(*relativeAdapter);
      int dims[3];
      dims[0] = 3;
      dims[1] = 3;
      dims[2] = temp.size();
      plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      for( int i = 0; i < temp.size(); i++ )
      {
        void * targetAddress = ((char*) mxGetData(plhs[0])) + i*9*sizeof(double);
        memcpy(targetAddress, temp[i].data(), 9*sizeof(double));
      }
      break;
    }
    case FIVEPT_KNEIP:
    {
      opengv::rotations_t temp;
      if(useIndices)
        temp = opengv::relative_pose::fivept_kneip(*relativeAdapter,indices);
      else
      {
        mexPrintf("opengv: Bad input to mex function opengv\n");
        mexPrintf("Assuming method: ");
        mexPrintf(methods[caseNumber]);
        mexPrintf("\n");
        mexPrintf("You must provide an indices vector\n");
        break;
      }
      int dims[3];
      dims[0] = 3;
      dims[1] = 3;
      dims[2] = temp.size();
      plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      for( int i = 0; i < temp.size(); i++ )
      {
        void * targetAddress = ((char*) mxGetData(plhs[0])) + i*9*sizeof(double);
        memcpy(targetAddress, temp[i].data(), 9*sizeof(double));
      }
      break;
    }
    case SEVENPT:
    {
      opengv::essentials_t temp;
      if(useIndices)
        temp = opengv::relative_pose::sevenpt(*relativeAdapter,indices);
      else
        temp = opengv::relative_pose::sevenpt(*relativeAdapter);
      int dims[3];
      dims[0] = 3;
      dims[1] = 3;
      dims[2] = temp.size();
      plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      for( int i = 0; i < temp.size(); i++ )
      {
        void * targetAddress = ((char*) mxGetData(plhs[0])) + i*9*sizeof(double);
        memcpy(targetAddress, temp[i].data(), 9*sizeof(double));
      }
      break;
    }
    case EIGHTPT:
    {
      opengv::essential_t temp;
      if(useIndices)
        temp = opengv::relative_pose::eightpt(*relativeAdapter,indices);
      else
        temp = opengv::relative_pose::eightpt(*relativeAdapter);
      int dims[2];
      dims[0] = 3;
      dims[1] = 3;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), temp.data(), 9*sizeof(double));
      break;
    }
    case EIGENSOLVER:
    {
      opengv::rotation_t temp;
      if(useIndices)
        temp = opengv::relative_pose::eigensolver(*relativeAdapter,indices);
      else
        temp = opengv::relative_pose::eigensolver(*relativeAdapter);
      int dims[2];
      dims[0] = 3;
      dims[1] = 3;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), temp.data(), 9*sizeof(double));
      break;
    }
    case ROTATIONONLY_RANSAC:
    {
      rotRansacPtr problem;
      if(useIndices)
        problem = rotRansacPtr( new rotRansac( *relativeAdapter, indices ) );
      else
        problem = rotRansacPtr( new rotRansac( *relativeAdapter ) );
      opengv::sac::Ransac<rotRansac> ransac;
      ransac.sac_model_ = problem;
      ransac.threshold_ = 2.0*(1.0 - cos(atan(sqrt(2.0)*0.5/800.0)));
      ransac.max_iterations_ = 50;
      ransac.computeModel();
      int dims[2];
      dims[0] = 3;
      dims[1] = 3;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), ransac.model_coefficients_.data(), 9*sizeof(double));
      break;
    }
    case FIVEPT_STEWENIUS_RANSAC:
    {
      relRansacPtr problem;
      if(useIndices)
        problem = relRansacPtr( new relRansac( *relativeAdapter, relRansac::STEWENIUS, indices ) );
      else
        problem = relRansacPtr( new relRansac( *relativeAdapter, relRansac::STEWENIUS ) );
      opengv::sac::Ransac<relRansac> ransac;
      ransac.sac_model_ = problem;
      ransac.threshold_ = 2.0*(1.0 - cos(atan(sqrt(2.0)*0.5/800.0)));
      ransac.max_iterations_ = 50;
      ransac.computeModel();
      
      opengv::transformation_t optimizedModel;
      problem->optimizeModelCoefficients(ransac.inliers_,ransac.model_coefficients_,optimizedModel);
      
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      //memcpy(mxGetData(plhs[0]), ransac.model_coefficients_.data(), 12*sizeof(double));
      memcpy(mxGetData(plhs[0]), optimizedModel.data(), 12*sizeof(double));
      break;
    }
    case FIVEPT_NISTER_RANSAC:
    {
      relRansacPtr problem;
      if(useIndices)
        problem = relRansacPtr( new relRansac( *relativeAdapter, relRansac::NISTER, indices ) );
      else
        problem = relRansacPtr( new relRansac( *relativeAdapter, relRansac::NISTER ) );
      opengv::sac::Ransac<relRansac> ransac;
      ransac.sac_model_ = problem;
      ransac.threshold_ = 2.0*(1.0 - cos(atan(sqrt(2.0)*0.5/800.0)));
      ransac.max_iterations_ = 50;
      ransac.computeModel();
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), ransac.model_coefficients_.data(), 12*sizeof(double));
      break;
    }
    case SEVENPT_RANSAC:
    {
      relRansacPtr problem;
      if(useIndices)
        problem = relRansacPtr( new relRansac( *relativeAdapter, relRansac::SEVENPT, indices ) );
      else
        problem = relRansacPtr( new relRansac( *relativeAdapter, relRansac::SEVENPT ) );
      opengv::sac::Ransac<relRansac> ransac;
      ransac.sac_model_ = problem;
      ransac.threshold_ = 2.0*(1.0 - cos(atan(sqrt(2.0)*0.5/800.0)));
      ransac.max_iterations_ = 50;
      ransac.computeModel();
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), ransac.model_coefficients_.data(), 12*sizeof(double));
      break;
    }
    case EIGHTPT_RANSAC:
    {
      relRansacPtr problem;
      if(useIndices)
        problem = relRansacPtr( new relRansac( *relativeAdapter, relRansac::EIGHTPT, indices ) );
      else
        problem = relRansacPtr( new relRansac( *relativeAdapter, relRansac::EIGHTPT ) );
      opengv::sac::Ransac<relRansac> ransac;
      ransac.sac_model_ = problem;
      ransac.threshold_ = 2.0*(1.0 - cos(atan(sqrt(2.0)*0.5/800.0)));
      ransac.max_iterations_ = 50;
      ransac.computeModel();
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), ransac.model_coefficients_.data(), 12*sizeof(double));
      break;
    }
    case EIGENSOLVER_RANSAC:
    {
      eigRansacPtr problem;
      if(useIndices)
        problem = eigRansacPtr( new eigRansac( *relativeAdapter, 10, indices ) );
      else
        problem = eigRansacPtr( new eigRansac( *relativeAdapter, 10 ) );
      opengv::sac::Ransac<eigRansac> ransac;
      ransac.sac_model_ = problem;
      ransac.threshold_ = 1.0;
      ransac.max_iterations_ = 50;
      ransac.computeModel();
      opengv::transformation_t temp;
      temp.block<3,3>(0,0) = ransac.model_coefficients_.rotation;
      temp.block<3,1>(0,3) = ransac.model_coefficients_.translation;
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), temp.data(), 12*sizeof(double));
      break;
    }
    case REL_NONLIN_CENTRAL:
    {
      opengv::transformation_t temp;
      if(useIndices)
        temp = opengv::relative_pose::optimize_nonlinear(*relativeAdapter,indices);
      else
        temp = opengv::relative_pose::optimize_nonlinear(*relativeAdapter);
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), temp.data(), 12*sizeof(double));
      break;
    }
    case SIXPT:
    {
      opengv::rotations_t temp;
      if(useIndices)
        temp = opengv::relative_pose::sixpt(*relativeAdapter,indices);
      else
        temp = opengv::relative_pose::sixpt(*relativeAdapter);	  
      int dims[3];
      dims[0] = 3;
      dims[1] = 3;
      dims[2] = temp.size();
      plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      for( int i = 0; i < temp.size(); i++ )
      {
        void * targetAddress = ((char*) mxGetData(plhs[0])) + i*9*sizeof(double);
        memcpy(targetAddress, temp[i].data(), 9*sizeof(double));
      }
      break;
    }
    case SEVENTEENPT:
    {
      opengv::transformation_t temp;
      if(useIndices)
        temp = opengv::relative_pose::seventeenpt(*relativeAdapter,indices);
      else
        temp = opengv::relative_pose::seventeenpt(*relativeAdapter);
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), temp.data(), 12*sizeof(double));
      break;
    }
    case GE:
    {
      opengv::rotation_t temp;
      opengv::geOutput_t output;
      output.rotation = relativeAdapter->getR12();
      if(useIndices)
        temp = opengv::relative_pose::ge(*relativeAdapter,indices,output);
      else
        temp = opengv::relative_pose::ge(*relativeAdapter,output);
      int dims[2];
      dims[0] = 3;
      dims[1] = 3;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]),temp.data(), 9*sizeof(double));
      break;
    }
    case SIXPT_RANSAC:
    {
      nrelRansacPtr problem;
      if(useIndices)
        problem = nrelRansacPtr( new nrelRansac( *relativeAdapter, nrelRansac::SIXPT, indices ) );
      else
        problem = nrelRansacPtr( new nrelRansac( *relativeAdapter, nrelRansac::SIXPT ) );
      opengv::sac::Ransac<nrelRansac> ransac;
      ransac.sac_model_ = problem;
      ransac.threshold_ = 2.0*(1.0 - cos(atan(sqrt(2.0)*0.5/800.0)));
      ransac.max_iterations_ = 50;
      ransac.computeModel();
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), ransac.model_coefficients_.data(), 12*sizeof(double));
      break;
    }
    case SEVENTEENPT_RANSAC:
    {
      nrelRansacPtr problem;
      if(useIndices)
        problem = nrelRansacPtr( new nrelRansac( *relativeAdapter, nrelRansac::SEVENTEENPT, indices ) );
      else
        problem = nrelRansacPtr( new nrelRansac( *relativeAdapter, nrelRansac::SEVENTEENPT ) );
      opengv::sac::Ransac<nrelRansac> ransac;
      ransac.sac_model_ = problem;
      ransac.threshold_ = 2.0*(1.0 - cos(atan(sqrt(2.0)*0.5/800.0)));
      ransac.max_iterations_ = 50;
      ransac.computeModel();
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), ransac.model_coefficients_.data(), 12*sizeof(double));
      break;
    }
    case GE_RANSAC:
    {
      nrelRansacPtr problem;
      if(useIndices)
        problem = nrelRansacPtr( new nrelRansac( *relativeAdapter, nrelRansac::GE, indices ) );
      else
        problem = nrelRansacPtr( new nrelRansac( *relativeAdapter, nrelRansac::GE ) );
      opengv::sac::Ransac<nrelRansac> ransac;
      ransac.sac_model_ = problem;
      ransac.threshold_ = 2.0*(1.0 - cos(atan(sqrt(2.0)*0.5/800.0)));
      ransac.max_iterations_ = 50;
      ransac.computeModel();
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), ransac.model_coefficients_.data(), 12*sizeof(double));
      break;
    }
    case REL_NONLIN_NONCENTRAL:
    {
      opengv::transformation_t temp;
      if(useIndices)
        temp = opengv::relative_pose::optimize_nonlinear(*relativeAdapter,indices);
      else
        temp = opengv::relative_pose::optimize_nonlinear(*relativeAdapter);
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), temp.data(), 12*sizeof(double));
      break;
    }
    case THREEPT_ARUN:
    {
      opengv::transformation_t temp;
      if(useIndices)
        temp = opengv::point_cloud::threept_arun(*pointCloudAdapter,indices);
      else
        temp = opengv::point_cloud::threept_arun(*pointCloudAdapter);
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), temp.data(), 12*sizeof(double));
      break;
    }
    case THREEPT_ARUN_RANSAC:
    {
      ptRansacPtr problem;
      if(useIndices)
        problem = ptRansacPtr( new ptRansac( *pointCloudAdapter, indices ) );
      else
        problem = ptRansacPtr( new ptRansac( *pointCloudAdapter ) );
      opengv::sac::Ransac<ptRansac> ransac;
      ransac.sac_model_ = problem;
      ransac.threshold_ = 0.1;
      ransac.max_iterations_ = 50;
      ransac.computeModel();
      int dims[2];
      dims[0] = 3;
      dims[1] = 4;
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      memcpy(mxGetData(plhs[0]), ransac.model_coefficients_.data(), 12*sizeof(double));
      break;
    }
    default: //-1
    {
      // impossible
      break;
    }
  }
}
