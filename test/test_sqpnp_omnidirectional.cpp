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

// Test for SQPnP with omnidirectional cameras (backward-facing bearing vectors)

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>
#include <opengv/math/cayley.hpp>
#include <sstream>
#include <fstream>

#include "random_generators.hpp"
#include "experiment_helpers.hpp"
#include "time_measurement.hpp"


using namespace std;
using namespace Eigen;
using namespace opengv;

int main( int argc, char** argv )
{
  //initialize random seed
  initializeRandomSeed();
  
  //set experiment parameters
  double noise = 0.0;
  double outlierFraction = 0.0;
  size_t numberPoints = 100;

  //create a random viewpoint pose
  translation_t position = generateRandomTranslation(2.0);
  rotation_t rotation = generateRandomRotation(0.5);
  
  //create a fake central camera
  translations_t camOffsets;
  rotations_t camRotations;
  generateCentralCameraSystem( camOffsets, camRotations );
  
  //derive correspondences based on random point-cloud
  bearingVectors_t bearingVectors;
  points_t points;
  std::vector<int> camCorrespondences;
  Eigen::MatrixXd gt(3,numberPoints);
  generateRandom2D3DCorrespondences(
      position, rotation, camOffsets, camRotations, numberPoints, noise, outlierFraction,
      bearingVectors, points, camCorrespondences, gt );

  std::cout << "==================================================" << std::endl;
  std::cout << "Test 1: Forward-facing bearing vectors (standard)" << std::endl;
  std::cout << "==================================================" << std::endl;
  
  //print the experiment characteristics
  printExperimentCharacteristics(
      position, rotation, noise, outlierFraction );

  //create a central absolute adapter
  absolute_pose::CentralAbsoluteAdapter adapter(
      bearingVectors,
      points,
      rotation );

  //run SQPnP
  std::cout << "running sqpnp with forward-facing bearing vectors" << std::endl;
  transformation_t sqpnp_transformation_forward = absolute_pose::sqpnp(adapter);
  
  std::cout << "Ground truth transformation:" << std::endl;
  transformation_t gt_transformation;
  gt_transformation.col(3) = position;
  gt_transformation.block<3,3>(0,0) = rotation;
  std::cout << gt_transformation << std::endl << std::endl;
  
  std::cout << "SQPnP result:" << std::endl;
  std::cout << sqpnp_transformation_forward << std::endl << std::endl;
  
  // Compute error
  double position_error = (sqpnp_transformation_forward.col(3) - position).norm();
  std::cout << "Position error: " << position_error << std::endl;
  
  std::cout << std::endl << "==================================================" << std::endl;
  std::cout << "Test 2: Backward-facing bearing vectors (omnidirectional)" << std::endl;
  std::cout << "==================================================" << std::endl;
  
  // Create a mix of forward and backward facing bearing vectors
  bearingVectors_t mixedBearingVectors;
  points_t mixedPoints;
  
  // Use first half as forward, second half as backward
  for(size_t i = 0; i < numberPoints/2; i++)
  {
    mixedBearingVectors.push_back(bearingVectors[i]);
    mixedPoints.push_back(points[i]);
  }
  
  // For second half, flip the bearing vectors to simulate backward-facing
  for(size_t i = numberPoints/2; i < numberPoints; i++)
  {
    bearingVector_t flipped_bearing = -bearingVectors[i];
    mixedBearingVectors.push_back(flipped_bearing);
    mixedPoints.push_back(points[i]);
  }
  
  std::cout << "Using " << numberPoints/2 << " forward-facing and " 
            << numberPoints - numberPoints/2 << " backward-facing bearing vectors" << std::endl;
  
  //create adapter with mixed bearing vectors
  absolute_pose::CentralAbsoluteAdapter adapter_mixed(
      mixedBearingVectors,
      mixedPoints,
      rotation );

  //run SQPnP with mixed bearing vectors
  std::cout << "running sqpnp with mixed forward and backward-facing bearing vectors" << std::endl;
  transformation_t sqpnp_transformation_mixed = absolute_pose::sqpnp(adapter_mixed);
  
  std::cout << "SQPnP result with mixed bearings:" << std::endl;
  std::cout << sqpnp_transformation_mixed << std::endl << std::endl;
  
  // Compute error for mixed case
  double position_error_mixed = (sqpnp_transformation_mixed.col(3) - position).norm();
  std::cout << "Position error with mixed bearings: " << position_error_mixed << std::endl;
  
  std::cout << std::endl << "==================================================" << std::endl;
  std::cout << "Test 3: All backward-facing bearing vectors" << std::endl;
  std::cout << "==================================================" << std::endl;
  
  // Create all backward-facing bearing vectors
  bearingVectors_t backwardBearingVectors;
  for(size_t i = 0; i < numberPoints; i++)
  {
    bearingVector_t flipped_bearing = -bearingVectors[i];
    backwardBearingVectors.push_back(flipped_bearing);
  }
  
  std::cout << "Using all backward-facing bearing vectors" << std::endl;
  
  //create adapter with all backward bearing vectors
  absolute_pose::CentralAbsoluteAdapter adapter_backward(
      backwardBearingVectors,
      points,
      rotation );

  //run SQPnP with backward bearing vectors
  std::cout << "running sqpnp with all backward-facing bearing vectors" << std::endl;
  transformation_t sqpnp_transformation_backward = absolute_pose::sqpnp(adapter_backward);
  
  std::cout << "SQPnP result with backward bearings:" << std::endl;
  std::cout << sqpnp_transformation_backward << std::endl << std::endl;
  
  // Compute error for backward case
  double position_error_backward = (sqpnp_transformation_backward.col(3) - position).norm();
  std::cout << "Position error with backward bearings: " << position_error_backward << std::endl;
  
  std::cout << std::endl << "==================================================" << std::endl;
  std::cout << "Summary: SQPnP successfully handles omnidirectional cameras" << std::endl;
  std::cout << "Forward-facing error: " << position_error << std::endl;
  std::cout << "Mixed bearings error: " << position_error_mixed << std::endl;
  std::cout << "Backward-facing error: " << position_error_backward << std::endl;
  std::cout << "==================================================" << std::endl;
  
  return 0;
}
