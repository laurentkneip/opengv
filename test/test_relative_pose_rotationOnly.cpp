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

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/relative_pose/CentralRelativeAdapter.hpp>
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
  // initialize random seed
  initializeRandomSeed();

  //set experiment parameters
  double noise = 0.0;
  double outlierFraction = 0.0;
  size_t numberPoints = 10;

  //generate a random pose for viewpoint 1
  translation_t position1 = Eigen::Vector3d::Zero();
  rotation_t rotation1 = Eigen::Matrix3d::Identity();

  //generate a random pose for viewpoint 2
  translation_t position2 = Eigen::Vector3d::Zero();
  rotation_t rotation2 = generateRandomRotation(0.5);

  //create a fake central camera
  translations_t camOffsets;
  rotations_t camRotations;
  generateCentralCameraSystem( camOffsets, camRotations );

  //derive correspondences based on random point-cloud
  bearingVectors_t bearingVectors1;
  bearingVectors_t bearingVectors2;
  std::vector<int> camCorrespondences1; //unused in the central case
  std::vector<int> camCorrespondences2; //unused in the central case
  Eigen::MatrixXd gt(3,numberPoints);
  generateRandom2D2DCorrespondences(
      position1, rotation1, position2, rotation2,
      camOffsets, camRotations, numberPoints, noise, outlierFraction,
      bearingVectors1, bearingVectors2,
      camCorrespondences1, camCorrespondences2, gt );

  //Extract the relative pose
  translation_t position; rotation_t rotation;
  extractRelativePose(
      position1, position2, rotation1, rotation2, position, rotation, false );

  //print experiment characteristics
  printExperimentCharacteristics( position, rotation, noise, outlierFraction );

  //create a central relative adapter
  relative_pose::CentralRelativeAdapter adapter(
      bearingVectors1,
      bearingVectors2,
      rotation);

  //timer
  struct timeval tic;
  struct timeval toc;
  size_t iterations = 100;

  //run experiment
  std::cout << "running twopt_rotationOnly (using first 2 correpondences";
  std::cout << std::endl;
  rotation_t twopt_rotationOnly_rotation;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    twopt_rotationOnly_rotation = relative_pose::twopt_rotationOnly(adapter);
  gettimeofday( &toc, 0 );
  double twopt_rotationOnly_time =
      TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "running rotationOnly (using all correspondences)" << std::endl;
  rotation_t rotationOnly_rotation;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    rotationOnly_rotation = relative_pose::rotationOnly(adapter);
  gettimeofday( &toc, 0 );
  double rotationOnly_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "running rotationOnly with 10 correspondences" << std::endl;
  std::vector<int> indices10 = getNindices(10);
  rotation_t rotationOnly10_rotation;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    rotationOnly10_rotation = relative_pose::rotationOnly(adapter,indices10);
  gettimeofday( &toc, 0 );
  double rotationOnly10_time =
      TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  //print results
  std::cout << "results from two-points rotation algorithm:" << std::endl;
  std::cout << twopt_rotationOnly_rotation << std::endl << std::endl;
  std::cout << "results from rotation algorithm:" << std::endl;
  std::cout << rotationOnly_rotation << std::endl << std::endl;
  std::cout << "results from rotation algorithm with 10 points:" << std::endl;
  std::cout << rotationOnly10_rotation << std::endl << std::endl;

  std::cout << "timings from two-points rotation algorithm: ";
  std::cout << twopt_rotationOnly_time << std::endl;
  std::cout << "timings from rotation algorithm: ";
  std::cout << rotationOnly_time << std::endl;
  std::cout << "timings from rotation algorithm with 10 points: ";
  std::cout << rotationOnly10_time << std::endl;
}
