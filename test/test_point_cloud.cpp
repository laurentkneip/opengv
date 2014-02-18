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
#include <opengv/point_cloud/methods.hpp>
#include <opengv/point_cloud/PointCloudAdapter.hpp>
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
  double noise = 0.05;
  double outlierFraction = 0.0;
  size_t numberPoints = 10;

  //generate a random pose for viewpoint 1
  translation_t position1 = Eigen::Vector3d::Zero();
  rotation_t rotation1 = Eigen::Matrix3d::Identity();

  //generate a random pose for viewpoint 2
  translation_t position2 = generateRandomTranslation(2.0);
  rotation_t rotation2 = generateRandomRotation(0.5);

  //derive the correspondences based on random point-cloud
  points_t points1;
  points_t points2;
  Eigen::MatrixXd gt(3,numberPoints);
  generateRandom3D3DCorrespondences(
      position1, rotation1, position2, rotation2,
      numberPoints, noise, outlierFraction, points1, points2, gt );
    
  //Extract the relative pose
  translation_t position; rotation_t rotation;
  extractRelativePose(
      position1, position2, rotation1, rotation2, position, rotation, false );

  //print experiment characteristics
  printExperimentCharacteristics( position, rotation, noise, outlierFraction );

  //create the point-cloud adapter
  point_cloud::PointCloudAdapter adapter(
      points1, points2, position, rotation);

  //timer
  struct timeval tic;
  struct timeval toc;
  size_t iterations = 100;

  //run experiments
  std::cout << "running threept with all the points" << std::endl;
  transformation_t threept_transformation;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    threept_transformation = point_cloud::threept_arun(adapter);
  gettimeofday( &toc, 0 );
  double threept_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "running threept with three points only" << std::endl;
  std::vector<int> indices3 = getNindices(3);
  transformation_t threept_transformation_3;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    threept_transformation_3 = point_cloud::threept_arun(adapter,indices3);
  gettimeofday( &toc, 0 );
  double threept_time_3 = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "setting perturbed pose and ";
  std::cout << "performing nonlinear optimization" << std::endl;
  //add a small perturbation to the rotation
  translation_t t_perturbed; rotation_t R_perturbed;
  getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.1);
  adapter.sett12(t_perturbed);
  adapter.setR12(R_perturbed);
  transformation_t nonlinear_transformation;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    nonlinear_transformation = point_cloud::optimize_nonlinear(adapter);
  gettimeofday( &toc, 0 );
  double nonlinear_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "setting perturbed pose and ";
  std::cout << "performing nonlinear optimization with three points only";
  std::cout << std::endl;
  //add a small perturbation to the rotation
  getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.01);
  adapter.sett12(t_perturbed);
  adapter.setR12(R_perturbed);
  transformation_t nonlinear_transformation_3 =
      point_cloud::optimize_nonlinear(adapter,indices3);

  //print results
  std::cout << "results from threept_arun algorithm:" << std::endl;
  std::cout << threept_transformation << std::endl << std::endl;
  std::cout << "results from threept_arun with 3 points only:" << std::endl;
  std::cout << threept_transformation_3 << std::endl << std::endl;
  std::cout << "results of nonlinear optimization:" << std::endl;
  std::cout << nonlinear_transformation << std::endl << std::endl;
  std::cout << "results of nonlinear optimization with 3 points only:" << std::endl;
  std::cout << nonlinear_transformation_3 << std::endl << std::endl;

  std::cout << "timings from threept_arun algorithm: ";
  std::cout << threept_time << std::endl;
  std::cout << "timings from threept_arun algorithm with 3 points only: ";
  std::cout << threept_time_3 << std::endl;
  std::cout << "timing of nonlinear optimization: ";
  std::cout << nonlinear_time << std::endl;
}
