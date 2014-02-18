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
  translation_t position2 = generateRandomTranslation(2.0);
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
      position1, position2, rotation1, rotation2, position, rotation );

  //print experiment characteristics
  printExperimentCharacteristics( position, rotation, noise, outlierFraction );
  
  //compute and print the essential-matrix
  printEssentialMatrix( position, rotation );

  //create a central relative adapter
  relative_pose::CentralRelativeAdapter adapter(
      bearingVectors1,
      bearingVectors2,
      rotation);

  //timer
  struct timeval tic;
  struct timeval toc;
  size_t iterations = 50;

  //running experiments
  std::cout << "running twopt" << std::endl;
  translation_t twopt_translation;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    twopt_translation = relative_pose::twopt(adapter,true);
  gettimeofday( &toc, 0 );
  double twopt_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "running fivept_stewenius" << std::endl;
  complexEssentials_t fivept_stewenius_essentials;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    fivept_stewenius_essentials = relative_pose::fivept_stewenius(adapter);
  gettimeofday( &toc, 0 );
  double fivept_stewenius_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "running fivept_nister" << std::endl;
  essentials_t fivept_nister_essentials;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    fivept_nister_essentials = relative_pose::fivept_nister(adapter);
  gettimeofday( &toc, 0 );
  double fivept_nister_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "running fivept_kneip" << std::endl;
  rotations_t fivept_kneip_rotations;
  gettimeofday( &tic, 0 );
  std::vector<int> indices5 = getNindices(5);
  for(size_t i = 0; i < iterations; i++)
    fivept_kneip_rotations = relative_pose::fivept_kneip(adapter,indices5);
  gettimeofday( &toc, 0 );
  double fivept_kneip_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "running sevenpt" << std::endl;
  essentials_t sevenpt_essentials;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    sevenpt_essentials = relative_pose::sevenpt(adapter);
  gettimeofday( &toc, 0 );
  double sevenpt_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "running eightpt" << std::endl;
  essential_t eightpt_essential;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    eightpt_essential = relative_pose::eightpt(adapter);
  gettimeofday( &toc, 0 );
  double eightpt_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "setting perturbed rotation and ";
  std::cout << "running eigensolver" << std::endl;
  translation_t t_perturbed; rotation_t R_perturbed;
  getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.01);
  rotation_t eigensolver_rotation;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
  {
    adapter.setR12(R_perturbed);
    eigensolver_rotation = relative_pose::eigensolver(adapter);
  }
  gettimeofday( &toc, 0 );
  double eigensolver_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "setting perturbed pose and ";
  std::cout << "performing nonlinear optimization" << std::endl;
  getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.1);
  transformation_t nonlinear_transformation;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
  {
    adapter.sett12(t_perturbed);
    adapter.setR12(R_perturbed);
    nonlinear_transformation = relative_pose::optimize_nonlinear(adapter);
  }
  gettimeofday( &toc, 0 );
  double nonlinear_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "setting perturbed pose and ";
  std::cout << "performing nonlinear optimization with 10 indices" << std::endl;
  std::vector<int> indices10 = getNindices(10);
  getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.1);
  adapter.sett12(t_perturbed);
  adapter.setR12(R_perturbed);
  transformation_t nonlinear_transformation_10 =
      relative_pose::optimize_nonlinear(adapter,indices10);

  //print results
  std::cout << "results from two-points algorithm:" << std::endl;
  std::cout << twopt_translation << std::endl << std::endl;
  std::cout << "results from stewenius' five-point algorithm:" << std::endl;
  for( size_t i = 0; i < fivept_stewenius_essentials.size(); i++ )
    std::cout << fivept_stewenius_essentials.at(i) << std::endl << std::endl;
  std::cout << "results from nisters' five-point algorithm:" << std::endl;
  for( size_t i = 0; i < fivept_nister_essentials.size(); i++ )
    std::cout << fivept_nister_essentials.at(i) << std::endl << std::endl;
  std::cout << "results from kneip's five-point algorithm:" << std::endl;
  for( size_t i = 0; i < fivept_kneip_rotations.size(); i++ )
    std::cout << fivept_kneip_rotations.at(i) << std::endl << std::endl;
  std::cout << "results from seven-point algorithm:" << std::endl;
  for( size_t i = 0; i < sevenpt_essentials.size(); i++ )
    std::cout << sevenpt_essentials.at(i) << std::endl << std::endl;
  std::cout << "results from eight-point algorithm:" << std::endl;
  std::cout << eightpt_essential << std::endl << std::endl;
  std::cout << "results from eigensystem based rotation solver:" << std::endl;
  std::cout << eigensolver_rotation << std::endl << std::endl << std::endl;
  std::cout << "results from nonlinear algorithm:" << std::endl;
  std::cout << nonlinear_transformation << std::endl << std::endl;
  std::cout << "results from nonlinear algorithm with only few correspondences:";
  std::cout << std::endl;
  std::cout << nonlinear_transformation_10 << std::endl << std::endl;

  std::cout << "timings from two-points algorithm: ";
  std::cout << twopt_time << std::endl;
  std::cout << "timings from stewenius' five-point algorithm: ";
  std::cout << fivept_stewenius_time << std::endl;
  std::cout << "timings from nisters' five-point algorithm: ";
  std::cout << fivept_nister_time << std::endl;
  std::cout << "timings from kneip's five-point algorithm: ";
  std::cout << fivept_kneip_time << std::endl;
  std::cout << "timings from seven-point algorithm: ";
  std::cout << sevenpt_time << std::endl;
  std::cout << "timings from eight-point algorithm: ";
  std::cout << eightpt_time << std::endl;
  std::cout << "timings from eigensystem based rotation solver: ";
  std::cout << eigensolver_time << std::endl;
  std::cout << "timings from nonlinear algorithm: ";
  std::cout << nonlinear_time << std::endl;
}
