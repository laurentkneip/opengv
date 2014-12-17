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
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/absolute_pose/NoncentralAbsoluteAdapter.hpp>
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
  int numberCameras = 4;

  //create a random viewpoint pose
  translation_t position = generateRandomTranslation(2.0);
  rotation_t rotation = generateRandomRotation(0.5);
  
  //create a random camera-system
  translations_t camOffsets;
  rotations_t camRotations;
  generateRandomCameraSystem( numberCameras, camOffsets, camRotations );
  
  //derive correspondences based on random point-cloud
  bearingVectors_t bearingVectors;
  points_t points;
  std::vector<int> camCorrespondences;
  Eigen::MatrixXd gt(3,numberPoints);
  generateRandom2D3DCorrespondences(
      position, rotation, camOffsets, camRotations, numberPoints, noise, outlierFraction,
      bearingVectors, points, camCorrespondences, gt );

  //print the experiment characteristics
  printExperimentCharacteristics(
      position, rotation, noise, outlierFraction );

  //create a non-central absolute adapter
  absolute_pose::NoncentralAbsoluteAdapter adapter(
      bearingVectors,
      camCorrespondences,
      points,
      camOffsets,
      camRotations,
      rotation);

  //timer
  struct timeval tic;
  struct timeval toc;
  size_t iterations = 50;

  //run the experiments
  std::cout << "running Kneip's GP3P (using first three correspondences/";
  std::cout << std::endl;
  transformations_t gp3p_transformations;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    gp3p_transformations = absolute_pose::gp3p(adapter);
  gettimeofday( &toc, 0 );
  double gp3p_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "running gpnp over all correspondences" << std::endl;
  transformation_t gpnp_transformation;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    gpnp_transformation = absolute_pose::gpnp(adapter);
  gettimeofday( &toc, 0 );
  double gpnp_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "running gpnp over 6 correspondences" << std::endl;
  std::vector<int> indices6 = getNindices(6);
  transformation_t gpnp_transformation_6 =
      absolute_pose::gpnp(adapter,indices6);
  
  std::cout << "running upnp over all correspondences" << std::endl;
  transformations_t upnp_transformations;
  gettimeofday( &tic, 0 );
  for( size_t i = 0; i < iterations; i++ )
    upnp_transformations = absolute_pose::upnp(adapter);
  gettimeofday( &toc, 0);
  double upnp_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
  
  //std::cout << "running upnp over 6 correspondences" << std::endl;
  //transformations_t upnp_transformations_6 =
  //    absolute_pose::upnp(adapter,indices6);

  std::cout << "setting perturbed pose and ";
  std::cout << "performing nonlinear optimization" << std::endl;
  //add a small perturbation to the rotation
  translation_t t_perturbed; rotation_t R_perturbed;
  getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.1 );
  transformation_t nonlinear_transformation;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
  {
    adapter.sett(t_perturbed);
    adapter.setR(R_perturbed);
    nonlinear_transformation = absolute_pose::optimize_nonlinear(adapter);
  }
  gettimeofday( &toc, 0 );
  double nonlinear_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "setting perturbed pose and ";
  std::cout << "performing nonlinear optimization with 10 correspondences";
  std::cout << std::endl;
  std::vector<int> indices10 = getNindices(10);
  //add a small perturbation to the rotation
  getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.1 );
  adapter.sett(t_perturbed);
  adapter.setR(R_perturbed);
  transformation_t nonlinear_transformation_10 =
      absolute_pose::optimize_nonlinear(adapter,indices10);

  //print the results
  std::cout << "results from gp3p algorithm:" << std::endl;
  for(size_t i = 0; i < gp3p_transformations.size(); i++)
    std::cout << gp3p_transformations[i] << std::endl << std::endl;
  std::cout << "results from gpnp algorithm:" << std::endl;
  std::cout << gpnp_transformation << std::endl << std::endl;
  std::cout << "results from gpnp algorithm with only 6 correspondences:";
  std::cout << std::endl;
  std::cout << gpnp_transformation_6 << std::endl << std::endl;
  std::cout << "results from upnp algorithm:" << std::endl;
  for(size_t i = 0; i < upnp_transformations.size(); i++)
    std::cout << upnp_transformations[i] << std::endl << std::endl;
  std::cout << "results from nonlinear algorithm:" << std::endl;
  std::cout << nonlinear_transformation << std::endl << std::endl;
  std::cout << "results from nonlinear algorithm with only 10 correspondences:";
  std::cout << std::endl;
  std::cout << nonlinear_transformation_10 << std::endl << std::endl;

  std::cout << "timings from gp3p algorithm: ";
  std::cout << gp3p_time << std::endl;
  std::cout << "timings from gpnp algorithm: ";
  std::cout << gpnp_time << std::endl;
  std::cout << "timings from upnp algorithm: ";
  std::cout << upnp_time << std::endl;
  std::cout << "timings from nonlinear algorithm: ";
  std::cout << nonlinear_time << std::endl;
}
