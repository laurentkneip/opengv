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
  std::vector<int> camCorrespondences; //unused in the central case!
  Eigen::MatrixXd gt(3,numberPoints);
  generateRandom2D3DCorrespondences(
      position, rotation, camOffsets, camRotations, numberPoints, noise, outlierFraction,
      bearingVectors, points, camCorrespondences, gt );

  //print the experiment characteristics
  printExperimentCharacteristics(
      position, rotation, noise, outlierFraction );

  //create a central absolute adapter
  absolute_pose::CentralAbsoluteAdapter adapter(
      bearingVectors,
      points,
      rotation );

  //timer
  struct timeval tic;
  struct timeval toc;
  size_t iterations = 50;

  //run the experiments
  std::cout << "running Kneip's P2P (first two correspondences)" << std::endl;
  translation_t p2p_translation;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    p2p_translation = absolute_pose::p2p(adapter);
  gettimeofday( &toc, 0 );
  double p2p_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "running Kneip's P3P (first three correspondences)" << std::endl;
  transformations_t p3p_kneip_transformations;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    p3p_kneip_transformations = absolute_pose::p3p_kneip(adapter);
  gettimeofday( &toc, 0 );
  double p3p_kneip_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "running Gao's P3P (first three correspondences)" << std::endl;
  transformations_t p3p_gao_transformations;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    p3p_gao_transformations = absolute_pose::p3p_gao(adapter);
  gettimeofday( &toc, 0 );
  double p3p_gao_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "running epnp (all correspondences)" << std::endl;
  transformation_t epnp_transformation;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    epnp_transformation = absolute_pose::epnp(adapter);
  gettimeofday( &toc, 0 );
  double epnp_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "running epnp with 6 correspondences" << std::endl;
  std::vector<int> indices6 = getNindices(6);
  transformation_t epnp_transformation_6 =
      absolute_pose::epnp( adapter, indices6 );

  std::cout << "running upnp with all correspondences" << std::endl;
  transformations_t upnp_transformations;
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    upnp_transformations = absolute_pose::upnp(adapter);
  gettimeofday( &toc, 0 );
  double upnp_time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;

  std::cout << "running upnp with 3 correspondences" << std::endl;
  std::vector<int> indices3 = getNindices(3);
  transformations_t upnp_transformations_3 =
      absolute_pose::upnp( adapter, indices3 );

  std::cout << "setting perturbed pose";
  std::cout << "and performing nonlinear optimization" << std::endl;
  //add a small perturbation to the pose
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

  std::cout << "setting perturbed pose";
  std::cout << "and performing nonlinear optimization with 10 correspondences";
  std::cout << std::endl;
  std::vector<int> indices10 = getNindices(10);
  //add a small perturbation to the pose
  getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.1 );
  adapter.sett(t_perturbed);
  adapter.setR(R_perturbed);
  transformation_t nonlinear_transformation_10 =
      absolute_pose::optimize_nonlinear(adapter,indices10);

  //print the results
  std::cout << "results from P2P algorithm:" << std::endl;
  std::cout << p2p_translation << std::endl << std::endl;
  std::cout << "results from Kneip's P3P algorithm:" << std::endl;
  for(size_t i = 0; i < p3p_kneip_transformations.size(); i++)
    std::cout << p3p_kneip_transformations[i] << std::endl << std::endl;
  std::cout << "results from Gao's P3P algorithm:" << std::endl;
  for(size_t i = 0; i < p3p_gao_transformations.size(); i++)
    std::cout << p3p_gao_transformations[i] << std::endl << std::endl;
  std::cout << "results from epnp algorithm:" << std::endl;
  std::cout << epnp_transformation << std::endl << std::endl;
  std::cout << "results from epnp algorithm with only 6 correspondences:";
  std::cout << std::endl;
  std::cout << epnp_transformation_6 << std::endl << std::endl;
  std::cout << "results from upnp:" << std::endl;
  for(size_t i = 0; i < upnp_transformations.size(); i++)
    std::cout << upnp_transformations[i] << std::endl << std::endl;
  std::cout << "results form upnp algorithm with only 3 correspondences:";
  std::cout << std::endl;
  for(size_t i = 0; i < upnp_transformations_3.size(); i++)
    std::cout << upnp_transformations_3[i] << std::endl << std::endl;
  std::cout << "results from nonlinear algorithm:" << std::endl;
  std::cout << nonlinear_transformation << std::endl << std::endl;
  std::cout << "results from nonlinear algorithm with only 10 correspondences:";
  std::cout << std::endl;
  std::cout << nonlinear_transformation_10 << std::endl << std::endl;

  std::cout << "timings from P2P algorithm: ";
  std::cout << p2p_time << std::endl;
  std::cout << "timings from Kneip's P3P algorithm: ";
  std::cout << p3p_kneip_time << std::endl;
  std::cout << "timings from Gao's P3P algorithm: ";
  std::cout << p3p_gao_time << std::endl;
  std::cout << "timings from epnp algorithm: ";
  std::cout << epnp_time << std::endl;
  std::cout << "timings for the upnp algorithm: ";
  std::cout << upnp_time << std::endl;
  std::cout << "timings from nonlinear algorithm: ";
  std::cout << nonlinear_time << std::endl;
}
