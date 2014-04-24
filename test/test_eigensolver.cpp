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
  double noise = 0.5;
  double outlierFraction = 0.0;
  size_t numberPoints = 10;

  //generate a random pose for viewpoint 1
  translation_t position1 = Eigen::Vector3d::Zero();
  rotation_t rotation1 = Eigen::Matrix3d::Identity();

  //generate a random pose for viewpoint 2
  translation_t position2 = generateRandomDirectionTranslation(0.0005);
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

  //create a central relative adapter
  relative_pose::CentralRelativeAdapter adapter(
      bearingVectors1,
      bearingVectors2,
      rotation);

  //run experiments
  std::cout << "running eigensolver with perturbed rotation" << std::endl;
  translation_t t_perturbed; rotation_t R_perturbed;
  getPerturbedPose( position, rotation, t_perturbed, R_perturbed, 0.01);
  adapter.setR12(R_perturbed);
  rotation_t eigensolver_rotation;
  eigensolverOutput_t output;
  output.rotation = adapter.getR12(); //transferring the initial value
  eigensolver_rotation = relative_pose::eigensolver(adapter,output);

  //print results
  std::cout << "results from eigensystem based rotation solver:" << std::endl;
  std::cout << eigensolver_rotation << std::endl << std::endl;
  std::cout << "the eigenvectors are: " << std::endl;
  std::cout << output.eigenvectors << std::endl << std::endl;
  std::cout << "the eigenvalues are: " << std::endl;
  std::cout << output.eigenvalues << std::endl << std::endl;
  std::cout << "the norms of the eigenvectors are: " << std::endl;
  std::cout << output.eigenvectors.col(0).norm() << " ";
  std::cout << output.eigenvectors.col(1).norm() << " ";
  std::cout << output.eigenvectors.col(2).norm() << std::endl << std::endl;
  std::cout << "the orthogongality is: " << std::endl;
  std::cout << output.eigenvectors.col(0).transpose()*output.eigenvectors.col(1);
  std::cout << " ";
  std::cout << output.eigenvectors.col(0).transpose()*output.eigenvectors.col(2);
  std::cout << " ";
  std::cout << output.eigenvectors.col(1).transpose()*output.eigenvectors.col(2);
  std::cout << std::endl << std::endl;
  std::cout << "the norm of the translation is: " << std::endl;
  std::cout << output.translation.norm() << std::endl;

}
