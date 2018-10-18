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
#include <limits.h>
#include <Eigen/Eigen>
#include <opengv/point_cloud/methods.hpp>
#include <opengv/point_cloud/PointCloudAdapter.hpp>
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac/Lmeds.hpp>
#include <opengv/sac_problems/point_cloud/PointCloudSacProblem.hpp>
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
  double outlierFraction = 0.1;
  size_t numberPoints = 100;

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

  //Create a PointCloudSacProblem and Ransac
  sac::Ransac<
      sac_problems::point_cloud::PointCloudSacProblem> ransac;
  std::shared_ptr<
      sac_problems::point_cloud::PointCloudSacProblem> relposeproblem_ptr(
      new sac_problems::point_cloud::PointCloudSacProblem(adapter));
  ransac.sac_model_ = relposeproblem_ptr;
  ransac.threshold_ = 0.1;
  ransac.max_iterations_ = 50;

  //Run the experiment
  struct timeval tic;
  struct timeval toc;
  gettimeofday( &tic, 0 );
  ransac.computeModel(0);
  gettimeofday( &toc, 0 );
  double ransac_time = TIMETODOUBLE(timeval_minus(toc,tic));

  //print the results
  std::cout << "the ransac threshold is: " << ransac.threshold_ << std::endl;
  std::cout << "the ransac results is: " << std::endl;
  std::cout << ransac.model_coefficients_ << std::endl << std::endl;
  std::cout << "Ransac needed " << ransac.iterations_ << " iterations and ";
  std::cout << ransac_time << " seconds" << std::endl << std::endl;
  std::cout << "the number of inliers is: " << ransac.inliers_.size();
  std::cout << std::endl << std::endl;
  std::cout << "the found inliers are: " << std::endl;
  for(size_t i = 0; i < ransac.inliers_.size(); i++)
    std::cout << ransac.inliers_[i] << " ";
  std::cout << std::endl << std::endl;

  // Create Lmeds
  sac::Lmeds<
      sac_problems::point_cloud::PointCloudSacProblem> lmeds;
  lmeds.sac_model_ = relposeproblem_ptr;
  lmeds.threshold_ = 0.1;
  lmeds.max_iterations_ = 50;

  //Run the experiment
  gettimeofday( &tic, 0 );
  lmeds.computeModel(0);
  gettimeofday( &toc, 0 );
  double lmeds_time = TIMETODOUBLE(timeval_minus(toc,tic));

  //print the results
  std::cout << "the lmeds threshold is: " << lmeds.threshold_ << std::endl;
  std::cout << "the lmeds results is: " << std::endl;
  std::cout << lmeds.model_coefficients_ << std::endl << std::endl;
  std::cout << "Lmeds needed " << lmeds.iterations_ << " iterations and ";
  std::cout << lmeds_time << " seconds" << std::endl << std::endl;
  std::cout << "the number of inliers is: " << lmeds.inliers_.size();
  std::cout << std::endl << std::endl;
  std::cout << "the found inliers are: " << std::endl;
  for(size_t i = 0; i < lmeds.inliers_.size(); i++)
    std::cout << lmeds.inliers_[i] << " ";
  std::cout << std::endl << std::endl;
}
