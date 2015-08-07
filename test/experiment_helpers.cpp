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
 
#include "experiment_helpers.hpp"
#include "random_generators.hpp"
#include <opengv/math/cayley.hpp>
#include <iostream>
#include <iomanip>


void
opengv::generateCentralCameraSystem(
    translations_t & camOffsets,
    rotations_t & camRotations )
{
  //create a fake camera system for a central camera
  camOffsets.push_back(Eigen::Vector3d::Zero());
  camRotations.push_back(Eigen::Matrix3d::Identity());
}

void
opengv::generateRandomCameraSystem(
    int numberCameras,
    translations_t & camOffsets,
    rotations_t & camRotations )
{
  double offset = 0.5; //this is the distance from the viewpoint origin
  
  for( int i = 0; i < numberCameras; i++ )
  {
    translation_t camOffset = generateRandomDirectionTranslation(offset);
    rotation_t camRotation = generateRandomRotation();
    camOffsets.push_back(camOffset);
    camRotations.push_back(camRotation);
  }
}

void
opengv::extractRelativePose(
    const translation_t & position1,
    const translation_t & position2,
    const rotation_t & rotation1,
    const rotation_t & rotation2,
    translation_t & relativePosition,
    rotation_t & relativeRotation,
    bool normalize )
{
  relativeRotation = rotation1.transpose() * rotation2;
  relativePosition = rotation1.transpose() * (position2 - position1);
  if(normalize)
    relativePosition = relativePosition / relativePosition.norm();
}

void
opengv::printExperimentCharacteristics(
    const translation_t & position,
    const rotation_t & rotation,
    double noise,
    double outlierFraction )
{
  std::cout << "the random position is:" << std::endl;
  std::cout << position << std::endl << std::endl;
  std::cout << "the random rotation is:" << std::endl;
  std::cout << rotation << std::endl << std::endl;
  std::cout << "the noise in the data is:" << std::endl;
  std::cout << noise << std::endl;
  std::cout << "the outlier fraction is:" << std::endl;
  std::cout << outlierFraction << std::endl;
}

void
opengv::printBearingVectorArraysMatlab(
    const bearingVectors_t & bearingVectors1,
    const bearingVectors_t & bearingVectors2 )
{
  size_t numberBearingVectors = bearingVectors1.size();
  
  //temp: print the vectors in matlab format (for debugging purposes)
  size_t precision = 10;
  std::cout << "ov1 = [";
  for( size_t i = 0; i < numberBearingVectors; i++ )
    std::cout << std::setprecision(precision) << bearingVectors1[i](0,0) << " ";
  std::cout << ";" << std::endl;
  for( size_t i = 0; i < 8; i++ )
    std::cout << std::setprecision(precision) << bearingVectors1[i](1,0) << " ";
  std::cout << ";" << std::endl;
  for( size_t i = 0; i < 8; i++ )
    std::cout << std::setprecision(precision) << bearingVectors1[i](2,0) << " ";
  std::cout << "];" << std::endl;
  std::cout << "ov2 = [";
  for( size_t i = 0; i < 8; i++ )
    std::cout << std::setprecision(precision) << bearingVectors2[i](0,0) << " ";
  std::cout << ";" << std::endl;
  for( size_t i = 0; i < 8; i++ )
    std::cout << std::setprecision(precision) << bearingVectors2[i](1,0) << " ";
  std::cout << ";" << std::endl;
  for( size_t i = 0; i < 8; i++ )
    std::cout << std::setprecision(precision) << bearingVectors2[i](2,0) << " ";
  std::cout << "];" << std::endl;
}

void
opengv::printEssentialMatrix(
    const translation_t & position,
    const rotation_t & rotation)
{
  //E transforms vectors from vp 2 to 1: x_1^T * E * x_2 = 0
  //and E = (t)_skew*R
  Eigen::Matrix3d t_skew = Eigen::Matrix3d::Zero();  
  t_skew(0,1) = -position[2];
  t_skew(0,2) = position[1];
  t_skew(1,0) = position[2];
  t_skew(1,2) = -position[0];
  t_skew(2,0) = -position[1];
  t_skew(2,1) = position[0];

  essential_t E = t_skew * rotation;
  double matrixNorm = 0.0;
  for(size_t r = 0; r < 3; r++)
  {
    for(size_t c = 0; c < 3; c++)
      matrixNorm += pow(E(r,c),2);
  }
  matrixNorm = sqrt(matrixNorm);
  for(size_t r = 0; r < 3; r++)
  {
    for(size_t c = 0; c < 3; c++)
      E(r,c) = E(r,c) / matrixNorm;
  }

  std::cout << "the random essential matrix is:" << std::endl;
  std::cout << E << std::endl;
}

void
opengv::getPerturbedPose(
    const translation_t & position,
    const rotation_t & rotation,
    translation_t & perturbedPosition,
    rotation_t & perturbedRotation,
    double amplitude )
{
  perturbedPosition = position;
  cayley_t cayley = math::rot2cayley(rotation);
  for( size_t i = 0; i < 3; i++ )
  {
    perturbedPosition[i] =
        perturbedPosition[i] +
        (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0*amplitude;
    cayley[i] =
        cayley[i] + (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0*amplitude;
  }
  perturbedRotation = math::cayley2rot(cayley);
}

std::vector<int>
opengv::getNindices( int n )
{
  std::vector<int> indices;
  for(int i = 0; i < n; i++)
    indices.push_back(i);
  return indices;
}





void
opengv::generateRandom2D3DCorrespondences(
    const translation_t & position,
    const rotation_t & rotation,
    const translations_t & camOffsets,
    const rotations_t & camRotations,
    size_t numberPoints,
    double noise,
    double outlierFraction,
    bearingVectors_t & bearingVectors,
    points_t & points,
    std::vector<int> & camCorrespondences,
    Eigen::MatrixXd & gt )
{
  //initialize point-cloud
  double minDepth = 4;
  double maxDepth = 8;
  
  for( size_t i = 0; i < (size_t) gt.cols(); i++ )
    gt.col(i) = generateRandomPoint( maxDepth, minDepth );
    
  //create the 2D3D-correspondences by looping through the cameras
  size_t numberCams = camOffsets.size();
  size_t camCorrespondence = 0;
  
  for( size_t i = 0; i < (size_t) gt.cols(); i++ )
  {
    //get the camera transformation
    translation_t camOffset = camOffsets[camCorrespondence];
    rotation_t camRotation = camRotations[camCorrespondence];

    //store the point
    points.push_back(gt.col(i));
    
    //project the point into the viewpoint frame
    point_t bodyPoint = rotation.transpose()*(gt.col(i) - position);
    
    //project the point into the camera frame
    bearingVectors.push_back(camRotation.transpose()*(bodyPoint - camOffset));

    //normalize the bearing-vector to 1
    bearingVectors[i] = bearingVectors[i] / bearingVectors[i].norm();

    //add noise
    if( noise > 0.0 )
      bearingVectors[i] = addNoise(noise,bearingVectors[i]);
    
    //push back the camera correspondence
    camCorrespondences.push_back(camCorrespondence++);
    if(camCorrespondence > (numberCams-1) )
      camCorrespondence = 0;  
  }
  
  //add outliers
  //compute the number of outliers based on fraction
  size_t numberOutliers = (size_t) floor(outlierFraction*numberPoints);
  //make the first numberOutliers points be outliers
  for(size_t i = 0; i < numberOutliers; i++)
  {
    //extract the camera transformation
    translation_t camOffset = camOffsets[camCorrespondences[i]];
    rotation_t camRotation = camRotations[camCorrespondences[i]];

    //create random point
    point_t p = generateRandomPoint(8,4);
    
    //project into viewpoint frame
    point_t bodyPoint = rotation.transpose()*(p - position);
    
    //project into camera-frame and use as outlier measurement
    bearingVectors[i] = camRotation.transpose()*(bodyPoint - camOffset);
    
    //normalize the bearing vector
    bearingVectors[i] = bearingVectors[i] / bearingVectors[i].norm();
  }
}

void
opengv::generateRandom2D2DCorrespondences(
    const translation_t & position1,
    const rotation_t & rotation1,
    const translation_t & position2,
    const rotation_t & rotation2,
    const translations_t & camOffsets,
    const rotations_t & camRotations,
    size_t numberPoints,
    double noise,
    double outlierFraction,
    bearingVectors_t & bearingVectors1,
    bearingVectors_t & bearingVectors2,
    std::vector<int> & camCorrespondences1,
    std::vector<int> & camCorrespondences2,
    Eigen::MatrixXd & gt )
{
  //initialize point-cloud
  double minDepth = 4;
  double maxDepth = 8;
  
  for( size_t i = 0; i < (size_t) gt.cols(); i++ )
    gt.col(i) = generateRandomPoint( maxDepth, minDepth );
  
  //create the 2D3D-correspondences by looping through the cameras
  size_t numberCams = camOffsets.size();
  size_t camCorrespondence = 0;
  
  for( size_t i = 0; i < (size_t) gt.cols(); i++ )
  {
    //get the camera transformation
    translation_t camOffset = camOffsets[camCorrespondence];
    rotation_t camRotation = camRotations[camCorrespondence];

    //get the point in viewpoint 1
    point_t bodyPoint1 = rotation1.transpose()*(gt.col(i) - position1);
    
    //get the point in viewpoint 2
    point_t bodyPoint2 = rotation2.transpose()*(gt.col(i) - position2);
    
    //get the point in the camera in viewpoint 1
    bearingVectors1.push_back(camRotation.transpose()*(bodyPoint1 - camOffset));
    
    //get the point in the camera in viewpoint 2
    bearingVectors2.push_back(camRotation.transpose()*(bodyPoint2 - camOffset));
    
    //normalize the bearing-vectors
    bearingVectors1[i] = bearingVectors1[i] / bearingVectors1[i].norm();
    bearingVectors2[i] = bearingVectors2[i] / bearingVectors2[i].norm();

    //add noise
    if( noise > 0.0 )
    {
      bearingVectors1[i] = addNoise(noise,bearingVectors1[i]);
      bearingVectors2[i] = addNoise(noise,bearingVectors2[i]);
    }

    //push back the camera correspondences
    camCorrespondences1.push_back(camCorrespondence);
    camCorrespondences2.push_back(camCorrespondence++);
    if(camCorrespondence > (numberCams - 1) )
      camCorrespondence = 0;
  }

  //add outliers
  size_t numberOutliers = (size_t) floor(outlierFraction*numberPoints);
  for(size_t i = 0; i < numberOutliers; i++)
  {
    //get the corresponding camera transformation
    translation_t camOffset = camOffsets[camCorrespondence];
    rotation_t camRotation = camRotations[camCorrespondence];

    //create random point
    point_t p = generateRandomPoint(8,4);
    
    //project this point into viewpoint 2
    point_t bodyPoint = rotation2.transpose()*(p - position2);
    
    //project this point into the corresponding camera in viewpoint 2
    //and use as outlier
    bearingVectors2[i] = camRotation.transpose()*(bodyPoint - camOffset);
    
    //normalize the bearing vector
    bearingVectors2[i] = bearingVectors2[i] / bearingVectors2[i].norm();
  }
}

void
opengv::generateRandom3D3DCorrespondences(
    const translation_t & position1,
    const rotation_t & rotation1,
    const translation_t & position2,
    const rotation_t & rotation2,
    size_t numberPoints,
    double noise,
    double outlierFraction,
    bearingVectors_t & points1,
    bearingVectors_t & points2,
    Eigen::MatrixXd & gt )
{
  //initialize point-cloud
  double minDepth = 4;
  double maxDepth = 8;

  for( size_t i = 0; i < (size_t) gt.cols(); i++ )
    gt.col(i) = generateRandomPoint( maxDepth, minDepth );
  
  //create the 3D-3D correspondences
  for( size_t i = 0; i < (size_t) gt.cols(); i++ )
  {
    //transform the points into the frames and store
    points1.push_back(rotation1.transpose()*(gt.col(i) - position1));
    points2.push_back(rotation2.transpose()*(gt.col(i) - position2));

    //add noise
    if( noise > 0.0 )
    {
      points1[i] = points1[i] + generateRandomTranslation(noise);
      points2[i] = points2[i] + generateRandomTranslation(noise);
    }
  }

  //add outliers
  size_t numberOutliers = (size_t) floor(outlierFraction*numberPoints);
  for(size_t i = 0; i < numberOutliers; i++)
  {
    //generate a random point
    point_t p = generateRandomPoint(8,4);
    
    //push-back in frame 2 only to create outlier
    points2[i] = rotation2.transpose()*(p - position2);
  }
}

void
opengv::generateMulti2D2DCorrespondences(
    const translation_t & position1,
    const rotation_t & rotation1,
    const translation_t & position2,
    const rotation_t & rotation2,
    const translations_t & camOffsets,
    const rotations_t & camRotations,
    size_t pointsPerCam,
    double noise,
    double outlierFraction,
    std::vector<boost::shared_ptr<bearingVectors_t> > & multiBearingVectors1,
    std::vector<boost::shared_ptr<bearingVectors_t> > & multiBearingVectors2,
    std::vector<boost::shared_ptr<Eigen::MatrixXd> > & gt )
{
  //initialize point-cloud
  double minDepth = 4;
  double maxDepth = 8;
  
  for( size_t cam = 0; cam < camOffsets.size(); cam++ )
  {
    boost::shared_ptr<Eigen::MatrixXd> gt_sub(new Eigen::MatrixXd(3,pointsPerCam));
    for( size_t i = 0; i < pointsPerCam; i++ )
      gt_sub->col(i) = generateRandomPoint( maxDepth, minDepth );
    gt.push_back(gt_sub);
  }

  //iterate through the cameras (pairs)
  for( size_t cam = 0; cam < camOffsets.size(); cam++ )
  {
    //create the bearing-vector arrays for this camera
    boost::shared_ptr<bearingVectors_t> bearingVectors1(new bearingVectors_t());
    boost::shared_ptr<bearingVectors_t> bearingVectors2(new bearingVectors_t());
    
    //get the offset and rotation of this camera
    translation_t camOffset = camOffsets[cam];
    rotation_t camRotation = camRotations[cam];

    //now iterate through the points of that camera
    for( size_t i = 0; i < (size_t) pointsPerCam; i++ )
    {
      //project a point into both viewpoint frames
      point_t bodyPoint1 = rotation1.transpose()*(gt[cam]->col(i) - position1);
      point_t bodyPoint2 = rotation2.transpose()*(gt[cam]->col(i) - position2);
      
      //project that point into the cameras
      bearingVectors1->push_back( camRotation.transpose()*(bodyPoint1 - camOffset) );
      bearingVectors2->push_back( camRotation.transpose()*(bodyPoint2 - camOffset) );
      
      //normalize the vectors
      (*bearingVectors1)[i] = (*bearingVectors1)[i] / (*bearingVectors1)[i].norm();
      (*bearingVectors2)[i] = (*bearingVectors2)[i] / (*bearingVectors2)[i].norm();

      //add noise
      if( noise > 0.0 )
      {
        (*bearingVectors1)[i] = addNoise(noise,(*bearingVectors1)[i]);
        (*bearingVectors2)[i] = addNoise(noise,(*bearingVectors2)[i]);
      }
    }

    //push back the stuff for this camera
    multiBearingVectors1.push_back(bearingVectors1);
    multiBearingVectors2.push_back(bearingVectors2);
  }

  //add outliers
  size_t outliersPerCam = (size_t) floor(outlierFraction*pointsPerCam);

  //iterate through the cameras
  for(size_t cam = 0; cam < camOffsets.size(); cam++)
  {
    //get the offset and rotation of this camera
    translation_t camOffset = camOffsets[cam];
    rotation_t camRotation = camRotations[cam];

    //add outliers
    for(size_t i = 0; i < outliersPerCam; i++)
    {
      //generate a random point
      point_t p = generateRandomPoint(8,4);
      
      //transform that point into viewpoint 2 only
      point_t bodyPoint = rotation2.transpose()*(p - position2);
      
      //use as measurement (outlier)
      (*multiBearingVectors2[cam])[i] =
          camRotation.transpose()*(bodyPoint - camOffset);
          
      //normalize
      (*multiBearingVectors2[cam])[i] =
          (*multiBearingVectors2[cam])[i] / (*multiBearingVectors2[cam])[i].norm();
    }
  }
}
