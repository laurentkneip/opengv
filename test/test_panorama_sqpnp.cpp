/******************************************************************************
 * Author:   ForwardSlasher                                                  *
 * Contact:  YonatanSimson@users.noreply.github.com                          *
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

// Example: Panorama pixel to ENU ray conversion for SQPnP
// This demonstrates how to:
// 1. Convert equirectangular panorama pixels (u, v) to spherical coordinates
// 2. Convert spherical coordinates to ENU (East-North-Up) bearing vectors
// 3. Use these bearing vectors with SQPnP for pose estimation
#define NOMINMAX
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
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

/**
 * Convert equirectangular panorama pixel coordinates to spherical angles
 * @param u Pixel x-coordinate (0 to width-1)
 * @param v Pixel y-coordinate (0 to height-1)
 * @param width Panorama image width
 * @param height Panorama image height
 * @param azimuth Output azimuth angle in radians (-π to 2π, 0 = North)
 * @param elevation Output elevation angle in radians (-π/2 to π/2, 0 = horizon)
 */
void panoramaPixelToSpherical(
    double u, double v,
    int width, int height,
    double &azimuth, double &elevation)
{
  // Normalize pixel coordinates to [0, 1]
  double u_norm = u / (width - 1.0);
  double v_norm = v / (height - 1.0);

  // Convert to spherical coordinates
  // Equirectangular panorama mapping:
  // u (longitude): maps linearly from 0 to 2π
  //   u=0.5 corresponds to azimuth=0 (North ENU convention)
  //   Standard: u=0 is left edge, u=width-1 is right edge
  //   We use: u=0 → azimuth=0 (North), increasing counterclockwise (right-handed)
  azimuth = 2.0 * M_PI * (u_norm - 0.5);
  if (azimuth >= M_PI) azimuth = M_PI - 1e-10;  // Ensure [-pi, pi)
  
  // v (latitude): maps from top (north pole, elevation=π/2) to bottom (south pole, elevation=-π/2)
  //   v=0 corresponds to elevation=π/2 (looking straight up)
  //   v=height-1 corresponds to elevation=-π/2 (looking straight down)
  elevation = M_PI * (-v_norm + 0.5);  // Range: [π/2, -π/2]
}

/**
 * Convert spherical coordinates to ENU (East-North-Up) bearing vector
 * @param azimuth Azimuth angle in radians (0 = East, increases counterclockwise)
 * @param elevation Elevation angle in radians (0 = horizon, positive = up)
 * @return Normalized bearing vector in ENU format [East, North, Up]
 */
bearingVector_t sphericalToENU(double azimuth, double elevation)
{
  bearingVector_t ray;
  
  // ENU coordinate system:
  // East (x): positive to the right
  // North (y): positive forward
  // Up (z): positive upward
  
  // Convert spherical to Cartesian in ENU frame
  // x (East) = cos(elevation) * sin(azimuth)
  // y (North) = cos(elevation) * cos(azimuth)
  // z (Up) = sin(elevation)
  
  double cos_elev = cos(elevation);
  ray[0] = cos_elev * sin(azimuth);      // East
  ray[1] = cos_elev * cos(azimuth);      // North
  ray[2] = sin(elevation);                // Up
  
  // Normalize (should already be normalized, but ensure it)
  ray.normalize();
  
  return ray;
}

/**
 * Convert panorama pixel directly to ENU bearing vector
 * Convenience function combining panoramaPixelToSpherical and sphericalToENU
 */
bearingVector_t panoramaPixelToENURay(
    double u, double v,
    int width, int height)
{
  double azimuth, elevation;
  panoramaPixelToSpherical(u, v, width, height, azimuth, elevation);
  return sphericalToENU(azimuth, elevation);
}

int main( int argc, char** argv )
{
  //initialize random seed
  initializeRandomSeed();
  
  std::cout << "==================================================" << std::endl;
  std::cout << "Panorama Pixel to ENU Ray Conversion Example" << std::endl;
  std::cout << "Using SQPnP for pose estimation" << std::endl;
  std::cout << "==================================================" << std::endl;
  
  //set experiment parameters
  double noise = 0.0;
  double outlierFraction = 0.0;
  size_t numberPoints = 200;  // Increased from 100 to reduce ambiguity
  
  // Panorama image dimensions (typical 360-degree camera)
  int panoramaWidth = 1920;
  int panoramaHeight = 960;
  
  //create a random viewpoint pose
  translation_t position = generateRandomTranslation(2.0);
  
  // For panorama example: assume camera frame = ENU frame (identity rotation)
  // This allows us to demonstrate panorama->ENU conversion correctly
  // In practice, you would transform ENU rays to camera frame using camera rotation
  rotation_t rotation = Eigen::Matrix3d::Identity();  // Camera frame = ENU frame
  
  // Generate random 3D points uniformly distributed across the full 360-degree panorama sphere
  // This ensures uniform coverage in azimuth and elevation (equirectangular space)
  points_t points;
  
  double minDepth = 4.0;
  double maxDepth = 8.0;
  
  for(size_t i = 0; i < numberPoints; i++)
  {
    // Generate uniform distribution in spherical coordinates
    // Azimuth: uniform in [-π, π] (full 360-degree)
    double azimuth = ((double)rand() / RAND_MAX) * 2.0 * M_PI - M_PI;
    
    // Elevation: uniform in [-π/2, π/2] (full vertical range)
    // For equirectangular panorama, uniform elevation gives uniform pixel distribution
    double elevation = ((double)rand() / RAND_MAX) * M_PI - M_PI / 2.0;
    
    // Random depth between minDepth and maxDepth
    double depth = minDepth + ((double)rand() / RAND_MAX) * (maxDepth - minDepth);
    
    // Convert spherical to Cartesian in ENU frame
    point_t point = sphericalToENU(azimuth, elevation) * depth;
    
    // Transform from camera frame (ENU) to world frame
    // Since rotation = Identity, camera frame = world frame, but we still need to add position
    point_t worldPoint = rotation * point + position;
    points.push_back(worldPoint);
  }
  
  std::cout << std::endl << "Generated " << numberPoints 
            << " random 3D points uniformly distributed across full 360-degree panorama sphere" << std::endl;
  std::cout << "Panorama dimensions: " << panoramaWidth 
            << " x " << panoramaHeight << " (equirectangular projection)" << std::endl;
  
  // Convert 3D points to panorama pixels, then to ENU bearing vectors
  bearingVectors_t bearingVectors;
  std::vector<std::pair<double, double> > pixelCoordinates; // Store (u, v) for reference
  
  // Statistics for distribution verification
  double minAzimuth = 2.0 * M_PI, maxAzimuth = 0.0;
  double minElevation = M_PI / 2.0, maxElevation = -M_PI / 2.0;
  
  std::cout << std::endl << "Converting 3D points to panorama pixels and ENU rays..." << std::endl;
  
  for(size_t i = 0; i < numberPoints; i++)
  {
    // Transform point to camera frame (ENU coordinates)
    point_t bodyPoint = rotation.transpose() * (points[i] - position);
    
    // Convert 3D point to spherical coordinates (for panorama projection)
    // In ENU: x=East, y=North, z=Up
    double x = bodyPoint[0]; // East
    double y = bodyPoint[1]; // North
    double z = bodyPoint[2]; // Up
    
    double range = sqrt(x*x + y*y + z*z);
    if(range < 1e-6) continue; // Skip points at origin
    
    // Calculate azimuth and elevation in ENU frame
    // Azimuth: angle from North (y-axis) towards East (x-axis)
    // atan2(East, North) = atan2(x, y) gives angle from North axis
    double azimuth = atan2(x, y);
    if(azimuth < 0) azimuth += 2.0 * M_PI;  // Normalize to [0, 2π]
    
    // Elevation: angle from horizontal plane (positive = up)
    double elevation = asin(z / range);
    
    // Track distribution statistics
    if(azimuth < minAzimuth) minAzimuth = azimuth;
    if(azimuth > maxAzimuth) maxAzimuth = azimuth;
    if(elevation < minElevation) minElevation = elevation;
    if(elevation > maxElevation) maxElevation = elevation;
    
    // Convert spherical to panorama pixel coordinates
    // Equirectangular mapping:
    // u (longitude): 0 at left edge (azimuth = -π), width-1 at right edge (azimuth = π)
    // v (latitude): 0 at top (elevation = π/2), height-1 at bottom (elevation = -π/2)

    double u_norm = (azimuth / (2.0 * M_PI)) + 0.5; // Normalize to [0, 1]
    double v_norm = 0.5 - (elevation / M_PI);       // Normalize to [0, 1]
    double u = u_norm * (panoramaWidth - 1.0);
    if(u < 0) u = 0;
    if(u >= panoramaWidth) u = panoramaWidth - 1.0;

    double v = v_norm * (panoramaHeight - 1.0);
    if(v < 0) v = 0;
    if(v >= panoramaHeight) v = panoramaHeight - 1.0;
    
    pixelCoordinates.push_back(std::make_pair(u, v));
    
    // Convert panorama pixel back to ENU bearing vector
    // (This simulates what you would do with real panorama images)
    bearingVector_t enuRay = panoramaPixelToENURay(u, v, panoramaWidth, panoramaHeight);
    
    // IMPORTANT: For SQPnP to work correctly, bearing vectors must be in the camera frame.
    // The camera frame is defined by the rotation matrix. Since we're demonstrating
    // panorama->ENU conversion, we assume the camera frame IS ENU (identity rotation).
    // In practice, you would transform: cameraRay = cameraRotation.transpose() * enuRay
    
    // For this example: camera frame = ENU frame, so we can use ENU ray directly
    // But we need to verify it matches bodyPoint (which is in camera frame)
    bearingVector_t expectedRay = bodyPoint / bodyPoint.norm();
    double angleError = acos(std::max(-1.0, std::min(1.0, enuRay.dot(expectedRay))));
    
    if(i < 5) // Print first 5 for verification
    {
      std::cout << "Point " << i << ":" << std::endl;
      std::cout << "  3D point (camera frame): [" << bodyPoint[0] << ", " 
                << bodyPoint[1] << ", " << bodyPoint[2] << "]" << std::endl;
      std::cout << "  Panorama pixel: (" << (int)u << ", " << (int)v << ")" << std::endl;
      std::cout << "  Spherical: azimuth=" << azimuth*180.0/M_PI 
                << " deg, elevation=" << elevation*180.0/M_PI << " deg" << std::endl;
      std::cout << "  ENU ray: [" << enuRay[0] << ", " << enuRay[1] << ", " << enuRay[2] << "]" << std::endl;
      std::cout << "  Expected ray (camera frame): [" << expectedRay[0] << ", " 
                << expectedRay[1] << ", " << expectedRay[2] << "]" << std::endl;
      std::cout << "  Conversion angle error: " << angleError*180.0/M_PI << " deg" << std::endl;
    }
    
    // Use the ENU ray as bearing vector (assuming camera frame = ENU frame)
    // If camera frame != ENU, transform: bearingVector = cameraRotation.transpose() * enuRay
    bearingVectors.push_back(enuRay);
  }
  
  // Print distribution statistics
  std::cout << std::endl << "Point distribution statistics:" << std::endl;
  std::cout << "  Azimuth range: " << minAzimuth*180.0/M_PI << " deg to " 
            << maxAzimuth*180.0/M_PI << " deg (span: " 
            << (maxAzimuth - minAzimuth)*180.0/M_PI << " deg)" << std::endl;
  std::cout << "  Elevation range: " << minElevation*180.0/M_PI << " deg to " 
            << maxElevation*180.0/M_PI << " deg (span: " 
            << (maxElevation - minElevation)*180.0/M_PI << " deg)" << std::endl;
  
  std::cout << std::endl << "==================================================" << std::endl;
  std::cout << "Running SQPnP with ENU bearing vectors from panorama" << std::endl;
  std::cout << "==================================================" << std::endl;
  
  //print the experiment characteristics
  printExperimentCharacteristics(
      position, rotation, noise, outlierFraction );
  
  // IMPORTANT: The adapter expects bearing vectors in the camera frame.
  // Since we're simulating a panorama camera with ENU coordinate system,
  // and we've set rotation = Identity (camera frame = ENU frame),
  // the ENU rays can be used directly as bearing vectors.
  //
  // In practice with a real camera that has a different orientation:
  //   cameraRay = cameraRotation.transpose() * enuRay
  
  //create a central absolute adapter
  absolute_pose::CentralAbsoluteAdapter adapter(
      bearingVectors,
      points,
      rotation );  // rotation is Identity, so camera frame = ENU frame
  
  //run SQPnP
  std::cout << "running sqpnp with panorama-derived ENU bearing vectors" << std::endl;
  transformation_t sqpnp_transformation = absolute_pose::sqpnp(adapter);
  
  std::cout << "Ground truth transformation:" << std::endl;
  transformation_t gt_transformation;
  gt_transformation.col(3) = position;
  gt_transformation.block<3,3>(0,0) = rotation;
  std::cout << gt_transformation << std::endl << std::endl;
  
  std::cout << "SQPnP result:" << std::endl;
  std::cout << sqpnp_transformation << std::endl << std::endl;
  
  // Compute errors
  double position_error = (sqpnp_transformation.col(3) - position).norm();
  
  // Compute angular error using Rodrigues vector norm: norm(rodrigues(R_est.T @ R_GT))
  // R_rel = R_est.T * R_GT
  rotation_t R_est = sqpnp_transformation.block<3,3>(0,0);
  rotation_t R_rel = R_est.transpose() * rotation;
  // Convert to angle-axis (Rodrigues) representation
  Eigen::AngleAxisd aa(R_rel);
  double angular_error = aa.angle();  // This is the norm of the Rodrigues vector
  
  std::cout << "Position error: " << position_error << std::endl;
  std::cout << "Angular error: " << angular_error << " rad (" << angular_error * 180.0 / M_PI << " deg)" << std::endl;
  
  std::cout << std::endl << "==================================================" << std::endl;
  std::cout << "Summary: SQPnP with panorama 360-degree bearing vectors" << std::endl;
  std::cout << "Position error: " << position_error << std::endl;
  std::cout << "Angular error: " << angular_error << " rad" << std::endl;
  std::cout << "==================================================" << std::endl;
  
  return 0;
}

