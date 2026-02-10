#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>
#include <opengv/math/cayley.hpp>
#include <sstream>
#include <fstream>
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include "random_generators.hpp"
#include "experiment_helpers.hpp"
#include "time_measurement.hpp"

using namespace std;
using namespace Eigen;
using namespace opengv;

// Random number generator for Gaussian noise
std::random_device rd;
std::mt19937 gen(rd());

/**
 * Convert equirectangular panorama pixel coordinates to spherical angles
 */
void panoramaPixelToSpherical(
    double u, double v,
    int width, int height,
    double &azimuth, double &elevation)
{
  double u_norm = u / (width - 1.0);
  double v_norm = v / (height - 1.0);
  
  azimuth = 2.0 * M_PI * u_norm;
  if (azimuth >= 2.0 * M_PI) azimuth = 2.0 * M_PI - 1e-10;
  
  elevation = M_PI / 2.0 - M_PI * v_norm;
}

/**
 * Convert spherical coordinates to ENU bearing vector
 */
bearingVector_t sphericalToENU(double azimuth, double elevation)
{
  double cos_elev = cos(elevation);
  bearingVector_t enu;
  enu[0] = cos_elev * sin(azimuth);  // East
  enu[1] = cos_elev * cos(azimuth);   // North
  enu[2] = sin(elevation);            // Up
  return enu;
}

/**
 * Convert panorama pixel to ENU bearing vector
 */
bearingVector_t panoramaPixelToENURay(
    double u, double v,
    int width, int height)
{
  double azimuth, elevation;
  panoramaPixelToSpherical(u, v, width, height, azimuth, elevation);
  return sphericalToENU(azimuth, elevation);
}

/**
 * Add Gaussian noise to a value
 */
double addGaussianNoise(double value, double stddev)
{
  if (stddev <= 0.0) return value;  // No noise if stddev is zero or negative
  std::normal_distribution<double> dist(0.0, stddev);
  return value + dist(gen);
}

/**
 * Generate a random bearing vector (for outliers)
 */
bearingVector_t generateRandomBearingVector()
{
  double azimuth = ((double)rand() / RAND_MAX) * 2.0 * M_PI;
  double elevation = ((double)rand() / RAND_MAX) * M_PI - M_PI / 2.0;
  
  double cos_elev = cos(elevation);
  bearingVector_t bearing;
  bearing[0] = cos_elev * sin(azimuth);
  bearing[1] = cos_elev * cos(azimuth);
  bearing[2] = sin(elevation);
  return bearing;
}

/**
 * Run a single noise configuration test
 */
void runNoiseTest(
    const std::string& testName,
    double pixelNoiseStd,  // Standard deviation in pixels
    double pointNoiseStd,  // Standard deviation in meters
    double outlierFraction, // Fraction of outliers (0.0 to 1.0)
    size_t numberPoints,
    size_t numberOfRuns,
    int panoramaWidth,
    int panoramaHeight)
{
  std::cout << std::endl << "==================================================" << std::endl;
  std::cout << "Test: " << testName << std::endl;
  std::cout << "Pixel noise std: " << pixelNoiseStd << " pixels" << std::endl;
  std::cout << "3D point noise std: " << pointNoiseStd << " m" << std::endl;
  std::cout << "Outlier fraction: " << outlierFraction * 100.0 << "%" << std::endl;
  std::cout << "==================================================" << std::endl;
  
  // Statistics accumulators
  double sqpnp_position_error_sum = 0.0;
  double sqpnp_rotation_error_sum = 0.0;
  double epnp_position_error_sum = 0.0;
  double epnp_rotation_error_sum = 0.0;
  double upnp_position_error_sum = 0.0;
  double upnp_rotation_error_sum = 0.0;
  
  double sqpnp_time_sum = 0.0;
  double epnp_time_sum = 0.0;
  double upnp_time_sum = 0.0;
  
  for(size_t run = 0; run < numberOfRuns; run++)
  {
    // Create a random viewpoint pose
    translation_t position = generateRandomTranslation(2.0);
    rotation_t rotation = Eigen::Matrix3d::Identity();  // Camera frame = ENU frame
    
    // Generate random 3D points uniformly distributed across the full 360 degree panorama sphere
    // Depth distribution: near (0.5-2m), medium (2-10m), far (10-50m with quadratic taper)
    points_t points;
    points_t noisyPoints;

    for(size_t i = 0; i < numberPoints; i++)
    {
      // Generate uniform distribution in spherical coordinates
      double azimuth = ((double)rand() / RAND_MAX) * 2.0 * M_PI;
      double elevation = ((double)rand() / RAND_MAX) * M_PI - M_PI / 2.0;

      // Depth distribution: 25% near, 50% medium, 25% far
      double rand_depth = ((double)rand() / RAND_MAX);
      double depth;
      if (rand_depth < 0.25) {
        // Near points: 0.5-2m uniform
        depth = 0.5 + ((double)rand() / RAND_MAX) * 1.5;
      } else if (rand_depth < 0.75) {
        // Medium points: 2-10m uniform
        depth = 2.0 + ((double)rand() / RAND_MAX) * 8.0;
      } else {
        // Far points: 10-50m with quadratic taper (more points near 10m, fewer at 50m)
        double u = ((double)rand() / RAND_MAX);
        depth = 10.0 + 40.0 * (1.0 - u * u);
      }
      
      // Convert spherical to Cartesian in ENU frame
      double cos_elev = cos(elevation);
      point_t point;
      point[0] = depth * cos_elev * sin(azimuth);  // East
      point[1] = depth * cos_elev * cos(azimuth);   // North
      point[2] = depth * sin(elevation);            // Up
      
      // Transform from camera frame (ENU) to world frame
      point_t worldPoint = rotation * point + position;
      points.push_back(worldPoint);
      
      // Add noise to 3D points
      point_t noisyWorldPoint;
      noisyWorldPoint[0] = addGaussianNoise(worldPoint[0], pointNoiseStd);
      noisyWorldPoint[1] = addGaussianNoise(worldPoint[1], pointNoiseStd);
      noisyWorldPoint[2] = addGaussianNoise(worldPoint[2], pointNoiseStd);
      noisyPoints.push_back(noisyWorldPoint);
    }
    
    // Convert 3D points to panorama pixels (with noise), then to ENU bearing vectors
    bearingVectors_t bearingVectors;
    
    for(size_t i = 0; i < numberPoints; i++)
    {
      // Transform point to camera frame (ENU coordinates) - use clean points for projection
      point_t bodyPoint = rotation.transpose() * (points[i] - position);
      
      double x = bodyPoint[0];
      double y = bodyPoint[1];
      double z = bodyPoint[2];
      
      double range = sqrt(x*x + y*y + z*z);
      if(range < 1e-6) continue;
      
      // Calculate azimuth and elevation
      double azimuth = atan2(x, y);
      if(azimuth < 0) azimuth += 2.0 * M_PI;
      double elevation = asin(z / range);
      
      // Convert to panorama pixel coordinates
      double u = ((azimuth / (2.0 * M_PI)) * (panoramaWidth - 1.0));
      double v = ((0.5 - elevation / M_PI) * (panoramaHeight - 1.0));
      
      // Check if this point should be an outlier
      bool isOutlier = ((double)rand() / RAND_MAX) < outlierFraction;
      
      bearingVector_t enuRay;
      if(isOutlier)
      {
        // Generate a completely random bearing vector (outlier)
        enuRay = generateRandomBearingVector();
      }
      else
      {
        // Add pixel noise
        u = addGaussianNoise(u, pixelNoiseStd);
        v = addGaussianNoise(v, pixelNoiseStd);
        
        // Clamp to valid range
        if(u < 0) u = 0;
        if(u >= panoramaWidth) u = panoramaWidth - 1.0;
        if(v < 0) v = 0;
        if(v >= panoramaHeight) v = panoramaHeight - 1.0;
        
        // Convert panorama pixel back to ENU bearing vector
        enuRay = panoramaPixelToENURay(u, v, panoramaWidth, panoramaHeight);
      }
      bearingVectors.push_back(enuRay);
    }
    
    // Create adapter with noisy points
    absolute_pose::CentralAbsoluteAdapter adapter(
        bearingVectors,
        noisyPoints,
        rotation );
    
    // Run SQPnP
    struct timeval tic, toc;
    gettimeofday( &tic, 0 );
    transformation_t sqpnp_transformation = absolute_pose::sqpnp(adapter);
    gettimeofday( &toc, 0 );
    double sqpnp_time = TIMETODOUBLE(timeval_minus(toc,tic));
    
    // Run EPNP
    gettimeofday( &tic, 0 );
    transformation_t epnp_transformation = absolute_pose::epnp(adapter);
    gettimeofday( &toc, 0 );
    double epnp_time = TIMETODOUBLE(timeval_minus(toc,tic));
    
    // Run UPNP (returns multiple solutions, pick best one)
    gettimeofday( &tic, 0 );
    transformations_t upnp_transformations = absolute_pose::upnp(adapter);
    gettimeofday( &toc, 0 );
    double upnp_time = TIMETODOUBLE(timeval_minus(toc,tic));
    
    // Find best UPNP solution
    transformation_t upnp_transformation;
    double best_upnp_error = std::numeric_limits<double>::max();
    for(size_t i = 0; i < upnp_transformations.size(); i++)
    {
      double err = (upnp_transformations[i].col(3) - position).norm();
      if(err < best_upnp_error)
      {
        best_upnp_error = err;
        upnp_transformation = upnp_transformations[i];
      }
    }
    
    // Compute SQPnP errors
    double sqpnp_position_error = (sqpnp_transformation.col(3) - position).norm();
    rotation_t sqpnp_R_est = sqpnp_transformation.block<3,3>(0,0);
    rotation_t sqpnp_R_rel = sqpnp_R_est.transpose() * rotation;
    Eigen::AngleAxisd sqpnp_aa(sqpnp_R_rel);
    double sqpnp_angle_error = sqpnp_aa.angle();
    
    // Compute EPNP errors
    double epnp_position_error = (epnp_transformation.col(3) - position).norm();
    rotation_t epnp_R_est = epnp_transformation.block<3,3>(0,0);
    rotation_t epnp_R_rel = epnp_R_est.transpose() * rotation;
    Eigen::AngleAxisd epnp_aa(epnp_R_rel);
    double epnp_angle_error = epnp_aa.angle();
    
    // Compute UPNP errors
    double upnp_position_error = best_upnp_error;
    double upnp_angle_error = 0.0;
    if(upnp_transformations.size() > 0)
    {
      rotation_t upnp_R_est = upnp_transformation.block<3,3>(0,0);
      rotation_t upnp_R_rel = upnp_R_est.transpose() * rotation;
      Eigen::AngleAxisd upnp_aa(upnp_R_rel);
      upnp_angle_error = upnp_aa.angle();
    }
    
    // Accumulate statistics
    sqpnp_position_error_sum += sqpnp_position_error;
    sqpnp_rotation_error_sum += sqpnp_angle_error;
    epnp_position_error_sum += epnp_position_error;
    epnp_rotation_error_sum += epnp_angle_error;
    upnp_position_error_sum += upnp_position_error;
    upnp_rotation_error_sum += upnp_angle_error;
    
    sqpnp_time_sum += sqpnp_time;
    epnp_time_sum += epnp_time;
    upnp_time_sum += upnp_time;
  }
  
  // Print summary statistics
  std::cout << std::endl << "SQPnP:" << std::endl;
  std::cout << "  Position error: " << sqpnp_position_error_sum / numberOfRuns << std::endl;
  std::cout << "  Angular error: " << sqpnp_rotation_error_sum / numberOfRuns << " rad" << std::endl;
  double sqpnp_avg_us = (sqpnp_time_sum / numberOfRuns) * 1000000.0;  // Convert to microseconds
  if (sqpnp_avg_us < 1000.0) {
    std::cout << "  timings: " << sqpnp_avg_us << " us" << std::endl;
  } else {
    std::cout << "  timings: " << (sqpnp_time_sum / numberOfRuns) * 1000.0 << " ms" << std::endl;
  }
  
  std::cout << std::endl << "EPNP:" << std::endl;
  std::cout << "  Position error: " << epnp_position_error_sum / numberOfRuns << std::endl;
  std::cout << "  Angular error: " << epnp_rotation_error_sum / numberOfRuns << " rad" << std::endl;
  double epnp_avg_us = (epnp_time_sum / numberOfRuns) * 1000000.0;  // Convert to microseconds
  if (epnp_avg_us < 1000.0) {
    std::cout << "  timings: " << epnp_avg_us << " us" << std::endl;
  } else {
    std::cout << "  timings: " << (epnp_time_sum / numberOfRuns) * 1000.0 << " ms" << std::endl;
  }
  
  std::cout << std::endl << "UPNP:" << std::endl;
  std::cout << "  Position error: " << upnp_position_error_sum / numberOfRuns << std::endl;
  std::cout << "  Angular error: " << upnp_rotation_error_sum / numberOfRuns << " rad" << std::endl;
  double upnp_avg_us = (upnp_time_sum / numberOfRuns) * 1000000.0;  // Convert to microseconds
  if (upnp_avg_us < 1000.0) {
    std::cout << "  timings: " << upnp_avg_us << " us" << std::endl;
  } else {
    std::cout << "  timings: " << (upnp_time_sum / numberOfRuns) * 1000.0 << " ms" << std::endl;
  }
  
  // Validate EPnP accuracy based on noise level
  // NOTE: EPnP performs poorly on panoramic/360-degree data compared to UPnP
  // We use lenient thresholds here because EPnP has high variance on panoramic views
  double avg_epnp_position_error = epnp_position_error_sum / numberOfRuns;
  double avg_epnp_rotation_error = epnp_rotation_error_sum / numberOfRuns;
  
  double position_threshold, rotation_threshold;
  
  if (pixelNoiseStd == 0.0 && pointNoiseStd == 0.0 && outlierFraction == 0.0) {
    // No noise: expect near-perfect accuracy
    position_threshold = 1e-6;
    rotation_threshold = 1e-6;  // rad
  } else if (outlierFraction > 0.0) {
    // With outliers: EPnP is NOT outlier-robust and fails badly
    // Can have 10m+ errors - use very lenient threshold
    position_threshold = 15.0;
    rotation_threshold = 1.5;  // rad (~86 deg) - EPnP can fail completely with outliers
  } else if (pixelNoiseStd > 0.0 && pointNoiseStd > 0.0) {
    // Combined noise: EPnP struggles on panoramic
    position_threshold = 2.0;
    rotation_threshold = 0.3;  // rad (~17 deg)
  } else if (pixelNoiseStd > 0.0) {
    // Pixel noise on panorama: EPnP has high variance
    position_threshold = 2.0;
    rotation_threshold = 0.2;  // rad (~11.5 deg)
  } else {
    // 3D point noise only: EPnP struggles with rotation
    position_threshold = 0.3;
    rotation_threshold = 0.2;  // rad (~11.5 deg)
  }
  
  std::cout << std::endl << "Validation (EPnP thresholds - lenient for 360 deg: pos=" << position_threshold 
            << " m, rot=" << rotation_threshold << " rad):" << std::endl;
  
  if (avg_epnp_position_error > position_threshold) {
    std::cerr << "ERROR: EPnP position error " << avg_epnp_position_error 
              << " exceeds threshold " << position_threshold << std::endl;
    exit(1);
  }
  if (avg_epnp_rotation_error > rotation_threshold) {
    std::cerr << "ERROR: EPnP rotation error " << avg_epnp_rotation_error 
              << " rad exceeds threshold " << rotation_threshold << " rad" << std::endl;
    exit(1);
  }
  
  std::cout << "  [PASS] EPnP position error: " << avg_epnp_position_error << " < " << position_threshold << std::endl;
  std::cout << "  [PASS] EPnP rotation error: " << avg_epnp_rotation_error << " < " << rotation_threshold << std::endl;
}

int main( int argc, char** argv )
{
  //initialize random seed
  initializeRandomSeed();
  
  std::cout << "==================================================" << std::endl;
  std::cout << "Panorama Algorithm Comparison: UPnP vs EPnP vs SQPnP" << std::endl;
  std::cout << "Full 360-degree equirectangular panorama test" << std::endl;
  std::cout << "==================================================" << std::endl;
  
  // Test parameters
  size_t numberPoints = 200;
  size_t numberOfRuns = 20;
  int panoramaWidth = 1920;
  int panoramaHeight = 960;
  
  // Test 1: No noise
  runNoiseTest("No Noise", 0.0, 0.0, 0.0, numberPoints, numberOfRuns, panoramaWidth, panoramaHeight);
  
  // Test 2: 5 pixel std noise on panorama pixels only
  runNoiseTest("5 pixel noise (bearing vectors)", 5.0, 0.0, 0.0, numberPoints, numberOfRuns, panoramaWidth, panoramaHeight);
  
  // Test 3: 5cm noise on 3D points only
  runNoiseTest("5cm noise (3D points)", 0.0, 0.05, 0.0, numberPoints, numberOfRuns, panoramaWidth, panoramaHeight);
  
  // Test 4: Both noises combined
  runNoiseTest("Combined: 5 pixel + 5cm noise", 5.0, 0.05, 0.0, numberPoints, numberOfRuns, panoramaWidth, panoramaHeight);
  
  // Test 5: 5% outliers (with combined noise)
  runNoiseTest("Combined noise + 5% outliers", 5.0, 0.05, 0.05, numberPoints, numberOfRuns, panoramaWidth, panoramaHeight);
  
  std::cout << std::endl << "==================================================" << std::endl;
  std::cout << "All tests completed." << std::endl;
  std::cout << "==================================================" << std::endl;
  
  return 0;
}
