#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>
#include <opengv/types.hpp>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;
using namespace opengv;

double computePositionError(const transformation_t& T, const translation_t& gt_position) {
  return (T.col(3) - gt_position).norm();
}

double computeRotationError(const rotation_t& R, const rotation_t& gt_rotation) {
  rotation_t R_rel = R.transpose() * gt_rotation;
  double trace = R_rel.trace();
  // Clamp trace to valid range for acos
  trace = std::min(3.0, std::max(-1.0, trace));
  return acos((trace - 1.0) / 2.0);
}

int main(int argc, char** argv) {
  cout << "================================================================================" << endl;
  cout << "Testing P3P Minimal Solvers with Backward-Facing Bearing Vectors" << endl;
  cout << "================================================================================" << endl;

  // Ground truth camera pose (Identity rotation, random position)
  translation_t gt_position;
  gt_position << 1.0, 0.5, 0.3;
  rotation_t gt_rotation = rotation_t::Identity();

  cout << "\nGround truth:" << endl;
  cout << "  Position: " << gt_position.transpose() << endl;
  cout << "  Rotation: Identity" << endl;

  // Create 3D points: 2 forward, 1 backward
  bearingVectors_t bearing_vectors;
  points_t points;

  // Point 0: Forward (z > 0 in camera frame)
  point_t p0;
  p0 << 2.0, 1.0, 5.0;
  points.push_back(p0);

  // Point 1: Forward (z > 0 in camera frame)
  point_t p1;
  p1 << 3.0, -1.0, 4.0;
  points.push_back(p1);

  // Point 2: Backward (z < 0 in camera frame)
  point_t p2;
  p2 << 1.5, 0.0, -3.0;
  points.push_back(p2);

  // Compute bearing vectors in camera frame
  cout << "\n3D Points (world frame) and bearing vectors:" << endl;
  for(size_t i = 0; i < 3; i++) {
    point_t p_cam = gt_rotation.transpose() * (points[i] - gt_position);
    bearingVector_t bv = p_cam.normalized();
    bearing_vectors.push_back(bv);

    cout << "  Point " << i << ": " << points[i].transpose()
         << " -> camera: " << p_cam.transpose()
         << " (z=" << (p_cam[2] > 0 ? "forward" : "BACKWARD") << ")" << endl;
    cout << "           bearing: " << bv.transpose()
         << " (z=" << fixed << setprecision(3) << bv[2] << ")" << endl;
  }

  // Create adapter
  absolute_pose::CentralAbsoluteAdapter adapter(
      bearing_vectors,
      points);

  // Test P3P Kneip
  cout << "\n" << string(80, '-') << endl;
  cout << "Testing P3P Kneip:" << endl;
  cout << string(80, '-') << endl;

  transformations_t kneip_solutions = absolute_pose::p3p_kneip(adapter);
  cout << "Found " << kneip_solutions.size() << " solutions" << endl;

  bool kneip_found_correct = false;
  for(size_t i = 0; i < kneip_solutions.size(); i++) {
    double pos_err = computePositionError(kneip_solutions[i], gt_position);
    double rot_err = computeRotationError(kneip_solutions[i].block<3,3>(0,0), gt_rotation);

    cout << "  Solution " << i << ": "
         << "pos_err=" << scientific << setprecision(2) << pos_err << " m, "
         << "rot_err=" << scientific << setprecision(2) << rot_err << " rad ("
         << fixed << setprecision(3) << (rot_err * 180.0 / M_PI) << " deg)" << endl;

    if(pos_err < 1e-6 && rot_err < 1e-6) {
      cout << "    -> CORRECT SOLUTION" << endl;
      kneip_found_correct = true;
    }
  }

  if(!kneip_found_correct) {
    cout << "\n[FAIL] P3P Kneip did not find correct solution for backward-facing vectors!" << endl;
    return 1;
  }
  cout << "[PASS] P3P Kneip correctly handles backward-facing vectors" << endl;

  // Test P3P Gao
  cout << "\n" << string(80, '-') << endl;
  cout << "Testing P3P Gao:" << endl;
  cout << string(80, '-') << endl;

  transformations_t gao_solutions = absolute_pose::p3p_gao(adapter);
  cout << "Found " << gao_solutions.size() << " solutions" << endl;

  bool gao_found_correct = false;
  for(size_t i = 0; i < gao_solutions.size(); i++) {
    double pos_err = computePositionError(gao_solutions[i], gt_position);
    double rot_err = computeRotationError(gao_solutions[i].block<3,3>(0,0), gt_rotation);

    cout << "  Solution " << i << ": "
         << "pos_err=" << scientific << setprecision(2) << pos_err << " m, "
         << "rot_err=" << scientific << setprecision(2) << rot_err << " rad ("
         << fixed << setprecision(3) << (rot_err * 180.0 / M_PI) << " deg)" << endl;

    if(pos_err < 1e-6 && rot_err < 1e-6) {
      cout << "    -> CORRECT SOLUTION" << endl;
      gao_found_correct = true;
    }
  }

  if(!gao_found_correct) {
    cout << "\n[FAIL] P3P Gao did not find correct solution for backward-facing vectors!" << endl;
    return 1;
  }
  cout << "[PASS] P3P Gao correctly handles backward-facing vectors" << endl;

  cout << "\n" << string(80, '=') << endl;
  cout << "Conclusion:" << endl;
  cout << string(80, '=') << endl;
  cout << "Both P3P Kneip and P3P Gao can handle bearing vectors in all hemispheres!" << endl;
  cout << "These minimal solvers are suitable for omnidirectional/panoramic cameras." << endl;
  cout << "\n[SUCCESS] All tests passed!" << endl;

  return 0;
}
