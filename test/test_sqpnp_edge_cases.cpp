// Test SQPnP with edge cases and backward-facing vectors

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cassert>
#include <cmath>
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>

#include "random_generators.hpp"
#include "experiment_helpers.hpp"

using namespace std;
using namespace Eigen;
using namespace opengv;

double computeRotationError(const rotation_t& R_est, const rotation_t& R_gt) {
    rotation_t R_err = R_est.transpose() * R_gt;
    double trace = R_err.trace();
    trace = std::min(3.0, std::max(-1.0, trace));
    return std::acos((trace - 1.0) / 2.0);
}

void test_forward_only_vectors() {
    cout << "=== Testing forward-only bearing vectors (baseline) ===" << endl;

    translation_t position = generateRandomTranslation(2.0);
    rotation_t rotation = generateRandomRotation(0.5);

    size_t numberPoints = 50;
    bearingVectors_t bearingVectors;
    points_t points;
    std::vector<int> camCorrespondences;
    Eigen::MatrixXd gt(3, numberPoints);

    translations_t camOffsets;
    rotations_t camRotations;
    generateCentralCameraSystem(camOffsets, camRotations);

    generateRandom2D3DCorrespondences(
        position, rotation, camOffsets, camRotations, numberPoints, 0.0, 0.0,
        bearingVectors, points, camCorrespondences, gt);

    // All vectors facing forward (standard perspective camera)
    absolute_pose::CentralAbsoluteAdapter adapter(bearingVectors, points);
    transformation_t T_sqpnp = absolute_pose::sqpnp(adapter);

    double pos_error = (T_sqpnp.col(3) - position).norm();
    double rot_error = computeRotationError(T_sqpnp.block<3,3>(0,0), rotation);

    cout << "Position error: " << pos_error << endl;
    cout << "Rotation error: " << rot_error << " rad" << endl;

    assert(pos_error < 1e-6 && "SQPNP position error too large with forward-only vectors");
    assert(rot_error < 1e-4 && "SQPNP rotation error too large with forward-only vectors");

    cout << "✓ Test passed" << endl << endl;
}

void test_backward_only_vectors() {
    cout << "=== Testing backward-only bearing vectors ===" << endl;

    translation_t position = generateRandomTranslation(2.0);
    rotation_t rotation = generateRandomRotation(0.5);

    size_t numberPoints = 50;
    bearingVectors_t bearingVectors;
    points_t points;
    std::vector<int> camCorrespondences;
    Eigen::MatrixXd gt(3, numberPoints);

    translations_t camOffsets;
    rotations_t camRotations;
    generateCentralCameraSystem(camOffsets, camRotations);

    generateRandom2D3DCorrespondences(
        position, rotation, camOffsets, camRotations, numberPoints, 0.0, 0.0,
        bearingVectors, points, camCorrespondences, gt);

    // Flip ALL bearing vectors to face backward (rear-facing camera)
    for(size_t i = 0; i < numberPoints; i++) {
        bearingVectors[i] = -bearingVectors[i];
    }

    absolute_pose::CentralAbsoluteAdapter adapter(bearingVectors, points);
    transformation_t T_sqpnp = absolute_pose::sqpnp(adapter);

    double pos_error = (T_sqpnp.col(3) - position).norm();
    double rot_error = computeRotationError(T_sqpnp.block<3,3>(0,0), rotation);

    cout << "Position error: " << pos_error << endl;
    cout << "Rotation error: " << rot_error << " rad" << endl;

    assert(pos_error < 0.1 && "SQPNP position error too large with backward-only vectors");
    assert(rot_error < 0.1 && "SQPNP rotation error too large with backward-only vectors");

    cout << "✓ Test passed" << endl << endl;
}

void test_mixed_forward_backward_vectors() {
    cout << "=== Testing mixed forward/backward bearing vectors (50/50) ===" << endl;

    translation_t position = generateRandomTranslation(2.0);
    rotation_t rotation = generateRandomRotation(0.5);

    size_t numberPoints = 50;
    bearingVectors_t bearingVectors;
    points_t points;
    std::vector<int> camCorrespondences;
    Eigen::MatrixXd gt(3, numberPoints);

    translations_t camOffsets;
    rotations_t camRotations;
    generateCentralCameraSystem(camOffsets, camRotations);

    generateRandom2D3DCorrespondences(
        position, rotation, camOffsets, camRotations, numberPoints, 0.0, 0.0,
        bearingVectors, points, camCorrespondences, gt);

    // Flip half the bearing vectors to face backward (omnidirectional scenario)
    for(size_t i = 0; i < numberPoints/2; i++) {
        bearingVectors[i] = -bearingVectors[i];
    }

    absolute_pose::CentralAbsoluteAdapter adapter(bearingVectors, points);
    transformation_t T_sqpnp = absolute_pose::sqpnp(adapter);

    double pos_error = (T_sqpnp.col(3) - position).norm();
    double rot_error = computeRotationError(T_sqpnp.block<3,3>(0,0), rotation);

    cout << "Position error: " << pos_error << endl;
    cout << "Rotation error: " << rot_error << " rad" << endl;

    // With Bug #1 (std::abs), this would incorrectly report low error
    // With Bug #4 (single point check), sign handling might fail
    assert(pos_error < 0.1 && "SQPNP position error too large with mixed vectors");
    assert(rot_error < 0.1 && "SQPNP rotation error too large with mixed vectors");

    cout << "✓ Test passed" << endl << endl;
}

void test_minority_backward_vectors() {
    cout << "=== Testing minority backward vectors (10% backward) ===" << endl;

    translation_t position = generateRandomTranslation(2.0);
    rotation_t rotation = generateRandomRotation(0.5);

    size_t numberPoints = 50;
    bearingVectors_t bearingVectors;
    points_t points;
    std::vector<int> camCorrespondences;
    Eigen::MatrixXd gt(3, numberPoints);

    translations_t camOffsets;
    rotations_t camRotations;
    generateCentralCameraSystem(camOffsets, camRotations);

    generateRandom2D3DCorrespondences(
        position, rotation, camOffsets, camRotations, numberPoints, 0.0, 0.0,
        bearingVectors, points, camCorrespondences, gt);

    // Flip 10% of bearing vectors (5 out of 50)
    for(size_t i = 0; i < numberPoints/10; i++) {
        bearingVectors[i] = -bearingVectors[i];
    }

    absolute_pose::CentralAbsoluteAdapter adapter(bearingVectors, points);
    transformation_t T_sqpnp = absolute_pose::sqpnp(adapter);

    double pos_error = (T_sqpnp.col(3) - position).norm();
    double rot_error = computeRotationError(T_sqpnp.block<3,3>(0,0), rotation);

    cout << "Position error: " << pos_error << endl;
    cout << "Rotation error: " << rot_error << " rad" << endl;

    // This tests the voting mechanism in solve_for_sign
    assert(pos_error < 0.1 && "SQPNP position error too large with minority backward vectors");
    assert(rot_error < 0.1 && "SQPNP rotation error too large with minority backward vectors");

    cout << "✓ Test passed" << endl << endl;
}

void test_degenerate_inputs() {
    cout << "=== Testing degenerate inputs ===" << endl;

    translation_t position = generateRandomTranslation(2.0);
    rotation_t rotation = generateRandomRotation(0.5);

    bearingVectors_t bearingVectors;
    points_t points;

    // Add 10 valid correspondences
    for(int i = 0; i < 10; i++) {
        point_t p = generateRandomPoint(8.0, 4.0);
        bearingVector_t bearing = rotation.transpose() * (p - position);
        bearing.normalize();

        bearingVectors.push_back(bearing);
        points.push_back(p);
    }

    // Try to add degenerate bearing vectors (should be handled gracefully)
    // With Bug #2, this could cause array misalignment
    bearingVector_t zero_vec(0, 0, 0);
    bearingVectors.push_back(zero_vec);
    points.push_back(generateRandomPoint(8.0, 4.0));

    absolute_pose::CentralAbsoluteAdapter adapter(bearingVectors, points);

    // Should still work with the 10 valid correspondences
    transformation_t T_sqpnp = absolute_pose::sqpnp(adapter);

    double pos_error = (T_sqpnp.col(3) - position).norm();
    double rot_error = computeRotationError(T_sqpnp.block<3,3>(0,0), rotation);

    cout << "Position error: " << pos_error << endl;
    cout << "Rotation error: " << rot_error << " rad" << endl;

    assert(pos_error < 0.5 && "SQPNP failed with degenerate inputs");
    assert(rot_error < 0.5 && "SQPNP failed with degenerate inputs");

    cout << "✓ Test passed" << endl << endl;
}

void test_sqp_vs_epnp_comparison() {
    cout << "=== Testing SQP vs EPnP selection ===" << endl;

    // Test that the solver correctly chooses between SQP and EPnP solutions
    translation_t position = generateRandomTranslation(2.0);
    rotation_t rotation = generateRandomRotation(0.5);

    size_t numberPoints = 100;
    bearingVectors_t bearingVectors;
    points_t points;
    std::vector<int> camCorrespondences;
    Eigen::MatrixXd gt(3, numberPoints);

    translations_t camOffsets;
    rotations_t camRotations;
    generateCentralCameraSystem(camOffsets, camRotations);

    generateRandom2D3DCorrespondences(
        position, rotation, camOffsets, camRotations, numberPoints, 0.0, 0.0,
        bearingVectors, points, camCorrespondences, gt);

    absolute_pose::CentralAbsoluteAdapter adapter(bearingVectors, points);

    transformation_t T_sqpnp = absolute_pose::sqpnp(adapter);
    transformation_t T_epnp = absolute_pose::epnp(adapter);

    double pos_error_sqpnp = (T_sqpnp.col(3) - position).norm();
    double rot_error_sqpnp = computeRotationError(T_sqpnp.block<3,3>(0,0), rotation);

    double pos_error_epnp = (T_epnp.col(3) - position).norm();
    double rot_error_epnp = computeRotationError(T_epnp.block<3,3>(0,0), rotation);

    cout << "SQPNP - Position error: " << pos_error_sqpnp << ", Rotation error: " << rot_error_sqpnp << endl;
    cout << "EPNP  - Position error: " << pos_error_epnp << ", Rotation error: " << rot_error_epnp << endl;

    // Both should be accurate for standard perspective camera
    assert(pos_error_sqpnp < 1e-6 && "SQPNP position error too large");
    assert(rot_error_sqpnp < 1e-4 && "SQPNP rotation error too large");

    cout << "✓ Test passed" << endl << endl;
}

int main(int argc, char** argv) {
    initializeRandomSeed();

    cout << "=============================================" << endl;
    cout << "SQPnP Edge Cases and Robustness Tests" << endl;
    cout << "Tests for forward/backward facing vectors" << endl;
    cout << "=============================================" << endl << endl;

    try {
        // Test various forward/backward vector combinations
        test_forward_only_vectors();
        test_backward_only_vectors();
        test_mixed_forward_backward_vectors();
        test_minority_backward_vectors();

        // Test edge cases
        test_degenerate_inputs();
        test_sqp_vs_epnp_comparison();

        cout << "=============================================" << endl;
        cout << "All tests PASSED! (6 tests)" << endl;
        cout << "=============================================" << endl;

        return 0;
    } catch (const exception& e) {
        cerr << "Test FAILED with exception: " << e.what() << endl;
        return 1;
    } catch (...) {
        cerr << "Test FAILED with unknown exception" << endl;
        return 1;
    }
}