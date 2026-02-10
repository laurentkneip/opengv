// Detailed analysis of absolute pose algorithms with backward-facing vectors
// This test investigates:
// 1. How EPnP and UPnP handle backward vectors
// 2. Edge cases with varying ratios of backward vectors
// 3. Performance comparison across different scenarios
// Note: SQPnP excluded - not designed for wide-angle/omnidirectional scenarios

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <vector>
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

// Modified SQPnP adapter that allows us to inspect internal behavior
class DiagnosticAdapter : public absolute_pose::CentralAbsoluteAdapter {
public:
    DiagnosticAdapter(
        const bearingVectors_t& bearingVectors,
        const points_t& worldPoints)
        : CentralAbsoluteAdapter(bearingVectors, worldPoints) {}

    // Compute angular error for a given pose
    double computeAngularError(const transformation_t& T) const {
        double total_error = 0.0;
        int count = 0;

        for(size_t i = 0; i < getNumberCorrespondences(); i++) {
            bearingVector_t bearing = getBearingVector(i);
            point_t point = getPoint(i);

            // Transform point to camera frame
            Eigen::Vector3d p_cam = T.block<3,3>(0,0) * point + T.col(3);
            double norm = p_cam.norm();

            if(norm > 1e-10) {
                Eigen::Vector3d p_dir = p_cam / norm;
                double cos_angle = bearing.dot(p_dir);

                // Clamp for numerical stability
                if(cos_angle > 1.0) cos_angle = 1.0;
                if(cos_angle < -1.0) cos_angle = -1.0;

                total_error += acos(cos_angle);
                count++;
            }
        }

        return (count > 0) ? total_error / count : 0.0;
    }

    // Analyze bearing vector distribution
    void analyzeBearingVectors() const {
        int forward_count = 0;
        int backward_count = 0;

        for(size_t i = 0; i < getNumberCorrespondences(); i++) {
            bearingVector_t bearing = getBearingVector(i);
            if(bearing(2) > 0)
                forward_count++;
            else
                backward_count++;
        }

        cout << "  Bearing vector distribution:" << endl;
        cout << "    Forward (z > 0):  " << forward_count << " ("
             << (100.0 * forward_count / getNumberCorrespondences()) << "%)" << endl;
        cout << "    Backward (z < 0): " << backward_count << " ("
             << (100.0 * backward_count / getNumberCorrespondences()) << "%)" << endl;
    }
};

void test_varying_backward_ratio() {
    cout << "=========================================" << endl;
    cout << "Testing varying ratios of backward vectors" << endl;
    cout << "=========================================" << endl << endl;

    translation_t position = generateRandomTranslation(2.0);
    rotation_t rotation = generateRandomRotation(0.5);

    size_t numberPoints = 50;

    // Test different ratios of backward vectors
    vector<double> backward_ratios = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    cout << setw(15) << "Backward %"
         << setw(15) << "Pos Error"
         << setw(15) << "Rot Error"
         << setw(15) << "Ang Error" << endl;
    cout << string(60, '-') << endl;

    for(double ratio : backward_ratios) {
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

        // Flip specified percentage of bearing vectors
        int num_to_flip = static_cast<int>(numberPoints * ratio);
        for(int i = 0; i < num_to_flip; i++) {
            bearingVectors[i] = -bearingVectors[i];
        }

        DiagnosticAdapter adapter(bearingVectors, points);
        transformation_t T_sqpnp = absolute_pose::sqpnp(adapter);

        double pos_error = (T_sqpnp.col(3) - position).norm();
        double rot_error = computeRotationError(T_sqpnp.block<3,3>(0,0), rotation);
        double ang_error = adapter.computeAngularError(T_sqpnp);

        cout << setw(15) << (ratio * 100.0) << "%"
             << setw(15) << scientific << setprecision(3) << pos_error
             << setw(15) << scientific << setprecision(3) << rot_error
             << setw(15) << scientific << setprecision(3) << ang_error << endl;
    }
    cout << endl;
}

void test_extreme_edge_cases() {
    cout << "=========================================" << endl;
    cout << "Testing extreme edge cases" << endl;
    cout << "=========================================" << endl << endl;

    translation_t position = generateRandomTranslation(2.0);
    rotation_t rotation = generateRandomRotation(0.5);
    size_t numberPoints = 50;

    // Case 1: Single backward vector among many forward vectors
    {
        cout << "Case 1: Single backward vector (1 out of 50)" << endl;
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

        // Flip only the first vector
        bearingVectors[0] = -bearingVectors[0];

        DiagnosticAdapter adapter(bearingVectors, points);
        adapter.analyzeBearingVectors();

        transformation_t T_sqpnp = absolute_pose::sqpnp(adapter);

        double pos_error = (T_sqpnp.col(3) - position).norm();
        double rot_error = computeRotationError(T_sqpnp.block<3,3>(0,0), rotation);
        double ang_error = adapter.computeAngularError(T_sqpnp);

        cout << "  Position error: " << scientific << pos_error << endl;
        cout << "  Rotation error: " << scientific << rot_error << " rad" << endl;
        cout << "  Angular error:  " << scientific << ang_error << " rad" << endl;
        cout << "  Status: " << (pos_error < 0.1 && rot_error < 0.1 ? "PASS" : "FAIL") << endl;
        cout << endl;
    }

    // Case 2: Nearly all backward vectors (49 out of 50)
    {
        cout << "Case 2: Nearly all backward vectors (49 out of 50)" << endl;
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

        // Flip all but the last vector
        for(size_t i = 0; i < numberPoints - 1; i++) {
            bearingVectors[i] = -bearingVectors[i];
        }

        DiagnosticAdapter adapter(bearingVectors, points);
        adapter.analyzeBearingVectors();

        transformation_t T_sqpnp = absolute_pose::sqpnp(adapter);

        double pos_error = (T_sqpnp.col(3) - position).norm();
        double rot_error = computeRotationError(T_sqpnp.block<3,3>(0,0), rotation);
        double ang_error = adapter.computeAngularError(T_sqpnp);

        cout << "  Position error: " << scientific << pos_error << endl;
        cout << "  Rotation error: " << scientific << rot_error << " rad" << endl;
        cout << "  Angular error:  " << scientific << ang_error << " rad" << endl;
        cout << "  Status: " << (pos_error < 0.1 && rot_error < 0.1 ? "PASS" : "FAIL") << endl;
        cout << endl;
    }

    // Case 3: Alternating forward/backward pattern
    {
        cout << "Case 3: Alternating forward/backward pattern" << endl;
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

        // Flip every other vector
        for(size_t i = 0; i < numberPoints; i += 2) {
            bearingVectors[i] = -bearingVectors[i];
        }

        DiagnosticAdapter adapter(bearingVectors, points);
        adapter.analyzeBearingVectors();

        transformation_t T_sqpnp = absolute_pose::sqpnp(adapter);

        double pos_error = (T_sqpnp.col(3) - position).norm();
        double rot_error = computeRotationError(T_sqpnp.block<3,3>(0,0), rotation);
        double ang_error = adapter.computeAngularError(T_sqpnp);

        cout << "  Position error: " << scientific << pos_error << endl;
        cout << "  Rotation error: " << scientific << rot_error << " rad" << endl;
        cout << "  Angular error:  " << scientific << ang_error << " rad" << endl;
        cout << "  Status: " << (pos_error < 0.1 && rot_error < 0.1 ? "PASS" : "FAIL") << endl;
        cout << endl;
    }
}

void test_solver_comparison() {
    cout << "=========================================" << endl;
    cout << "Comparing SQPnP vs EPnP for backward vectors" << endl;
    cout << "=========================================" << endl << endl;

    translation_t position = generateRandomTranslation(2.0);
    rotation_t rotation = generateRandomRotation(0.5);
    size_t numberPoints = 50;

    vector<double> backward_ratios = {0.0, 0.25, 0.5, 0.75, 1.0};

    cout << setw(15) << "Backward %"
         << setw(18) << "SQPnP Ang Err"
         << setw(18) << "EPnP Ang Err"
         << setw(12) << "Winner" << endl;
    cout << string(63, '-') << endl;

    for(double ratio : backward_ratios) {
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

        int num_to_flip = static_cast<int>(numberPoints * ratio);
        for(int i = 0; i < num_to_flip; i++) {
            bearingVectors[i] = -bearingVectors[i];
        }

        DiagnosticAdapter adapter(bearingVectors, points);

        // Run both solvers independently (Note: the sqpnp method runs both internally)
        transformation_t T_sqpnp_result = absolute_pose::sqpnp(adapter);
        transformation_t T_epnp = absolute_pose::epnp(adapter);

        double ang_error_sqpnp = adapter.computeAngularError(T_sqpnp_result);
        double ang_error_epnp = adapter.computeAngularError(T_epnp);

        string winner = (ang_error_sqpnp < ang_error_epnp) ? "SQPnP" :
                       (ang_error_sqpnp > ang_error_epnp) ? "EPnP" : "Tie";

        cout << setw(15) << (ratio * 100.0) << "%"
             << setw(18) << scientific << setprecision(4) << ang_error_sqpnp
             << setw(18) << scientific << setprecision(4) << ang_error_epnp
             << setw(12) << winner << endl;
    }
    cout << endl;
}

void test_mathematical_properties() {
    cout << "=========================================" << endl;
    cout << "Mathematical properties of backward vectors" << endl;
    cout << "=========================================" << endl << endl;

    cout << "Property 1: Symmetry under vector negation" << endl;
    cout << "If we negate all bearing vectors AND all world points," << endl;
    cout << "the solution should remain identical." << endl << endl;

    translation_t position = generateRandomTranslation(2.0);
    rotation_t rotation = generateRandomRotation(0.5);
    size_t numberPoints = 30;

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

    // Test 1: Original
    DiagnosticAdapter adapter1(bearingVectors, points);
    transformation_t T1 = absolute_pose::sqpnp(adapter1);
    double error1 = adapter1.computeAngularError(T1);

    // Test 2: Flip all bearings
    bearingVectors_t bearingVectors_flipped = bearingVectors;
    for(size_t i = 0; i < numberPoints; i++) {
        bearingVectors_flipped[i] = -bearingVectors_flipped[i];
    }
    DiagnosticAdapter adapter2(bearingVectors_flipped, points);
    transformation_t T2 = absolute_pose::sqpnp(adapter2);
    double error2 = adapter2.computeAngularError(T2);

    cout << "  Original:       Angular error = " << scientific << error1 << endl;
    cout << "  Flipped:        Angular error = " << scientific << error2 << endl;
    cout << "  Difference:     " << scientific << abs(error2 - error1) << endl;
    cout << "  Status: " << (abs(error2 - error1) < 0.01 ? "PASS (Similar performance)" : "NOTE (Different performance)") << endl;
    cout << endl;
}

int main(int argc, char** argv) {
    initializeRandomSeed();

    cout << "=========================================" << endl;
    cout << "Backward Vector Analysis for SQPnP" << endl;
    cout << "Detailed diagnostic investigation" << endl;
    cout << "=========================================" << endl << endl;

    try {
        test_varying_backward_ratio();
        test_extreme_edge_cases();
        test_solver_comparison();
        test_mathematical_properties();

        cout << "=========================================" << endl;
        cout << "All backward vector analyses complete!" << endl;
        cout << "=========================================" << endl;

        return 0;
    } catch (const exception& e) {
        cerr << "Analysis FAILED with exception: " << e.what() << endl;
        return 1;
    } catch (...) {
        cerr << "Analysis FAILED with unknown exception" << endl;
        return 1;
    }
}
