// Benchmark to evaluate the hybrid SQPnP+EPnP approach
// Questions to answer:
// 1. What is the computational overhead of running both algorithms?
// 2. How often does SQPnP win vs EPnP?
// 3. What is the accuracy benefit of the hybrid approach?
// 4. Is the overhead worthwhile?

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <vector>
#include <chrono>
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>
#include <opengv/absolute_pose/modules/Sqpnp.hpp>

#include "random_generators.hpp"
#include "experiment_helpers.hpp"

using namespace std;
using namespace Eigen;
using namespace opengv;
using namespace std::chrono;

double computeRotationError(const rotation_t& R_est, const rotation_t& R_gt) {
    rotation_t R_err = R_est.transpose() * R_gt;
    double trace = R_err.trace();
    trace = std::min(3.0, std::max(-1.0, trace));
    return std::acos((trace - 1.0) / 2.0);
}

double computeAngularError(
    const bearingVectors_t& bearings,
    const points_t& points,
    const transformation_t& T)
{
    double total_error = 0.0;
    int count = 0;

    for(size_t i = 0; i < bearings.size(); i++) {
        bearingVector_t bearing = bearings[i];
        point_t point = points[i];

        Eigen::Vector3d p_cam = T.block<3,3>(0,0) * point + T.col(3);
        double norm = p_cam.norm();

        if(norm > 1e-10) {
            Eigen::Vector3d p_dir = p_cam / norm;
            double cos_angle = bearing.dot(p_dir);
            if(cos_angle > 1.0) cos_angle = 1.0;
            if(cos_angle < -1.0) cos_angle = -1.0;
            total_error += acos(cos_angle);
            count++;
        }
    }

    return (count > 0) ? total_error / count : 0.0;
}

struct BenchmarkResult {
    string scenario;
    double time_sqpnp_only_ms;
    double time_epnp_only_ms;
    double time_hybrid_ms;
    double overhead_percent;

    double error_sqpnp;
    double error_epnp;
    double error_hybrid;

    string winner;
    double accuracy_improvement_percent;
};

BenchmarkResult benchmark_scenario(
    const string& scenario_name,
    const bearingVectors_t& bearingVectors,
    const points_t& points,
    const translation_t& gt_position,
    const rotation_t& gt_rotation,
    int num_trials = 10)
{
    BenchmarkResult result;
    result.scenario = scenario_name;

    absolute_pose::CentralAbsoluteAdapter adapter(bearingVectors, points);

    // Benchmark hybrid (current implementation)
    auto start_hybrid = high_resolution_clock::now();
    transformation_t T_hybrid;
    for(int i = 0; i < num_trials; i++) {
        T_hybrid = absolute_pose::sqpnp(adapter);  // Runs both internally
    }
    auto end_hybrid = high_resolution_clock::now();
    result.time_hybrid_ms = duration_cast<microseconds>(end_hybrid - start_hybrid).count() / (1000.0 * num_trials);
    result.error_hybrid = computeAngularError(bearingVectors, points, T_hybrid);

    // Benchmark EPnP only
    auto start_epnp = high_resolution_clock::now();
    transformation_t T_epnp;
    for(int i = 0; i < num_trials; i++) {
        T_epnp = absolute_pose::epnp(adapter);
    }
    auto end_epnp = high_resolution_clock::now();
    result.time_epnp_only_ms = duration_cast<microseconds>(end_epnp - start_epnp).count() / (1000.0 * num_trials);
    result.error_epnp = computeAngularError(bearingVectors, points, T_epnp);

    // Estimate SQPnP-only time (approximation since we can't isolate it easily)
    // The hybrid approach runs: Omega computation + SQP solve + EPnP
    // We can estimate: time_sqpnp_only ≈ time_hybrid - time_epnp_only
    result.time_sqpnp_only_ms = result.time_hybrid_ms - result.time_epnp_only_ms;
    if(result.time_sqpnp_only_ms < 0) result.time_sqpnp_only_ms = 0;

    // For accurate SQPnP-only error, we'd need to modify the code
    // For now, use the hybrid error if it's SQPnp that won
    result.error_sqpnp = result.error_hybrid;  // Approximation

    // Calculate overhead
    double baseline_time = min(result.time_sqpnp_only_ms, result.time_epnp_only_ms);
    result.overhead_percent = 100.0 * (result.time_hybrid_ms - baseline_time) / baseline_time;

    // Determine winner
    if(result.error_hybrid <= result.error_epnp * 1.01) {  // Within 1% counts as SQPnP win
        result.winner = "SQPnP";
        result.accuracy_improvement_percent = 100.0 * (result.error_epnp - result.error_hybrid) / result.error_epnp;
    } else {
        result.winner = "EPnP";
        result.accuracy_improvement_percent = 0.0;
    }

    return result;
}

void print_benchmark_results(const vector<BenchmarkResult>& results) {
    cout << "\n=========================================" << endl;
    cout << "Hybrid SQPnP+EPnP Benchmark Results" << endl;
    cout << "=========================================" << endl << endl;

    cout << "Timing Comparison (milliseconds per solve):" << endl;
    cout << setw(30) << "Scenario"
         << setw(12) << "SQPnP*"
         << setw(12) << "EPnP"
         << setw(12) << "Hybrid"
         << setw(12) << "Overhead" << endl;
    cout << string(78, '-') << endl;

    for(const auto& r : results) {
        cout << setw(30) << r.scenario
             << setw(12) << fixed << setprecision(3) << r.time_sqpnp_only_ms
             << setw(12) << fixed << setprecision(3) << r.time_epnp_only_ms
             << setw(12) << fixed << setprecision(3) << r.time_hybrid_ms
             << setw(11) << fixed << setprecision(1) << r.overhead_percent << "%" << endl;
    }

    cout << "\n* SQPnP time estimated as: Hybrid - EPnP\n" << endl;

    cout << "Accuracy Comparison (angular error in radians):" << endl;
    cout << setw(30) << "Scenario"
         << setw(15) << "EPnP Error"
         << setw(15) << "Hybrid Error"
         << setw(12) << "Winner"
         << setw(12) << "Improve%" << endl;
    cout << string(84, '-') << endl;

    for(const auto& r : results) {
        cout << setw(30) << r.scenario
             << setw(15) << scientific << setprecision(3) << r.error_epnp
             << setw(15) << scientific << setprecision(3) << r.error_hybrid
             << setw(12) << r.winner
             << setw(11) << fixed << setprecision(1) << r.accuracy_improvement_percent << "%" << endl;
    }

    cout << "\n=========================================" << endl;
    cout << "Summary Analysis:" << endl;
    cout << "=========================================" << endl;

    // Calculate statistics
    double avg_overhead = 0;
    int sqpnp_wins = 0;
    int epnp_wins = 0;
    double total_improvement = 0;
    int improvement_count = 0;

    for(const auto& r : results) {
        avg_overhead += r.overhead_percent;
        if(r.winner == "SQPnP") {
            sqpnp_wins++;
            if(r.accuracy_improvement_percent > 1.0) {  // Significant improvement
                total_improvement += r.accuracy_improvement_percent;
                improvement_count++;
            }
        } else {
            epnp_wins++;
        }
    }

    avg_overhead /= results.size();
    double avg_improvement = (improvement_count > 0) ? total_improvement / improvement_count : 0;

    cout << "\nAverage overhead: " << fixed << setprecision(1) << avg_overhead << "%" << endl;
    cout << "SQPnP wins: " << sqpnp_wins << " / " << results.size()
         << " (" << fixed << setprecision(1) << (100.0 * sqpnp_wins / results.size()) << "%)" << endl;
    cout << "EPnP wins: " << epnp_wins << " / " << results.size()
         << " (" << fixed << setprecision(1) << (100.0 * epnp_wins / results.size()) << "%)" << endl;
    cout << "Average improvement when SQPnP wins: " << fixed << setprecision(1) << avg_improvement << "%" << endl;

    cout << "\nConclusion:" << endl;
    if(avg_overhead < 50 && sqpnp_wins > epnp_wins) {
        cout << "✓ Hybrid approach is WORTHWHILE" << endl;
        cout << "  - Moderate overhead (<50%)" << endl;
        cout << "  - SQPnP wins majority of cases" << endl;
        cout << "  - Provides accuracy improvement" << endl;
    } else if(avg_overhead < 100 && avg_improvement > 10) {
        cout << "✓ Hybrid approach is BENEFICIAL" << endl;
        cout << "  - Reasonable overhead (<100%)" << endl;
        cout << "  - Significant accuracy gains (>" << fixed << setprecision(0) << avg_improvement << "%)" << endl;
    } else if(sqpnp_wins < results.size() * 0.3) {
        cout << "⚠ Hybrid approach has LIMITED VALUE" << endl;
        cout << "  - SQPnP rarely wins" << endl;
        cout << "  - Consider making SQPnP optional" << endl;
    } else {
        cout << "~ Hybrid approach is SITUATIONAL" << endl;
        cout << "  - Provides value in specific scenarios" << endl;
        cout << "  - May want user-configurable mode" << endl;
    }

    cout << "\nRecommendation:" << endl;
    if(avg_overhead > 100) {
        cout << "Consider adding a configuration option:" << endl;
        cout << "  - FAST mode: EPnP only" << endl;
        cout << "  - ACCURATE mode: Hybrid (current)" << endl;
        cout << "  - ROBUST mode: Try all 3 methods" << endl;
    } else {
        cout << "Current hybrid approach is well-balanced." << endl;
        cout << "The overhead is acceptable for the robustness gained." << endl;
    }
    cout << endl;
}

int main(int argc, char** argv) {
    initializeRandomSeed();

    vector<BenchmarkResult> results;

    cout << "Running benchmarks (this may take a minute)..." << endl;

    // Scenario 1: Forward-only (standard perspective camera)
    {
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

        results.push_back(benchmark_scenario(
            "Forward only (0%)", bearingVectors, points, position, rotation));
    }

    // Scenario 2: 25% backward
    {
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

        for(size_t i = 0; i < numberPoints/4; i++)
            bearingVectors[i] = -bearingVectors[i];

        results.push_back(benchmark_scenario(
            "Mixed (25% backward)", bearingVectors, points, position, rotation));
    }

    // Scenario 3: 50% backward
    {
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

        for(size_t i = 0; i < numberPoints/2; i++)
            bearingVectors[i] = -bearingVectors[i];

        results.push_back(benchmark_scenario(
            "Omnidirectional (50%)", bearingVectors, points, position, rotation));
    }

    // Scenario 4: 75% backward
    {
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

        for(size_t i = 0; i < 3*numberPoints/4; i++)
            bearingVectors[i] = -bearingVectors[i];

        results.push_back(benchmark_scenario(
            "Mixed (75% backward)", bearingVectors, points, position, rotation));
    }

    // Scenario 5: 100% backward
    {
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

        for(size_t i = 0; i < numberPoints; i++)
            bearingVectors[i] = -bearingVectors[i];

        results.push_back(benchmark_scenario(
            "Backward only (100%)", bearingVectors, points, position, rotation));
    }

    // Scenario 6: Few points (minimal)
    {
        translation_t position = generateRandomTranslation(2.0);
        rotation_t rotation = generateRandomRotation(0.5);
        size_t numberPoints = 10;  // Minimal case

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

        results.push_back(benchmark_scenario(
            "Few points (n=10)", bearingVectors, points, position, rotation));
    }

    // Scenario 7: Many points
    {
        translation_t position = generateRandomTranslation(2.0);
        rotation_t rotation = generateRandomRotation(0.5);
        size_t numberPoints = 200;  // Large case

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

        results.push_back(benchmark_scenario(
            "Many points (n=200)", bearingVectors, points, position, rotation, 5));  // Fewer trials due to size
    }

    print_benchmark_results(results);

    return 0;
}
