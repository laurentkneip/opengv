// Proof-of-concept: Parallel Hybrid SQPnP Implementation
// This file demonstrates how to implement parallel execution
// To use: Add compute_pose_parallel() to Sqpnp class

#include <future>
#include <opengv/absolute_pose/modules/Sqpnp.hpp>

namespace opengv {
namespace absolute_pose {
namespace modules {

// Helper struct to pass results between threads
struct SqpnpPoseResult {
    double R[3][3];
    double t[3];
    double error;

    SqpnpPoseResult() : error(std::numeric_limits<double>::max()) {
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++)
                R[i][j] = 0.0;
            t[i] = 0.0;
        }
    }
};

// Isolated SQPnP computation (thread-safe)
SqpnpPoseResult compute_sqpnp_isolated(
    const double* pws,
    const double* us,
    const int* signs,
    int n)
{
    // Create temporary Sqpnp instance for this thread
    Sqpnp solver;
    solver.set_maximum_number_of_correspondences(n);
    solver.reset_correspondences();

    // Copy correspondences
    for(int i = 0; i < n; i++) {
        solver.add_correspondence(
            pws[3*i], pws[3*i+1], pws[3*i+2],
            us[3*i], us[3*i+1], us[3*i+2]);
    }

    // Compute Omega matrix and run SQP solver
    Eigen::MatrixXd Omega;
    solver.compute_omega_matrix(Omega);

    Eigen::Matrix3d R_sqp;
    Eigen::Vector3d t_sqp;
    double error = solver.sqp_solve(Omega, R_sqp, t_sqp);

    // Package result
    SqpnpPoseResult result;
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++)
            result.R[i][j] = R_sqp(i, j);
        result.t[i] = t_sqp(i);
    }
    result.error = error;

    return result;
}

// Isolated EPnP computation (thread-safe)
SqpnpPoseResult compute_epnp_isolated(
    const double* pws,
    const double* us,
    const int* signs,
    int n)
{
    // Create temporary Sqpnp instance for this thread
    Sqpnp solver;
    solver.set_maximum_number_of_correspondences(n);
    solver.reset_correspondences();

    // Copy correspondences
    for(int i = 0; i < n; i++) {
        solver.add_correspondence(
            pws[3*i], pws[3*i+1], pws[3*i+2],
            us[3*i], us[3*i+1], us[3*i+2]);
    }

    // Run EPnP-based algorithm
    solver.choose_control_points();
    solver.compute_barycentric_coordinates();

    Eigen::MatrixXd M(2*n, 12);
    for(int i = 0; i < n; i++) {
        solver.fill_M(M, 2*i,
                     solver.get_alphas() + 4*i,
                     us[3*i], us[3*i+1], us[3*i+2]);
    }

    Eigen::MatrixXd MtM = M.transpose() * M;
    Eigen::JacobiSVD<Eigen::MatrixXd> SVD(
        MtM, Eigen::ComputeFullV | Eigen::ComputeFullU);
    Eigen::MatrixXd Ut = SVD.matrixU().transpose();

    // Apply Gram-Schmidt
    solver.gram_schmidt_orthogonalize(Ut, 6, 6);

    Eigen::Matrix<double,6,10> L_6x10;
    Eigen::Matrix<double,6,1> Rho;
    solver.compute_L_6x10(Ut, L_6x10);
    solver.compute_rho(Rho);

    // Try different beta approximations
    double Betas[4][4], rep_errors[4];
    double Rs[4][3][3], ts[4][3];

    solver.find_betas_approx_1(L_6x10, Rho, Betas[1]);
    solver.gauss_newton(L_6x10, Rho, Betas[1]);
    rep_errors[1] = solver.compute_R_and_t(Ut, Betas[1], Rs[1], ts[1]);

    solver.find_betas_approx_2(L_6x10, Rho, Betas[2]);
    solver.gauss_newton(L_6x10, Rho, Betas[2]);
    rep_errors[2] = solver.compute_R_and_t(Ut, Betas[2], Rs[2], ts[2]);

    solver.find_betas_approx_3(L_6x10, Rho, Betas[3]);
    solver.gauss_newton(L_6x10, Rho, Betas[3]);
    rep_errors[3] = solver.compute_R_and_t(Ut, Betas[3], Rs[3], ts[3]);

    // Find best EPnP solution
    int best_idx = 1;
    for(int i = 2; i <= 3; i++) {
        if(rep_errors[i] < rep_errors[best_idx])
            best_idx = i;
    }

    // Package result
    SqpnpPoseResult result;
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++)
            result.R[i][j] = Rs[best_idx][i][j];
        result.t[i] = ts[best_idx][i];
    }
    result.error = rep_errors[best_idx];

    return result;
}

// Add this method to Sqpnp class
double Sqpnp::compute_pose_parallel(double R[3][3], double t[3])
{
    // NOTE: This requires exposing pws, us, signs as const accessors
    // or making copies here

    // For now, create copies (thread-safe approach)
    int n = number_of_correspondences;
    double* pws_copy = new double[3*n];
    double* us_copy = new double[3*n];
    int* signs_copy = new int[n];

    memcpy(pws_copy, pws, 3*n*sizeof(double));
    memcpy(us_copy, us, 3*n*sizeof(double));
    memcpy(signs_copy, signs, n*sizeof(int));

    // Launch SQPnP in parallel thread
    auto sqpnp_future = std::async(std::launch::async,
        compute_sqpnp_isolated, pws_copy, us_copy, signs_copy, n);

    // Run EPnP in current thread
    SqpnpPoseResult epnp_result = compute_epnp_isolated(
        pws_copy, us_copy, signs_copy, n);

    // Wait for SQPnP to complete
    SqpnpPoseResult sqpnp_result = sqpnp_future.get();

    // Clean up copies
    delete[] pws_copy;
    delete[] us_copy;
    delete[] signs_copy;

    // Select best solution
    if(sqpnp_result.error < epnp_result.error) {
        copy_R_and_t(sqpnp_result.R, sqpnp_result.t, R, t);
        return sqpnp_result.error;
    } else {
        copy_R_and_t(epnp_result.R, epnp_result.t, R, t);
        return epnp_result.error;
    }
}

// Alternative: Zero-copy version using const methods
// Requires refactoring Sqpnp to have const compute methods

double Sqpnp::compute_pose_parallel_zerocopy(double R[3][3], double t[3])
{
    // This version avoids copies by making compute methods const-correct
    // and using separate temporary state for each algorithm

    // TODO: Requires refactoring Sqpnp class to separate:
    // - Input data (const, shareable)
    // - Working memory (thread-local)
    // - Output results

    // Design pattern:
    // struct SqpnpInput {
    //     const double* pws;
    //     const double* us;
    //     const int* signs;
    //     int n;
    // };
    //
    // struct SqpnpWorkspace {
    //     // Temporary arrays for computation
    // };
    //
    // SqpnpPoseResult compute_sqpnp(
    //     const SqpnpInput& input,
    //     SqpnpWorkspace& workspace);

    // For now, use the copy-based version above
    return compute_pose_parallel(R, t);
}

}  // namespace modules
}  // namespace absolute_pose
}  // namespace opengv


// ========================================================================
// Example: Adding to public API
// ========================================================================

namespace opengv {
namespace absolute_pose {

// In methods.hpp:
/**
 * \brief Compute camera pose using hybrid SQPnP+EPnP with parallel execution.
 *
 * Runs both SQPnP and EPnP algorithms in parallel threads and selects
 * the result with lower angular error. This provides the robustness of
 * the hybrid approach with minimal timing overhead.
 *
 * \param[in] adapter Visitor holding bearing vectors and world points.
 * \return The transformation from the world to the camera frame.
 *
 * \note Requires C++11 and multi-core CPU for performance benefit.
 * \note Uses ~2× memory during computation (released immediately).
 * \note Wall-clock time ≈ max(T_sqpnp, T_epnp) instead of sum.
 */
transformation_t sqpnp_parallel(const AbsoluteAdapterBase & adapter);

// In methods.cpp:
transformation_t sqpnp_parallel(const AbsoluteAdapterBase & adapter)
{
    modules::Sqpnp Sqpnp_obj;
    Sqpnp_obj.set_maximum_number_of_correspondences(
        adapter.getNumberCorrespondences());
    Sqpnp_obj.reset_correspondences();

    for(size_t i = 0; i < adapter.getNumberCorrespondences(); i++) {
        bearingVector_t f = adapter.getBearingVector(i);
        point_t p = adapter.getPoint(i);
        Sqpnp_obj.add_correspondence(p[0], p[1], p[2], f[0], f[1], f[2]);
    }

    transformation_t result = transformation_t::Identity();
    Sqpnp_obj.compute_pose_parallel(
        result.topLeftCorner<3,3>().data(),
        result.col(3).data());

    return result;
}

}  // namespace absolute_pose
}  // namespace opengv


// ========================================================================
// Example: Usage from Python bindings
// ========================================================================

// In pyopengv.cpp:
/*
    m.def("sqpnp_parallel",
          &absolute_pose::sqpnp_parallel,
          "adapter"_a,
          "Compute absolute pose using parallel hybrid SQPnP+EPnP.\n"
          "Faster than sqpnp() by running both algorithms simultaneously.\n\n"
          "Args:\n"
          "    adapter: Adapter with bearing vectors and world points\n\n"
          "Returns:\n"
          "    4×4 transformation matrix [R | t]\n\n"
          "Note:\n"
          "    Uses 2 CPU cores during computation for ~1.5-2× speedup.\n");
*/


// ========================================================================
// Benchmarking Code
// ========================================================================

void benchmark_sequential_vs_parallel()
{
    using namespace std::chrono;

    bearingVectors_t bearings;
    points_t points;
    // ... generate test data ...

    absolute_pose::CentralAbsoluteAdapter adapter(bearings, points);

    // Benchmark sequential
    auto t1 = high_resolution_clock::now();
    transformation_t result_seq = absolute_pose::sqpnp(adapter);
    auto t2 = high_resolution_clock::now();
    double time_seq_ms = duration_cast<microseconds>(t2 - t1).count() / 1000.0;

    // Benchmark parallel
    auto t3 = high_resolution_clock::now();
    transformation_t result_par = absolute_pose::sqpnp_parallel(adapter);
    auto t4 = high_resolution_clock::now();
    double time_par_ms = duration_cast<microseconds>(t4 - t3).count() / 1000.0;

    std::cout << "Sequential: " << time_seq_ms << " ms\n";
    std::cout << "Parallel:   " << time_par_ms << " ms\n";
    std::cout << "Speedup:    " << (time_seq_ms / time_par_ms) << "×\n";
}
