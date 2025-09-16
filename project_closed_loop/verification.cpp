#include "../solver.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>

#ifdef __cpp_lib_math_special_functions
    // C++17 special functions are available
    #define HAS_STD_BESSEL 1
#else
    // Fall back to custom implementation
    #define HAS_STD_BESSEL 0
#endif

// Use C++17 standard library Bessel functions
// std::cyl_bessel_k(nu, x) computes K_nu(x)
// K0(x) = std::cyl_bessel_k(0, x)
// K1(x) = std::cyl_bessel_k(1, x)
mpf_class modified_bessel_k0(const mpf_class& x);
mpf_class modified_bessel_k1(const mpf_class& x);

// Constants from require.txt
const double R = 0.2;          // m
const double Q = 1.0;          // m^3/d
const double Df = 1.0;         // m^2/d  
const double Db = 1e-9;        // m^2/d (Dr in require.txt)
const double x_pos = 1000.0;   // m
const double h0 = 30.0;        // K
const int n = 300;             // Gaver-Stehfest method terms

class VerificationFunctions {
private:
    // High-precision exp function using Taylor series
    static mpf_class exp_mpf(const mpf_class& x, int terms = 50) {
        mpf_class result(1.0);
        mpf_class term(1.0);
        
        for (int i = 1; i <= terms; i++) {
            term = term * x / mpf_class(i);
            result += term;
        }
        return result;
    }

public:
    // H(p) = 30/p (Laplace transform of step function h(t) = 30)
    static mpf_class H_laplace(const mpf_class& p) {
        return mpf_class(h0) / p;
    }
    
    // gamma1 = Q / (pi * R^2 * Df)
    static mpf_class gamma1() {
        const mpf_class pi(3.141592653589793238462643383279502884197);
        return mpf_class(Q) / (pi * mpf_class(R * R) * mpf_class(Df));
    }
    
    // gamma2 = (Db/Df * 1/R * sqrt(p/Db) * K1(sqrt(p/Db)*R) / K0(sqrt(p/Db)*R) + p/Df)
    static mpf_class gamma2(const mpf_class& p) {
        const mpf_class sqrt_p_over_Db = sqrt(p / mpf_class(Db));
        const mpf_class arg = sqrt_p_over_Db * mpf_class(R);
        
        const mpf_class k1_val = modified_bessel_k1(arg);
        const mpf_class k0_val = modified_bessel_k0(arg);
        
        // Check for division by zero
        if (abs(k0_val) < mpf_class(1e-100)) {
            std::cerr << "Warning: K0 is very small, potential division by zero" << std::endl;
            return mpf_class(1e10); // Return large value to avoid division by zero
        }
        
        const mpf_class term1 = (mpf_class(Db) / mpf_class(Df)) * 
                               (mpf_class(1.0) / mpf_class(R)) * 
                               sqrt_p_over_Db * 
                               (k1_val / k0_val);
        
        const mpf_class term2 = p / mpf_class(Df);
        
        return term1 + term2;
    }
    
    // Laplace space solution: T_rf_bar = H(p) * exp((gamma1 - sqrt(gamma1^2 + 4*gamma2))/2 * x)
    static mpf_class laplace_solution(const mpf_class& p) {
        const mpf_class g1 = gamma1();
        const mpf_class g2 = gamma2(p);
        
        const mpf_class discriminant = g1 * g1 + 4.0 * g2;
        const mpf_class sqrt_discriminant = sqrt(discriminant);
        const mpf_class exponent = (g1 - sqrt_discriminant) / 2.0 * mpf_class(x_pos);
        
        return H_laplace(p) * exp_mpf(exponent);
    }
    
    // Time domain analytical solution: T_rf = h(t)/2 * erfc((x - Q/(pi*R^2)*t) / (2*sqrt(Df*t)))
    static mpf_class time_domain_solution(const mpf_class& t) {
        const mpf_class pi(3.141592653589793238462643383279502884197);
        const mpf_class velocity = mpf_class(Q) / (pi * mpf_class(R * R));
        const mpf_class advection_term = mpf_class(x_pos) - velocity * t;
        const mpf_class diffusion_term = 2.0 * sqrt(mpf_class(Df) * t);
        
        // Check for division by zero
        if (abs(diffusion_term) < mpf_class(1e-100)) {
            std::cerr << "Warning: diffusion term is very small for t=" << t.get_d() << std::endl;
            return mpf_class(0.0);
        }
        
        const mpf_class erfc_arg = advection_term / diffusion_term;
        
        // Check for reasonable erfc argument range
        double erfc_arg_double = erfc_arg.get_d();
        if (abs(erfc_arg_double) > 100.0) {
            // For very large positive arguments, erfc ≈ 0
            // For very large negative arguments, erfc ≈ 2
            if (erfc_arg_double > 100.0) {
                return mpf_class(0.0);
            } else {
                return mpf_class(h0);
            }
        }
        
        const mpf_class erfc_val = complementary_error_function(erfc_arg);
        
        return mpf_class(h0) / 2.0 * erfc_val;
    }
    
private:
    // Complementary error function implementation for mpf_class
    static mpf_class complementary_error_function(const mpf_class& x) {
        // For now, use a series approximation or convert to double for std::erfc
        // This is a placeholder - should be implemented with high precision
        double x_double = x.get_d();
        return mpf_class(std::erfc(x_double));
    }
};

// Helper function for high-precision exponential
mpf_class exp_mpf(const mpf_class& x, int terms = 50) {
    mpf_class result(1.0);
    mpf_class term(1.0);
    
    for (int i = 1; i <= terms; i++) {
        term = term * x / mpf_class(i);
        result += term;
    }
    return result;
}

// GSL-based Bessel function implementations 
// Using GSL's well-tested Bessel functions with error handling
mpf_class modified_bessel_k0(const mpf_class& x) {
    double x_val = x.get_d();
    
    // Check for reasonable input range
    if (x_val <= 0.0) {
        std::cerr << "Warning: K0 called with non-positive argument: " << x_val << std::endl;
        return mpf_class(1e100); // Large value to indicate divergence
    }
    
    // For very large x, K0 becomes extremely small (underflow)
    if (x_val > 700.0) {
        return mpf_class(0.0); // Return 0 for very large arguments
    }
    
    // Use GSL's modified Bessel function K0
    gsl_sf_result result;
    int status = gsl_sf_bessel_K0_e(x_val, &result);
    
    if (status != GSL_SUCCESS) {
        // Handle GSL error gracefully
        if (x_val > 100.0) {
            // For large x, use asymptotic approximation
            const double pi = 3.141592653589793;
            return mpf_class(sqrt(pi / (2.0 * x_val)) * exp(-x_val));
        } else {
            std::cerr << "GSL K0 error for x=" << x_val << ", status=" << status << std::endl;
            return mpf_class(0.0);
        }
    }
    
    return mpf_class(result.val);
}

mpf_class modified_bessel_k1(const mpf_class& x) {
    double x_val = x.get_d();
    
    // Check for reasonable input range
    if (x_val <= 0.0) {
        std::cerr << "Warning: K1 called with non-positive argument: " << x_val << std::endl;
        return mpf_class(1e100); // Large value to indicate divergence
    }
    
    // For very large x, K1 becomes extremely small (underflow)
    if (x_val > 700.0) {
        return mpf_class(0.0); // Return 0 for very large arguments
    }
    
    // Use GSL's modified Bessel function K1
    gsl_sf_result result;
    int status = gsl_sf_bessel_K1_e(x_val, &result);
    
    if (status != GSL_SUCCESS) {
        // Handle GSL error gracefully
        if (x_val > 100.0) {
            // For large x, use asymptotic approximation
            const double pi = 3.141592653589793;
            return mpf_class(sqrt(pi / (2.0 * x_val)) * exp(-x_val));
        } else {
            std::cerr << "GSL K1 error for x=" << x_val << ", status=" << status << std::endl;
            return mpf_class(0.0);
        }
    }
    
    return mpf_class(result.val);
}

int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    // Set GSL error handler to not abort on errors
    gsl_set_error_handler_off();
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (rank == 0) {
        std::cout << "=== Verification of Inverse Laplace Transform ===" << std::endl;
        std::cout << "Comparing Laplace space solution with time domain analytical solution" << std::endl;
        std::cout << "Running on " << size << " MPI process(es)" << std::endl;
        std::cout << "Parameters:" << std::endl;
        std::cout << "  R = " << R << " m" << std::endl;
        std::cout << "  Q = " << Q << " m³/d" << std::endl;
        std::cout << "  Df = " << Df << " m²/d" << std::endl;
        std::cout << "  Db = " << Db << " m²/d" << std::endl;
        std::cout << "  x = " << x_pos << " m" << std::endl;
        std::cout << "  h(t) = " << h0 << " K" << std::endl;
        std::cout << "  n = " << n << " (Gaver-Stehfest terms)" << std::endl;
        std::cout << std::endl;
    }
    
    // Initialize solver
    InverseLaplaceSolver solver(100); // High precision
    
    // Set the Laplace space function
    solver.set_target_function(VerificationFunctions::laplace_solution);
    
    // Set the analytical solution for comparison
    solver.set_analytical_solution(VerificationFunctions::time_domain_solution);
    
    // Generate t values from 100 to 1000 (log scale, step 0.2 for faster computation)
    std::vector<mpf_class> t_values = solver.generate_t_values(2.0, 3.0, 0.2);
    
    if (rank == 0) {
        std::cout << "Generated " << t_values.size() << " t values from 100 to 1000" << std::endl;
        std::cout << "Running inverse Laplace transform with n = " << n << " terms..." << std::endl;
    }
    
    // Solve the inverse Laplace transform
    SolverResults results = solver.solve(n, t_values);
    
    if (rank == 0) {
        // Print some results
        solver.print_results(results, 20); // Show 20 sample points
        
        // Write results to file for Python plotting
        std::string filename = "verification_results.txt";
        solver.write_results_to_file(results, filename);
        
        std::cout << std::endl;
        std::cout << "Results written to: " << filename << std::endl;
        std::cout << "RMSE Error: " << results.rmse_error << std::endl;
        std::cout << "Calculation time: " << results.calculation_time << " seconds" << std::endl;
        std::cout << "Time per point: " << results.time_per_point << " seconds" << std::endl;
    }
    
    // Finalize MPI
    MPI_Finalize();
    return 0;
}