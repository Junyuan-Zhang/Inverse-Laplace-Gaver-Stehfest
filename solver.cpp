#include "solver.h"
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <cmath>

InverseLaplaceSolver::InverseLaplaceSolver(int prec) : precision(prec) {
    // Set GMP precision for floating point operations
    mpf_set_default_prec(precision * 3.32);  // Convert decimal places to bits
    
    // Initialize MPI
    int provided;
    MPI_Query_thread(&provided);
    if (provided < MPI_THREAD_MULTIPLE) {
        // MPI doesn't support multithreading, but we can still use it
    }
    
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
}

InverseLaplaceSolver::~InverseLaplaceSolver() {
    // MPI_Finalize is called in main
}

void InverseLaplaceSolver::set_target_function(std::function<mpf_class(const mpf_class&)> func) {
    target_function = func;
}

void InverseLaplaceSolver::set_analytical_solution(std::function<mpf_class(const mpf_class&)> func) {
    analytical_function = func;
    has_analytical = true;
}

mpf_class InverseLaplaceSolver::Vi(int n, int i) {
    mpf_class result(0.0);
    
    // Calculate sign: (-1)^((n/2)+i)
    int sign_exp = (n/2) + i;
    mpf_class sign = (sign_exp % 2 == 0) ? mpf_class(1.0) : mpf_class(-1.0);
    
    int min_term = std::min(i, n/2);
    int kin = (i + 1) / 2;
    
    mpf_class sum(0.0);
    
    for (int k = kin; k <= min_term; k++) {
        // Calculate k^(n/2)
        mpf_class k_power;
        mpf_pow_ui(k_power.get_mpf_t(), mpf_class(k).get_mpf_t(), n/2);
        
        // Calculate factorials
        mpz_class fact_2k, fact_n2_k, fact_k, fact_k_1, fact_i_k, fact_2k_i;
        
        mpz_fac_ui(fact_2k.get_mpz_t(), 2*k);
        mpz_fac_ui(fact_n2_k.get_mpz_t(), (n/2) - k);
        mpz_fac_ui(fact_k.get_mpz_t(), k);
        mpz_fac_ui(fact_k_1.get_mpz_t(), k - 1);
        mpz_fac_ui(fact_i_k.get_mpz_t(), i - k);
        mpz_fac_ui(fact_2k_i.get_mpz_t(), 2*k - i);
        
        // Convert factorials to floating point
        mpf_class dumm1 = k_power * mpf_class(fact_2k);
        mpf_class dumm2a = mpf_class(fact_n2_k) * mpf_class(fact_k) * 
                          mpf_class(fact_k_1) * mpf_class(fact_i_k) * mpf_class(fact_2k_i);
        
        sum += dumm1 / dumm2a;
    }
    
    return sign * sum;
}

mpf_class InverseLaplaceSolver::inverse_laplace_transform(int n, const mpf_class& t) {
    mpf_class sum(0.0);
    
    // Pre-computed ln(2) for better accuracy
    mpf_class ln2("0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186875");
    mpf_class const_term = ln2 / t;
    
    // Distribute work among MPI processes
    int terms_per_process = n / mpi_size;
    int remainder = n % mpi_size;
    
    int start_i = mpi_rank * terms_per_process + 1;
    int end_i = (mpi_rank + 1) * terms_per_process;
    
    // Handle remainder terms
    if (mpi_rank < remainder) {
        start_i += mpi_rank;
        end_i += mpi_rank + 1;
    } else {
        start_i += remainder;
        end_i += remainder;
    }
    
    mpf_class local_sum(0.0);
    
    // Calculate local portion
    for (int i = start_i; i <= end_i && i <= n; i++) {
        mpf_class arg = mpf_class(i) * ln2 / t;
        mpf_class vi_val = Vi(n, i);
        mpf_class target_val = target_function(arg);
        
        local_sum += vi_val * target_val;
    }
    
    // Reduce all partial sums - FIXED: Use master-only computation to avoid precision loss
    mpf_class global_sum(0.0);
    
    if (mpi_rank == 0) {
        // Master process calculates everything to maintain precision
        for (int i = 1; i <= n; i++) {
            mpf_class arg = mpf_class(i) * ln2 / t;
            mpf_class vi_val = Vi(n, i);
            mpf_class target_val = target_function(arg);
            
            global_sum += vi_val * target_val;
        }
    }
    
    // Note: For now we sacrifice parallelization to maintain numerical accuracy
    // Future enhancement: implement high-precision MPI communication
    
    return const_term * global_sum;
}

std::vector<mpf_class> InverseLaplaceSolver::generate_t_values(double log_start, double log_end, double step) {
    std::vector<mpf_class> t_values;
    
    for (double log_t = log_start; log_t < log_end; log_t += step) {
        double t_double = std::pow(10.0, log_t);
        mpf_class t_val(t_double);
        t_values.push_back(t_val);
    }
    
    return t_values;
}

SolverResults InverseLaplaceSolver::solve(int n, const std::vector<mpf_class>& t_values) {
    SolverResults results;
    results.numerical_results.reserve(t_values.size());
    results.analytical_results.reserve(t_values.size());
    results.t_values = t_values;
    results.mpi_processes_used = mpi_size;
    
    // Start timing
    double start_time = MPI_Wtime();
    
    if (mpi_rank == 0) {
        std::cout << "Starting inverse Laplace transform calculation..." << std::endl;
        std::cout << "Using " << mpi_size << " MPI processes" << std::endl;
        std::cout << "Processing " << t_values.size() << " t values with n = " << n << " terms" << std::endl;
        std::cout << "Progress: ";
    }
    
    for (size_t i = 0; i < t_values.size(); i++) {
        if (mpi_rank == 0 && i % (t_values.size() / 20) == 0) {
            std::cout << "*" << std::flush;
        }
        
        mpf_class t_val = t_values[i];
        mpf_class numerical = inverse_laplace_transform(n, t_val);
        results.numerical_results.push_back(numerical);
        
        if (has_analytical) {
            mpf_class analytical = analytical_function(t_val);
            results.analytical_results.push_back(analytical);
        }
    }
    
    // End timing
    double end_time = MPI_Wtime();
    results.calculation_time = end_time - start_time;
    results.time_per_point = results.calculation_time / t_values.size();
    
    if (mpi_rank == 0) {
        std::cout << " Done!" << std::endl;
    }
    
    // Calculate RMSE if analytical solution is available
    if (has_analytical) {
        results.rmse_error = calculate_rmse(results.numerical_results, results.analytical_results);
    }
    
    return results;
}

double InverseLaplaceSolver::calculate_rmse(const std::vector<mpf_class>& numerical, 
                                          const std::vector<mpf_class>& analytical) {
    if (numerical.size() != analytical.size()) {
        throw std::invalid_argument("Vector sizes must match");
    }
    
    mpf_class sum_squared_error(0.0);
    size_t n = numerical.size();
    
    for (size_t i = 0; i < n; i++) {
        mpf_class diff = numerical[i] - analytical[i];
        mpf_class squared_diff;
        mpf_pow_ui(squared_diff.get_mpf_t(), diff.get_mpf_t(), 2);
        sum_squared_error += squared_diff;
    }
    
    mpf_class mean_squared_error = sum_squared_error / mpf_class(n);
    mpf_class rmse;
    mpf_sqrt(rmse.get_mpf_t(), mean_squared_error.get_mpf_t());
    
    return rmse.get_d();
}

void InverseLaplaceSolver::print_mpf(const mpf_class& value, int decimal_places) {
    mp_exp_t exp;
    std::string str = value.get_str(exp, 10, decimal_places);
    
    if (exp > 0 && exp <= (int)str.length()) {
        std::cout << str.substr(0, exp) << "." << str.substr(exp);
    } else if (exp <= 0) {
        std::cout << "0.";
        for (int i = 0; i < -exp; i++) std::cout << "0";
        std::cout << str;
    } else {
        std::cout << str;
        for (int i = str.length(); i < exp; i++) std::cout << "0";
    }
}

void InverseLaplaceSolver::print_results(const SolverResults& results, int num_samples) {
    if (mpi_rank != 0) return;  // Only master process prints
    
    std::cout << "\n=== Results ===" << std::endl;
    
    if (results.rmse_error > 0) {
        std::cout << "RMSE Error: " << std::scientific << results.rmse_error << std::endl;
    }
    
    std::cout << "\nSample Results (first " << num_samples << " points):" << std::endl;
    std::cout << "t\t\tNumerical";
    if (has_analytical) {
        std::cout << "\t\tAnalytical\t\tError";
    }
    std::cout << std::endl;
    std::cout << "-------------------------------------------------------------------" << std::endl;
    
    int samples = std::min(num_samples, (int)results.t_values.size());
    for (int i = 0; i < samples; i++) {
        print_mpf(results.t_values[i], 6);
        std::cout << "\t";
        print_mpf(results.numerical_results[i], 8);
        
        if (has_analytical && i < (int)results.analytical_results.size()) {
            std::cout << "\t";
            print_mpf(results.analytical_results[i], 8);
            std::cout << "\t";
            
            mpf_class error = results.numerical_results[i] - results.analytical_results[i];
            print_mpf(error, 6);
        }
        std::cout << std::endl;
    }
}

void InverseLaplaceSolver::write_results_to_file(const SolverResults& results, const std::string& filename) {
    if (mpi_rank != 0) return;  // Only master process writes files
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << " for writing" << std::endl;
        return;
    }
    
    // Write header
    file << "# Inverse Laplace Transform Results" << std::endl;
    file << "# Generated by Gaver-Stehfest Algorithm" << std::endl;
    file << "# Column 1: Time (t)" << std::endl;
    file << "# Column 2: Numerical Result f(t)" << std::endl;
    if (has_analytical && !results.analytical_results.empty()) {
        file << "# Column 3: Analytical Result (if available)" << std::endl;
        file << "# Column 4: Absolute Error" << std::endl;
    }
    file << "# " << std::endl;
    
    // Set precision for output
    file << std::fixed << std::setprecision(12);
    
    // Write data
    for (size_t i = 0; i < results.t_values.size(); i++) {
        // Convert mpf_class to double for file output
        double t_val = results.t_values[i].get_d();
        double num_val = results.numerical_results[i].get_d();
        
        file << t_val << "\t" << num_val;
        
        if (has_analytical && i < results.analytical_results.size()) {
            double ana_val = results.analytical_results[i].get_d();
            double error = num_val - ana_val;
            file << "\t" << ana_val << "\t" << error;
        }
        
        file << std::endl;
    }
    
    file.close();
    std::cout << "Results written to file: " << filename << std::endl;
}

// Example Functions Implementation
namespace ExampleFunctions {
    
    mpf_class default_target_function(const mpf_class& s) {
        mpf_class s2;
        mpf_pow_ui(s2.get_mpf_t(), s.get_mpf_t(), 2);  // s^2
        
        mpf_class s_minus_2 = s - mpf_class(2.0);
        mpf_class s3;
        mpf_pow_ui(s3.get_mpf_t(), s_minus_2.get_mpf_t(), 3);  // (s-2)^3
        
        mpf_class numerator = mpf_class(5.0) * s2 - mpf_class(15.0) * s + mpf_class(7.0);
        mpf_class denominator = (s + mpf_class(1.0)) * s3;
        
        return numerator / denominator;
    }
    
    mpf_class default_analytical_solution(const mpf_class& t) {
        // Using Taylor series approximation for exp function
        auto exp_approx = [](const mpf_class& x, int terms = 50) -> mpf_class {
            mpf_class result(1.0);
            mpf_class term(1.0);
            
            for (int i = 1; i <= terms; i++) {
                term = term * x / mpf_class(i);
                result += term;
            }
            return result;
        };
        
        mpf_class exp_neg_t = exp_approx(-t);
        mpf_class exp_2t = exp_approx(mpf_class(2.0) * t);
        
        mpf_class t_squared;
        mpf_pow_ui(t_squared.get_mpf_t(), t.get_mpf_t(), 2);
        
        mpf_class term1 = -exp_neg_t;
        mpf_class term2 = -mpf_class(0.5) * t_squared * exp_2t;
        mpf_class term3 = mpf_class(2.0) * t * exp_2t;
        mpf_class term4 = exp_2t;
        
        return term1 + term2 + term3 + term4;
    }
    
    std::function<mpf_class(const mpf_class&)> exponential_decay_target(double a) {
        return [a](const mpf_class& s) -> mpf_class {
            return mpf_class(1.0) / (s + mpf_class(a));
        };
    }
    
    std::function<mpf_class(const mpf_class&)> exponential_decay_analytical(double a) {
        return [a](const mpf_class& t) -> mpf_class {
            // exp(-a*t) using Taylor series
            mpf_class arg = -mpf_class(a) * t;
            mpf_class result(1.0);
            mpf_class term(1.0);
            
            for (int i = 1; i <= 50; i++) {
                term = term * arg / mpf_class(i);
                result += term;
            }
            return result;
        };
    }
    
    mpf_class step_function_target(const mpf_class& s) {
        return mpf_class(1.0) / s;
    }
    
    mpf_class step_function_analytical(const mpf_class& t) {
        // Step function is simply 1 for t > 0
        return mpf_class(1.0);
    }
}
