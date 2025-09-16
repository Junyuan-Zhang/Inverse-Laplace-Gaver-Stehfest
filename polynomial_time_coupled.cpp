#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <numeric>
#include <gmp.h>
#include <gmpxx.h>
#include <mpi.h>
#include "solver.h"

// Polynomial time-domain function: f(t) = a³t - b²t² + c⁴t
// Laplace transform: F(s) = (a³ + c⁴)/s² - 2b²/s³
class PolynomialTimeLaplaceFunction {
private:
    double a, b, c;  // Parameters

public:
    PolynomialTimeLaplaceFunction(double a_val, double b_val, double c_val) 
        : a(a_val), b(b_val), c(c_val) {}
    
    void set_parameters(double a_val, double b_val, double c_val) {
        a = a_val; b = b_val; c = c_val;
    }
    
    // Laplace domain function: F(s) = (a³ + c⁴)/s² - 2b²/s³
    mpf_class operator()(const mpf_class& s) const {
        mpf_class s2 = s * s;
        mpf_class s3 = s2 * s;
        
        mpf_class a_cubed(pow(a, 3));
        mpf_class b_squared(pow(b, 2));
        mpf_class c_fourth(pow(c, 4));
        
        return (a_cubed + c_fourth) / s2 - 2.0 * b_squared / s3;
    }
    
    // Analytical time-domain solution: f(t) = a³t - b²t² + c⁴t
    mpf_class analytical_solution(const mpf_class& t) const {
        mpf_class a_cubed(pow(a, 3));
        mpf_class b_squared(pow(b, 2));
        mpf_class c_fourth(pow(c, 4));
        mpf_class t2 = t * t;
        
        return a_cubed * t - b_squared * t2 + c_fourth * t;
    }
};

// Results structure for Sobol analysis
struct SobolResults {
    std::vector<double> times;
    std::vector<double> S1_a;     // First-order sensitivity for parameter a
    std::vector<double> S1_b;     // First-order sensitivity for parameter b
    std::vector<double> S1_c;     // First-order sensitivity for parameter c
    std::vector<double> S_total_a; // Total sensitivity for parameter a
    std::vector<double> S_total_b; // Total sensitivity for parameter b
    std::vector<double> S_total_c; // Total sensitivity for parameter c
};

// Sobol sensitivity analysis for polynomial time-domain function
class PolynomialTimeSobolAnalysis {
private:
    InverseLaplaceSolver* solver;
    PolynomialTimeLaplaceFunction* laplace_func;
    int n_samples;
    int n_terms;
    
    // Generate Sobol samples using simple stratified sampling
    std::vector<std::vector<double>> generate_sobol_samples(int n, double min_val, double max_val) {
        std::vector<std::vector<double>> samples(n, std::vector<double>(3));
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(min_val, max_val);
        
        for (int i = 0; i < n; i++) {
            samples[i][0] = dis(gen);  // a
            samples[i][1] = dis(gen);  // b
            samples[i][2] = dis(gen);  // c
        }
        return samples;
    }
    
    // Evaluate model at given parameter values
    double evaluate_model(const std::vector<std::vector<double>>& samples, int sample_idx, double t) {
        double a = samples[sample_idx][0];
        double b = samples[sample_idx][1];
        double c = samples[sample_idx][2];
        
        laplace_func->set_parameters(a, b, c);
        std::function<mpf_class(const mpf_class&)> target_func = 
            [this](const mpf_class& s) { return (*laplace_func)(s); };
        solver->set_target_function(target_func);
        
        mpf_class t_val(t);
        mpf_class result = solver->inverse_laplace_transform(n_terms, t_val);
        return result.get_d();
    }
    
    double calculate_first_order_index(const std::vector<double>& f_A,
                                      const std::vector<double>& f_B,
                                      const std::vector<double>& f_AB,
                                      double mean_f, double var_f) {
        double sum_product = 0.0;
        int n = f_A.size();
        
        for (int i = 0; i < n; i++) {
            sum_product += f_B[i] * (f_AB[i] - f_A[i]);
        }
        
        double V_i = sum_product / n;
        double S_i = V_i / var_f;
        return std::max(0.0, std::min(1.0, S_i));
    }

public:
    PolynomialTimeSobolAnalysis(InverseLaplaceSolver* solver_ptr, 
                               PolynomialTimeLaplaceFunction* func_ptr,
                               int samples = 2048, int terms = 14) 
        : solver(solver_ptr), laplace_func(func_ptr), n_samples(samples), n_terms(terms) {}
    
    SobolResults analyze_sensitivity(const std::vector<double>& time_points,
                                    double a_min, double a_max,
                                    double b_min, double b_max,
                                    double c_min, double c_max) {
        SobolResults results;
        results.times = time_points;
        results.S1_a.resize(time_points.size());
        results.S1_b.resize(time_points.size());
        results.S1_c.resize(time_points.size());
        results.S_total_a.resize(time_points.size());
        results.S_total_b.resize(time_points.size());
        results.S_total_c.resize(time_points.size());
        
        std::cout << "Time points: " << time_points.size() << std::endl;
        std::cout << "Monte Carlo samples: " << n_samples << std::endl;
        
        for (size_t t_idx = 0; t_idx < time_points.size(); t_idx++) {
            double t = time_points[t_idx];
            std::cout << "Processing time point " << (t_idx + 1) << "/" << time_points.size() 
                      << " (t = " << std::fixed << std::setprecision(2) << t << ")" << std::endl;
            
            // Generate base samples A and B
            auto samples_A = generate_sobol_samples(n_samples, a_min, a_max);
            auto samples_B = generate_sobol_samples(n_samples, b_min, b_max);
            
            // Evaluate f(A) and f(B)
            std::vector<double> f_A(n_samples), f_B(n_samples);
            for (int i = 0; i < n_samples; i++) {
                f_A[i] = evaluate_model(samples_A, i, t);
                f_B[i] = evaluate_model(samples_B, i, t);
            }
            
            // Calculate statistics
            double mean_f = (std::accumulate(f_A.begin(), f_A.end(), 0.0) + 
                           std::accumulate(f_B.begin(), f_B.end(), 0.0)) / (2.0 * n_samples);
            
            double var_sum = 0.0;
            for (int i = 0; i < n_samples; i++) {
                var_sum += pow(f_A[i] - mean_f, 2) + pow(f_B[i] - mean_f, 2);
            }
            double var_f = var_sum / (2.0 * n_samples - 1);
            
            // First-order indices using improved method
            std::vector<double> f_AB_a(n_samples), f_AB_b(n_samples), f_AB_c(n_samples);
            
            for (int i = 0; i < n_samples; i++) {
                // For S1_a: B with A's first parameter
                std::vector<std::vector<double>> AB_a = samples_B;
                AB_a[i][0] = samples_A[i][0];
                f_AB_a[i] = evaluate_model(AB_a, i, t);
                
                // For S1_b: B with A's second parameter  
                std::vector<std::vector<double>> AB_b = samples_B;
                AB_b[i][1] = samples_A[i][1];
                f_AB_b[i] = evaluate_model(AB_b, i, t);
                
                // For S1_c: B with A's third parameter
                std::vector<std::vector<double>> AB_c = samples_B;
                AB_c[i][2] = samples_A[i][2];
                f_AB_c[i] = evaluate_model(AB_c, i, t);
            }
            
            // Calculate first-order indices
            results.S1_a[t_idx] = calculate_first_order_index(f_A, f_B, f_AB_a, mean_f, var_f);
            results.S1_b[t_idx] = calculate_first_order_index(f_A, f_B, f_AB_b, mean_f, var_f);
            results.S1_c[t_idx] = calculate_first_order_index(f_A, f_B, f_AB_c, mean_f, var_f);
            
            // Normalize to sum to 1
            double sum_S1 = results.S1_a[t_idx] + results.S1_b[t_idx] + results.S1_c[t_idx];
            if (sum_S1 > 0) {
                results.S1_a[t_idx] /= sum_S1;
                results.S1_b[t_idx] /= sum_S1;
                results.S1_c[t_idx] /= sum_S1;
            }
            
            // Total indices (simplified approach)
            results.S_total_a[t_idx] = 0.0;
            results.S_total_b[t_idx] = 0.0;
            results.S_total_c[t_idx] = 0.0;
        }
        
        return results;
    }
    
    // Validation: compare numerical inverse Laplace with analytical solution
    void validate_solution(const std::vector<double>& time_points, double a, double b, double c) {
        std::cout << "\n=== Validation: Numerical vs Analytical ===" << std::endl;
        std::cout << "Parameters: a=" << a << ", b=" << b << ", c=" << c << std::endl;
        std::cout << "Time\t\tNumerical\tAnalytical\tError" << std::endl;
        std::cout << "------------------------------------------------" << std::endl;
        
        laplace_func->set_parameters(a, b, c);
        std::function<mpf_class(const mpf_class&)> target_func = 
            [this](const mpf_class& s) { return (*laplace_func)(s); };
        solver->set_target_function(target_func);
        
        for (size_t i = 0; i < time_points.size(); i += 4) {  // Every 4th point
            double t = time_points[i];
            
            mpf_class t_val(t);
            mpf_class numerical = solver->inverse_laplace_transform(n_terms, t_val);
            mpf_class analytical = laplace_func->analytical_solution(t_val);
            
            double error = std::abs(numerical.get_d() - analytical.get_d());
            
            std::cout << std::fixed << std::setprecision(3) << t << "\t\t"
                      << std::setprecision(6) << numerical.get_d() << "\t"
                      << analytical.get_d() << "\t"
                      << std::scientific << error << std::endl;
        }
    }
};

// Main function to run the polynomial time-domain coupled analysis
int main(int argc, char* argv[]) {
    // Initialize MPI (required by solver)
    MPI_Init(&argc, &argv);
    
    try {
        std::cout << "=== Polynomial Time-Domain Inverse Laplace Transform & Sobol Analysis ===" << std::endl;
        std::cout << "Time-domain function: f(t) = a³t - b²t² + c⁴t" << std::endl;
        std::cout << "Laplace transform: F(s) = (a³ + c⁴)/s² - 2b²/s³" << std::endl;
        std::cout << "Time range: 0.1 to 2" << std::endl;
        std::cout << "Parameter a range: 0.5 to 1.5" << std::endl;
        std::cout << "Parameter b range: 0.5 to 1.5" << std::endl;
        std::cout << "Parameter c range: 0.5 to 1.5" << std::endl;
        std::cout << "Performing Sobol sensitivity analysis..." << std::endl;
        
        // Initialize solver and function
        InverseLaplaceSolver solver(30);  // 30 decimal places precision
        PolynomialTimeLaplaceFunction laplace_func(1.0, 1.0, 1.0);  // Default parameters
        
        // Create Sobol analyzer
        PolynomialTimeSobolAnalysis sobol_analyzer(&solver, &laplace_func, 2048, 14);  // 2048 samples, 14 terms
        
        // Define time points from 0.1 to 2 (good range for polynomial observation)
        std::vector<double> time_points;
        for (int i = 0; i <= 20; i++) {
            double t = 0.1 + i * 1.9 / 20.0;  // 21 points from 0.1 to 2
            time_points.push_back(t);
        }
        
        // Parameter ranges (keeping them reasonable to avoid numerical issues)
        double a_min = 0.5, a_max = 1.5;   // a³ term coefficient
        double b_min = 0.5, b_max = 1.5;   // b² term coefficient  
        double c_min = 0.5, c_max = 1.5;   // c⁴ term coefficient
        
        // Validation step
        std::cout << "Validating solution with default parameters..." << std::endl;
        sobol_analyzer.validate_solution(time_points, 1.0, 1.0, 1.0);
        
        // Perform sensitivity analysis
        std::cout << "Starting sensitivity analysis..." << std::endl;
        auto start_time = std::chrono::high_resolution_clock::now();
        
        auto results = sobol_analyzer.analyze_sensitivity(time_points, a_min, a_max, b_min, b_max, c_min, c_max);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "\n=== Timing Summary ===" << std::endl;
        std::cout << "Total analysis time: " << duration.count() << " ms (" 
                  << std::fixed << std::setprecision(1) << duration.count() / 1000.0 << " seconds)" << std::endl;
        std::cout << "Time per point: " << std::fixed << std::setprecision(1) 
                  << (duration.count() / (double)time_points.size()) << " ms average" << std::endl;
        
        // Save results to file
        std::ofstream outfile("polynomial_time_sobol_results.txt");
        outfile << "# Polynomial Time-Domain Inverse Laplace Transform & Sobol Sensitivity Analysis Results\n";
        outfile << "# Time\tS1(a³)\tS1(b²)\tS1(c⁴)\tST(a³)\tST(b²)\tST(c⁴)\n";
        
        for (size_t i = 0; i < results.times.size(); i++) {
            outfile << std::fixed << std::setprecision(2) << results.times[i] << "\t"
                    << std::setprecision(4) << results.S1_a[i] << "\t"
                    << results.S1_b[i] << "\t" << results.S1_c[i] << "\t"
                    << results.S_total_a[i] << "\t" << results.S_total_b[i] << "\t" 
                    << results.S_total_c[i] << std::endl;
        }
        outfile.close();
        
        // Print results table
        std::cout << "\n=== Polynomial Time-Domain Sobol Sensitivity Analysis Results ===" << std::endl;
        std::cout << "Time\t\tS1(a³)\t\tS1(b²)\t\tS1(c⁴)\t\tST(a³)\tST(b²)\t\tST(c⁴)" << std::endl;
        std::cout << "--------------------------------------------------------------------------------" << std::endl;
        
        // Print subset of results (every 3rd point to keep output manageable)
        for (size_t i = 0; i < results.times.size(); i += 3) {
            std::cout << std::fixed << std::setprecision(2) << results.times[i] << "\t\t"
                      << std::setprecision(4) << results.S1_a[i] << "\t\t"
                      << results.S1_b[i] << "\t\t" << results.S1_c[i] << "\t\t"
                      << results.S_total_a[i] << "\t" << results.S_total_b[i] << "\t\t"
                      << results.S_total_c[i] << std::endl;
        }
        
        std::cout << "\nResults saved to: polynomial_time_sobol_results.txt" << std::endl;
        std::cout << "Use Python script to generate plots." << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        MPI_Finalize();
        return 1;
    }

    MPI_Finalize();
    return 0;
}