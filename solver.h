#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <vector>
#include <functional>
#include <string>
#include <cmath>
#include <gmp.h>
#include <gmpxx.h>

struct SolverResults {
    std::vector<mpf_class> t_values;
    std::vector<mpf_class> numerical_results;
    std::vector<mpf_class> analytical_results;
    double rmse_error = -1.0;
};

class InverseLaplaceSolver {
private:
    int precision;
    int mpi_rank;
    int mpi_size;
    bool has_analytical;
    
    std::function<mpf_class(const mpf_class&)> target_function;
    std::function<mpf_class(const mpf_class&)> analytical_function;
    
    // Taylor series approximation for exp(x)
    mpf_class exp_approx(const mpf_class& x, int terms = 50) {
        mpf_class result(1.0);
        mpf_class term(1.0);
        
        for (int i = 1; i <= terms; i++) {
            term = term * x / mpf_class(i);
            result += term;
        }
        return result;
    }
    
public:
    InverseLaplaceSolver(int prec = 100);
    ~InverseLaplaceSolver();
    
    // Set the target function F(s) to be inverted
    void set_target_function(std::function<mpf_class(const mpf_class&)> func);
    
    // Set analytical solution for verification (optional)
    void set_analytical_solution(std::function<mpf_class(const mpf_class&)> func);
    
    // Calculate Vi(n, i) coefficient function
    mpf_class Vi(int n, int i);
    
    // Main inverse Laplace transform calculation with MPI support
    mpf_class inverse_laplace_transform(int n, const mpf_class& t);
    
    // Generate logarithmic range of t values
    std::vector<mpf_class> generate_t_values(double log_start = -4.0, double log_end = 2.0, double step = 0.1);
    
    // Solve the inverse Laplace transform for multiple t values
    SolverResults solve(int n, const std::vector<mpf_class>& t_values);
    
    // Calculate RMSE error between numerical and analytical solutions
    double calculate_rmse(const std::vector<mpf_class>& numerical, const std::vector<mpf_class>& analytical);
    
    // Print high-precision floating point number
    void print_mpf(const mpf_class& value, int decimal_places = 10);
    
    // Print results in formatted table
    void print_results(const SolverResults& results, int num_samples = 10);
    
    // Get MPI rank and size
    int get_mpi_rank() const { return mpi_rank; }
    int get_mpi_size() const { return mpi_size; }
};

// Predefined example functions
namespace ExampleFunctions {
    // Example target function: (5*s^2 - 15*s + 7) / ((s+1)*(s-2)^3)
    mpf_class default_target_function(const mpf_class& s);
    
    // Corresponding analytical solution: -exp(-t) - 0.5*t^2*exp(2*t) + 2*t*exp(2*t) + exp(2*t)
    mpf_class default_analytical_solution(const mpf_class& t);
    
    // Example: Simple exponential decay F(s) = 1/(s+a), f(t) = exp(-a*t)
    std::function<mpf_class(const mpf_class&)> exponential_decay_target(double a);
    std::function<mpf_class(const mpf_class&)> exponential_decay_analytical(double a);
    
    // Example: Step function F(s) = 1/s, f(t) = 1
    mpf_class step_function_target(const mpf_class& s);
    mpf_class step_function_analytical(const mpf_class& t);
}

#endif // SOLVER_H
