#include "../solver.h"
#include <cmath>
#include <iostream>
#include <mpi.h>

// Simple test function that should work
mpf_class simple_test_function(const mpf_class& s) {
    return mpf_class(1.0) / (s + mpf_class(1.0));
}

mpf_class simple_analytical(const mpf_class& t) {
    // exp(-t) approximation
    double t_val = t.get_d();
    return mpf_class(exp(-t_val));
}

int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (rank == 0) {
        std::cout << "=== Simple Verification Test ===" << std::endl;
        std::cout << "Testing basic inverse Laplace transform" << std::endl;
    }
    
    try {
        // Initialize solver
        InverseLaplaceSolver solver(50); // Lower precision to avoid issues
        
        // Set simple functions
        solver.set_target_function(simple_test_function);
        solver.set_analytical_solution(simple_analytical);
        
        // Generate a few t values
        std::vector<mpf_class> t_values;
        for (int i = 0; i < 5; i++) {
            t_values.push_back(mpf_class(0.1 * (i + 1)));
        }
        
        if (rank == 0) {
            std::cout << "Generated " << t_values.size() << " test values" << std::endl;
            std::cout << "Running inverse Laplace transform..." << std::endl;
        }
        
        // Solve with fewer terms
        SolverResults results = solver.solve(50, t_values);
        
        if (rank == 0) {
            std::cout << "Test completed successfully!" << std::endl;
            std::cout << "RMSE Error: " << results.rmse_error << std::endl;
            
            // Write simple results
            std::ofstream file("verification_results.txt");
            file << "# t_value numerical analytical\n";
            for (size_t i = 0; i < results.t_values.size(); i++) {
                file << results.t_values[i].get_d() << " " 
                     << results.numerical_results[i].get_d() << " "
                     << results.analytical_results[i].get_d() << "\n";
            }
            file.close();
            std::cout << "Results written to verification_results.txt" << std::endl;
        }
        
    } catch (const std::exception& e) {
        if (rank == 0) {
            std::cout << "Error: " << e.what() << std::endl;
        }
    }
    
    // Finalize MPI
    MPI_Finalize();
    return 0;
}