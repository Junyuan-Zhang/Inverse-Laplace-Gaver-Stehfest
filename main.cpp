#include "solver.h"
#include "custom_function.h"
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>

void print_banner() {
    std::cout << "=====================================================" << std::endl;
    std::cout << "    Inverse Laplace Transform Solver with MPI      " << std::endl;
    std::cout << "         High-Precision GMP Implementation          " << std::endl;
    std::cout << "=====================================================" << std::endl;
}

void print_help() {
    std::cout << "\n=== Available Functions ===" << std::endl;
    std::cout << "1. Default Example: (5*s^2 - 15*s + 7) / ((s+1)*(s-2)^3)" << std::endl;
    std::cout << "2. Exponential Decay: 1/(s+a) -> exp(-a*t)" << std::endl;
    std::cout << "3. Step Function: 1/s -> 1" << std::endl;
    std::cout << "4. Custom Function (defined in custom_function.cpp)" << std::endl;
    std::cout << "\n=== Parameters ===" << std::endl;
    std::cout << "- n: Number of terms in summation (higher = more accurate, slower)" << std::endl;
    std::cout << "- t_range: Range of t values to compute (log scale)" << std::endl;
    std::cout << "- precision: Decimal precision for GMP calculations" << std::endl;
}

int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    try {
        if (rank == 0) {
            print_banner();
            std::cout << "Running on " << size << " MPI processes" << std::endl;
        }
        
        // Get user input (only from master process)
        int function_choice = 1;
        int n = 100;
        int precision = 100;
        double log_start = -4.0, log_end = 2.0, log_step = 0.1;
        bool use_analytical = true;
        
        if (rank == 0) {
            print_help();
            
            std::cout << "\n=== Configuration ===" << std::endl;
            
            std::cout << "Select function type (1-4): ";
            std::cin >> function_choice;
            
            std::cout << "Enter number of terms (n) [default: 100]: ";
            std::cin >> n;
            
            std::cout << "Enter precision (decimal places) [default: 100]: ";
            std::cin >> precision;
            
            std::cout << "Enter t-value range:" << std::endl;
            std::cout << "  Log start [default: -4]: ";
            std::cin >> log_start;
            std::cout << "  Log end [default: 2]: ";
            std::cin >> log_end;
            std::cout << "  Step size [default: 0.1]: ";
            std::cin >> log_step;
            
            std::cout << "Use analytical solution for verification? (1=yes, 0=no): ";
            std::cin >> use_analytical;
            
            std::cin.ignore(); // Clear input buffer
        }
        
        // Broadcast parameters to all processes
        MPI_Bcast(&function_choice, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&precision, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&log_start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&log_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&log_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&use_analytical, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        
        // Create solver
        InverseLaplaceSolver solver(precision);
        
        // Set up functions based on user choice
        switch (function_choice) {
            case 1: // Default example
                solver.set_target_function(ExampleFunctions::default_target_function);
                if (use_analytical) {
                    solver.set_analytical_solution(ExampleFunctions::default_analytical_solution);
                }
                if (rank == 0) {
                    std::cout << "Using default function: (5*s^2 - 15*s + 7) / ((s+1)*(s-2)^3)" << std::endl;
                }
                break;
                
            case 2: // Exponential decay
                {
                    double a = 1.0;
                    if (rank == 0) {
                        std::cout << "Enter decay constant a for 1/(s+a): ";
                        std::cin >> a;
                    }
                    MPI_Bcast(&a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    
                    solver.set_target_function(ExampleFunctions::exponential_decay_target(a));
                    if (use_analytical) {
                        solver.set_analytical_solution(ExampleFunctions::exponential_decay_analytical(a));
                    }
                    if (rank == 0) {
                        std::cout << "Using exponential decay: 1/(s+" << a << ")" << std::endl;
                    }
                }
                break;
                
            case 3: // Step function
                solver.set_target_function(ExampleFunctions::step_function_target);
                if (use_analytical) {
                    solver.set_analytical_solution(ExampleFunctions::step_function_analytical);
                }
                if (rank == 0) {
                    std::cout << "Using step function: 1/s" << std::endl;
                }
                break;
                
            case 4: // Custom function
                solver.set_target_function(custom_target_function);
                if (HAS_ANALYTICAL_SOLUTION && use_analytical) {
                    solver.set_analytical_solution(custom_analytical_function);
                }
                if (rank == 0) {
                    std::cout << "Using custom function from custom_function.cpp" << std::endl;
                    if (HAS_ANALYTICAL_SOLUTION) {
                        std::cout << "Analytical solution: Available for verification" << std::endl;
                    } else {
                        std::cout << "Analytical solution: Not provided" << std::endl;
                    }
                }
                break;
                
            default:
                if (rank == 0) {
                    std::cout << "Invalid choice, using default function." << std::endl;
                }
                solver.set_target_function(ExampleFunctions::default_target_function);
                if (use_analytical) {
                    solver.set_analytical_solution(ExampleFunctions::default_analytical_solution);
                }
                break;
        }
        
        // Generate t values
        std::vector<mpf_class> t_values = solver.generate_t_values(log_start, log_end, log_step);
        
        if (rank == 0) {
            std::cout << "\n=== Starting Calculation ===" << std::endl;
            std::cout << "Parameters:" << std::endl;
            std::cout << "  Terms (n): " << n << std::endl;
            std::cout << "  Precision: " << precision << " decimal places" << std::endl;
            std::cout << "  T-values: " << t_values.size() << " points from 10^" << log_start 
                      << " to 10^" << log_end << std::endl;
            std::cout << "  MPI processes: " << size << std::endl;
        }
        
        // Solve
        auto start_time = MPI_Wtime();
        SolverResults results = solver.solve(n, t_values);
        auto end_time = MPI_Wtime();
        
        // Print results (only master process)
        if (rank == 0) {
            std::cout << "\n=== Calculation Complete ===" << std::endl;
            std::cout << "Total time: " << std::fixed << std::setprecision(3) 
                      << (end_time - start_time) << " seconds" << std::endl;
            
            solver.print_results(results, 15);
            
            // Ask user if they want to save results to a file
            std::cout << "\nDo you want to save results to a text file? (y/n): ";
            char save_choice;
            std::cin >> save_choice;
            
            if (save_choice == 'y' || save_choice == 'Y') {
                // Generate default filename with timestamp
                std::string filename = "inverse_laplace_results.txt";
                std::cout << "Saving results to: " << filename << std::endl;
                solver.write_results_to_file(results, filename);
            }
            
            std::cout << "\n=== Summary ===" << std::endl;
            std::cout << "Successfully computed " << t_values.size() << " points" << std::endl;
            std::cout << "Used " << n << " terms in the summation" << std::endl;
            std::cout << "Precision: " << precision << " decimal places" << std::endl;
            
            // Performance metrics
            std::cout << "\n=== Performance Metrics ===" << std::endl;
            std::cout << "MPI processes used: " << results.mpi_processes_used << std::endl;
            std::cout << "Calculation time: " << std::fixed << std::setprecision(4) 
                      << results.calculation_time << " seconds" << std::endl;
            std::cout << "Time per point: " << std::fixed << std::setprecision(6) 
                      << results.time_per_point << " seconds" << std::endl;
            std::cout << "Processing rate: " << std::fixed << std::setprecision(2) 
                      << (1.0 / results.time_per_point) << " points/second" << std::endl;
            
            // MPI efficiency suggestion
            if (results.mpi_processes_used > 1) {
                std::cout << "Note: Currently using master-only calculation for precision." << std::endl;
                std::cout << "      Try different process counts to find optimal performance." << std::endl;
            } else {
                std::cout << "Tip: Try using multiple MPI processes: mpirun -np 2 or -np 4" << std::endl;
            }
            
            if (results.rmse_error > 0) {
                std::cout << "\nRMSE Error vs Analytical: " << std::scientific 
                          << results.rmse_error << std::endl;
            }
        }
        
    } catch (const std::exception& e) {
        if (rank == 0) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
        MPI_Finalize();
        return 1;
    }
    
    MPI_Finalize();
    return 0;
}
