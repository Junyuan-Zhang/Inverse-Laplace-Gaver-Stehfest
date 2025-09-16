#!/usr/bin/env python3
"""
Verification Plot Script for Inverse Laplace Transform
Compares the numerical inverse Laplace solution with analytical time domain solution.
Usage: conda activate ml_env && python verification_plot.py
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import subprocess
import sys

def activate_conda_env_and_run():
    """Run the verification.cpp executable if needed"""
    verification_exe = "./verification"
    results_file = "verification_results.txt"
    
    # Check if results file exists and is recent
    if not os.path.exists(results_file):
        print("Results file not found. Please compile and run verification.cpp first:")
        print("cd /Users/junyuanzhang/Desktop/Solver/project_closed_loop")
        print("g++ -std=c++17 -I.. -I/opt/homebrew/include -L/opt/homebrew/lib -lmpfr -lgmp verification.cpp ../solver.o ../mpfr_float.o -o verification")
        print("./verification")
        return False
    
    return True

def read_verification_results(filename="verification_results.txt"):
    """Read the verification results from the output file"""
    try:
        # Read the data file
        data = []
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        # Skip header lines and parse data
        for line in lines:
            line = line.strip()
            if line and not line.startswith('#') and not line.startswith('='):
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        t = float(parts[0])
                        numerical = float(parts[1])
                        analytical = float(parts[2])
                        data.append([t, numerical, analytical])
                    except ValueError:
                        continue
        
        if not data:
            print("No valid data found in results file")
            return None
            
        # Convert to numpy array
        data = np.array(data)
        return {
            't_values': data[:, 0],
            'numerical_solution': data[:, 1], 
            'analytical_solution': data[:, 2]
        }
        
    except FileNotFoundError:
        print(f"Results file {filename} not found!")
        return None
    except Exception as e:
        print(f"Error reading results: {e}")
        return None

def create_verification_plots(results):
    """Create comprehensive verification plots"""
    if results is None:
        return
        
    t = results['t_values']
    numerical = results['numerical_solution']
    analytical = results['analytical_solution']
    
    # Calculate relative error
    relative_error = np.abs((numerical - analytical) / (analytical + 1e-15)) * 100
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Verification of Inverse Laplace Transform\nLaplace Space Solution vs Time Domain Analytical Solution', 
                 fontsize=16, fontweight='bold')
    
    # Plot 1: Solutions comparison (log-log scale)
    axes[0, 0].loglog(t, numerical, 'b-', linewidth=2, label='Numerical (Inverse Laplace)', alpha=0.8)
    axes[0, 0].loglog(t, analytical, 'r--', linewidth=2, label='Analytical (Time Domain)', alpha=0.8)
    axes[0, 0].set_xlabel('Time t [d]')
    axes[0, 0].set_ylabel('Temperature T [K]')
    axes[0, 0].set_title('Solutions Comparison (Log-Log Scale)')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Plot 2: Solutions comparison (semi-log scale)
    axes[0, 1].semilogx(t, numerical, 'b-', linewidth=2, label='Numerical (Inverse Laplace)', alpha=0.8)
    axes[0, 1].semilogx(t, analytical, 'r--', linewidth=2, label='Analytical (Time Domain)', alpha=0.8)
    axes[0, 1].set_xlabel('Time t [d]')
    axes[0, 1].set_ylabel('Temperature T [K]')
    axes[0, 1].set_title('Solutions Comparison (Semi-Log Scale)')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Plot 3: Relative error
    axes[1, 0].loglog(t, relative_error, 'g-', linewidth=2, alpha=0.8)
    axes[1, 0].set_xlabel('Time t [d]')
    axes[1, 0].set_ylabel('Relative Error [%]')
    axes[1, 0].set_title('Relative Error between Solutions')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Plot 4: Residuals
    residuals = numerical - analytical
    axes[1, 1].semilogx(t, residuals, 'purple', linewidth=2, alpha=0.8)
    axes[1, 1].axhline(y=0, color='k', linestyle='--', alpha=0.5)
    axes[1, 1].set_xlabel('Time t [d]')
    axes[1, 1].set_ylabel('Residuals (Numerical - Analytical) [K]')
    axes[1, 1].set_title('Residuals')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save the plot
    output_file = 'verification_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Verification plot saved as: {output_file}")
    
    # Show statistics
    print("\n=== Verification Statistics ===")
    print(f"Number of data points: {len(t)}")
    print(f"Time range: {t.min():.1e} to {t.max():.1e} days")
    print(f"Temperature range (numerical): {numerical.min():.3f} to {numerical.max():.3f} K")
    print(f"Temperature range (analytical): {analytical.min():.3f} to {analytical.max():.3f} K")
    print(f"Mean relative error: {np.mean(relative_error):.3f}%")
    print(f"Max relative error: {np.max(relative_error):.3f}%")
    print(f"RMSE: {np.sqrt(np.mean(residuals**2)):.6f} K")
    
    plt.show()

def create_parameter_info_plot():
    """Create a plot showing the problem parameters"""
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    # Problem parameters from require.txt
    params_text = """
    Problem Parameters:
    
    R = 0.2 m (radius)
    Q = 1 m³/d (flow rate)
    Df = 1 m²/d (fracture diffusivity)
    Db = 1×10⁻⁹ m²/d (matrix diffusivity)
    h(t) = 30 K (step function)
    x = 1000 m (distance)
    n = 300 (Gaver-Stehfest terms)
    
    Verification Approach:
    1. Laplace space solution with Bessel functions:
       T̄rf = H(p) × exp((γ₁ - √(γ₁² + 4γ₂))/2 × x)
       
    2. Time domain analytical solution:
       Trf = h(t)/2 × erfc((x - Qt/(πR²))/(2√(Dft)))
       
    3. Inverse Laplace transform of (1) should match (2)
    """
    
    ax.text(0.05, 0.95, params_text, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    ax.axis('off')
    ax.set_title('Verification Problem Setup', fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('verification_parameters.png', dpi=300, bbox_inches='tight')
    print("Parameter info plot saved as: verification_parameters.png")
    plt.show()

def main():
    """Main function"""
    print("=== Inverse Laplace Transform Verification Plotter ===")
    print("Activating ml_env environment and creating plots...")
    
    # Check if we can read the results
    if not activate_conda_env_and_run():
        return
    
    # Read verification results
    print("Reading verification results...")
    results = read_verification_results()
    
    if results is None:
        print("Could not read verification results. Exiting.")
        return
    
    # Create plots
    print("Creating verification plots...")
    create_verification_plots(results)
    
    print("Creating parameter information plot...")
    create_parameter_info_plot()
    
    print("\nVerification plots completed successfully!")
    print("Files generated:")
    print("  - verification_comparison.png")
    print("  - verification_parameters.png")

if __name__ == "__main__":
    main()