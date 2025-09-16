import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

def load_and_analyze_data(filename):
    """Load and analyze polynomial time-domain Sobol results"""
    try:
        # Read the data file
        data = pd.read_csv(filename, sep='\t', comment='#', 
                          names=['Time', 'S1_a', 'S1_b', 'S1_c', 'ST_a', 'ST_b', 'ST_c'])
        
        print("Polynomial Time-Domain Inverse Laplace Transform & Sobol Sensitivity Analysis - Visualization")
        print("=" * 100)
        print(f"Loaded {len(data)} data points from {filename}")
        print(f"Time range: {data['Time'].min():.2f} to {data['Time'].max():.2f}")
        
        return data
    except Exception as e:
        print(f"Error loading data: {e}")
        return None

def print_summary_statistics(data):
    """Print summary statistics for the sensitivity analysis"""
    print("\n=== Summary Statistics ===")
    
    # First-order indices statistics
    for param in ['a', 'b', 'c']:
        col = f'S1_{param}'
        mean_val = data[col].mean()
        std_val = data[col].std()
        print(f"Parameter {param}³ - Mean S1: {mean_val:.4f}, Std: {std_val:.4f}")
    
    # Total indices statistics
    for param in ['a', 'b', 'c']:
        col = f'ST_{param}'
        mean_val = data[col].mean()
        std_val = data[col].std()
        print(f"Parameter {param}³ - Mean ST: {mean_val:.4f}, Std: {std_val:.4f}")
    
    # Check sum of first-order indices
    data['S1_sum'] = data['S1_a'] + data['S1_b'] + data['S1_c']
    print(f"\nSum of first-order indices - Mean: {data['S1_sum'].mean():.4f}, Max: {data['S1_sum'].max():.4f}")
    if data['S1_sum'].max() > 1.1:
        print("Warning: Sum of first-order indices exceeds 1, indicating potential numerical issues")
    
    # Find maximum sensitivities
    for param in ['a', 'b', 'c']:
        col = f'S1_{param}'
        max_idx = data[col].idxmax()
        max_val = data[col].max()
        max_time = data.loc[max_idx, 'Time']
        print(f"\nMaximum sensitivity for parameter {param}³: {max_val:.4f} at time {max_time:.2f}")

def create_comprehensive_plot(data):
    """Create a comprehensive 4-panel plot showing all aspects of the analysis"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Polynomial Time-Domain f(t) = a³t - b²t² + c⁴t\nInverse Laplace Transform & Sobol Sensitivity Analysis', 
                 fontsize=14, fontweight='bold')
    
    # Colors for each parameter
    colors = {'a': '#FF6B6B', 'b': '#4ECDC4', 'c': '#45B7D1'}
    
    # Panel 1: First-order sensitivity indices
    ax1 = axes[0, 0]
    ax1.plot(data['Time'], data['S1_a'], 'o-', color=colors['a'], linewidth=2, markersize=4, label='S₁(a³)')
    ax1.plot(data['Time'], data['S1_b'], 's-', color=colors['b'], linewidth=2, markersize=4, label='S₁(b²)')
    ax1.plot(data['Time'], data['S1_c'], '^-', color=colors['c'], linewidth=2, markersize=4, label='S₁(c⁴)')
    ax1.set_xlabel('Time t')
    ax1.set_ylabel('First-Order Sensitivity Index')
    ax1.set_title('First-Order Sobol Indices')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 1)
    
    # Panel 2: Total sensitivity indices
    ax2 = axes[0, 1]
    ax2.plot(data['Time'], data['ST_a'], 'o-', color=colors['a'], linewidth=2, markersize=4, label='Sₜ(a³)')
    ax2.plot(data['Time'], data['ST_b'], 's-', color=colors['b'], linewidth=2, markersize=4, label='Sₜ(b²)')
    ax2.plot(data['Time'], data['ST_c'], '^-', color=colors['c'], linewidth=2, markersize=4, label='Sₜ(c⁴)')
    ax2.set_xlabel('Time t')
    ax2.set_ylabel('Total Sensitivity Index')
    ax2.set_title('Total Sobol Indices')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, 1)
    
    # Panel 3: Comparison of first-order vs total for each parameter
    ax3 = axes[1, 0]
    width = 0.25
    x = np.arange(len(data))
    
    ax3.bar(x - width, data['S1_a'], width, label='S₁(a³)', color=colors['a'], alpha=0.7)
    ax3.bar(x, data['S1_b'], width, label='S₁(b²)', color=colors['b'], alpha=0.7)
    ax3.bar(x + width, data['S1_c'], width, label='S₁(c⁴)', color=colors['c'], alpha=0.7)
    
    # Only show every 4th time label to avoid overcrowding
    time_labels = [f"{t:.1f}" if i % 4 == 0 else "" for i, t in enumerate(data['Time'])]
    ax3.set_xticks(x)
    ax3.set_xticklabels(time_labels, rotation=45)
    ax3.set_xlabel('Time t')
    ax3.set_ylabel('First-Order Sensitivity')
    ax3.set_title('First-Order Sensitivity Distribution')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Panel 4: Interaction effects (Total - First-order)
    ax4 = axes[1, 1]
    interaction_a = data['ST_a'] - data['S1_a']
    interaction_b = data['ST_b'] - data['S1_b']
    interaction_c = data['ST_c'] - data['S1_c']
    
    ax4.plot(data['Time'], interaction_a, 'o-', color=colors['a'], linewidth=2, markersize=4, label='Interaction(a³)')
    ax4.plot(data['Time'], interaction_b, 's-', color=colors['b'], linewidth=2, markersize=4, label='Interaction(b²)')
    ax4.plot(data['Time'], interaction_c, '^-', color=colors['c'], linewidth=2, markersize=4, label='Interaction(c⁴)')
    ax4.set_xlabel('Time t')
    ax4.set_ylabel('Interaction Effect (Sₜ - S₁)')
    ax4.set_title('Parameter Interaction Effects')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('polynomial_time_sobol_analysis.png', dpi=300, bbox_inches='tight')
    print(f"Comprehensive analysis plot saved as: polynomial_time_sobol_analysis.png")
    
    return fig

def create_first_order_plot(data):
    """Create the main requested plot: first-order Sobol indices vs time"""
    plt.figure(figsize=(12, 8))
    
    # Colors and styles for each parameter
    colors = {'a': '#FF6B6B', 'b': '#4ECDC4', 'c': '#45B7D1', 'sum': '#9B59B6'}
    markers = {'a': 'o', 'b': 's', 'c': '^', 'sum': 'D'}
    
    # Calculate sum of first-order indices
    data['S1_sum'] = data['S1_a'] + data['S1_b'] + data['S1_c']
    
    # Plot first-order sensitivity indices
    plt.plot(data['Time'], data['S1_a'], marker=markers['a'], color=colors['a'], 
             linewidth=3, markersize=8, label='S₁(a³) - Linear term coeff.', markerfacecolor='white', markeredgewidth=2)
    plt.plot(data['Time'], data['S1_b'], marker=markers['b'], color=colors['b'], 
             linewidth=3, markersize=8, label='S₁(b²) - Quadratic term coeff.', markerfacecolor='white', markeredgewidth=2)
    plt.plot(data['Time'], data['S1_c'], marker=markers['c'], color=colors['c'], 
             linewidth=3, markersize=8, label='S₁(c⁴) - Additional linear coeff.', markerfacecolor='white', markeredgewidth=2)
    
    # Plot sum of first-order indices
    plt.plot(data['Time'], data['S1_sum'], marker=markers['sum'], color=colors['sum'], 
             linewidth=3, markersize=8, label='∑S₁ - Sum of all first-order indices', 
             markerfacecolor='white', markeredgewidth=2, linestyle='--', alpha=0.8)
    
    plt.xlabel('Time t', fontsize=14, fontweight='bold')
    plt.ylabel('First-Order Sobol Sensitivity Index', fontsize=14, fontweight='bold')
    plt.title('First-Order Sobol Indices for Polynomial Time Function\nf(t) = a³t - b²t² + c⁴t', 
              fontsize=16, fontweight='bold', pad=20)
    
    plt.legend(fontsize=12, loc='best', frameon=True, fancybox=True, shadow=True)
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.ylim(0, 1.05)
    
    # Add annotations for interesting points
    max_a_idx = data['S1_a'].idxmax()
    max_b_idx = data['S1_b'].idxmax()
    max_c_idx = data['S1_c'].idxmax()
    max_sum_idx = data['S1_sum'].idxmax()
    
    plt.annotate(f'Max a³: {data["S1_a"].max():.3f}', 
                xy=(data.loc[max_a_idx, 'Time'], data.loc[max_a_idx, 'S1_a']),
                xytext=(10, 10), textcoords='offset points',
                bbox=dict(boxstyle='round,pad=0.3', facecolor=colors['a'], alpha=0.7),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    
    plt.annotate(f'Max ∑S₁: {data["S1_sum"].max():.3f}', 
                xy=(data.loc[max_sum_idx, 'Time'], data.loc[max_sum_idx, 'S1_sum']),
                xytext=(-10, -20), textcoords='offset points',
                bbox=dict(boxstyle='round,pad=0.3', facecolor=colors['sum'], alpha=0.7),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    
    # Add horizontal line at y=1 for reference
    plt.axhline(y=1.0, color='red', linestyle=':', alpha=0.7, linewidth=2, label='Perfect additivity (∑S₁ = 1)')
    
    # Add text box with summary statistics
    textstr = f'Mean ∑S₁: {data["S1_sum"].mean():.3f}\nMax ∑S₁: {data["S1_sum"].max():.3f}\nMin ∑S₁: {data["S1_sum"].min():.3f}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.savefig('polynomial_time_first_order_sobol_indices.png', dpi=300, bbox_inches='tight')
    print(f"First-order Sobol indices plot saved as: polynomial_time_first_order_sobol_indices.png")

def create_3d_visualization(data):
    """Create 3D surface plots showing sensitivity evolution"""
    try:
        fig = plt.figure(figsize=(16, 12))
        
        # Create time meshgrid for 3D plotting
        t = data['Time'].values
        params = ['a', 'b', 'c']
        param_labels = ['a³', 'b²', 'c⁴']
        
        for i, (param, label) in enumerate(zip(params, param_labels)):
            # First-order sensitivities
            ax1 = fig.add_subplot(2, 3, i+1, projection='3d')
            
            T, P = np.meshgrid(t, [0, 1])
            S1 = np.vstack([data[f'S1_{param}'].values, data[f'S1_{param}'].values])
            
            ax1.plot_surface(T, P, S1, alpha=0.7, cmap='viridis')
            ax1.plot(t, np.zeros_like(t), data[f'S1_{param}'].values, 'r-', linewidth=3)
            ax1.set_xlabel('Time t')
            ax1.set_ylabel('Parameter Range')
            ax1.set_zlabel(f'S₁({label})')
            ax1.set_title(f'First-Order Sensitivity: {label}')
            
            # Total sensitivities
            ax2 = fig.add_subplot(2, 3, i+4, projection='3d')
            
            ST = np.vstack([data[f'ST_{param}'].values, data[f'ST_{param}'].values])
            
            ax2.plot_surface(T, P, ST, alpha=0.7, cmap='plasma')
            ax2.plot(t, np.zeros_like(t), data[f'ST_{param}'].values, 'r-', linewidth=3)
            ax2.set_xlabel('Time t')
            ax2.set_ylabel('Parameter Range')
            ax2.set_zlabel(f'Sₜ({label})')
            ax2.set_title(f'Total Sensitivity: {label}')
        
        plt.suptitle('3D Sensitivity Analysis for Polynomial Time Function f(t) = a³t - b²t² + c⁴t', 
                     fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig('polynomial_time_3d_sobol_surfaces.png', dpi=300, bbox_inches='tight')
        print(f"3D surface plots saved as: polynomial_time_3d_sobol_surfaces.png")
        
    except Exception as e:
        print(f"Could not create 3D visualization: {e}")

def main():
    # Load and analyze data
    data = load_and_analyze_data('polynomial_time_sobol_results.txt')
    if data is None:
        return
    
    # Print summary statistics
    print_summary_statistics(data)
    
    # Create plots
    print(f"\nGenerating visualization plots...")
    
    # Main comprehensive plot
    create_comprehensive_plot(data)
    
    # First-order Sobol indices plot (main requested plot)
    create_first_order_plot(data)
    
    # 3D visualization
    print(f"\nGenerating additional 3D visualization...")
    create_3d_visualization(data)
    
    print(f"\nAll plots have been generated and saved:")
    print(f"- polynomial_time_sobol_analysis.png: Comprehensive 4-panel analysis")
    print(f"- polynomial_time_first_order_sobol_indices.png: Main requested plot (first-order vs time)")
    print(f"- polynomial_time_3d_sobol_surfaces.png: 3D surface plots (if available)")

if __name__ == "__main__":
    main()