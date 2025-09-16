# Inverse Laplace Transform with Sobol Sensitivity Analysis
# Main Makefile - delegates to polynomial example

# Default target - build polynomial analysis
all:
	@echo "Building Inverse Laplace Transform with Sobol Sensitivity Analysis..."
	@echo "Using polynomial time-domain example: f(t) = a³t - b²t² + c⁴t"
	$(MAKE) -f Makefile_polynomial

# Run the polynomial analysis
run:
	$(MAKE) -f Makefile_polynomial run

# Generate plots
plot:
	$(MAKE) -f Makefile_polynomial plot

# Complete workflow
full-analysis:
	$(MAKE) -f Makefile_polynomial full-analysis

# Clean build files
clean:
	$(MAKE) -f Makefile_polynomial clean

# Check dependencies
check-deps:
	$(MAKE) -f Makefile_polynomial check-deps

# Help
help:
	@echo "Available targets:"
	@echo "  all          - Build the polynomial analysis program"
	@echo "  run          - Run the analysis"
	@echo "  plot         - Generate plots from results"
	@echo "  full-analysis- Complete workflow: build, run, and plot"
	@echo "  clean        - Remove build files and results"
	@echo "  check-deps   - Check required dependencies"
	@echo "  help         - Show this help message"
	@echo ""
	@echo "Example workflow:"
	@echo "  make check-deps    # Verify dependencies"
	@echo "  make full-analysis # Complete analysis"

.PHONY: all run plot full-analysis clean check-deps help