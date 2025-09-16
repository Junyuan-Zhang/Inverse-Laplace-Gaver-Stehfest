# Makefile for Verification Program
# Compares Laplace space solution with time domain analytical solution

# Compiler settings
CXX = /opt/homebrew/bin/mpic++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra
INCLUDES = -I.. -I/opt/homebrew/include -I/opt/homebrew/opt/mpfr/include
LIBS = -L/opt/homebrew/lib -L/opt/homebrew/opt/mpfr/lib -lmpfr -lgmp -lmpi -lgsl -lgslcblas

# Source files
VERIFICATION_SRC = verification.cpp
SOLVER_SRC = ../solver.cpp
MPFR_FLOAT_SRC = ../mpfr_float.cpp
SOLVER_HDR = ../solver.h
MPFR_FLOAT_HDR = ../mpfr_float.h

# Object files
VERIFICATION_OBJ = verification.o
SOLVER_OBJ = ../solver.o
MPFR_FLOAT_OBJ = ../mpfr_float.o

# Executable
VERIFICATION_EXE = verification

# Python script
PLOT_SCRIPT = verification_plot.py

# Default target
all: $(VERIFICATION_EXE)

# Build the verification executable
$(VERIFICATION_EXE): $(VERIFICATION_OBJ) $(SOLVER_OBJ) $(MPFR_FLOAT_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# Build verification object file
$(VERIFICATION_OBJ): $(VERIFICATION_SRC) $(SOLVER_HDR) $(MPFR_FLOAT_HDR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Build solver object file (if not already built)
$(SOLVER_OBJ): $(SOLVER_SRC) $(SOLVER_HDR)
	cd .. && $(CXX) $(CXXFLAGS) $(INCLUDES) -c solver.cpp -o solver.o

# Build MPFR float object file (if not already built)
$(MPFR_FLOAT_OBJ): $(MPFR_FLOAT_SRC) $(MPFR_FLOAT_HDR)
	cd .. && $(CXX) $(CXXFLAGS) $(INCLUDES) -c mpfr_float.cpp -o mpfr_float.o

# Run the verification
run: $(VERIFICATION_EXE)
	@echo "Running verification analysis..."
	./$(VERIFICATION_EXE)

# Generate plots using Python
plot: verification_results.txt
	@echo "Generating verification plots..."
	@if command -v conda >/dev/null 2>&1; then \
		echo "Activating ml_env and running Python plotting script..."; \
		conda activate ml_env && python $(PLOT_SCRIPT); \
	else \
		echo "Conda not found, running with system Python..."; \
		python3 $(PLOT_SCRIPT); \
	fi

# Complete verification workflow
verify: $(VERIFICATION_EXE)
	@echo "=== Complete Verification Workflow ==="
	@echo "1. Running inverse Laplace transform verification..."
	./$(VERIFICATION_EXE)
	@echo ""
	@echo "2. Generating comparison plots..."
	@if command -v conda >/dev/null 2>&1; then \
		conda activate ml_env && python $(PLOT_SCRIPT); \
	else \
		python3 $(PLOT_SCRIPT); \
	fi
	@echo ""
	@echo "Verification complete! Check verification_comparison.png and verification_parameters.png"

# Clean build files and results
clean:
	rm -f *.o $(VERIFICATION_EXE) verification_results.txt verification_comparison.png verification_parameters.png

# Check dependencies
check-deps:
	@echo "Checking dependencies..."
	@echo -n "MPI compiler: "
	@if command -v $(CXX) >/dev/null 2>&1; then echo "✓ Found"; else echo "✗ Missing"; fi
	@echo -n "MPFR library: "
	@if [ -f /opt/homebrew/lib/libmpfr.dylib ] || [ -f /usr/local/lib/libmpfr.so ]; then echo "✓ Found"; else echo "✗ Missing"; fi
	@echo -n "GMP library: "
	@if [ -f /opt/homebrew/lib/libgmp.dylib ] || [ -f /usr/local/lib/libgmp.so ]; then echo "✓ Found"; else echo "✗ Missing"; fi
	@echo -n "Python: "
	@if command -v python3 >/dev/null 2>&1; then echo "✓ Found"; else echo "✗ Missing"; fi
	@echo -n "Conda: "
	@if command -v conda >/dev/null 2>&1; then echo "✓ Found"; else echo "⚠ Missing (optional)"; fi

# Help
help:
	@echo "Available targets:"
	@echo "  all          - Build the verification program"
	@echo "  run          - Run the verification analysis"
	@echo "  plot         - Generate verification plots (requires Python)"
	@echo "  verify       - Complete workflow: build, run, and plot"
	@echo "  clean        - Remove build files and results"
	@echo "  check-deps   - Check required dependencies"
	@echo "  help         - Show this help message"
	@echo ""
	@echo "Example usage:"
	@echo "  make check-deps    # Verify dependencies"
	@echo "  make verify        # Complete verification"

.PHONY: all run plot verify clean check-deps help