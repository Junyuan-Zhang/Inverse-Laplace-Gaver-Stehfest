# Makefile for Inverse Laplace Transform Solver with MPI and GMP

CXX = mpicxx
CXXFLAGS = -std=c++17 -Wall -Wextra -O2
GMP_INCLUDE = /opt/homebrew/Cellar/gmp/6.3.0/include
GMP_LIB = /opt/homebrew/Cellar/gmp/6.3.0/lib
LIBS = -lgmp -lgmpxx

# Include and library flags
INCLUDES = -I$(GMP_INCLUDE)
LDFLAGS = -L$(GMP_LIB) $(LIBS)

# Target executable
TARGET = inverse_laplace_solver
OBJECTS = main.o solver.o custom_function.o

.PHONY: all clean run run-single run-parallel help

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS) -o $@

main.o: main.cpp solver.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

solver.o: solver.cpp solver.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

custom_function.o: custom_function.cpp custom_function.h solver.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

run: $(TARGET)
	mpirun -np 4 ./$(TARGET)

run-single: $(TARGET)
	mpirun -np 1 ./$(TARGET)

clean:
	rm -f $(TARGET) $(OBJECTS)

help:
	@echo "Available targets:"
	@echo "  all           - Build the solver"
	@echo "  run           - Run with 4 MPI processes"
	@echo "  run-single    - Run with 1 process"
	@echo "  clean         - Remove all build files"
	@echo "  help          - Show this help message"
