# QuantumSim

A C++ implementation of a quantum system simulator using the **Quantum Approximate Optimization Algorithm (QAOA)**.  
This project supports parallel execution with **MPI** to distribute sampling across multiple CPUs.

## Features
- Quantum state representation with a custom `Qbit` class.
- Gate operations (Hadamard, Pauli-X, etc.) defined in `Gate.h`.
- Energy evaluation for a simple Hamiltonian.
- QAOA loop with parameter sweeps (`gamma`, `alpha`).
- Parallel sampling with **MPI** to leverage multiple CPUs.

## Requirements
- C++17 or newer  
- [MPI library](https://www.open-mpi.org/) (e.g., OpenMPI or MPICH)  

## Compilation
Use `mpic++` (MPI C++ compiler wrapper) to build:

```bash
mpic++ -std=c++17 -O2 QuantumSim.cpp -o QuantumSim
```

## Execution
Run on N CPUs using mpirun or mpiexec:
```bash
mpirun -np 4 ./QuantumSim
```

This launches the simulation across 4 processes. Each process computes a subset of the QAOA samples, and results are gathered at the root rank.

## File Structure
QuantumSim.cpp      # Main program
Gate.h              # Quantum gate definitions
Qbit.h              # Qubit state representation

## Example Workflow
```bash
# Compile
mpic++ -std=c++17 -O2 QuantumSim.cpp -o QuantumSim

# Run with 2 CPUs
mpirun -np 2 ./QuantumSim

# Run with 8 CPUs and save results
mpirun -np 8 ./QuantumSim > results.txt
```