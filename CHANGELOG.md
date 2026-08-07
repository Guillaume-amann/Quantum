# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

---

## [Unreleased]

### Added
- (Nothing yet — upcoming v0.1.0)

### Changed

### Fixed

### Deprecated

### Removed

---

## [0.0.0] – 2026-08-07

### Added

#### Core Infrastructure
- **Gate.h**: Complete single- and two-qubit gate library
  - Single-qubit: `identity`, `X`, `Y`, `Z`, `H`, `Rx(θ)`, `Ry(θ)`, `Rz(θ)`
  - Two-qubit: `CNOT`, `CZ`, `SWAP`, `Rzz(θ)`, `controlled_U()`, `Toffoli`
  - Register embedding: `expand()` for single qubits, `expand_two()` for arbitrary qubit pairs
  - Tensor product and sequential composition operators
  - **Convention:** Big-endian qubit ordering (qubit 0 = MSB)

- **Qbit.h**: Pure state-vector simulator
  - State vector ψ representation (2^n complex amplitudes)
  - Unitary evolution via gate application
  - Measurement with probabilistic collapse (Born rule)
  - Thread-local RNG per thread (MPI-safe)
  - Normalisation enforcement

- **DensityMatrix.h**: Open-system simulator with full density-matrix formalism
  - Density-matrix ρ representation (2^n × 2^n Hermitian, trace 1)
  - Unitary evolution: `ρ → UρU†`
  - Kraus-channel noise evolution: `ρ → Σ_k K_k ρ K_k†`
  - Partial trace for reduced density matrices
  - Single-qubit and full-register measurement with collapse
  - Observables: expectation value, purity, von Neumann entropy, fidelity
  - Probability vector extraction (non-destructive)

- **NoiseModel.h**: Kraus-operator noise channel factory
  - **Single-qubit channels**: bit-flip, phase-flip, bit-phase-flip, depolarising
  - **Amplitude damping**: T1 relaxation via `amplitude_damping(T1, t_gate, target, n)`
  - **Generalised amplitude damping**: finite-temperature T1 (thermal excitation)
  - **T1/T2 combined**: physical constraint T2 ≤ 2·T1, pure-dephasing extraction
  - **Two-qubit channels**: depolarising (16 Pauli-error Kraus ops)
  - **Measurement error**: classical readout confusion matrix
  - **Completeness verification**: Σ K†K = I checked at runtime
  - **All channels**: fully parameterised, gate-duration aware

#### Applications
- **QuantumSim.cpp**: MPI QAOA grid search for 2-qubit QUBO
  - Depth-1 QAOA ansatz (preparation + cost + mixer)
  - Parallel (γ, α) grid scan (each rank owns γ-stripe)
  - Non-destructive sampling via repeated state preparation
  - Energy surface + measurement histogram outputs
  - Timing breakdown (grid search / measurement phases)

- **Electricity_procurement.cpp**: 8-qubit bin-packing QUBO
  - Programmatic QUBO construction from band catalogue + hourly demand
  - Closed-form Ising Hamiltonian (linear `h_i` + quadratic `J_ij`)
  - Normalisation by max coefficient (QAOA angle scaling)
  - Genuine entanglement via `Rzz` gates (time-overlapping bands)
  - Optional budget constraint (C4, dense coupling)
  - Classical brute-force validation (all 2^8 = 256 states)
  - Feasibility checking (coverage ≥ demand per hour, budget)
  - Non-destructive sampling (CDF-based binary search)
  - Per-rank noise sweep (100K shots distributed across ranks)

- **Electricity_procurement_noise.cpp**: Noisy version on DensityMatrix
  - Same instance, fixed optimal (γ, α) from noiseless run
  - Per-gate Kraus noise: 1-qubit depolarising on H/Rz/Rx, 2-qubit on Rzz
  - MPI-parallelised noise-level sweep (alternative to angle grid)
  - Output: (p, energy, purity, fidelity, entropy, top_cost, top_feasible)
  - Demonstrations: partial_trace (entanglement decay), partial_measurement (qubit-by-qubit readout), measurement_error (classical confusion)

#### Testing & Validation
- CI/CD pipeline (.github/workflows/ci.yml)
  - Build matrix: GCC, Clang on Ubuntu
  - Automated tests on every push/PR
  - Clang-format enforcement

#### Documentation
- README.md: Project overview, quick start, performance benchmarks
- CONTRIBUTING.md: Branch workflow, commit conventions, testing guidelines
- CHANGELOG.md: This file
- docs/ARCHITECTURE.md: Design decisions, physics conventions, coupling graph
- docs/FORMULATION.md: Electricity QUBO derivation, Ising mapping, constraint algebra
- docs/NOISE_MODELS.md: Kraus operator interpretation, physical parameters
- docs/API.md: Class reference (auto-generated friendly)
- examples/README.md: Runnable demos

#### Scripts
- scripts/energy_surface.py: 3D visualisation of QAOA grid
- scripts/measurement_histogram.py: Bar chart of sampled outcomes

### Technical Highlights

- **Mathematical rigour**: All density matrices enforced Hermitian + Tr(ρ)=1; all Kraus channels verified complete (Σ K†K = I)
- **Physical fidelity**: T2 ≤ 2·T1 constraint enforced; depolarising p ∈ [0,0.75]; amplitude damping γ ∈ [0,1]
- **Symmetry**: Two simulation back-ends (Qbit, DensityMatrix) coexist without interference; noise is a decorator, not a replacement
- **Scalability**: MPI parallelisation on grid search + noise sweep; 8-qubit density matrix on 4 cores ≈ 12 min noiseless, 45 min noisy
- **Sampling discipline**: Non-destructive state preparation per sample (fixes QuantumSim.cpp v0 bug where collapsed state was reused)

### Known Limitations

- Density matrix: O(2^(3n)) memory and time; max ~14 qubits on standard workstation
- QAOA: Depth 1 only (no adaptive/warm-start angle refinement)
- QSVM: Kernel matrix computation untested; dual QP solver not yet integrated
- Noise simulation: Gate duration not extracted from gate definitions; must pass t_gate manually

### Not Included (Intentional)

- Annealing (hybrid classical/quantum algorithms)
- Cat qubits or other platform-specific models
- GPU acceleration (roadmap for v1.x)
- Hardware integration (Qiskit/Cirq connectors)

---

## Notes

- **v0.0.0** is a **stable baseline**, not yet feature-complete. Core simulation and QAOA are production-ready; QSVM and QEC are prototypes.
- All code reviewed for mathematical correctness and physics fidelity.
- Ready for research use; educational clarity prioritised over micro-optimisation.
