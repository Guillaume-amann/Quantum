# Quantum — a from-scratch gate-model quantum simulator

A small C++ quantum-computing simulator built on [Eigen](https://eigen.tuxfamily.org)
and MPI, written to support the accompanying paper (`Quantum.pdf`). 
Implements the gate model (not annealing, not cat qubits) with a pure state-vector backend (`Qbit`) and a density-matrix backend (`DensityMatrix`) for noise and mixed-state simulation.
It uses it to explore QAOA for combinatorial optimisation and an electricity-procurement QUBO, with a roadmap toward open-system (noisy) simulation and quantum machine learning (QSVM) for anomaly detection.

The codebase has two simulation back-ends that sit side by side:

- **`Qbit`** — a pure state-vector simulator for the ideal, noise-free case.
- **`DensityMatrix`** — a density-matrix simulator for mixed states and noise, kept as
  a parallel class rather than a replacement for `Qbit`.

## Features

### Core
- **Gate library**: 1- and 2-qubit gates (Pauli, Hadamard, rotation families); arbitrary-qubit embedding via `expand_two()`
- **Pure state simulation** (`Qbit`): state-vector backend, fast ideal-case QAOA
- **Open-system simulation** (`DensityMatrix`): density-matrix formalism, Kraus-channel noise evolution, partial trace, measurement collapse, entropy/purity/fidelity diagnostics
- **Noise models**: 9 Kraus channels (bit-flip, phase-flip, depolarising, amplitude damping, T1/T2, two-qubit depolarising, measurement error)
- **QAOA solver**: depth-1 ansatz, MPI-parallelised grid search over (γ, α), non-destructive sampling
- **QUBO construction**: closed-form Ising Hamiltonian from band-packing constraints, normalisation, quadratic coupling sparsity

### Applications
- **Electricity procurement**: 8-qubit bin-packing QUBO, hourly demand coverage, optional budget constraint, feasibility checking
- **QSVM (prototype)**: quantum feature-map kernel, kernel matrix via circuit fidelity, classical dual QP solver (planned)

### Physics
- **Big-endian qubit convention** (qubit 0 = MSB) maintained consistently across all files
- **Hermiticity enforcement** on density matrices post-operation
- **Completeness verification** on Kraus operator sets (Σ K†K = I)
- **Physical constraints** validated (T2 ≤ 2·T1, depolarising p ∈ [0,0.75], etc.)

## Key Concepts

### Big-Endian Convention
Qubit 0 is the most significant bit. For a 3-qubit state |q₀q₁q₂⟩, the index in the state vector is:
```
index = q₀·2² + q₁·2¹ + q₂·2⁰
```
This convention is enforced across `Gate::expand_two()`, `expand()`, and `DensityMatrix::partial_trace()`.

### Noise Model Design
Each Kraus channel is a **closed-form factory** in `NoiseModel.h`, decoupled from the simulator back-end. A channel is applied as:
```cpp
rho.apply_kraus(NoiseModel::depolarising(p, qubit, n_total));
```
The completeness relation Σ K†K = I is verified at call time.

### QAOA Cost Hamiltonian Diagonality
For classical QUBO problems, H_C contains only Z operators (no X or Y), so it is diagonal in the computational basis and the cost unitary exp(−iγH_C) factors exactly without Trotterisation. Off-diagonal terms (X or Y) would indicate a qualitatively different problem class and would require approximate time evolution.

### QSVM Kernel Geometry
The quantum advantage in QSVM is **geometric, not computational**: a quantum feature map |φ(x)⟩ can embed data into a richer Hilbert space than classical features allow. The kernel matrix K_Q(xᵢ, xⱼ) = |⟨φ(xⱼ) | φ(xᵢ)⟩|² is computed via circuit fidelity measurements. The dual QP solver and α vector are entirely classical and identical to standard SVM; the quantum component is only the kernel.

## Repository layout

| Path | Purpose |
|------|---------|
| `Gate.h` | Gate library (single- and two-qubit gates, register embedding). |
| `Qbit.h` | Pure state-vector simulator. |
| `DensityMatrix.h` | Density-matrix simulator (mixed states, noise, entropy, fidelity). |
| `NoiseModel.h` | Factory of Kraus-operator noise channels for `DensityMatrix`. |
| `QuantumSim.cpp` | MPI QAOA grid-search for a 2-qubit bin-packing QUBO. |
| `QuantumSim_V1.cpp` | Earlier version of the simulator, kept for reference. |
| `Electricity_procurement.cpp` | MPI QAOA for the electricity-procurement QUBO (8 qubits, entangling `Rzz`, optional budget). |
| `Electricity_procurement_noise.cpp` | The same, on `DensityMatrix` with Kraus noise: sweeps the noise level and tracks purity/fidelity/entropy. |
| `Results/` | Python plotting scripts; CSV outputs land here (git-ignored). |
| `Books/` | Reference texts (quantum computing introductions, related papers). |
| `Quantum.pdf` | The paper this simulator backs (git-ignored). |

to be. restructured as 

```
quantum-simulator/
├── README.md
├── CONTRIBUTING.md               # Dev setup, branch workflow
├── LICENSE                       # MIT
├── VERSION                       # 0.0.0
├── CHANGELOG.md                  # Release history
├── CMakeLists.txt
├── .clang-format                 # Code style
├── .github/
│   ├── workflows/
│   │   ├── ci.yml               # Build + test on push/PR
│   │   └── clang-format.yml     # Style enforcement
│   └── PULL_REQUEST_TEMPLATE.md
├── src/
│   ├── core/                    # Gate, Qbit, DensityMatrix, NoiseModel
│   │   ├── Gate.h
│   │   ├── Qbit.h
│   │   ├── DensityMatrix.h
│   │   ├── NoiseModel.h
│   │   └── CMakeLists.txt
│   ├── applications/            # QAOA, QSVM
│   │   ├── QuantumSimQAOA.cpp   # Noiseless 2-qubit sandbox
│   │   ├── ElectricityQAOA.cpp  # 8-qubit electricity procurement
│   │   ├── ElectricityNoise.cpp # Noisy version
│   │   └── CMakeLists.txt
│   └── CMakeLists.txt
├── include/
│   └── quantum_sim/             # Public headers (if building as library)
│       ├── gate.h
│       ├── simulator.h
│       └── noise.h
├── tests/
│   ├── CMakeLists.txt
│   ├── test_gate.cpp            # Unitary verification, tensor products
│   ├── test_density_matrix.cpp  # Partial trace, Kraus ops, entropy
│   ├── test_qsvm.cpp            # Kernel matrix, fidelity
│   └── fixtures/
│       └── bell_state.h         # Precomputed reference states
├── examples/
│   ├── CMakeLists.txt
│   ├── simple_bell.cpp          # Minimal entanglement demo
│   ├── qaoa_electricity.cpp     # Link to src/applications
│   └── README.md                # Usage walkthrough
├── docs/
│   ├── ARCHITECTURE.md          # Design, physics background, conventions
│   ├── FORMULATION.md           # QUBO → Ising mapping (electricity instance)
│   ├── NOISE_MODELS.md          # Kraus channels, physical interpretation
│   ├── API.md                   # Class reference (auto-gen friendly)
│   └── images/
│       ├── energy_surface.png
│       └── coupling_graph.svg
├── scripts/
│   ├── energy_surface.py        # 3D plot from CSV
│   ├── measurement_histogram.py # Bar chart
│   └── requirements.txt         # matplotlib, pandas
├── benchmarks/
│   ├── CMakeLists.txt
│   ├── mpi_scaling.cpp          # Strong scaling on dense coupling
│   └── density_matrix_vs_qbit.cpp
├── Results/                     # Generated outputs (in .gitignore)
│   └── .gitkeep
└── .gitignore
```

## Components

### `Gate.h`
Gates are `Eigen::MatrixXcd` wrappers carrying a qubit count.

- **Single-qubit:** `identity`, `X`, `Y`, `Z`, `H`, `Rx(θ)`, `Ry(θ)`, `Rz(θ)`.
- **Two-qubit:** `CNOT`, `CZ`, `SWAP`, `Rzz(θ)` (native ZZ-interaction gate),
  `controlled_U(U)` (a controlled-anything factory), and `Toffoli`.
- **Embedding into an n-qubit register:**
  - `expand(total, target)` lifts a single-qubit gate to the full `2ⁿ` space.
  - `expand_two(G, ctrl, tgt, total)` is the single general-purpose mechanism that
    lifts *any* 4×4 two-qubit gate onto an arbitrary qubit pair, including
    non-adjacent and reversed-order pairs. Convenience wrappers `cnot`, `cz`, `swap`,
    `rzz`, and `controlled` route through it.
- `tensor(A, B)` and `operator*` provide tensor products and sequential composition.

### `Qbit.h`
Pure state-vector simulator. Applies gates via an Eigen map over the state buffer
(no copy), supports `measure()` (collapsing, returns a bitstring) and exposes the raw
state through `get_state()` / `access_state()`.

### `DensityMatrix.h`
Density-matrix formalism, constructible from a state vector or a raw matrix:

- `apply(U)` — unitary evolution `ρ → UρU†`.
- `apply_kraus(ops)` — general CPTP channel, with a completeness check `Σ KᵢKᵢ = I`.
- `partial_trace(qubit)` — trace out a qubit.
- `partial_measurement(qubit)` and `measure()` — single-qubit and full-register
  measurement with collapse.
- `expectation(O)`, `purity()`, `fidelity(ψ)`, `probabilities()`.
- `entropy()` — von Neumann entropy via Eigen's `SelfAdjointEigenSolver`.

### `NoiseModel.h`
A factory of Kraus-operator sets, decoupled from `DensityMatrix` so the channels can be
proofread and extended independently. Channels: `bit_flip`, `phase_flip`,
`bit_phase_flip`, `depolarising`, `amplitude_damping` (from T1 and gate duration),
`generalised_amplitude_damping` (finite temperature), `t1_t2` (combined T1/T2 with
pure-dephasing separated out, returned as a `T1T2Channels` struct),
`depolarising_2q` (16 two-qubit Pauli operators), and `measurement_error` — a classical
readout confusion matrix rather than a Kraus channel on ρ.

## The QAOA simulator (`QuantumSim.cpp`)

A QAOA solver for a 2-qubit bin-packing QUBO, parallelised with MPI:

1. **Grid search** over the QAOA angles `(γ, α)` (`STEPS` per axis), each rank scanning a
   slice of the grid and writing its `(γ, α, energy)` rows to `Results/energy_surface.csv`.
2. **Best-parameter selection** via `MPI_Allreduce` with `MINLOC`.
3. **Sampling** of the best state across ranks (`SAMPLES` shots), reduced into
   `Results/measurement_histogram.csv`.
4. A timing summary (grid / measurement / total).

The cost is currently a hand-coded `compute_energy()` (`0.45·Z₀ + 0.15·Z₁`) with an
Rz-only cost unitary built on `Qbit`. **This is intentional and correct for this
instance:** the classical cost function has no interaction term between the two
variables, so the Hamiltonian is diagonal and the simulator does not yet produce
entanglement. That limitation is documented in the paper, not a bug.

`QuantumSim_V1.cpp` is the same program with the earlier symmetric coefficients
(`0.4·Z₀ + 0.4·Z₁`) and CSV output to the working directory instead of `Results/`.

## The electricity-procurement QAOA (`Electricity_procurement.cpp`)

The same MPI grid-search skeleton, "stepped up" to a real QUBO derived from
`electricity_qubo_formulation.md`. It selects 1-MW generation bands (nuclear / wind /
gas) to meet an hourly demand profile at least cost.

- **`build_problem()`** turns the band catalogue + demand into the closed-form Ising
  Hamiltonian of the formulation (Sections 4–9): one qubit per band placement, a linear
  `hᵢ`, and a quadratic `Jᵢⱼ` from per-hour coverage overlaps plus the optional budget
  term (C4, `B_max`). The budget term makes the coupling graph **dense** — every pair of
  bands is coupled through total spend.
- **Entanglement** is now genuine: the `Rzz` couplings entangle the qubits, unlike
  `QuantumSim.cpp` whose diagonal cost left the state a product state.
- The raw penalty weights push the Ising coefficients into the thousands, so `H_C` is
  **normalised by its largest coefficient** before the circuit — this keeps the QAOA
  angles `2γ·coef` O(1) over `γ∈[0,2π]` and never changes the optimum (a positive
  rescaling preserves `argmin`). Decoded costs are reported in real €.
- **Validation:** a brute-force `classical_optimum()` over all `2ⁿ` states, plus a
  decode + per-hour feasibility (`cov(t) ≥ d_t`) and budget check. The grid search is
  compared against it; the cheapest feasible sampled shot is the recommended procurement.
- **Measurement is non-destructive:** the optimised state is sampled from its Born
  distribution without collapsing (`QuantumSim.cpp` re-measured one collapsed copy, so
  every shot returned the same outcome — fixed here).
- Optional `argv[1]` overrides the grid resolution (default `STEPS = 1000`) for quick
  runs; output goes to the same `Results/*.csv` so the existing plot scripts apply.

> **Header fix.** Building this surfaced a latent bug in `Gate::Rx`: its lower-left entry
> had `+i·sin(θ/2)` (copied from `Ry`'s antisymmetric pattern) instead of `−i·sin(θ/2)`,
> making the gate Hermitian rather than unitary and inflating the state norm. Now
> corrected — `Rx` is a proper unitary, which also fixes `QuantumSim.cpp`'s mixer.

## The noisy version (`Electricity_procurement_noise.cpp`)

The open-system counterpart, built on `DensityMatrix.h` + `NoiseModel.h`. It solves the
**same instance at the same optimal `(γ*,α*)`** the noiseless run found, so the two are
directly comparable — the only change is the back-end and what it can express.

- **State:** a density matrix `ρ` (2ⁿ×2ⁿ) instead of a state vector, so mixed states and
  decoherence are representable. Each gate is `ρ → UρU†`; every gate is followed by a
  **Kraus channel** (`apply_kraus`): single-qubit `depolarising` on 1-qubit gates, the
  2-qubit `depolarising_2q` on each entangling `Rzz` (2-qubit gates are 5× noisier).
- **What it sweeps:** because `ρ` costs `O((2ⁿ)³)` per gate, it does *not* re-scan the
  angle grid; it fixes `(γ*,α*)` and **sweeps the physical noise strength `p`** (MPI-
  parallelised). That curve — quality vs. noise — is the open-system analogue of the
  noiseless energy surface, written to `Results/electricity_noise_sweep.csv`
  (`p, energy, purity, fidelity, entropy, …`).
- **`partial_trace`:** traces out all but one qubit to watch entanglement wash out. At
  `p=0` the reduced 1-qubit state is mixed *because of entanglement* (purity 0.665); as
  noise grows the whole state decoheres (global purity → 1/2ⁿ) and the marginal goes
  maximally mixed for a different reason — the two are distinguished by the global purity.
- **`partial_measurement` + `measurement_error`:** projects the register out qubit-by-
  qubit, then optionally distorts the classical bit with a readout confusion matrix.
- **Findings:** `p=0` reproduces the noiseless `⟨H_C⟩` exactly; purity → `1/256`,
  fidelity → 0, entropy → `ln 256` as `p` grows. Degradation is steep — 28 entangling
  gates at elevated 2-qubit error accumulate quickly — a realistic caution about deep
  QAOA on noisy hardware.

## Design decisions

- Stay on the **gate model** throughout (not annealing, not cat qubits).
- Build `DensityMatrix` as a **parallel** class to `Qbit`, keeping the pure, noise-free
  path intact.
- Keep all noise channels in a **dedicated factory** (`NoiseModel.h`), decoupled from the
  density matrix.
- Use **`expand_two()` as the one mechanism** for lifting two-qubit gates to the full
  register — no more ad-hoc CNOT/Rzz decompositions.

## Roadmap

**Done**
- Two-qubit gates and `expand_two()` in `Gate.h`.
- Full `DensityMatrix.h` (unitary + Kraus evolution, partial trace, measurement,
  expectation, entropy, purity, fidelity).
- `NoiseModel.h` channel factory.
- **Electricity-procurement QAOA** (`Electricity_procurement.cpp`): a `QuboHamiltonian`
  (linear `h` + quadratic `J`) built from a band catalogue and a per-hour coverage
  constraint, with `Rzz` couplings that produce genuine entanglement, an optional budget
  term, a brute-force classical check, and non-destructive sampling.
- Fixed the non-unitary `Gate::Rx` (see the header-fix note above).
- **Noisy density-matrix QAOA** (`Electricity_procurement_noise.cpp`): the same QAOA on
  `DensityMatrix` with per-gate Kraus noise, a noise-level sweep, and purity / fidelity /
  entropy / `partial_trace` / `partial_measurement` diagnostics.

**Next**
- **Longer term:** a `QEC.h` for simulated error correction (encoding circuits, syndrome
  measurement on ancillas via `partial_measurement`, correction lookup), reusing the
  infrastructure already in place — a natural next step now that noisy evolution is in.

### v0.1.0 (planned)
- [ ] QSVM kernel matrix computation (two-qubit feature maps, fidelity circuit)
- [ ] Classical dual QP solver (SciPy / CVXPY binding)
- [ ] Anomaly detection pipeline on synthetic billing data

### v0.2.0 (planned)
- [ ] Variational Quantum Eigensolver (VQE) for Hamiltonian simulation
- [ ] Quantum Error Correction (QEC) simulation framework
- [ ] Realistic gate calibration (pulse-level constraints)

## Building and running

Requires an MPI compiler and Eigen (here, Homebrew's `eigen3`):

```sh
mpic++ -std=c++23 -O2 \
  -I/opt/homebrew/Cellar/eigen/5.0.1/include/eigen3 \
  QuantumSim.cpp -o QuantumSim

mpirun -np 4 ./QuantumSim
```

The electricity-procurement solver builds and runs the same way (give it more ranks —
8 qubits with the dense budget coupling and a full 1000×1000 grid is heavy):

```sh
mpic++ -std=c++23 -O2 \
  -I/opt/homebrew/Cellar/eigen/5.0.1/include/eigen3 \
  Electricity_procurement.cpp -o Electricity_procurement

mpirun -np 8 ./Electricity_procurement        # full grid (≈20 min on 8 cores)
mpirun -np 8 ./Electricity_procurement 100     # quick 100×100 grid for a smoke test
```

Then plot the outputs:
```sh
python3 Results/energy_surface.py
python3 Results/measurement_histogram.py
```

---

## Physics References

- Kaye, P., Laflamme, R., & Mosca, M. (n.d.). An introduction to quantum computing. Oxford University Press on Demand.
- Grynberg, G., Aspect, A., & Fabre, C. (1997). Introduction aux lasers et à l’optique quantique. Ellipses Marketing.
- LaPierre, R. (2021). Introduction to Quantum computing. In The Materials Research Society series. 