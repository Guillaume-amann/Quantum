# Quantum — a from-scratch gate-model quantum simulator

A small C++ quantum-computing simulator built on [Eigen](https://eigen.tuxfamily.org)
and MPI, written to support the accompanying paper (`Quantum.pdf`). It takes the
gate-model approach (not annealing, not cat qubits) and is used to explore QAOA for
combinatorial optimisation, with a roadmap toward open-system (noisy) simulation and
an electricity-procurement QUBO.
à
The codebase has two simulation back-ends that sit side by side:

- **`Qbit`** — a pure state-vector simulator for the ideal, noise-free case.
- **`DensityMatrix`** — a density-matrix simulator for mixed states and noise, kept as
  a parallel class rather than a replacement for `Qbit`.

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

## Conventions

- **Big-endian qubit ordering** (qubit 0 = most significant bit) is used consistently
  across `Gate.h`, `Qbit.h`, `DensityMatrix.h`, and `NoiseModel.h`.
- Noise is expressed in the **Kraus / operator-sum** formalism; unitary evolution is the
  single-operator special case.

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
