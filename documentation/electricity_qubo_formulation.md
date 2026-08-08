# Electricity Procurement via Bin Packing — Complete QUBO/QAOA Formulation

## 1. Problem Statement

Electricity markets require matching a continuous, hourly demand profile with discrete
purchases of supply. Supply is bought in **bands**: each band is a contract that
delivers exactly **1 MW continuously** over a fixed number of consecutive hours,
starting at a given hour, at a given cost. Different generation types (nuclear, wind,
coal, gas, ...) offer bands of different durations and costs.

The goal is to select a subset of available bands that **meets demand at every hour**,
at **minimum total cost**, while allowing controlled overflow (buying more than
strictly needed) when doing so is cheaper than exact coverage.

This is modelled as a Quadratic Unconstrained Binary Optimization (QUBO) problem and
solved via QAOA.

> **Note on a discarded earlier model.** An earlier version of this formulation
> represented demand as a 2D grid (hour × MW-row) and used one binary variable per
> grid cell. This was abandoned: bands are continuous MW purchases over time, not
> stackable physical layers, so a "MW-layer ordering constraint" between rows is
> physically meaningless. The corrected model below uses one binary variable per
> **band placement** and a single **per-hour** coverage constraint — both simpler
> and physically accurate.

---

## 2. Sets

| Symbol | Definition |
|---|---|
| `T = {0, 1, ..., T-1}` | Time slots (hours) in the procurement horizon |
| `B` | Band *types* (e.g. nuclear, wind, coal, gas), each with a fixed duration `δ_b` and per-MWh cost `c_b` |
| `P = { (b,s) : b∈B, s∈T, s+δ_b ≤ T }` | **Derived set** of all valid band *placements* — type `b` starting at hour `s`. Generated programmatically from `B` and `T`, never hand-enumerated. |

Each placement `i = (b,s) ∈ P` is treated as a single index for the rest of this
document.

---

## 3. Parameters

| Symbol | Meaning | Source |
|---|---|---|
| `d_t` | Demand at hour `t`, in MW (integer) | given, exogenous |
| `δ_i` | Duration of placement `i` (hours), inherited from its band type `b` | given |
| `c_i` | Total cost of placement `i` = `c_b · δ_b` | given |
| `λ_cov` | Coverage penalty weight | chosen, must dominate total achievable cost |
| `B_max` | Budget cap (only if the optional budget constraint is enabled) | chosen |
| `λ_bud` | Budget penalty weight (only if budget constraint enabled) | chosen |

---

## 4. Derived Quantities

These are not new variables — they are fixed numbers, computable once `B`, `T`, and
`d_t` are known, **before** any optimization runs.

**Coverage set of an hour** — placements active at hour `t`:
```
S_t = { i ∈ P : i covers hour t }
```

**Total demand under a placement**:
```
D_i = Σ_{t ∈ T_i} d_t      where T_i = the set of hours placement i covers
```

**Pairwise temporal overlap** (hours during which both placements i and j are active):
```
overlap_ij = |T_i ∩ T_j| = Σ_{t∈T} A_{i,t}·A_{j,t}
           = max(0,  min(s_i+δ_i, s_j+δ_j) − max(s_i, s_j))      (closed form, contiguous bands)
```
where `A_{i,t} ∈ {0,1}` is the fixed indicator "placement i covers hour t".

---

## 5. Decision Variables

```
x_i ∈ {0,1}     for every i ∈ P
```
`x_i = 1` means placement `i` is purchased. There is exactly one qubit per valid band
placement — not per grid cell.

**Derived (not a separate variable) — hourly coverage**:
```
cov(t) = Σ_{i ∈ S_t} x_i
```
the total MW supplied at hour `t` given the current choice of `x`.

---

## 6. Constraints

| ID | Statement | Status |
|---|---|---|
| C1 | MW-layer ordering | **Dropped** — artefact of the discarded grid model; bands are continuous MW quantities, not stackable rows, so no such constraint exists physically. |
| C2 | No same-type temporal overlap | **Dropped** — each band type represents an unlimited, identical contract; buying two overlapping placements of the same type is valid bin-packing behaviour, not a violation. |
| C3 | Grid boundary: `s + δ_b ≤ T` | **Kept** — enforced for free when generating `P`; no penalty term needed. |
| C4 | Budget cap: `Σ_i c_i x_i ≤ B_max` | **Kept, optional** — encoded as a soft QUBO penalty (Section 7.3). |

The only *hard* requirement left to enforce via penalty is demand coverage, handled
below as part of the objective (Section 7.2), using the symmetric penalty (**Option
B**, see below) rather than a separate hard/soft split.

---

## 7. Objective Function — Derivation

### 7.1 Cost term

```
Σ_{i∈P} c_i · x_i
```
Linear in `x`; no derivation needed.

### 7.2 Coverage term (Option B — symmetric penalty, no slack qubits)

```
λ_cov · Σ_{t∈T} (d_t − cov(t))²
```
Penalizes both deficit and surplus equally per unit; the asymmetric "overflow is
cheap" behaviour emerges from the interaction with the cost term, not from a separate
overflow weight. (Earlier considered: an asymmetric split with `λ_over` and exact slack
qubits — Option A — rejected here in favour of the qubit-free Option B.)

**Expansion.** Let `cov(t) = Σ_{i∈S_t} x_i`. Expanding the square:
```
(d_t − cov(t))² = d_t² − 2d_t·cov(t) + cov(t)²
```
Squaring the sum `cov(t)` (multiply every term in the sum by every term):
```
cov(t)² = Σ_{i∈S_t} x_i²  +  2 Σ_{i<j ∈ S_t} x_i x_j
```
Since `x_i ∈ {0,1}`, `x_i² = x_i`. Substituting back and summing over all `t`:
```
λ_cov Σ_t (d_t−cov(t))²
  = λ_cov Σ_t d_t²                                   (constant)
  + λ_cov Σ_i (δ_i − 2D_i) x_i                        (linear)
  + 2λ_cov Σ_{i<j} overlap_ij · x_i x_j               (quadratic)
```
The linear coefficient comes from `Σ_{t∈T_i}(1−2d_t) = δ_i − 2D_i`. The quadratic
coefficient comes from counting, for each pair `(i,j)`, how many hours they jointly
appear in `S_t` — exactly `overlap_ij`.

### 7.3 Budget term (optional, only if C4 enabled)

```
λ_bud · (Σ_i c_i x_i − B_max)²
```
**Expansion.** Let `C = Σ_i c_i x_i`. Binomial expansion:
```
(C − B_max)² = C² − 2B_max·C + B_max²
```
Squaring `C` (sum with coefficients, same rule as 7.2 but with weights `c_i`):
```
C² = Σ_i c_i² x_i²  +  2 Σ_{i<j} c_i c_j x_i x_j  = Σ_i c_i² x_i + 2Σ_{i<j} c_ic_j x_ix_j
```
(again using `x_i²=x_i`). And `−2B_max·C` distributes directly: `Σ_i (−2B_max c_i) x_i`.
Combining and multiplying by `λ_bud`:
```
λ_bud(C−B_max)²
  = λ_bud B_max²                                      (constant)
  + Σ_i λ_bud(c_i² − 2B_max c_i) x_i                   (linear)
  + 2λ_bud Σ_{i<j} c_i c_j x_i x_j                     (quadratic)
```
**Structural note**: unlike the coverage term, this quadratic coupling is **dense** —
nonzero for every pair `(i,j)` regardless of whether they overlap in time, since total
spend is a global sum. Enabling C4 turns a sparse, time-local coupling graph into a
fully connected one.

---

## 8. QUBO Canonical Form

Collecting Sections 7.1–7.3:
```
E(x) = const  +  Σ_i a_i x_i  +  Σ_{i<j} b_ij x_i x_j
```

```
const = λ_cov · Σ_{t∈T} d_t²                              [+ λ_bud · B_max²                if C4]

a_i   = c_i  +  λ_cov·(δ_i − 2D_i)                          [+ λ_bud·(c_i² − 2B_max·c_i)     if C4]

b_ij  = 2λ_cov · overlap_ij                                 [+ 2λ_bud · c_i · c_j            if C4]
```

All three quantities are closed-form functions of problem data (band catalogue, demand
profile, penalty weights) — no search or optimization is required to compute them.

---

## 9. Ising Mapping

Substitute `x_i = (1 − Z_i)/2` (standard QUBO → Ising transformation):
```
a_i x_i        = a_i/2 − (a_i/2) Z_i
b_ij x_i x_j   = (b_ij/4)(1 − Z_i − Z_j + Z_i Z_j)
```

**Collect by operator type** (constant / single-Z / ZZ):

```
H_C = E_0  +  Σ_i h_i Z_i  +  Σ_{i<j} J_ij Z_i Z_j
```

```
E_0  = const + Σ_i a_i/2 + Σ_{i<j} b_ij/4        (dropped for QAOA — shifts all energies equally, never changes the optimum)

h_i  = −a_i/2  −  Σ_{j≠i} b_ij/4

J_ij =  b_ij/4
```

**Substituting the closed forms from Section 8:**
```
h_i  = −[c_i + λ_cov(δ_i−2D_i)] / 2   −   Σ_{j≠i} [2λ_cov·overlap_ij + 2λ_bud·c_ic_j] / 4
                                                                  (drop the λ_bud term if C4 disabled)

J_ij =  λ_cov·overlap_ij / 2   +   λ_bud·c_i·c_j / 2
                                                                  (drop the λ_bud term if C4 disabled)
```

These are exactly the `H.h[i]` and `H.J[i][j]` entries to populate in the
`QuboHamiltonian` struct.

---

## 10. QAOA Circuit

`H_C` is a sum of **mutually commuting, diagonal** Pauli terms (every term built purely
from `Z`), so the cost-unitary exponential factors *exactly* — no Trotterization or
approximation needed:
```
exp(−iγ H_C) = [ Π_i exp(−iγh_i Z_i) ] · [ Π_{i<j} exp(−iγJ_ij Z_iZ_j) ]
```
Matching to the gate definitions:
```
Rz(θ)  := exp(−iθ/2 · Z)        ⟹   exp(−iγh_i Z_i)     = Rz(2γh_i)        on qubit i
Rzz(θ) := exp(−iθ/2 · Z⊗Z)      ⟹   exp(−iγJ_ij Z_iZ_j)  = Rzz(2γJ_ij)      on qubits i,j
```

**Full depth-1 QAOA circuit:**
```
1. Prepare |+⟩^⊗|P| via H on every qubit
2. Cost unitary U_C(γ):
     for each i with h_i ≠ 0:        apply Rz(2γh_i) on qubit i
     for each (i,j) with J_ij ≠ 0:   apply Rzz(2γJ_ij) on qubits i,j
3. Mixer unitary U_M(α):
     for each qubit i:                apply Rx(2α) on qubit i
```
Coupling sparsity follows directly from `overlap_ij`: pairs of bands that never share an
hour have `J_ij = 0` (and, with C4 disabled, need no `Rzz` gate at all between them).

---

## 11. Implementation Checklist

1. Generate `P` programmatically: loop over `(b,s)` pairs satisfying C3 (`s+δ_b≤T`).
2. Compute `δ_i, D_i` per placement, and `overlap_ij` per pair (matrix or interval
   formula).
3. Compute `a_i, b_ij, const` from Section 8.
4. Apply the Ising substitution to get `h_i, J_ij, E_0` from Section 9.
5. Populate `QuboHamiltonian.h` and `QuboHamiltonian.J`.
6. Choose `λ_cov` (≫ max achievable total cost) and decide whether C4 is active
   (and if so, `B_max`, `λ_bud`).
7. Run classical brute-force (`classical_optimum()`) on the new Hamiltonian for
   validation.
8. Run the existing MPI QAOA grid search (`apply_qaoa`, unchanged) and compare.
9. Sample, decode, and check feasibility (`Σ_{i∈S_t} x_i ≥ d_t` for all `t`).

No changes are required to `apply_qaoa()`, the MPI grid-search loop, or the
measurement/sampling logic — they are already fully generic with respect to any
`QuboHamiltonian`. All remaining work is in `build_problem()`.
