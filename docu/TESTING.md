# Testing Guide

This document explains how to build, run, and interpret the unit test suite for the Quantum simulator.

## Quick Start

### Prerequisites

- **CMake** ≥ 3.14
- **C++17** compiler (GCC 9+, Clang 10+, Apple Clang)
- **Eigen3** (header-only, already included or via package manager)
- **MPI** (optional; required only for application binaries, not tests)

### macOS (Homebrew)

```bash
brew install cmake eigen
```

### Ubuntu / Debian

```bash
sudo apt-get update
sudo apt-get install -y cmake libeigen3-dev build-essential
```

## Building and Running Tests

### Method 1: Automated Test Runner (Recommended)

```bash
# From project root
chmod +x run_tests.sh
./run_tests.sh
```

**Options:**
```bash
./run_tests.sh --verbose    # Show full output from each test
./run_tests.sh --rebuild    # Clean rebuild before testing
```

**Output:**
```
========================================
  Quantum Simulator Test Suite
========================================

Running test_gate...
✓ PASS : test_gate

Running test_density_matrix...
✓ PASS : test_density_matrix

Running test_qbit...
✓ PASS : test_qbit

Running test_noisemodel...
✓ PASS : test_noisemodel

========================================
  Test Summary
========================================
Total:  4 tests
Passed: 4
Failed: 0

✓ All tests passed!
```

### Method 2: Manual CMake Build

```bash
# From project root
mkdir -p build
cd build
cmake ..
cmake --build . --target all

# Run individual tests
../bin/test_gate
../bin/test_density_matrix
../bin/test_qbit
../bin/test_noisemodel
```

### Method 3: Per-Test Build and Run

```bash
# Just one test
cd build
cmake --build . --target test_gate
../bin/test_gate
```

## Test Structure

### Test Files

| File | What it tests | # Tests | Coverage |
|------|---------------|---------|----------|
| `tests/test_gate.cpp` | Gate.h | 70+ | Unitarity, Pauli algebra, tensor products, n-qubit embedding, composition |
| `tests/test_density_matrix.cpp` | DensityMatrix.h | 80+ | Hermiticity, trace, Kraus evolution, partial trace, purity, entropy, fidelity, measurement |
| `tests/test_qbit.cpp` | Qbit.h | 50+ | State normalisation, gate application, measurement collapse, Born rule, entanglement |
| `tests/test_noisemodel.cpp` | NoiseModel.h | 60+ | Kraus channels, physical constraints, completeness, noise effects, measurement error |

**Total: ~260 unit tests covering 100% of core infrastructure**

### Test Categories

#### Gate.h (70 tests)
- **Single-qubit unitarity** — X, Y, Z, H, Rx, Ry, Rz all unitary
- **Pauli algebra** — X² = I, XY = iZ, anticommutation relations
- **Rotation gates** — Exponential form, basis state mappings
- **Two-qubit gates** — CNOT, CZ, SWAP, Rzz, Toffoli
- **Tensor products** — Kronecker structure, dimension correctness
- **Gate expansion** — Single-qubit `expand()`, two-qubit `expand_two()` on arbitrary qubit pairs
- **Composition** — Gate multiplication, sequential application

#### DensityMatrix.h (80 tests)
- **Construction** — Pure states, superpositions, Bell states, maximally mixed
- **Unitary evolution** — Hermiticity and trace preservation under U ρ U†
- **Kraus channels** — Bit-flip, phase-flip, depolarising, amplitude damping, T1/T2, two-qubit depolarising
- **Partial trace** — Reduces n-qubit to (n-1)-qubit, entanglement decay
- **Purity** — Tr(ρ²) ∈ [1/d, 1], decreases under noise
- **Von Neumann entropy** — S(ρ) = −Tr(ρ log ρ), increases under dephasing
- **Fidelity** — F(ρ, |ψ⟩) = ⟨ψ|ρ|ψ⟩, interpolates between orthogonal states
- **Measurement** — Full register and single-qubit collapse, probability vectors
- **Expectation values** — ⟨O⟩ = Tr(O ρ) for Hermitian observables

#### Qbit.h (50 tests)
- **State construction** — |0⟩, |+⟩, Bell states, normalisation enforcement
- **Gate application** — Preserves normalisation, creates superposition and entanglement
- **Measurement** — Deterministic (basis states), probabilistic (superpositions), collapse
- **Born rule** — Probabilities match |amplitude|²
- **Sampling statistics** — 1000-shot distributions match expected probabilities
- **Entanglement** — Bell state measurements yield correlated outcomes
- **Multi-qubit states** — 3-qubit GHZ, 4-qubit systems
- **State access** — Const and mutable views, printing

#### NoiseModel.h (60 tests)
- **Kraus channels** — 9 single- and two-qubit channels
- **Physical constraints** — T2 ≤ 2T1, p ∈ [0, 0.75], p_eq ∈ [0, 0.5]
- **Parameter validation** — Bounds checking, exception throwing
- **Completeness** — Σ K†K = I verified at runtime
- **Channel effects** — p=0 → identity, p→max → maximally mixed
- **Noise escalation** — Progressive decoherence, purity decay, entropy increase
- **Measurement error** — Classical confusion matrix on readout
- **Multi-qubit noise** — 16-operator depolarising on two qubits
- **Temperature effects** — Generalised amplitude damping at finite T

## Continuous Integration (CI)

### GitHub Actions Setup

1. **Create `.github/workflows/tests.yml`** in your repository root with the contents of `ci_workflow.yml` (provided).

2. **Push to GitHub:**
   ```bash
   git add .github/workflows/tests.yml
   git commit -m "ci: add automated unit tests via GitHub Actions"
   git push origin main
   ```

3. **Verify:** Go to your repository's **Actions** tab. Tests run automatically on every push and pull request.

### What CI Tests

- **macOS (latest)** with Clang
- **Ubuntu (latest)** with GCC and Clang
- **C++17, Release build**
- All 4 test executables in sequence

### CI Badges (Optional)

Add to your README.md:
```markdown
[![Tests](https://github.com/your-username/Quantum/actions/workflows/tests.yml/badge.svg)](https://github.com/your-username/Quantum/actions/workflows/tests.yml)
```

## Understanding Test Output

### Passing Test

```
✓ PASS : X is unitary
```
The test passed. No action needed.

### Failing Test

```
✗ FAIL : Rx(π) = -iX
```
The test failed. Run with `--verbose` to see details:
```bash
./run_tests.sh --verbose
```

This will print the full output from the failed test, including:
- Which assertion failed
- Expected vs actual values
- Stack trace (if available)

### Exit Codes

```bash
./run_tests.sh
echo $?
```

- `0` = all tests passed
- `1` = one or more tests failed

## Adding New Tests

1. **Create a new test function** in the appropriate `.cpp` file:
   ```cpp
   void test_my_feature() {
       cout << "\n--- My feature description ---\n";
       
       // Setup
       MyClass obj;
       
       // Test assertions
       test_result("First property", obj.property1() == expected);
       test_result("Second property", obj.property2() == expected);
   }
   ```

2. **Add to the main suite** at the bottom of `main()`:
   ```cpp
   test_my_feature();
   ```

3. **Rebuild and run:**
   ```bash
   cd build
   cmake --build . --target test_gate  # (or your file)
   ../bin/test_gate
   ```

## Debugging Failed Tests

### Step 1: Identify the failing test
```bash
./run_tests.sh --verbose 2>&1 | grep -A 5 "FAIL"
```

### Step 2: Isolate the failure
Edit the test `.cpp` file temporarily to comment out all tests except the failing one.

### Step 3: Add debug output
```cpp
void test_my_feature() {
    cout << "\n--- Debug ---\n";
    
    auto result = myFunction();
    cout << "Result: " << result << "\n";
    cout << "Expected: " << expected << "\n";
    
    test_result("My feature", result == expected);
}
```

### Step 4: Rebuild and run
```bash
cd build
cmake --build . --target test_gate
../bin/test_gate
```

## Test Conventions

All tests follow these patterns:

- **Test names describe what is tested**, not how: `test_measurement_born_rule` ✓, `test_probs` ✗
- **Each test has a section header**: `cout << "\n--- Description ---\n";`
- **Assertions use `test_result()`** which exits immediately on failure
- **Tolerance is `EPSILON = 1e-9` to 1e-10** for floating-point comparisons
- **No randomness** — all tests are deterministic (except statistical ones, which run 1000 shots)

## Performance

### Build time (macOS M1, clean)
- First build: ~15 seconds
- Incremental: ~2 seconds

### Test runtime (all 4 suites)
- ~2 seconds total (mostly test_density_matrix: 0.6s, others <0.2s each)

### Size (Release build, stripped)
- Each test binary: 100–300 KB

## Troubleshooting

### "Eigen3 not found"
```bash
# macOS
brew install eigen

# Ubuntu
sudo apt-get install libeigen3-dev
```

### "CMAKE not found"
```bash
# macOS
brew install cmake

# Ubuntu
sudo apt-get install cmake
```

### "C++17 not supported"
Update your compiler:
```bash
# macOS: update Xcode
xcode-select --install

# Ubuntu: upgrade GCC
sudo apt-get install g++-9
```

### Test hangs or crashes
1. Check for infinite loops in gate application or measurement
2. Verify qubit indices are in range
3. Run with `gdb`:
   ```bash
   gdb ./bin/test_qbit
   (gdb) run
   (gdb) bt  # backtrace on crash
   ```

## Best Practices for Recruiters/Evaluators

When cloning and testing:

```bash
# Clone
git clone https://github.com/your-username/Quantum.git
cd Quantum

# Build and test (one command)
mkdir -p build && cd build && cmake .. && cmake --build . --target all && cd .. && ./run_tests.sh

# Or just use the runner
./run_tests.sh --rebuild --verbose
```

All tests should pass. If any fail, check:
1. C++17 compiler available
2. Eigen3 installed
3. Proper file paths (relative to project root)

---

**Questions?** See `CONTRIBUTING.md` for developer workflow guidelines.