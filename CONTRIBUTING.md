# Contributing to Quantum Simulator

## Development Workflow

### Branches
- **`main`** — production-ready code, all tests pass, tagged releases
- **`develop`** — integration branch, latest development (may be unstable between releases)
- **Feature branches** — branched from `develop`, merged back via pull request

### Naming Convention
```
feat/<feature-name>        # New feature (e.g. feat/qsvm-kernel)
fix/<bug-description>      # Bug fix (e.g. fix/measurement-collapse)
refactor/<scope>           # Code refactoring (e.g. refactor/gate-docs)
docs/<section>             # Documentation only (e.g. docs/architecture)
perf/<improvement>         # Performance (e.g. perf/mpi-scaling)
test/<coverage>            # Test addition (e.g. test/kraus-completeness)
```

### Workflow Steps

1. **Create branch from `develop`:**
   ```bash
   git checkout develop
   git pull origin develop
   git checkout -b feat/your-feature-name
   ```

2. **Work incrementally:** Commit early, commit often (2–3 times per day is ideal).
   ```bash
   git add src/core/MyFeature.h
   git commit -m "feat(core): add QuantumKernel class skeleton
   
   - Public methods: compute(), fidelity_matrix()
   - Private: feature_map_circuit()
   - [WIP] Integration with dual QP solver pending"
   ```

3. **Keep feature branch up-to-date:**
   ```bash
   git fetch origin
   git rebase origin/develop
   ```
   (Rebase, not merge, keeps feature history linear.)

4. **Push regularly:**
   ```bash
   git push origin feat/your-feature-name
   ```

5. **Open a Pull Request** when ready for review:
   - Target: `develop` (not `main`)
   - Use the PR template (auto-generated)
   - Link any related issues: `Resolves #42`
   - Wait for CI to pass, then request review

6. **Address review feedback:**
   ```bash
   git add <files>
   git commit -m "Address review: clarify Kraus operator contract"
   git push origin feat/your-feature-name
   ```
   (Do NOT force-push during review.)

7. **Merge:** Once approved and CI passes, maintainer squash-merges into `develop`:
   ```bash
   git checkout develop
   git pull origin develop
   git merge --squash feat/your-feature-name
   git commit -m "feat(qsvm): implement quantum kernel circuit (#45)"
   git push origin develop
   ```

---

## Local Setup

### Prerequisites
- CMake ≥ 3.12
- C++17 compiler (GCC 9+, Clang 10+)
- Eigen 3.3+
- OpenMPI 3.0+
- Python 3.8+ (for scripts + plotting)

### macOS (Homebrew)
```bash
brew install cmake eigen openmpi python3
pip3 install matplotlib pandas numpy
```

### Ubuntu / Debian
```bash
sudo apt-get install -y cmake libeigen3-dev libopenmpi-dev
pip3 install matplotlib pandas numpy
```

### First Build
```bash
git clone https://github.com/<your-account>/quantum-simulator.git
cd quantum-simulator
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
ctest --output-on-failure
```

### Incremental Development
```bash
cmake --build build -j4 && ctest --output-on-failure
```

---

## Code Style

### Automatic Formatting
This project uses **Clang 14+** for automatic formatting. Install it:
```bash
brew install clang-format  # or apt-get install clang-format-14
```

### Run Before Commit
```bash
clang-format -i src/**/*.{h,cpp}
```

Or set up a pre-commit hook:
```bash
ln -s ../../.githooks/pre-commit .git/hooks/pre-commit
chmod +x .git/hooks/pre-commit
```

### Style Rules
- **Indentation:** 4 spaces (no tabs)
- **Column limit:** 100 characters
- **Braces:** Allman style (opening brace on new line for functions, same line for inline)
- **Naming:** `snake_case` for functions/variables, `CamelCase` for classes, `UPPER_CASE` for constants

### Example
```cpp
class QuantumKernel {
private:
    int n_qubits;
    vector<double> params;
    
    Gate feature_map_circuit(const vector<double>& x) const {
        // Implementation
    }

public:
    double compute_kernel(const vector<double>& x, const vector<double>& y) {
        Qbit psi = apply_feature_map(x);
        Qbit phi = apply_feature_map(y);
        return fidelity(psi, phi);
    }
};
```

---

## Commit Message Format

Follow the [Conventional Commits](https://www.conventionalcommits.org/) standard:

```
<type>(<scope>): <subject>

<body>

<footer>
```

### Type
- `feat` — new feature
- `fix` — bug fix
- `refactor` — code refactoring (no new feature, no bug fix)
- `perf` — performance improvement
- `docs` — documentation only
- `test` — test addition / modification
- `chore` — build, CI, dependency updates

### Scope
The area affected (e.g. `gate`, `noise-model`, `qsvm`, `mpi`)

### Subject
- Imperative mood ("add" not "added")
- No period
- Max 50 characters
- Lowercase

### Body
- Wrap at 72 characters
- Explain *why*, not *what*
- Reference any related issues: `Refs #42`

### Footer
- Breaking changes: `BREAKING CHANGE: description`
- Closes issues: `Closes #42` or `Resolves #42`

### Example
```
feat(density-matrix): add partial_trace method

Implement partial trace over arbitrary qubit k by constructing
the reduced density matrix via the sum of outer products over
the traced-out qubit's basis states.

This enables single-qubit marginal analysis and entanglement
diagnostics. Tested against known Bell state marginals.

Resolves #8
```

---

## Testing

### Write Tests First (TDD)
For any new feature, add a test in `tests/` *before* implementation:

```cpp
// tests/test_qsvm.cpp
#include <gtest/gtest.h>
#include "quantum_sim/gate.h"
#include "quantum_sim/qsvm.h"

TEST(QuantumKernel, IdenticalPointsReturnUnitKernel) {
    QuantumKernel kernel(4);  // 4 qubits
    vector<double> x = {0.1, 0.2, 0.3, 0.4};
    double K = kernel.compute(x, x);
    EXPECT_NEAR(K, 1.0, 1e-10);
}

TEST(QuantumKernel, KernelSymmetry) {
    vector<double> x = {0.1, 0.2};
    vector<double> y = {0.5, 0.6};
    double Kxy = kernel.compute(x, y);
    double Kyx = kernel.compute(y, x);
    EXPECT_NEAR(Kxy, Kyx, 1e-10);
}
```

### Run Tests
```bash
cmake --build build && ctest --output-on-failure
```

### Test Coverage (Optional)
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DENABLE_COVERAGE=ON
cmake --build build
ctest
gcovr --root . --filter 'src/.*' --print-summary
```

---

## Documentation

### When You Write Code, Document It
- **Headers**: Doxygen-style comments on public classes/methods
  ```cpp
  /// Compute the quantum kernel matrix between two data points.
  /// @param x Input vector (feature dimension)
  /// @param y Input vector (feature dimension)
  /// @return K_Q(x,y) = |⟨φ(y)|φ(x)⟩|², in [0,1]
  double compute(const vector<double>& x, const vector<double>& y);
  ```

- **Complex logic**: Inline comments explaining *why*
  ```cpp
  // Big-endian indexing: qubit 0 is MSB, so bit position = n-1-i
  int k_shift = n - 1 - qubit;
  ```

- **Public APIs**: Update [docs/API.md](docs/API.md) when you add a class or method

- **Design decisions**: Update [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) if you change how components interact

### Docs Checklist
- [ ] Doxygen comments on public methods
- [ ] Inline comments on non-obvious logic
- [ ] Examples in docstrings for complex APIs
- [ ] ARCHITECTURE.md updated if design changes
- [ ] CHANGELOG.md entry under "[Unreleased]"

---

## Pull Request Checklist

Before submitting a PR, ensure:

- [ ] Code builds without warnings: `cmake --build build`
- [ ] All tests pass: `ctest --output-on-failure`
- [ ] Code formatted: `clang-format -i src/**/*.{h,cpp}`
- [ ] Commit messages follow Conventional Commits
- [ ] CHANGELOG.md updated
- [ ] Docstrings added for public APIs
- [ ] No hardcoded debug prints or TODOs in final code
- [ ] Feature is linked to an issue (or issue created)

---

## Review Expectations

- **Maintainer reviews** all PRs targeting `develop` or `main`
- **Turnaround:** 2–3 business days
- **Scope:** code correctness, style, tests, documentation
- **Feedback:** constructive, specific, suggests fixes
- **Approval:** requires 1 approval before merge to `develop`; 2 for `main`

---

## Release Process

### Creating a Release
1. Update `VERSION` file and `CHANGELOG.md`
2. Create a tag:
   ```bash
   git tag -a v0.1.0 -m "First stable release: QSVM kernel, MPI scaling verified"
   git push origin v0.1.0
   ```
3. Create a GitHub Release from the tag (auto-fills from tag message)
4. Update `develop` version to next `-dev`:
   ```bash
   echo "0.2.0-dev" > VERSION
   git commit -am "chore: bump to v0.2.0-dev"
   git push origin develop
   ```

---

## Questions?

Refer to:
- [ARCHITECTURE.md](docs/ARCHITECTURE.md) — design & conventions
- [API.md](docs/API.md) — class reference
- GitHub Issues — existing problems & discussions
- Guillaume (physics/algorithm questions)

Thank you for contributing!
