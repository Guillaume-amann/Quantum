# Contributing to Quantum Simulator

## Development Workflow

### Branches
- **`main`** — production-ready code, all tests pass, tagged releases
- **`dev`** — integration branch, latest development (may be unstable between releases)
- **Feature branches** — branched from `dev`, merged back via pull request

### Naming Convention
```
feature-<feature-name>        # New feature (e.g. feat-qsvm)
bug-<bug-description>         # Bug fix (e.g. bug-measurement-collapse)
refactor-<scope>              # Code refactoring (e.g. refactor-gates)
```

### Workflow Steps

1. **Create branch from `dev`:**
   ```bash
   git checkout dev
   git pull origin dev
   git checkout -b feature-feature-name
   ```

2. **Work incrementally:** Commit early.
   ```bash
   git add MyFeature.h
   git commit -m "feat(core): add QuantumKernel class skeleton
   
   - Public methods: compute(), fidelity_matrix()
   - Private: feature_map_circuit()
   - [WIP] Integration with dual QP solver pending"
   ```

3. **Keep feature branch up-to-date:**
   ```bash
   git fetch origin
   git rebase origin/dev
   ```
   (Rebase, not merge, keeps feature history linear.)

4. **Push regularly:**
   ```bash
   git push origin feature-feature-name
   ```

5. **Open a Pull Request** when ready for review:
   - Target: `dev` (not `main`)
   - Use the PR template (auto-generated)
   - Link any related issues: `Resolves #42`
   - Wait for CI to pass, then request review

6. **Address review feedback:**
   ```bash
   git add <files>
   git commit -m "Address review: clarify Kraus operator contract"
   git push origin feature-feature-name
   ```
   (Do NOT force-push during review.)

7. **Merge:** Once approved and CI passes, maintainer squash-merges into `dev`:
   ```bash
   git checkout dev
   git pull origin dev
   git merge --squash feature-feature-name
   git commit -m "feat(qsvm): implement quantum kernel circuit (#45)"
   git push origin dev
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
git clone https://github.com/Guillaume-amann/Quantum.git
cd Quantum
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
ctest --output-on-failure
```

### Incremental Development
```bash
cmake --build build -j4 && ctest --output-on-failure
```

---

## Review Expectations

- **Maintainer reviews** all PRs targeting `dev` or `main`
- **Scope:** code correctness, style, tests
- **Feedback:** constructive, specific, suggests fixes
- **Approval:** requires approval before merge to `dev` or `main`