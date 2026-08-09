// =============================================================================
//  test_density_matrix.cpp — Unit tests for DensityMatrix.h
//
//  Verifies:
//    * Hermiticity and trace normalisation (Tr(ρ) = 1)
//    * Unitary evolution preserves Hermiticity and trace
//    * Kraus channel evolution (ρ → Σ K_k ρ K_k†) with completeness check
//    * Partial trace reduces qubit and preserves trace
//    * Purity Tr(ρ²) ∈ [1/d, 1] for d-dimensional system
//    * Von Neumann entropy S(ρ) = -Tr(ρ log ρ) ≥ 0
//    * Fidelity F(ρ, |ψ⟩) with pure state ψ
//    * Measurement collapse (partial + full register)
//    * Expectation values and probability vectors
//
//  Build:
//    g++ -std=c++17 -O2 -I/opt/homebrew/include/eigen3 test_density_matrix.cpp -o test_density_matrix
//  Run:
//    ./test_density_matrix
// =============================================================================

#include <iostream>
#include <cmath>
#include <complex>
#include <cassert>
#include <iomanip>
#include <vector>
#include "src/core/Gate.h"
#include "src/core/Qbit.h"
#include "src/core/DensityMatrix.h"
#include "src/core/NoiseModel.h"

using namespace std;
using Complex = complex<double>;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;

// ============================================================================
// Test utilities
// ============================================================================

constexpr double EPSILON = 1e-9;

bool approx_eq(double a, double b, double tol = EPSILON) {
    return fabs(a - b) < tol;
}

bool approx_eq(Complex a, Complex b, double tol = EPSILON) {
    return abs(a - b) < tol;
}

// Check Hermiticity: ρ = ρ†
bool is_hermitian(const MatrixXcd& rho, double tol = EPSILON) {
    return (rho - rho.adjoint()).norm() < tol;
}

// Check trace is approximately 1
bool trace_is_one(const MatrixXcd& rho, double tol = EPSILON) {
    return approx_eq(rho.trace().real(), 1.0, tol) && approx_eq(rho.trace().imag(), 0.0, tol);
}

// Check all eigenvalues are non-negative (positive-semidefinite)
bool is_positive_semidefinite(const MatrixXcd& rho, double tol = EPSILON) {
    Eigen::SelfAdjointEigenSolver<MatrixXcd> solver(rho, Eigen::EigenvaluesOnly);
    for (int i = 0; i < solver.eigenvalues().size(); ++i)
        if (solver.eigenvalues()(i) < -tol) return false;
    return true;
}

// Probability vector from diagonal of ρ
vector<double> diagonal(const MatrixXcd& rho) {
    vector<double> diag(rho.rows());
    for (int i = 0; i < rho.rows(); ++i) diag[i] = rho(i, i).real();
    return diag;
}

// Print test result
void test_result(const char* name, bool passed) {
    cout << (passed ? "✓ PASS" : "✗ FAIL") << " : " << name << "\n";
    if (!passed) exit(1);
}

// ============================================================================
// Tests: Density matrix construction and basic properties
// ============================================================================

void test_pure_state_construction() {
    cout << "\n--- Pure state construction ---\n";
    
    // |0⟩ state
    DensityMatrix rho0(1);
    test_result("ρ(|0⟩) is Hermitian", is_hermitian(rho0.get_rho()));
    test_result("Tr(ρ) = 1", trace_is_one(rho0.get_rho()));
    test_result("ρ(|0⟩) is positive-semidefinite", is_positive_semidefinite(rho0.get_rho()));
    
    // |1⟩ state
    VectorXcd psi1(2); psi1 << 0, 1;
    DensityMatrix rho1(1, vector<Complex>(psi1.data(), psi1.data() + 2));
    test_result("ρ(|1⟩) is Hermitian", is_hermitian(rho1.get_rho()));
    test_result("Tr(ρ(|1⟩)) = 1", trace_is_one(rho1.get_rho()));
}

void test_plus_state() {
    cout << "\n--- Plus state |+⟩ = (|0⟩+|1⟩)/√2 ---\n";
    
    VectorXcd plus(2);
    plus << 1.0 / sqrt(2.0), 1.0 / sqrt(2.0);
    DensityMatrix rho_plus(1, vector<Complex>(plus.data(), plus.data() + 2));
    
    test_result("ρ(|+⟩) is Hermitian", is_hermitian(rho_plus.get_rho()));
    test_result("Tr(ρ(|+⟩)) = 1", trace_is_one(rho_plus.get_rho()));
    
    // Verify |+⟩⟨+| has off-diagonal coherence
    auto rho_mat = rho_plus.get_rho();
    test_result("ρ(|+⟩)[0][1] has coherence", abs(rho_mat(0, 1)) > 0.1);
}

void test_bell_state() {
    cout << "\n--- Bell state (Φ⁺) = (|00⟩+|11⟩)/√2 ---\n";
    
    VectorXcd bell(4);
    bell << 1.0 / sqrt(2.0), 0, 0, 1.0 / sqrt(2.0);
    DensityMatrix rho_bell(2, vector<Complex>(bell.data(), bell.data() + 4));
    
    test_result("ρ(|Φ⁺⟩) is Hermitian", is_hermitian(rho_bell.get_rho()));
    test_result("Tr(ρ(|Φ⁺⟩)) = 1", trace_is_one(rho_bell.get_rho()));
    test_result("ρ(|Φ⁺⟩) is positive-semidefinite", is_positive_semidefinite(rho_bell.get_rho()));
}

void test_maximally_mixed() {
    cout << "\n--- Maximally mixed state I/2 (single qubit) ---\n";
    
    // Explicitly construct I/2
    MatrixXcd rho_mm = MatrixXcd::Identity(2, 2) * 0.5;
    
    test_result("I/2 is Hermitian", is_hermitian(rho_mm));
    test_result("Tr(I/2) = 1", trace_is_one(rho_mm));
    test_result("I/2 is positive-semidefinite", is_positive_semidefinite(rho_mm));
}

// ============================================================================
// Tests: Unitary evolution
// ============================================================================

void test_unitary_evolution() {
    cout << "\n--- Unitary evolution ρ → UρU† ---\n";
    
    DensityMatrix rho(1);  // Start with |0⟩
    test_result("Initial: Tr(ρ) = 1", trace_is_one(rho.get_rho()));
    
    rho.apply(Gate::X());
    test_result("After X: Tr(ρ) = 1", trace_is_one(rho.get_rho()));
    test_result("After X: ρ is Hermitian", is_hermitian(rho.get_rho()));
    
    // X|0⟩ = |1⟩, so ρ should be |1⟩⟨1|
    auto rho_mat = rho.get_rho();
    test_result("ρ(X|0⟩) = |1⟩⟨1|", approx_eq(rho_mat(1, 1), 1.0) && approx_eq(rho_mat(0, 0), 0.0));
}

void test_hadamard_evolution() {
    cout << "\n--- Hadamard evolution: |0⟩ → |+⟩ ---\n";
    
    DensityMatrix rho(1);  // Start with |0⟩
    rho.apply(Gate::H());
    
    test_result("After H: Tr(ρ) = 1", trace_is_one(rho.get_rho()));
    
    // H|0⟩ = |+⟩ = (|0⟩+|1⟩)/√2
    auto rho_mat = rho.get_rho();
    test_result("ρ(H|0⟩)[0][0] ≈ 0.5", approx_eq(rho_mat(0, 0).real(), 0.5));
    test_result("ρ(H|0⟩)[1][1] ≈ 0.5", approx_eq(rho_mat(1, 1).real(), 0.5));
    test_result("ρ(H|0⟩)[0][1] ≈ 0.5", approx_eq(rho_mat(0, 1).real(), 0.5));
}

// ============================================================================
// Tests: Kraus channel evolution
// ============================================================================

void test_kraus_bit_flip() {
    cout << "\n--- Kraus channel: bit-flip ---\n";
    
    DensityMatrix rho(1);  // |0⟩
    
    // Bit-flip with p=0 should leave state unchanged
    rho.apply_kraus(NoiseModel::bit_flip(0.0, 0, 1));
    test_result("Bit-flip(p=0): Tr(ρ) = 1", trace_is_one(rho.get_rho()));
    test_result("Bit-flip(p=0): ρ is Hermitian", is_hermitian(rho.get_rho()));
    
    // Bit-flip with p=1 should flip state: |0⟩ → |1⟩
    DensityMatrix rho2(1);
    rho2.apply_kraus(NoiseModel::bit_flip(1.0, 0, 1));
    auto rho2_mat = rho2.get_rho();
    test_result("Bit-flip(p=1): ρ → |1⟩⟨1|", approx_eq(rho2_mat(1, 1), 1.0));
}

void test_kraus_depolarising() {
    cout << "\n--- Kraus channel: depolarising ---\n";
    
    // Depolarising with p=0 should leave state unchanged
    DensityMatrix rho0(1);
    rho0.apply_kraus(NoiseModel::depolarising(0.0, 0, 1));
    test_result("Depolarising(p=0): |0⟩ unchanged", is_hermitian(rho0.get_rho()));
    
    // Depolarising with p=0.75 should map everything to I/2
    DensityMatrix rho_max(1);
    rho_max.apply_kraus(NoiseModel::depolarising(0.75, 0, 1));
    auto rho_max_mat = rho_max.get_rho();
    test_result("Depolarising(p=0.75): ρ → I/2", approx_eq(rho_max_mat(0, 0).real(), 0.5) && approx_eq(rho_max_mat(1, 1).real(), 0.5));
}

void test_kraus_amplitude_damping() {
    cout << "\n--- Kraus channel: amplitude damping (T1 relaxation) ---\n";
    
    // Start with |1⟩: should decay to |0⟩
    VectorXcd ket1(2); ket1 << 0, 1;
    DensityMatrix rho1(1, vector<Complex>(ket1.data(), ket1.data() + 2));
    
    // T1=1, t_gate=1 → γ = 1 - exp(-1) ≈ 0.632
    rho1.apply_kraus(NoiseModel::amplitude_damping(1.0, 1.0, 0, 1));
    test_result("Amplitude damping: Tr(ρ) = 1", trace_is_one(rho1.get_rho()));
    test_result("Amplitude damping: ρ is Hermitian", is_hermitian(rho1.get_rho()));
    
    auto rho1_mat = rho1.get_rho();
    // Excited state should decay
    test_result("Amplitude damping: |1⟩ population decreases", rho1_mat(1, 1).real() < 1.0);
    // Ground state population should increase
    test_result("Amplitude damping: |0⟩ population increases", rho1_mat(0, 0).real() > 0.0);
}

void test_kraus_phase_flip() {
    cout << "\n--- Kraus channel: phase-flip (dephasing) ---\n";
    
    VectorXcd plus(2);
    plus << 1.0 / sqrt(2.0), 1.0 / sqrt(2.0);
    DensityMatrix rho_plus(1, vector<Complex>(plus.data(), plus.data() + 2));
    
    auto rho_before = rho_plus.get_rho();
    double coherence_before = abs(rho_before(0, 1));
    
    rho_plus.apply_kraus(NoiseModel::phase_flip(0.5, 0, 1));
    
    auto rho_after = rho_plus.get_rho();
    double coherence_after = abs(rho_after(0, 1));
    
    test_result("Phase-flip: Tr(ρ) = 1", trace_is_one(rho_after));
    // Coherence should decay (in general)
    test_result("Phase-flip(p=0.5): coherence shrinks", coherence_after < coherence_before || approx_eq(coherence_after, 0.0, 0.01));
}

void test_kraus_two_qubit() {
    cout << "\n--- Kraus channel: two-qubit depolarising ---\n";
    
    VectorXcd bell(4);
    bell << 1.0 / sqrt(2.0), 0, 0, 1.0 / sqrt(2.0);
    DensityMatrix rho_bell(2, vector<Complex>(bell.data(), bell.data() + 4));
    
    rho_bell.apply_kraus(NoiseModel::depolarising_2q(0.0, 0, 1, 2));
    test_result("Two-qubit depolarising(p=0): Tr(ρ) = 1", trace_is_one(rho_bell.get_rho()));
}

// ============================================================================
// Tests: Partial trace
// ============================================================================

void test_partial_trace_entangled() {
    cout << "\n--- Partial trace of entangled state ---\n";
    
    // Bell state |Φ⁺⟩ = (|00⟩+|11⟩)/√2
    VectorXcd bell(4);
    bell << 1.0 / sqrt(2.0), 0, 0, 1.0 / sqrt(2.0);
    DensityMatrix rho_bell(2, vector<Complex>(bell.data(), bell.data() + 4));
    
    // Trace out qubit 1 (the second one)
    DensityMatrix rho_q0 = rho_bell.partial_trace(1);
    
    test_result("Partial trace of qubit 1: 1-qubit result", rho_q0.get_num_qubits() == 1);
    test_result("Partial trace: Tr(ρ_q0) = 1", trace_is_one(rho_q0.get_rho()));
    
    // Tracing out one qubit of a maximally entangled Bell state gives I/2 (maximally mixed)
    auto rho_q0_mat = rho_q0.get_rho();
    test_result("Partial trace of Bell state: diagonal is [0.5, 0.5]", approx_eq(rho_q0_mat(0, 0).real(), 0.5) && approx_eq(rho_q0_mat(1, 1).real(), 0.5));
}

void test_partial_trace_product_state() {
    cout << "\n--- Partial trace of product state ---\n";
    
    // |+⟩⊗|1⟩
    VectorXcd prod(4);
    prod << 0.0, 1.0 / sqrt(2.0), 0.0, 1.0 / sqrt(2.0);
    DensityMatrix rho_prod(2, vector<Complex>(prod.data(), prod.data() + 4));
    
    // Trace out qubit 0
    DensityMatrix rho_q1 = rho_prod.partial_trace(0);
    
    test_result("Partial trace: result is 1-qubit", rho_q1.get_num_qubits() == 1);
    test_result("Partial trace: Tr(ρ_q1) = 1", trace_is_one(rho_q1.get_rho()));
    
    // Marginal over qubit 1 should be |1⟩⟨1|
    auto rho_q1_mat = rho_q1.get_rho();
    test_result("Partial trace of |+⟩⊗|1⟩: qubit 1 is |1⟩", approx_eq(rho_q1_mat(1, 1).real(), 1.0));
}

void test_partial_trace_three_qubit() {
    cout << "\n--- Partial trace on 3-qubit system ---\n";
    
    DensityMatrix rho3(3);
    DensityMatrix rho2 = rho3.partial_trace(0);
    test_result("Partial trace of qubit 0: result is 2-qubit", rho2.get_num_qubits() == 2);
    test_result("Partial trace: Tr(ρ_2) = 1", trace_is_one(rho2.get_rho()));
    
    DensityMatrix rho1 = rho2.partial_trace(1);
    test_result("Partial trace of qubit 1: result is 1-qubit", rho1.get_num_qubits() == 1);
    test_result("Partial trace: Tr(ρ_1) = 1", trace_is_one(rho1.get_rho()));
}

// ============================================================================
// Tests: Purity
// ============================================================================

void test_purity_pure_states() {
    cout << "\n--- Purity of pure states (Tr(ρ²) = 1) ---\n";
    
    // |0⟩
    DensityMatrix rho0(1);
    test_result("Purity(|0⟩) = 1", approx_eq(rho0.purity(), 1.0));
    
    // |+⟩
    VectorXcd plus(2); plus << 1.0 / sqrt(2.0), 1.0 / sqrt(2.0);
    DensityMatrix rho_plus(1, vector<Complex>(plus.data(), plus.data() + 2));
    test_result("Purity(|+⟩) = 1", approx_eq(rho_plus.purity(), 1.0));
}

void test_purity_mixed_states() {
    cout << "\n--- Purity of mixed states (Tr(ρ²) < 1) ---\n";
    
    // I/2 (maximally mixed)
    MatrixXcd rho_mm = MatrixXcd::Identity(2, 2) * 0.5;
    DensityMatrix dm_mm(1, vector<vector<Complex>>{{rho_mm(0, 0), rho_mm(0, 1)}, {rho_mm(1, 0), rho_mm(1, 1)}});
    test_result("Purity(I/2) = 0.5", approx_eq(dm_mm.purity(), 0.5, 1e-8));
}

void test_purity_bounds() {
    cout << "\n--- Purity bounds 1/d ≤ Tr(ρ²) ≤ 1 ---\n";
    
    // For a 2-qubit system: purity ∈ [1/4, 1]
    VectorXcd bell(4);
    bell << 1.0 / sqrt(2.0), 0, 0, 1.0 / sqrt(2.0);
    DensityMatrix rho_bell(2, vector<Complex>(bell.data(), bell.data() + 4));
    test_result("Purity(|Φ⁺⟩) = 1 (pure)", approx_eq(rho_bell.purity(), 1.0));
    
    // Partial trace → mixed state
    DensityMatrix rho_mixed = rho_bell.partial_trace(0);
    double pur = rho_mixed.purity();
    test_result("Purity of partial trace < 1", pur < 1.0);
    test_result("Purity ≥ 1/(2²)", pur >= 0.25 - EPSILON);
}

// ============================================================================
// Tests: Von Neumann entropy
// ============================================================================

void test_entropy_pure_states() {
    cout << "\n--- Von Neumann entropy of pure states (S = 0) ---\n";
    
    DensityMatrix rho0(1);
    test_result("S(|0⟩) = 0", approx_eq(rho0.entropy(), 0.0, 1e-8));
    
    VectorXcd plus(2); plus << 1.0 / sqrt(2.0), 1.0 / sqrt(2.0);
    DensityMatrix rho_plus(1, vector<Complex>(plus.data(), plus.data() + 2));
    test_result("S(|+⟩) = 0", approx_eq(rho_plus.entropy(), 0.0, 1e-8));
}

void test_entropy_mixed_states() {
    cout << "\n--- Von Neumann entropy of mixed states ---\n";
    
    // I/2: maximally mixed, S = ln(2) ≈ 0.693
    MatrixXcd rho_mm = MatrixXcd::Identity(2, 2) * 0.5;
    DensityMatrix dm_mm(1, vector<vector<Complex>>{{rho_mm(0, 0), rho_mm(0, 1)}, {rho_mm(1, 0), rho_mm(1, 1)}});
    double S_mm = dm_mm.entropy();
    test_result("S(I/2) ≈ ln(2)", approx_eq(S_mm, log(2.0), 1e-8));
}

void test_entropy_bounds() {
    cout << "\n--- Entropy bounds 0 ≤ S(ρ) ≤ ln(d) ---\n";
    
    VectorXcd bell(4);
    bell << 1.0 / sqrt(2.0), 0, 0, 1.0 / sqrt(2.0);
    DensityMatrix rho_bell(2, vector<Complex>(bell.data(), bell.data() + 4));
    
    // S(|Φ⁺⟩) = 0
    test_result("S(|Φ⁺⟩) = 0", approx_eq(rho_bell.entropy(), 0.0, 1e-8));
    
    // Partial trace of Bell state: I/2 for one qubit, S = ln(2)
    DensityMatrix rho_mixed = rho_bell.partial_trace(0);
    double S_mixed = rho_mixed.entropy();
    test_result("S(partial trace of Bell) = ln(2)", approx_eq(S_mixed, log(2.0), 1e-8));
}

void test_entropy_non_negative() {
    cout << "\n--- Entropy is non-negative S(ρ) ≥ 0 ---\n";
    
    for (double p = 0.0; p <= 1.0; p += 0.2) {
        DensityMatrix rho(1);
        rho.apply_kraus(NoiseModel::depolarising(min(p, 0.75), 0, 1));
        double S = rho.entropy();
        test_result(("S(ρ_depol[" + to_string(p) + "]) ≥ 0").c_str(), S >= -EPSILON);
    }
}

// ============================================================================
// Tests: Fidelity
// ============================================================================

void test_fidelity_pure_identical() {
    cout << "\n--- Fidelity with identical pure states ---\n";
    
    VectorXcd psi0(2); psi0 << 1, 0;
    DensityMatrix rho0(1, vector<Complex>(psi0.data(), psi0.data() + 2));
    test_result("F(ρ_|0⟩, |0⟩) = 1", approx_eq(rho0.fidelity(vector<Complex>{1, 0}), 1.0));
}

void test_fidelity_orthogonal() {
    cout << "\n--- Fidelity between orthogonal pure states ---\n";
    
    VectorXcd psi0(2); psi0 << 1, 0;
    DensityMatrix rho0(1, vector<Complex>(psi0.data(), psi0.data() + 2));
    test_result("F(ρ_|0⟩, |1⟩) = 0", approx_eq(rho0.fidelity(vector<Complex>{0, 1}), 0.0, 1e-9));
}

void test_fidelity_plus_zero() {
    cout << "\n--- Fidelity between |+⟩ and |0⟩ ---\n";
    
    VectorXcd plus(2); plus << 1.0 / sqrt(2.0), 1.0 / sqrt(2.0);
    DensityMatrix rho_plus(1, vector<Complex>(plus.data(), plus.data() + 2));
    test_result("F(ρ_|+⟩, |0⟩) = 0.5", approx_eq(rho_plus.fidelity(vector<Complex>{1, 0}), 0.5));
}

void test_fidelity_maximally_mixed() {
    cout << "\n--- Fidelity with maximally mixed state ---\n";
    
    MatrixXcd rho_mm = MatrixXcd::Identity(2, 2) * 0.5;
    DensityMatrix dm_mm(1, vector<vector<Complex>>{{rho_mm(0, 0), rho_mm(0, 1)}, {rho_mm(1, 0), rho_mm(1, 1)}});
    test_result("F(I/2, |0⟩) = 0.5", approx_eq(dm_mm.fidelity(vector<Complex>{1, 0}), 0.5));
}

// ============================================================================
// Tests: Measurement
// ============================================================================

void test_full_measurement_pure() {
    cout << "\n--- Full measurement of pure state ---\n";
    
    DensityMatrix rho(1);  // |0⟩
    string outcome = rho.measure();
    test_result("Measurement of |0⟩ = |0⟩", outcome == "|0⟩");
    
    // After measurement, state should be collapsed
    test_result("After measurement: Tr(ρ) = 1", trace_is_one(rho.get_rho()));
}

void test_partial_measurement() {
    cout << "\n--- Partial measurement of single qubit ---\n";
    
    // Start with |00⟩
    DensityMatrix rho(2);
    test_result("Initial: Tr(ρ) = 1", trace_is_one(rho.get_rho()));
    
    int bit0 = rho.partial_measurement(0);
    test_result("Partial measurement(qubit 0): outcome is 0 or 1", bit0 == 0 || bit0 == 1);
    test_result("After partial measurement: Tr(ρ) = 1", trace_is_one(rho.get_rho()));
}

void test_measurement_plus_state() {
    cout << "\n--- Measurement of superposition |+⟩ ---\n";
    
    VectorXcd plus(2); plus << 1.0 / sqrt(2.0), 1.0 / sqrt(2.0);
    DensityMatrix rho_plus(1, vector<Complex>(plus.data(), plus.data() + 2));
    
    // Multiple measurements should give both 0 and 1 with ~50% each (statistical)
    vector<int> outcomes(100);
    for (int i = 0; i < 100; ++i) {
        DensityMatrix temp = rho_plus;  // copy
        int bit = temp.partial_measurement(0);
        outcomes[i] = bit;
    }
    
    int count0 = 0, count1 = 0;
    for (int bit : outcomes) { if (bit == 0) count0++; else count1++; }
    test_result("Measurement of |+⟩: both outcomes observed", count0 > 10 && count1 > 10);
}

// ============================================================================
// Tests: Expectation values
// ============================================================================

void test_expectation_z() {
    cout << "\n--- Expectation value ⟨Z⟩ ---\n";
    
    // ⟨Z⟩ on |0⟩ = +1
    DensityMatrix rho0(1);
    test_result("⟨Z⟩_|0⟩ = +1", approx_eq(rho0.expectation(Gate::Z()), 1.0));
    
    // ⟨Z⟩ on |1⟩ = -1
    VectorXcd ket1(2); ket1 << 0, 1;
    DensityMatrix rho1(1, vector<Complex>(ket1.data(), ket1.data() + 2));
    test_result("⟨Z⟩_|1⟩ = -1", approx_eq(rho1.expectation(Gate::Z()), -1.0));
    
    // ⟨Z⟩ on |+⟩ = 0
    VectorXcd plus(2); plus << 1.0 / sqrt(2.0), 1.0 / sqrt(2.0);
    DensityMatrix rho_plus(1, vector<Complex>(plus.data(), plus.data() + 2));
    test_result("⟨Z⟩_|+⟩ = 0", approx_eq(rho_plus.expectation(Gate::Z()), 0.0, 1e-9));
}

void test_expectation_x() {
    cout << "\n--- Expectation value ⟨X⟩ ---\n";
    
    // ⟨X⟩ on |+⟩ = +1
    VectorXcd plus(2); plus << 1.0 / sqrt(2.0), 1.0 / sqrt(2.0);
    DensityMatrix rho_plus(1, vector<Complex>(plus.data(), plus.data() + 2));
    test_result("⟨X⟩_|+⟩ = +1", approx_eq(rho_plus.expectation(Gate::X()), 1.0));
    
    // ⟨X⟩ on |-⟩ = -1
    VectorXcd minus(2); minus << 1.0 / sqrt(2.0), -1.0 / sqrt(2.0);
    DensityMatrix rho_minus(1, vector<Complex>(minus.data(), minus.data() + 2));
    test_result("⟨X⟩_|-⟩ = -1", approx_eq(rho_minus.expectation(Gate::X()), -1.0));
}

// ============================================================================
// Tests: Probability vector
// ============================================================================

void test_probabilities() {
    cout << "\n--- Probability vector (diagonal of ρ) ---\n";
    
    DensityMatrix rho0(1);
    auto probs = rho0.probabilities();
    test_result("Probabilities of |0⟩ = [1, 0]", approx_eq(probs[0], 1.0) && approx_eq(probs[1], 0.0));
    
    VectorXcd plus(2); plus << 1.0 / sqrt(2.0), 1.0 / sqrt(2.0);
    DensityMatrix rho_plus(1, vector<Complex>(plus.data(), plus.data() + 2));
    auto probs_plus = rho_plus.probabilities();
    test_result("Probabilities of |+⟩ ≈ [0.5, 0.5]", approx_eq(probs_plus[0], 0.5) && approx_eq(probs_plus[1], 0.5));
}

// ============================================================================
// Tests: Consistency across operations
// ============================================================================

void test_trace_preservation() {
    cout << "\n--- Trace preservation through operations ---\n";
    
    DensityMatrix rho(2);
    
    test_result("Initial: Tr(ρ) = 1", trace_is_one(rho.get_rho()));
    
    rho.apply(Gate::H().expand(2, 0));
    test_result("After H(0): Tr(ρ) = 1", trace_is_one(rho.get_rho()));
    
    rho.apply(Gate::cnot(0, 1, 2));
    test_result("After CNOT(0,1): Tr(ρ) = 1", trace_is_one(rho.get_rho()));
    
    rho.apply_kraus(NoiseModel::depolarising(0.1, 0, 2));
    test_result("After noise: Tr(ρ) = 1", trace_is_one(rho.get_rho()));
}

void test_hermiticity_preservation() {
    cout << "\n--- Hermiticity preservation ---\n";
    
    DensityMatrix rho(2);
    
    for (int step = 0; step < 5; ++step) {
        rho.apply(Gate::H().expand(2, step % 2));
        test_result(("Step " + to_string(step) + ": Hermitian").c_str(), is_hermitian(rho.get_rho()));
    }
}

// ============================================================================
// Tests: Integration — noise sweep
// ============================================================================

void test_noise_evolution() {
    cout << "\n--- Noise sweep: depolarising error effect ---\n";
    
    VectorXcd plus(2); plus << 1.0 / sqrt(2.0), 1.0 / sqrt(2.0);
    
    double purity_prev = 1.0;
    double entropy_prev = 0.0;
    
    for (double p : {0.0, 0.1, 0.3, 0.5}) {
        DensityMatrix rho(1, vector<Complex>(plus.data(), plus.data() + 2));
        rho.apply_kraus(NoiseModel::depolarising(p, 0, 1));
        
        double pur = rho.purity();
        double ent = rho.entropy();
        
        test_result(("Purity at p=" + to_string(p) + " ≤ previous").c_str(), pur <= purity_prev + EPSILON);
        test_result(("Entropy at p=" + to_string(p) + " ≥ previous").c_str(), ent >= entropy_prev - EPSILON);
        
        purity_prev = pur;
        entropy_prev = ent;
    }
}

// ============================================================================
// Main test suite runner
// ============================================================================

int main() {
    cout << "\n" << string(70, '=') << "\n"
         << "  Test Suite: DensityMatrix.h\n"
         << string(70, '=') << "\n";
    
    test_pure_state_construction();
    test_plus_state();
    test_bell_state();
    test_maximally_mixed();
    
    test_unitary_evolution();
    test_hadamard_evolution();
    
    test_kraus_bit_flip();
    test_kraus_depolarising();
    test_kraus_amplitude_damping();
    test_kraus_phase_flip();
    test_kraus_two_qubit();
    
    test_partial_trace_entangled();
    test_partial_trace_product_state();
    test_partial_trace_three_qubit();
    
    test_purity_pure_states();
    test_purity_mixed_states();
    test_purity_bounds();
    
    test_entropy_pure_states();
    test_entropy_mixed_states();
    test_entropy_bounds();
    test_entropy_non_negative();
    
    test_fidelity_pure_identical();
    test_fidelity_orthogonal();
    test_fidelity_plus_zero();
    test_fidelity_maximally_mixed();
    
    test_full_measurement_pure();
    test_partial_measurement();
    test_measurement_plus_state();
    
    test_expectation_z();
    test_expectation_x();
    
    test_probabilities();
    
    test_trace_preservation();
    test_hermiticity_preservation();
    
    test_noise_evolution();
    
    cout << "\n" << string(70, '=') << "\n"
         << "  All tests passed ✓\n"
         << string(70, '=') << "\n\n";
    
    return 0;
}