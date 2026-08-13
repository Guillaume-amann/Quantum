// =============================================================================
//  test_noisemodel.cpp — Unit tests for NoiseModel.h
//
//  Verifies:
//    * All Kraus channel factories produce valid operators
//    * Completeness relation: Σ K†K = I (verified at runtime)
//    * Physical parameter bounds enforced (T2 ≤ 2T1, p ∈ [0,0.75], etc.)
//    * Channel effects on density matrices (unitary or dissipative)
//    * Noise escalation: p=0 → identity, p→max → maximally mixed
//    * Two-qubit depolarising 16-operator structure
//    * Measurement error confusion matrix application
//    * T1/T2 combined channel decomposition
//    * Generalised amplitude damping at finite temperature
//
//  Build:
//    g++ -std=c++17 -O2 -I/opt/homebrew/include/eigen3 test_noisemodel.cpp -o test_noisemodel
//  Run:
//    ./test_noisemodel
// =============================================================================

#include <iostream>
#include <cmath>
#include <complex>
#include <cassert>
#include <iomanip>
#include <vector>
#include "src/core/Gate.h"
#include "src/core/DensityMatrix.h"
#include "src/core/NoiseModel.h"

using namespace std;
using Complex = complex<double>;

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

// Check Hermiticity: M = M†
bool is_hermitian(const Eigen::MatrixXcd& M, double tol = EPSILON) {
    return (M - M.adjoint()).norm() < tol;
}

// Check trace is 1
bool trace_is_one(const Eigen::MatrixXcd& M, double tol = EPSILON) {
    return approx_eq(M.trace().real(), 1.0, tol) && approx_eq(M.trace().imag(), 0.0, tol);
}

// Print test result
void test_result(const char* name, bool passed) {
    cout << (passed ? "✓ PASS" : "✗ FAIL") << " : " << name << "\n";
    if (!passed) exit(1);
}

// ============================================================================
// Tests: Bit-flip channel
// ============================================================================

void test_bit_flip_identity() {
    cout << "\n--- Bit-flip channel with p=0 ---\n";
    
    auto kraus_ops = NoiseModel::bit_flip(0.0, 0, 1);
    test_result("bit_flip(0) returns 2 Kraus operators", kraus_ops.size() == 2);
    
    DensityMatrix rho(1);  // |0⟩
    rho.apply_kraus(kraus_ops);
    
    auto rho_mat = rho.get_rho();
    test_result("After bit_flip(0): Tr(ρ) = 1", trace_is_one(rho_mat));
    test_result("After bit_flip(0): ρ[0][0] = 1", approx_eq(rho_mat(0, 0).real(), 1.0));
}

void test_bit_flip_total_flip() {
    cout << "\n--- Bit-flip channel with p=1 ---\n";
    
    auto kraus_ops = NoiseModel::bit_flip(1.0, 0, 1);
    
    DensityMatrix rho(1);  // |0⟩
    rho.apply_kraus(kraus_ops);
    
    auto rho_mat = rho.get_rho();
    test_result("After bit_flip(1): |0⟩→|1⟩", approx_eq(rho_mat(1, 1).real(), 1.0));
}

void test_bit_flip_superposition() {
    cout << "\n--- Bit-flip on superposition ---\n";
    
    vector<Complex> plus_state(2);
    plus_state[0] = 1.0 / sqrt(2.0);
    plus_state[1] = 1.0 / sqrt(2.0);
    DensityMatrix rho(1, plus_state);
    
    auto kraus_ops = NoiseModel::bit_flip(0.5, 0, 1);
    rho.apply_kraus(kraus_ops);
    
    test_result("After bit_flip(0.5)|+⟩: Tr(ρ) = 1", trace_is_one(rho.get_rho()));
    test_result("After bit_flip(0.5)|+⟩: ρ is Hermitian", is_hermitian(rho.get_rho()));
}

// ============================================================================
// Tests: Phase-flip (dephasing) channel
// ============================================================================

void test_phase_flip_coherence_decay() {
    cout << "\n--- Phase-flip channel (dephasing) ---\n";
    
    vector<Complex> plus_state(2);
    plus_state[0] = 1.0 / sqrt(2.0);
    plus_state[1] = 1.0 / sqrt(2.0);
    DensityMatrix rho(1, plus_state);
    
    auto rho_before = rho.get_rho();
    double coherence_before = abs(rho_before(0, 1));
    
    // Apply phase-flip: p=0.5 → partial dephasing
    rho.apply_kraus(NoiseModel::phase_flip(0.5, 0, 1));
    
    auto rho_after = rho.get_rho();
    double coherence_after = abs(rho_after(0, 1));
    
    test_result("After phase_flip(0.5): coherence shrinks", coherence_after < coherence_before);
    test_result("After phase_flip(0.5): diagonal unchanged", 
                approx_eq(rho_after(0, 0).real(), 0.5) && approx_eq(rho_after(1, 1).real(), 0.5));
}

void test_phase_flip_total_dephasing() {
    cout << "\n--- Phase-flip complete dephasing (p=0.5) ---\n";
    
    vector<Complex> plus_state(2);
    plus_state[0] = 1.0 / sqrt(2.0);
    plus_state[1] = 1.0 / sqrt(2.0);
    DensityMatrix rho(1, plus_state);
    
    // At p=0.5, phase-flip completely dephases (off-diagonals → 0)
    rho.apply_kraus(NoiseModel::phase_flip(0.5, 0, 1));
    
    auto rho_mat = rho.get_rho();
    test_result("After phase_flip(0.5): off-diagonals ≈ 0", 
                abs(rho_mat(0, 1)) < EPSILON && abs(rho_mat(1, 0)) < EPSILON);
}

// ============================================================================
// Tests: Depolarising channel
// ============================================================================

void test_depolarising_identity() {
    cout << "\n--- Depolarising channel with p=0 ---\n";
    
    auto kraus_ops = NoiseModel::depolarising(0.0, 0, 1);
    test_result("depolarising(0) returns 4 Kraus operators", kraus_ops.size() == 4);
    
    DensityMatrix rho(1);  // |0⟩
    rho.apply_kraus(kraus_ops);
    
    auto rho_mat = rho.get_rho();
    test_result("After depolarising(0): ρ unchanged", approx_eq(rho_mat(0, 0).real(), 1.0));
}

void test_depolarising_maximally_mixed() {
    cout << "\n--- Depolarising channel with p=0.75 (maximum) ---\n";
    
    DensityMatrix rho(1);  // |0⟩
    rho.apply_kraus(NoiseModel::depolarising(0.75, 0, 1));
    
    auto rho_mat = rho.get_rho();
    test_result("After depolarising(0.75): ρ → I/2", 
                approx_eq(rho_mat(0, 0).real(), 0.5) && approx_eq(rho_mat(1, 1).real(), 0.5));
    test_result("After depolarising(0.75): off-diagonals = 0", 
                approx_eq(rho_mat(0, 1).real(), 0.0));
}

void test_depolarising_superposition() {
    cout << "\n--- Depolarising on superposition ---\n";
    
    vector<Complex> plus_state(2);
    plus_state[0] = 1.0 / sqrt(2.0);
    plus_state[1] = 1.0 / sqrt(2.0);
    DensityMatrix rho(1, plus_state);
    
    rho.apply_kraus(NoiseModel::depolarising(0.3, 0, 1));
    
    auto rho_mat = rho.get_rho();
    test_result("After depolarising(0.3)|+⟩: Tr(ρ) = 1", trace_is_one(rho_mat));
    // Diagonal should stay symmetric at 0.5, but coherence shrinks
    double coherence = abs(rho_mat(0, 1));
    test_result("After depolarising(0.3)|+⟩: coherence shrinks", coherence < 0.5);
}

// ============================================================================
// Tests: Amplitude damping (T1 relaxation)
// ============================================================================

void test_amplitude_damping_excited_state() {
    cout << "\n--- Amplitude damping: |1⟩ decays to |0⟩ ---\n";
    
    vector<Complex> ket1 = {0.0, 1.0};
    DensityMatrix rho(1, vector<Complex>(ket1.data(), ket1.data() + 2));
    
    // T1=1, t_gate=1 → γ = 1 - exp(-1) ≈ 0.632
    rho.apply_kraus(NoiseModel::amplitude_damping(1.0, 1.0, 0, 1));
    
    auto rho_mat = rho.get_rho();
    test_result("After amplitude_damping: |1⟩ population decreases", rho_mat(1, 1).real() < 1.0);
    test_result("After amplitude_damping: |0⟩ population increases", rho_mat(0, 0).real() > 0.0);
    test_result("After amplitude_damping: Tr(ρ) = 1", trace_is_one(rho_mat));
}

void test_amplitude_damping_ground_state() {
    cout << "\n--- Amplitude damping: |0⟩ unchanged ---\n";
    
    DensityMatrix rho(1);  // |0⟩
    rho.apply_kraus(NoiseModel::amplitude_damping(1.0, 1.0, 0, 1));
    
    auto rho_mat = rho.get_rho();
    test_result("After amplitude_damping |0⟩: ρ[0][0] = 1", approx_eq(rho_mat(0, 0).real(), 1.0));
}

void test_amplitude_damping_no_time() {
    cout << "\n--- Amplitude damping with t_gate=0 ---\n";
    
    vector<Complex> ket1 = {0.0, 1.0};
    DensityMatrix rho(1, vector<Complex>(ket1.data(), ket1.data() + 2));
    
    // t_gate=0 → γ = 0 (no time for decay)
    rho.apply_kraus(NoiseModel::amplitude_damping(1.0, 0.0, 0, 1));
    
    auto rho_mat = rho.get_rho();
    test_result("After amplitude_damping(t=0): |1⟩ unchanged", approx_eq(rho_mat(1, 1).real(), 1.0));
}

// ============================================================================
// Tests: Generalised amplitude damping (finite temperature)
// ============================================================================

void test_generalised_amplitude_damping_cold() {
    cout << "\n--- Generalised amplitude damping at T=0 (p_eq=0) ---\n";
    
    vector<Complex> ket1 = {0.0, 1.0};
    DensityMatrix rho(1, vector<Complex>(ket1.data(), ket1.data() + 2));
    
    // p_eq=0: pure decay (T=0 limit)
    rho.apply_kraus(NoiseModel::generalised_amplitude_damping(1.0, 1.0, 0.0, 0, 1));
    
    auto rho_mat = rho.get_rho();
    test_result("Generalised AD (T=0): |1⟩ decays", rho_mat(1, 1).real() < 1.0);
}

void test_generalised_amplitude_damping_hot() {
    cout << "\n--- Generalised amplitude damping at finite temperature (p_eq=0.3) ---\n";
    
    DensityMatrix rho(1);  // |0⟩
    
    // p_eq=0.3: finite thermal excitation
    rho.apply_kraus(NoiseModel::generalised_amplitude_damping(1.0, 1.0, 0.3, 0, 1));
    
    auto rho_mat = rho.get_rho();
    // |0⟩ can now be thermally excited to |1⟩
    test_result("Generalised AD (T>0): |1⟩ population grows", rho_mat(1, 1).real() > 0.0);
    test_result("Generalised AD: Tr(ρ) = 1", trace_is_one(rho_mat));
}

// ============================================================================
// Tests: T1/T2 combined channel
// ============================================================================

void test_t1_t2_physical_constraint() {
    cout << "\n--- T1/T2 physical constraint (T2 ≤ 2T1) ---\n";
    
    // Valid: T1=1, T2=2
    auto [ad, pd] = NoiseModel::t1_t2(1.0, 2.0, 1.0, 0, 1);
    test_result("T1/T2(T2=2T1): returns valid channels", ad.size() > 0 && pd.size() > 0);
    
    // T2 = 2T1 (limiting case) should work
    DensityMatrix rho(1);
    rho.apply_kraus(ad);
    test_result("T1/T2: amplitude damping applied without crash", trace_is_one(rho.get_rho()));
}

void test_t1_t2_invalid_constraint() {
    cout << "\n--- T1/T2 constraint violation detection ---\n";
    
    // Invalid: T2 > 2T1
    try {
        auto [ad, pd] = NoiseModel::t1_t2(1.0, 2.1, 1.0, 0, 1);  // T2 > 2T1
        test_result("T1/T2(T2 > 2T1): should throw", false);  // Should not reach here
    } catch (const invalid_argument&) {
        test_result("T1/T2(T2 > 2T1): throws invalid_argument", true);
    }
}

void test_t1_t2_combined_effect() {
    cout << "\n--- T1/T2 combined effect on |+⟩ ---\n";
    
    vector<Complex> plus_state(2);
    plus_state[0] = 1.0 / sqrt(2.0);
    plus_state[1] = 1.0 / sqrt(2.0);
    DensityMatrix rho(1, vector<Complex>(plus_state.data(), plus_state.data() + 2));
    
    auto [ad, pd] = NoiseModel::t1_t2(1.0, 1.5, 1.0, 0, 1);  // T2 < 2T1
    rho.apply_kraus(ad);
    rho.apply_kraus(pd);
    
    auto rho_mat = rho.get_rho();
    test_result("T1/T2 combined: coherence decays", abs(rho_mat(0, 1)) < 0.5);
    test_result("T1/T2 combined: Tr(ρ) = 1", trace_is_one(rho_mat));
}

// ============================================================================
// Tests: Two-qubit depolarising channel
// ============================================================================

void test_two_qubit_depolarising_identity() {
    cout << "\n--- Two-qubit depolarising with p=0 ---\n";
    
    auto kraus_ops = NoiseModel::depolarising_2q(0.0, 0, 1, 2);
    test_result("depolarising_2q(0) returns 16 operators", kraus_ops.size() == 16);
    
    VectorXcd bell(4); bell << 1.0 / sqrt(2.0), 0, 0, 1.0 / sqrt(2.0);
    DensityMatrix rho(2, vector<Complex>(bell.data(), bell.data() + 4));
    
    rho.apply_kraus(kraus_ops);
    test_result("After depolarising_2q(0): state unchanged", approx_eq(rho.purity(), 1.0));
}

void test_two_qubit_depolarising_mixed() {
    cout << "\n--- Two-qubit depolarising with p>0 ---\n";
    
    VectorXcd bell(4); bell << 1.0 / sqrt(2.0), 0, 0, 1.0 / sqrt(2.0);
    DensityMatrix rho(2, vector<Complex>(bell.data(), bell.data() + 4));
    
    auto kraus_ops = NoiseModel::depolarising_2q(0.3, 0, 1, 2);
    rho.apply_kraus(kraus_ops);
    
    double pur = rho.purity();
    test_result("After depolarising_2q(0.3): purity decreases", pur < 1.0);
    test_result("After depolarising_2q(0.3): Tr(ρ) = 1", trace_is_one(rho.get_rho()));
}

void test_two_qubit_different_qubits() {
    cout << "\n--- Two-qubit depolarising on different qubit pairs ---\n";
    
    DensityMatrix rho(3);  // Start with |000⟩
    
    // Apply noise to qubits (0,1) and (1,2)
    auto kraus1 = NoiseModel::depolarising_2q(0.1, 0, 1, 3);
    auto kraus2 = NoiseModel::depolarising_2q(0.1, 1, 2, 3);
    
    rho.apply_kraus(kraus1);
    rho.apply_kraus(kraus2);
    
    test_result("Two 2-qubit channels on 3-qubit system: Tr(ρ) = 1", trace_is_one(rho.get_rho()));
}

// ============================================================================
// Tests: Measurement error (classical confusion matrix)
// ============================================================================

void test_measurement_error_identity() {
    cout << "\n--- Measurement error with p=0 (no error) ---\n";
    
    auto me = NoiseModel::measurement_error(0.0, 0.0);
    vector<double> probs = {1.0, 0.0};
    
    auto noisy = me.apply(probs, 0, 1);
    test_result("measurement_error(0,0): probs unchanged", approx_eq(noisy[0], 1.0) && approx_eq(noisy[1], 0.0));
}

void test_measurement_error_bit_flip() {
    cout << "\n--- Measurement error with p_01=1 (report 1 when true is 0) ---\n";
    
    auto me = NoiseModel::measurement_error(1.0, 0.0);
    vector<double> probs = {1.0, 0.0};  // True outcome: |0⟩
    
    auto noisy = me.apply(probs, 0, 1);
    // With p_01=1: always report 1 when true is 0
    test_result("measurement_error(1,0)|0⟩: reported as |1⟩", approx_eq(noisy[1], 1.0));
}

void test_measurement_error_symmetric() {
    cout << "\n--- Measurement error symmetric confusion (p_01=p_10=0.05) ---\n";
    
    auto me = NoiseModel::measurement_error(0.05, 0.05);
    vector<double> probs = {0.5, 0.5};  // Equal superposition
    
    auto noisy = me.apply(probs, 0, 1);
    // With symmetric 5% error, distribution should remain approximately symmetric
    test_result("measurement_error symmetric: probs still ~[0.5, 0.5]", 
                approx_eq(noisy[0], 0.5, 0.1) && approx_eq(noisy[1], 0.5, 0.1));
}

void test_measurement_error_sample() {
    cout << "\n--- Measurement error sampling ---\n";
    
    auto me = NoiseModel::measurement_error(0.1, 0.1);
    vector<double> probs = {1.0, 0.0};  // |0⟩
    
    // Sample from noisy distribution
    int bit = me.sample(probs, 0, 1);
    test_result("measurement_error.sample: returns 0 or 1", bit == 0 || bit == 1);
}

// ============================================================================
// Tests: Physical parameter validation
// ============================================================================

void test_depolarising_p_bounds() {
    cout << "\n--- Depolarising parameter bounds ---\n";
    
    // Valid: p ∈ [0, 0.75]
    auto valid0 = NoiseModel::depolarising(0.0, 0, 1);
    test_result("depolarising(0): valid", valid0.size() == 4);
    
    auto valid75 = NoiseModel::depolarising(0.75, 0, 1);
    test_result("depolarising(0.75): valid", valid75.size() == 4);
    
    // Invalid: p > 0.75
    try {
        NoiseModel::depolarising(0.8, 0, 1);
        test_result("depolarising(0.8): should throw", false);
    } catch (const invalid_argument&) {
        test_result("depolarising(0.8): throws invalid_argument", true);
    }
}

void test_amplitude_damping_T1_positive() {
    cout << "\n--- Amplitude damping T1 > 0 ---\n";
    
    // Valid: T1 > 0
    auto valid = NoiseModel::amplitude_damping(1.0, 1.0, 0, 1);
    test_result("amplitude_damping(T1=1.0): valid", valid.size() == 2);
    
    // Invalid: T1 ≤ 0
    try {
        NoiseModel::amplitude_damping(0.0, 1.0, 0, 1);
        test_result("amplitude_damping(T1=0): should throw", false);
    } catch (const invalid_argument&) {
        test_result("amplitude_damping(T1=0): throws invalid_argument", true);
    }
}

void test_generalised_ad_p_eq_bounds() {
    cout << "\n--- Generalised amplitude damping p_eq ∈ [0, 0.5] ---\n";
    
    // Valid: p_eq ∈ [0, 0.5]
    auto valid0 = NoiseModel::generalised_amplitude_damping(1.0, 1.0, 0.0, 0, 1);
    test_result("generalised_ad(p_eq=0): valid", valid0.size() == 4);
    
    auto valid50 = NoiseModel::generalised_amplitude_damping(1.0, 1.0, 0.5, 0, 1);
    test_result("generalised_ad(p_eq=0.5): valid", valid50.size() == 4);
    
    // Invalid: p_eq > 0.5 (thermodynamically impossible)
    try {
        NoiseModel::generalised_amplitude_damping(1.0, 1.0, 0.6, 0, 1);
        test_result("generalised_ad(p_eq=0.6): should throw", false);
    } catch (const invalid_argument&) {
        test_result("generalised_ad(p_eq=0.6): throws invalid_argument", true);
    }
}

void test_measurement_error_p_bounds() {
    cout << "\n--- Measurement error p ∈ [0, 1] ---\n";
    
    // Valid
    auto valid = NoiseModel::measurement_error(0.5, 0.5);
    test_result("measurement_error(0.5, 0.5): valid", true);
    
    // Invalid: p_01 > 1
    try {
        NoiseModel::measurement_error(1.1, 0.5);
        test_result("measurement_error(1.1, 0.5): should throw", false);
    } catch (const invalid_argument&) {
        test_result("measurement_error(1.1, 0.5): throws invalid_argument", true);
    }
}

// ============================================================================
// Tests: Qubit index validation
// ============================================================================

void test_qubit_index_bounds() {
    cout << "\n--- Qubit index bounds checking ---\n";
    
    // Valid: qubit 0 in 2-qubit system
    auto valid = NoiseModel::bit_flip(0.1, 0, 2);
    test_result("bit_flip(target=0, n=2): valid", valid.size() > 0);
    
    // Invalid: qubit 2 in 2-qubit system
    try {
        NoiseModel::bit_flip(0.1, 2, 2);  // Out of bounds
        test_result("bit_flip(target=2, n=2): should throw", false);
    } catch (const invalid_argument&) {
        test_result("bit_flip(target=2, n=2): throws invalid_argument", true);
    }
}

void test_two_qubit_same_qubit() {
    cout << "\n--- Two-qubit gate on same qubit (invalid) ---\n";
    
    try {
        NoiseModel::depolarising_2q(0.1, 0, 0, 2);  // Same qubit
        test_result("depolarising_2q(0, 0, 2): should throw", false);
    } catch (const invalid_argument&) {
        test_result("depolarising_2q(0, 0, 2): throws invalid_argument", true);
    }
}

// ============================================================================
// Main test suite runner
// ============================================================================

int main() {
    cout << "\n" << string(70, '=') << "\n"
         << "  Test Suite: NoiseModel.h\n"
         << string(70, '=') << "\n";
    
    test_bit_flip_identity();
    test_bit_flip_total_flip();
    test_bit_flip_superposition();
    
    test_phase_flip_coherence_decay();
    test_phase_flip_total_dephasing();
    
    test_depolarising_identity();
    test_depolarising_maximally_mixed();
    test_depolarising_superposition();
    
    test_amplitude_damping_excited_state();
    test_amplitude_damping_ground_state();
    test_amplitude_damping_no_time();
    
    test_generalised_amplitude_damping_cold();
    test_generalised_amplitude_damping_hot();
    
    test_t1_t2_physical_constraint();
    test_t1_t2_invalid_constraint();
    test_t1_t2_combined_effect();
    
    test_two_qubit_depolarising_identity();
    test_two_qubit_depolarising_mixed();
    test_two_qubit_different_qubits();
    
    test_measurement_error_identity();
    test_measurement_error_bit_flip();
    test_measurement_error_symmetric();
    test_measurement_error_sample();
    
    test_depolarising_p_bounds();
    test_amplitude_damping_T1_positive();
    test_generalised_ad_p_eq_bounds();
    test_measurement_error_p_bounds();
    
    test_qubit_index_bounds();
    test_two_qubit_same_qubit();
    
    cout << "\n" << string(70, '=') << "\n"
         << "  All tests passed ✓\n"
         << string(70, '=') << "\n\n";
    
    return 0;
}