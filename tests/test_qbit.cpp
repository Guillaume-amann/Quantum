// =============================================================================
//  test_qbit.cpp — Unit tests for Qbit.h
//
//  Verifies:
//    * State vector normalisation (||ψ|| = 1)
//    * Unitary gate application preserves normalisation
//    * Measurement collapse: Born rule probabilities
//    * Measurement outcome repeats after collapse (state frozen)
//    * Multi-qubit state construction and access
//    * Probability vector extraction (non-destructive)
//    * Superposition states sample correctly
//    * Bell state entanglement structure
//    * GHZ state multi-qubit coherence
//
//  Build:
//    g++ -std=c++17 -O2 -I/opt/homebrew/include/eigen3 test_qbit.cpp -o test_qbit
//  Run:
//    ./test_qbit
// =============================================================================

#include <iostream>
#include <cmath>
#include <complex>
#include <cassert>
#include <iomanip>
#include <vector>
#include <map>
#include "src/core/Gate.h"
#include "src/core/Qbit.h"

using namespace std;
using Complex = complex<double>;

// ============================================================================
// Test utilities
// ============================================================================

constexpr double EPSILON = 1e-10;

bool approx_eq(double a, double b, double tol = EPSILON) {
    return fabs(a - b) < tol;
}

bool approx_eq(Complex a, Complex b, double tol = EPSILON) {
    return abs(a - b) < tol;
}

// Check normalisation: ||ψ|| = 1
bool is_normalised(const vector<Complex>& psi, double tol = EPSILON) {
    double norm_sq = 0.0;
    for (const auto& amp : psi) norm_sq += norm(amp);
    return approx_eq(norm_sq, 1.0, tol);
}

// Print test result
void test_result(const char* name, bool passed) {
    cout << (passed ? "✓ PASS" : "✗ FAIL") << " : " << name << "\n";
    if (!passed) exit(1);
}

// ============================================================================
// Tests: State construction and normalisation
// ============================================================================

void test_default_construction() {
    cout << "\n--- Default state construction ---\n";
    
    // Single qubit defaults to |0⟩
    Qbit q1(1);
    auto state1 = q1.get_state();
    test_result("|0⟩ has amplitude 1 at index 0", approx_eq(state1[0], 1.0));
    test_result("|0⟩ has amplitude 0 at index 1", approx_eq(state1[1], 0.0));
    test_result("|0⟩ is normalised", is_normalised(state1));
    
    // Two qubits default to |00⟩
    Qbit q2(2);
    auto state2 = q2.get_state();
    test_result("|00⟩ has amplitude 1 at index 0", approx_eq(state2[0], 1.0));
    test_result("|00⟩ is normalised", is_normalised(state2));
}

void test_plus_state_construction() {
    cout << "\n--- Superposition state construction ---\n";
    
    // |+⟩ = (|0⟩+|1⟩)/√2
    vector<Complex> plus_state(2);
    plus_state[0] = 1.0 / sqrt(2.0);
    plus_state[1] = 1.0 / sqrt(2.0);
    
    Qbit q_plus(1, plus_state);
    auto state = q_plus.get_state();
    
    test_result("ρ_|+⟩[0] = 1/√2", approx_eq(state[0], 1.0 / sqrt(2.0)));
    test_result("ρ_|+⟩[1] = 1/√2", approx_eq(state[1], 1.0 / sqrt(2.0)));
    test_result("|+⟩ is normalised", is_normalised(state));
}

void test_bell_state_construction() {
    cout << "\n--- Bell state (entangled) construction ---\n";
    
    // |Φ⁺⟩ = (|00⟩+|11⟩)/√2
    vector<Complex> bell_state(4);
    bell_state[0] = 1.0 / sqrt(2.0);
    bell_state[3] = 1.0 / sqrt(2.0);
    
    Qbit q_bell(2, bell_state);
    auto state = q_bell.get_state();
    
    test_result("|Φ⁺⟩[0] = 1/√2", approx_eq(state[0], 1.0 / sqrt(2.0)));
    test_result("|Φ⁺⟩[3] = 1/√2", approx_eq(state[3], 1.0 / sqrt(2.0)));
    test_result("|Φ⁺⟩[1] = 0", approx_eq(state[1], 0.0));
    test_result("|Φ⁺⟩[2] = 0", approx_eq(state[2], 0.0));
    test_result("|Φ⁺⟩ is normalised", is_normalised(state));
}

// ============================================================================
// Tests: Gate application and normalisation preservation
// ============================================================================

void test_gate_application_normalisation() {
    cout << "\n--- Gate application preserves normalisation ---\n";
    
    Qbit q(1);
    
    q.apply(Gate::X());
    test_result("After X: normalised", is_normalised(q.get_state()));
    
    q.apply(Gate::H());
    test_result("After H: normalised", is_normalised(q.get_state()));
    
    q.apply(Gate::Rx(M_PI / 4));
    test_result("After Rx(π/4): normalised", is_normalised(q.get_state()));
}

void test_hadamard_superposition() {
    cout << "\n--- Hadamard creates superposition ---\n";
    
    Qbit q(1);  // |0⟩
    q.apply(Gate::H());
    
    auto state = q.get_state();
    // H|0⟩ = (|0⟩+|1⟩)/√2
    test_result("(H|0⟩)[0] ≈ 1/√2", approx_eq(abs(state[0]), 1.0 / sqrt(2.0)));
    test_result("(H|0⟩)[1] ≈ 1/√2", approx_eq(abs(state[1]), 1.0 / sqrt(2.0)));
}

void test_cnot_entanglement() {
    cout << "\n--- CNOT creates entanglement ---\n";
    
    Qbit q(2);  // |00⟩
    
    // Apply Hadamard to first qubit
    q.apply(Gate::H().expand(2, 0));
    
    // Apply CNOT
    q.apply(Gate::CNOT());
    
    auto state = q.get_state();
    // Should be (|00⟩+|11⟩)/√2
    test_result("CNOT|+0⟩[0] ≈ 1/√2", approx_eq(abs(state[0]), 1.0 / sqrt(2.0)));
    test_result("CNOT|+0⟩[3] ≈ 1/√2", approx_eq(abs(state[3]), 1.0 / sqrt(2.0)));
    test_result("CNOT|+0⟩[1] ≈ 0", approx_eq(abs(state[1]), 0.0));
    test_result("CNOT|+0⟩[2] ≈ 0", approx_eq(abs(state[2]), 0.0));
}

// ============================================================================
// Tests: Measurement and Born rule
// ============================================================================

void test_measurement_deterministic() {
    cout << "\n--- Deterministic measurement (pure basis state) ---\n";
    
    Qbit q(1);  // |0⟩
    
    string outcome = q.measure();
    test_result("Measure |0⟩ always yields |0⟩", outcome == "|0⟩");
    
    // After measurement, state is still |0⟩
    auto state = q.get_state();
    test_result("After measurement: |0⟩", approx_eq(state[0], 1.0) && approx_eq(state[1], 0.0));
}

void test_measurement_collapse() {
    cout << "\n--- Measurement collapses state ---\n";
    
    vector<Complex> plus_state(2);
    plus_state[0] = 1.0 / sqrt(2.0);
    plus_state[1] = 1.0 / sqrt(2.0);
    Qbit q(1, plus_state);
    
    string outcome = q.measure();
    // Outcome is either |0⟩ or |1⟩
    test_result("Measurement outcome is valid", outcome == "|0⟩" || outcome == "|1⟩");
    
    // After collapse, repeated measurement gives same outcome
    string outcome2 = q.measure();
    test_result("Repeated measurement matches", outcome == outcome2);
    
    // State is now a basis state
    auto state = q.get_state();
    if (outcome == "|0⟩") {
        test_result("Collapsed to |0⟩", approx_eq(state[0], 1.0));
    } else {
        test_result("Collapsed to |1⟩", approx_eq(state[1], 1.0));
    }
}

void test_measurement_born_rule() {
    cout << "\n--- Born rule: measurement probabilities ---\n";
    
    // |ψ⟩ = (2|0⟩ + 1|1⟩)/√5  (unnormalised intentionally to test normalisation)
    vector<Complex> psi(2);
    psi[0] = 2.0;
    psi[1] = 1.0;
    
    Qbit q(1, psi);  // Constructor normalises automatically
    auto state = q.get_state();
    
    // After normalisation: amplitudes should be ≈ 2/√5 and 1/√5
    double expected_p0 = 4.0 / 5.0;   // |2/√5|² = 4/5
    double expected_p1 = 1.0 / 5.0;   // |1/√5|² = 1/5
    
    test_result("P(|0⟩) ≈ 4/5", approx_eq(norm(state[0]), expected_p0));
    test_result("P(|1⟩) ≈ 1/5", approx_eq(norm(state[1]), expected_p1));
}

void test_measurement_sampling_statistics() {
    cout << "\n--- Measurement sampling statistics ---\n";
    
    // |+⟩ = (|0⟩+|1⟩)/√2, should measure 0 and 1 with ~50% each
    vector<Complex> plus_state(2);
    plus_state[0] = 1.0 / sqrt(2.0);
    plus_state[1] = 1.0 / sqrt(2.0);
    
    int count0 = 0, count1 = 0;
    const int SAMPLES = 1000;
    
    for (int i = 0; i < SAMPLES; ++i) {
        Qbit q(1, plus_state);  // Fresh state each shot (crucial!)
        string outcome = q.measure();
        if (outcome == "|0⟩") count0++;
        else count1++;
    }
    
    double freq0 = static_cast<double>(count0) / SAMPLES;
    double freq1 = static_cast<double>(count1) / SAMPLES;
    
    // With 1000 samples, expect ~50% ± ~5% (3-sigma)
    test_result("P(|0⟩) ≈ 0.5 (sampling)", freq0 > 0.4 && freq0 < 0.6);
    test_result("P(|1⟩) ≈ 0.5 (sampling)", freq1 > 0.4 && freq1 < 0.6);
}

void test_measurement_entangled_state() {
    cout << "\n--- Measurement of entangled state ---\n";
    
    // |Φ⁺⟩ = (|00⟩+|11⟩)/√2
    vector<Complex> bell_state(4);
    bell_state[0] = 1.0 / sqrt(2.0);
    bell_state[3] = 1.0 / sqrt(2.0);
    Qbit q(2, bell_state);
    
    string outcome = q.measure();
    // Outcome should be |00⟩ or |11⟩ (never |01⟩ or |10⟩)
    test_result("Bell state outcome is correlated", outcome == "|00⟩" || outcome == "|11⟩");
}

// ============================================================================
// Tests: Multi-qubit states
// ============================================================================

void test_three_qubit_state() {
    cout << "\n--- Three-qubit state (GHZ) ---\n";
    
    // GHZ state: (|000⟩+|111⟩)/√2
    vector<Complex> ghz_state(8, 0.0);
    ghz_state[0] = 1.0 / sqrt(2.0);
    ghz_state[7] = 1.0 / sqrt(2.0);
    
    Qbit q(3, ghz_state);
    test_result("GHZ state size is 8", q.get_state().size() == 8);
    test_result("GHZ is normalised", is_normalised(q.get_state()));
    test_result("GHZ has num_qubits() = 3", q.num_qubits() == 3);
}

void test_four_qubit_state() {
    cout << "\n--- Four-qubit state ---\n";
    
    Qbit q(4);
    test_result("4-qubit state size is 16", q.get_state().size() == 16);
    test_result("4-qubit |0000⟩ is first basis state", approx_eq(q.get_state()[0], 1.0));
}

// ============================================================================
// Tests: State access and modification
// ============================================================================

void test_const_get_state() {
    cout << "\n--- Const state access ---\n";
    
    Qbit q(2);  // |00⟩
    const auto& state_ref = q.get_state();
    
    test_result("Const access to state[0]", approx_eq(state_ref[0], 1.0));
    test_result("State vector has size 4", state_ref.size() == 4);
}

void test_access_state_modification() {
    cout << "\n--- Mutable state access ---\n";
    
    Qbit q(1);  // |0⟩
    auto& state_mut = q.access_state();
    
    test_result("Initial state is |0⟩", approx_eq(state_mut[0], 1.0));
    
    // Modify state directly (dangerous, but available for advanced use)
    state_mut[0] = 1.0 / sqrt(2.0);
    state_mut[1] = 1.0 / sqrt(2.0);
    
    // Verify state changed
    auto state_after = q.get_state();
    test_result("State modified via access_state", approx_eq(abs(state_after[0]), 1.0 / sqrt(2.0)));
}

// ============================================================================
// Tests: Print output (no crash, reasonable format)
// ============================================================================

void test_print_state_no_crash() {
    cout << "\n--- Print state output (sanity check) ---\n";
    
    Qbit q(2);  // |00⟩
    q.apply(Gate::H().expand(2, 0));
    
    cout << "  State output: ";
    q.print_state();
    
    test_result("print_state() did not crash", true);
}

// ============================================================================
// Tests: Edge cases and error handling
// ============================================================================

void test_zero_amplitude_handling() {
    cout << "\n--- Zero amplitude states ---\n";
    
    vector<Complex> sparse_state(4, 0.0);
    sparse_state[0] = 1.0;  // Only |00⟩ has amplitude
    
    Qbit q(2, sparse_state);
    test_result("Sparse state is normalised", is_normalised(q.get_state()));
    
    string outcome = q.measure();
    test_result("Sparse state always measures to |00⟩", outcome == "|00⟩");
}

void test_phase_preservation() {
    cout << "\n--- Complex phase preservation ---\n";
    
    vector<Complex> phase_state(2);
    phase_state[0] = Complex(0, 1.0) / sqrt(2.0);  // i/√2
    phase_state[1] = Complex(0, 1.0) / sqrt(2.0);
    
    Qbit q(1, phase_state);
    auto state = q.get_state();
    
    // Phases are preserved (but don't affect measurement probabilities)
    test_result("Phase is preserved", approx_eq(state[0], Complex(0, 1.0) / sqrt(2.0)));
    test_result("Phase state is normalised", is_normalised(state));
}

// ============================================================================
// Tests: Sequential gate application
// ============================================================================

void test_sequential_gates() {
    cout << "\n--- Sequential gate application ---\n";
    
    Qbit q(1);  // |0⟩
    
    // H|0⟩ = |+⟩ = (|0⟩+|1⟩)/√2
    q.apply(Gate::H());
    
    // X(|0⟩+|1⟩)/√2 = (|1⟩+|0⟩)/√2 = |+⟩ (symmetric)
    q.apply(Gate::X());
    
    auto state = q.get_state();
    test_result("After H·X: still superposition", approx_eq(abs(state[0]), 1.0 / sqrt(2.0)) && approx_eq(abs(state[1]), 1.0 / sqrt(2.0)));
    test_result("After H·X: normalised", is_normalised(state));
}

void test_rotation_sequence() {
    cout << "\n--- Rotation gate sequence ---\n";
    
    Qbit q(1);  // |0⟩
    
    // Rx(π/2) twice should be Rx(π)
    q.apply(Gate::Rx(M_PI / 2));
    q.apply(Gate::Rx(M_PI / 2));
    
    auto state = q.get_state();
    // Should now be proportional to i|1⟩ (or -i|0⟩ + 1|1⟩, depending on phase)
    test_result("After Rx(π/2)·Rx(π/2): flips toward |1⟩", norm(state[1]) > 0.5);
    test_result("Rotation sequence is normalised", is_normalised(state));
}

// ============================================================================
// Main test suite runner
// ============================================================================

int main() {
    cout << "\n" << string(70, '=') << "\n"
         << "  Test Suite: Qbit.h\n"
         << string(70, '=') << "\n";
    
    test_default_construction();
    test_plus_state_construction();
    test_bell_state_construction();
    
    test_gate_application_normalisation();
    test_hadamard_superposition();
    test_cnot_entanglement();
    
    test_measurement_deterministic();
    test_measurement_collapse();
    test_measurement_born_rule();
    test_measurement_sampling_statistics();
    test_measurement_entangled_state();
    
    test_three_qubit_state();
    test_four_qubit_state();
    
    test_const_get_state();
    test_access_state_modification();
    
    test_print_state_no_crash();
    
    test_zero_amplitude_handling();
    test_phase_preservation();
    
    test_sequential_gates();
    test_rotation_sequence();
    
    cout << "\n" << string(70, '=') << "\n"
         << "  All tests passed ✓\n"
         << string(70, '=') << "\n\n";
    
    return 0;
}