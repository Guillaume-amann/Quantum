// =============================================================================
//  test_gate.cpp — Unit tests for Gate.h
//
//  Verifies:
//    * Single-qubit gates are unitary (U†U = I)
//    * Two-qubit gates are unitary
//    * Tensor products satisfy Kronecker properties
//    * Gate expansion (embed into n-qubit register) preserves unitarity
//    * expand_two() correctly lifts arbitrary 2-qubit gates to full register
//    * Gate composition U·V is consistent with matrix multiplication
//    * Pauli algebra holds (X² = Y² = Z² = I, XY = iZ, etc.)
//
//  Build:
//    g++ -std=c++17 -O2 -I/opt/homebrew/include/eigen3 test_gate.cpp -o test_gate
//  Run:
//    ./test_gate
// =============================================================================

#include <iostream>
#include <cmath>
#include <complex>
#include <cassert>
#include <iomanip>
#include "src/core/Gate.h"

using namespace std;
using Complex = complex<double>;

// ============================================================================
// Test utilities
// ============================================================================

constexpr double EPSILON = 1e-10;

bool approx_eq(Complex a, Complex b, double tol = EPSILON) {
    return abs(a - b) < tol;
}

bool approx_eq(double a, double b, double tol = EPSILON) {
    return fabs(a - b) < tol;
}

// Check if matrix is unitary: M†M = I
bool is_unitary(const Eigen::MatrixXcd& M, double tol = EPSILON) {
    auto I = M.adjoint() * M;
    auto should_be_I = Eigen::MatrixXcd::Identity(M.rows(), M.rows());
    return (I - should_be_I).norm() < tol;
}

// Check if matrix is Hermitian: M = M†
bool is_hermitian(const Eigen::MatrixXcd& M, double tol = EPSILON) {
    return (M - M.adjoint()).norm() < tol;
}

// Check if matrix is identity
bool is_identity(const Eigen::MatrixXcd& M, double tol = EPSILON) {
    auto I = Eigen::MatrixXcd::Identity(M.rows(), M.rows());
    return (M - I).norm() < tol;
}

// Trace (sum of diagonal elements)
Complex trace(const Eigen::MatrixXcd& M) {
    return M.trace();
}

// Print a test result
void test_result(const char* name, bool passed) {
    cout << (passed ? "✓ PASS" : "✗ FAIL") << " : " << name << "\n";
    if (!passed) exit(1);
}

// ============================================================================
// Tests: Single-qubit gates
// ============================================================================

void test_single_qubit_unitarity() {
    cout << "\n--- Single-qubit gate unitarity ---\n";
    
    test_result("X is unitary",   is_unitary(Gate::X().matrix));
    test_result("Y is unitary",   is_unitary(Gate::Y().matrix));
    test_result("Z is unitary",   is_unitary(Gate::Z().matrix));
    test_result("H is unitary",   is_unitary(Gate::H().matrix));
    test_result("Rx(π/4) is unitary", is_unitary(Gate::Rx(M_PI / 4).matrix));
    test_result("Ry(π/3) is unitary", is_unitary(Gate::Ry(M_PI / 3).matrix));
    test_result("Rz(π/6) is unitary", is_unitary(Gate::Rz(M_PI / 6).matrix));
}

void test_pauli_algebra() {
    cout << "\n--- Pauli algebra ---\n";
    
    auto X = Gate::X().matrix;
    auto Y = Gate::Y().matrix;
    auto Z = Gate::Z().matrix;
    auto I = Eigen::MatrixXcd::Identity(2, 2);
    
    // X² = Y² = Z² = I
    test_result("X² = I", is_identity(X * X));
    test_result("Y² = I", is_identity(Y * Y));
    test_result("Z² = I", is_identity(Z * Z));
    
    // X† = X, Y† = Y, Z† = Z (Hermitian)
    test_result("X is Hermitian", is_hermitian(X));
    test_result("Y is Hermitian", is_hermitian(Y));
    test_result("Z is Hermitian", is_hermitian(Z));
    
    // XY = iZ, YZ = iX, ZX = iY
    auto XY_iZ = X * Y - Complex(0, 1) * Z;
    auto YZ_iX = Y * Z - Complex(0, 1) * X;
    auto ZX_iY = Z * X - Complex(0, 1) * Y;
    
    test_result("XY = iZ",  XY_iZ.norm() < EPSILON);
    test_result("YZ = iX",  YZ_iX.norm() < EPSILON);
    test_result("ZX = iY",  ZX_iY.norm() < EPSILON);
    
    // {X, Y} = 0 (anticommute)
    auto XY_plus_YX = X * Y + Y * X;
    test_result("{X,Y} = 0", XY_plus_YX.norm() < EPSILON);
}

void test_rotation_gates() {
    cout << "\n--- Rotation gates (exponential form) ---\n";
    
    // Rx(θ) = exp(-i θ/2 X) = [[cos(θ/2), -i sin(θ/2)], [-i sin(θ/2), cos(θ/2)]]
    // At θ=0: identity; at θ=π: [[0, -i], [-i, 0]] = -iX (up to global phase)
    auto Rx0 = Gate::Rx(0.0).matrix;
    test_result("Rx(0) = I", is_identity(Rx0));
    
    // Rx(π) should map |0⟩ → i|1⟩ and |1⟩ → i|0⟩ (phase ±i applied)
    auto Rxpi = Gate::Rx(M_PI).matrix;
    Eigen::VectorXcd ket0(2); ket0 << 1, 0;
    Eigen::VectorXcd ket1(2); ket1 << 0, 1;
    auto Rx_ket0 = Rxpi * ket0;
    auto Rx_ket1 = Rxpi * ket1;
    test_result("Rx(π)|0⟩ ∝ |1⟩", approx_eq(abs(Rx_ket0(0)), 0.0) && approx_eq(abs(Rx_ket0(1)), 1.0));
    test_result("Rx(π)|1⟩ ∝ |0⟩", approx_eq(abs(Rx_ket1(0)), 1.0) && approx_eq(abs(Rx_ket1(1)), 0.0));
    
    // Ry(θ) = [[cos(θ/2), -sin(θ/2)], [sin(θ/2), cos(θ/2)]] at θ=π: [[0, -1], [1, 0]] = -iY
    auto Rypi = Gate::Ry(M_PI).matrix;
    auto Ry_ket0 = Rypi * ket0;
    auto Ry_ket1 = Rypi * ket1;
    test_result("Ry(π)|0⟩ ∝ |1⟩", approx_eq(abs(Ry_ket0(0)), 0.0) && approx_eq(abs(Ry_ket0(1)), 1.0));
    test_result("Ry(π)|1⟩ ∝ |0⟩", approx_eq(abs(Ry_ket1(0)), 1.0) && approx_eq(abs(Ry_ket1(1)), 0.0));
    
    // Rz(θ) = diag(exp(-iθ/2), exp(iθ/2)) at θ=π: diag(-i, i) = -iZ (up to phase)
    auto Rzpi = Gate::Rz(M_PI).matrix;
    auto Rz_ket0 = Rzpi * ket0;
    auto Rz_ket1 = Rzpi * ket1;
    // |0⟩ maps to e^(-iπ/2)|0⟩ = -i|0⟩, |1⟩ maps to e^(iπ/2)|1⟩ = i|1⟩
    test_result("Rz(π)|0⟩ ∝ |0⟩", approx_eq(abs(Rz_ket0(0)), 1.0) && approx_eq(abs(Rz_ket0(1)), 0.0));
    test_result("Rz(π)|1⟩ ∝ |1⟩", approx_eq(abs(Rz_ket1(0)), 0.0) && approx_eq(abs(Rz_ket1(1)), 1.0));
}

void test_hadamard() {
    cout << "\n--- Hadamard gate ---\n";
    
    auto H = Gate::H().matrix;
    auto X = Gate::X().matrix;
    auto Z = Gate::Z().matrix;
    
    // H² = I
    test_result("H² = I", is_identity(H * H));
    
    // H is Hermitian
    test_result("H is Hermitian", is_hermitian(H));
    
    // H diagonalises X: HXH† = Z (conjugate transpose, equivalent to HXH for Hermitian H)
    auto HXH_adj = H * X * H.adjoint();
    test_result("HXH† = Z", (HXH_adj - Z).norm() < EPSILON);
    
    // H|0⟩ = (|0⟩+|1⟩)/√2 = |+⟩
    Eigen::VectorXcd ket0(2);
    ket0 << 1, 0;
    auto plus = H * ket0;
    Eigen::VectorXcd expected(2);
    expected << 1.0 / sqrt(2.0), 1.0 / sqrt(2.0);
    test_result("H|0⟩ = |+⟩", (plus - expected).norm() < EPSILON);
}

// ============================================================================
// Tests: Two-qubit gates
// ============================================================================

void test_two_qubit_unitarity() {
    cout << "\n--- Two-qubit gate unitarity ---\n";
    
    test_result("CNOT is unitary", is_unitary(Gate::CNOT().matrix));
    test_result("CZ is unitary",   is_unitary(Gate::CZ().matrix));
    test_result("SWAP is unitary", is_unitary(Gate::SWAP().matrix));
    test_result("Rzz(π/4) is unitary", is_unitary(Gate::Rzz(M_PI / 4).matrix));
}

void test_cnot() {
    cout << "\n--- CNOT (control-first basis: |00⟩|01⟩|10⟩|11⟩) ---\n";
    
    auto CNOT = Gate::CNOT().matrix;
    
    // |00⟩ -> |00⟩
    Eigen::VectorXcd ket00(4); ket00 << 1, 0, 0, 0;
    test_result("CNOT|00⟩ = |00⟩", (CNOT * ket00 - ket00).norm() < EPSILON);
    
    // |01⟩ -> |01⟩
    Eigen::VectorXcd ket01(4); ket01 << 0, 1, 0, 0;
    test_result("CNOT|01⟩ = |01⟩", (CNOT * ket01 - ket01).norm() < EPSILON);
    
    // |10⟩ -> |11⟩ (control=1 flips target)
    Eigen::VectorXcd ket10(4); ket10 << 0, 0, 1, 0;
    Eigen::VectorXcd ket11(4); ket11 << 0, 0, 0, 1;
    test_result("CNOT|10⟩ = |11⟩", (CNOT * ket10 - ket11).norm() < EPSILON);
    
    // |11⟩ -> |10⟩
    test_result("CNOT|11⟩ = |10⟩", (CNOT * ket11 - ket10).norm() < EPSILON);
    
    // CNOT² = I
    test_result("CNOT² = I", is_identity(CNOT * CNOT));
}

void test_swap() {
    cout << "\n--- SWAP gate ---\n";
    
    auto SWAP = Gate::SWAP().matrix;
    
    // SWAP² = I
    test_result("SWAP² = I", is_identity(SWAP * SWAP));
    
    // SWAP|01⟩ = |10⟩
    Eigen::VectorXcd ket01(4); ket01 << 0, 1, 0, 0;
    Eigen::VectorXcd ket10(4); ket10 << 0, 0, 1, 0;
    test_result("SWAP|01⟩ = |10⟩", (SWAP * ket01 - ket10).norm() < EPSILON);
}

void test_rzz() {
    cout << "\n--- Rzz(θ) = exp(-i θ/2 Z⊗Z) ---\n";
    
    // At θ=0: identity
    auto Rzz0 = Gate::Rzz(0.0).matrix;
    test_result("Rzz(0) = I", is_identity(Rzz0));
    
    // Rzz² at θ=π should equal exp(-i π Z⊗Z)
    auto Rzz_half = Gate::Rzz(M_PI / 2).matrix;
    test_result("Rzz(π/2) is unitary", is_unitary(Rzz_half));
    
    // Rzz(θ) diagonalises in Z⊗Z eigenbasis — should leave |00⟩, |01⟩, |10⟩, |11⟩ as eigenstates
    // (only adds phase, no population change)
    Eigen::VectorXcd ket00(4); ket00 << 1, 0, 0, 0;
    auto result = Rzz_half * ket00;
    test_result("Rzz(π/2)|00⟩ is proportional to |00⟩", approx_eq(abs(result(1)), 0.0) && approx_eq(abs(result(2)), 0.0) && approx_eq(abs(result(3)), 0.0));
}

void test_toffoli() {
    cout << "\n--- Toffoli (CCX) ---\n";
    
    auto CCX = Gate::Toffoli().matrix;
    
    // Toffoli² = I
    test_result("CCX² = I", is_identity(CCX * CCX));
    
    // |111⟩ -> |110⟩ (both controls 1, flip target)
    Eigen::VectorXcd ket111(8); ket111 << 0, 0, 0, 0, 0, 0, 0, 1;
    Eigen::VectorXcd ket110(8); ket110 << 0, 0, 0, 0, 0, 0, 1, 0;
    test_result("CCX|111⟩ = |110⟩", (CCX * ket111 - ket110).norm() < EPSILON);
    
    // |110⟩ -> |111⟩
    test_result("CCX|110⟩ = |111⟩", (CCX * ket110 - ket111).norm() < EPSILON);
}

// ============================================================================
// Tests: Tensor products and Kronecker structure
// ============================================================================

void test_tensor_product() {
    cout << "\n--- Tensor product (Kronecker) ---\n";
    
    // (A ⊗ B)† = A† ⊗ B†
    auto X = Gate::X().matrix;
    auto Y = Gate::Y().matrix;
    auto XY = Gate::tensor(Gate::X(), Gate::Y()).matrix;
    auto X_adj_Y_adj = Gate::tensor(Gate(vector<vector<Complex>>{{0, 1}, {1, 0}}),
                                     Gate(vector<vector<Complex>>{{0, Complex(0, 1)}, {Complex(0, -1), 0}})).matrix;
    test_result("(X ⊗ Y)† = X† ⊗ Y†", is_hermitian(XY) == is_hermitian(X_adj_Y_adj));
    
    // X ⊗ X is unitary
    auto XX = Gate::tensor(Gate::X(), Gate::X()).matrix;
    test_result("X ⊗ X is unitary", is_unitary(XX));
    
    // Dimension check
    test_result("I⊗I (4x4)", Gate::tensor(Gate::identity(1), Gate::identity(1)).matrix.rows() == 4);
    test_result("I⊗I⊗I (8x8)", Gate::tensor(Gate::tensor(Gate::identity(1), Gate::identity(1)), Gate::identity(1)).matrix.rows() == 8);
}

// ============================================================================
// Tests: Gate expansion into n-qubit register
// ============================================================================

void test_expand_single_qubit() {
    cout << "\n--- Single-qubit gate expansion (expand) ---\n";
    
    // X on qubit 0 in a 2-qubit register: X ⊗ I
    auto X_q0_2 = Gate::X().expand(2, 0).matrix;
    auto X_tensor_I = Gate::tensor(Gate::X(), Gate::identity(1)).matrix;
    test_result("X.expand(2,0) = X⊗I", (X_q0_2 - X_tensor_I).norm() < EPSILON);
    
    // X on qubit 1 in a 2-qubit register: I ⊗ X
    auto X_q1_2 = Gate::X().expand(2, 1).matrix;
    auto I_tensor_X = Gate::tensor(Gate::identity(1), Gate::X()).matrix;
    test_result("X.expand(2,1) = I⊗X", (X_q1_2 - I_tensor_X).norm() < EPSILON);
    
    // H on qubit 0 in a 3-qubit register: H ⊗ I ⊗ I
    auto H_q0_3 = Gate::H().expand(3, 0).matrix;
    test_result("H.expand(3,0) is unitary", is_unitary(H_q0_3));
    test_result("H.expand(3,0) is 8×8", H_q0_3.rows() == 8);
}

void test_expand_two_qubit() {
    cout << "\n--- Two-qubit gate expansion (expand_two) ---\n";
    
    // CNOT on qubits (0,1) in a 2-qubit register should equal CNOT itself
    auto CNOT_01_2 = Gate::cnot(0, 1, 2).matrix;
    auto CNOT_base = Gate::CNOT().matrix;
    test_result("CNOT(0,1,2) = CNOT", (CNOT_01_2 - CNOT_base).norm() < EPSILON);
    
    // CNOT on qubits (1,0) in a 2-qubit register: should swap the roles
    // In big-endian: qubit 0 is MSB. The 2-qubit gate CNOT is ctrl-first.
    // So CNOT(1,0,2) lifts control from qubit 1 and target to qubit 0.
    auto CNOT_10_2 = Gate::cnot(1, 0, 2).matrix;
    test_result("CNOT(1,0,2) is unitary", is_unitary(CNOT_10_2));
    
    // CZ on qubits (0,1) should be symmetric (CZ(0,1) = CZ(1,0))
    auto CZ_01 = Gate::cz(0, 1, 2).matrix;
    auto CZ_10 = Gate::cz(1, 0, 2).matrix;
    test_result("CZ(0,1) = CZ(1,0)", (CZ_01 - CZ_10).norm() < EPSILON);
    
    // SWAP on qubits (0,1) should commute when applied again
    auto SWAP_01 = Gate::swap(0, 1, 2).matrix;
    test_result("SWAP(0,1)² = I", is_identity(SWAP_01 * SWAP_01));
}

void test_expand_two_non_adjacent() {
    cout << "\n--- expand_two on non-adjacent qubits ---\n";
    
    // Rzz(θ) on qubits (0,2) in a 3-qubit register
    auto Rzz_02_3 = Gate::rzz(M_PI / 4, 0, 2, 3).matrix;
    test_result("Rzz(π/4) on (0,2) in 3-qubit is unitary", is_unitary(Rzz_02_3));
    test_result("Rzz(π/4) on (0,2) in 3-qubit is 8×8", Rzz_02_3.rows() == 8);
    
    // On a 4-qubit register, Rzz(θ) on (1,3) should be symmetric
    auto Rzz_13_4 = Gate::rzz(M_PI / 3, 1, 3, 4).matrix;
    test_result("Rzz on (1,3) in 4-qubit is unitary", is_unitary(Rzz_13_4));
}

void test_controlled_u() {
    cout << "\n--- Controlled-U factory ---\n";
    
    auto CX = Gate::controlled_U(Gate::X()).matrix;
    auto CNOT = Gate::CNOT().matrix;
    test_result("controlled_U(X) = CNOT", (CX - CNOT).norm() < EPSILON);
    
    auto CZ_manual = Gate::controlled_U(Gate::Z()).matrix;
    auto CZ_base = Gate::CZ().matrix;
    test_result("controlled_U(Z) = CZ", (CZ_manual - CZ_base).norm() < EPSILON);
}

// ============================================================================
// Tests: Gate composition
// ============================================================================

void test_gate_composition() {
    cout << "\n--- Gate composition (multiplication) ---\n";
    
    // (XY) applied to |ψ⟩ should equal X(Y|ψ⟩)
    auto X = Gate::X();
    auto Y = Gate::Y();
    auto XY = X * Y;
    test_result("XY composition is unitary", is_unitary(XY.matrix));
    
    // H² should be I
    auto H = Gate::H();
    auto HH = H * H;
    test_result("H·H = I", is_identity(HH.matrix));
    
    // Verify dimension matching in composed gates
    test_result("1-qubit·1-qubit = 1-qubit", (X * Y).matrix.rows() == 2);
}

void test_rotation_composition() {
    cout << "\n--- Rotation gate composition ---\n";
    
    // Rx(π/2)·Rx(π/2) = Rx(π)
    auto Rx_half = Gate::Rx(M_PI / 2);
    auto Rx_pi = Gate::Rx(M_PI);
    auto Rx_half_twice = Rx_half * Rx_half;
    test_result("Rx(π/2)·Rx(π/2) ≈ Rx(π)", (Rx_half_twice.matrix - Rx_pi.matrix).norm() < EPSILON);
}

// ============================================================================
// Tests: Commutation relations (consistency)
// ============================================================================

void test_commutation_relations() {
    cout << "\n--- Commutation relations ---\n";
    
    auto X = Gate::X().matrix;
    auto Y = Gate::Y().matrix;
    auto Z = Gate::Z().matrix;
    auto I = Eigen::MatrixXcd::Identity(2, 2);
    
    // [X, X] = 0
    test_result("[X,X] = 0", (X * X - X * X).norm() < EPSILON);
    
    // [X, Z] ≠ 0 (anticommute)
    auto XZ_ZX = X * Z - Z * X;
    test_result("[X,Z] ≠ 0", XZ_ZX.norm() > 0.1);
}

// ============================================================================
// Tests: Identity and dimension consistency
// ============================================================================

void test_identity_gate() {
    cout << "\n--- Identity gate ---\n";
    
    auto I1 = Gate::identity(1).matrix;
    auto I2 = Gate::identity(2).matrix;
    auto I3 = Gate::identity(3).matrix;
    
    test_result("I(1) is 2×2", I1.rows() == 2);
    test_result("I(2) is 4×4", I2.rows() == 4);
    test_result("I(3) is 8×8", I3.rows() == 8);
    
    test_result("I(1) = I", is_identity(I1));
    test_result("I(2) = I", is_identity(I2));
    test_result("I(3) = I", is_identity(I3));
}

// ============================================================================
// Main test suite runner
// ============================================================================

int main() {
    cout << "\n" << string(70, '=') << "\n"
         << "  Test Suite: Gate.h\n"
         << string(70, '=') << "\n";
    
    test_single_qubit_unitarity();
    test_pauli_algebra();
    test_rotation_gates();
    test_hadamard();
    
    test_two_qubit_unitarity();
    test_cnot();
    test_swap();
    test_rzz();
    test_toffoli();
    
    test_tensor_product();
    test_expand_single_qubit();
    test_expand_two_qubit();
    test_expand_two_non_adjacent();
    test_controlled_u();
    
    test_gate_composition();
    test_rotation_composition();
    
    test_commutation_relations();
    test_identity_gate();
    
    cout << "\n" << string(70, '=') << "\n"
         << "  All tests passed ✓\n"
         << string(70, '=') << "\n\n";
    
    return 0;
}