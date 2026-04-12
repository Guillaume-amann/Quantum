#pragma once
#include <vector>
#include <complex>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
using Complex = complex<double>;

class Gate {
public:
    MatrixXcd matrix;
    int size;  // number of qubits the gate acts on

    explicit Gate(int num_qubits) : matrix(MatrixXcd::Zero(1 << num_qubits, 1 << num_qubits)), size(num_qubits) {}

    Gate(const vector<vector<Complex>>& mat) : matrix(mat.size(), mat.size()), size(0) {
        int d = static_cast<int>(mat.size());
        if ((d & (d - 1)) != 0) throw invalid_argument("Matrix size must be a power of 2");
        size = static_cast<int>(log2(d));
        for (int i = 0; i < d; ++i)
            for (int j = 0; j < d; ++j)
                matrix(i, j) = mat[i][j];
    }

    // =========================================================================
    // Single-qubit gates

    static Gate identity(int n) {
        Gate g(n);
        g.matrix = MatrixXcd::Identity(1 << n, 1 << n);
        return g;
    }

    static Gate X() {
        Gate g(1);
        g.matrix << 0, 1,
                    1, 0;
        return g;
    }

    static Gate Y() {
        Gate g(1);
        g.matrix << Complex(0, 0), Complex(0,-1),
                    Complex(0, 1), Complex(0, 0);
        return g;
    }

    static Gate Z() {
        Gate g(1);
        g.matrix << 1,  0,
                    0, -1;
        return g;
    }

    static Gate H() {
        Gate g(1);
        g.matrix << 1.0 / sqrt(2.0),  1.0 / sqrt(2.0),
                    1.0 / sqrt(2.0), -1.0 / sqrt(2.0);
        return g;
    }

    static Gate Rx(double theta) {
        Gate g(1);
        g.matrix << cos(theta / 2),  Complex(0, -sin(theta / 2)),
                    Complex(0, sin(theta / 2)),  cos(theta / 2);
        return g;
    }

    static Gate Ry(double theta) {
        Gate g(1);
        g.matrix << cos(theta / 2), -sin(theta / 2),
                    sin(theta / 2),  cos(theta / 2);
        return g;
    }

    static Gate Rz(double theta) {
        Gate g(1);
        g.matrix << exp(Complex(0, -theta / 2)), 0,
                    0, exp(Complex(0,  theta / 2));
        return g;
    }

    // =========================================================================
    // Two-qubit gates (4×4, ctrl-first basis: |00⟩ |01⟩ |10⟩ |11⟩)

    static Gate CNOT() {
        Gate g(2);
        g.matrix << 1,0,0,0,
                    0,1,0,0,
                    0,0,0,1,
                    0,0,1,0;
        return g;
    }

    static Gate CZ() {
        Gate g(2);
        g.matrix << 1, 0, 0,  0,
                    0, 1, 0,  0,
                    0, 0, 1,  0,
                    0, 0, 0, -1;
        return g;
    }

    static Gate SWAP() {
        Gate g(2);
        g.matrix << 1,0,0,0,
                    0,0,1,0,
                    0,1,0,0,
                    0,0,0,1;
        return g;
    }

    // Rzz(θ) = exp(-i θ/2  Z⊗Z)
    static Gate Rzz(double theta) {
        Gate g(2);
        g.matrix << exp(Complex(0, -theta / 2.0)), 0, 0,  0,
                    0, exp(Complex(0,  theta / 2.0)), 0,  0,
                    0, 0, exp(Complex(0,  theta / 2.0)),  0,
                    0, 0, 0, exp(Complex(0, -theta / 2.0));
        return g;
    }

    // Controlled-U: applies single-qubit U to target when control = |1⟩
    static Gate controlled_U(const Gate& U) {
        if (U.size != 1) throw invalid_argument("controlled_U(): U must be a single-qubit gate");
        Gate g(2);
        g.matrix(0,0) = 1;
        g.matrix(1,1) = 1;
        g.matrix(2,2) = U.matrix(0,0);  g.matrix(2,3) = U.matrix(0,1);
        g.matrix(3,2) = U.matrix(1,0);  g.matrix(3,3) = U.matrix(1,1);
        return g;
    }

    // =========================================================================
    // Three-qubit gate
    
    // Toffoli (CCX): flips target when both controls are |1⟩
    static Gate Toffoli() {
        Gate g(3);
        g.matrix = MatrixXcd::Identity(8, 8);
        g.matrix(6,6) = 0;  g.matrix(6,7) = 1;
        g.matrix(7,7) = 0;  g.matrix(7,6) = 1;
        return g;
    }

    // =========================================================================


    // Tensor product:  A ⊗ B
    static Gate tensor(const Gate& A, const Gate& B) {
        int dimA = 1 << A.size, dimB = 1 << B.size;
        Gate result(A.size + B.size);
        for (int i = 0; i < dimA; ++i)
            for (int j = 0; j < dimA; ++j)
                for (int k = 0; k < dimB; ++k)
                    for (int l = 0; l < dimB; ++l)
                        result.matrix(i * dimB + k, j * dimB + l) =
                            A.matrix(i, j) * B.matrix(k, l);
        return result;
    }

    Gate expand(int total_qubits, int target) const {
        if (size != 1) throw invalid_argument("expand(): only single-qubit gates can be expanded");
        Gate result = identity(0);
        for (int i = 0; i < total_qubits; ++i)
            result = tensor(result, (i == target) ? *this : identity(1));
        return result;
    }

    static Gate rotationX(int target, int total_qubits, double theta) {
        return Rx(theta).expand(total_qubits, target);
    }


    static Gate expand_two(const Gate& G, int ctrl, int tgt, int total_qubits) {
        if (G.size != 2) throw invalid_argument("expand_two(): G must be a 2-qubit gate");
        if (ctrl == tgt || ctrl < 0 || tgt < 0 || ctrl >= total_qubits || tgt >= total_qubits) throw invalid_argument("expand_two(): invalid ctrl/tgt indices");

        int n   = total_qubits;
        int dim = 1 << n;
        // Bit positions of ctrl and tgt within an n-bit integer (MSB = qubit 0).
        int ctrl_shift = n - 1 - ctrl;
        int tgt_shift  = n - 1 - tgt;
        Gate result(n);

        for (int row = 0; row < dim; ++row) {
            // Extract the ctrl and tgt bits from the row basis state and
            // pack them into the 2-bit index expected by G (ctrl is MSB).
            int row_ctrl = (row >> ctrl_shift) & 1;
            int row_tgt  = (row >> tgt_shift)  & 1;
            int row_pair = row_ctrl * 2 + row_tgt;

            for (int col = 0; col < dim; ++col) {
                // Mask covering every qubit except ctrl and tgt.
                // If the spectator bits differ between row and col the gate
                // acts as identity on them, so the matrix element is zero.
                int other_mask = ~((1 << ctrl_shift) | (1 << tgt_shift));
                if ((row & other_mask) != (col & other_mask)) continue;

                // Same 2-bit packing for the column basis state.
                int col_ctrl = (col >> ctrl_shift) & 1;
                int col_tgt  = (col >> tgt_shift)  & 1;
                int col_pair = col_ctrl * 2 + col_tgt;

                // Copy the corresponding element of the 2-qubit gate G.
                result.matrix(row, col) = G.matrix(row_pair, col_pair);
            }
        }
        return result;
    }
    
    static Gate cnot (int ctrl, int tgt, int n) { return expand_two(CNOT(), ctrl, tgt, n); }
    static Gate cz (int ctrl, int tgt, int n) { return expand_two(CZ(), ctrl, tgt, n); }
    static Gate swap (int q0,   int q1,  int n) { return expand_two(SWAP(), q0, q1, n); }
    static Gate rzz (double t, int q0,  int q1, int n) { return expand_two(Rzz(t), q0, q1, n); }
    static Gate controlled(const Gate& U, int ctrl, int tgt, int n) { return expand_two(controlled_U(U), ctrl, tgt, n); }
};

// =============================================================================
// Gate multiplication (sequential composition):  C = A·B  →  C|ψ⟩ = A(B|ψ⟩)
// =============================================================================

inline Gate operator*(const Gate& A, const Gate& B) {
    if (A.matrix.rows() != B.matrix.rows())
        throw invalid_argument("Gate size mismatch for multiplication");
    Gate result(A.size);
    result.matrix.noalias() = A.matrix * B.matrix;
    return result;
}