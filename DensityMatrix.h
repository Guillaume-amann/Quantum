#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <string>
#include <Eigen/Dense>
#include "Gate.h"

using namespace std;
using namespace Eigen;
using Complex = complex<double>;

// =============================================================================
//   Pure state definition          ->  ρ  = |ψ⟩⟨ψ|,                Tr(ρ^2) = 1 (purity condition)
//   Unitary pure state evolution   ->  ρ' = U ρ U+
//   Mixed state definition         ->  ρ  = Σ_k p_k |ψ_k⟩⟨ψ_k|,    purity < 1
//   Kraus Noisy evolution          ->  ρ' = Σ_k K_k ρ K_k+
//
//   - Gate noise via Kraus operators:  ρ' = Σ_k K_k ρ K_k+
//   - Partial trace produces mixed states even from a pure input
//   - Von Neumann entropy is only meaningful on ρ, not on state vectors
//
// Basis convention: big-endian (same as Qbit and Gate)
//   |q0 q1 ... q_{n-1}⟩  ->  index = q0*2^(n-1) + q1*2^(n-2) + ... + q_{n-1}
// =============================================================================

class DensityMatrix {
private:
    int n;    // number of qubits
    int dim;  // 2^n
    MatrixXcd rho;

    // Enforce Hermiticity and Tr(ρ)=1 (norm. condition) after any in-place operation
    void renormalise() {
        rho = (rho + rho.adjoint()) * 0.5;
        double tr = rho.trace().real();
        if (tr < 1e-15) {
            rho.setZero();
            rho(0, 0) = 1.0;
        } else {
            rho /= tr;
        }
    }

    // Private constructor used by partial_trace (cf. L121) (skips size checks)
    DensityMatrix(int num_qubits, MatrixXcd rho_mat) : n(num_qubits), dim(1 << num_qubits), rho(std::move(rho_mat)){
        renormalise();
    }

public:
    explicit DensityMatrix(int num_qubits) : n(num_qubits), dim(1 << num_qubits), rho(MatrixXcd::Zero(1 << num_qubits, 1 << num_qubits)) {
        rho(0, 0) = 1.0;
    }

    // From a state vector |ψ⟩  ->  ρ = |ψ⟩⟨ψ|
    DensityMatrix(int num_qubits, const vector<Complex>& psi) : n(num_qubits), dim(1 << num_qubits), rho(MatrixXcd::Zero(1 << num_qubits, 1 << num_qubits)){
        if ((int)psi.size() != dim) throw invalid_argument("DensityMatrix: state vector size mismatch");
        Map<const VectorXcd> v(psi.data(), dim);
        VectorXcd vn = v.normalized();
        rho = vn * vn.adjoint();
    }

    // From a raw ρ matrix (e.g. returned by partial_trace)
    DensityMatrix(int num_qubits, const vector<vector<Complex>>& rho_in) : n(num_qubits), dim(1 << num_qubits), rho(dim, dim) {
        if ((int)rho_in.size() != dim || (int)rho_in[0].size() != dim) throw invalid_argument("DensityMatrix: matrix size mismatch");
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                rho(i, j) = rho_in[i][j];
        renormalise();
    }

    // Unitary evolution:  ρ' = U ρ U+
    void apply(const Gate& U) {
        if (U.matrix.rows() != dim) throw invalid_argument("apply(): gate dimension does not match register size");
        rho = U.matrix * rho * U.matrix.adjoint();
        renormalise();
    }

    // Kraus channel Noisy evolution:  ρ' = Σ_k K_k ρ K_k+
    //      pass the output of NoiseModel::kraus_*() directly
    void apply_kraus(const vector<Gate>& kraus_ops) {
        // Accumulate Σ_k K_k† K_k — must equal I for the channel to be trace-preserving
        // (If < I the channel is trace-decreasing, i.e. it can "lose" probability)
        MatrixXcd completeness = MatrixXcd::Zero(dim, dim);
        MatrixXcd result = MatrixXcd::Zero(dim, dim);
        for (const Gate& K : kraus_ops) {
            if (K.matrix.rows() != dim) throw invalid_argument("apply_kraus(): Kraus operator size mismatch");
            completeness.noalias() += K.matrix.adjoint() * K.matrix;
            result.noalias() += K.matrix * rho * K.matrix.adjoint();
        }
        // Frobenius norm of the deviation from I; anything above numerical noise means
        // the operator set is incomplete and the evolution would not preserve Tr(ρ) = 1.
        if ((completeness - MatrixXcd::Identity(dim, dim)).norm() > 1e-9) throw invalid_argument("apply_kraus(): Kraus operators do not satisfy completeness relation Σ K†K = I");
        rho = result;
        renormalise();
    }

    // Partial trace over qubit k, producing a reduced density matrix for the remaining n-1 qubits
    DensityMatrix partial_trace(int qubit) const {
        if (qubit < 0 || qubit >= n) throw invalid_argument("partial_trace(): qubit index out of range");

        // Output system has one fewer qubit
        int n_out   = n - 1;
        int dim_out = 1 << n_out;

        // Bit position of the traced-out qubit in big-endian indexing   (  full index bits:  [ q0 | q1 | q2 ]  )
        // and a mask for the bits that sit below it                     (  bit positions:       2    1    0    )
        int k_shift = n - 1 - qubit;
        int lo_mask = (1 << k_shift) - 1;

        MatrixXcd rho_out = MatrixXcd::Zero(dim_out, dim_out);

        // reduced matrix of the system after tracing out qubit k 
        // |q0 q1 q2⟩ becomes |q0 q2⟩ after tracing out q1
        for (int i_out = 0; i_out < dim_out; ++i_out)
            for (int j_out = 0; j_out < dim_out; ++j_out)
                for (int b = 0; b < 2; ++b) {
                    int i_full = ((i_out >> k_shift) << (k_shift + 1)) | (b << k_shift) | (i_out & lo_mask);
                    int j_full = ((j_out >> k_shift) << (k_shift + 1)) | (b << k_shift) | (j_out & lo_mask);
                    rho_out(i_out, j_out) += rho(i_full, j_full);
                }

        return DensityMatrix(n_out, std::move(rho_out));
    }

    // Measure a single qubit k  ->  returns outcome (0 or 1) and collapses ρ
    int partial_measurement(int qubit) {
        if (qubit < 0 || qubit >= n) throw invalid_argument("partial_measurement(): qubit index out of range");

        int k_shift = n - 1 - qubit;

        double p0 = 0.0;
        for (int i = 0; i < dim; ++i)
            if (((i >> k_shift) & 1) == 0)
                p0 += rho(i, i).real();

        static thread_local mt19937 gen(random_device{}());
        uniform_real_distribution<double> dist(0.0, 1.0);
        int outcome = (dist(gen) < p0) ? 0 : 1;
        double prob = (outcome == 0)   ? p0 : 1.0 - p0;

        if (prob < 1e-15) throw runtime_error("partial_measurement(): outcome probability ~0");

        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                if (((i >> k_shift) & 1) != outcome || ((j >> k_shift) & 1) != outcome)
                    rho(i, j) = 0;

        renormalise();
        return outcome;
    }

    // Measure (and collapse) the whole system
    string measure() {
        vector<double> cumulative(dim);
        cumulative[0] = rho(0, 0).real();
        for (int i = 1; i < dim; ++i)
            cumulative[i] = cumulative[i-1] + rho(i, i).real();

        static thread_local mt19937 gen(random_device{}());
        uniform_real_distribution<double> dist(0.0, 1.0);
        double r = dist(gen);
        int idx = 0;
        while (idx < dim - 1 && r > cumulative[idx]) ++idx;

        rho.setZero();
        rho(idx, idx) = 1.0;

        string s = "|";
        for (int i = n - 1; i >= 0; --i)
            s += ((idx >> i) & 1) ? '1' : '0';
        return s + "⟩";
    }

    // Displays amplitudes without collapsing — a "cheat view" impossible on real hardware.
    string print_state() const {
        vector<double> cumulative(dim);
        cumulative[0] = rho(0, 0).real();
        for (int i = 1; i < dim; ++i)
            cumulative[i] = cumulative[i-1] + rho(i, i).real();

        static thread_local mt19937 gen(random_device{}());
        uniform_real_distribution<double> dist(0.0, 1.0);
        double r = dist(gen);
        int idx = 0;
        while (idx < dim - 1 && r > cumulative[idx]) ++idx;

        string s = "|";
        for (int i = n - 1; i >= 0; --i)
            s += ((idx >> i) & 1) ? '1' : '0';
        return s + "⟩";
    }

    // Expectation value of a Hermitian observable O:  <O> = Tr(O ρ)
    double expectation(const Gate& O) const {
        if (O.matrix.rows() != dim) throw invalid_argument("expectation(): observable dimension mismatch");
        return (O.matrix * rho).trace().real();
    }

    // Purity:  Tr(ρ^2)  -- 1 for pure, 1/dim for maximally mixed
    double purity() const {
        return (rho * rho).trace().real();
    }

    // Von Neumann entropy:  S(ρ) = -Tr(ρ log ρ) = -Σ_k λ_k log(λ_k)
    double entropy() const {
        SelfAdjointEigenSolver<MatrixXcd> solver(rho, EigenvaluesOnly);
        double S = 0.0;
        for (int k = 0; k < solver.eigenvalues().size(); ++k) {
            double lam = solver.eigenvalues()(k);
            if (lam > 1e-15)
                S -= lam * log(lam);
        }
        return S;
    }

    // Fidelity with a pure state |ψ⟩:  F = ⟨ψ|ρ|ψ⟩
    double fidelity(const vector<Complex>& psi) const {
        if ((int)psi.size() != dim) throw invalid_argument("fidelity(): state vector size mismatch");
        Map<const VectorXcd> v(psi.data(), dim);
        return v.dot(rho * v).real();
    }

    // Probability vector: diagonal of ρ
    vector<double> probabilities() const {
        vector<double> p(dim);
        for (int i = 0; i < dim; ++i) p[i] = rho(i, i).real();
        return p;
    }

    // Prints measurement probabilities (diagonal of ρ) and purity
    // Physically meaningful — equivalent to what repeated measurements would reconstruct statistically
    void print_probabilities() const {
        cout << "ρ (non-zero diagonals):\n";
        for (int i = 0; i < dim; ++i) {
            if (rho(i, i).real() > 1e-10) {
                cout << "  p(|";
                for (int b = n - 1; b >= 0; --b)
                    cout << ((i >> b) & 1);
                cout << "⟩) = " << rho(i, i).real() << "\n";
            }
        }
        cout << "  purity = " << purity() << "\n";
    }

    // Prints the full ρ matrix including off-diagonal coherences
    // Cheat view — coherences are invisible to any real measurement
    void print_rho() const {
        cout << "rho (" << dim << "x" << dim << "):\n";
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                if (abs(rho(i, j)) > 1e-10) {
                    cout << "  [" << i << "][" << j << "] = " << rho(i, j).real();
                    if (rho(i, j).imag() >  1e-10) cout << "+" << rho(i, j).imag() << "i";
                    else if (rho(i, j).imag() < -1e-10) cout << rho(i, j).imag() << "i";
                    cout << "\n";
                }
            }
        }
    }

    int get_num_qubits() const { return n; }
    int get_dimension() const { return dim; }
    const MatrixXcd& get_rho() const { return rho; }
};