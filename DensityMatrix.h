#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <stdexcept>
#include <string>
#include <Eigen/Dense>
#include "Gate.h"

using namespace std;
using Complex = complex<double>;

// =============================================================================
// DensityMatrix — mixed-state description of an n-qubit register.
//
// State: rho, a 2^n x 2^n complex Hermitian matrix with Tr(rho) = 1.
//
//   Pure state |psi>  ->  rho = |psi><psi|,  purity Tr(rho^2) = 1
//   Mixed state       ->  rho = sum_k p_k |psi_k><psi_k|,  purity < 1
//
// Why over Qbit:
//   - Gate noise via Kraus operators:  rho' = sum_k K_k rho K_k+
//   - Partial trace produces mixed states even from a pure input
//   - Von Neumann entropy is only meaningful on rho
//
// Basis convention: big-endian (same as Qbit and Gate).
//   |q0 q1 ... q_{n-1}>  ->  index = q0*2^(n-1) + q1*2^(n-2) + ... + q_{n-1}
// =============================================================================

class DensityMatrix {
private:
    int n;    // number of qubits
    int dim;  // 2^n
    Eigen::MatrixXcd rho;

    // Enforce Hermiticity and Tr(rho)=1 after any in-place operation
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

    // Private constructor used by partial_trace (skips size checks)
    DensityMatrix(int num_qubits, Eigen::MatrixXcd rho_mat)
        : n(num_qubits), dim(1 << num_qubits), rho(std::move(rho_mat))
    {
        renormalise();
    }

public:
    // =========================================================================
    // Constructors
    // =========================================================================

    // Default: |0...0><0...0|
    explicit DensityMatrix(int num_qubits)
        : n(num_qubits), dim(1 << num_qubits),
          rho(Eigen::MatrixXcd::Zero(1 << num_qubits, 1 << num_qubits))
    {
        rho(0, 0) = 1.0;
    }

    // From a state vector |psi>  ->  rho = |psi><psi|
    DensityMatrix(int num_qubits, const vector<Complex>& psi)
        : n(num_qubits), dim(1 << num_qubits),
          rho(Eigen::MatrixXcd::Zero(1 << num_qubits, 1 << num_qubits))
    {
        if ((int)psi.size() != dim)
            throw invalid_argument("DensityMatrix: state vector size mismatch");
        Eigen::Map<const Eigen::VectorXcd> v(psi.data(), dim);
        Eigen::VectorXcd vn = v.normalized();
        rho = vn * vn.adjoint();
    }

    // From a raw rho matrix (e.g. returned by partial_trace)
    DensityMatrix(int num_qubits, const vector<vector<Complex>>& rho_in)
        : n(num_qubits), dim(1 << num_qubits), rho(dim, dim)
    {
        if ((int)rho_in.size() != dim || (int)rho_in[0].size() != dim)
            throw invalid_argument("DensityMatrix: matrix size mismatch");
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                rho(i, j) = rho_in[i][j];
        renormalise();
    }

    // =========================================================================
    // Unitary evolution:  rho' = U rho U+
    // =========================================================================
    void apply(const Gate& U) {
        if (U.matrix.rows() != dim)
            throw invalid_argument("apply(): gate dimension does not match register size");
        rho = U.matrix * rho * U.matrix.adjoint();
        renormalise();
    }

    // =========================================================================
    // Kraus channel:  rho' = sum_k K_k rho K_k+
    //
    // Kraus operators must satisfy sum_k K_k+ K_k = I (not checked).
    // Pass the output of NoiseModel::kraus_*() directly.
    // =========================================================================
    void apply_kraus(const vector<Gate>& kraus_ops) {
        Eigen::MatrixXcd result = Eigen::MatrixXcd::Zero(dim, dim);
        for (const Gate& K : kraus_ops) {
            if (K.matrix.rows() != dim)
                throw invalid_argument("apply_kraus(): Kraus operator size mismatch");
            result.noalias() += K.matrix * rho * K.matrix.adjoint();
        }
        rho = result;
        renormalise();
    }

    // =========================================================================
    // Partial trace over qubit k  ->  new DensityMatrix on (n-1) qubits
    //
    // rho_reduced[i][j] = sum_{b=0,1}  rho[insert(i,k,b)][insert(j,k,b)]
    // =========================================================================
    DensityMatrix partial_trace(int qubit) const {
        if (qubit < 0 || qubit >= n)
            throw invalid_argument("partial_trace(): qubit index out of range");

        int n_out   = n - 1;
        int dim_out = 1 << n_out;
        int k_shift = n - 1 - qubit;
        int lo_mask = (1 << k_shift) - 1;
        Eigen::MatrixXcd rho_out = Eigen::MatrixXcd::Zero(dim_out, dim_out);

        for (int i_out = 0; i_out < dim_out; ++i_out)
            for (int j_out = 0; j_out < dim_out; ++j_out)
                for (int b = 0; b < 2; ++b) {
                    int i_full = ((i_out >> k_shift) << (k_shift + 1)) | (b << k_shift) | (i_out & lo_mask);
                    int j_full = ((j_out >> k_shift) << (k_shift + 1)) | (b << k_shift) | (j_out & lo_mask);
                    rho_out(i_out, j_out) += rho(i_full, j_full);
                }

        return DensityMatrix(n_out, std::move(rho_out));
    }

    // =========================================================================
    // Measure a single qubit k  ->  returns outcome (0 or 1) and collapses rho
    // =========================================================================
    int measure_qubit(int qubit) {
        if (qubit < 0 || qubit >= n)
            throw invalid_argument("measure_qubit(): qubit index out of range");

        int k_shift = n - 1 - qubit;

        double p0 = 0.0;
        for (int i = 0; i < dim; ++i)
            if (((i >> k_shift) & 1) == 0)
                p0 += rho(i, i).real();

        static thread_local mt19937 gen(random_device{}());
        uniform_real_distribution<double> dist(0.0, 1.0);
        int outcome = (dist(gen) < p0) ? 0 : 1;
        double prob = (outcome == 0)   ? p0 : 1.0 - p0;

        if (prob < 1e-15)
            throw runtime_error("measure_qubit(): outcome probability ~0");

        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                if (((i >> k_shift) & 1) != outcome || ((j >> k_shift) & 1) != outcome)
                    rho(i, j) = 0;

        renormalise();
        return outcome;
    }

    // =========================================================================
    // Full-register measurement  ->  collapses to a basis state, returns bitstring
    // =========================================================================
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
        return s + ">";
    }

    // =========================================================================
    // Non-destructive sample (does not collapse rho)
    // =========================================================================
    string sample() const {
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
        return s + ">";
    }

    // =========================================================================
    // Expectation value of a Hermitian observable O:  <O> = Tr(O rho)
    // =========================================================================
    double expectation(const Gate& O) const {
        if (O.matrix.rows() != dim)
            throw invalid_argument("expectation(): observable dimension mismatch");
        return (O.matrix * rho).trace().real();
    }

    // =========================================================================
    // Purity:  Tr(rho^2)  -- 1 for pure, 1/dim for maximally mixed
    // =========================================================================
    double purity() const {
        return (rho * rho).trace().real();
    }

    // =========================================================================
    // Von Neumann entropy:  S(rho) = -Tr(rho log rho) = -sum_k lambda_k log(lambda_k)
    // =========================================================================
    double entropy() const {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(rho, Eigen::EigenvaluesOnly);
        double S = 0.0;
        for (int k = 0; k < solver.eigenvalues().size(); ++k) {
            double lam = solver.eigenvalues()(k);
            if (lam > 1e-15)
                S -= lam * log(lam);
        }
        return S;
    }

    // =========================================================================
    // Fidelity with a pure state |psi>:  F = <psi|rho|psi>
    // =========================================================================
    double fidelity(const vector<Complex>& psi) const {
        if ((int)psi.size() != dim)
            throw invalid_argument("fidelity(): state vector size mismatch");
        Eigen::Map<const Eigen::VectorXcd> v(psi.data(), dim);
        return v.dot(rho * v).real();
    }

    // =========================================================================
    // Probability vector: diagonal of rho
    // =========================================================================
    vector<double> probabilities() const {
        vector<double> p(dim);
        for (int i = 0; i < dim; ++i) p[i] = rho(i, i).real();
        return p;
    }

    // =========================================================================
    // Accessors
    // =========================================================================

    int num_qubits()                  const { return n; }
    int dimension()                   const { return dim; }
    const Eigen::MatrixXcd& get_rho() const { return rho; }

    void print_state() const {
        cout << "rho (non-zero diagonals):\n";
        for (int i = 0; i < dim; ++i) {
            if (rho(i, i).real() > 1e-10) {
                cout << "  p(|";
                for (int b = n - 1; b >= 0; --b)
                    cout << ((i >> b) & 1);
                cout << ">) = " << rho(i, i).real() << "\n";
            }
        }
        cout << "  purity = " << purity() << "\n";
    }

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
};