#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <Eigen/Dense>
#include "Gate.h"
#include "DensityMatrix.h"

using namespace std;
using namespace Eigen;
using Complex = complex<double>;

// =============================================================================
//  Quantum Kernel — QSVM feature map and kernel computation
//
//  Purpose:
//    Encode classical data x ∈ ℝ^d into quantum feature map |ψ(x)⟩,
//    then compute quantum kernel K_Q(xᵢ, xⱼ) = |⟨ψ(xᵢ)|ψ(xⱼ)⟩|²
//
//  Architecture:
//    1. Feature map: RealAmplitudes ansatz (RY rotations + CNOT entanglement)
//    2. Kernel computation via fidelity (DensityMatrix backend)
//    3. Kernel matrix assembly for SVM dual problem
//
//  Usage:
//    QKernel qk(n_qubits=8, n_features=11, depth=2);
//    auto K_Q = qk.compute_kernel_matrix(X_train);
//    // Pass K_Q to classical SVM solver
//
//  Mathematical background:
//    φ_Q(x) = U(θ(x))|0...0⟩  where U = parametrized quantum circuit
//    K_Q(xᵢ, xⱼ) = |⟨φ_Q(xᵢ)|φ_Q(xⱼ)⟩|² = fidelity squared
//
// =============================================================================

class QKernel {
private:
    int n_qubits;           // number of qubits
    int n_features;         // input dimension (data features)
    int depth;              // circuit depth (number of layers)
    
    // Circuit parameters: weights for RY rotations in each layer
    // Structure: params[layer][qubit][feature] = weight for RY(θ·x)
    // Total: depth × n_qubits × (n_features + 1) parameters
    // (+1 for bias term in each RY angle)
    vector<vector<vector<double>>> params;
    
    constexpr static double PI = 3.14159265358979323846;

public:
    // =========================================================================
    // Constructor: initialize with random parameters
    // =========================================================================
    QKernel(int num_qubits = 8, int num_features = 11, int circuit_depth = 2)
        : n_qubits(num_qubits), n_features(num_features), depth(circuit_depth)
    {
        if (n_qubits < 2 || n_qubits > 15) 
            throw invalid_argument("QKernel: n_qubits must be in [2, 15]");
        if (n_features < 1 || n_features > 100)
            throw invalid_argument("QKernel: n_features must be in [1, 100]");
        if (circuit_depth < 1 || circuit_depth > 5)
            throw invalid_argument("QKernel: depth must be in [1, 5]");
        
        // Initialize parameters (small random values)
        init_random_params();
    }

    // =========================================================================
    // Initialize circuit parameters randomly
    // =========================================================================
    void init_random_params() {
        params.clear();
        static thread_local mt19937 gen(random_device{}());
        uniform_real_distribution<double> dist(-0.1, 0.1);
        
        for (int l = 0; l < depth; ++l) {
            vector<vector<double>> layer;
            for (int q = 0; q < n_qubits; ++q) {
                vector<double> qubit_params;
                for (int f = 0; f <= n_features; ++f) {  // +1 for bias
                    qubit_params.push_back(dist(gen));
                }
                layer.push_back(qubit_params);
            }
            params.push_back(layer);
        }
    }

    // =========================================================================
    // Prepare quantum feature map: |ψ(x)⟩ = U(θ(x))|0...0⟩
    // =========================================================================
    DensityMatrix prepare_feature_map(const vector<double>& x) const {
        if ((int)x.size() != n_features)
            throw invalid_argument("prepare_feature_map: input dimension mismatch");
        
        // Start with |0...0⟩ in density matrix form
        DensityMatrix rho(n_qubits);  // |0⟩⟨0| on each qubit
        
        // Apply RealAmplitudes circuit: L layers of (RY rotations + entanglement)
        for (int layer = 0; layer < depth; ++layer) {
            // Layer: RY rotations on each qubit, parameterized by x
            for (int q = 0; q < n_qubits; ++q) {
                // Compute angle: θ = Σⱼ αⱼ xⱼ + β (linear in x + bias)
                double angle = 0.0;
                for (int f = 0; f < n_features; ++f) {
                    angle += params[layer][q][f] * x[f];
                }
                angle += params[layer][q][n_features];  // bias term
                
                // Apply RY(angle) on qubit q
                rho.apply(Gate::Ry(angle).expand(n_qubits, q));
            }
            
            // Entanglement: chain of CNOT gates (linear connectivity)
            // CNOT(0,1), CNOT(2,3), ... for even n_qubits
            // This creates a simple entanglement pattern
            for (int q = 0; q < n_qubits - 1; ++q) {
                rho.apply(Gate::cnot(q, q+1, n_qubits));
            }
        }
        
        return rho;
    }

    // =========================================================================
    // Compute quantum kernel: K_Q(xᵢ, xⱼ) = |⟨ψ(xᵢ)|ψ(xⱼ)⟩|²
    // =========================================================================
    double kernel(const vector<double>& x_i, const vector<double>& x_j) const {
        if ((int)x_i.size() != n_features || (int)x_j.size() != n_features)
            throw invalid_argument("kernel: input dimension mismatch");
        
        // Prepare feature maps
        DensityMatrix psi_i = prepare_feature_map(x_i);  // |ψ(xᵢ)⟩⟨ψ(xᵢ)|
        
        // Extract state vector from density matrix (assuming pure state)
        // For a pure state |ψ⟩, ρ = |ψ⟩⟨ψ|, so ρ diagonal entries are |⟨k|ψ⟩|²
        // We reconstruct ψ from ρ (assuming it's pure)
        const auto& rho_i_matrix = psi_i.get_rho();
        vector<Complex> psi_i_vec(1 << n_qubits);
        for (int idx = 0; idx < (1 << n_qubits); ++idx) {
            // Extract amplitude from density matrix
            // For pure |ψ⟩, ⟨idx|ψ⟩ can be extracted via careful measurement
            // Simplified: take diagonal as probability and phase from off-diagonals
            double prob = rho_i_matrix(idx, idx).real();
            psi_i_vec[idx] = Complex(sqrt(max(0.0, prob)), 0.0);
        }
        
        // Compute fidelity F = |⟨ψⱼ|ψᵢ⟩|² using DensityMatrix.fidelity()
        // which expects a state vector
        double fidelity_val = psi_i.fidelity(psi_i_vec);
        
        // Actually, fidelity expects pure state. Better approach:
        // Compute overlap directly using DensityMatrix of xⱼ
        DensityMatrix psi_j = prepare_feature_map(x_j);
        
        // K_Q = Tr(ρᵢ · ρⱼ) for pure states = |⟨ψᵢ|ψⱼ⟩|²
        // But simpler: use fidelity of the two density matrices
        auto rho_j_matrix = psi_j.get_rho();
        
        // Compute fidelity: F(ρ,σ) = Tr(√(√ρ σ √ρ))²
        // For pure states: F = |⟨ψ|φ⟩|²
        // Direct computation: Tr(ρᵢ · ρⱼ)
        Complex overlap_trace = (rho_i_matrix * rho_j_matrix).trace();
        double K_Q = abs(overlap_trace) * abs(overlap_trace);  // squared
        
        return K_Q;
    }

    // =========================================================================
    // Compute full kernel matrix K_Q[i][j] for training set
    // =========================================================================
    MatrixXd compute_kernel_matrix(const vector<vector<double>>& X) const {
        int n = X.size();
        MatrixXd K_Q(n, n);
        
        for (int i = 0; i < n; ++i) {
            for (int j = i; j < n; ++j) {
                double k_ij = kernel(X[i], X[j]);
                K_Q(i, j) = k_ij;
                K_Q(j, i) = k_ij;  // symmetric
            }
        }
        
        return K_Q;
    }

    // =========================================================================
    // Compute kernel row: K_Q[i, :] for a single test point
    // =========================================================================
    VectorXd kernel_row(const vector<double>& x_test, 
                       const vector<vector<double>>& X_train) const {
        int n = X_train.size();
        VectorXd k_row(n);
        
        for (int i = 0; i < n; ++i) {
            k_row(i) = kernel(X_train[i], x_test);
        }
        
        return k_row;
    }

    // =========================================================================
    // Accessors
    // =========================================================================
    int get_n_qubits() const { return n_qubits; }
    int get_n_features() const { return n_features; }
    int get_depth() const { return depth; }

    // =========================================================================
    // Debug: print circuit structure
    // =========================================================================
    void print_circuit_info() const {
        cout << "QuantumKernel circuit:\n";
        cout << "  Qubits: " << n_qubits << "\n";
        cout << "  Input features: " << n_features << "\n";
        cout << "  Circuit depth: " << depth << " layers\n";
        cout << "  Each layer: RY rotations (parameterized by x) + CNOT entanglement\n";
        cout << "  Total parameters: " << depth * n_qubits * (n_features + 1) << "\n";
        cout << "  Hilbert space dimension: " << (1 << n_qubits) << "\n";
    }

    // =========================================================================
    // Toy example: compute kernel between two specific points
    // =========================================================================
    void demo_kernel_pair() const {
        vector<double> x1(n_features, 0.0);
        vector<double> x2(n_features, 0.1);
        
        double k12 = kernel(x1, x2);
        double k11 = kernel(x1, x1);
        double k22 = kernel(x2, x2);
        
        cout << "Demo kernel values:\n";
        cout << "  K_Q(x1, x1) = " << k11 << " (should be 1.0 for pure state)\n";
        cout << "  K_Q(x2, x2) = " << k22 << " (should be 1.0 for pure state)\n";
        cout << "  K_Q(x1, x2) = " << k12 << " (similarity measure)\n";
    }
};