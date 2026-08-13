#pragma once
#include <vector>
#include <complex>
#include <random>
#include <Eigen/Dense>
#include "Gate.h"

using namespace std;
using namespace Eigen;
using Complex = complex<double>;
 
// =============================================================================
// NoiseModel — ready to pass to:
//   DensityMatrix::apply_kraus(NoiseModel::xxx(...))
//
// All channels are returned as full 2^n × 2^n operators acting on the
// target qubit inside an n-qubit register (identity on all other qubits)
//
// Physical background
// -------------------
// A noise channel is a completely positive, trace-preserving (CPTP) map:
//
//   ε(ρ) = Σ_k  K_k . ρ . K_k†       with  Σ_k K_k† K_k = I  (completeness)
//
// Each K_k is a Kraus operator. Unitary evolution is the special case where
// there is only one Kraus operator (K_0 = U).
//
// Naming convention for physical parameters
// -----------------------------------------
//   p          — error probability per gate application (dimensionless, 0 ≤ p ≤ 1)
//   T1         — energy relaxation time (same unit as t_gate, typically ns)
//   T2         — dephasing time, T2 ≤ 2·T1  always  (same unit as t_gate)
//   t_gate     — gate duration (same unit as T1, T2)
//   target     — qubit index (0-based, big-endian convention from Gate.h)
//   n          — total number of qubits in the register
//
// =============================================================================
 
class NoiseModel {
private:
    // builds I ⊗ ... ⊗ K ⊗ ... ⊗ I.
    static Gate expand1(const Gate& K, int target, int n) { return K.expand(n, target); }
 
    static Gate scale(const Gate& G, double s) {
        Gate out = G;
        out.matrix *= Complex(s, 0.0);
        return out;
    }

    static Gate gate2(Complex a, Complex b, Complex c, Complex d) {
        Gate g(1);
        g.matrix(0,0) = a;  g.matrix(0,1) = b;
        g.matrix(1,0) = c;  g.matrix(1,1) = d;
        return g;
    }
 
public:
    // =========================================================================
    // 1. Bit-flip channel
    //
    // Physical meaning: with probability p, a qubit spontaneously flips |0⟩↔|1⟩.
    // This is the quantum analogue of a classical bit-flip error.
    //
    // Kraus operators:
    //   K_0 = √(1-p) · I    (no error)
    //   K_1 = √p     · X    (bit flip)
    //
    // Completeness: K_0†K_0 + K_1†K_1 = (1-p)I + p·X†X = (1-p)I + pI = I  ✓
    //
    // Effect on Bloch sphere: shrinks the X and Y components, Z unchanged.
    //   ρ' = (1-p)ρ + p·XρX†
    // =========================================================================
    static vector<Gate> bit_flip(double p, int target, int n) {
        if (target < 0 || target >= n) {
            throw invalid_argument("Qubit index " + to_string(target) + 
                                " out of bounds for " + to_string(n) + 
                                "-qubit system");
        }
        if (p < 0.0 || p > 1.0) throw invalid_argument("Bit-flip probability p must be in [0, 1]");
        if (p < 0 || p > 1) throw invalid_argument("Bit_flip: p must be in [0,1]");
        return {
            expand1(scale(Gate::identity(1), sqrt(1.0 - p)), target, n),
            expand1(scale(Gate::X(),         sqrt(p)),       target, n)
        };
    }
 
    // =========================================================================
    // 2. Phase-flip (dephasing) channel
    //
    // Physical meaning: with probability p, the qubit acquires a random π
    // phase kick, destroying coherence between |0⟩ and |1⟩.
    // This is the most common error in real superconducting hardware.
    //
    // Kraus operators:
    //   K_0 = √(1-p) · I
    //   K_1 = √p     · Z
    //
    // Completeness: K_0†K_0 + K_1†K_1 = (1-p)I + p·Z†Z = I  ✓
    //
    // Effect on density matrix:
    //   ρ' = (1-p)ρ + p·ZρZ†
    //   Off-diagonal elements: ρ'[0][1] = (1-2p) · ρ[0][1]
    //   → At p=0.5: complete dephasing (off-diag → 0), diagonal unchanged.
    //   → At p=1:   coherence flips sign (not destroyed — this is a pure Z gate).
    // =========================================================================
    static vector<Gate> phase_flip(double p, int target, int n) {
        if (target < 0 || target >= n) {
            throw invalid_argument("Qubit index " + to_string(target) + 
                                " out of bounds for " + to_string(n) + 
                                "-qubit system");
        }
        if (p < 0.0 || p > 1.0) throw invalid_argument("Phase-flip probability p must be in [0, 1]");
        if (p < 0 || p > 1) throw invalid_argument("Phase_flip: p must be in [0,1]");
        return {
            expand1(scale(Gate::identity(1), sqrt(1.0 - p)), target, n),
            expand1(scale(Gate::Z(),         sqrt(p)),       target, n)
        };
    }
 
    // =========================================================================
    // 3. Bit-phase-flip channel (Y error)
    //
    // Physical meaning: simultaneous bit and phase flip. Rarely modelled alone,
    // but included for completeness and for building depolarising below.
    //
    // Kraus operators:
    //   K_0 = √(1-p) · I
    //   K_1 = √p     · Y
    // =========================================================================
    static vector<Gate> bit_phase_flip(double p, int target, int n) {
        if (p < 0 || p > 1) throw invalid_argument("Bit_phase_flip: p must be in [0,1]");
        return {
            expand1(scale(Gate::identity(1), sqrt(1.0 - p)), target, n),
            expand1(scale(Gate::Y(),         sqrt(p)),       target, n)
        };
    }
 
    // =========================================================================
    // 4. Depolarising channel
    //
    // Physical meaning: the standard model for gate infidelity in real hardware.
    // With probability p, the qubit is replaced by the maximally mixed state I/2.
    // Equivalently: with probability p/3 each, a random X, Y, or Z error occurs.
    //
    // Kraus operators (4 operators):
    //   K_0 = √(1 - p)    · I
    //   K_1 = √(p/3)      · X
    //   K_2 = √(p/3)      · Y
    //   K_3 = √(p/3)      · Z
    //
    // Completeness: (1-p)I + (p/3)(X†X + Y†Y + Z†Z) = (1-p)I + pI = I  ✓
    //
    // Effect on Bloch sphere: uniform shrinkage toward the origin.
    //   Bloch vector r → (1 - 4p/3) · r
    //   → At p=0.75: any state maps to I/2 (maximally mixed).
    //   → p ∈ [0, 0.75] for this to be a valid (contractive) channel.
    //
    // This is the channel reported in hardware specs as "gate fidelity F = 1 - p".
    // =========================================================================
    static vector<Gate> depolarising(double p, int target, int n) {
        if (target < 0 || target >= n) {
            throw invalid_argument("Qubit index " + to_string(target) + 
                                " out of bounds for " + to_string(n) + 
                                "-qubit system");
        }
        if (p < 0.0 || p > 1.0) throw invalid_argument("Depolarising: p must be in [0, 1]");
        if (p < 0 || p > 0.75) throw invalid_argument("Depolarising: p must be in [0, 0.75]");
        double s0 = sqrt(1.0 - p);
        double s1 = sqrt(p / 3.0);
        return {
            expand1(scale(Gate::identity(1), s0), target, n),
            expand1(scale(Gate::X(),         s1), target, n),
            expand1(scale(Gate::Y(),         s1), target, n),
            expand1(scale(Gate::Z(),         s1), target, n)
        };
    }
 
    // =========================================================================
    // 5. Amplitude damping channel
    //
    // Physical meaning: energy relaxation — the qubit spontaneously decays from
    // |1⟩ to |0⟩ (like an excited atom emitting a photon). This is the primary
    // source of error in superconducting qubits and is parameterised by T1.
    //
    // The decay probability over a gate of duration t_gate is:
    //   γ = 1 - exp(-t_gate / T1)
    //
    // Kraus operators:
    //   K_0 = [ 1        0      ]    (no decay, ground state preserved)
    //         [ 0   √(1-γ)     ]
    //
    //   K_1 = [ 0   √γ         ]    (decay |1⟩ → |0⟩)
    //         [ 0    0         ]
    //
    // Completeness:
    //   K_0†K_0 + K_1†K_1
    //   = diag(1, 1-γ) + diag(0, γ)
    //   = diag(1, 1) = I  ✓
    //
    // Effect on density matrix:
    //   ρ[0][0]' = ρ[0][0] + γ·ρ[1][1]      (ground pop. grows)
    //   ρ[1][1]' = (1-γ)·ρ[1][1]            (excited pop. decays)
    //   ρ[0][1]' = √(1-γ) · ρ[0][1]         (coherence damps)
    // =========================================================================
    static vector<Gate> amplitude_damping(double T1, double t_gate, int target, int n) {
        if (target < 0 || target >= n) {
            throw invalid_argument("Qubit index " + to_string(target) + 
                                " out of bounds for " + to_string(n) + 
                                "-qubit system");
        }
        if (T1 <= 0) throw invalid_argument("Amplitude_damping: T1 must be > 0");
        if (t_gate < 0) throw invalid_argument("Amplitude_damping: t_gate must be >= 0");
 
        double gamma = 1.0 - exp(-t_gate / T1);   // decay probability
 
        Gate K0 = gate2(1.0,         0.0,
                        0.0, sqrt(1.0 - gamma));
        Gate K1 = gate2(0.0, sqrt(gamma),
                        0.0,         0.0);
        return {
            expand1(K0, target, n),
            expand1(K1, target, n)
        };
    }
 
    // =========================================================================
    // 6. Generalised amplitude damping channel
    //
    // Physical meaning: amplitude damping at finite temperature. At T=0, the
    // qubit can only decay (|1⟩→|0⟩). At finite temperature, it can also be
    // thermally excited (|0⟩→|1⟩). This is the physically correct model for
    // any hardware operating above absolute zero.
    //
    // Parameters:
    //   T1     — energy relaxation time
    //   t_gate — gate duration
    //   p_eq   — equilibrium excited-state population  p_eq = 1/(exp(ℏω/kT)+1)
    //            → 0 at T=0 (pure decay), → 0.5 at T=∞ (maximally mixed)
    //
    // Decay probability:   γ = 1 - exp(-t_gate / T1)
    //
    // Kraus operators (4 operators):
    //   K_0 = √(1-p_eq) · [ 1        0    ]
    //                      [ 0   √(1-γ)   ]
    //
    //   K_1 = √(1-p_eq) · [ 0   √γ       ]
    //                      [ 0    0       ]
    //
    //   K_2 = √p_eq     · [ √(1-γ)   0   ]
    //                      [ 0        1   ]
    //
    //   K_3 = √p_eq     · [ 0     0      ]
    //                      [ √γ    0      ]
    //
    // Completeness: verified by (1-p_eq)(K_0†K_0 + K_1†K_1) + p_eq(K_2†K_2 + K_3†K_3) = I ✓
    //
    // Steady state: ρ_ss = diag(1-p_eq, p_eq)  — the thermal state.
    // =========================================================================
    static vector<Gate> generalised_amplitude_damping(double T1, double t_gate, double p_eq, int target, int n) {
        if (target < 0 || target >= n) {
            throw invalid_argument("Qubit index " + to_string(target) + 
                                " out of bounds for " + to_string(n) + 
                                "-qubit system");
        }
        if (p_eq < 0.0 || p_eq > 1.0) throw invalid_argument("Generalised amplitude damping probability p_eq must be in [0, 1]");
        if (T1 <= 0)           throw invalid_argument("Generalised_amplitude_damping: T1 must be > 0");
        if (t_gate < 0)        throw invalid_argument("Generalised_amplitude_damping: t_gate must be >= 0");
        if (p_eq < 0 || p_eq > 0.5)
            throw invalid_argument("Generalised_amplitude_damping: p_eq must be in [0, 0.5]");
 
        double gamma  = 1.0 - exp(-t_gate / T1);
        double s0 = sqrt(1.0 - p_eq);
        double s1 = sqrt(p_eq);
        double sg = sqrt(gamma);
        double sc = sqrt(1.0 - gamma);
 
        Gate K0 = scale(gate2(1.0,  0.0,
                               0.0,  sc),  s0);
        Gate K1 = scale(gate2(0.0,  sg,
                               0.0,  0.0), s0);
        Gate K2 = scale(gate2(sc,   0.0,
                               0.0,  1.0), s1);
        Gate K3 = scale(gate2(0.0,  0.0,
                               sg,   0.0), s1);
        return {
            expand1(K0, target, n),
            expand1(K1, target, n),
            expand1(K2, target, n),
            expand1(K3, target, n)
        };
    }
 
    // =========================================================================
    // 7. Combined T1/T2 channel  (the standard hardware noise model)
    //
    // Physical meaning: real qubits suffer both energy relaxation (T1) and
    // dephasing (T2). These are NOT independent: dephasing includes a
    // contribution from relaxation, so T2 ≤ 2·T1 always.
    //
    // The pure dephasing time T_phi is extracted as:
    //   1/T2 = 1/(2·T1) + 1/T_phi
    //   T_phi = 1 / (1/T2 - 1/(2·T1))
    //
    // Implementation: compose amplitude damping (captures T1) and pure
    // dephasing (captures the remaining T_phi contribution).
    //
    //   Pure dephasing Kraus operators:
    //     K_0 = [ 1            0        ]
    //           [ 0   √(1 - λ)          ]       λ = 1 - exp(-t_gate / T_phi)
    //
    //     K_1 = [ 0       0             ]
    //           [ 0   √λ                ]
    //
    // These two channels are applied sequentially; their composition gives
    // the correct Lindblad evolution for both T1 and T2 simultaneously.
    //
    // Completeness of each sub-channel: verified independently above.
    //
    // Returned as a struct so the caller can apply them in sequence:
    //   auto [ad, pd] = NoiseModel::t1_t2(T1, T2, t_gate, q, n);
    //   dm.apply_kraus(ad);
    //   dm.apply_kraus(pd);
    // =========================================================================
    struct T1T2Channels {
        vector<Gate> amplitude_damping_ops;
        vector<Gate> pure_dephasing_ops;
    };
 
    static T1T2Channels t1_t2(double T1, double T2, double t_gate, int target, int n) {
        if (T1 <= 0) throw invalid_argument("T1T2Channels::t1_t2: T1 must be > 0");
        if (T2 <= 0) throw invalid_argument("T1T2Channels::t1_t2: T2 must be > 0");
        if (T2 > 2.0 * T1) throw invalid_argument("T1T2Channels::t1_t2: T2 > 2·T1 violates the physical constraint");
        if (t_gate < 0) throw invalid_argument("T1T2Channels::t1_t2: t_gate must be >= 0");
 
        // Amplitude damping (T1 process)
        auto ad_ops = amplitude_damping(T1, t_gate, target, n);
 
        // Pure dephasing rate: 1/T_phi = 1/T2 - 1/(2·T1)
        double rate_phi = 1.0 / T2 - 1.0 / (2.0 * T1);
        vector<Gate> pd_ops;
 
        if (rate_phi < 1e-15) {
            // T2 = 2·T1: no pure dephasing beyond what T1 already causes.
            // Pure dephasing channel = identity (one Kraus op = I).
            pd_ops = { Gate::identity(n) };
        } else {
            double T_phi = 1.0 / rate_phi;
            double lam   = 1.0 - exp(-t_gate / T_phi);   // pure dephasing probability
 
            Gate Kpd0 = gate2(1.0,          0.0,
                               0.0, sqrt(1.0 - lam));
            Gate Kpd1 = gate2(0.0,  0.0,
                               0.0, sqrt(lam));
            pd_ops = {
                expand1(Kpd0, target, n),
                expand1(Kpd1, target, n)
            };
        }
 
        return { ad_ops, pd_ops };
    }
 
    // =========================================================================
    // 8. Two-qubit depolarising channel
    //
    // Physical meaning: the standard noise model for two-qubit gates (CNOT, CZ).
    // Two-qubit gates are typically 10-100× noisier than single-qubit gates in
    // real hardware. Models uniform Pauli error on both qubits.
    //
    // Kraus operators: 16 operators  { √(1-p)·I⊗I } ∪ { √(p/15)·P_i⊗P_j }
    // where P ∈ {I, X, Y, Z} and the identity⊗identity term is excluded from
    // the error set (it is absorbed into the (1-p) term).
    //
    // Completeness: (1-p)·I + (p/15)·15·I = I  ✓  (since each P_i⊗P_j is unitary)
    //
    // Note: the two target qubits are specified as (q0, q1); the gate acts on
    // the 4-dimensional subspace of the n-qubit register.
    // =========================================================================
    static vector<Gate> depolarising_2q(double p, int q0, int q1, int n) {
        if (q0 < 0 || q0 >= n) {
            throw invalid_argument("Qubit index " + to_string(q0) + 
                                " out of bounds for " + to_string(n) + 
                                "-qubit system");
        }
        if (q1 < 0 || q1 >= n) {
            throw invalid_argument("Qubit index " + to_string(q1) + 
                                " out of bounds for " + to_string(n) + 
                                "-qubit system");
        }
        if (p < 0.0 || p > 1.0) throw invalid_argument("Depolarising_2q: p must be in [0, 1]");
        if (p < 0 || p > 1) throw invalid_argument("Depolarising_2q: p must be in [0,1]");
        if (q0 == q1) throw invalid_argument("Depolarising_2q: q0 and q1 must differ");
 
        vector<Gate> paulis_1q = {
            Gate::identity(1), Gate::X(), Gate::Y(), Gate::Z()
        };
 
        vector<Gate> result;
        result.reserve(16);
 
        double s_id  = sqrt(1.0 - p);
        double s_err = sqrt(p / 15.0);
 
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                bool is_II = (i == 0 && j == 0);
                double s = is_II ? s_id : s_err;
 
                // Build the two-qubit Pauli P_i ⊗ P_j in the full n-qubit space.
                // Apply P_i to q0 and P_j to q1 independently (they commute on
                // different qubits), then compose.
                Gate full_i = paulis_1q[i].expand(n, q0);
                Gate full_j = paulis_1q[j].expand(n, q1);
                Gate combined = scale(full_i * full_j, s);
                result.push_back(combined);
            }
        }
        return result;
    }
 
    // =========================================================================
    // 9. Measurement error channel
    //
    // Physical meaning: readout errors — the classical outcome is flipped with
    // probability p_01 (reporting 1 when the qubit is |0⟩) or p_10 (reporting
    // 0 when the qubit is |1⟩). Models imperfect photon detection in
    // superconducting readout circuits.
    //
    // This is NOT a Kraus channel on ρ — it is a classical confusion matrix
    // applied to the probability distribution after measurement.
    //
    // Returned as a 2×2 matrix M where M[b][a] = P(report b | true outcome a).
    //
    //   M = [ 1-p_01    p_10  ]
    //       [  p_01   1-p_10  ]
    //
    // Usage:
    //   auto M = NoiseModel::measurement_error(p_01, p_10);
    //   auto probs = dm.probabilities();
    //   // For a single qubit: p_noisy[b] = M[b][0]*probs[0] + M[b][1]*probs[1]
    // =========================================================================
    struct MeasurementError {
        double p_01;   // P(report 1 | qubit is |0⟩)
        double p_10;   // P(report 0 | qubit is |1⟩)
 
        // Apply confusion matrix to a probability vector (full register).
        // For a single target qubit: marginalise over all other qubits.
        // Returns the distorted probability distribution — does NOT collapse ρ.
        vector<double> apply(const vector<double>& probs, int target, int n) const {
            int dim = 1 << n;
            int k_shift = n - 1 - target;
            vector<double> out(dim, 0.0);
            for (int i = 0; i < dim; ++i) {
                int b = (i >> k_shift) & 1;   // true outcome bit
                for (int b_rep = 0; b_rep < 2; ++b_rep) {
                    // M[b_rep][b]
                    double m = (b_rep == b) ? (b == 0 ? 1.0 - p_01 : 1.0 - p_10)
                                            : (b == 0 ? p_01       : p_10);
                    // Flip target bit in i to get the reported index
                    int i_rep = (i & ~(1 << k_shift)) | (b_rep << k_shift);
                    out[i_rep] += m * probs[i];
                }
            }
            return out;
        }

        // Apply the confusion matrix then draw a classical outcome.
        // Returns the reported bit (0 or 1) for the target qubit.
        // Marginalises over all other qubits before sampling.
        //
        // Usage in QuantumSim.cpp:
        //   auto me = NoiseModel::measurement_error(0.01, 0.02);
        //   int bit = me.sample(dm.probabilities(), target, n);
        int sample(const vector<double>& probs, int target, int n) const {
            vector<double> noisy = apply(probs, target, n);

            // Marginalise over all qubits except target
            int k_shift = n - 1 - target;
            double p0 = 0.0;
            for (int i = 0; i < (int)noisy.size(); ++i)
                if (((i >> k_shift) & 1) == 0) p0 += noisy[i];

            static thread_local mt19937 gen(random_device{}());
            uniform_real_distribution<double> dist(0.0, 1.0);
            return (dist(gen) < p0) ? 0 : 1;
        }
    };
 
    static MeasurementError measurement_error(double p_01, double p_10) {
        if (p_01 < 0 || p_01 > 1) throw invalid_argument("measurement_error: p_01 must be in [0,1]");
        if (p_10 < 0 || p_10 > 1) throw invalid_argument("measurement_error: p_10 must be in [0,1]");
        return { p_01, p_10 };
    }
};