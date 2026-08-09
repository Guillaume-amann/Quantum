#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <bitset>
#include <string>
#include <Eigen/Dense>
#include "Gate.h"

using namespace std;
using namespace Eigen;
using Complex = complex<double>;

class Qbit {
private:
    int n;                          // Number of qubits
    vector<complex<double>> state;  // State vector ψ : size = 2^n

    void normalize() {
        double norm_sq = 0.0;
        for (const auto& amp : state)
            norm_sq += norm(amp);

        if (norm_sq == 0.0) {
            state.assign(state.size(), {0.0, 0.0});
            state[0] = {1.0, 0.0};
        } else {
            double norm_factor = sqrt(norm_sq);
            for (auto& amp : state)
                amp /= norm_factor;
        }
    }

public:
    Qbit(int num_qubits = 1) : n(num_qubits), state(1 << n, {0.0, 0.0}) {
        state[0] = {1.0, 0.0};
    }

    Qbit(int num_qubits, const vector<complex<double>>& initial_state) : n(num_qubits), state(initial_state) {
        if (state.size() != (1 << n)) throw invalid_argument("Initial state vector size mismatch");
        normalize();
    }

    // Measure and collapse the state
    string measure() {
        // Thread-local RNG seeded once per thread to avoid contention
        static thread_local mt19937 gen(random_device{}());
        uniform_real_distribution<double> dist(0.0, 1.0);

        // Build a cumulative probability distribution from the squared amplitudes (Born rule)
        vector<double> cumulative(state.size());
        cumulative[0] = norm(state[0]);
        for (size_t i = 1; i < state.size(); ++i) {
            cumulative[i] = cumulative[i - 1] + norm(state[i]);
        }

        // Sample a random value and find which basis state it falls into (inverse CDF sampling)
        double r = dist(gen);
        size_t idx = 0;
        while (idx < cumulative.size() && r > cumulative[idx]) ++idx;
        if (idx == cumulative.size()) idx = cumulative.size() - 1;  // clamp in case of floating-point rounding

        // Collapse the state: zero all amplitudes then set the measured basis state to 1
        state.assign(state.size(), {0.0, 0.0});
        state[idx] = {1.0, 0.0};

        // Format the measured index in a Dirac ket string format : |XX...X⟩>
        string bitstring = "|";
        for (int i = n - 1; i >= 0; --i) {bitstring += ((idx >> i) & 1) ? '1' : '0';}
        return bitstring += "⟩";
    }

    // Apply a gate to Qbit state
    void apply(const Gate& g) {
        if (g.matrix.rows() != (int)state.size()) throw invalid_argument("Gate dimension does not match Qbit state size");
        // Map state vector memory into Eigen without copying, so gate matrix multiplication
        // can operate directly on the std::vector data using Eigen's linear algebra ops.
        Map<VectorXcd> sv(state.data(), state.size());
        VectorXcd result = g.matrix * sv;
        sv = result;
    }

    // Displays amplitudes without collapsing — a "cheat view" impossible on real hardware.
    void print_state() const {
        for (size_t i = 0; i < state.size(); ++i) {
            if (norm(state[i]) > 1e-10) {
                cout << "(" << state[i] << ") |"
                     << bitset<64>(i).to_string().substr(64 - n) << "⟩ + ";
            }
        }
        cout << "\b\b \n";
    }

    int num_qubits() const { return n; }
    const vector<complex<double>>& get_state() const { return state; }
    vector<complex<double>>& access_state() { return state; }
};