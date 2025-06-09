#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <numeric>  // for accumulate

using namespace std;

class Qbit {
private:
    int n;  // number of qubits
    vector<complex<double>> state; // 2^n complex amplitudes

    void normalize() {
        double norm_sq = 0.0;
        for (auto& amp : state)
            norm_sq += norm(amp);

        if (norm_sq == 0.0) {
            // Reset to |0...0>
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
    state[0] = {1.0, 0.0}; // set initial amplitude of |0...0> state to 1
}

    Qbit(int num_qubits, const vector<complex<double>>& initial_state)
        : n(num_qubits), state(initial_state) {
        if (state.size() != (1 << n)) {
            throw invalid_argument("Initial state vector size mismatch");
        }
        normalize();
    }

    // Measure the multi-qubit state, returning the classical bitstring as an integer
    string measure() {
        static thread_local std::mt19937 gen(std::random_device{}());
        uniform_real_distribution<double> dist(0.0, 1.0);

        // Compute cumulative probabilities
        vector<double> cumulative(state.size(), 0.0);
        cumulative[0] = std::norm(state[0]);
        for (size_t i = 1; i < state.size(); ++i) {
            cumulative[i] = cumulative[i - 1] + std::norm(state[i]);
        }

        double r = dist(gen);
        size_t idx = 0;
        while (idx < cumulative.size() && r > cumulative[idx]) ++idx;
        if (idx == cumulative.size()) idx = cumulative.size() - 1;

        state.assign(state.size(), {0.0, 0.0});
        state[idx] = {1.0, 0.0};

        // Convert index to ket string
        string ket = "|";
        for (int i = n - 1; i >= 0; --i) {
            ket += ((idx >> i) & 1) ? '1' : '0';
        }
        ket += ">";
        return ket;
    }

    static int ketToIndex(const string& ket);

    // Print full quantum state
    void print_state() const {
        for (size_t i = 0; i < state.size(); ++i) {
            if (norm(state[i]) > 1e-10) {
                cout << "(" << state[i] << ") |" << bitset<16>(i).to_string().substr(16-n) << "⟩ + ";
            }
        }
        cout << "\b\b \n";  // Erase last "+ "
    }

    int num_qubits() const { return n; }

    const vector<complex<double>>& getState() const { return state; }
};