#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <bitset>
#include <stdexcept>
#include <string>
#include <Eigen/Dense>
#include "Gate.h"

using namespace std;
using Complex = complex<double>;

class Qbit {
private:
    int n;  // Number of qubits
    vector<complex<double>> state;  // State vector: size = 2^n

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

    Qbit(int num_qubits, const vector<complex<double>>& initial_state)
        : n(num_qubits), state(initial_state) {
        if (state.size() != (1 << n)) {
            throw invalid_argument("Initial state vector size mismatch");
        }
        normalize();
    }

    // Measure and collapse the state
    string measure() {
        static thread_local mt19937 gen(random_device{}());
        uniform_real_distribution<double> dist(0.0, 1.0);

        vector<double> cumulative(state.size());
        cumulative[0] = norm(state[0]);
        for (size_t i = 1; i < state.size(); ++i) {
            cumulative[i] = cumulative[i - 1] + norm(state[i]);
        }

        double r = dist(gen);
        size_t idx = 0;
        while (idx < cumulative.size() && r > cumulative[idx]) ++idx;
        if (idx == cumulative.size()) idx = cumulative.size() - 1;

        state.assign(state.size(), {0.0, 0.0});
        state[idx] = {1.0, 0.0};

        string bitstring = "|";
        for (int i = n - 1; i >= 0; --i) {
            bitstring += ((idx >> i) & 1) ? '1' : '0';
        }
        bitstring += ">";
        return bitstring;
    }

    // Apply a gate to Qbit state
    void apply(const Gate& g) {
        if (g.matrix.rows() != (int)state.size())
            throw invalid_argument("Gate dimension does not match Qbit state size");
        Eigen::Map<Eigen::VectorXcd> sv(state.data(), state.size());
        Eigen::VectorXcd result = g.matrix * sv;
        sv = result;
    }

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
    vector<complex<double>>& access_state() { return state; } // allows in-place gate ops
};