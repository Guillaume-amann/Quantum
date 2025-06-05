#include <iostream>
#include <cmath>
#include <vector>
#include <complex>

using namespace std;

class Qbit {
public:

    vector<complex<double>> state; // State vector ( alpha )
                                                // ( beta  )  a,b belongs to C

    // Constructor to initialise the Quantum state's basis |0> at 1             (and |1> at 0)                                            
    Qbit() : state({{1.0, 0.0}, 
                    {0.0, 0.0}}) {}

    // Apply a gate to the qubit
    void applyGate(const vector<vector<complex<double>>>& gate) {
        vector<complex<double>> newState(2, {0.0, 0.0});
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                newState[i] += gate[i][j] * state[j];
            }
        }
        state = newState;
    }

    // Measure the Qbit (collapse the superposition and return one basis state)
    int measure() {
        double P0 = norm(state[0]);
        double P1 = norm(state[1]);

        double randomValue = static_cast<double>(rand()) / RAND_MAX;

        if (randomValue < P0) {
            state = {{1.0, 0.0}, {0.0, 0.0}}; // Collapse to |0⟩
            return 0;
        } else {
            state = {{0.0, 0.0}, {1.0, 0.0}}; // Collapse to |1⟩
            return 1;
        }
    }

    void getQstate() const {cout << "|0⟩: " << state[0] << ", |1⟩: " << state[1] << endl;}
};