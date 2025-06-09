#pragma once
#include <iostream>
#include <cmath>
#include <random>
#include <complex>

using namespace std;

class Qbit {

private:
    // Ensure normalization: |α|² + |β|² = 1
    void normalize() {
        double norma = sqrt(norm(alpha) + norm(beta));
        if (norma == 0) {
            alpha = 1.0;
            beta = 0.0;
        } else {
            alpha /= norma;
            beta /= norma;
        }
    }

public:
    complex<double> alpha;  // amplitude for |0⟩
    complex<double> beta;   // amplitude for |1⟩

    // Constructor to initialise the Quantum state's basis |0> at 1 (and |1> at 0)
    // This allow for Qbit() or Qbit(z1,z2) instantiation
    Qbit(complex<double> a = complex<double>(1.0, 0.0),
         complex<double> b = complex<double>(0.0, 0.0))
            : alpha(a), beta(b) { normalize(); }

                                         

    // Measure the Qbit (collapse the superposition and return one basis state)
    int measure() {
        static thread_local mt19937 gen(random_device{}());
        uniform_real_distribution<double> dist(0.0, 1.0);

        double P0 = std::norm(alpha);
        double randomValue = dist(gen);

        if (randomValue < P0) {
            alpha = {1.0, 0.0};
            beta = {0.0, 0.0};
            return 0;
        } else {
            alpha = {0.0, 0.0};
            beta = {1.0, 0.0};
            return 1;
        }
    }

    // Print quantum state in classical form: a|0⟩ + b|1⟩
    void getQstate() const {
        cout << alpha << " |0⟩ + " << beta << " |1⟩" << endl;
    }
};