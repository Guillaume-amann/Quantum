#pragma once
#include <complex>
#include <cmath>
#include <optional>
#include "Qbit.h"

using namespace std;

class Gate {
public:
    // 2x2 matrix elements: m[row][col]
    complex<double> m[2][2];

    // Constructor with matrix elements (row-major)
    Gate(complex<double> m00, complex<double> m01,
         complex<double> m10, complex<double> m11) {
        m[0][0] = m00; m[0][1] = m01;
        m[1][0] = m10; m[1][1] = m11;
    }

    // Apply this gate to a Qbit and return the transformed Qbit
    Qbit apply(const Qbit& q, optional<int> target = nullopt) const {
        int n = q.num_qubits();
        
        int actual_target = 0;
        if (target.has_value()) {
            actual_target = target.value();
        } else if (n == 1) {
            actual_target = 0;
        } else {
            throw invalid_argument("Target qubit index must be specified for multi-qubit systems.");
        }

        int dim = 1 << n;
        vector<complex<double>> new_state(dim, {0.0, 0.0});

        for (int i = 0; i < dim; ++i) {
            int bit = (i >> (n - 1 - actual_target)) & 1;
            int j = i ^ (1 << (n - 1 - actual_target));
            new_state[i] += m[bit][0] * q.getState()[i] + m[bit][1] * q.getState()[j];
        }

        return Qbit(n, new_state);
    }

    // Common 1-qubit gates:
    static Gate I() {
        return Gate(complex<double>(1,0), complex<double>(0,0),
                    complex<double>(0,0), complex<double>(1,0));
    }
    static Gate X() {
        return Gate(complex<double>(0,0), complex<double>(1,0),
                    complex<double>(1,0), complex<double>(0,0));
    }
    static Gate Y() {
        return Gate(complex<double>(0,0), complex<double>(0,-1),
                    complex<double>(0,1), complex<double>(0,0));
    }
    static Gate Z() {
        return Gate(complex<double>(1,0), complex<double>(0,0),
                    complex<double>(0,0), complex<double>(-1,0));
    }
    static Gate H() {
        complex<double> invSqrt2 = 1.0 / sqrt(2);
        return Gate(invSqrt2, invSqrt2,
                    invSqrt2, -invSqrt2);
    }
    static Gate S() {
        return Gate(complex<double>(1,0), complex<double>(0,0),
                    complex<double>(0,0), complex<double>(0,1));
    }
    static Gate T() {
        complex<double> exp_i_pi_4 = exp(complex<double>(0, M_PI/4));
        return Gate(complex<double>(1,0), complex<double>(0,0),
                    complex<double>(0,0), exp_i_pi_4);
    }

    // Rotation around X,Y,Z axis by angle theta (radians)
    static Gate Rx(double theta) {
        complex<double> c = cos(theta/2);
        complex<double> is = complex<double>(0, -sin(theta/2));
        return Gate(c, is,
                    is, c);
    }
    static Gate Ry(double theta) {
        double c = cos(theta/2);
        double s = sin(theta/2);
        return Gate(c, complex<double>(-s, 0),
                    complex<double>(s, 0), c);
    }
    static Gate Rz(double theta) {
        complex<double> exp_minus = exp(complex<double>(0, -theta/2));
        complex<double> exp_plus = exp(complex<double>(0, theta/2));
        return Gate(exp_minus, 0, 0, exp_plus);
    }
};