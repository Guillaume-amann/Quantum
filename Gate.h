#pragma once
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <optional>

using namespace std;
using Complex = complex<double>;

class Gate {
public:
    vector<vector<complex<double>>> matrix;
    int size;  // number of qubits the gate acts on

    Gate(int num_qubits) : size(num_qubits) {
        int dim = 1 << num_qubits;
        matrix.resize(dim, vector<complex<double>>(dim, {0.0, 0.0}));
    }

    Gate(const vector<vector<complex<double>>>& mat) : matrix(mat) {
        int dim = mat.size();
        if ((dim & (dim - 1)) != 0) throw invalid_argument("Matrix size must be a power of 2");
        size = log2(dim);
    }

    static Gate identity(int n) {
        int dim = 1 << n;
        Gate g(n);
        for (int i = 0; i < dim; ++i)
            g.matrix[i][i] = {1.0, 0.0};
        return g;
    }

    static Gate tensor(const Gate& A, const Gate& B) {
        int dimA = 1 << A.size;
        int dimB = 1 << B.size;
        int dim = dimA * dimB;
        Gate result(A.size + B.size);

        for (int i = 0; i < dimA; ++i) {
            for (int j = 0; j < dimA; ++j) {
                for (int k = 0; k < dimB; ++k) {
                    for (int l = 0; l < dimB; ++l) {
                        result.matrix[i * dimB + k][j * dimB + l] = A.matrix[i][j] * B.matrix[k][l];
                    }
                }
            }
        }
        return result;
    }

    static Gate X() {
        return Gate({
            {0, 1},
            {1, 0}
        });
    }

    static Gate Y() {
        return Gate({
            {0, {-0, 1}},
            {{0, -1}, 0}
        });
    }

    static Gate Z() {
        return Gate({
            {1, 0},
            {0, -1}
        });
    }

    static Gate H() {
        complex<double> inv_sqrt2 = 1.0 / sqrt(2);
        return Gate({
            {inv_sqrt2, inv_sqrt2},
            {inv_sqrt2, -inv_sqrt2}
        });
    }

    static Gate Rx(double theta) {
        complex<double> c = cos(theta / 2);
        complex<double> is = {0, -sin(theta / 2)};
        return Gate({
            {c, is},
            {is, c}
        });
    }

    static Gate Ry(double theta) {
        double c = cos(theta / 2);
        double s = sin(theta / 2);
        return Gate({
            {c, -s},
            {s, c}
        });
    }

    static Gate Rz(double theta) {
        complex<double> minus = exp(complex<double>(0, -theta / 2));
        complex<double> plus = exp(complex<double>(0, theta / 2));
        return Gate({
            {minus, 0},
            {0, plus}
        });
    }

    // Extend Gate class with static method rotationX acting on target qubit
    static Gate rotationX(int target, int total_qubits, double theta) {
        Gate rx = Gate::Rx(theta);
        return rx.expand(total_qubits, target);
    }

    // Extend this gate to n qubits, acting on target qubit index
    Gate expand(int total_qubits, int target) const {
        if (size != 1) throw invalid_argument("Only single-qubit gates can be expanded");

        Gate result = identity(0);
        for (int i = 0; i < total_qubits; ++i) {
            if (i == target)
                result = tensor(result, *this);
            else
                result = tensor(result, identity(1));
        }
        return result;
    }
};

// Matrix multiplication helper for Gate multiplication
Gate operator*(const Gate& A, const Gate& B) {
    if (A.matrix.size() != B.matrix.size())
        throw invalid_argument("Gate size mismatch for multiplication");
    int dim = A.matrix.size();
    Gate result(A.size);
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j) {
            Complex sum = 0;
            for (int k = 0; k < dim; ++k)
                sum += A.matrix[i][k] * B.matrix[k][j];
            result.matrix[i][j] = sum;
        }
    return result;
}