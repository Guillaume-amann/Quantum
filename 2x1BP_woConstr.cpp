#include <iostream>
#include <limits>
#include <cmath>
#include "Qbit.h"
#include "Gate.h"

using namespace std;


// If cost = sum of independent item costs, no constraints, yes — no entanglement strictly needed


// Run one QAOA step with given gamma and beta
string run_qaoa_step(double gamma, double beta) {
    Qbit q(2);  // 2 qubits

    // Apply Hadamard to all qubits to create superposition
    Gate h = Gate::H();
    q = h.apply(q, 0);  // qubit 0
    q = h.apply(q, 1);  // qubit 1

    // Cost unitary: Rz on each qubit with angle 2*gamma
    Gate rz = Gate::Rz(2 * gamma);
    q = rz.apply(q, 0);
    q = rz.apply(q, 1);

    // Mixer unitary: Rx on each qubit with angle 2*beta
    Gate rx = Gate::Rx(2 * beta);
    q = rx.apply(q, 0);
    q = rx.apply(q, 1);

    string measured_ket = q.measure();  // returns e.g. "|01>"
    return measured_ket;  // change your function return type accordingly
}

int count_packed_items(const std::string& ket) {
    int count = 0;
    for (char c : ket) {
        if (c == '1') ++count;
    }
    return count;
}

double estimate_expectation(double gamma, double beta, int shots = 1000) {
    int sum = 0;

    for (int i = 0; i < shots; ++i) {
        string result = run_qaoa_step(gamma, beta);
        int packed = count_packed_items(result);
        sum += packed;  // reward per packed item
    }

    return static_cast<double>(sum) / shots;
}

int main() {
    cout << "QAOA for 2x1 Bin Packing (2 qubits)" << endl;

    double best_gamma = 0.0;
    double best_beta = 0.0;
    double best_expectation = -1.0;

    for (double gamma = 0.0; gamma <= M_PI; gamma += 0.1) {
        for (double beta = 0.0; beta <= M_PI / 2; beta += 0.1) {
            double expectation = estimate_expectation(gamma, beta);

            if (expectation > best_expectation) {
                best_expectation = expectation;
                best_gamma = gamma;
                best_beta = beta;
            }
        }
    }

    cout << "\nOptimal parameters:" << endl;
    cout << "γ = " << best_gamma << ", β = " << best_beta
         << " → max ⟨#packed⟩ = " << best_expectation << endl;  

    // Run a single measurement at optimal γ, β
    string final_result = run_qaoa_step(best_gamma, best_beta);
    cout << "Measurement for optimal (γ, β): " << final_result << endl;

    return 0;
}