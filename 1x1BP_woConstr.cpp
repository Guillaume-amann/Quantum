#include <iostream>
#include <limits>
#include <cmath>
#include "Qbit.h"
#include "Gate.h"

using namespace std;


// If cost = sum of independent item costs, no constraints, yes — no entanglement strictly needed


// Run one QAOA step with given gamma and beta
int run_qaoa_step(double gamma, double beta) {
    Qbit q;  // 1-qubit system initialized to |0>

    q = Gate::H().apply(q);  // Hadamard on qubit 0

    // Cost Hamiltonian: H_C = -2.5·I - 7.5·Z → U_C = Rz(-15γ)
    Gate cost_unitary = Gate::Rz(15 * gamma);
    q = cost_unitary.apply(q);

    // Mixer: Rx(2β) retrieval of THE solution
    Gate mixer = Gate::Rx(2 * beta);
    q = mixer.apply(q);

    // 0 = empty, 1 = filled (low energy = high profit)
    string measured_ket = q.measure();           // e.g. "|0>" or "|1>"
    return Qbit::ketToIndex(measured_ket);       // convert to int 0 or 1
}

double estimate_expectation(double gamma, double beta, int shots = 1000) {
    int sum = 0;

    for (int i = 0; i < shots; ++i) {
        int result = run_qaoa_step(gamma, beta);
        // Energy: +10 for empty (0)(this "cost much"), -5 for filled (1)
        sum += (result == 0 ? +100 : -5);
    }

    return static_cast<double>(sum) / shots;
}

int Qbit::ketToIndex(const string& ket) {
    int result = 0;
    for (char c : ket.substr(1, ket.size() - 2)) {  // remove '|' and '>'
        result = (result << 1) | (c == '1' ? 1 : 0);
    }
    return result;
}

int main() {
    cout << "QAOA for 1x1 Bin Packing (1 qubit)" << endl;

    double best_gamma = 0.0;
    double best_beta = 0.0;
    double best_energy = numeric_limits<double>::infinity();

    for (double gamma = 0.0; gamma <= M_PI; gamma += 0.1) {
        for (double beta = 0.0; beta <= M_PI / 2; beta += 0.1) {
            double energy = estimate_expectation(gamma, beta);
            // cout << "γ = " << gamma << ", β = " << beta
            //      << " → ⟨H⟩ = " << energy << endl;

            if (energy < best_energy) {
                best_energy = energy;
                best_gamma = gamma;
                best_beta = beta;
            }
        }
    }

    cout << "\nOptimal parameters:" << endl;
    cout << "γ = " << best_gamma << ", β = " << best_beta
         << " → min ⟨H⟩ = " << best_energy << endl;

    // Run a single measurement at optimal γ, β
    int final_result = run_qaoa_step(best_gamma, best_beta);
    cout << "Measurement for optimal (γ, β): "
         << (final_result == 0 ? "empty (0)" : "filled (1)") << endl;

    return 0;
}