#include <iostream>
#include "Qbit.h"
#include "Gate.h"

using namespace std;

int main() {
    Qbit q(1);  // 1-qubit system initialized to |0⟩
    cout << "Initial state q:" << endl;
    q.getState();

    Gate h = Gate::H();
    Qbit q2 = h.apply(q, 0);  // Apply Hadamard on qubit 0
    cout << "After H on qubit 0 (q2):" << endl;
    q2.getState();  // Should show (|0⟩ + |1⟩)/√2 superposition

    string measurement = q2.measure();  // Returns ket string like "|0>" or "|1>"
    cout << "Measurement of q2: " << measurement << endl;

    Gate x = Gate::X();
    Qbit q3 = x.apply(q2, 0);  // Bit-flip on qubit 0
    cout << "After X on qubit 0 (q3):" << endl;
    q3.getState();

    // Rotation example:
    Gate rx = Gate::Rx(M_PI / 2);
    Qbit q4 = rx.apply(q, 0);
    cout << "After Rx(pi/2) on qubit 0 (q4):" << endl;
    q4.getState();

    return 0;
}