#include "Qbit.h"
#include "Gate.h"

int main() {
    Qbit q(2); // 2-qubit system initialized to |00>
    Gate h = Gate::H();
    q = h.apply(q, 0);  // Apply Hadamard to qubit 0
    q = h.apply(q, 1);  // Apply Hadamard to qubit 1
    q.print_state();
    string measure = q.measure();
    cout << measure << endl;
}