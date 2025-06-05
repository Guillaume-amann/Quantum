#include <iostream>
#include <cmath>

#include "Qbit.h"
#include "Gate.h"
#include "Circuit.h"

int main() {
    Qbit q;
    q.printState(); // Should print |0⟩: 1, |1⟩: 0

    HadamardGate hadamard; // Create a Hadamard gate
    hadamard.apply(q); // Apply the gate
    q.printState(); // Prints the new superposition state

    int measurement = q.measure(); // Measure the qubit
    cout << "Measurement result: " << measurement << endl;
    q.printState(); // Prints the collapsed state
}