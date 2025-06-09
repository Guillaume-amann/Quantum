#include "Qbit.h"
#include "Gate.h"

int main() {
    Qbit q;  // Initialized to |0⟩
    q.getQstate();

    Gate h = Gate::H();
    Qbit q2 = h.apply(q);
    q2.getQstate();  // Should show superposition (|0> + |1>)/√2
    int q2_value = q2.measure();
    cout << q2_value << endl;

    Gate x = Gate::X();
    Qbit q3 = x.apply(q2);
    q3.getQstate();  // Bit-flip of q2 state

    // Rotation example:
    Gate rx = Gate::Rx(M_PI/2);
    Qbit q4 = rx.apply(q);
    q4.getQstate();

    return 0;
}