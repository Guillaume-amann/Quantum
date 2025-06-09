#pragma once
#include <vector>
#include "Qbit.h"
#include "Gate.h"

using namespace std;

class Circuit {
private:
    vector<Qbit> qbits;
    vector<Gate*> gates;

public:
    void addQbit(Qbit q) { qbits.push_back(q); }

    void addGate(Gate* g) { gates.push_back(g); }

    void run() {
        for (auto& gate : gates) {
            for (auto& qbit : qbits) {
                gate->apply(qbit);
            }
        }
    }

    void printCircuitState() {
        for (auto& qbit : qbits) {
            qbit.getQstate();
        }
    }
};