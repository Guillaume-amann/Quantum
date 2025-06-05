class Gate {
public:
    virtual void apply(Qbit& q) = 0; // Virtual method to apply gate
};

class HadamardGate : public Gate {
public:
    void apply(Qbit& q) override {
        // Hadamard gate matrix
        std::vector<std::vector<std::complex<double>>> hadamard = {
            {1.0 / std::sqrt(2), 1.0 / std::sqrt(2)},
            {1.0 / std::sqrt(2), -1.0 / std::sqrt(2)}
        };
        q.applyGate(hadamard);
    }
};