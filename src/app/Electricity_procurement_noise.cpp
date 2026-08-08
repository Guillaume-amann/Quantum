// =============================================================================
//  Electricity Procurement via Bin Packing — NOISY (open-system) QAOA
//
//  The density-matrix counterpart of Electricity_procurement.cpp. Identical
//  problem, identical parameters, evaluated at the SAME optimal QAOA angles the
//  noiseless run found — so the two programs are directly comparable on the same
//  input. The difference is the simulation back-end and what it can express:
//
//    Electricity_procurement.cpp        this file
//    -----------------------------      --------------------------------------
//    Qbit (pure state vector |ψ⟩)       DensityMatrix (mixed state ρ)
//    unitary evolution only             unitary + Kraus-channel noise (ρ→ΣKρK†)
//    scans the (γ,α) energy surface     scans the NOISE level p (the analogue)
//    entanglement, no decoherence       decoherence: purity/fidelity/entropy decay
//
//  Concepts exercised (all from DensityMatrix.h / NoiseModel.h):
//    * Kraus-channel noise evolution — apply_kraus() after every gate, with
//      single-qubit depolarising on 1-qubit gates and the 2-qubit depolarising
//      channel on each entangling Rzz (2-qubit gates are the noisy ones).
//    * partial_trace — reduce to one qubit and watch entanglement wash out.
//    * partial_measurement — project the register out qubit-by-qubit at readout.
//    * measurement_error — a classical readout confusion matrix on top of ρ.
//
//  Because ρ is 2ⁿ×2ⁿ and every gate costs O((2ⁿ)³), we do NOT re-scan the angle
//  grid here; we fix (γ*,α*) and sweep the physical noise strength instead. That
//  sweep is what MPI parallelises, and it is the open-system analogue of the
//  noiseless energy surface: quality vs. noise instead of energy vs. angles.
//
//  Build:
//    mpic++ -std=c++23 -O2 -I/opt/homebrew/include/eigen3 \
//        Electricity_procurement_noise.cpp -o Electricity_procurement_noise
//  Run (optional argv[1] overrides the number of noise samples):
//    mpirun -np 8 ./Electricity_procurement_noise
// =============================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <mpi.h>
#include "src/core/Gate.h"
#include "src/core/Qbit.h"
#include "src/core/DensityMatrix.h"
#include "src/core/NoiseModel.h"

using namespace std;

// ============================ Tunable constants ==============================
constexpr double PI = 3.14159265358979323846;

// Optimal QAOA angles found by Electricity_procurement.cpp's full 1000×1000
// noiseless grid search. Hard-coded here so both programs operate at exactly the
// same point (fair comparison); re-derive them there if the instance changes.
constexpr double BEST_GAMMA = 0.534605;
constexpr double BEST_ALPHA = 2.868;

// Noise sweep: physical error probability p on 1-qubit gates. 2-qubit gates are
// TWO_Q_FACTOR× noisier (as on real hardware). p is swept over [0, P_MAX].
constexpr int    NOISE_STEPS   = 120;
constexpr double P_MAX         = 0.05;      // 5% single-qubit gate error at the top
constexpr double TWO_Q_FACTOR  = 5.0;       // 2-qubit gates 5× noisier
constexpr double P_DEMO        = 0.01;      // "moderate noise" point for the demos
constexpr double READOUT_ERR   = 0.03;      // measurement_error probability (demo)
constexpr int    SAMPLES       = 100000;    // shots represented in the histograms

// ---- Problem parameters — identical to Electricity_procurement.cpp ----
constexpr double LAMBDA_COV    = 3000.0;
constexpr bool   ENABLE_BUDGET = true;
constexpr double B_MAX         = 400.0;
constexpr double LAMBDA_BUD    = 0.002;

// ======================== Problem data structures ============================
struct Placement {
    const char* type;
    int    start, duration;
    double cost;
    int    hours_begin, hours_end;
};

struct QuboHamiltonian {
    int n = 0;
    vector<double> h;
    vector<vector<double>> J;
    double E0 = 0.0;
    double scale = 1.0;
    vector<Placement> P;
    vector<int> demand;
};

// ======================= QUBO construction (identical instance) ==============
QuboHamiltonian build_problem() {
    const vector<int> demand = {2, 3, 3, 3, 2, 1};
    const int         T      = static_cast<int>(demand.size());

    struct Spec { const char* type; int start, duration; double cost; };
    const vector<Spec> catalogue = {
        {"nuclear", 0, 5,  70.0}, {"nuclear", 1, 5,  70.0},
        {"wind",    0, 3,  65.0}, {"wind",    1, 3,  65.0}, {"wind", 3, 3, 65.0},
        {"gas",     0, 3, 170.0}, {"gas",     2, 3, 170.0}, {"gas",  3, 3, 170.0},
    };

    vector<Placement> P;
    for (const Spec& s : catalogue) {
        if (s.start + s.duration > T) throw invalid_argument("band window exceeds horizon");
        P.push_back({s.type, s.start, s.duration, s.cost, s.start, s.start + s.duration});
    }
    const int n = static_cast<int>(P.size());

    auto demand_under = [&](const Placement& p) {
        double D = 0; for (int t = p.hours_begin; t < p.hours_end; ++t) D += demand[t]; return D;
    };
    auto overlap = [&](const Placement& a, const Placement& c) {
        return max(0, min(a.hours_end, c.hours_end) - max(a.hours_begin, c.hours_begin));
    };

    double konst = 0.0;
    for (int t = 0; t < T; ++t) konst += LAMBDA_COV * demand[t] * demand[t];
    if (ENABLE_BUDGET) konst += LAMBDA_BUD * B_MAX * B_MAX;

    vector<double> a(n);
    for (int i = 0; i < n; ++i) {
        double Di = demand_under(P[i]);
        a[i] = P[i].cost + LAMBDA_COV * (P[i].duration - 2.0 * Di);
        if (ENABLE_BUDGET) a[i] += LAMBDA_BUD * (P[i].cost * P[i].cost - 2.0 * B_MAX * P[i].cost);
    }

    vector<vector<double>> b(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) if (i != j) {
            double bij = 2.0 * LAMBDA_COV * overlap(P[i], P[j]);
            if (ENABLE_BUDGET) bij += 2.0 * LAMBDA_BUD * P[i].cost * P[j].cost;
            b[i][j] = bij;
        }

    QuboHamiltonian H;
    H.n = n; H.P = std::move(P); H.demand = demand;
    H.h.assign(n, 0.0);
    H.J.assign(n, vector<double>(n, 0.0));

    H.E0 = konst;
    for (int i = 0; i < n; ++i)               H.E0 += a[i] / 2.0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)       H.E0 += b[i][j] / 4.0;

    for (int i = 0; i < n; ++i) {
        double sum_b = 0.0;
        for (int j = 0; j < n; ++j) if (j != i) sum_b += b[i][j];
        H.h[i] = -a[i] / 2.0 - sum_b / 4.0;
    }
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) H.J[i][j] = b[i][j] / 4.0;

    double maxabs = 0.0;
    for (int i = 0; i < n; ++i) maxabs = max(maxabs, fabs(H.h[i]));
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) maxabs = max(maxabs, fabs(H.J[i][j]));
    H.scale = (maxabs > 0.0) ? 2.0 * maxabs : 1.0;
    for (int i = 0; i < n; ++i) H.h[i] /= H.scale;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) H.J[i][j] /= H.scale;
    H.E0 /= H.scale;
    return H;
}

// ============================ Shared helpers =================================
double basis_energy(int idx, const QuboHamiltonian& H) {
    double E = 0.0;
    for (int i = 0; i < H.n; ++i) {
        int zi = ((idx >> (H.n - 1 - i)) & 1) ? -1 : 1;
        E += H.h[i] * zi;
        for (int j = i + 1; j < H.n; ++j) {
            if (H.J[i][j] == 0.0) continue;
            int zj = ((idx >> (H.n - 1 - j)) & 1) ? -1 : 1;
            E += H.J[i][j] * zi * zj;
        }
    }
    return E;
}

struct Decoded {
    vector<int> chosen; double cost = 0.0; vector<int> coverage;
    bool feasible = true, within_budget = true;
};
Decoded decode(int idx, const QuboHamiltonian& H) {
    Decoded d;
    d.coverage.assign(H.demand.size(), 0);
    for (int i = 0; i < H.n; ++i)
        if ((idx >> (H.n - 1 - i)) & 1) {
            d.chosen.push_back(i);
            d.cost += H.P[i].cost;
            for (int t = H.P[i].hours_begin; t < H.P[i].hours_end; ++t) d.coverage[t]++;
        }
    for (size_t t = 0; t < H.demand.size(); ++t)
        if (d.coverage[t] < H.demand[t]) d.feasible = false;
    d.within_budget = !ENABLE_BUDGET || d.cost <= B_MAX;
    return d;
}

string ket(int idx, int n) {
    string s = "|";
    for (int i = n - 1; i >= 0; --i) s += ((idx >> i) & 1) ? '1' : '0';
    return s + ">";
}

// Ideal (noiseless) QAOA state vector at the operating point — the reference for
// fidelity(ρ) and the p=0 baseline. Uses the pure Qbit path (fast, exact).
Qbit ideal_state(const QuboHamiltonian& H) {
    Qbit q(H.n);
    { Gate Hg = Gate::H(); for (int i = 0; i < H.n; ++i) q.apply(Hg.expand(H.n, i)); }
    for (int i = 0; i < H.n; ++i)
        q.apply(Gate::Rz(2 * BEST_GAMMA * H.h[i]).expand(H.n, i));
    for (int i = 0; i < H.n; ++i)
        for (int j = i + 1; j < H.n; ++j)
            if (H.J[i][j] != 0.0)
                q.apply(Gate::rzz(2 * BEST_GAMMA * H.J[i][j], i, j, H.n));
    for (int i = 0; i < H.n; ++i)
        q.apply(Gate::Rx(2 * BEST_ALPHA).expand(H.n, i));
    return q;
}

// ============================ Noisy QAOA on ρ ================================
// Same depth-1 QAOA, but on a density matrix, with a Kraus channel applied after
// every gate: single-qubit depolarising on the 1-qubit gates, the 2-qubit
// depolarising channel on each entangling Rzz. p is the 1-qubit error probability.
DensityMatrix run_noisy_qaoa(const QuboHamiltonian& H, double p) {
    const int n = H.n;
    const double p1 = min(p, 0.75);                       // depolarising: p ∈ [0,0.75]
    const double p2 = min(TWO_Q_FACTOR * p, 1.0);         // depolarising_2q: p ∈ [0,1]
    DensityMatrix rho(n);                                 // starts at |0…0⟩⟨0…0|
    Gate Hg = Gate::H();

    auto noise1 = [&](int q) { if (p1 > 0) rho.apply_kraus(NoiseModel::depolarising(p1, q, n)); };
    auto noise2 = [&](int a, int b) { if (p2 > 0) rho.apply_kraus(NoiseModel::depolarising_2q(p2, a, b, n)); };

    // step 1: prepare |+⟩^⊗n
    for (int i = 0; i < n; ++i) { rho.apply(Hg.expand(n, i)); noise1(i); }
    // step 2: cost unitary (single-Z phases, then entangling ZZ couplings)
    for (int i = 0; i < n; ++i) { rho.apply(Gate::Rz(2 * BEST_GAMMA * H.h[i]).expand(n, i)); noise1(i); }
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (H.J[i][j] != 0.0) { rho.apply(Gate::rzz(2 * BEST_GAMMA * H.J[i][j], i, j, n)); noise2(i, j); }
    // step 3: mixer
    for (int i = 0; i < n; ++i) { rho.apply(Gate::Rx(2 * BEST_ALPHA).expand(n, i)); noise1(i); }
    return rho;
}

// ⟨H_C⟩ = Tr(H_C ρ) — H_C is diagonal, so it is Σ_i ρ_ii · E(i).
double energy_of(const DensityMatrix& rho, const vector<double>& E_basis) {
    vector<double> pr = rho.probabilities();
    double E = 0.0;
    for (size_t i = 0; i < pr.size(); ++i) E += pr[i] * E_basis[i];
    return E;
}

// Reduce ρ to qubit 0 by tracing out every other qubit (repeatedly dropping the
// current highest index, which leaves qubit 0). Returns that 1-qubit ρ.
DensityMatrix reduce_to_qubit0(DensityMatrix rho) {
    while (rho.get_num_qubits() > 1) rho = rho.partial_trace(rho.get_num_qubits() - 1);
    return rho;
}

void write_hist(const string& path, const vector<double>& probs, int n) {
    ofstream f(path);
    f << "state,count\n";
    for (int i = 0; i < (int)probs.size(); ++i)
        f << ket(i, n) << "," << (long long)llround(probs[i] * SAMPLES) << "\n";
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double t_start = MPI_Wtime();

    const int steps = (argc > 1) ? max(2, atoi(argv[1])) : NOISE_STEPS;

    QuboHamiltonian H = build_problem();
    const int dim = 1 << H.n;
    vector<double> E_basis(dim);
    for (int i = 0; i < dim; ++i) E_basis[i] = basis_energy(i, H);
    int classical_best = 0;
    for (int i = 1; i < dim; ++i) if (E_basis[i] < E_basis[classical_best]) classical_best = i;

    Qbit ideal = ideal_state(H);
    const vector<complex<double>> psi_ideal = ideal.get_state();

    if (rank == 0) {
        cout << "=== Noisy electricity-procurement QAOA (density matrix) ===\n";
        cout << "Qubits = " << H.n << " (rho is " << dim << "x" << dim << ")   "
             << "operating point (gamma,alpha) = (" << BEST_GAMMA << "," << BEST_ALPHA << ")\n";
        cout << "Noise: 1-qubit depolarising p in [0," << P_MAX << "], 2-qubit "
             << TWO_Q_FACTOR << "x, " << steps << " samples on " << size << " rank(s).\n";
        cout << "Classical optimum: " << ket(classical_best, H.n)
             << "  cost " << decode(classical_best, H).cost << "  (Ising E " << E_basis[classical_best] << ")\n";
    }

    // ---- Noise sweep (MPI-parallelised): one ρ evolution per noise level ----
    // Row layout: p, energy, purity, fidelity, entropy, top_cost, top_feasible
    constexpr int NCOL = 7;
    vector<double> local;
    for (int si = rank; si < steps; si += size) {
        double p = P_MAX * si / (steps - 1);
        DensityMatrix rho = run_noisy_qaoa(H, p);
        vector<double> pr = rho.probabilities();
        int top = 0; for (int i = 1; i < dim; ++i) if (pr[i] > pr[top]) top = i;
        Decoded d = decode(top, H);

        local.push_back(p);
        local.push_back(energy_of(rho, E_basis));
        local.push_back(rho.purity());
        local.push_back(rho.fidelity(psi_ideal));
        local.push_back(rho.entropy());
        local.push_back(d.cost);
        local.push_back((d.feasible && d.within_budget) ? 1.0 : 0.0);
    }

    // ---- Gather rows to rank 0 and write the sweep CSV ----
    if (rank == 0) {
        ofstream f("Results/electricity_noise_sweep.csv");
        f << "noise_p,energy,purity,fidelity,entropy,top_cost,top_feasible\n";
        auto dump = [&](const vector<double>& buf) {
            for (size_t i = 0; i + NCOL <= buf.size(); i += NCOL) {
                for (int c = 0; c < NCOL; ++c) f << buf[i + c] << (c + 1 < NCOL ? "," : "\n");
            }
        };
        dump(local);
        for (int r = 1; r < size; ++r) {
            int cnt; MPI_Recv(&cnt, 1, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            vector<double> buf(cnt);
            MPI_Recv(buf.data(), cnt, MPI_DOUBLE, r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            dump(buf);
        }
        f.close();
    } else {
        int cnt = (int)local.size();
        MPI_Send(&cnt, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(local.data(), cnt, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }

    double t_sweep_end = MPI_Wtime();

    // ---- Rank-0 demonstrations at representative noise levels ----
    if (rank == 0) {
        cout << "\n--- partial_trace: entanglement of qubit 0 with the rest, vs noise ---\n";
        cout << "  (pure global state => mixed 1-qubit marginal <=> entangled; noise makes ρ globally mixed)\n";
        for (double p : {0.0, P_DEMO, P_MAX}) {
            DensityMatrix rho = run_noisy_qaoa(H, p);
            DensityMatrix q0 = reduce_to_qubit0(rho);
            cout << "  p=" << p << ":  global purity=" << rho.purity()
                 << "  fidelity_vs_ideal=" << rho.fidelity(psi_ideal)
                 << "  |  reduced-q0 purity=" << q0.purity()
                 << "  entropy=" << q0.entropy() << " nats\n";
        }

        // Distribution comparison: p=0 vs p=P_DEMO (exact ρ-diagonal readout).
        DensityMatrix rho0    = run_noisy_qaoa(H, 0.0);
        DensityMatrix rhoNsy  = run_noisy_qaoa(H, P_DEMO);
        write_hist("Results/electricity_noiseless_hist.csv", rho0.probabilities(),   H.n);
        write_hist("Results/electricity_noisy_hist.csv",     rhoNsy.probabilities(), H.n);

        auto report_top = [&](const char* tag, const DensityMatrix& rho) {
            vector<double> pr = rho.probabilities();
            int top = 0; for (int i = 1; i < dim; ++i) if (pr[i] > pr[top]) top = i;
            Decoded d = decode(top, H);
            cout << "  " << tag << " most-likely " << ket(top, H.n) << " (p=" << pr[top]
                 << ")  cost " << d.cost << "  "
                 << (d.feasible ? "feasible" : "infeasible")
                 << (d.within_budget ? ", within budget" : ", OVER budget") << "\n";
        };
        cout << "\n--- measurement distribution: noiseless vs p=" << P_DEMO << " ---\n";
        report_top("noiseless:", rho0);
        report_top("noisy    :", rhoNsy);

        // partial_measurement: project the register out one qubit at a time,
        // then optionally distort the classical outcome with a readout confusion
        // matrix (measurement_error) — a full noisy readout of the noisy state.
        cout << "\n--- partial_measurement: qubit-by-qubit readout of the noisy state (p=" << P_DEMO << ") ---\n";
        DensityMatrix shot = run_noisy_qaoa(H, P_DEMO);
        vector<double> pre = shot.probabilities();          // for the readout-error demo
        int measured = 0;
        for (int q = 0; q < H.n; ++q) {
            int bit = shot.partial_measurement(q);           // collapses ρ on qubit q
            measured = (measured << 1) | bit;
            cout << "  measured qubit " << q << " -> " << bit << "\n";
        }
        Decoded dm = decode(measured, H);
        cout << "  collapsed outcome " << ket(measured, H.n) << "  cost " << dm.cost
             << (dm.feasible ? "  feasible" : "  infeasible")
             << (dm.within_budget ? ", within budget\n" : ", OVER budget\n");

        auto me = NoiseModel::measurement_error(READOUT_ERR, READOUT_ERR);
        int noisy_bit0 = me.sample(pre, 0, H.n);
        cout << "  measurement_error(" << READOUT_ERR << ") applied to qubit 0 readout -> reported "
             << noisy_bit0 << " (classical readout confusion on top of ρ)\n";
    }

    double t_end = MPI_Wtime();
    MPI_Finalize();
    if (rank == 0) {
        cout << "\nWrote Results/electricity_noise_sweep.csv, electricity_noiseless_hist.csv, "
                "electricity_noisy_hist.csv\n";
        cout << "Timing: noise sweep " << (t_sweep_end - t_start)
             << " s, demos " << (t_end - t_sweep_end) << " s, total " << (t_end - t_start) << " s\n";
    }
    return 0;
}
