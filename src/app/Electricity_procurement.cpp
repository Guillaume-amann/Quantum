// =============================================================================
//  Electricity Procurement via Bin Packing — QUBO/QAOA solver
//
//  Implements the formulation in electricity_qubo_formulation.md on the same
//  gate-model stack as QuantumSim.cpp (Gate.h + Qbit.h, parallelised with MPI).
//
//  What is "stepped up" relative to QuantumSim.cpp:
//    * |P| qubits instead of 2 — one per valid band placement (Section 2/5).
//    * A real QUBO built programmatically in build_problem() (Sections 4–9),
//      not a hand-coded 2-term cost.
//    * Genuine entanglement: the coverage term produces ZᵢZⱼ couplings for every
//      pair of placements that overlap in time, applied as native Rzz gates
//      (Section 10). The state is no longer a product state.
//    * Brute-force classical_optimum() for validation (checklist step 7).
//    * Non-destructive sampling of the optimised state (QuantumSim.cpp re-used a
//      collapsed Qbit, which returns the same outcome every shot — fixed here).
//    * Decode + feasibility check of the sampled solution (checklist step 9).
//
//  The MPI grid-search / CSV-output / timing skeleton is intentionally identical
//  to QuantumSim.cpp so the existing Results/*.py plotting scripts still apply.
//
//  Build:
//    mpic++ -std=c++23 -O2 -I/opt/homebrew/include/eigen3 \
//        Electricity_procurement.cpp -o Electricity_procurement
//  Run (optional argv[1] overrides the grid resolution, default STEPS):
//    mpirun -np 4 ./Electricity_procurement
// =============================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <cmath>
#include <limits>
#include <string>
#include <algorithm>
#include <random>
#include <mpi.h>
#include "src/core/Gate.h"
#include "src/core/Qbit.h"
#include "src/core/DensityMatrix.h"
#include "src/core/NoiseModel.h"

using namespace std;

// ============================ Tunable constants ==============================
constexpr double PI        = 3.14159265358979323846;
constexpr int    STEPS     = 1000;          // (γ,α) grid resolution per axis
constexpr double GAMMA_MAX = 2 * PI;
constexpr double ALPHA_MAX = 2 * PI;
constexpr int    SAMPLES   = 100000;        // measurement shots for the histogram

// Penalty weights (Section 3). λ_cov must dominate the max achievable total cost
// (here Σ all bands = €845) so coverage behaves as a hard constraint; λ_bud is the
// soft budget weight. With the budget on (C4 / Section 7.3) the coupling graph
// becomes dense — every pair of bands is coupled through total spend, not just the
// time-overlapping ones — which is the extra entanglement asked for.
//
// λ_bud is kept small on purpose: the .md budget term is the symmetric penalty
// λ_bud·(Σcᵢxᵢ − B_max)², which targets spend ≈ B_max rather than "≤ B_max". Choosing
// λ_bud so its parabola's vertex (B_max − 1/2λ_bud) sits below the cheapest feasible
// spend keeps the optimum at the cheapest feasible plan while still discouraging
// overspend and supplying the dense ZᵢZⱼ coupling.
constexpr double LAMBDA_COV    = 3000.0;    // ≫ €845  → coverage effectively hard
constexpr bool   ENABLE_BUDGET = true;      // C4 on → dense all-pairs ZᵢZⱼ coupling
constexpr double B_MAX         = 400.0;     // budget cap (€)        (if ENABLE_BUDGET)
constexpr double LAMBDA_BUD    = 0.002;     // budget penalty weight (if ENABLE_BUDGET)

// ======================== Problem data structures ============================
// One band placement i: a generation contract of a given type, started at hour
// `start`, delivering 1 MW over the contiguous half-open range [begin, end).
// In this instance the catalogue is given explicitly (a fixed list of 8 bands)
// rather than auto-enumerated, so each Placement is one qubit. (Section 2/5)
struct Placement {
    const char* type;            // "nuclear" / "wind" / "gas"
    int    start, duration;      // start hour s, duration δ_i
    double cost;                 // c_i — total cost of one band (€)
    int    hours_begin, hours_end;   // covers [begin, end)
};

// Closed-form Ising Hamiltonian H_C = E0 + Σ hᵢ Zᵢ + Σ_{i<j} Jᵢⱼ ZᵢZⱼ  (Section 9),
// plus the metadata needed to decode a bitstring back to a procurement plan.
// h and J are stored already normalised (divided by `scale`) so the QAOA rotation
// angles 2γ·coef stay O(1) over γ∈[0,2π]; a positive rescaling leaves argmin — and
// therefore the optimal plan — unchanged. Decoded costs below use the real €.
struct QuboHamiltonian {
    int n = 0;                          // qubits = |P|
    vector<double> h;                   // hᵢ (normalised)
    vector<vector<double>> J;           // Jᵢⱼ (normalised, upper-triangular i<j)
    double E0 = 0.0;                    // constant offset (dropped in the QAOA energy)
    double scale = 1.0;                 // normalisation divisor applied to h, J, E0
    vector<Placement> P;
    vector<int> demand;
};

// ======================= QUBO construction (Sections 2–9) ====================
QuboHamiltonian build_problem() {
    // ---- Instance: a realistic day-part procurement over T = 6 hours.
    //      Demand peaks mid-horizon. Three generation types, cheap→dear per MWh:
    //      nuclear (€70 / 5 MWh = €14/MWh, baseload), wind (€65 / 3 MWh ≈ €21.7/MWh),
    //      gas (€170 / 3 MWh ≈ €56.7/MWh, peaker). Eight bands (qubits) total: the
    //      two valid nuclear windows plus three wind and three gas at spread starts.
    const vector<int> demand = {2, 3, 3, 3, 2, 1};       // d_t, MW per hour
    const int         T      = static_cast<int>(demand.size());

    struct Spec { const char* type; int start, duration; double cost; };
    const vector<Spec> catalogue = {
        {"nuclear", 0, 5,  70.0},   // q0: covers [0,5)
        {"nuclear", 1, 5,  70.0},   // q1: covers [1,6)
        {"wind",    0, 3,  65.0},   // q2: covers [0,3)
        {"wind",    1, 3,  65.0},   // q3: covers [1,4)
        {"wind",    3, 3,  65.0},   // q4: covers [3,6)
        {"gas",     0, 3, 170.0},   // q5: covers [0,3)
        {"gas",     2, 3, 170.0},   // q6: covers [2,5)
        {"gas",     3, 3, 170.0},   // q7: covers [3,6)
    };

    vector<Placement> P;
    for (const Spec& s : catalogue) {
        if (s.start + s.duration > T)                    // C3: window must fit horizon
            throw invalid_argument("build_problem(): band window exceeds horizon");
        P.push_back({s.type, s.start, s.duration, s.cost, s.start, s.start + s.duration});
    }
    const int n = static_cast<int>(P.size());

    // ---- Derived quantities (Section 4) ----
    auto demand_under = [&](const Placement& p) {            // Dᵢ = Σ_{t∈Tᵢ} d_t
        double D = 0; for (int t = p.hours_begin; t < p.hours_end; ++t) D += demand[t]; return D;
    };
    auto overlap = [&](const Placement& a, const Placement& c) {  // |Tᵢ ∩ Tⱼ|, closed form
        return max(0, min(a.hours_end, c.hours_end) - max(a.hours_begin, c.hours_begin));
    };

    // ---- QUBO canonical coefficients a_i, b_ij, const (Section 8) ----
    double konst = 0.0;
    for (int t = 0; t < T; ++t) konst += LAMBDA_COV * demand[t] * demand[t];
    if (ENABLE_BUDGET) konst += LAMBDA_BUD * B_MAX * B_MAX;

    vector<double> a(n);
    for (int i = 0; i < n; ++i) {
        double Di = demand_under(P[i]);
        a[i] = P[i].cost + LAMBDA_COV * (P[i].duration - 2.0 * Di);
        if (ENABLE_BUDGET) a[i] += LAMBDA_BUD * (P[i].cost * P[i].cost - 2.0 * B_MAX * P[i].cost);
    }

    vector<vector<double>> b(n, vector<double>(n, 0.0));     // symmetric, i≠j
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) if (i != j) {
            double bij = 2.0 * LAMBDA_COV * overlap(P[i], P[j]);
            if (ENABLE_BUDGET) bij += 2.0 * LAMBDA_BUD * P[i].cost * P[j].cost;
            b[i][j] = bij;
        }

    // ---- Ising mapping  hᵢ = −aᵢ/2 − Σ_{j≠i} bᵢⱼ/4,  Jᵢⱼ = bᵢⱼ/4  (Section 9) ----
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
        for (int j = i + 1; j < n; ++j)
            H.J[i][j] = b[i][j] / 4.0;

    // ---- Normalise so max|coefficient| = 0.5 (scale = 2·max|coef|). Keeps the
    //      QAOA angles 2γ·coef ∈ O(1) over γ∈[0,2π]; argmin (the optimal plan) is
    //      invariant under any positive rescaling. ----
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

// =============================== Energy ======================================
// Ising energy of a computational basis state (E0 dropped — it shifts every
// state equally and never changes the optimum, Section 9). zᵢ = +1 for bit 0,
// −1 for bit 1; big-endian, so qubit i is bit (n-1-i) (matches Gate.h/Qbit.h).
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

// ⟨H_C⟩ = Σ_b |ψ_b|² · E(b) — H_C is diagonal, so this is an exact expectation.
double compute_energy(const Qbit& q, const vector<double>& E_basis) {
    const auto& s = q.get_state();
    double E = 0.0;
    for (size_t i = 0; i < s.size(); ++i) E += norm(s[i]) * E_basis[i];
    return E;
}

// =============================== QAOA circuit (Section 10) ====================
// Depth-1: a precomputed |+⟩^⊗n, then the cost unitary (Rz single-Z phases +
// Rzz entangling couplings), then the transverse-field mixer (Rx on every qubit).
Qbit apply_qaoa(double gamma, double alpha, const QuboHamiltonian& H, const Qbit& plus) {
    Qbit q = plus;                                              // step 1: |+⟩^⊗n

    // step 2: U_C(γ) = Πᵢ exp(−iγhᵢZᵢ) · Π_{i<j} exp(−iγJᵢⱼZᵢZⱼ)
    for (int i = 0; i < H.n; ++i)
        if (H.h[i] != 0.0)
            q.apply(Gate::Rz(2 * gamma * H.h[i]).expand(H.n, i));
    for (int i = 0; i < H.n; ++i)
        for (int j = i + 1; j < H.n; ++j)
            if (H.J[i][j] != 0.0)
                q.apply(Gate::rzz(2 * gamma * H.J[i][j], i, j, H.n));  // the entangler

    // step 3: U_M(α) = Πᵢ Rx(2α)
    for (int i = 0; i < H.n; ++i)
        q.apply(Gate::Rx(2 * alpha).expand(H.n, i));

    return q;
}

// =========================== Decode & feasibility (step 9) ===================
struct Decoded {
    vector<int> chosen;        // selected placement indices
    double      cost = 0.0;
    vector<int> coverage;      // cov(t)
    bool        feasible = true;     // cov(t) ≥ d_t  for every hour
    bool        within_budget = true;
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

void print_solution(const char* tag, int idx, const QuboHamiltonian& H, const vector<double>& E_basis) {
    Decoded d = decode(idx, H);
    cout << "  " << tag << " " << ket(idx, H.n)
         << "  Ising E = " << E_basis[idx]
         << "  cost = "    << d.cost
         << "  coverage = [";
    for (size_t t = 0; t < d.coverage.size(); ++t) cout << d.coverage[t] << (t + 1 < d.coverage.size() ? "," : "");
    cout << "]  vs demand [";
    for (size_t t = 0; t < H.demand.size(); ++t)  cout << H.demand[t]   << (t + 1 < H.demand.size()  ? "," : "");
    cout << "]  -> " << (d.feasible ? "FEASIBLE" : "infeasible");
    if (ENABLE_BUDGET) cout << ", " << (d.within_budget ? "within budget" : "OVER budget");
    cout << "\n";
}

// ============================ Reporting (rank 0) =============================
void print_problem(const QuboHamiltonian& H) {
    cout << "=== Electricity procurement QUBO ===\n";
    cout << "Qubits |P| = " << H.n << "   (Hilbert space dim = " << (1 << H.n) << ")\n";
    cout << "Demand d_t = [";
    for (size_t t = 0; t < H.demand.size(); ++t) cout << H.demand[t] << (t + 1 < H.demand.size() ? "," : "");
    cout << "]   lambda_cov = " << LAMBDA_COV;
    if (ENABLE_BUDGET) cout << "   [budget C4 ON: B_max = " << B_MAX << ", lambda_bud = " << LAMBDA_BUD << "]";
    else               cout << "   [budget C4 off]";
    cout << "\n";
    cout << "Bands (1 MW each, cost in EUR per band):\n";
    for (int i = 0; i < H.n; ++i)
        cout << "  q" << i << " = " << H.P[i].type << " @ h" << H.P[i].start
             << "  covers [" << H.P[i].hours_begin << "," << H.P[i].hours_end << ")"
             << "  cost " << H.P[i].cost << "\n";
    cout << "Ising H_C normalised by scale = " << H.scale << " (QAOA angle scaling; optimum unchanged)\n";
    cout << "  h_i: ";
    for (int i = 0; i < H.n; ++i) cout << "h" << i << "=" << H.h[i] << "  ";
    cout << "\n  J_ij" << (ENABLE_BUDGET ? " (dense — every pair coupled via total spend):\n" : " (nonzero = overlapping bands):\n");
    for (int i = 0; i < H.n; ++i)
        for (int j = i + 1; j < H.n; ++j)
            if (H.J[i][j] != 0.0)
                cout << "    J(" << i << "," << j << ") = " << H.J[i][j] << "\n";
    cout << "====================================\n";
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double t_start = MPI_Wtime();

    const int steps = (argc > 1) ? max(2, atoi(argv[1])) : STEPS;  // optional override

    // ---- Build the Hamiltonian and precompute angle-independent data ----
    QuboHamiltonian H = build_problem();
    const int dim = 1 << H.n;

    vector<double> E_basis(dim);
    for (int idx = 0; idx < dim; ++idx) E_basis[idx] = basis_energy(idx, H);

    Qbit plus(H.n);                                   // |+⟩^⊗n, built once, reused every grid point
    { Gate Hg = Gate::H(); for (int i = 0; i < H.n; ++i) plus.apply(Hg.expand(H.n, i)); }

    // Classical brute-force optimum (validation, checklist step 7)
    int classical_best = 0;
    for (int idx = 1; idx < dim; ++idx) if (E_basis[idx] < E_basis[classical_best]) classical_best = idx;

    if (rank == 0) {
        print_problem(H);
        cout << "Classical brute-force optimum:\n";
        print_solution("argmin =", classical_best, H, E_basis);
        cout << "Running QAOA grid search on " << size << " rank(s), "
             << steps << "x" << steps << " (gamma,alpha) grid...\n";
    }

    // ---- QAOA grid search over (γ, α): each rank scans a γ-stripe ----
    double local_best_energy = numeric_limits<double>::max();
    double local_best_gamma = 0.0, local_best_alpha = 0.0;
    vector<tuple<double, double, double>> local_data;

    for (int gi = rank; gi < steps; gi += size) {
        double gamma = GAMMA_MAX * gi / (steps - 1);
        for (int ai = 0; ai < steps; ++ai) {
            double alpha = ALPHA_MAX * ai / (steps - 1);
            Qbit q = apply_qaoa(gamma, alpha, H, plus);
            double E = compute_energy(q, E_basis);
            local_data.emplace_back(gamma, alpha, E);

            if (E < local_best_energy) {
                local_best_energy = E;
                local_best_gamma  = gamma;
                local_best_alpha  = alpha;
            }
        }
    }

    double t_grid_end = MPI_Wtime();

    // ---- Best (γ,α) via MINLOC reduction, then broadcast it ----
    struct { double energy; int rank; } local_result{local_best_energy, rank}, global_result;
    MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

    double best_gamma = 0.0, best_alpha = 0.0;
    if (rank == global_result.rank) { best_gamma = local_best_gamma; best_alpha = local_best_alpha; }
    MPI_Bcast(&best_gamma, 1, MPI_DOUBLE, global_result.rank, MPI_COMM_WORLD);
    MPI_Bcast(&best_alpha, 1, MPI_DOUBLE, global_result.rank, MPI_COMM_WORLD);

    // ---- Rank 0 gathers the full surface and writes energy_surface.csv ----
    if (rank == 0) {
        ofstream energy_file("Results/energy_surface.csv");
        energy_file << "gamma,alpha,energy\n";
        for (int r = 0; r < size; ++r) {
            if (r == 0) {
                for (const auto& [gamma, alpha, E] : local_data)
                    energy_file << gamma << "," << alpha << "," << E << "\n";
            } else {
                int count;
                MPI_Recv(&count, 1, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                vector<double> buf(count * 3);
                MPI_Recv(buf.data(), count * 3, MPI_DOUBLE, r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int i = 0; i < count; ++i)
                    energy_file << buf[i * 3] << "," << buf[i * 3 + 1] << "," << buf[i * 3 + 2] << "\n";
            }
        }
        energy_file.close();
    } else {
        int count = local_data.size();
        MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        vector<double> buf(count * 3);
        for (int i = 0; i < count; ++i) {
            auto [gamma, alpha, E] = local_data[i];
            buf[i * 3] = gamma; buf[i * 3 + 1] = alpha; buf[i * 3 + 2] = E;
        }
        MPI_Send(buf.data(), count * 3, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }

    // ---- Measurement: non-destructive sampling of the optimised state ----
    // Build the best state once, form its CDF, and draw shots without collapsing
    // (so every shot is an independent draw from |ψ|² — the actual distribution).
    Qbit best_q = apply_qaoa(best_gamma, best_alpha, H, plus);
    const auto& amp = best_q.get_state();
    vector<double> cdf(dim);
    cdf[0] = norm(amp[0]);
    for (int i = 1; i < dim; ++i) cdf[i] = cdf[i - 1] + norm(amp[i]);
    double total = cdf[dim - 1];

    mt19937 rng(1234 + rank);                          // per-rank seed: reproducible, decorrelated
    uniform_real_distribution<double> dist(0.0, 1.0);
    int local_samples = SAMPLES / size + (rank < SAMPLES % size ? 1 : 0);
    vector<int> local_counts(dim, 0), global_counts(dim, 0);
    for (int s = 0; s < local_samples; ++s) {
        double r = dist(rng) * total;
        int idx = static_cast<int>(lower_bound(cdf.begin(), cdf.end(), r) - cdf.begin());
        if (idx >= dim) idx = dim - 1;                 // clamp against FP rounding
        local_counts[idx]++;
    }
    MPI_Reduce(local_counts.data(), global_counts.data(), dim, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "\nBest gamma: " << best_gamma << "\n";
        cout << "Best alpha: " << best_alpha << "\n";
        cout << "Minimum <H_C> (normalised units): " << global_result.energy << "\n";

        // (a) Distribution mode — the single most-probable basis state.
        int mode = 0;
        for (int idx = 1; idx < dim; ++idx) if (global_counts[idx] > global_counts[mode]) mode = idx;
        cout << "QAOA distribution mode (" << global_counts[mode] << "/" << SAMPLES << " shots):\n";
        print_solution("mode        =", mode, H, E_basis);

        // (b) Best feasible, within-budget shot = QAOA's recommended procurement.
        //     This is how QAOA is used for optimisation: sample, post-select feasible,
        //     keep the cheapest. Ties broken by how often it was sampled.
        int best_shot = -1;
        for (int idx = 0; idx < dim; ++idx) {
            if (global_counts[idx] == 0) continue;
            Decoded d = decode(idx, H);
            if (!d.feasible || !d.within_budget) continue;
            if (best_shot < 0 || decode(idx, H).cost < decode(best_shot, H).cost) best_shot = idx;
        }
        if (best_shot < 0) {
            cout << "QAOA recommended procurement: none of the sampled states was feasible+within-budget.\n";
        } else {
            cout << "QAOA recommended procurement (cheapest feasible shot, "
                 << global_counts[best_shot] << "/" << SAMPLES << " shots):\n";
            print_solution("best shot   =", best_shot, H, E_basis);
            cout << (best_shot == classical_best
                         ? "  -> matches the classical brute-force optimum.\n"
                         : "  -> feasible & within budget, above the classical optimum cost (depth-1 limit).\n");
        }

        ofstream hist_file("Results/measurement_histogram.csv");
        hist_file << "state,count\n";
        for (int idx = 0; idx < dim; ++idx)
            hist_file << ket(idx, H.n) << "," << global_counts[idx] << "\n";
        hist_file.close();
    }

    double t_measure_end = MPI_Wtime();
    double total_runtime = t_measure_end - t_start;
    double grid_time     = t_grid_end - t_start;
    double measure_time  = t_measure_end - t_grid_end;

    MPI_Finalize();

    if (rank == 0) {
        cout << "Timing summary:\n";
        cout << "  Grid search time: " << grid_time   << " s\n";
        cout << "  Measurement time: " << measure_time << " s\n";
        cout << "  Total time: "       << total_runtime << " s\n";
    }
    return 0;
}
