#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <limits>
#include <mpi.h>
#include "src/core/Gate.h"
#include "src/core/Qbit.h"

using namespace std;

mt19937 gen(42); //setting seed fro reproductability

constexpr double PI = 3.14159265358979323846;
constexpr int STEPS = 1000;
constexpr double GAMMA_MAX = 2 * PI;
constexpr double ALPHA_MAX = 2 * PI;
constexpr int SAMPLES = 1000;

double compute_energy(const Qbit& q) {
    const auto& state = q.get_state();
    // Z0 = +1 for |0x⟩, -1 for |1x⟩
    // Z1 = +1 for |x0⟩, -1 for |x1⟩
    double E = 0.0;
    for (int i = 0; i < 4; ++i) {
        int b0 = (i >> 1) & 1;
        int b1 = i & 1;
        int z0 = b0 == 0 ? 1 : -1;
        int z1 = b1 == 0 ? 1 : -1;
        E += norm(state[i]) * (0.45 * z0 + 0.15 * z1);
    }
    return E;
}

Qbit apply_qaoa(double gamma, double alpha) {
    Qbit q(2);

    // Prepare |+⟩ state
    Gate H = Gate::H();
    q.apply(H.expand(2, 0));
    q.apply(H.expand(2, 1));

    // Apply cost unitary: exp(-i * gamma * H_C) = Rz(2 * gamma * 1)
    q.apply(Gate::Rz(2 * gamma * 0.45).expand(2, 0));
    q.apply(Gate::Rz(2 * gamma * 0.15).expand(2, 1));

    // Apply mixer unitary: Rx(2 * alpha)
    q.apply(Gate::Rx(2 * alpha).expand(2, 0));
    q.apply(Gate::Rx(2 * alpha).expand(2, 1));

    return q;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double t_start = MPI_Wtime();

    double local_best_energy = numeric_limits<double>::max();
    double local_best_gamma = 0.0;
    double local_best_alpha = 0.0;
    vector<tuple<double, double, double>> local_data;

    for (int gi = rank; gi < STEPS; gi += size) {
        double gamma = GAMMA_MAX * gi / (STEPS - 1);
        for (int ai = 0; ai < STEPS; ++ai) {
            double alpha = ALPHA_MAX * ai / (STEPS - 1);
            Qbit q = apply_qaoa(gamma, alpha);
            double E = compute_energy(q);
            local_data.emplace_back(gamma, alpha, E);

            if (E < local_best_energy) {
                local_best_energy = E;
                local_best_gamma = gamma;
                local_best_alpha = alpha;
            }
        }
    }

    double t_grid_end = MPI_Wtime();

    struct {
        double energy;
        int rank;
    } local_result{local_best_energy, rank}, global_result;

    MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

    double best_gamma = 0.0;
    double best_alpha = 0.0;

    if (rank == global_result.rank) {
        best_gamma = local_best_gamma;
        best_alpha = local_best_alpha;
    }

    // Broadcast best (gamma, alpha) to all ranks
    MPI_Bcast(&best_gamma, 1, MPI_DOUBLE, global_result.rank, MPI_COMM_WORLD);
    MPI_Bcast(&best_alpha, 1, MPI_DOUBLE, global_result.rank, MPI_COMM_WORLD);

    // Rank 0 writes energy surface
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
            buf[i * 3] = gamma;
            buf[i * 3 + 1] = alpha;
            buf[i * 3 + 2] = E;
        }
        MPI_Send(buf.data(), count * 3, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }

    // Parallel measurement phase
    Qbit best_q = apply_qaoa(best_gamma, best_alpha);

    // const auto& amplitudes = best_q.get_state();
    // cout << "Probabilities:\n";
    // for (int i = 0; i < (1 << best_q.num_qubits()); ++i) {
    //     double p = norm(amplitudes[i]);
    //     cout << "|" << bitset<2>(i) << ">: " << p << "\n";
    // }

    map<string, int> local_hist = {{"|00>", 0}, {"|01>", 0}, {"|10>", 0}, {"|11>", 0}};
    int local_samples = SAMPLES / size;
    for (int i = 0; i < local_samples; ++i) {
        string outcome = best_q.measure();
        local_hist[outcome]++;
    }

    // Collect counts to rank 0
    map<string, int> global_hist = {{"|00>", 0}, {"|01>", 0}, {"|10>", 0}, {"|11>", 0}};
    for (auto& [state, count] : local_hist) {
        int total;
        MPI_Reduce(&count, &total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) global_hist[state] = total;
    }

    if (rank == 0) {
        cout << "Best gamma: " << best_gamma << "\n";
        cout << "Best alpha: " << best_alpha << "\n";
        cout << "Minimum energy: " << global_result.energy << "\n";

        ofstream hist_file("Results/measurement_histogram.csv");
        hist_file << "state,count\n";
        for (const auto& [state, count] : global_hist) {
            hist_file << state << "," << count << "\n";
        }
        hist_file.close();
    }

    double t_measure_end = MPI_Wtime();

    double total_runtime = t_measure_end - t_start;
    double grid_time = t_grid_end - t_start;
    double measure_time = t_measure_end - t_grid_end;

    MPI_Finalize();

    if (rank == 0) {
        cout << "Timing summary:\n";
        cout << "  Grid search time: " << grid_time << " s\n";
        cout << "  Measurement time: " << measure_time << " s\n";
        cout << "  Total time: " << total_runtime << " s\n";
    }       
    return 0;
}