// =============================================================================
//  QuantumSVM.cpp — QSVM anomaly detection for H29 Petite Enfance
//
//  Purpose:
//    Train QSVM on historical H29 crèche cost data (24 months, 50 crèches).
//    Compute quantum kernel matrix K_Q.
//    Score monthly data for anomalies (deviations from normal profile).
//    Output anomaly flags, confidence scores, and decision boundary distances.
//
//  Workflow:
//    1. Load H29 historical data (11 features per crèche-month)
//    2. Compute quantum kernel matrix K_Q (all pairs)
//    3. Export K_Q for classical SVM solver (Python/scikit-learn)
//    4. Load monthly data, compute kernel row, score anomalies
//    5. Output CSV and summary report
//
//  Data format (11 features):
//    [presences_planned, presences_recorded, presences_billed,
//     cost_direct, cost_locaux, cost_personnel, cost_energie,
//     tarif_moyen, ecart_prev_vs_fact, age_mix, absenteeisme]
//
//  Build:
//    g++ -std=c++17 -O2 -I/opt/homebrew/include/eigen3 \
//        -I/path/to/Quantum/src/core \
//        QuantumSVM.cpp -o QuantumSVM
//
//  Run:
//    ./QuantumSVM  [optional: n_qubits (default 8)]
//
// =============================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <random>
#include <iomanip>
#include "src/core/Gate.h"
#include "src/core/DensityMatrix.h"
#include "src/core/QKernel.h"

using namespace std;
using namespace Eigen;

constexpr int N_FEATURES = 11;
constexpr int N_TRAIN_creche = 50;
constexpr int N_TRAIN_MONTHS = 24;
constexpr int N_TEST_MONTHS = 1;  // score current month

// ============================================================================
//  Synthetic data generation (for PoC: replace with real H29 load)
// ============================================================================

struct H29Record {
    string creche_id;
    int month;
    vector<double> features;  // 11 features
};

// Generate synthetic "normal" data: 50 crèches × 24 months = 1200 points
vector<H29Record> generate_training_data(int n_creche = 50, int n_months = 24) {
    vector<H29Record> data;
    static mt19937 gen(42);  // seed for reproducibility
    
    normal_distribution<double> presence_dist(10.0, 2.0);     // avg 10 kids/day
    normal_distribution<double> cost_dist(100.0, 20.0);       // costs
    normal_distribution<double> tarif_dist(35.0, 5.0);        // tarif moyen
    
    for (int c = 0; c < n_creche; ++c) {
        string crèche_name = "Crèche_" + to_string(c);
        
        for (int m = 0; m < n_months; ++m) {
            H29Record rec;
            rec.creche_id
     = crèche_name;
            rec.month = m;
            
            // Synthetic features (11 dims)
            double presences_planned = max(1.0, presence_dist(gen));
            double presences_recorded = presences_planned * (0.95 + 0.10 * (gen() % 100) / 100.0);
            double presences_billed = presences_recorded * (0.98 + 0.04 * (gen() % 100) / 100.0);
            
            rec.features.push_back(presences_planned);
            rec.features.push_back(presences_recorded);
            rec.features.push_back(presences_billed);
            rec.features.push_back(max(10.0, cost_dist(gen)));          // cost_direct
            rec.features.push_back(max(10.0, cost_dist(gen)));          // cost_locaux
            rec.features.push_back(max(10.0, cost_dist(gen)));          // cost_personnel
            rec.features.push_back(max(5.0, cost_dist(gen) * 0.5));     // cost_energie
            rec.features.push_back(max(20.0, tarif_dist(gen)));         // tarif_moyen
            rec.features.push_back(fabs(presences_planned - presences_billed) / max(1.0, presences_planned));  // écart
            rec.features.push_back(0.5 + 0.2 * (gen() % 100) / 100.0);  // age_mix (0.5 to 0.7)
            rec.features.push_back(0.05 + 0.05 * (gen() % 100) / 100.0); // absentéisme
            
            data.push_back(rec);
        }
    }
    
    return data;
}

// Generate test data with some anomalies injected
vector<H29Record> generate_test_data(int n_creche = 50, double anomaly_rate = 0.1) {
    vector<H29Record> data = generate_training_data(n_creche, 1);  // 1 month
    
    // Inject anomalies in ~10% of crèches
    static mt19937 gen(123);
    uniform_real_distribution<double> anomaly_dist(0.0, 1.0);
    
    for (auto& rec : data) {
        if (anomaly_dist(gen) < anomaly_rate) {
            // Inject anomaly: spike in tarif_moyen and cost_energie
            rec.features[7] *= 1.5;  // tarif_moyen spike
            rec.features[6] *= 2.0;  // cost_energie spike
        }
    }
    
    return data;
}

// ============================================================================
//  Normalize data (z-score)
// ============================================================================

pair<vector<vector<double>>, pair<VectorXd, VectorXd>> 
normalize_data(const vector<H29Record>& data) {
    vector<vector<double>> X;
    for (const auto& rec : data) X.push_back(rec.features);
    
    VectorXd mean = VectorXd::Zero(N_FEATURES);
    VectorXd std_dev = VectorXd::Zero(N_FEATURES);
    
    // Compute mean
    for (const auto& x : X) {
        for (int f = 0; f < N_FEATURES; ++f) {
            mean(f) += x[f];
        }
    }
    mean /= X.size();
    
    // Compute std dev
    for (const auto& x : X) {
        for (int f = 0; f < N_FEATURES; ++f) {
            std_dev(f) += (x[f] - mean(f)) * (x[f] - mean(f));
        }
    }
    std_dev = std_dev.array() / X.size();
    std_dev = std_dev.array().sqrt();
    
    // Normalize
    vector<vector<double>> X_norm;
    for (const auto& x : X) {
        vector<double> x_norm;
        for (int f = 0; f < N_FEATURES; ++f) {
            double z = (x[f] - mean(f)) / max(1e-9, std_dev(f));
            x_norm.push_back(z);
        }
        X_norm.push_back(x_norm);
    }
    
    return {X_norm, {mean, std_dev}};
}

// ============================================================================
//  Main QSVM workflow
// ============================================================================

int main(int argc, char** argv) {
    int n_qubits = (argc > 1) ? atoi(argv[1]) : 8;
    
    if (n_qubits < 2 || n_qubits > 15) {
        cerr << "Usage: " << argv[0] << " [n_qubits (2-15, default 8)]\n";
        return 1;
    }
    
    cout << "=== Quantum SVM for H29 anomaly detection ===\n";
    cout << "Qubits: " << n_qubits << "\n";
    cout << "Features: " << N_FEATURES << "\n";
    cout << "Training set: " << N_TRAIN_creche << " crèches × " << N_TRAIN_MONTHS << " months = " 
         << (N_TRAIN_creche * N_TRAIN_MONTHS) << " points\n\n";
    
    // ========== STEP 1: Load and normalize training data ==========
    cout << "Step 1: Load training data...\n";
    auto train_data = generate_training_data(N_TRAIN_creche, N_TRAIN_MONTHS);
    auto [X_train_norm, norm_stats] = normalize_data(train_data);
    auto [mean, std_dev] = norm_stats;
    
    cout << "  ✓ Loaded " << train_data.size() << " training records\n";
    cout << "  ✓ Normalized (z-score)\n";
    
    // ========== STEP 2: Create quantum kernel and compute K_Q matrix ==========
    cout << "\nStep 2: Initialize quantum kernel and compute K_Q matrix...\n";
    cout << "  (This may take a few moments...)\n";
    
    QKernel qk(n_qubits, N_FEATURES, 2);  // 2-layer circuit
    qk.print_circuit_info();
    
    auto start = chrono::high_resolution_clock::now();
    auto K_Q = qk.compute_kernel_matrix(X_train_norm);
    auto end = chrono::high_resolution_clock::now();
    
    double kernel_time = chrono::duration<double>(end - start).count();
    cout << "  ✓ K_Q computed in " << kernel_time << " s\n";
    cout << "  Kernel matrix: " << K_Q.rows() << "×" << K_Q.cols() << "\n";
    
    // ========== STEP 3: Export K_Q for classical SVM solver ==========
    cout << "\nStep 3: Export kernel matrix for classical SVM...\n";
    
    ofstream kernel_file("Results/quantum_kernel_matrix.csv");
    kernel_file << "# Quantum kernel matrix K_Q for H29 training set\n";
    kernel_file << "# Size: " << K_Q.rows() << "x" << K_Q.cols() << "\n";
    kernel_file << "# Format: CSV (space or comma separated)\n";
    for (int i = 0; i < K_Q.rows(); ++i) {
        for (int j = 0; j < K_Q.cols(); ++j) {
            kernel_file << K_Q(i, j);
            if (j < K_Q.cols() - 1) kernel_file << ",";
        }
        kernel_file << "\n";
    }
    kernel_file.close();
    
    cout << "  ✓ K_Q exported to Results/quantum_kernel_matrix.csv\n";
    cout << "  → Next: Use this with sklearn.svm.SVC(kernel='precomputed') in Python\n";
    
    // ========== STEP 4: Load and score test data ==========
    cout << "\nStep 4: Load monthly test data and score anomalies...\n";
    
    auto test_data = generate_test_data(N_TRAIN_creche, 0.1);  // 10% anomalies
    auto [X_test_norm, _] = normalize_data(test_data);
    
    cout << "  ✓ Loaded " << test_data.size() << " test records\n";
    
    // ========== STEP 5: Compute anomaly scores (distance to hyperplane) ==========
    // This is a placeholder: real scores come from SVM solver
    // For now, compute distance to mean in kernel space as proxy
    
    cout << "\nStep 5: Compute anomaly scores...\n";
    
    ofstream scores_file("Results/anomaly_scores.csv");
    scores_file << "creche_id,month,kernel_distance,anomaly_flag,confidence\n";
    
    vector<pair<string, double>> anomaly_list;
    
    for (size_t i = 0; i < X_test_norm.size(); ++i) {
        // Compute kernel row (distance to all training points)
        VectorXd k_row = qk.kernel_row(X_test_norm[i], X_train_norm);
        
        // Proxy score: average kernel value (similarity to training set)
        double avg_kernel = k_row.mean();
        double std_kernel = sqrt((k_row.array() - avg_kernel).square().mean());
        
        // Z-score: how far from typical similarity?
        double z_score = (avg_kernel - 0.5) / max(0.1, std_kernel);  // 0.5 is nominal
        
        bool is_anomaly = abs(z_score) > 1.5;  // threshold
        double confidence = min(1.0, abs(z_score) / 3.0);
        
        string crèche = test_data[i].creche_id
;
        scores_file << crèche << "," << test_data[i].month << ","
                    << avg_kernel << "," << (is_anomaly ? "YES" : "NO") << ","
                    << confidence << "\n";
        
        if (is_anomaly) {
            anomaly_list.push_back({crèche, confidence});
        }
    }
    scores_file.close();
    
    cout << "  ✓ Anomaly scores computed\n";
    cout << "  → Exported to Results/anomaly_scores.csv\n";
    
    // ========== STEP 6: Summary report ==========
    cout << "\n" << string(60, '=') << "\n";
    cout << "SUMMARY REPORT\n";
    cout << string(60, '=') << "\n";
    
    cout << "\nTraining set: " << train_data.size() << " records (normal data)\n";
    cout << "Test set: " << test_data.size() << " records (July 2026)\n";
    cout << "Quantum kernel computation time: " << kernel_time << " s\n";
    cout << "Hilbert space dimension: 2^" << n_qubits << " = " << (1 << n_qubits) << "\n";
    
    cout << "\nAnomalies detected (threshold: z-score > 1.5):\n";
    if (anomaly_list.empty()) {
        cout << "  None (all crèches within normal range)\n";
    } else {
        cout << "  Total: " << anomaly_list.size() << " crèches\n";
        for (const auto& [crèche, conf] : anomaly_list) {
            cout << "    • " << crèche << " (confidence: " << fixed << setprecision(2) 
                 << (conf * 100) << "%)\n";
        }
    }
    
    cout << "\nOutput files:\n";
    cout << "  • Results/quantum_kernel_matrix.csv (K_Q for SVM solver)\n";
    cout << "  • Results/anomaly_scores.csv (per-crèche scores)\n";
    
    cout << "\nNext steps:\n";
    cout << "  1. Train classical SVM on K_Q with sklearn\n";
    cout << "  2. Fine-tune anomaly threshold based on known ground truth\n";
    cout << "  3. Deploy monthly monitoring pipeline\n";
    cout << "  4. Compare QSVM vs classical RBF kernel performance\n";
    
    cout << "\n" << string(60, '=') << "\n";
    cout << "✓ QSVM anomaly detection complete\n";
    cout << string(60, '=') << "\n\n";
    
    return 0;
}