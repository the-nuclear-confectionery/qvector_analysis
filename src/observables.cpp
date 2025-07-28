#include "observables.h"
#include <numeric>
#include <cmath>

// --- Registration ---
void Observables::register_integrated(const std::string& name, IntegratedObservable fn) {
    integrated_observables[name] = fn;
}

void Observables::register_differential(const std::string& name, DifferentialObservable fn) {
    differential_observables[name] = fn;
}

void Observables::evaluate_all_integrated(const std::vector<event>& events, cuts& cut,
                                          const std::vector<int>& pids, int n_max) {
    // Temporary storage
    std::unordered_map<std::string, std::unordered_map<int, std::vector<double>>> value_acc;
    std::unordered_map<std::string, std::unordered_map<int, std::vector<double>>> weight_acc;

    for (const event& ev : events) {
        for (const auto& [name, fn] : integrated_observables) {
            for (int pid : pids) {
                for (int n = 1; n <= n_max; ++n) {
                    auto [val, w] = fn(const_cast<event&>(ev), cut, pid, n);
                    value_acc[name][pid].push_back(val);
                    weight_acc[name][pid].push_back(w);
                }
            }
        }
    }

    // Finalize
    for (const auto& [name, pid_map] : value_acc) {
        for (const auto& [pid, vals] : pid_map) {
            const auto& wgts = weight_acc[name][pid];
            std::string key = name + "_pid" + std::to_string(pid);
            scalar_values[key] = weighted_average(vals, wgts);
            scalar_errors[key] = compute_scalar_jackknife_error(vals, wgts);
        }
    }
}

void Observables::evaluate_all_differential(const std::vector<event>& events, cuts& cut,
                                            const std::vector<int>& pids, int n_max) {
    // Temporary storage
    std::unordered_map<std::string, std::unordered_map<int, std::vector<std::vector<double>>>> val_acc;
    std::unordered_map<std::string, std::unordered_map<int, std::vector<std::vector<double>>>> wgt_acc;

    for (const event& ev : events) {
        for (const auto& [name, fn] : differential_observables) {
            for (int pid : pids) {
                for (int n = 1; n <= n_max; ++n) {
                    auto vec = fn(const_cast<event&>(ev), cut, pid, n);
                    std::vector<double> vals, wgts;
                    for (const auto& [v, w] : vec) {
                        vals.push_back(v);
                        wgts.push_back(w);
                    }
                    val_acc[name][pid].push_back(vals);
                    wgt_acc[name][pid].push_back(wgts);
                }
            }
        }
    }

    // Finalize
    for (const auto& [name, pid_map] : val_acc) {
        for (const auto& [pid, vals] : pid_map) {
            const auto& wgts = wgt_acc[name][pid];
            std::string key = name + "_pid" + std::to_string(pid);

            // Collapse weight vectors into per-bin sums
            std::vector<double> weight_sums(wgts[0].size(), 0.0);
            for (const auto& w : wgts)
                for (size_t i = 0; i < w.size(); ++i)
                    weight_sums[i] += w[i];

            vector_values[key] = weighted_average(vals, weight_sums);
            vector_errors[key] = compute_vector_jackknife_error(vals, weight_sums);
        }
    }
}


// --- Weighted Averages ---
double Observables::weighted_average(const std::vector<double>& values, const std::vector<double>& weights) {
    double wsum = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (wsum == 0.0) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < values.size(); ++i)
        sum += values[i] * weights[i];
    return sum / wsum;
}

std::vector<double> Observables::weighted_average(const std::vector<std::vector<double>>& values, const std::vector<double>& weights) {
    size_t N = values.size();
    size_t B = values[0].size();
    std::vector<double> result(B, 0.0);

    double wsum = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (wsum == 0.0) return result;

    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < B; ++j)
            result[j] += values[i][j] * weights[i];

    for (double& x : result)
        x /= wsum;

    return result;
}

// --- Jackknife Errors ---
double Observables::compute_scalar_jackknife_error(const std::vector<double>& vals, const std::vector<double>& weights) {
    const size_t N = vals.size();
    if (N <= 1 || vals.empty()){
        std::cerr << "[jackknife] Warning: Not enough samples (N = " << N
                  << ") to compute jackknife error. Returning 0.0.\n";
        return 0.0;
    } 

    // Precompute total sums
    double total_weight = 0.0;
    double total_weighted_sum = 0.0;
    for (size_t i = 0; i < N; ++i) {
        total_weight += weights[i];
        total_weighted_sum += weights[i] * vals[i];
    }

    // Compute jackknife leave-one-out means efficiently
    std::vector<double> jk_means(N);
    for (size_t i = 0; i < N; ++i) {
        double w_i = weights[i];
        double x_i = vals[i];

        double sub_weight = total_weight - w_i;
        if (sub_weight <= 0.0) {
            jk_means[i] = 0.0; // fallback, should rarely happen
        } else {
            double sub_sum = total_weighted_sum - w_i * x_i;
            jk_means[i] = sub_sum / sub_weight;
        }
    }

    // Compute mean of jackknife estimates
    double mean = std::accumulate(jk_means.begin(), jk_means.end(), 0.0) / N;

    // Compute jackknife variance
    double err2 = 0.0;
    for (size_t i = 0; i < N; ++i) {
        double diff = jk_means[i] - mean;
        err2 += diff * diff;
    }

    return std::sqrt((N - 1.0) / N * err2);
}


std::vector<double> Observables::compute_vector_jackknife_error(
    const std::vector<std::vector<double>>& vals,
    const std::vector<double>& weights) {

    const size_t N = vals.size();
    if (N <= 1 || vals.empty() || vals[0].empty()) {
        std::cerr << "[jackknife] Warning: Not enough samples (N = " << N
                  << ") to compute jackknife error. Returning zeros.\n";
        return std::vector<double>(vals.empty() ? 0 : vals[0].size(), 0.0);
    }

    const size_t B = vals[0].size();  // number of bins/components
    std::vector<double> total_weights(B, 0.0);
    std::vector<double> total_weighted_sum(B, 0.0);

    // Precompute total weighted sums per bin
    for (size_t i = 0; i < N; ++i) {
        double wi = weights[i];
        for (size_t b = 0; b < B; ++b) {
            total_weights[b] += wi;
            total_weighted_sum[b] += wi * vals[i][b];
        }
    }

    // Compute leave-one-out means
    std::vector<std::vector<double>> jk_means(N, std::vector<double>(B, 0.0));
    for (size_t i = 0; i < N; ++i) {
        double wi = weights[i];
        for (size_t b = 0; b < B; ++b) {
            double w_tot = total_weights[b] - wi;
            if (w_tot <= 0.0) {
                jk_means[i][b] = 0.0;
            } else {
                double sum_excl = total_weighted_sum[b] - wi * vals[i][b];
                jk_means[i][b] = sum_excl / w_tot;
            }
        }
    }

    // Compute jackknife mean and error
    std::vector<double> mean(B, 0.0), err(B, 0.0);
    for (size_t b = 0; b < B; ++b) {
        for (size_t i = 0; i < N; ++i)
            mean[b] += jk_means[i][b];
        mean[b] /= N;

        for (size_t i = 0; i < N; ++i) {
            double diff = jk_means[i][b] - mean[b];
            err[b] += diff * diff;
        }

        err[b] = std::sqrt((N - 1.0) / N * err[b]);
    }

    return err;
}
