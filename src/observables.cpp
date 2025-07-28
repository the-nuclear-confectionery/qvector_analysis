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


bool Observables::requires_n(const std::string& name) const {
    return name.find("vn") != std::string::npos;  
}

void Observables::evaluate_all_integrated(const std::vector<event>& events, cuts& cut,
                                          const std::vector<int>& pids, int n_max) {
    std::unordered_map<std::string, std::vector<double>> value_acc;
    std::unordered_map<std::string, std::vector<double>> weight_acc;

    for (const event& ev : events) {
        for (const auto& [name, fn] : integrated_observables) {
            for (int pid : pids) {
                if (requires_n(name)) {
                    for (int n = 1; n <= n_max; ++n) {
                        auto [val, w] = fn(const_cast<event&>(ev), cut, pid, n);
                        std::string key = name + "_pid" + std::to_string(pid) + "_n" + std::to_string(n);
                        value_acc[key].push_back(val);
                        weight_acc[key].push_back(w);
                    }
                } else {
                    auto [val, w] = fn(const_cast<event&>(ev), cut, pid, 0);
                    std::string key = name + "_pid" + std::to_string(pid);
                    value_acc[key].push_back(val);
                    weight_acc[key].push_back(w);
                }
            }
        }
    }

    for (const auto& [key, vals] : value_acc) {
        const auto& wgts = weight_acc[key];
        scalar_values[key] = weighted_average(vals, wgts);
        scalar_errors[key] = compute_scalar_jackknife_error(vals, wgts);
    }
}

void Observables::evaluate_all_differential(const std::vector<event>& events, cuts& cut,
                                            const std::vector<int>& pids, int n_max) {
    std::unordered_map<std::string, std::vector<std::vector<double>>> val_acc;
    std::unordered_map<std::string, std::vector<std::vector<double>>> wgt_acc;

    for (const event& ev : events) {
        for (const auto& [name, fn] : differential_observables) {
            for (int pid : pids) {
                if (requires_n(name)) {
                    for (int n = 1; n <= n_max; ++n) {
                        auto vec = fn(const_cast<event&>(ev), cut, pid, n);
                        std::vector<double> vals, wgts;
                        for (const auto& [v, w] : vec) {
                            vals.push_back(v);
                            wgts.push_back(w);
                        }

                        std::string key = name + "_pid" + std::to_string(pid) + "_n" + std::to_string(n);
                        val_acc[key].push_back(vals);
                        wgt_acc[key].push_back(wgts);
                    }
                } else {
                    auto vec = fn(const_cast<event&>(ev), cut, pid, 0);
                    std::vector<double> vals, wgts;
                    for (const auto& [v, w] : vec) {
                        vals.push_back(v);
                        wgts.push_back(w);
                    }

                    std::string key = name + "_pid" + std::to_string(pid);
                    val_acc[key].push_back(vals);
                    wgt_acc[key].push_back(wgts);
                }
            }
        }
    }

    for (const auto& [key, vals] : val_acc) {
        const auto& wgts = wgt_acc[key];
        vector_values[key] = weighted_average(vals, wgts);
        vector_errors[key] = compute_vector_jackknife_error(vals, wgts);
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

std::vector<double> Observables::weighted_average(
    const std::vector<std::vector<double>>& values,
    const std::vector<std::vector<double>>& weights) {

    size_t N = values.size();    // number of events
    size_t B = values[0].size(); // number of bins

    std::vector<double> result(B, 0.0);
    std::vector<double> weight_sums(B, 0.0);

    for (size_t i = 0; i < N; ++i) {
        for (size_t b = 0; b < B; ++b) {
            result[b] += values[i][b] * weights[i][b];
            weight_sums[b] += weights[i][b];
        }
    }

    for (size_t b = 0; b < B; ++b) {
        if (weight_sums[b] > 0.0)
            result[b] /= weight_sums[b];
        else
            result[b] = 0.0;
    }

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
    const std::vector<std::vector<double>>& weights) {

    const size_t N = vals.size();
    if (N <= 1 || vals.empty() || vals[0].empty()) {
        std::cerr << "[jackknife] Warning: Not enough samples (N = " << N
                  << ") to compute jackknife error. Returning zeros.\n";
        return std::vector<double>(vals.empty() ? 0 : vals[0].size(), 0.0);
    }

    const size_t B = vals[0].size();  // number of bins
    std::vector<double> total_weights(B, 0.0);
    std::vector<double> total_weighted_sum(B, 0.0);

    // Precompute total weighted sums per bin
    for (size_t i = 0; i < N; ++i) {
        for (size_t b = 0; b < B; ++b) {
            total_weights[b] += weights[i][b];
            total_weighted_sum[b] += weights[i][b] * vals[i][b];
        }
    }

    // Compute jackknife leave-one-out means
    std::vector<std::vector<double>> jk_means(N, std::vector<double>(B, 0.0));
    for (size_t i = 0; i < N; ++i) {
        for (size_t b = 0; b < B; ++b) {
            double w_tot = total_weights[b] - weights[i][b];
            if (w_tot <= 0.0) {
                jk_means[i][b] = 0.0;
            } else {
                double sum_excl = total_weighted_sum[b] - weights[i][b] * vals[i][b];
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
