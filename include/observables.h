#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include "event.h"
#include "cuts.h"
#include "config.h"
#include <cmath>
#include <iostream>
#include <complex>
#include <string>
#include <vector>
#include <unordered_map>
#include <functional>
#include <utility>

using IntegratedObservable = std::function<std::pair<double, double>(event&, cuts&, int, int)>;
using DifferentialObservable = std::function<std::vector<std::pair<double, double>>(event&, cuts&, int, int)>;

class Observables {
public:
    void register_integrated(const std::string& name, IntegratedObservable fn);
    void register_differential(const std::string& name, DifferentialObservable fn);

    void evaluate_all_integrated(const std::vector<event>& events, cuts& cut,
                                 const std::vector<int>& pids, int n_max);

    void evaluate_all_differential(const std::vector<event>& events, cuts& cut,
                                    const std::vector<int>& pids, int n_max);


    std::unordered_map<std::string, double> scalar_values;
    std::unordered_map<std::string, double> scalar_errors;
    std::unordered_map<std::string, std::vector<double>> vector_values;
    std::unordered_map<std::string, std::vector<double>> vector_errors;

private:
    std::unordered_map<std::string, IntegratedObservable> integrated_observables;
    std::unordered_map<std::string, DifferentialObservable> differential_observables;

    // helper functions
    double weighted_average(const std::vector<double>& values, const std::vector<double>& weights);
    std::vector<double> weighted_average(const std::vector<std::vector<double>>& values, const std::vector<double>& weights);
    double compute_scalar_jackknife_error(const std::vector<double>& values, const std::vector<double>& weights);
    std::vector<double> compute_vector_jackknife_error(const std::vector<std::vector<double>>& values, const std::vector<double>& weights);
};

#endif
