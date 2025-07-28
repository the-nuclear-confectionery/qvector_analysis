
#include "event.h"
#include "cuts.h"
#include "observables.h"
#include "integrated_observables.h"
#include "differential_observables.h"
#include "config.h"

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " config.yaml file1.root [file2.root ...]\n";
        return 1;
    }

    // === Load configuration ===
    config cfg = config::load(argv[1]);

    // Add charged PID (0) if enabled
    if (cfg.calculate_charged && std::find(cfg.pids.begin(), cfg.pids.end(), 0) == cfg.pids.end()) {
        cfg.pids.insert(cfg.pids.begin(), 0);  // add to beginning (optional)
    }

    // === Load events ===
    std::vector<event> events;
    for (int i = 2; i < argc; ++i) {
        event ev(argv[i], cfg);
        events.push_back(std::move(ev));
    }

    // === Setup observables ===
    Observables obs;

    // --- Register integrated observables ---
    for (const std::string& name : cfg.integrated_observables) {
        if (name == "M")             obs.register_integrated(name, M);
        else if (name == "mean_pT")  obs.register_integrated(name, mean_pT);
        else if (name == "vn{EP}")   obs.register_integrated(name, vnEP);
        else if (name == "vn{2}")    obs.register_integrated(name, vn2);
        else if (name == "vn{4}")    obs.register_integrated(name, vn4);  // will be removed if unsupported
        else std::cerr << "Warning: Unknown integrated observable '" << name << "'\n";
    }

    // --- Register differential observables ---
    for (const std::string& name : cfg.differential_observables) {
        if (name == "dN/deta")            obs.register_differential(name, dN_deta);
        else if (name == "dN/2piptdptdy") obs.register_differential(name, dN_2piptdpTdy);
        else if (name == "vn{2}_pt")      obs.register_differential(name, vn2_pt);
        else if (name == "vn{2}_eta")     obs.register_differential(name, vn2_eta);
        else std::cerr << "Warning: Unknown differential observable '" << name << "'\n";
    }

    // === Run observable evaluation ===
    obs.evaluate_all_integrated(events, cfg.cut, cfg.pids, cfg.max_n);
    obs.evaluate_all_differential(events, cfg.cut, cfg.pids, cfg.max_n);

    // === Output ===
    std::cout << "=== Integrated Observables ===\n";
    for (const auto& [name, value] : obs.scalar_values) {
        double err = obs.scalar_errors.at(name);
        std::cout << "  " << name << " = " << value << " ± " << err << "\n";
    }

    std::cout << "\n=== Differential Observables ===\n";
    for (const auto& [name, vec] : obs.vector_values) {
        const auto& err = obs.vector_errors.at(name);
        std::cout << "  " << name << ":\n";

        std::vector<double> bin_centers;

        if (name.find("pt") != std::string::npos) {
            bin_centers = events[0].pt_centers;
        } else if (name.find("eta") != std::string::npos) {
            bin_centers = events[0].eta_centers;
        } else {
            bin_centers.resize(vec.size());
            std::iota(bin_centers.begin(), bin_centers.end(), 0); // fallback index
        }


        for (size_t i = 0; i < vec.size(); ++i) {
            std::cout << " " << bin_centers[i]
                      << " " << vec[i] << " " << err[i] << "\n";
        }
    }


    return 0;
}
