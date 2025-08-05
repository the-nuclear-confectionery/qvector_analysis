#include "event.h"
#include "cuts.h"
#include "observables.h"
#include "integrated_observables.h"
#include "differential_observables.h"
#include "output.h"
#include "config.h"

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

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
    std::vector<std::string> integrated_registered;
    std::vector<std::string> differential_registered;

    // --- Register integrated observables ---
    for (const std::string& name : cfg.integrated_observables) {
        if (name == "M")             obs.register_integrated(name, M);
        else if (name == "mean_pT")  obs.register_integrated(name, mean_pT);
        else if (name == "vn{EP}")   obs.register_integrated(name, vnEP);
        else if (name == "vn{2}")    obs.register_integrated(name, vn2);
        else if (name == "vn{4}")    obs.register_integrated(name, vn4);
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

    obs.evaluate_all_integrated(events, cfg.cut, cfg.pids, cfg.max_n);
    obs.evaluate_all_differential(events, cfg.cut, cfg.pids, cfg.max_n);


    Output output(cfg, obs, events);
    std::string output_dir = output.build_output_path();

    for (int pid : cfg.pids) {
        output.write_integrated(pid, output_dir);
        output.write_all_differential(pid, output_dir);
    }


    return 0;
}