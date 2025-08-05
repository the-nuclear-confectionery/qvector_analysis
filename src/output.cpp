#include "output.h"
#include <fstream>
#include <iomanip>
#include <regex>
#include <map>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <filesystem>
namespace fs = std::filesystem;

Output::Output(const config& cfg, const Observables& obs, const std::vector<event>& events)
    : cfg_(cfg), obs_(obs), events_(events) {}

std::string Output::build_output_path() const {
    char folder[256];
    snprintf(folder, sizeof(folder), "%s_pt%.1f-%.1f_eta%.1f",
             cfg_.output_prefix.c_str(),
             cfg_.cut.pt_min, cfg_.cut.pt_max,
             cfg_.cut.eta_cut);

    fs::path full_path = fs::path(cfg_.output_directory) / folder;
    if (!fs::exists(full_path))
        fs::create_directories(full_path);

    return full_path.string();
}

void Output::write_all() const {
    std::string outdir = build_output_path();
    for (int pid : cfg_.pids) {
        write_integrated(pid, outdir);
        write_all_differential(pid, outdir);
    }
}

void Output::write_integrated(int pid, const std::string& outdir) const {
    std::string filename = outdir + "/integrated_pid" + std::to_string(pid) + ".dat";
    std::ofstream out(filename);
    out << std::setprecision(8) << std::scientific;

    std::vector<std::string> order = {"M", "mean_pT", "vn{EP}", "vn{2}", "vn{4}"};

    std::map<std::string, std::pair<double, double>> sorted_obs;
    for (const auto& [name, val] : obs_.scalar_values) {
        if (name.find("pid" + std::to_string(pid)) != std::string::npos) {
            double err = obs_.scalar_errors.count(name) ? obs_.scalar_errors.at(name) : 0.0;
            sorted_obs[name] = std::make_pair(val, err);
        }
    }

    out << "centrality";
    for (const auto& type : order) {
        for (const auto& [name, _] : sorted_obs) {
            if (name.find(type + "_") == 0) {
                std::string label = std::regex_replace(name, std::regex("vn\\{(\\d+)\\}"), "v$1");
                out << " " << label << " " << label << "_err";
            }
        }
    }
    out << "\n";

    out << "0";
    for (const auto& type : order) {
        for (const auto& [name, valerr] : sorted_obs) {
            if (name.find(type + "_") == 0) {
                out << " " << valerr.first << " " << valerr.second;
            }
        }
    }
    out << "\n";
}

void Output::write_all_differential(int pid, const std::string& outdir) const {
    std::map<std::string, std::map<int, std::string>> harmonic_groups;
    std::vector<std::string> dndx_names;

    std::regex harmonic_pattern("(vn\\{\\d+\\}_(pt|eta)_pid" + std::to_string(pid) + ")_n(\\d+)");
    for (const auto& [name, _] : obs_.vector_values) {
        std::smatch match;
        if (std::regex_match(name, match, harmonic_pattern)) {
            std::string base = match[1];
            int n = std::stoi(match[3]);
            harmonic_groups[base][n] = name;
        } else if ((name.find("dN/deta_pid" + std::to_string(pid)) != std::string::npos) ||
                   (name.find("dN/2piptdptdy_pid" + std::to_string(pid)) != std::string::npos)) {
            dndx_names.push_back(name);
        }
    }

    // Write vn groups
    for (const auto& [base, harmonics] : harmonic_groups) {
        std::string safe_base = base;
        std::replace(safe_base.begin(), safe_base.end(), '/', '_');
        std::string filename = outdir + "/" + safe_base + ".dat";

        std::ofstream out(filename);
        out << std::setprecision(8) << std::scientific;

        std::vector<double> bin_centers;
        if (base.find("_pt_") != std::string::npos)
            bin_centers = events_[0].pt_centers;
        else if (base.find("_eta_") != std::string::npos)
            bin_centers = events_[0].eta_centers;
        else {
            size_t n = obs_.vector_values.begin()->second.size();
            bin_centers.resize(n);
            std::iota(bin_centers.begin(), bin_centers.end(), 0);
        }

        out << "# bin";
        for (const auto& [n, name] : harmonics) {
            std::string label = "v" + std::to_string(n);
            out << " " << label << " " << label << "_err";
        }
        out << "\n";

        for (size_t i = 0; i < bin_centers.size(); ++i) {
            out << bin_centers[i];
            for (const auto& [n, name] : harmonics) {
                double val = obs_.vector_values.at(name)[i];
                double err = obs_.vector_errors.at(name)[i];
                out << " " << val << " " << err;
            }
            out << "\n";
        }
    }

    // Write dN/deta and dN/2piptdptdy
    for (const auto& name : dndx_names) {
        std::string safe_name = name;
        std::replace(safe_name.begin(), safe_name.end(), '/', '_');
        std::string filename = outdir + "/" + safe_name + ".dat";

        std::ofstream out(filename);
        out << std::setprecision(8) << std::scientific;

        const auto& vec = obs_.vector_values.at(name);
        const auto& err = obs_.vector_errors.at(name);

        std::vector<double> bin_centers;
        if (name.find("pt") != std::string::npos)
            bin_centers = events_[0].pt_centers;
        else if (name.find("eta") != std::string::npos)
            bin_centers = events_[0].eta_centers;
        else {
            bin_centers.resize(vec.size());
            std::iota(bin_centers.begin(), bin_centers.end(), 0);
        }

        out << "# bin val err\n";
        for (size_t i = 0; i < vec.size(); ++i) {
            out << bin_centers[i] << " " << vec[i] << " " << err[i] << "\n";
        }
    }
}
