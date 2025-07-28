#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <string>
#include <yaml-cpp/yaml.h>
#include "cuts.h"

struct config {
    std::vector<int> pids;
    int max_n;
    cuts cut;
    bool calculate_charged = false;
    std::vector<std::string> integrated_observables;
    std::vector<std::string> differential_observables;

    static config load(const std::string& filename) {
        config cfg;
        YAML::Node node = YAML::LoadFile(filename);

        cfg.pids = node["pids"].as<std::vector<int>>();
        cfg.max_n = node["max_n"].as<int>();
        cfg.cut.eta_cut = node["cuts"]["eta_cut"].as<double>();
        cfg.cut.pt_min  = node["cuts"]["pt_min"].as<double>();
        cfg.cut.pt_max  = node["cuts"]["pt_max"].as<double>();
        cfg.calculate_charged = node["calculate_charged"].as<bool>();

        cfg.integrated_observables = node["observables"]["integrated"].as<std::vector<std::string>>();
        cfg.differential_observables = node["observables"]["differential"].as<std::vector<std::string>>();

        return cfg;
    }
};


#endif
