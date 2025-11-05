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

    std::vector<double> centrality_edges;  

    std::string database_path;
    std::string centrality_dict_path;

    std::string output_directory = "output";
    std::string output_prefix = "";
    std::string output_format = "dat";
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

        
        cfg.output_directory = node["output"]["directory"] ? node["output"]["directory"].as<std::string>() : "output";
        cfg.output_prefix    = node["output"]["prefix"] ? node["output"]["prefix"].as<std::string>() : "";
        cfg.output_format    = node["output"]["format"] ? node["output"]["format"].as<std::string>() : "dat";

        cfg.database_path = node["centrality"]["database_path"].as<std::string>();
        cfg.centrality_dict_path = node["centrality"]["centrality_dict_path"].as<std::string>();
        if (node["centrality"]["bins"]) {
            cfg.centrality_edges = node["centrality"]["bins"].as<std::vector<double>>();
        }

        return cfg;
    }
};


#endif
