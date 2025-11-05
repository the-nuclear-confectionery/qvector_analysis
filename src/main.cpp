#include "event.h"
#include "cuts.h"
#include "observables.h"
#include "integrated_observables.h"
#include "differential_observables.h"
#include "output.h"
#include "config.h"
#include "centrality.h"

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <optional>
#include <unordered_map>
#include <fstream>
#include <filesystem>
#include <sqlite3.h>

namespace fs = std::filesystem;

struct DbHandle {
    sqlite3* db = nullptr;
    ~DbHandle() { if (db) sqlite3_close(db); }
};

static std::optional<double> fetch_centrality(sqlite3* db, int event_id) {
    const char* SQL = "SELECT centrality_estimator FROM events WHERE event_id = ?1;";
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db, SQL, -1, &stmt, nullptr) != SQLITE_OK) {
        std::cerr << "[sqlite] prepare failed: " << sqlite3_errmsg(db) << "\n";
        return {};
    }
    sqlite3_bind_int(stmt, 1, event_id);
    std::optional<double> result;
    if (sqlite3_step(stmt) == SQLITE_ROW &&
        sqlite3_column_type(stmt, 0) != SQLITE_NULL)
    {
        result = sqlite3_column_double(stmt, 0);
    }
    sqlite3_finalize(stmt);
    return result;
}

static int bin_from_edges(double x, const std::vector<double>& edges) {
    if (edges.size() < 2) return -1;
    for (size_t i = 0; i + 1 < edges.size(); ++i) {
        if (x >= edges[i] && x < edges[i+1]) return static_cast<int>(i);
    }
    if (x == edges.back()) return static_cast<int>(edges.size() - 2);
    return -1;
}

static std::string bin_folder_label(int idx, const std::vector<double>& edges) {
    if (idx < 0 || static_cast<size_t>(idx+1) >= edges.size()) return "unknown";
    int a = static_cast<int>(edges[idx]);
    int b = static_cast<int>(edges[idx+1]);
    char buf[16];
    snprintf(buf, sizeof(buf), "%d%d", a, b);
    return std::string(buf);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " config.yaml file1.root [...]\n";
        return 1;
    }

    config cfg = config::load(argv[1]);

    if (cfg.calculate_charged &&
        std::find(cfg.pids.begin(), cfg.pids.end(), 0) == cfg.pids.end()) {
        cfg.pids.insert(cfg.pids.begin(), 0);
    }

    if (cfg.database_path.empty() || cfg.centrality_dict_path.empty()) {
        std::cerr << "Error: database_path or centrality_dict_path missing in YAML\n";
        return 2;
    }

    std::vector<double> user_edges = cfg.centrality_edges;
    if (user_edges.size() < 2) user_edges = {0.0, 100.0};

    std::vector<CentBin> cent_bins = load_centrality_dict(cfg.centrality_dict_path);
    if (cent_bins.empty()) return 3;

    DbHandle db;
    if (sqlite3_open(cfg.database_path.c_str(), &db.db) != SQLITE_OK) {
        std::cerr << "Error opening DB: " << cfg.database_path << "\n";
        return 4;
    }

    std::vector<event> events;
    events.reserve(argc - 2);

    for (int i = 2; i < argc; ++i) {
        std::string path = argv[i];
        event ev(path.c_str(), cfg);

        auto eid = extract_event_id_from_path(path);
        if (eid) {
            if (auto est = fetch_centrality(db.db, *eid)) {
                int cbin = find_centrality(*est, cent_bins, user_edges);
                if (cbin >= 0) {
                    ev.centrality = 0.5 * (user_edges[cbin] + user_edges[cbin + 1]);
                } else {
                    ev.centrality = 0.0;
                }
            } else {
                ev.centrality = 0.0;
            }
        } else {
            ev.centrality = 0.0;
        }

        events.push_back(std::move(ev));
    }

    Observables obs_base;
    Output base_out(cfg, obs_base, events);
    std::string base_dir = base_out.build_output_path();
    fs::create_directories(base_dir);

    const size_t NB = user_edges.size() - 1;
    std::vector<std::vector<size_t>> bin_indices(NB);
    size_t skipped = 0;

    for (size_t i = 0; i < events.size(); ++i) {
        double c = events[i].centrality;
        int b = bin_from_edges(c, user_edges);
        if (b >= 0 && (size_t)b < NB) bin_indices[b].push_back(i);
        else ++skipped;
    }

    std::cerr << "[centrality] skipped (no/invalid centrality): " << skipped << "\n";

    for (size_t b = 0; b < NB; ++b) {
        if (bin_indices[b].empty()) continue;

        std::vector<event> ev_subset;
        ev_subset.reserve(bin_indices[b].size());
        for (size_t idx : bin_indices[b]) ev_subset.push_back(events[idx]);

        Observables obs_bin;

        for (const std::string& name : cfg.integrated_observables) {
            if (name == "M")             obs_bin.register_integrated(name, M);
            else if (name == "mean_pT")  obs_bin.register_integrated(name, mean_pT);
            else if (name == "vn{EP}")   obs_bin.register_integrated(name, vnEP);
            else if (name == "vn{2}")    obs_bin.register_integrated(name, vn2);
            else if (name == "vn{4}")    obs_bin.register_integrated(name, vn4);
            else std::cerr << "Warning: Unknown integrated observable '" << name << "'\n";
        }

        for (const std::string& name : cfg.differential_observables) {
            if (name == "dN/deta")            obs_bin.register_differential(name, dN_deta);
            else if (name == "dN/2piptdptdy") obs_bin.register_differential(name, dN_2piptdpTdy);
            else if (name == "vn{2}_pt")      obs_bin.register_differential(name, vn2_pt);
            else if (name == "vn{2}_eta")     obs_bin.register_differential(name, vn2_eta);
            else std::cerr << "Warning: Unknown differential observable '" << name << "'\n";
        }

        obs_bin.evaluate_all_integrated(ev_subset, cfg.cut, cfg.pids, cfg.max_n);
        obs_bin.evaluate_all_differential(ev_subset, cfg.cut, cfg.pids, cfg.max_n);

        Output out_bin(cfg, obs_bin, ev_subset);

        std::string clabel = bin_folder_label((int)b, user_edges);
        double cmid = 0.5 * (user_edges[b] + user_edges[b+1]);

        for (int pid : cfg.pids) {
            out_bin.append_integrated_row(pid, base_dir, cmid);
        }

        for (int pid : cfg.pids) {
            out_bin.write_all_differential_for_bin(pid, base_dir, clabel);
        }
    }

    return 0;
}
