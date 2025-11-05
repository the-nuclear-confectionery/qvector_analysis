#ifndef CENTRALITY_H
#define CENTRALITY_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <optional>
#include <limits>
#include <iostream>
#include <regex>  

struct CentBin {
    double p_lo, p_hi;   // fine percentile bounds (e.g., 0,1)
    double e_lo, e_hi;   // estimator edges for that fine bin
};

inline std::vector<CentBin> load_centrality_dict(const std::string& path) {
    std::vector<CentBin> bins;
    std::ifstream f(path);
    if (!f) { std::cerr << "[centrality] cannot open dict: " << path << "\n"; return bins; }
    std::string line; size_t ln=0; 
    while (std::getline(f, line)) {
        ++ln; if (line.empty() || line[0]=='#') continue;
        CentBin b{};
        std::istringstream iss(line);
        if (iss >> b.p_lo >> b.p_hi >> b.e_lo >> b.e_hi) bins.push_back(b);
        else std::cerr << "[centrality] bad line " << ln << ": " << line << "\n";
    }
    return bins;
}

// ---- Single-shot mapper ----
// Returns the user coarse-bin index for the given estimator `v`.
// Steps: (1) find fine dict bin; (2) take pct = p_hi; (3) map into [edges[i],edges[i+1)).
// If v falls outside all dict bins, it clamps to the nearest fine bin first.
// Returns -1 on failure (e.g., empty dict or <2 edges).
inline int find_centrality(double v,
                           const std::vector<CentBin>& dict,
                           const std::vector<double>& edges)
{
    if (dict.empty() || edges.size() < 2) return -1;

    // 1) find fine bin (support ascending/descending estimator edges)
    const CentBin* pick = nullptr;
    for (const auto& b : dict) {
        const double lo = std::min(b.e_lo, b.e_hi);
        const double hi = std::max(b.e_lo, b.e_hi);
        if (v >= lo && v < hi) { pick = &b; break; }
    }

    // If not found, clamp to nearest interval
    if (!pick) {
        double best = std::numeric_limits<double>::infinity();
        for (const auto& b : dict) {
            const double lo = std::min(b.e_lo, b.e_hi);
            const double hi = std::max(b.e_lo, b.e_hi);
            double d = (v < lo) ? (lo - v) : (v >= hi ? (v - hi) : 0.0);
            if (d < best) { best = d; pick = &b; }
        }
        if (!pick) return -1;
    }

    // 2) pick a single percentile representative for fine bin
    const double pct = pick->p_hi;  // e.g., fine [3–4] → use 4

    // 3) map to user coarse bin via edges
    for (size_t i = 0; i + 1 < edges.size(); ++i) {
        if (pct >= edges[i] && pct < edges[i + 1]) return static_cast<int>(i);
    }
    if (pct == edges.back()) return static_cast<int>(edges.size() - 2); // inclusive top edge

    return -1;
}

// (nice-to-have) label helper: idx->"A-B"
inline std::string centrality_label(int idx, const std::vector<double>& edges) {
    if (idx < 0 || static_cast<size_t>(idx+1) >= edges.size()) return "";
    int a = static_cast<int>(edges[idx]);
    int b = static_cast<int>(edges[idx+1]);
    return std::to_string(a) + "-" + std::to_string(b);
}


// ======================
// Event ID discovery
// ======================

// Extract last integer found in filename (e.g. "event_18.root" → 18)
inline std::optional<int> extract_event_id_from_path(const std::string& path) {
    std::regex re("(\\d+)(?!.*\\d)");
    std::smatch m;
    if (std::regex_search(path, m, re)) {
        try { return std::stoi(m.str(1)); }
        catch (...) { return std::nullopt; }
    }
    return std::nullopt;
}

#endif // CENTRALITY_H
