#ifndef OUTPUT_H
#define OUTPUT_H

#include "config.h"
#include "observables.h"
#include "event.h"
#include <vector>
#include <string>
#include <vector>
#include <string>
#include <numeric>

class Output {
public:
    Output(const config& cfg, const Observables& obs, const std::vector<event>& events);

    std::string build_output_path() const;

    void write_all() const;
    void write_integrated(int pid, const std::string& outdir) const;
    void write_all_differential(int pid, const std::string& outdir) const;

    // ------------------------------------------------------------
    // New API (keeps old API intact)
    // ------------------------------------------------------------

    // Append ONE row (this centrality bin) into integrated_pidX.dat at cut folder.
    // The header is created automatically if the file does not exist yet.
    void append_integrated_row(int pid, const std::string& cut_dir, double centrality_mid) const;

    // Write all differential vectors for this bin into:
    //    <cut_dir>/differential/<cent_label>/
    // where cent_label is e.g. "05", "510", "1020", ...
    void write_all_differential_for_bin(int pid, const std::string& cut_dir, const std::string& cent_label) const;

private:
    const config& cfg_;
    const Observables& obs_;
    const std::vector<event>& events_;

    // Ensure header exists in integrated_pidX.dat (idempotent)
    void ensure_integrated_header_(int pid, const std::string& path) const;
};

#endif // OUTPUT_H
