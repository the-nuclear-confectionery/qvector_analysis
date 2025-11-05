#ifndef EVENT_H
#define EVENT_H

#include "config.h"
#include <string>
#include <vector>
#include <complex>
#include <map>
#include <TH2D.h>

class event {
public:
    event(const std::string& filename, config& cfg);

    std::string event_id;
    int nsamples = 1;
    int detected_n_max = 0;
    std::vector<int> available_pids;
    std::vector<double> eta_centers;
    std::vector<double> pt_centers;
    double centrality = -1.0;

    std::map<int, std::vector<std::vector<std::vector<std::complex<double>>>>> Qn; // [pid][n][eta][pt]
    std::map<int, std::vector<std::vector<double>>> Q0;                           // [pid][eta][pt]

    std::vector<std::vector<std::complex<double>>> get_Qn(int pid, int n) const;

    std::complex<double> get_Qn_integrated_over_eta_and_pt(int pid, int n,
        double eta_min, double eta_max, double pt_min, double pt_max) const;

    std::vector<std::complex<double>> get_Qn_integrated_over_eta(int pid, int n,
        double eta_min, double eta_max, double pt_min, double pt_max) const;

    std::vector<std::complex<double>> get_Qn_integrated_over_pt(int pid, int n,
        double eta_min, double eta_max, double pt_min, double pt_max) const;

    std::vector<std::vector<double>> get_Q0(int pid) const;

    double get_Q0_integrated_over_eta_and_pt(int pid,
        double eta_min, double eta_max, double pt_min, double pt_max) const;

    std::vector<double> get_Q0_integrated_over_eta(int pid,
        double eta_min, double eta_max, double pt_min, double pt_max) const;

    std::vector<double> get_Q0_integrated_over_pt(int pid,
        double eta_min, double eta_max, double pt_min, double pt_max) const;

private:
    void read_file(const std::string& filename, config& cfg);
    int read_qvectors(int pid);
    void init_binning(const TH2D* hist);
    bool binning_initialized = false;
};

#endif
