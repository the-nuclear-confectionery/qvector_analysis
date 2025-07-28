#ifndef DIFFERENTIAL_OBSERVABLES_H
#define DIFFERENTIAL_OBSERVABLES_H

#include "event.h"
#include "cuts.h"
#include <cmath>
#include <complex>
#include <utility>
#include <vector>

// === dN/deta ===
inline std::vector<std::pair<double, double>> dN_deta(event& ev, cuts& cut, int pid, int dummy = 0) {
    std::vector<double> dNdeta = ev.get_Q0_integrated_over_pt(pid, -cut.eta_cut, cut.eta_cut, cut.pt_min, cut.pt_max);
    //normalize by event number of samples
    double nsamples = static_cast<double>(ev.nsamples);
    std::vector<std::pair<double, double>> result;
    for (size_t i = 0; i < dNdeta.size(); ++i) {
        if (std::abs(ev.eta_centers[i]) > cut.eta_cut) continue;
        result.emplace_back(dNdeta[i] / nsamples, 1.0); // weight is 1.0
    }
    return result;
}

// === dN/2pi pT dpT dy ===
inline std::vector<std::pair<double, double>> dN_2piptdpTdy(event& ev, cuts& cut, int pid,  int dummy = 0) {
    std::vector<double> dNdpt = ev.get_Q0_integrated_over_eta(pid, -cut.eta_cut, cut.eta_cut, cut.pt_min, cut.pt_max);
    std::vector<std::pair<double, double>> result;
    double nsamples = static_cast<double>(ev.nsamples);
    double dy = 2.0 * cut.eta_cut; // rapidity range
    for (size_t j = 0; j < dNdpt.size(); ++j) {
        double pt = ev.pt_centers[j];
        if (pt < cut.pt_min || pt > cut.pt_max) continue;
        double val = dNdpt[j] / (2 * M_PI * pt * dy) / nsamples; // normalize by event number of samples
        result.emplace_back(val, 1.0);
    }
    return result;
}

// === vn{2}(pt) ===
inline std::vector<std::pair<double, double>> vn2_pt(event& ev, cuts& cut, int pid, int n) {
    const auto& Qn_pt = ev.get_Qn_integrated_over_eta(pid, n, -cut.eta_cut, cut.eta_cut, cut.pt_min, cut.pt_max);
    const auto& Q0_pt = ev.get_Q0_integrated_over_eta(pid, -cut.eta_cut, cut.eta_cut, cut.pt_min, cut.pt_max);
    const auto& Qn_all = ev.get_Qn(pid, n);
    const auto& Q0_all = ev.get_Q0(pid);

    std::complex<double> Qn_ref = {0.0, 0.0};
    double Mq = 0.0;
    for (size_t i = 0; i < Qn_all.size(); ++i) {
        for (size_t j = 0; j < Qn_all[i].size(); ++j) {
            double eta = ev.eta_centers[i];
            double pt = ev.pt_centers[j];
            if (pt < cut.pt_min || pt > cut.pt_max) continue;
            if (std::abs(eta) > cut.eta_cut) continue;
            Qn_ref += Qn_all[i][j];
            Mq += Q0_all[i][j];
        }
    }

    double c2_ref = 0.0;
    if (Mq > 1.0)
        c2_ref = (std::norm(Qn_ref) - Mq) / (Mq * (Mq - 1.0));

    std::vector<std::pair<double, double>> result;
    for (size_t j = 0; j < Qn_pt.size(); ++j) {
        double pt = ev.pt_centers[j];
        if (pt < cut.pt_min || pt > cut.pt_max) continue;

        std::complex<double> Pn = Qn_pt[j];
        double Mp = Q0_pt[j];
        double Mpq = Q0_pt[j]; // perfect overlap assumed

        double denom = Mp * (Mq - Mpq);
        double two_prime = 0.0;
        if (denom > 0.0)
            two_prime = (std::real(Pn * std::conj(Qn_ref)) - Mpq) / denom;

        double vn2_val = (c2_ref > 0.0) ? two_prime / std::sqrt(c2_ref) : 0.0;
        result.emplace_back(vn2_val, Mq * (Mq - 1.0)); // w2 is total
    }
    return result;
}

// === vn{2}(eta) ===
inline std::vector<std::pair<double, double>> vn2_eta(event& ev, cuts& cut, int pid, int n) {
    const auto& Qn_eta = ev.get_Qn_integrated_over_pt(pid, n, -cut.eta_cut, cut.eta_cut, cut.pt_min, cut.pt_max);
    const auto& Q0_eta = ev.get_Q0_integrated_over_pt(pid, -cut.eta_cut, cut.eta_cut, cut.pt_min, cut.pt_max);
    const auto& Qn_all = ev.get_Qn(pid, n);
    const auto& Q0_all = ev.get_Q0(pid);

    std::complex<double> Qn_ref = {0.0, 0.0};
    double Mq = 0.0;
    for (size_t i = 0; i < Qn_all.size(); ++i) {
        for (size_t j = 0; j < Qn_all[i].size(); ++j) {
            double eta = ev.eta_centers[i];
            double pt = ev.pt_centers[j];
            if (pt < cut.pt_min || pt > cut.pt_max) continue;
            if (std::abs(eta) > cut.eta_cut) continue;
            Qn_ref += Qn_all[i][j];
            Mq += Q0_all[i][j];
        }
    }

    double c2_ref = 0.0;
    if (Mq > 1.0)
        c2_ref = (std::norm(Qn_ref) - Mq) / (Mq * (Mq - 1.0));

    std::vector<std::pair<double, double>> result;
    for (size_t i = 0; i < Qn_eta.size(); ++i) {
        double eta = ev.eta_centers[i];
        if (std::abs(eta) > cut.eta_cut) continue;

        std::complex<double> Pn = Qn_eta[i];
        double Mp = Q0_eta[i];
        double Mpq = Q0_eta[i]; // perfect overlap assumed

        double denom = Mp * (Mq - Mpq);
        double two_prime = 0.0;
        if (denom > 0.0)
            two_prime = (std::real(Pn * std::conj(Qn_ref)) - Mpq) / denom;

        double vn2_val = (c2_ref > 0.0) ? two_prime / std::sqrt(c2_ref) : 0.0;
        result.emplace_back(vn2_val, Mq * (Mq - 1.0));
    }
    return result;
}

#endif
