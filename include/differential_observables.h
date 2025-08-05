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
    std::vector<double> dN = ev.get_Q0_integrated_over_pt(pid, -cut.eta_cut, cut.eta_cut, cut.pt_min, cut.pt_max);
    double nsamples = static_cast<double>(ev.nsamples);
    std::vector<std::pair<double, double>> result;
    result.reserve(ev.eta_centers.size());
    for (size_t i = 0; i < ev.eta_centers.size(); ++i) {
        double deta = ev.eta_centers[i + 1] - ev.eta_centers[i];
        if (std::abs(ev.eta_centers[i]) > cut.eta_cut)
            result.emplace_back(0.0, 0.0);
        else
            result.emplace_back(dN[i]/deta / nsamples, 1.0);
    }
    return result;
}

// === dN/2pi pT dpT dy ===
inline std::vector<std::pair<double, double>> dN_2piptdpTdy(event& ev, cuts& cut, int pid, int dummy = 0) {
    std::vector<double> dNdpt = ev.get_Q0_integrated_over_eta(pid, -cut.eta_cut, cut.eta_cut, cut.pt_min, cut.pt_max);
    double nsamples = static_cast<double>(ev.nsamples);
    double dy = 2.0 * cut.eta_cut;
    std::vector<std::pair<double, double>> result;
    result.reserve(ev.pt_centers.size());
    for (size_t j = 0; j < ev.pt_centers.size(); ++j) {
        double pt = ev.pt_centers[j];
        double dpt = ev.pt_centers[j + 1] - ev.pt_centers[j]; // assuming pt_centers is sorted
        if (pt < cut.pt_min || pt > cut.pt_max)
            result.emplace_back(0.0, 0.0);
        else
            result.emplace_back(dNdpt[j] / (2 * M_PI * pt * dy * dpt) / nsamples, 1.0);
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
    double M = 0.0;
    for (size_t i = 0; i < ev.eta_centers.size(); ++i) {
        if (std::abs(ev.eta_centers[i]) > cut.eta_cut) continue;
        for (size_t j = 0; j < ev.pt_centers.size(); ++j) {
            double pt = ev.pt_centers[j];
            if (pt < cut.pt_min || pt > cut.pt_max) continue;
            Qn_ref += Qn_all[i][j];
            M += Q0_all[i][j];
        }
    }

    double c2_ref = (M > 1.0) ? (std::norm(Qn_ref) - M) / (M * (M - 1.0)) : 0.0;

    std::vector<std::pair<double, double>> result;
    result.reserve(ev.pt_centers.size());
    for (size_t j = 0; j < ev.pt_centers.size(); ++j) {
        double pt = ev.pt_centers[j];
        if (pt < cut.pt_min || pt > cut.pt_max) {
            result.emplace_back(0.0, 0.0);
            continue;
        }

        std::complex<double> Pn = Qn_pt[j];
        double mp = Q0_pt[j];
        double mq = mp;

        double denom = mp*M - mq;
        double two_prime = (denom > 0.0) ? (std::real(Pn * std::conj(Qn_ref)) - mq) / denom : 0.0;
        double vn2_val = (c2_ref > 0.0) ? two_prime / std::sqrt(c2_ref) : 0.0;
        result.emplace_back(vn2_val, denom);
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
    double M = 0.0;
    for (size_t i = 0; i < ev.eta_centers.size(); ++i) {
        if (std::abs(ev.eta_centers[i]) > cut.eta_cut) continue;
        for (size_t j = 0; j < ev.pt_centers.size(); ++j) {
            double pt = ev.pt_centers[j];
            if (pt < cut.pt_min || pt > cut.pt_max) continue;
            Qn_ref += Qn_all[i][j];
            M += Q0_all[i][j];
        }
    }

    double c2_ref = (M > 1.0) ? (std::norm(Qn_ref) - M) / (M * (M - 1.0)) : 0.0;

    std::vector<std::pair<double, double>> result;
    result.reserve(ev.eta_centers.size());
    for (size_t i = 0; i < ev.eta_centers.size(); ++i) {
        if (std::abs(ev.eta_centers[i]) > cut.eta_cut) {
            result.emplace_back(0.0, 0.0);
            continue;
        }

        std::complex<double> Pn = Qn_eta[i];
        double mp = Q0_eta[i];
        double mq = mp;

        double denom = mp*M - mq;
        double two_prime = (denom > 0.0) ? (std::real(Pn * std::conj(Qn_ref)) - mq) / denom : 0.0;
        double vn2_val = (c2_ref > 0.0) ? two_prime / std::sqrt(c2_ref) : 0.0;
        result.emplace_back(vn2_val, denom);
    }
    return result;
}

#endif
