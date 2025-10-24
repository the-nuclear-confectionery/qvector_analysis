#ifndef INTEGRATED_OBSERVABLES_H
#define INTEGRATED_OBSERVABLES_H

#include "event.h"
#include "cuts.h"
#include <cmath>
#include <complex>
#include <utility>

// === Multiplicity: M ===
inline std::pair<double, double> M(event& ev, cuts& cut, int pid, int dummy = 0) {
    double mult = 0.0;
    const auto& Q0 = ev.get_Q0(pid);
    // event number of samples converted to double 
    double nsamples = static_cast<double>(ev.nsamples);


    for (size_t i = 0; i < ev.eta_centers.size(); ++i) {\
        if (std::abs(ev.eta_centers[i]) > cut.eta_cut) continue;
        for (size_t j = 0; j < ev.pt_centers.size(); ++j) {
            if (ev.pt_centers[j] < cut.pt_min || ev.pt_centers[j] > cut.pt_max) continue;
            mult += Q0[i][j]/ nsamples; // normalize by number of samples
        }
    }
    return {mult, 1.0};
}

// === Mean pT: mean_pT ===
inline std::pair<double, double> mean_pT(event& ev, cuts& cut, int pid,  int dummy = 0) {
    const auto& Q0 = ev.get_Q0(pid);
    const auto& pt = ev.pt_centers;

    double numerator = 0.0;
    double denominator = 0.0;

    for (size_t i = 0; i < ev.eta_centers.size(); ++i) {
        if (std::abs(ev.eta_centers[i]) > cut.eta_cut) continue;
        for (size_t j = 0; j < pt.size(); ++j) {
            if (pt[j] < cut.pt_min || pt[j] > cut.pt_max) continue;
            double w = Q0[i][j];
            numerator += pt[j] * w;
            denominator += w;
        }
    }

    double result = (denominator > 0.0) ? numerator / denominator : 0.0;
    return {result, denominator};
}

// === v_n{2}: vn{2} ===
inline std::pair<double, double> vn2(event& ev, cuts& cut, int pid, int n) {
    std::complex<double> Qn = {0.0, 0.0};
    double Mval = 0.0;

    const auto& Qvec = ev.get_Qn(pid, n);
    const auto& Q0 = ev.get_Q0(pid);

    for (size_t i = 0; i < ev.eta_centers.size(); ++i) {
        if (std::abs(ev.eta_centers[i]) > cut.eta_cut) continue;
        for (size_t j = 0; j < ev.pt_centers.size(); ++j) {
            if (ev.pt_centers[j] < cut.pt_min || ev.pt_centers[j] > cut.pt_max) continue;
            Qn += Qvec[i][j];
            Mval += Q0[i][j];
        }
    }

    double weight = Mval * (Mval - 1);
    double c2 = (weight > 0.0) ? (std::norm(Qn) - Mval) / weight : 0.0;
    double v2 = sqrt(c2);
    return {v2, weight};
}

// === v_n{EP}: vn{EP} ===
inline std::pair<double, double> vnEP(event& ev, cuts& cut, int pid, int n) {
    std::complex<double> Qn = {0.0, 0.0};
    double Mval = 0.0;

    const auto& Qvec = ev.get_Qn(pid, n);
    const auto& Q0 = ev.get_Q0(pid);

    for (size_t i = 0; i < ev.eta_centers.size(); ++i) {
        if (std::abs(ev.eta_centers[i]) > cut.eta_cut) continue;
        for (size_t j = 0; j < ev.pt_centers.size(); ++j) {
            if (ev.pt_centers[j] < cut.pt_min || ev.pt_centers[j] > cut.pt_max) continue;
            Qn += Qvec[i][j];
            Mval += Q0[i][j];
        }
    }

    std::complex<double> qn_norm = (Mval > 0.0) ? Qn / Mval : 0.0;
    double psi_n = std::arg(qn_norm) / n;

    double result = 0.0;
    if (Mval > 0.0)
        result = std::real(qn_norm * std::exp(std::complex<double>(0, -n * psi_n)));

    return {result, Mval};
}

// === v_n{4}: vn{4} ===
inline std::pair<double, double> vn4(event& ev, cuts& cut, int pid, int n) {
    std::complex<double> Qn = {0.0, 0.0};
    std::complex<double> Q2n = {0.0, 0.0};
    double Mval = 0.0;

    const auto& Qn_vec = ev.get_Qn(pid, n);
    const auto& Q2n_vec = ev.get_Qn(pid, 2 * n);
    const auto& Q0 = ev.get_Q0(pid);

    for (size_t i = 0; i < ev.eta_centers.size(); ++i) {
        if (std::abs(ev.eta_centers[i]) > cut.eta_cut) continue;
        for (size_t j = 0; j < ev.pt_centers.size(); ++j) {
            if (ev.pt_centers[j] < cut.pt_min || ev.pt_centers[j] > cut.pt_max) continue;
            Qn += Qn_vec[i][j];
            Q2n += Q2n_vec[i][j];
            Mval += Q0[i][j];
        }
    }

    double w = Mval * (Mval - 1) * (Mval - 2) * (Mval - 3);
    double result = 0.0;
    if (w > 0.0) {
        double Qn2 = std::norm(Qn);
        double Qn4 = Qn2 * Qn2;
        double Q2n2 = std::norm(Q2n);
        double Re_term = std::real(Q2n * std::pow(std::conj(Qn), 2));

        result = (Qn4 + Q2n2 - 2.0 * Re_term
                 - 2.0 * (2.0 * (Mval - 2) * Qn2 - Mval * (Mval - 3))) / w;
    }
    double v4 = sqrt(result)
    return {result, w};
}

#endif
