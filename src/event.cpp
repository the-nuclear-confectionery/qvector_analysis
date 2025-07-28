#include "event.h"
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <algorithm>

event::event(const std::string& filename, config& cfg) {
    read_file(filename, cfg);
}

void event::read_file(const std::string& filename, config& cfg) {
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Failed to open file: " << filename << "\n";
        return;
    }

    TString base = gSystem->BaseName(filename.c_str());
    base.ReplaceAll(".root", "");
    event_id = base.Data();

    TH1D* hSample = dynamic_cast<TH1D*>(file->Get("hSampleCounter"));
    nsamples = hSample ? static_cast<int>(hSample->GetBinContent(1)) : 1;

    detected_n_max = 0;

    if (cfg.calculate_charged) {
        TString test = "ReQ_charged_n0";
        if (file->GetListOfKeys()->Contains(test)) {
            int local_n = read_qvectors(0);
            detected_n_max = std::max(detected_n_max, local_n);
            available_pids.push_back(0);
        } else {
            std::cerr << "Warning: Requested charged particles but Q_charged not found in file.\n";
        }
    }

    for (int pid : cfg.pids) {
        if (pid == 0) continue;
        std::string prefix = "Q_pid" + std::to_string(pid);
        TString test = TString::Format("Re%s_n0", prefix.c_str());

        if (file->GetListOfKeys()->Contains(test)) {
            int local_n = read_qvectors(pid);
            detected_n_max = std::max(detected_n_max, local_n);
            available_pids.push_back(pid);
        } else {
            std::cerr << "Warning: Requested PID " << pid << " not found in file.\n";
        }
    }

    file->Close();

    if (cfg.max_n > detected_n_max) {
        std::cerr << "Warning: Requested max_n = " << cfg.max_n
                  << " but only n ≤ " << detected_n_max << " found in file.\n";
    }

    bool has_vn4 = std::find(cfg.integrated_observables.begin(),
                             cfg.integrated_observables.end(), "vn{4}") != cfg.integrated_observables.end();

    if (has_vn4 && detected_n_max < 2 * cfg.max_n) {
        std::cerr << "Warning: vn{4} requires Q-vectors up to 2×n = " << 2 * cfg.max_n
                  << " but only n ≤ " << detected_n_max << " found in file. vn{4} will be disabled.\n";

        cfg.integrated_observables.erase(
            std::remove(cfg.integrated_observables.begin(),
                        cfg.integrated_observables.end(), "vn{4}"),
            cfg.integrated_observables.end()
        );
    }
}

int event::read_qvectors(int pid) {
    std::string prefix = (pid == 0) ? "Q_charged" : ("Q_pid" + std::to_string(pid));
    int found_n_max = -1;

    std::vector<std::vector<std::vector<std::complex<double>>>> qn_list;
    std::vector<std::vector<double>> q0;

    for (int n = 0; ; ++n) {
        TString name_re = TString::Format("Re%s_n%d", prefix.c_str(), n);
        TString name_im = TString::Format("Im%s_n%d", prefix.c_str(), n);

        TH2D* hre = dynamic_cast<TH2D*>(gDirectory->Get(name_re));
        TH2D* him = dynamic_cast<TH2D*>(gDirectory->Get(name_im));
        if (!hre || !him) break;

        found_n_max = n;

        if (!binning_initialized) {
            init_binning(hre);
            binning_initialized = true;
        }

        int nx = hre->GetNbinsX();
        int ny = hre->GetNbinsY();
        std::vector<std::vector<std::complex<double>>> qn(nx, std::vector<std::complex<double>>(ny));

        for (int ix = 0; ix < nx; ++ix)
            for (int iy = 0; iy < ny; ++iy)
                qn[ix][iy] = {hre->GetBinContent(ix + 1, iy + 1), him->GetBinContent(ix + 1, iy + 1)};
        
        qn_list.push_back(std::move(qn));

        if (n == 0) {
            std::vector<std::vector<double>> m0(nx, std::vector<double>(ny));
            for (int ix = 0; ix < nx; ++ix)
                for (int iy = 0; iy < ny; ++iy)
                    m0[ix][iy] = hre->GetBinContent(ix + 1, iy + 1);
            q0 = std::move(m0);
        }
    }

    Qn[pid] = std::move(qn_list);
    Q0[pid] = std::move(q0);
    return found_n_max;
}

void event::init_binning(const TH2D* hist) {
    int nx = hist->GetNbinsX();
    int ny = hist->GetNbinsY();
    eta_centers.resize(nx);
    pt_centers.resize(ny);

    for (int ix = 1; ix <= nx; ++ix)
        eta_centers[ix - 1] = hist->GetXaxis()->GetBinCenter(ix);
    for (int iy = 1; iy <= ny; ++iy)
        pt_centers[iy - 1] = hist->GetYaxis()->GetBinCenter(iy);
}

std::vector<std::vector<std::complex<double>>> event::get_Qn(int pid, int n) const {
    return Qn.at(pid).at(n);
}

std::complex<double> event::get_Qn_integrated_over_eta_and_pt(int pid, int n,
    double eta_min, double eta_max, double pt_min, double pt_max) const {
    const auto& q = Qn.at(pid).at(n);
    std::complex<double> result = {0.0, 0.0};
    for (size_t ix = 0; ix < eta_centers.size(); ++ix)
        if (eta_centers[ix] >= eta_min && eta_centers[ix] <= eta_max)
            for (size_t iy = 0; iy < pt_centers.size(); ++iy)
                if (pt_centers[iy] >= pt_min && pt_centers[iy] <= pt_max)
                    result += q[ix][iy];
    return result;
}

std::vector<std::complex<double>> event::get_Qn_integrated_over_eta(int pid, int n,
    double eta_min, double eta_max, double pt_min, double pt_max) const {
    const auto& q = Qn.at(pid).at(n);
    std::vector<std::complex<double>> result(pt_centers.size(), {0.0, 0.0});
    for (size_t iy = 0; iy < pt_centers.size(); ++iy)
        if (pt_centers[iy] >= pt_min && pt_centers[iy] <= pt_max)
            for (size_t ix = 0; ix < eta_centers.size(); ++ix)
                if (eta_centers[ix] >= eta_min && eta_centers[ix] <= eta_max)
                    result[iy] += q[ix][iy];
    return result;
}

std::vector<std::complex<double>> event::get_Qn_integrated_over_pt(int pid, int n,
    double eta_min, double eta_max, double pt_min, double pt_max) const {
    const auto& q = Qn.at(pid).at(n);
    std::vector<std::complex<double>> result(eta_centers.size(), {0.0, 0.0});
    for (size_t ix = 0; ix < eta_centers.size(); ++ix)
        if (eta_centers[ix] >= eta_min && eta_centers[ix] <= eta_max)
            for (size_t iy = 0; iy < pt_centers.size(); ++iy)
                if (pt_centers[iy] >= pt_min && pt_centers[iy] <= pt_max)
                    result[ix] += q[ix][iy];
    return result;
}

std::vector<std::vector<double>> event::get_Q0(int pid) const {
    return Q0.at(pid);
}

double event::get_Q0_integrated_over_eta_and_pt(int pid,
    double eta_min, double eta_max, double pt_min, double pt_max) const {
    const auto& q0 = Q0.at(pid);
    double result = 0.0;
    for (size_t ix = 0; ix < eta_centers.size(); ++ix)
        if (eta_centers[ix] >= eta_min && eta_centers[ix] <= eta_max)
            for (size_t iy = 0; iy < pt_centers.size(); ++iy)
                if (pt_centers[iy] >= pt_min && pt_centers[iy] <= pt_max)
                    result += q0[ix][iy];
    return result;
}

std::vector<double> event::get_Q0_integrated_over_eta(int pid,
    double eta_min, double eta_max, double pt_min, double pt_max) const {
    const auto& q0 = Q0.at(pid);
    std::vector<double> result(pt_centers.size(), 0.0);
    for (size_t iy = 0; iy < pt_centers.size(); ++iy)
        if (pt_centers[iy] >= pt_min && pt_centers[iy] <= pt_max)
            for (size_t ix = 0; ix < eta_centers.size(); ++ix)
                if (eta_centers[ix] >= eta_min && eta_centers[ix] <= eta_max)
                    result[iy] += q0[ix][iy];
    return result;
}

std::vector<double> event::get_Q0_integrated_over_pt(int pid,
    double eta_min, double eta_max, double pt_min, double pt_max) const {
    const auto& q0 = Q0.at(pid);
    std::vector<double> result(eta_centers.size(), 0.0);
    for (size_t ix = 0; ix < eta_centers.size(); ++ix)
        if (eta_centers[ix] >= eta_min && eta_centers[ix] <= eta_max)
            for (size_t iy = 0; iy < pt_centers.size(); ++iy)
                if (pt_centers[iy] >= pt_min && pt_centers[iy] <= pt_max)
                    result[ix] += q0[ix][iy];
    return result;
}
