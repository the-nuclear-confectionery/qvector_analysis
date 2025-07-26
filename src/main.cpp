#include event.h

#include <iostream>
#include <functional>
#include <numeric>
#include "event.h"
#include "observable_caching.h" 



int main(int argc, char** argv) {
    if (argc < 2) return 1;

    event ev(argv[1]);

    int pid = 0;     // charged
    int n = 2;       // v2
    double eta_min = -1.0, eta_max = 1.0;
    double pt_min = 0.2, pt_max = 2.0;

    // Just initialize and test method calls (no output)
    auto qn      = ev.get_Qn(pid, n);
    auto qn_eta  = ev.get_Qn_integrated_over_pt(pid, n, eta_min, eta_max, pt_min, pt_max);
    auto qn_pt   = ev.get_Qn_integrated_over_eta(pid, n, eta_min, eta_max, pt_min, pt_max);
    auto qn_int  = ev.get_Qn_integrated_over_eta_and_pt(pid, n, eta_min, eta_max, pt_min, pt_max);

    auto q0      = ev.get_Q0(pid);
    auto q0_eta  = ev.get_Q0_integrated_over_pt(pid, eta_min, eta_max, pt_min, pt_max);
    auto q0_pt   = ev.get_Q0_integrated_over_eta(pid, eta_min, eta_max, pt_min, pt_max);
    auto q0_int  = ev.get_Q0_integrated_over_eta_and_pt(pid, eta_min, eta_max, pt_min, pt_max);


    // Test observable caching
    ObservableCache cache({"two_particle_correlation"}, 
                          {two_particle_correlation});

    
    return 0;
}




double two_particle_correlation(const event& ev, const cuts& cut, int pid){
    int n = 2; // v2

    // extract cuts
    double eta_min = cut.eta_min;
    double eta_max = cut.eta_max;
    double pt_min = cut.pt_min;
    double pt_max = cut.pt_max;

    std::complex<double> Qn = ev.get_Qn_integrated_over_eta_and_pt(pid, n, eta_min, eta_max, pt_min, pt_max);
    double Q0 = ev.get_Q0_integrated_over_eta_and_pt(pid, eta_min, eta_max, pt_min, pt_max);

    return (std::norm(Qn) - Q0) / (Q0 * (Q0 - 1.0));
}

double four_particle_correlation(const event& ev, int pid, int n,
    double eta_min, double eta_max, double pt_min, double pt_max) 
{
    using namespace std;

    std::complex<double> Qn = ev.get_Qn_integrated_over_eta_and_pt(pid, n, eta_min, eta_max, pt_min, pt_max);
    std::complex<double> Q2n = ev.get_Qn_integrated_over_eta_and_pt(pid, 2 * n, eta_min, eta_max, pt_min, pt_max);
    double M = ev.get_Q0_integrated_over_eta_and_pt(pid, eta_min, eta_max, pt_min, pt_max);

    if (M < 4.0) return 0.0; // Avoid division by zero or negative denominator

    double Qn2 = norm(Qn);               // |Qn|^2
    double Qn4 = Qn2 * Qn2;              // |Qn|^4
    double Q2n2 = norm(Q2n);             // |Q_{2n}|^2
    double Re_cross = real(Q2n * conj(Qn) * conj(Qn)); // Re[Q_{2n} (Qn^*)^2]

    double numerator = Qn4 + Q2n2 - 2.0 * Re_cross 
                     - 2.0 * (2.0 * (M - 2.0) * Qn2 - M * (M - 3.0));

    double denominator = M * (M - 1.0) * (M - 2.0) * (M - 3.0);

    return numerator / denominator;
}


double vn_two_particle(const event& ev, int pid, int n,
    double eta_min, double eta_max, double pt_min, double pt_max) {

    std::complex<double> Qn = ev.get_Qn_integrated_over_eta_and_pt(pid, n, eta_min, eta_max, pt_min, pt_max);
    double Q0 = ev.get_Q0_integrated_over_eta_and_pt(pid, eta_min, eta_max, pt_min, pt_max);

    
}













// Scalar average
std::vector<double> integrated_observable( const event& ev, int pid, 
    double eta_min, double eta_max, double pt_min, double pt_max,
    std::vector<std::function<double(const event& ev, int pid, double eta_min, double eta_max, double pt_min, double pt_max)>> funcs) {

    std::vector<double> results;
    for (const auto& func : funcs) {
        results.push_back(func(ev, pid, eta_min, eta_max, pt_min, pt_max));
    }
    return results;

}



double weighted_average(const std::vector<double>& values, const std::vector<double>& weights) {
    if (values.size() != weights.size() || values.empty()) return 0.0;

    double weighted_sum = 0.0;
    double total_weight = 0.0;

    for (size_t i = 0; i < values.size(); ++i) {
        weighted_sum += values[i] * weights[i];
        total_weight += weights[i];
    }

    return total_weight > 0.0 ? weighted_sum / total_weight : 0.0;
}