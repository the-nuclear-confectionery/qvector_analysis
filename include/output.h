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

    void write_all() const;
    void write_integrated(int pid, const std::string& outdir) const;
    void write_all_differential(int pid, const std::string& outdir) const;

    std::string build_output_path() const;
private:
    const config& cfg_;
    const Observables& obs_;
    const std::vector<event>& events_;

};

#endif // OUTPUT_H
