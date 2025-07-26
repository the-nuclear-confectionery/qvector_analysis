#include "ObservableCache.h"


// === Constructor ===
ObservableCache::ObservableCache() {
    // Add your function bindings here (must be declared somewhere)
    observable_functions["two_particle_correlation"] = two_particle_correlation;
    observable_functions["vn_two_particle"]          = vn_two_particle;
    observable_functions["four_particle_correlation"] = four_particle_correlation;
}

// === Store ===
void ObservableCache::store_observable(const std::string& key, double value) {
    cache[key] = value;
}

// === Get or Compute ===
std::optional<double> ObservableCache::get_observable(const std::string& key, event& ev, cuts& cut, int pid) {
    auto cached = cache.find(key);
    if (cached != cache.end())
        return cached->second;

    auto func_it = observable_functions.find(key);
    if (func_it != observable_functions.end()) {
        double value = func_it->second(ev, cut, pid);
        cache.emplace(key, value);
        return value;
    }

    std::cerr << "Observable function not found for key: " << key << std::endl;
    return std::nullopt;
}
