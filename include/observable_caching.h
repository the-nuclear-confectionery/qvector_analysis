#ifndef OBSERVABLE_CACHE_H
#define OBSERVABLE_CACHE_H

#include <string>
#include <unordered_map>
#include <optional>
#include <functional>
#include <iostream>

#include "event.h"
#include "cuts.h"

// Signature of any observable function
using ObservableFunction = std::function<double(event&, cuts&, ObservableCache&, int)>;

class ObservableCache {
public:
    ObservableCache();

    // Store computed value
    void store_observable(const std::string& key, double value);

    // Retrieve or compute observable
    std::optional<double> get_observable(const std::string& key,
                                         event& ev,
                                         cuts& cut,
                                         ObservableCache& cache,
                                         int pid);

private:
    std::unordered_map<std::string, double> cache;
    std::unordered_map<std::string, ObservableFunction> observable_functions;
};

#endif  // OBSERVABLE_CACHE_H
