// dyn_occ_data.h
// Model-specific response data for dynamic (multi-season) occupancy models
// MacKenzie et al. (2003) colonization-extinction dynamics

#ifndef TULPAOCC_DYN_OCC_DATA_H
#define TULPAOCC_DYN_OCC_DATA_H

#include <vector>

namespace tulpaObs {

// Multi-season occupancy response data
// N observations = n_sites (likelihood loops over seasons internally)
struct DynOccResponseData {
    int n_sites;
    int n_seasons;
    int max_visits;

    // Detection history: y[site * n_seasons * max_visits + season * max_visits + visit]
    // -1 = missing (no visit)
    std::vector<int> y;

    // Number of actual visits per site-season: n_visits[site * n_seasons + season]
    std::vector<int> n_visits;

    // Whether any detection occurred in each site-season
    std::vector<bool> any_detected;  // [n_sites * n_seasons]
};

} // namespace tulpaObs

#endif // TULPAOCC_DYN_OCC_DATA_H
