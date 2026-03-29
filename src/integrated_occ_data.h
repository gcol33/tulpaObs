// integrated_occ_data.h
// Response data for integrated (multi-source) occupancy models.
// Shared occupancy, per-source detection.

#ifndef TULPAOCC_INTEGRATED_OCC_DATA_H
#define TULPAOCC_INTEGRATED_OCC_DATA_H

#include <vector>

namespace tulpaOcc {

struct IntegratedOccResponseData {
    int n_sites;
    int n_sources;

    // Per-source detection histories and metadata
    // y_flat[source]: flattened y for source s, length n_sites_s * max_visits_s
    std::vector<std::vector<int>> y;
    std::vector<int> max_visits;        // [n_sources]
    std::vector<int> n_sites_per;       // [n_sources] sites observed by each source
    std::vector<std::vector<int>> n_visits;       // [source][site]
    std::vector<std::vector<bool>> any_detected;  // [source][site]

    // site_map[source][local_idx] = global site index (0-based)
    std::vector<std::vector<int>> site_map;
};

} // namespace tulpaOcc

#endif // TULPAOCC_INTEGRATED_OCC_DATA_H
