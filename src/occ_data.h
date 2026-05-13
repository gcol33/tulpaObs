// occ_data.h
// Model-specific response data for occupancy models
// Opaque to tulpa — passed via ModelData::model_response_data

#ifndef TULPAOCC_OCC_DATA_H
#define TULPAOCC_OCC_DATA_H

#include <vector>

namespace tulpaObs {

// Single-season occupancy response data
struct OccResponseData {
    int n_sites;
    int max_visits;

    // Detection history: y[site * max_visits + visit], -1 = missing (no visit)
    std::vector<int> y;

    // Number of actual visits per site
    std::vector<int> n_visits;

    // Visit-level detection covariates (optional)
    // X_det_visit[site * max_visits * p_det_visit + visit * p_det_visit + j]
    std::vector<double> X_det_visit;
    int p_det_visit = 0;

    // Whether site had at least one detection
    std::vector<bool> any_detected;

    // Precomputed: number of detections per site
    std::vector<int> n_detections;
};

} // namespace tulpaObs

#endif // TULPAOCC_OCC_DATA_H
