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

    // Random effects over visit rows on the detection logit. They are extra
    // parameters after the visit covariates, type-blocked the way the engine
    // lays out its own random effects: every term's log SDs, then every
    // correlated term's Cholesky raw values, then every term's latent z.
    struct VisitReTerm {
        int n_groups = 0;
        int n_coefs = 0;
        bool correlated = false;
        double sigma_scale = 1.0;
        // Per unit-major visit row: 0-based group, -1 where the visit is absent.
        std::vector<int> group;
        // Per visit row, n_coefs design values (intercept column first).
        std::vector<double> Z;
        // Absolute parameter positions, set once the layout is known.
        std::vector<int> log_sigma_idx;
        int chol_start = -1;
        int re_start = 0;
    };
    std::vector<VisitReTerm> visit_re;
};

} // namespace tulpaObs

#endif // TULPAOCC_OCC_DATA_H
