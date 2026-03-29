// integrated_occ_likelihood.h
// Integrated (multi-source) occupancy likelihood.
//
// n_processes = 1 + n_sources:
//   Process 0: occupancy — eta[0] = logit(psi_i)
//   Process 1..S: detection per source — eta[s] = logit(p_is)
//
// Shared occupancy: psi_i is common across all sources.
// Per-source detection: each source has its own p.
//
// Likelihood per site i:
//   If any source detected:
//     log(psi_i) + sum_s log P(y_is | p_is, z_i=1)
//   If no source detected:
//     log(psi_i * prod_s P(y_is=all0 | p_is) + (1-psi_i))

#ifndef TULPAOCC_INTEGRATED_OCC_LIKELIHOOD_H
#define TULPAOCC_INTEGRATED_OCC_LIKELIHOOD_H

#include <cmath>
#include <vector>
#include <type_traits>
#include "integrated_occ_data.h"
#include "occ_likelihood.h"  // For log_inv_logit, log1m_inv_logit, log_sum_exp
#include <tulpa/likelihood.h>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/autodiff_arena.h>
#include <tulpa/autodiff_fwd.h>

namespace tulpaOcc {

template<typename T>
T integrated_occ_log_likelihood(
    int i,                            // Site index (global)
    const T* eta,                     // eta[0]=logit(psi), eta[1..S]=logit(p_s)
    const T& logit_zi,
    const T& logit_oi,
    const std::vector<T>& params,
    const tulpa::ModelData& data,
    const tulpa::ParamLayout& layout,
    const void* model_data
) {
    if constexpr (std::is_same_v<T, tulpa::arena::Var>) {
        tulpa::arena::current_arena() = eta[0].arena_;
    }

    const auto* intd = static_cast<const IntegratedOccResponseData*>(model_data);

    T logit_psi = eta[0];
    T log_psi = log_inv_logit(logit_psi);
    T log1m_psi = log1m_inv_logit(logit_psi);

    bool any_source_detected = false;
    T log_det_given_occ = T(0.0);    // sum over sources of log P(y_is | z=1, p_is)
    T sum_log_prob_all0 = T(0.0);    // sum over sources of log P(y_is=all0 | z=1, p_is)

    for (int s = 0; s < intd->n_sources; s++) {
        // Check if this site was observed by source s
        int local_idx = -1;
        for (int j = 0; j < intd->n_sites_per[s]; j++) {
            if (intd->site_map[s][j] == i) {
                local_idx = j;
                break;
            }
        }
        if (local_idx < 0) continue;  // Source s didn't observe site i

        T logit_p = eta[1 + s];
        int K = intd->n_visits[s][local_idx];
        if (K == 0) continue;

        T sum_log1m_p = T(0.0);
        T log_p_data = T(0.0);

        int y_offset = local_idx * intd->max_visits[s];
        for (int j = 0; j < K; j++) {
            int y_ij = intd->y[s][y_offset + j];
            if (y_ij < 0) continue;

            T log_p = log_inv_logit(logit_p);
            T log1m_p = log1m_inv_logit(logit_p);

            if (y_ij == 1) {
                log_p_data = log_p_data + log_p;
                any_source_detected = true;
            } else {
                log_p_data = log_p_data + log1m_p;
            }
            sum_log1m_p = sum_log1m_p + log1m_p;
        }

        if (intd->any_detected[s][local_idx]) {
            log_det_given_occ = log_det_given_occ + log_p_data;
        } else {
            log_det_given_occ = log_det_given_occ + log_p_data;
            sum_log_prob_all0 = sum_log_prob_all0 + sum_log1m_p;
        }
    }

    if (any_source_detected) {
        return log_psi + log_det_given_occ;
    } else {
        T log_term1 = log_psi + sum_log_prob_all0;
        T log_term2 = log1m_psi;
        return log_sum_exp(log_term1, log_term2);
    }
}

} // namespace tulpaOcc

#endif // TULPAOCC_INTEGRATED_OCC_LIKELIHOOD_H
