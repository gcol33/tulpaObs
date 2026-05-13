// dyn_occ_likelihood.h
// Multi-season dynamic occupancy likelihood
//
// Four processes:
//   Process 0: initial occupancy — eta[0] = logit(psi1_i)
//   Process 1: detection — eta[1] = logit(p_i), site-level
//   Process 2: colonization — eta[2] = logit(gamma_i)
//   Process 3: extinction — eta[3] = logit(epsilon_i)
//
// HMM forward algorithm per site:
//   alpha[0] = psi1
//   alpha[t] = alpha[t-1]*(1-epsilon) + (1-alpha[t-1])*gamma  for t >= 1
//   L_site = prod_t [alpha[t]*P(y_t|1) + (1-alpha[t])*P(y_t|0)]
//
// Reference: MacKenzie et al. (2003) Ecology

#ifndef TULPAOCC_DYN_OCC_LIKELIHOOD_H
#define TULPAOCC_DYN_OCC_LIKELIHOOD_H

#include <cmath>
#include <vector>
#include <type_traits>
#include "dyn_occ_data.h"
#include "occ_likelihood.h"  // For safe_log1pexp, log_inv_logit, log1m_inv_logit, log_sum_exp
#include <tulpa/likelihood.h>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/autodiff_arena.h>
#include <tulpa/autodiff_fwd.h>

namespace TulpaObs {

// ============================================================================
// Dynamic occupancy log-likelihood (per site)
//
// N = n_sites. Each call computes the full HMM forward probability for site i
// across all seasons, using process-level linear predictors.
//
// eta[0] = logit(psi1), eta[1] = logit(p), eta[2] = logit(gamma), eta[3] = logit(epsilon)
// ============================================================================
template<typename T>
T dyn_occ_log_likelihood(
    int i,                            // Site index
    const T* eta,                     // eta[0..3] for this site
    const T& logit_zi,                // Unused
    const T& logit_oi,                // Unused
    const std::vector<T>& params,     // Full parameter vector
    const tulpa::ModelData& data,
    const tulpa::ParamLayout& layout,
    const void* model_data
) {
    // Cross-DLL arena sync (Windows/MinGW)
    if constexpr (std::is_same_v<T, tulpa::arena::Var>) {
        tulpa::arena::current_arena() = eta[0].arena_;
    }

    const auto* dyn = static_cast<const DynOccResponseData*>(model_data);
    const int T_seasons = dyn->n_seasons;
    const int K = dyn->max_visits;

    // Process linear predictors (constant across seasons for this site)
    T logit_psi1 = eta[0];
    T logit_p    = eta[1];
    T logit_gam  = eta[2];
    T logit_eps  = eta[3];

    // Pre-compute log probabilities
    T log_psi1     = log_inv_logit(logit_psi1);
    T log1m_psi1   = log1m_inv_logit(logit_psi1);
    T log_gam      = log_inv_logit(logit_gam);
    T log1m_gam    = log1m_inv_logit(logit_gam);
    T log_eps      = log_inv_logit(logit_eps);
    T log1m_eps    = log1m_inv_logit(logit_eps);

    // HMM forward algorithm in log space
    // log_alpha[0] = log(P(occupied at t=0))
    // log_alpha[1] = log(P(unoccupied at t=0))
    T log_alpha_occ = log_psi1;
    T log_alpha_unocc = log1m_psi1;

    T log_lik = T(0.0);

    for (int t = 0; t < T_seasons; t++) {
        // Compute detection contribution for this season
        int base_y = i * T_seasons * K + t * K;
        int base_nv = i * T_seasons + t;
        int nv = dyn->n_visits[base_nv];
        bool detected_this_season = dyn->any_detected[base_nv];

        T log_p_data_occ = T(0.0);    // log P(data_t | occupied)

        if (nv > 0) {
            T sum_log1m_p = T(0.0);
            for (int j = 0; j < nv; j++) {
                int y_ij = dyn->y[base_y + j];
                if (y_ij < 0) continue;

                T logit_p_ij = logit_p;
                // Add visit-level covariates if present
                if (dyn->p_det_visit > 0) {
                    int xbase = i * T_seasons * K * dyn->p_det_visit
                              + t * K * dyn->p_det_visit
                              + j * dyn->p_det_visit;
                    int beta_offset = layout.extra_offset;
                    for (int c = 0; c < dyn->p_det_visit; c++) {
                        logit_p_ij = logit_p_ij + T(dyn->X_det_visit[xbase + c]) * params[beta_offset + c];
                    }
                }

                T log_p = log_inv_logit(logit_p_ij);
                T log1m_p = log1m_inv_logit(logit_p_ij);

                if (y_ij == 1) {
                    log_p_data_occ = log_p_data_occ + log_p;
                } else {
                    log_p_data_occ = log_p_data_occ + log1m_p;
                }
                sum_log1m_p = sum_log1m_p + log1m_p;
            }

            if (detected_this_season) {
                // P(data | unoccupied) = 0 → log = -inf
                // Contribution: log(alpha_occ * P(data|occ))
                log_lik = log_lik + log_alpha_occ + log_p_data_occ;

                // After observing detection, we know z_t = 1
                // Update alpha for transition
                log_alpha_occ = T(0.0);       // log(1)
                log_alpha_unocc = T(-1e10);   // log(0)
            } else {
                // P(data | occ) = prod(1-p), P(data | unocc) = 1
                // Contribution: log(alpha_occ * prod(1-p) + alpha_unocc)
                T term1 = log_alpha_occ + sum_log1m_p;
                T term2 = log_alpha_unocc;
                log_lik = log_lik + log_sum_exp(term1, term2);

                // Posterior P(z_t = 1 | data_t) for transition
                // P(z=1|data) = alpha_occ * prod(1-p) / [alpha_occ * prod(1-p) + alpha_unocc]
                T log_norm = log_sum_exp(term1, term2);
                log_alpha_occ = term1 - log_norm;
                log_alpha_unocc = term2 - log_norm;
            }
        }
        // If nv == 0, no data this season, alpha unchanged, no lik contribution

        // Transition to next season (if not last)
        if (t < T_seasons - 1) {
            // P(z_{t+1}=1) = P(z_t=1)*(1-eps) + P(z_t=0)*gam
            T new_log_occ = log_sum_exp(
                log_alpha_occ + log1m_eps,
                log_alpha_unocc + log_gam
            );
            T new_log_unocc = log_sum_exp(
                log_alpha_occ + log_eps,
                log_alpha_unocc + log1m_gam
            );
            log_alpha_occ = new_log_occ;
            log_alpha_unocc = new_log_unocc;
        }
    }

    return log_lik;
}

} // namespace TulpaObs

#endif // TULPAOCC_DYN_OCC_LIKELIHOOD_H
