// occ_likelihood.h
// Single-season occupancy likelihood for tulpa engine
//
// Two processes:
//   Process 0: occupancy — eta[0] = logit(psi_i), site-level covariates
//   Process 1: detection — eta[1] = logit(p_i), site-level covariates
//
// Visit-level detection covariates (if any) are stored in OccResponseData
// and added to eta[1] inside the likelihood.
//
// Likelihood per site i:
//   If any detection: log(psi_i) + sum_j [y_ij*log(p_ij) + (1-y_ij)*log(1-p_ij)]
//   If all zeros:     log(psi_i * prod_j(1-p_ij) + (1-psi_i))
//
// Reference: MacKenzie et al. (2002) Ecology

#ifndef TULPAOCC_OCC_LIKELIHOOD_H
#define TULPAOCC_OCC_LIKELIHOOD_H

#include <cmath>
#include <vector>
#include <type_traits>
#include "occ_data.h"
#include <tulpa/likelihood.h>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/autodiff_arena.h>
#include <tulpa/autodiff_fwd.h>

namespace tulpaObs {

// AD-safe log(1 + exp(x)) that avoids overflow.
// Uses comparison operators (work for double, arena::Var, fwd::Dual)
// and AD-tracked math ops to preserve gradient flow.
//
// For T = double: std::exp/std::log1p via ADL on double args.
// For T = arena::Var: tulpa::arena::exp/log1p via ADL on Var args.
// For T = fwd::Dual: fwd::exp/log1p via ADL on Dual args.
template<typename T>
inline T safe_log1pexp(const T& x) {
    using std::exp; using std::log1p;  // fallback for T = double
    if (x > 35.0) return x;                   // exp(x) >> 1, so log(1+exp(x)) ~ x
    if (x < -10.0) return exp(x);             // exp(x) << 1, so log(1+exp(x)) ~ exp(x)
    return log1p(exp(x));                      // general case
}

// Safe log(exp(a) + exp(b)) = max + log(1 + exp(min - max))
template<typename T>
inline T log_sum_exp(const T& a, const T& b) {
    T mx = (a >= b) ? a : b;
    T mn = (a >= b) ? b : a;
    return mx + safe_log1pexp(mn - mx);
}

// log(inv_logit(x)) = x - log(1 + exp(x))
template<typename T>
inline T log_inv_logit(const T& x) {
    return x - safe_log1pexp(x);
}

// log(1 - inv_logit(x)) = -log(1 + exp(x))
template<typename T>
inline T log1m_inv_logit(const T& x) {
    return T(0.0) - safe_log1pexp(x);
}

// ============================================================================
// Single-season occupancy log-likelihood (per site)
// ============================================================================
template<typename T>
T occ_log_likelihood(
    int i,                            // Site index
    const T* eta,                     // eta[0]=logit(psi), eta[1]=logit(p) base
    const T& logit_zi,                // Unused
    const T& logit_oi,                // Unused
    const std::vector<T>& params,     // Full parameter vector
    const tulpa::ModelData& data,
    const tulpa::ParamLayout& layout,
    const void* model_data
) {
    // Cross-DLL arena sync: on Windows/MinGW, the thread-local current_arena()
    // is duplicated per DLL. Sync tulpaObs's copy from the incoming Var objects
    // so that T(0.0) constructors inside this function find the correct arena.
    if constexpr (std::is_same_v<T, tulpa::arena::Var>) {
        tulpa::arena::current_arena() = eta[0].arena_;
    }

    const auto* occ = static_cast<const OccResponseData*>(model_data);

    T logit_psi = eta[0];
    T log_psi = log_inv_logit(logit_psi);
    T log1m_psi = log1m_inv_logit(logit_psi);

    int K = occ->n_visits[i];
    if (K == 0) {
        // No visits — contributes nothing (or marginal over psi)
        return T(0.0);
    }

    // Sum log p(y_ij | p_ij) across visits, and accumulate log(1-p_ij) for
    // the all-zero case
    T log_p_det_given_occ = T(0.0);  // sum_j [y_ij*log(p_ij) + (1-y_ij)*log(1-p_ij)]
    T sum_log1m_p = T(0.0);          // sum_j log(1 - p_ij)

    for (int j = 0; j < K; j++) {
        int y_ij = occ->y[i * occ->max_visits + j];
        if (y_ij < 0) continue;  // Missing visit

        // Detection linear predictor for visit j
        T logit_p_ij = eta[1];

        // Add visit-level covariates if present
        if (occ->p_det_visit > 0) {
            int base = i * occ->max_visits * occ->p_det_visit + j * occ->p_det_visit;
            int beta_offset = layout.extra_offset;
            for (int c = 0; c < occ->p_det_visit; c++) {
                logit_p_ij = logit_p_ij + T(occ->X_det_visit[base + c]) * params[beta_offset + c];
            }
        }

        T log_p = log_inv_logit(logit_p_ij);
        T log1m_p = log1m_inv_logit(logit_p_ij);

        if (y_ij == 1) {
            log_p_det_given_occ = log_p_det_given_occ + log_p;
        } else {
            log_p_det_given_occ = log_p_det_given_occ + log1m_p;
        }

        sum_log1m_p = sum_log1m_p + log1m_p;
    }

    if (occ->any_detected[i]) {
        // Site detected: P(data|occupied) * P(occupied)
        // = psi * prod[p^y * (1-p)^(1-y)]
        return log_psi + log_p_det_given_occ;
    } else {
        // Site not detected: P(occupied, not detected) + P(not occupied)
        // = psi * prod(1-p) + (1-psi)
        T log_term1 = log_psi + sum_log1m_p;   // log(psi * prod(1-p))
        T log_term2 = log1m_psi;                // log(1-psi)
        return log_sum_exp(log_term1, log_term2);
    }
}

// ============================================================================
// Residual function for H-mode gradients (d log_lik / d eta[k])
// ============================================================================
inline void occ_residual(
    int i,
    const double* eta,
    double logit_zi,
    double logit_oi,
    const std::vector<double>& params,
    const tulpa::ModelData& data,
    const tulpa::ParamLayout& layout,
    const void* model_data,
    double* resid_out
) {
    const auto* occ = static_cast<const OccResponseData*>(model_data);

    double logit_psi = eta[0];
    double psi = 1.0 / (1.0 + std::exp(-logit_psi));

    int K = occ->n_visits[i];

    // d(ll)/d(logit_psi): requires careful derivation
    // d(ll)/d(logit_p): sum over visits

    if (K == 0) {
        resid_out[0] = 0.0;
        resid_out[1] = 0.0;
        return;
    }

    // Compute detection probabilities and sum_log1m_p
    double prod_1m_p = 1.0;
    double d_logit_p_sum = 0.0;

    for (int j = 0; j < K; j++) {
        int y_ij = occ->y[i * occ->max_visits + j];
        if (y_ij < 0) continue;

        double logit_p_ij = eta[1];
        if (occ->p_det_visit > 0) {
            int base = i * occ->max_visits * occ->p_det_visit + j * occ->p_det_visit;
            int beta_offset = layout.extra_offset;
            for (int c = 0; c < occ->p_det_visit; c++) {
                logit_p_ij += occ->X_det_visit[base + c] * params[beta_offset + c];
            }
        }

        double p_ij = 1.0 / (1.0 + std::exp(-logit_p_ij));
        prod_1m_p *= (1.0 - p_ij);

        // d(ll_det)/d(logit_p_ij) = y_ij - p_ij (for the part that goes through eta[1])
        // Only site-level part (visit-level betas have separate gradient)
        d_logit_p_sum += (y_ij - p_ij);
    }

    if (occ->any_detected[i]) {
        // d log(psi) / d logit_psi = 1 - psi
        resid_out[0] = 1.0 - psi;
        resid_out[1] = d_logit_p_sum;
    } else {
        // ll = log(psi * prod(1-p) + (1-psi))
        double denom = psi * prod_1m_p + (1.0 - psi);

        // d(ll)/d(logit_psi) = psi*(1-psi) * (prod(1-p) - 1) / denom
        resid_out[0] = psi * (1.0 - psi) * (prod_1m_p - 1.0) / denom;

        // d(ll)/d(logit_p) is more complex — via prod(1-p) derivative
        // d prod(1-p)/d logit_p_j = -p_j*(1-p_j) * prod_{k!=j}(1-p_k)
        // For site-level logit_p (all visits share eta[1]):
        // d(ll)/d eta[1] = psi / denom * d(prod(1-p))/d(eta[1])
        // = psi / denom * sum_j [-p_j * prod_{k!=j}(1-p_k)]
        // For uniform p: = psi / denom * K * (-p*(1-p)^{K-1})
        // General case: recalculate per-visit
        double d_prod_d_eta1 = 0.0;
        for (int j = 0; j < K; j++) {
            int y_ij = occ->y[i * occ->max_visits + j];
            if (y_ij < 0) continue;

            double logit_p_ij = eta[1];
            if (occ->p_det_visit > 0) {
                int base = i * occ->max_visits * occ->p_det_visit + j * occ->p_det_visit;
                int beta_offset = layout.extra_offset;
                for (int c = 0; c < occ->p_det_visit; c++) {
                    logit_p_ij += occ->X_det_visit[base + c] * params[beta_offset + c];
                }
            }
            double p_ij = 1.0 / (1.0 + std::exp(-logit_p_ij));
            double prod_others = (1.0 - p_ij) > 1e-300 ? prod_1m_p / (1.0 - p_ij) : 0.0;
            d_prod_d_eta1 += -p_ij * (1.0 - p_ij) * prod_others;
        }

        resid_out[1] = psi * d_prod_d_eta1 / denom;
    }
}

} // namespace tulpaObs

#endif // TULPAOCC_OCC_LIKELIHOOD_H
