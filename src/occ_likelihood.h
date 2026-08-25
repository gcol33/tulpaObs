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
#include "occu_coupling_shared.h"  // nodet_mixture_block (the canonical mixture)
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

// The detection logit for visit j at site i: the site-level eta[1] plus the
// visit-level design row. Templated because the log-likelihood assembles it in
// autodiff Var and the residual in plain double; T = double is the trivial
// instantiation and compiles to the same arithmetic.
template<typename T>
inline T occ_visit_logit_p(const OccResponseData* occ, const T* eta,
                           const std::vector<T>& params,
                           const tulpa::ParamLayout& layout, int i, int j) {
    T logit_p_ij = eta[1];
    if (occ->p_det_visit > 0) {
        const int base = i * occ->max_visits * occ->p_det_visit + j * occ->p_det_visit;
        const int beta_offset = layout.extra_offset;
        for (int c = 0; c < occ->p_det_visit; c++) {
            logit_p_ij = logit_p_ij + T(occ->X_det_visit[base + c]) * params[beta_offset + c];
        }
    }
    return logit_p_ij;
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

    for (int j = 0; j < occ->max_visits; j++) {
        int y_ij = occ->y[i * occ->max_visits + j];
        if (y_ij < 0) continue;  // Missing visit

        T logit_p_ij = occ_visit_logit_p(occ, eta, params, layout, i, j);

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

    // Detection logits for this site's valid visits, compact (one entry per
    // valid visit) so the undetected branch can hand them straight to the
    // shared no-detection mixture.
    double  eta_p_stack[64];
    std::vector<double> eta_p_heap;
    double* eta_p = eta_p_stack;
    if (occ->max_visits > 64) {
        eta_p_heap.assign(occ->max_visits, 0.0);
        eta_p = eta_p_heap.data();
    }

    int nv = 0;
    double d_logit_p_sum = 0.0;

    for (int j = 0; j < occ->max_visits; j++) {
        int y_ij = occ->y[i * occ->max_visits + j];
        if (y_ij < 0) continue;

        double logit_p_ij = occ_visit_logit_p(occ, eta, params, layout, i, j);
        eta_p[nv++] = logit_p_ij;

        // d(ll_det)/d(logit_p_ij) = y_ij - p_ij (for the part that goes through eta[1])
        // Only site-level part (visit-level betas have separate gradient)
        d_logit_p_sum += (y_ij - 1.0 / (1.0 + std::exp(-logit_p_ij)));
    }

    if (occ->any_detected[i]) {
        // d log(psi) / d logit_psi = 1 - psi
        resid_out[0] = 1.0 - psi;
        resid_out[1] = d_logit_p_sum;
        return;
    }

    // No detection here: the likelihood is the MacKenzie mixture
    // psi * prod(1 - p) + (1 - psi). That object belongs to
    // occu_coupling_shared.h, which accumulates prod(1 - p) in log space; this
    // branch used to re-derive it and recover each visit's contribution by
    // dividing (1 - p_j) back out of the product, hard-zeroing the quotient
    // below 1e-300. resid_out[1] is the site-level detection score, so the
    // block's per-visit scores are summed.
    double  g_psi = 0.0, nh_psi = 0.0;
    double  g_p_stack[64];
    std::vector<double> g_p_heap;
    double* g_p = g_p_stack;
    if (nv > 64) { g_p_heap.assign(nv, 0.0); g_p = g_p_heap.data(); }

    nodet_mixture_block(psi, eta_p, nv, /*want_hess=*/false, /*expected=*/false,
                        g_psi, nh_psi, g_p, nullptr, nullptr, nullptr);

    resid_out[0] = g_psi;
    double g_p_sum = 0.0;
    for (int v = 0; v < nv; v++) g_p_sum += g_p[v];
    resid_out[1] = g_p_sum;
}

} // namespace tulpaObs

#endif // TULPAOCC_OCC_LIKELIHOOD_H
