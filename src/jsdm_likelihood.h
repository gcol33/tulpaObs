// jsdm_likelihood.h
// Joint Species Distribution Model likelihood (no detection process).
// Single process: occupancy — eta[0] = logit(psi_i)
// Standard Bernoulli: y_i ~ Bernoulli(psi_i)

#ifndef TULPAOCC_JSDM_LIKELIHOOD_H
#define TULPAOCC_JSDM_LIKELIHOOD_H

#include <cmath>
#include <vector>
#include <type_traits>
#include "occ_likelihood.h"  // For log_inv_logit, log1m_inv_logit
#include <tulpa/likelihood.h>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/autodiff_arena.h>
#include <tulpa/autodiff_fwd.h>

namespace tulpaOcc {

// JSDM response: simple binary presence/absence vector
struct JSDMResponseData {
    int n_obs;                 // N = n_sites * n_species
    std::vector<int> y;        // Binary: 0/1
};

template<typename T>
T jsdm_log_likelihood(
    int i,
    const T* eta,
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

    const auto* jsdm = static_cast<const JSDMResponseData*>(model_data);
    int y_i = jsdm->y[i];

    if (y_i == 1) {
        return log_inv_logit(eta[0]);
    } else {
        return log1m_inv_logit(eta[0]);
    }
}

} // namespace tulpaOcc

#endif // TULPAOCC_JSDM_LIKELIHOOD_H
