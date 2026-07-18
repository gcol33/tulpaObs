// fp_occu_kernel.h
// Per-site marginal log-likelihood + gradients for the Miller et al. (2011)
// multistate false-positive occupancy model. Detections are classified into
// three states y_ij in {0, 1, 2}:
//   0 = no detection
//   1 = ambiguous detection (a true detection OR a false positive)
//   2 = certain / confirmed detection (only possible when the site is occupied)
//
//   z_i ~ Bernoulli(psi_i),  logit psi_i = X_psi beta_psi
//   per-visit cell probabilities given the latent state:
//     occupied  (z=1): P(0)=1-p11, P(1)=p11(1-b), P(2)=p11 b
//     unoccupied(z=0): P(0)=1-p10, P(1)=p10,      P(2)=0
//   logit p11 = X_p11 beta_p11   (true detection probability)
//   logit p10 = X_p10 beta_p10   (false-positive probability)
//   logit b   = X_b   beta_b     (prob a true detection is certain/confirmed)
//
// The latent occupancy z marginalises in closed form (two states):
//   L_i = psi_i prod_j g1_ij + (1 - psi_i) prod_j g0_ij
// A certain detection (y=2) makes the unoccupied term zero, which is what
// identifies the model (Miller et al. 2011; the unmarked::occuFP confirmation
// design). The arms are site-level here (one logit predictor per arm per site;
// p11 constant across a site's visits), so every gradient is a scalar per site
// and sandwiches with the arm design in the driver. The detection score's latent
// coefficient gives the clean occupancy-style gradient grad_eta_psi = w1 - psi
// with w1 = P(z=1 | y_i) the posterior occupancy.
//
// References:
//   Miller, D. A. W., Nichols, J. D., McClintock, B. T., Grant, E. H. C., Bailey,
//     L. L., Weir, L. A. (2011) Ecology 92: 1422-1428.
//   Royle, J. A., Link, W. A. (2006) Ecology 87: 835-841.

#ifndef TULPAOBS_FP_OCCU_KERNEL_H
#define TULPAOBS_FP_OCCU_KERNEL_H

#include <cmath>
#include <limits>

namespace tulpaObs {

struct FpOccuSiteResult {
    double log_lik;
    double grad_eta_psi, grad_eta_p11, grad_eta_p10, grad_eta_b;
    double w1;                 // posterior occupancy P(z = 1 | y_i)
};

// Stable inverse logit.
inline double fp_inv_logit(double e) {
    if (e > 0.0) return 1.0 / (1.0 + std::exp(-e));
    double ee = std::exp(e);
    return ee / (1.0 + ee);
}

// log p and log(1-p) under the logit link (numerically stable).
inline void fp_logit_log_probs(double eta, double& lp, double& l1mp) {
    if (eta > 0.0) {
        double sp = std::log1p(std::exp(-eta));
        lp = -sp; l1mp = -eta - sp;
    } else {
        double sp = std::log1p(std::exp(eta));
        lp = eta - sp; l1mp = -sp;
    }
}

// Per-site false-positive occupancy marginal. `y` are the J_i valid detection
// states (NA visits dropped upstream); the four eta are the site-level logit
// predictors. Returns the marginal log-lik, the per-arm gradients (eta scale),
// and the posterior occupancy.
inline FpOccuSiteResult compute_fp_occu_site(
    const int* y, int J,
    double eta_psi, double eta_p11, double eta_p10, double eta_b
) {
    FpOccuSiteResult res;
    const double psi = fp_inv_logit(eta_psi);
    const double p11 = fp_inv_logit(eta_p11);
    const double p10 = fp_inv_logit(eta_p10);
    const double b   = fp_inv_logit(eta_b);

    double lp11, l1mp11, lp10, l1mp10, lb, l1mb;
    fp_logit_log_probs(eta_p11, lp11, l1mp11);
    fp_logit_log_probs(eta_p10, lp10, l1mp10);
    fp_logit_log_probs(eta_b,   lb,   l1mb);

    double logA = 0.0, logB = 0.0;
    bool   B_zero = false;
    // Gradient accumulators of the within-state log-likelihoods.
    double dA_p11 = 0.0, dA_b = 0.0, dB_p10 = 0.0;
    for (int j = 0; j < J; ++j) {
        const int yj = y[j];
        if (yj == 0) {
            logA += l1mp11;                 // log(1 - p11)
            dA_p11 += -p11;                 // d log(1-p11) / d eta_p11
            logB += l1mp10;                 // log(1 - p10)
            dB_p10 += -p10;
        } else if (yj == 1) {
            logA += lp11 + l1mb;            // log p11 + log(1 - b)
            dA_p11 += (1.0 - p11);
            dA_b   += -b;
            logB += lp10;                   // log p10
            dB_p10 += (1.0 - p10);
        } else {                            // yj == 2 (certain detection)
            logA += lp11 + lb;              // log p11 + log b
            dA_p11 += (1.0 - p11);
            dA_b   += (1.0 - b);
            B_zero = true;                  // impossible when unoccupied
        }
    }

    double lpsi, l1mpsi;
    fp_logit_log_probs(eta_psi, lpsi, l1mpsi);  // stable log(psi), log(1-psi)
    const double t1 = lpsi + logA;          // log(psi * A)
    double log_lik, w1;
    if (B_zero) {
        log_lik = t1;
        w1 = 1.0;
        dB_p10 = 0.0;
    } else {
        const double t0 = l1mpsi + logB;    // log((1-psi) * B)
        const double m = (t1 > t0) ? t1 : t0;
        log_lik = m + std::log(std::exp(t1 - m) + std::exp(t0 - m));
        w1 = std::exp(t1 - log_lik);
    }
    const double w0 = 1.0 - w1;

    res.log_lik      = log_lik;
    res.w1           = w1;
    res.grad_eta_psi = w1 - psi;            // posterior occupancy minus prior
    res.grad_eta_p11 = w1 * dA_p11;
    res.grad_eta_b   = w1 * dA_b;
    res.grad_eta_p10 = w0 * dB_p10;
    return res;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_FP_OCCU_KERNEL_H
