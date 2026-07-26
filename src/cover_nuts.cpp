// cover_nuts.cpp - NUTS target for the non-spatial standalone cover hurdle
// family (cover()).
//
// The cover hurdle is occu_cover() minus the occupancy / detection latent
// mixture: there is no z to marginalise. The two arms are conditionally
// independent given the data, so the joint log-posterior over the packed
// coefficient vector
//   theta = c(beta_presence, beta_positive, log_dispersion)
// is a plain sum of per-row terms:
//
//   sum_i [ presence Bernoulli logpdf
//           + (present_i ? positive logpdf : 0) ]
//   - 0.5 ||beta||^2 / sigma_beta^2      # weak Gaussian coef priors
//   - 0.5  log_disp^2 / sigma_ld^2       # weak log-dispersion prior
//
// The presence arm is a Bernoulli on 1{cover > 0} with predictor
// X_presence beta_presence (logit). The positive arm is the beta / lognormal
// density on cover | cover > 0 with predictor X_pos beta_pos and one
// dispersion scalar log_disp (log phi for beta, log sigma for lognormal),
// evaluated only at present rows. The positive-arm density, eta-gradient, and
// log-dispersion score reuse the LognormalPositive / BetaPositive policies from
// occu_coupling_shared.h, so the sampler shares the cover Laplace path's
// likelihood math with no new derivation. The coefficient gradient is the
// design-sandwiched eta-gradient. The R oracle (.tobs_cover_nuts_logpost)
// mirrors this target and is cross-checked byte-for-byte before it drives
// tulpa's NUTS engine.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "occu_coupling_shared.h"
#include "nuts_engine.h"

namespace tulpaObs {

// Packed, NUTS-ready view of one bound non-spatial cover hurdle. The presence
// design has one row per observation; the positive design has one row per
// present observation (cover > 0), in the same order the R encoder slices it
// (enc$idx_pos). `y_pos` is the per-present-row response: the raw cover in
// (0, 1) for beta, the natural-scale cover for lognormal (the policy takes the
// log internally).
struct CoverNutsData {
    int n_obs = 0;          // presence-arm rows
    int n_pos = 0;          // positive-arm rows (present observations)
    int pos_code = 0;       // 0 lognormal, 3 beta, 4 gaussian (#112)

    Rcpp::IntegerVector present;   // [n_obs] 1{cover > 0}
    Rcpp::NumericVector y_pos;     // [n_pos] cover at present rows
    Rcpp::NumericMatrix X_pres;    // [n_obs x p_pres]
    Rcpp::NumericMatrix X_pos;     // [n_pos x p_pos]

    int p_pres = 0, p_pos = 0, total = 0;
};

inline CoverNutsData cover_nuts_build_data(const Rcpp::List& spec) {
    CoverNutsData d;
    d.pos_code = Rcpp::as<int>(spec["pos_code"]);
    d.present  = Rcpp::as<Rcpp::IntegerVector>(spec["present"]);
    d.y_pos    = Rcpp::as<Rcpp::NumericVector>(spec["y_pos"]);
    d.X_pres   = Rcpp::as<Rcpp::NumericMatrix>(spec["X_pres"]);
    d.X_pos    = Rcpp::as<Rcpp::NumericMatrix>(spec["X_pos"]);

    d.n_obs  = d.X_pres.nrow();
    d.n_pos  = d.X_pos.nrow();
    d.p_pres = d.X_pres.ncol();
    d.p_pos  = d.X_pos.ncol();
    d.total  = d.p_pres + d.p_pos + 1;
    return d;
}

// log-posterior + gradient over theta = [beta_presence | beta_pos | log_disp].
// Writes the full gradient into `grad` (length d.total). NUTS maximises, so
// this returns the (un-negated) log-posterior.
inline double cover_nuts_eval(const CoverNutsData& d, const double* theta,
                              double sigma_beta, double sigma_logdisp,
                              double* grad) {
    const int n_obs = d.n_obs, n_pos = d.n_pos;
    const int p_pres = d.p_pres, p_pos = d.p_pos, total = d.total;

    const double* b_pres   = theta;
    const double* b_pos    = theta + p_pres;
    const double  log_disp = theta[total - 1];
    const double  disp     = std::exp(log_disp);

    const int g_bpres = 0;
    const int g_bpos  = p_pres;
    const int g_ld    = total - 1;

    for (int k = 0; k < total; ++k) grad[k] = 0.0;
    double lp = 0.0;
    double g_logdisp = 0.0;

    // Presence arm: plain Bernoulli on 1{cover > 0}.
    for (int i = 0; i < n_obs; ++i) {
        double eta = 0.0;
        for (int k = 0; k < p_pres; ++k) eta += d.X_pres(i, k) * b_pres[k];
        const double pr = sigmoid_(eta);
        double g_eta;
        if (d.present[i] == 1) { lp += log_safe(pr);       g_eta = 1.0 - pr; }
        else                   { lp += log_safe(1.0 - pr); g_eta = -pr; }
        for (int k = 0; k < p_pres; ++k) grad[g_bpres + k] += g_eta * d.X_pres(i, k);
    }

    // Positive arm: beta / lognormal density at the present rows.
    for (int j = 0; j < n_pos; ++j) {
        double eta = 0.0;
        for (int k = 0; k < p_pos; ++k) eta += d.X_pos(j, k) * b_pos[k];
        const double yp = d.y_pos[j];
        lp        += pos_log_density(d.pos_code, yp, eta, disp);
        double g_eta = pos_grad_eta(d.pos_code, yp, eta, disp);
        g_logdisp += pos_grad_logdisp(d.pos_code, yp, eta, disp);
        for (int k = 0; k < p_pos; ++k) grad[g_bpos + k] += g_eta * d.X_pos(j, k);
    }
    grad[g_ld] = g_logdisp;

    // Weak Gaussian priors: N(0, sigma_beta^2) on every coefficient and a broad
    // N(0, sigma_logdisp^2) on log_disp to keep the dispersion proper.
    const double ib2 = 1.0 / (sigma_beta * sigma_beta);
    const int n_beta = total - 1;
    for (int k = 0; k < n_beta; ++k) {
        lp      -= 0.5 * ib2 * theta[k] * theta[k];
        grad[k] -= ib2 * theta[k];
    }
    const double ild2 = 1.0 / (sigma_logdisp * sigma_logdisp);
    lp         -= 0.5 * ild2 * log_disp * log_disp;
    grad[g_ld] -= ild2 * log_disp;

    return lp;
}

// Model wrapper handed to the shared NUTS engine through ModelData; the
// FullGradFn reaches it via ModelData.model_response_data.
struct CoverNutsModel {
    CoverNutsData d;
    double sigma_beta = 5.0;
    double sigma_logdisp = 5.0;
};

inline void cover_nuts_full_grad(const std::vector<double>& params,
                                 const tulpa::ModelData& data,
                                 const tulpa::ParamLayout& /*layout*/,
                                 std::vector<double>& grad, double* log_post_out) {
    const CoverNutsModel* m =
        static_cast<const CoverNutsModel*>(data.model_response_data);
    grad.assign((std::size_t) m->d.total, 0.0);
    const double lp = cover_nuts_eval(m->d, params.data(),
                                      m->sigma_beta, m->sigma_logdisp, grad.data());
    if (log_post_out) *log_post_out = lp;
}

}  // namespace tulpaObs

// Cross-check entry: the full-vector joint log-posterior + gradient, validated
// byte-for-byte against the R oracle .tobs_cover_nuts_logpost.
// [[Rcpp::export]]
Rcpp::List cpp_cover_nuts_logpost(Rcpp::List spec, Rcpp::NumericVector theta,
                                  double sigma_beta, double sigma_logdisp) {
    tulpaObs::CoverNutsData d = tulpaObs::cover_nuts_build_data(spec);
    if ((int) theta.size() != d.total)
        Rcpp::stop("theta length %d != expected %d", (int) theta.size(), d.total);
    Rcpp::NumericVector grad(d.total);
    const double lp = tulpaObs::cover_nuts_eval(
        d, theta.begin(), sigma_beta, sigma_logdisp, grad.begin());
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// Sample the exact non-spatial cover-hurdle coefficient posterior via tulpa's
// NUTS engine and the in-tree FullGradFn, warm-started at the Laplace mode with
// a diagonal Laplace metric. Returns draws + sampler diagnostics.
// [[Rcpp::export]]
Rcpp::List cpp_cover_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                          double sigma_beta, double sigma_logdisp,
                          Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                          int n_iter, int n_warmup, int max_treedepth,
                          double adapt_delta, int seed, bool verbose) {
    tulpaObs::CoverNutsModel m;
    m.d = tulpaObs::cover_nuts_build_data(spec);
    m.sigma_beta = sigma_beta;
    m.sigma_logdisp = sigma_logdisp;
    return tulpaObs::run_tulpa_nuts(
        &tulpaObs::cover_nuts_full_grad, &m, m.d.total, theta0, sigma_beta,
        inv_metric, n_iter, n_warmup, max_treedepth, adapt_delta, seed, verbose);
}
