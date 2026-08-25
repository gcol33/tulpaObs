// occu_cover_nuts.cpp - NUTS target for the non-spatial joint occupancy +
// cover-hurdle family (occu_cover()).
//
// The Laplace fit (R/occu_cover.R -> .tobs_fit_occu_cover) sums the latent
// occupancy state z out in closed form (two states per cell) and returns a
// Gaussian observed-Fisher posterior over the packed coefficient vector
//   theta = c(beta_psi, beta_p, beta_pos, log_dispersion).
// NUTS instead samples the exact marginal posterior of that vector, which gives
// calibrated (non-Gaussian) intervals and the per-draw pointwise likelihood
// WAIC / LOO need.
//
// The bare non-spatial target carries no latent structure at all, so its
// parameter vector is just the flat three-arm coefficient block plus one
// log-dispersion scalar. Two optional blocks extend it, each trailing that
// vector and each a no-op when absent: a coupled areal field on the latent state
// (nuts_field_hyper.h) and one observation-arm random intercept per grouping
// factor on the detection / positive-cover arms (nuts_re_block.h). The joint
// log-posterior of the bare vector is
//
//   log p(theta | y) = sum_i log m_i(theta)            # per-cell two-state marginal
//                      - 0.5 ||beta||^2 / sigma_beta^2  # weak Gaussian coef priors
//                      - 0.5  log_disp^2 / sigma_ld^2   # weak log-dispersion prior
//
// where m_i is the exact occu_cover cell marginal:
//
//   any detection : log psi + sum_v [ y log p + (1-y) log(1-p) ]
//                            + sum_{v: y=1} log f_pos(y_pos; eta_pos, disp)
//   no detection  : log( psi * prod_v (1 - p_v) + (1 - psi) )
//
// The per-cell occupancy / detection mixture reuses the canonical no-detection
// block (nodet_mixture_block) and the cover arm reuses the LognormalPositive /
// BetaPositive policies (log-density, eta-gradient, and the log-dispersion
// score) from occu_coupling_shared.h, so the sampler shares the Laplace path's
// likelihood math with no new derivation. The coefficient gradient is the
// design-sandwiched eta-gradient. The R oracle (.tobs_occu_cover_nuts_logpost)
// mirrors this target and is cross-checked byte-for-byte before it drives
// tulpa's NUTS engine.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "tobs_shape.h"
#include "occu_coupling_shared.h"
#include "nuts_field_hyper.h"
#include "nuts_re_block.h"
#include "nuts_engine.h"

namespace tulpaObs {

// Packed, NUTS-ready view of one bound non-spatial occu_cover model. The design
// matrices are held as light Rcpp handles into the (call-protected) spec list;
// every fit reaches its data through this struct. Visit-level design rows are in
// site-major order: row (i, v) = i * max_visits + v (0-based), matching the R
// builder .occu_cover_eta_from_par (matrix(..., byrow = TRUE)).
struct OccuCoverNutsData {
    int n_sites = 0;
    int max_visits = 0;
    int pos_code = 0;       // 0 lognormal, 3 beta, 4 gaussian (#112)

    Rcpp::IntegerMatrix y;        // [n_sites x max_visits] detection 0/1 (NA -> 0)
    Rcpp::NumericMatrix y_pos;    // [n_sites x max_visits] cover (0 where y != 1)
    Rcpp::IntegerMatrix valid;    // [n_sites x max_visits] visit observed 0/1
    Rcpp::NumericMatrix X_occ;        // [n_sites x p_occ]
    Rcpp::NumericMatrix X_det_site;   // [n_sites x p_det_site]
    Rcpp::NumericMatrix X_det_visit;  // [(n_sites*max_visits) x p_det_visit] (0 cols ok)
    Rcpp::NumericMatrix X_pos_site;   // [n_sites x p_pos_site]
    Rcpp::NumericMatrix X_pos_visit;  // [(n_sites*max_visits) x p_pos_visit] (0 cols ok)

    int p_occ = 0, p_det_site = 0, p_det_visit = 0, p_pos_site = 0, p_pos_visit = 0;
    int p_p = 0, p_pos = 0, total = 0;

    // Optional coupled areal field(s) on the latent state z (the spatial NUTS
    // path). Each non-centered field enters psi additively, weighted by its own
    // per-site design column, and the cover arm scaled by its own copy
    // amplitude alpha:
    //   eta_psi_c  += sum_b w_b(c) z_b[cell(c)],
    //   eta_pos_cv += sum_b alpha_b w_b(c) z_b[cell(c)].
    // An unweighted block is the shared intercept field; a weighted one is a
    // spatially-varying coefficient. nuts_field_hyper.h owns the field
    // parameterisation: a fixed basis with the field SD sigma, the mixing /
    // correlation rho and alpha either sampled as bounded coordinates or pinned
    // at a nested-Laplace estimate (the #74 / #113 fixed-hyper loading, which
    // marshals as the pinned configuration and is byte-identical to it). An
    // empty list leaves every field branch guarded, so the non-spatial sampler
    // is unchanged.
    std::vector<HyperFieldBlock> fb;

    // Optional observation-arm random intercepts. Each grouping factor on the
    // detection or positive-cover formula is one shared ReBlock
    // (nuts_re_block.h), non-centered with its own SAMPLED log_sigma_re,
    // carrying ONE group code per observation ROW (site-major, row = i * J + v)
    // rather than per site -- the observation arms are per-visit. Crossed /
    // nested groupings are simply several blocks. `re_det` / `re_pos` index into
    // `re` so a visit adds only the blocks that load on its own arm.
    std::vector<ReBlock> re;
    std::vector<int> re_det, re_pos;
};

inline OccuCoverNutsData occu_cover_nuts_build_data(const Rcpp::List& spec) {
    OccuCoverNutsData d;
    d.n_sites    = Rcpp::as<int>(spec["n_sites"]);
    d.max_visits = Rcpp::as<int>(spec["max_visits"]);
    d.pos_code   = Rcpp::as<int>(spec["pos_code"]);
    d.y          = Rcpp::as<Rcpp::IntegerMatrix>(spec["y"]);
    d.y_pos      = Rcpp::as<Rcpp::NumericMatrix>(spec["y_pos"]);
    d.valid      = Rcpp::as<Rcpp::IntegerMatrix>(spec["valid"]);
    d.X_occ      = Rcpp::as<Rcpp::NumericMatrix>(spec["X_occ"]);
    d.X_det_site = Rcpp::as<Rcpp::NumericMatrix>(spec["X_det_site"]);
    d.X_det_visit= Rcpp::as<Rcpp::NumericMatrix>(spec["X_det_visit"]);
    d.X_pos_site = Rcpp::as<Rcpp::NumericMatrix>(spec["X_pos_site"]);
    d.X_pos_visit= Rcpp::as<Rcpp::NumericMatrix>(spec["X_pos_visit"]);

    d.p_occ      = d.X_occ.ncol();
    d.p_det_site = d.X_det_site.ncol();
    d.p_det_visit= d.X_det_visit.ncol();
    d.p_pos_site = d.X_pos_site.ncol();
    d.p_pos_visit= d.X_pos_visit.ncol();
    d.p_p   = d.p_det_site + d.p_det_visit;
    d.p_pos = d.p_pos_site + d.p_pos_visit;
    int base = d.p_occ + d.p_p + d.p_pos + 1;   // +1 log_dispersion
    d.fb = hyper_field_build_list(spec, base, d.n_sites);
    base += hyper_field_list_size(d.fb);

    // RE blocks trail the field block, so a fit that carries no random effect
    // keeps the exact coefficient / field layout the sampler had before.
    // Arm codes: 1 = detection, 2 = positive cover (0, the latent state arm, is
    // the grid-integrated engine's own RE and never reaches this target).
    d.re = re_block_build_list(spec, "re_blocks", base,
                               d.n_sites * d.max_visits, /*max_arm=*/2);
    for (std::size_t b = 0; b < d.re.size(); ++b) {
        if (!d.re[b].active()) continue;
        if (d.re[b].arm == 1)      d.re_det.push_back((int) b);
        else if (d.re[b].arm == 2) d.re_pos.push_back((int) b);
        else Rcpp::stop("occu_cover NUTS: re_arm 0 (the latent state arm) is not "
                        "sampled by this target; use method = \"nested_laplace\".");
    }
    d.total = base + re_block_list_size(d.re);
    return d;
}

// log-posterior + gradient over theta = [beta_psi | beta_p | beta_pos | log_disp]
// (beta_p / beta_pos each = site-level block then visit-level block). Writes the
// full gradient into `grad` (length d.total). NUTS maximises, so this returns the
// (un-negated) log-posterior.
inline double occu_cover_nuts_eval(const OccuCoverNutsData& d, const double* theta,
                                   double sigma_beta, double sigma_logdisp,
                                   double* grad) {
    const int N = d.n_sites, J = d.max_visits;
    const int p_occ = d.p_occ;
    const int p_det_site = d.p_det_site, p_det_visit = d.p_det_visit;
    const int p_pos_site = d.p_pos_site, p_pos_visit = d.p_pos_visit;
    const int p_p = d.p_p, p_pos = d.p_pos, total = d.total;
    // The coefficient block runs to n_coef; log_dispersion sits at index n_coef
    // and any trailing raw-field block (spatial NUTS) follows it. log_disp /
    // g_ld / n_beta index into this leading block, never the field tail.
    const int n_coef     = p_occ + p_p + p_pos;

    const double* bo         = theta;
    const double* bp_site    = theta + p_occ;
    const double* bp_visit   = theta + p_occ + p_det_site;
    const double* bpos_site  = theta + p_occ + p_p;
    const double* bpos_visit = theta + p_occ + p_p + p_pos_site;
    const double  log_disp   = theta[n_coef];
    const double  disp       = std::exp(log_disp);

    // Gradient block base offsets into `grad`.
    const int g_bo        = 0;
    const int g_bp_site   = p_occ;
    const int g_bp_visit  = p_occ + p_det_site;
    const int g_bpos_site = p_occ + p_p;
    const int g_bpos_visit= p_occ + p_p + p_pos_site;
    const int g_ld        = n_coef;

    for (int k = 0; k < total; ++k) grad[k] = 0.0;
    double lp = 0.0;
    double g_logdisp = 0.0;

    // Coupled field(s): z per cell (nuts_field_hyper.h), each entering psi with
    // its own per-site weight and the cover arm scaled by its own copy amplitude
    // alpha. Per block the per-cell field score is accumulated in g_f[b] and the
    // alpha data-score in g_alpha_data[b]; both are handed back to that block
    // below.
    const std::size_t n_fb = d.fb.size();
    const bool has_field = n_fb > 0;
    std::vector<HyperFieldState> fs(n_fb);
    std::vector<double> g_alpha_data(n_fb, 0.0);
    std::vector<std::vector<double> > g_f(n_fb);
    for (std::size_t b = 0; b < n_fb; ++b) {
        hyper_field_forward(d.fb[b], theta, fs[b]);
        g_f[b].assign(d.fb[b].n_units, 0.0);
    }

    // Observation-arm random intercepts: each block's sampled sigma_re, and the
    // running data score of its log_sigma_re coordinate.
    const std::size_t n_re = d.re.size();
    std::vector<double> re_sig(n_re, 0.0), re_glsig(n_re, 0.0);
    for (std::size_t b = 0; b < n_re; ++b)
        re_sig[b] = re_block_sigma(d.re[b], theta);

    std::vector<double> eta_p(J), g_eta_p(J), g_eta_pos(J), eta_p_compact(J), g_p_compact(J);

    for (int i = 0; i < N; ++i) {
        // Field offsets this site carries on the two coupled arms.
        double f_psi = 0.0, f_pos = 0.0;
        if (has_field) hyper_field_site_offsets(d.fb, fs, i, f_psi, f_pos);
        // Occupancy linear predictor (+ the shared field(s) on psi).
        double eta_psi = f_psi;
        for (int k = 0; k < p_occ; ++k) eta_psi += d.X_occ(i, k) * bo[k];
        const double psi = sigmoid_(eta_psi);

        // Per-visit detection linear predictor (site block broadcast + visit block).
        double eta_p_site = 0.0;
        for (int k = 0; k < p_det_site; ++k) eta_p_site += d.X_det_site(i, k) * bp_site[k];
        bool any_det = false;
        int n_valid = 0;
        for (int v = 0; v < J; ++v) {
            g_eta_p[v] = 0.0; g_eta_pos[v] = 0.0;
            if (d.valid(i, v) == 0) { eta_p[v] = 0.0; continue; }
            double e = eta_p_site;
            if (p_det_visit > 0) {
                const int row = i * J + v;
                for (int k = 0; k < p_det_visit; ++k) e += d.X_det_visit(row, k) * bp_visit[k];
            }
            for (std::size_t t = 0; t < d.re_det.size(); ++t) {
                const int b = d.re_det[t];
                e += re_block_offset(d.re[b], re_sig[b], theta, i * J + v);
            }
            eta_p[v] = e;
            ++n_valid;
            if (d.y(i, v) == 1) any_det = true;
        }

        double g_eta_psi = 0.0;
        if (any_det) {
            lp += log_safe(psi);
            g_eta_psi = 1.0 - psi;
            for (int v = 0; v < J; ++v) {
                if (d.valid(i, v) == 0) continue;
                const double pv = sigmoid_(eta_p[v]);
                if (d.y(i, v) == 1) { lp += log_safe(pv);       g_eta_p[v] = 1.0 - pv; }
                else                { lp += log_safe(1.0 - pv); g_eta_p[v] = -pv; }
            }
            // Cover arm at detected visits (+ the shared field scaled by alpha).
            for (int v = 0; v < J; ++v) {
                if (d.valid(i, v) == 0 || d.y(i, v) != 1) continue;
                const double yp = d.y_pos(i, v);
                // Missing-at-random cover: a detected visit with no cover value
                // (NA -> non-finite) drops out of the cover factor; its cover-arm
                // score stays 0.
                if (!std::isfinite(yp)) continue;
                double eta_pos = f_pos;
                for (int k = 0; k < p_pos_site; ++k) eta_pos += d.X_pos_site(i, k) * bpos_site[k];
                if (p_pos_visit > 0) {
                    const int row = i * J + v;
                    for (int k = 0; k < p_pos_visit; ++k) eta_pos += d.X_pos_visit(row, k) * bpos_visit[k];
                }
                for (std::size_t t = 0; t < d.re_pos.size(); ++t) {
                    const int b = d.re_pos[t];
                    eta_pos += re_block_offset(d.re[b], re_sig[b], theta, i * J + v);
                }
                lp += pos_log_density(d.pos_code, yp, eta_pos, disp);
                g_eta_pos[v] = pos_grad_eta(d.pos_code, yp, eta_pos, disp);
                g_logdisp   += pos_grad_logdisp(d.pos_code, yp, eta_pos, disp);
            }
        } else {
            // No-detection mixture L = psi * prod_v(1-p_v) + (1-psi). Reuse the
            // canonical block over the cell's valid visits (compact buffers), then
            // scatter the per-visit scores back to the J-wide slots.
            int nv = 0;
            for (int v = 0; v < J; ++v) {
                if (d.valid(i, v) == 0) continue;
                eta_p_compact[nv] = eta_p[v];
                ++nv;
            }
            double g_w = 0.0, nh_w = 0.0;
            const double cell_ll = nodet_mixture_block(
                psi, eta_p_compact.data(), nv, /*want_hess=*/false, /*expected=*/false,
                g_w, nh_w, g_p_compact.data(), nullptr, nullptr, nullptr);
            lp += cell_ll;
            g_eta_psi = g_w;
            int j = 0;
            for (int v = 0; v < J; ++v) {
                if (d.valid(i, v) == 0) continue;
                g_eta_p[v] = g_p_compact[j];
                ++j;
            }
        }

        // Design-sandwich the eta-gradients onto the coefficient blocks.
        for (int k = 0; k < p_occ; ++k) grad[g_bo + k] += g_eta_psi * d.X_occ(i, k);

        double g_eta_p_sum = 0.0, g_eta_pos_sum = 0.0;
        for (int v = 0; v < J; ++v) { g_eta_p_sum += g_eta_p[v]; g_eta_pos_sum += g_eta_pos[v]; }
        for (int k = 0; k < p_det_site; ++k) grad[g_bp_site + k] += g_eta_p_sum * d.X_det_site(i, k);
        for (int k = 0; k < p_pos_site; ++k) grad[g_bpos_site + k] += g_eta_pos_sum * d.X_pos_site(i, k);

        // Field score for this cell: the psi-arm eta-score plus alpha times the
        // cover-arm eta-score sum (the cover arm sees alpha * f), both weighted
        // by the block's own per-site design column. Each copy amplitude's own
        // data-score is the mirror term, its loading times that sum.
        if (has_field)
            hyper_field_site_score(d.fb, fs, i, g_eta_psi, g_eta_pos_sum,
                                   g_f, g_alpha_data);
        if (p_det_visit > 0) {
            for (int v = 0; v < J; ++v) {
                if (g_eta_p[v] == 0.0) continue;
                const int row = i * J + v;
                for (int k = 0; k < p_det_visit; ++k)
                    grad[g_bp_visit + k] += g_eta_p[v] * d.X_det_visit(row, k);
            }
        }
        if (p_pos_visit > 0) {
            for (int v = 0; v < J; ++v) {
                if (g_eta_pos[v] == 0.0) continue;
                const int row = i * J + v;
                for (int k = 0; k < p_pos_visit; ++k)
                    grad[g_bpos_visit + k] += g_eta_pos[v] * d.X_pos_visit(row, k);
            }
        }

        // Chain the per-visit arm scores into each RE block's whitened z and its
        // log_sigma_re; an invalid visit carries a zero eta-score anyway.
        for (std::size_t t = 0; t < d.re_det.size(); ++t) {
            const int b = d.re_det[t];
            for (int v = 0; v < J; ++v) {
                if (d.valid(i, v) == 0) continue;
                re_block_accumulate(d.re[b], re_sig[b], g_eta_p[v], i * J + v,
                                    theta, grad, re_glsig[b]);
            }
        }
        for (std::size_t t = 0; t < d.re_pos.size(); ++t) {
            const int b = d.re_pos[t];
            for (int v = 0; v < J; ++v) {
                if (g_eta_pos[v] == 0.0) continue;
                re_block_accumulate(d.re[b], re_sig[b], g_eta_pos[v], i * J + v,
                                    theta, grad, re_glsig[b]);
            }
        }
    }

    grad[g_ld] = g_logdisp;

    // Weak Gaussian priors: N(0, sigma_beta^2) on every coefficient (the
    // data-dominated ridge, matching the Laplace path's default) and a broad
    // N(0, sigma_logdisp^2) on log_disp to keep the dispersion proper.
    const double ib2 = 1.0 / (sigma_beta * sigma_beta);
    const int n_beta = n_coef;   // log_disp + field tail keep their own priors
    for (int k = 0; k < n_beta; ++k) {
        lp        -= 0.5 * ib2 * theta[k] * theta[k];
        grad[k]   -= ib2 * theta[k];
    }
    const double ild2 = 1.0 / (sigma_logdisp * sigma_logdisp);
    lp           -= 0.5 * ild2 * log_disp * log_disp;
    grad[g_ld]   -= ild2 * log_disp;

    // Whitened field prior, the chain rule onto raw, and any sampled hyper's
    // bounded-transform log-density + gradient, per block.
    for (std::size_t b = 0; b < n_fb; ++b)
        lp += hyper_field_backward(d.fb[b], theta, fs[b], g_f[b],
                                   g_alpha_data[b], grad);

    // Whitened RE prior z ~ N(0, I) plus the sampled log_sigma_re's own
    // N(0, sigma_lsd^2) prior, and that coordinate's accumulated data score.
    for (std::size_t b = 0; b < n_re; ++b)
        lp += re_block_backward(d.re[b], theta, re_glsig[b], grad);

    return lp;
}

// Model wrapper handed to the shared NUTS engine through ModelData; the
// FullGradFn reaches it via ModelData.model_response_data.
struct OccuCoverNutsModel {
    OccuCoverNutsData d;
    double sigma_beta = 5.0;
    double sigma_logdisp = 5.0;
};

inline void occu_cover_nuts_full_grad(const std::vector<double>& params,
                                      const tulpa::ModelData& data,
                                      const tulpa::ParamLayout& /*layout*/,
                                      std::vector<double>& grad, double* log_post_out) {
    const OccuCoverNutsModel* m =
        static_cast<const OccuCoverNutsModel*>(data.model_response_data);
    grad.assign((std::size_t) m->d.total, 0.0);
    const double lp = occu_cover_nuts_eval(m->d, params.data(),
                                           m->sigma_beta, m->sigma_logdisp, grad.data());
    if (log_post_out) *log_post_out = lp;
}

}  // namespace tulpaObs

// Cross-check entry: the full-vector joint log-posterior + gradient, validated
// byte-for-byte against the R oracle .tobs_occu_cover_nuts_logpost.
// [[Rcpp::export]]
Rcpp::List cpp_occu_cover_nuts_joint_logpost(Rcpp::List spec, Rcpp::NumericVector theta,
                                             double sigma_beta, double sigma_logdisp) {
    tulpaObs::OccuCoverNutsData d = tulpaObs::occu_cover_nuts_build_data(spec);
    if ((int) theta.size() != d.total)
        Rcpp::stop("theta length %d != expected %d", (int) theta.size(), d.total);
    Rcpp::NumericVector grad(d.total);
    const double lp = tulpaObs::occu_cover_nuts_eval(
        d, theta.begin(), sigma_beta, sigma_logdisp, grad.begin());
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// Sample the exact non-spatial occu_cover coefficient posterior via tulpa's NUTS
// engine and the in-tree FullGradFn, warm-started at the Laplace mode with a
// diagonal Laplace metric. Returns draws + sampler diagnostics.
// [[Rcpp::export]]
Rcpp::List cpp_occu_cover_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                               double sigma_beta, double sigma_logdisp,
                               Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                               int n_iter, int n_warmup, int max_treedepth,
                               double adapt_delta, int seed, bool verbose) {
    tulpaObs::OccuCoverNutsModel m;
    m.d = tulpaObs::occu_cover_nuts_build_data(spec);
    m.sigma_beta = sigma_beta;
    m.sigma_logdisp = sigma_logdisp;
    return tulpaObs::run_tulpa_nuts(
        &tulpaObs::occu_cover_nuts_full_grad, &m, m.d.total, theta0, sigma_beta,
        tulpaObs::shape::optional_numeric(inv_metric.get(), "inv_metric"), n_iter, n_warmup, max_treedepth, adapt_delta, seed, verbose);
}
