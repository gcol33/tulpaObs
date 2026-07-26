// distance_laplace.cpp
// Non-spatial fixed-effects Laplace (maximum-likelihood) fit for binned
// distance-sampling abundance models with a Poisson OR negative-binomial
// abundance mixing distribution and a half-normal or hazard-rate detection key.
//
// Parameter vector: beta = (beta_lambda [p_lam], beta_sigma [p_sig],
// eta_b [1, hazard-rate only]), plus the outer NB dispersion theta = log r.
// The inner Newton maximises beta against the marginal observed information
// (distance_kernel.h: E[I_c|y] minus the rank-1 Var[N|y] coupling), falling back
// to the positive-definite Fisher-scoring matrix when the observed block is not
// PSD; the outer loop profiles theta exactly as the (lambda, logit-p) families do
// (marginal_count_laplace.h). The abundance arm and the NB dispersion row/col are
// the shared count-marginal math (nmix_kernel.h); only the detection arm (a
// site-level log-scale sigma and an optional scalar hazard shape, with bin
// probabilities integrated by quadrature) is distance-specific.

#include "distance_kernel.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

namespace {

using tulpaObs::DistSiteResult;
using tulpaObs::DistQuad;
using tulpaObs::compute_distance_site;

// Slack on the backtracking accept test. A Newton step at a converged mode moves
// the objective by less than the rounding error of the sweep that measures it, so
// a strict `ll_try >= ll` would reject the final step and exhaust the halving
// loop; accepting a decrease this small lets the step through instead.
constexpr double kLineSearchSlack = 1e-10;

// Floor on the negative-binomial size r during the profile. As r grows the NB
// tends to the Poisson, so the log-r axis is only weakly identified from below
// and an unbounded profile can walk toward zero; log(r) is optimised inside
// [log(kMinNbSize), theta_max].
constexpr double kMinNbSize = 1e-4;

// One per-site kernel pass at the current (beta, r); fills `out` and returns the
// total marginal log-likelihood.
inline double dist_sweep(const Rcpp::IntegerMatrix& y,
                         const Map<MatrixXd>& Xl, const Map<MatrixXd>& Xs,
                         const VectorXd& beta_lam, const VectorXd& beta_sig,
                         double eta_b, int key, const DistQuad& quad,
                         int K_max, double r,
                         std::vector<DistSiteResult>& out) {
    const int n_sites = y.nrow(), n_bins = y.ncol();
    VectorXd eta_lam = Xl * beta_lam;
    VectorXd eta_sig = Xs * beta_sig;
    double ll = 0.0;
    std::vector<int> y_site(n_bins);
    for (int s = 0; s < n_sites; ++s) {
        for (int b = 0; b < n_bins; ++b) y_site[b] = y(s, b);
        out[s] = compute_distance_site(y_site.data(), n_bins, eta_lam(s),
                                       eta_sig(s), eta_b, key, quad, K_max, r);
        ll += out[s].log_lik;
    }
    return ll;
}

// Cheap log-lik-only sweep for line-search trial points.
inline double dist_loglik(const Rcpp::IntegerMatrix& y,
                          const Map<MatrixXd>& Xl, const Map<MatrixXd>& Xs,
                          const VectorXd& beta_lam, const VectorXd& beta_sig,
                          double eta_b, int key, const DistQuad& quad,
                          int K_max, double r) {
    const int n_sites = y.nrow(), n_bins = y.ncol();
    VectorXd eta_lam = Xl * beta_lam;
    VectorXd eta_sig = Xs * beta_sig;
    double ll = 0.0;
    std::vector<int> y_site(n_bins);
    for (int s = 0; s < n_sites; ++s) {
        for (int b = 0; b < n_bins; ++b) y_site[b] = y(s, b);
        DistSiteResult res = compute_distance_site(y_site.data(), n_bins, eta_lam(s),
                                                   eta_sig(s), eta_b, key, quad, K_max, r);
        if (!R_finite(res.log_lik)) return res.log_lik;
        ll += res.log_lik;
    }
    return ll;
}

// Beta gradient (length pb) at the current sweep state.
inline VectorXd dist_grad_beta(const Map<MatrixXd>& Xl, const Map<MatrixXd>& Xs,
                               int p_lam, int p_sig, bool hazard,
                               const std::vector<DistSiteResult>& st) {
    const int n_sites = (int)st.size();
    const int pb = p_lam + p_sig + (hazard ? 1 : 0);
    VectorXd g = VectorXd::Zero(pb);
    for (int s = 0; s < n_sites; ++s) {
        g.segment(0, p_lam)     += st[s].grad_eta_lambda * Xl.row(s).transpose();
        g.segment(p_lam, p_sig) += st[s].grad_eta_d[0]   * Xs.row(s).transpose();
        if (hazard) g(p_lam + p_sig) += st[s].grad_eta_d[1];
    }
    return g;
}

// Beta-block Hessian. `observed` selects the marginal observed information
// (E[I_c|y] minus rank-1 Var[N|y]); otherwise the PSD Fisher-scoring matrix.
inline MatrixXd dist_hess_beta(const Map<MatrixXd>& Xl, const Map<MatrixXd>& Xs,
                               int p_lam, int p_sig, bool hazard,
                               const std::vector<DistSiteResult>& st, bool observed) {
    const int n_sites = (int)st.size();
    const int pb = p_lam + p_sig + (hazard ? 1 : 0);
    const int sb = p_lam;                        // sigma block offset
    const int bb = p_lam + p_sig;                // shape param offset
    MatrixXd H = MatrixXd::Zero(pb, pb);
    for (int s = 0; s < n_sites; ++s) {
        const DistSiteResult& r = st[s];
        VectorXd xl = Xl.row(s).transpose();
        VectorXd xs = Xs.row(s).transpose();
        // lambda-lambda (complete-data Fisher; PSD in both modes)
        H.block(0, 0, p_lam, p_lam).selfadjointView<Eigen::Lower>()
            .rankUpdate(xl, r.info_eta_lambda);
        const double (*Id)[2] = observed ? r.info_eta_d : r.info_eta_d_fs;
        // sigma-sigma
        H.block(sb, sb, p_sig, p_sig).selfadjointView<Eigen::Lower>()
            .rankUpdate(xs, Id[0][0]);
        if (hazard) {
            // sigma-shape and shape-shape (shape design is the unit intercept)
            H.block(bb, sb, 1, p_sig) += Id[0][1] * xs.transpose();
            H(bb, bb) += Id[1][1];
        }
        if (observed && r.var_N > 0.0) {
            VectorXd u = VectorXd::Zero(pb);
            u.segment(0, p_lam)     = r.score_wt_lambda * xl;
            u.segment(sb, p_sig)    = r.vN_d[0] * xs;
            if (hazard) u(bb) = r.vN_d[1];
            H.selfadjointView<Eigen::Lower>().rankUpdate(u, -r.var_N);
        }
    }
    return H.selfadjointView<Eigen::Lower>();
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List cpp_distance_laplace_fixed(
    Rcpp::IntegerMatrix y,               // n_sites x n_bins bin counts
    Rcpp::NumericMatrix X_lambda_R,      // n_sites x p_lam
    Rcpp::NumericMatrix X_sigma_R,       // n_sites x p_sig
    Rcpp::NumericVector cutpoints,       // length n_bins + 1
    int transect,                        // 0 line, 1 point
    int key,                             // 0 half-normal, 1 hazard-rate
    Rcpp::NumericVector beta_lambda_init,
    Rcpp::NumericVector beta_sigma_init,
    double eta_b_init,                   // log-shape init (hazard only)
    int K_max, int max_iter, double tol, bool verbose,
    bool nb, double log_r_init, double theta_max, int quad_order
) {
    const int n_sites = y.nrow(), n_bins = y.ncol();
    const int p_lam = X_lambda_R.ncol(), p_sig = X_sigma_R.ncol();
    const bool hazard = (key == tulpaObs::DIST_HAZARD);
    if (X_lambda_R.nrow() != n_sites) Rcpp::stop("nrow(X_lambda) must equal nrow(y).");
    if (X_sigma_R.nrow()  != n_sites) Rcpp::stop("nrow(X_sigma) must equal nrow(y).");
    if ((int)cutpoints.size() != n_bins + 1) Rcpp::stop("length(cutpoints) must equal ncol(y) + 1.");
    if ((int)beta_lambda_init.size() != p_lam) Rcpp::stop("beta_lambda_init length must equal ncol(X_lambda).");
    if ((int)beta_sigma_init.size()  != p_sig) Rcpp::stop("beta_sigma_init length must equal ncol(X_sigma).");

    std::vector<double> cut(cutpoints.begin(), cutpoints.end());
    DistQuad quad = tulpaObs::dist_build_quad(cut, transect, quad_order);

    Map<MatrixXd> Xl(REAL(X_lambda_R), n_sites, p_lam);
    Map<MatrixXd> Xs(REAL(X_sigma_R), n_sites, p_sig);
    VectorXd beta_lam = Map<VectorXd>(REAL(beta_lambda_init), p_lam);
    VectorXd beta_sig = Map<VectorXd>(REAL(beta_sigma_init), p_sig);
    double eta_b = hazard ? eta_b_init : 0.0;

    const int pb = p_lam + p_sig + (hazard ? 1 : 0);
    const int p_total = pb + (nb ? 1 : 0);
    const double theta_min = std::log(kMinNbSize);
    double theta = nb ? std::min(std::max(log_r_init, theta_min), theta_max) : R_PosInf;
    double r = nb ? std::exp(theta) : std::numeric_limits<double>::infinity();

    std::vector<DistSiteResult> st(n_sites);

    // --- Inner Newton on beta at fixed r, run inside the theta profile loop. ---
    auto inner_newton = [&](double& log_lik, double& grad_norm, bool& converged) {
        log_lik = R_NegInf; grad_norm = R_PosInf; converged = false;
        for (int iter = 0; iter < max_iter; ++iter) {
            log_lik = dist_sweep(y, Xl, Xs, beta_lam, beta_sig, eta_b, key, quad,
                                 K_max, r, st);
            VectorXd grad = dist_grad_beta(Xl, Xs, p_lam, p_sig, hazard, st);
            grad_norm = grad.norm();
            if (verbose) Rcpp::Rcout << "    [beta] iter " << iter << "  log_lik "
                                     << log_lik << "  grad_norm " << grad_norm << "\n";
            if (grad_norm < tol) { converged = true; break; }

            MatrixXd H = dist_hess_beta(Xl, Xs, p_lam, p_sig, hazard, st, true);
            VectorXd delta;
            Eigen::LLT<MatrixXd> chol(H);
            if (chol.info() == Eigen::Success) {
                delta = chol.solve(grad);
            } else {
                MatrixXd Hf = dist_hess_beta(Xl, Xs, p_lam, p_sig, hazard, st, false);
                Eigen::LLT<MatrixXd> cf(Hf);
                if (cf.info() != Eigen::Success) {
                    Rcpp::warning("Fisher fallback Cholesky failure at beta iter %d.", iter);
                    break;
                }
                delta = cf.solve(grad);
                if (verbose) Rcpp::Rcout << "      (Fisher-scoring fallback)\n";
            }
            VectorXd dl = delta.segment(0, p_lam);
            VectorXd ds = delta.segment(p_lam, p_sig);
            double db = hazard ? delta(p_lam + p_sig) : 0.0;

            double step = 1.0; bool stepped = false;
            for (int h = 0; h < 12; ++h) {
                VectorXd bl = beta_lam + step * dl;
                VectorXd bs = beta_sig + step * ds;
                double eb = eta_b + step * db;
                double ll_try = dist_loglik(y, Xl, Xs, bl, bs, eb, key, quad, K_max, r);
                if (R_finite(ll_try) && ll_try >= log_lik - kLineSearchSlack) {
                    beta_lam = bl; beta_sig = bs; eta_b = eb; stepped = true; break;
                }
                step *= 0.5;
            }
            if (!stepped) { Rcpp::warning("beta step halving exhausted at iter %d.", iter); break; }
        }
        log_lik = dist_sweep(y, Xl, Xs, beta_lam, beta_sig, eta_b, key, quad,
                             K_max, r, st);
    };

    double log_lik = R_NegInf, beta_grad_norm = R_PosInf, grad_theta = 0.0;
    bool beta_conv = false, theta_converged = !nb, dispersion_boundary = false;
    int outer_iter = 0;
    const int outer_max = nb ? max_iter : 1;
    for (outer_iter = 0; outer_iter < outer_max; ++outer_iter) {
        inner_newton(log_lik, beta_grad_norm, beta_conv);
        if (!nb) break;
        grad_theta = 0.0; double info_theta_tot = 0.0;
        for (int s = 0; s < n_sites; ++s) {
            grad_theta     += st[s].grad_theta;
            info_theta_tot += st[s].info_theta - st[s].var_stheta;
        }
        if (verbose) Rcpp::Rcout << "  [theta] outer " << outer_iter << "  log_r " << theta
                                 << "  grad_theta " << grad_theta << "\n";
        if (std::abs(grad_theta) < tol && beta_conv) { theta_converged = true; break; }
        double dth = (info_theta_tot > 1e-8) ? grad_theta / info_theta_tot
                                             : (grad_theta > 0 ? 0.5 : -0.5);
        if (dth > 1.5) dth = 1.5; if (dth < -1.5) dth = -1.5;
        double ll_cur = log_lik, step = 1.0; bool stepped = false;
        for (int h = 0; h < 25; ++h) {
            double th_try = std::min(std::max(theta + step * dth, theta_min), theta_max);
            if (th_try == theta) break;
            double r_try = std::exp(th_try);
            double ll_try = dist_loglik(y, Xl, Xs, beta_lam, beta_sig, eta_b, key, quad, K_max, r_try);
            if (R_finite(ll_try) && ll_try >= ll_cur - kLineSearchSlack) {
                theta = th_try; r = r_try; stepped = true; break;
            }
            step *= 0.5;
        }
        if (!stepped) {
            dispersion_boundary = (theta >= theta_max - 1e-6);
            theta_converged = dispersion_boundary || (std::abs(grad_theta) < 1e-2);
            break;
        }
    }

    // --- Final joint observed-information Hessian at the mode. ---
    double log_lik_final = dist_sweep(y, Xl, Xs, beta_lam, beta_sig, eta_b, key, quad,
                                      K_max, r, st);
    MatrixXd H_beta = dist_hess_beta(Xl, Xs, p_lam, p_sig, hazard, st, true);
    MatrixXd H_obs = MatrixXd::Zero(p_total, p_total);
    H_obs.topLeftCorner(pb, pb) = H_beta;
    if (nb) {
        const int th = pb, sb = p_lam, bb = p_lam + p_sig;
        grad_theta = 0.0;
        for (int s = 0; s < n_sites; ++s) {
            const DistSiteResult& rr = st[s];
            grad_theta += rr.grad_theta;
            double w_lam_th = rr.info_lambda_theta - rr.cov_N_stheta * rr.score_wt_lambda;
            for (int k = 0; k < p_lam; ++k) H_obs(th, k) += w_lam_th * Xl(s, k);
            double w_sig_th = -rr.cov_N_stheta * rr.vN_d[0];
            for (int k = 0; k < p_sig; ++k) H_obs(th, sb + k) += w_sig_th * Xs(s, k);
            if (hazard) H_obs(th, bb) += -rr.cov_N_stheta * rr.vN_d[1];
            H_obs(th, th) += rr.info_theta - rr.var_stheta;
        }
        for (int j = 0; j < pb; ++j) H_obs(j, th) = H_obs(th, j);
    }

    MatrixXd vcov; bool vcov_ok = true;
    {
        Eigen::LLT<MatrixXd> chol(H_obs);
        if (chol.info() == Eigen::Success)
            vcov = chol.solve(MatrixXd::Identity(p_total, p_total));
        else { vcov_ok = false; vcov = MatrixXd::Constant(p_total, p_total, NA_REAL); }
    }

    VectorXd mean_N(n_sites), var_N(n_sites), bw(n_sites), p_det(n_sites);
    for (int s = 0; s < n_sites; ++s) {
        mean_N(s) = st[s].mean_N; var_N(s) = st[s].var_N;
        bw(s) = st[s].boundary_weight; p_det(s) = st[s].p_det;
    }

    Rcpp::NumericVector beta_lam_out(beta_lam.data(), beta_lam.data() + p_lam);
    Rcpp::NumericVector beta_sig_out(beta_sig.data(), beta_sig.data() + p_sig);
    bool converged = beta_conv && theta_converged;
    double grad_norm_report = (nb && !dispersion_boundary)
        ? std::sqrt(beta_grad_norm * beta_grad_norm + grad_theta * grad_theta)
        : beta_grad_norm;

    return Rcpp::List::create(
        Rcpp::Named("beta_lambda")         = beta_lam_out,
        Rcpp::Named("beta_sigma")          = beta_sig_out,
        Rcpp::Named("eta_b")               = hazard ? eta_b : NA_REAL,
        Rcpp::Named("shape")               = hazard ? std::exp(eta_b) : NA_REAL,
        Rcpp::Named("log_r")               = nb ? theta : NA_REAL,
        Rcpp::Named("r")                   = nb ? r : NA_REAL,
        Rcpp::Named("grad_theta")          = nb ? grad_theta : NA_REAL,
        Rcpp::Named("dispersion_boundary") = nb ? dispersion_boundary : false,
        Rcpp::Named("log_lik")             = log_lik_final,
        Rcpp::Named("vcov")                = Rcpp::wrap(vcov),
        Rcpp::Named("vcov_ok")             = vcov_ok,
        Rcpp::Named("H_obs")               = Rcpp::wrap(H_obs),
        Rcpp::Named("n_iter")              = outer_iter + (converged ? 1 : 0),
        Rcpp::Named("converged")           = converged,
        Rcpp::Named("grad_norm")           = grad_norm_report,
        Rcpp::Named("mean_N")              = Rcpp::wrap(mean_N),
        Rcpp::Named("var_N")               = Rcpp::wrap(var_N),
        Rcpp::Named("boundary_weight")     = Rcpp::wrap(bw),
        Rcpp::Named("p_det")               = Rcpp::wrap(p_det)
    );
}

// Total marginal log-likelihood + per-site gradients / observed-info pieces for
// the distance model -- the per-site-marginal surface the NUTS oracle and any
// random-effect integrator read (mirrors cpp_nmix_total_log_lik). The detection
// gradients are returned in eta space (d/d eta_sigma per site, and the summed
// d/d eta_b for the scalar hazard shape).
// [[Rcpp::export]]
Rcpp::List cpp_distance_total_log_lik(
    Rcpp::IntegerMatrix y,
    Rcpp::NumericVector eta_lambda,
    Rcpp::NumericVector eta_sigma,
    double eta_b,
    Rcpp::NumericVector cutpoints,
    int transect, int key, int K_max, double r, int quad_order
) {
    const int n_sites = y.nrow(), n_bins = y.ncol();
    if ((int)eta_lambda.size() != n_sites) Rcpp::stop("length(eta_lambda) must equal nrow(y).");
    if ((int)eta_sigma.size()  != n_sites) Rcpp::stop("length(eta_sigma) must equal nrow(y).");
    if ((int)cutpoints.size() != n_bins + 1) Rcpp::stop("length(cutpoints) must equal ncol(y) + 1.");
    const bool hazard = (key == tulpaObs::DIST_HAZARD);
    std::vector<double> cut(cutpoints.begin(), cutpoints.end());
    DistQuad quad = tulpaObs::dist_build_quad(cut, transect, quad_order);

    double total_ll = 0.0, total_grad_theta = 0.0, grad_eta_b = 0.0;
    Rcpp::NumericVector log_lik_site(n_sites), grad_eta_lambda(n_sites),
        grad_eta_sigma(n_sites), mean_N(n_sites), var_N(n_sites),
        boundary_weight(n_sites), p_det(n_sites);
    std::vector<int> y_site(n_bins);
    int n_K_inadmissible = 0;
    for (int s = 0; s < n_sites; ++s) {
        for (int b = 0; b < n_bins; ++b) y_site[b] = y(s, b);
        DistSiteResult res = compute_distance_site(y_site.data(), n_bins, eta_lambda[s],
                                                   eta_sigma[s], eta_b, key, quad, K_max, r);
        if (!R_finite(res.log_lik)) ++n_K_inadmissible;
        total_ll += res.log_lik; total_grad_theta += res.grad_theta;
        log_lik_site[s] = res.log_lik;
        grad_eta_lambda[s] = res.grad_eta_lambda;
        grad_eta_sigma[s]  = res.grad_eta_d[0];
        if (hazard) grad_eta_b += res.grad_eta_d[1];
        mean_N[s] = res.mean_N; var_N[s] = res.var_N;
        boundary_weight[s] = res.boundary_weight; p_det[s] = res.p_det;
    }
    return Rcpp::List::create(
        Rcpp::Named("log_lik")          = total_ll,
        Rcpp::Named("log_lik_site")     = log_lik_site,
        Rcpp::Named("grad_eta_lambda")  = grad_eta_lambda,
        Rcpp::Named("grad_eta_sigma")   = grad_eta_sigma,
        Rcpp::Named("grad_eta_b")       = grad_eta_b,
        Rcpp::Named("grad_theta")       = total_grad_theta,
        Rcpp::Named("mean_N")           = mean_N,
        Rcpp::Named("var_N")            = var_N,
        Rcpp::Named("boundary_weight")  = boundary_weight,
        Rcpp::Named("p_det")            = p_det,
        Rcpp::Named("n_K_inadmissible") = n_K_inadmissible
    );
}


// Per-site sweep for the areal-spatial distance fit (gcol33/tulpaObs#51): returns,
// for each site, the abundance-arm marginal moments and the half-normal detection
// arm's gradient / Fisher / N-coupling, all from compute_distance_site. The R
// inner Newton assembles the joint (beta_lambda, beta_sigma, z) gradient and the
// marginal observed-information Hessian (per-site eta-space block
// diag(info_lam, info_sig) - var_N v v', v = (score_wt_lambda, vN_d), the form
// distance_kernel.h documents) scattered through X_lambda / X_sigma + the field
// loading, plus the CAR prior. `key` selects the detection key (0 half-normal,
// 1 hazard-rate); under the hazard key the scalar log-shape `eta_b` is a global
// coordinate and the sweep returns the summed shape score grad_eta_b (and its
// summed Fisher information info_b) so the driver can fold it into the fixed
// block. Poisson or NB (the NB size r is the outer-grid axis).
// [[Rcpp::export]]
Rcpp::List cpp_distance_site_sweep(
    Rcpp::IntegerMatrix y_bins,
    Rcpp::NumericVector eta_lambda, Rcpp::NumericVector eta_sigma,
    Rcpp::NumericVector cutpoints, int transect, int quad_order, int K_max,
    bool nb, double r, int key = 0, double eta_b = 0.0
) {
    const int n_sites = y_bins.nrow(), n_bins = y_bins.ncol();
    if (eta_lambda.size() != n_sites || eta_sigma.size() != n_sites)
        Rcpp::stop("eta_lambda / eta_sigma must have length nrow(y_bins).");
    std::vector<double> cut(cutpoints.begin(), cutpoints.end());
    tulpaObs::DistQuad quad = tulpaObs::dist_build_quad(cut, transect, quad_order);
    const double rr = nb ? r : std::numeric_limits<double>::infinity();
    const int key_code = (key == 1) ? tulpaObs::DIST_HAZARD : tulpaObs::DIST_HALFNORMAL;
    const bool hazard = (key_code == tulpaObs::DIST_HAZARD);

    Rcpp::NumericVector logL(n_sites), grad_lam(n_sites), info_lam(n_sites),
        var_N(n_sites), swl(n_sites), grad_sig(n_sites), info_sig_obs(n_sites),
        info_sig_fs(n_sites), vN_sig(n_sites), boundary(n_sites), p_det(n_sites);
    std::vector<int> yb(n_bins);
    int n_inadmissible = 0;
    double grad_logr = 0.0;       // summed NB dispersion score (0 under Poisson)
    double grad_b = 0.0, info_b = 0.0;  // summed hazard log-shape score / Fisher info
    for (int s = 0; s < n_sites; ++s) {
        for (int b = 0; b < n_bins; ++b) yb[b] = y_bins(s, b);
        tulpaObs::DistSiteResult d = tulpaObs::compute_distance_site(
            yb.data(), n_bins, eta_lambda[s], eta_sigma[s], eta_b,
            key_code, quad, K_max, rr);
        if (!R_finite(d.log_lik)) ++n_inadmissible;
        logL[s] = d.log_lik;
        grad_lam[s] = d.grad_eta_lambda; info_lam[s] = d.info_eta_lambda;
        var_N[s] = d.var_N; swl[s] = d.score_wt_lambda;
        grad_sig[s] = d.grad_eta_d[0];
        info_sig_obs[s] = d.info_eta_d[0][0]; info_sig_fs[s] = d.info_eta_d_fs[0][0];
        vN_sig[s] = d.vN_d[0]; boundary[s] = d.boundary_weight; p_det[s] = d.p_det;
        if (nb) grad_logr += d.grad_theta;
        if (hazard) { grad_b += d.grad_eta_d[1]; info_b += d.info_eta_d[1][1]; }
    }
    return Rcpp::List::create(
        Rcpp::Named("log_lik") = logL,
        Rcpp::Named("grad_lam") = grad_lam, Rcpp::Named("info_lam") = info_lam,
        Rcpp::Named("var_N") = var_N, Rcpp::Named("swl") = swl,
        Rcpp::Named("grad_sig") = grad_sig,
        Rcpp::Named("info_sig_obs") = info_sig_obs,
        Rcpp::Named("info_sig_fs") = info_sig_fs,
        Rcpp::Named("vN_sig") = vN_sig, Rcpp::Named("boundary") = boundary,
        Rcpp::Named("p_det") = p_det, Rcpp::Named("grad_logr") = grad_logr,
        Rcpp::Named("grad_b") = grad_b, Rcpp::Named("info_b") = info_b,
        Rcpp::Named("n_inadmissible") = n_inadmissible);
}
