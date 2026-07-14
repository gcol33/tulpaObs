// Community relative-abundance (count) NUTS: the non-spatial community Poisson
// GLMM sampled over the exact joint posterior (community means, per-species
// coefficient deviations, and the community covariance), via tulpa's in-tree
// FullGradFn engine. This is the reduced counterpart of ms_abun_nuts.cpp: the
// community count has NO detection arm and NO latent-N marginalisation, so the
// per-(species, site) contribution is a plain Poisson log-likelihood
//   log p = sum_{s,i} [ y_{s,i} eta_{s,i} - exp(eta_{s,i}) - lgamma(y_{s,i}+1) ]
// with eta_{s,i} = X_i . (mu + b_s). NON-CENTERED: b_s = C z_s, z_s ~ N(0, I),
// so the community covariance (log-Cholesky C) enters only the data term. The
// arm Cholesky helpers are the shared community_chol.h. Byte-exact vs the R
// oracle .tobs_ms_count_nuts_logpost (gcol33/tulpaObs#117).

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/likelihood.h>
#include <tulpa/nuts_api.h>
#include "community_chol.h"

using namespace Rcpp;

namespace tulpaObs {

struct MsCountPri { double logdiag_mean = 0.0, logdiag_sd = 1.5, offdiag_sd = 1.0; };

struct MsCountNutsData {
    int n_sites = 0, n_species = 0, p = 0;
    std::vector<double> y;          // n_sites x n_species, column-major (species-major)
    std::vector<double> lgy;        // lgamma(y + 1), same layout
    NumericMatrix X;                // n_sites x p
    // layout: mu (p), z species-major (S*p), chol (p*(p+1)/2).
    int q = 0, total = 0, mu_off = 0, b_off = 0, chol_off = 0;
};

inline void ms_count_nuts_layout(MsCountNutsData& d) {
    d.q = d.p * (d.p + 1) / 2;
    d.mu_off = 0; d.b_off = d.p; d.chol_off = d.p + d.n_species * d.p;
    d.total = d.chol_off + d.q;
}

inline MsCountNutsData ms_count_nuts_build_data(const Rcpp::List& spec) {
    MsCountNutsData d;
    d.X = Rcpp::as<NumericMatrix>(spec["X"]);
    NumericMatrix ym = Rcpp::as<NumericMatrix>(spec["y"]);   // n_sites x n_species
    d.n_sites = d.X.nrow(); d.p = d.X.ncol(); d.n_species = ym.ncol();
    if (ym.nrow() != d.n_sites) Rcpp::stop("y must have n_sites rows.");
    d.y.assign((std::size_t) d.n_sites * d.n_species, 0.0);
    d.lgy.assign((std::size_t) d.n_sites * d.n_species, 0.0);
    for (int s = 0; s < d.n_species; ++s)
        for (int i = 0; i < d.n_sites; ++i) {
            const double yy = ym(i, s);
            d.y[(std::size_t) s * d.n_sites + i]  = yy;
            d.lgy[(std::size_t) s * d.n_sites + i] = std::lgamma(yy + 1.0);
        }
    ms_count_nuts_layout(d);
    return d;
}

// Joint log-posterior + gradient over theta = (mu, {z_s}, chol). NUTS maximises.
inline double ms_count_nuts_eval(const MsCountNutsData& d, const double* th,
                                 double sigma_beta, const MsCountPri& pr, double* g) {
    const int P = d.p, S = d.n_species, N = d.n_sites;
    for (int j = 0; j < d.total; ++j) g[j] = 0.0;
    const double* mu = th + d.mu_off;
    const double* z  = th + d.b_off;
    double* g_mu = g + d.mu_off;
    double* g_z  = g + d.b_off;

    std::vector<double> C;
    chol_unpack_cpp(th + d.chol_off, P, C);
    std::vector<double> A((std::size_t) P * P, 0.0);   // A[i,j] = sum_s grad_b_{s,i} z_{s,j}

    std::vector<double> gmu_s((std::size_t) S * P, 0.0), lp_s(S, 0.0);
    std::vector<double> A_s((std::size_t) S * P * P, 0.0);

    #pragma omp parallel for schedule(static)
    for (int s = 0; s < S; ++s) {
        const double* z_s = z + s * P;
        std::vector<double> b(P), gb(P, 0.0);
        for (int i = 0; i < P; ++i) {                  // b = C z (lower-tri C)
            double v = 0.0;
            for (int j = 0; j <= i; ++j) v += C[(std::size_t) i * P + j] * z_s[j];
            b[i] = v;
        }
        const double* ys = &d.y[(std::size_t) s * N];
        const double* lg = &d.lgy[(std::size_t) s * N];
        double* gmu_loc = &gmu_s[(std::size_t) s * P];
        double lp_loc = 0.0;
        for (int i = 0; i < N; ++i) {
            double eta = 0.0;
            for (int k = 0; k < P; ++k) eta += d.X(i, k) * (mu[k] + b[k]);
            const double lam = std::exp(eta < 700.0 ? eta : 700.0);
            lp_loc += ys[i] * eta - lam - lg[i];
            const double ge = ys[i] - lam;             // d log p / d eta
            for (int k = 0; k < P; ++k) {
                const double gx = ge * d.X(i, k);
                gmu_loc[k] += gx; gb[k] += gx;
            }
        }
        lp_s[s] = lp_loc;
        double* gz_s = g_z + s * P;                    // z grad (data) = C' gb
        for (int v = 0; v < P; ++v) {
            double sg = 0.0;
            for (int i = v; i < P; ++i) sg += C[(std::size_t) i * P + v] * gb[i];
            gz_s[v] += sg;
        }
        double* As = &A_s[(std::size_t) s * P * P];    // A_s = grad_b z'
        for (int i = 0; i < P; ++i)
            for (int j = 0; j <= i; ++j) As[(std::size_t) i * P + j] = gb[i] * z_s[j];
    }

    double lp = 0.0;                                   // serial reduction (species order)
    for (int s = 0; s < S; ++s) {
        const double* gmu_loc = &gmu_s[(std::size_t) s * P];
        for (int k = 0; k < P; ++k) g_mu[k] += gmu_loc[k];
        lp += lp_s[s];
        const double* As = &A_s[(std::size_t) s * P * P];
        for (int t = 0; t < P * P; ++t) A[t] += As[t];
    }

    for (int j = 0; j < S * P; ++j) {                  // z prior N(0, I)
        const double zz = z[j]; g_z[j] -= zz; lp += -0.5 * zz * zz;
    }
    lp += chol_data_grad_noncentered(A, C, P, th + d.chol_off,
                                     pr.logdiag_mean, pr.logdiag_sd, pr.offdiag_sd,
                                     g + d.chol_off);
    const double inv_sb2 = 1.0 / (sigma_beta * sigma_beta);
    for (int k = 0; k < P; ++k) {                      // community-mean prior
        g_mu[k] -= inv_sb2 * mu[k]; lp += -0.5 * inv_sb2 * mu[k] * mu[k];
    }
    return lp;
}

struct MsCountNutsModel { MsCountNutsData d; double sigma_beta = 10.0; MsCountPri pr; };

inline void ms_count_nuts_full_grad(const std::vector<double>& params,
                                    const tulpa::ModelData& data,
                                    const tulpa::ParamLayout& /*layout*/,
                                    std::vector<double>& grad, double* log_post_out) {
    const MsCountNutsModel* m =
        static_cast<const MsCountNutsModel*>(data.model_response_data);
    grad.assign((std::size_t) m->d.total, 0.0);
    const double lp = ms_count_nuts_eval(m->d, params.data(), m->sigma_beta,
                                         m->pr, grad.data());
    if (log_post_out) *log_post_out = lp;
}

inline MsCountPri ms_count_pri_from_list(const Rcpp::List& pri) {
    MsCountPri pr;
    pr.logdiag_mean = Rcpp::as<double>(pri["chol_logdiag_mean"]);
    pr.logdiag_sd   = Rcpp::as<double>(pri["chol_logdiag_sd"]);
    pr.offdiag_sd   = Rcpp::as<double>(pri["chol_offdiag_sd"]);
    return pr;
}

} // namespace tulpaObs

// [[Rcpp::export]]
Rcpp::List cpp_ms_count_nuts_joint_logpost(Rcpp::List spec, Rcpp::NumericVector theta,
                                           Rcpp::List pri, double sigma_beta) {
    tulpaObs::MsCountNutsData d = tulpaObs::ms_count_nuts_build_data(spec);
    if ((int) theta.size() != d.total)
        Rcpp::stop("theta length %d != expected %d", (int) theta.size(), d.total);
    tulpaObs::MsCountPri pr = tulpaObs::ms_count_pri_from_list(pri);
    Rcpp::NumericVector grad(d.total);
    const double lp = tulpaObs::ms_count_nuts_eval(
        d, theta.begin(), sigma_beta, pr, grad.begin());
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// [[Rcpp::export]]
Rcpp::List cpp_ms_count_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                             Rcpp::List pri, double sigma_beta,
                             Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                             int n_iter, int n_warmup, int max_treedepth,
                             double adapt_delta, int seed, bool verbose) {
    tulpaObs::MsCountNutsModel m;
    m.d = tulpaObs::ms_count_nuts_build_data(spec);
    m.sigma_beta = sigma_beta;
    m.pr = tulpaObs::ms_count_pri_from_list(pri);
    if ((int) theta0.size() != m.d.total)
        Rcpp::stop("theta0 length %d != expected %d", (int) theta0.size(), m.d.total);

    tulpa::LikelihoodSpec lspec;
    lspec.name = "ms_count";
    lspec.n_processes = 1;
    lspec.gradient_fn = &tulpaObs::ms_count_nuts_full_grad;

    tulpa::ModelData data;
    data.N = m.d.n_sites;
    data.n_processes = 1;
    data.sigma_beta = sigma_beta;
    data.model_response_data = &m;
    data.likelihood_spec = &lspec;
    data.sharing.init(1);
    data.zi_type = tulpa::ZIType::NONE;
    data.p_zi = 0; data.p_oi = 0;

    tulpa::ParamLayout layout;
    layout.total_params = m.d.total;
    tulpa::set_gradient_mode_str("H");

    std::vector<double> init(theta0.begin(), theta0.end());
    std::vector<double> imv; const double* im = nullptr;
    if (inv_metric.isNotNull()) {
        Rcpp::NumericVector v(inv_metric); imv.assign(v.begin(), v.end()); im = imv.data();
    }

    tulpa::NUTSFn run_nuts = tulpa::get_nuts_fn();
    tulpa::NUTSResult result = {};
    run_nuts(&data, &layout, init.data(), m.d.total, n_iter, n_warmup,
             max_treedepth, adapt_delta, static_cast<unsigned int>(seed),
             verbose ? 1 : 0, im, &result);

    const int n_samples = result.n_sample, np = m.d.total;
    Rcpp::NumericMatrix draws(n_samples, np);
    Rcpp::NumericVector lp(n_samples), ap(n_samples);
    Rcpp::IntegerVector div(n_samples), td(n_samples);
    for (int s = 0; s < n_samples; ++s) {
        for (int j = 0; j < np; ++j) draws(s, j) = result.samples[s * np + j];
        lp[s] = result.log_prob[s]; ap[s] = result.accept_prob[s];
        div[s] = result.divergent[s]; td[s] = result.treedepth[s];
    }
    const double epsilon = result.epsilon;
    result.free_buffers();
    return Rcpp::List::create(
        Rcpp::Named("draws") = draws, Rcpp::Named("log_prob") = lp,
        Rcpp::Named("accept_prob") = ap, Rcpp::Named("divergent") = div,
        Rcpp::Named("treedepth") = td, Rcpp::Named("epsilon") = epsilon,
        Rcpp::Named("n_params") = np);
}
