// ms_occu_cover_spatial_nuts.cpp
// C++ marginal log-likelihood + (later) gradient for the reduced-rank
// spatial-factor community occu_cover NUTS target (gcol33/tulpa#67). The R
// reference .ms_ocs_joint_logpost (R/ms_occu_cover_spatial_nuts.R) is the oracle;
// this port mirrors it block by block and is cross-checked against it.
//
// Step 1 (this commit): the per-species per-cell occu_cover marginal log-lik with
// the shared latent field injected on the occupancy predictor (and, with a
// cover-arm factor, on the cover predictor). Cover granularity is per-visit
// ("none"): the cover term is one positive-arm log-density per detected visit.
// The occupancy state z is integrated out in closed form over its two states,
// reusing nodet_mixture_block / LognormalPositive / BetaPositive from
// occu_coupling_shared.h so the likelihood stays on one source of truth.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "occu_coupling_shared.h"

using namespace Rcpp;

namespace tulpaObs {

// Marshalled per-fit data for the spatial-factor community occu_cover marginal.
// y / y_pos / valid are length-S lists of n_sites x max_visits matrices (one per
// species); the designs are cell-level (site rows), broadcast across visits.
struct MsOcsData {
    int n_sites = 0, max_visits = 0, S = 0, K = 0;
    int P_occ = 0, P_p = 0, P_pos = 0;
    bool cover_factor = false;
    bool is_beta = false;
    std::vector<NumericMatrix> y, y_pos, valid;   // length S
    NumericMatrix X_occ, X_p, X_pos;              // n_sites x P_arm
    NumericMatrix W;                              // n_sites x K (filled per eval)
};

inline MsOcsData ms_ocs_build_data(const List& spec) {
    MsOcsData d;
    d.n_sites      = as<int>(spec["n_sites"]);
    d.max_visits   = as<int>(spec["max_visits"]);
    d.S            = as<int>(spec["S"]);
    d.K            = as<int>(spec["K"]);
    d.P_occ        = as<int>(spec["P_occ"]);
    d.P_p          = as<int>(spec["P_p"]);
    d.P_pos        = as<int>(spec["P_pos"]);
    d.cover_factor = as<bool>(spec["cover_factor"]);
    d.is_beta      = as<bool>(spec["is_beta"]);
    d.X_occ = as<NumericMatrix>(spec["X_occ"]);
    d.X_p   = as<NumericMatrix>(spec["X_p"]);
    d.X_pos = as<NumericMatrix>(spec["X_pos"]);
    List yl = spec["y"], ypl = spec["y_pos"], vl = spec["valid"];
    d.y.reserve(d.S); d.y_pos.reserve(d.S); d.valid.reserve(d.S);
    for (int s = 0; s < d.S; ++s) {
        d.y.push_back(as<NumericMatrix>(yl[s]));
        d.y_pos.push_back(as<NumericMatrix>(ypl[s]));
        d.valid.push_back(as<NumericMatrix>(vl[s]));
    }
    return d;
}

inline double clamp30(double e) { return e < -30.0 ? -30.0 : (e > 30.0 ? 30.0 : e); }

// X (n x p, column-major) times the length-p coefficient vector at row i.
inline double row_dot(const NumericMatrix& X, int i, const double* beta, int p) {
    double acc = 0.0;
    for (int j = 0; j < p; ++j) acc += X(i, j) * beta[j];
    return acc;
}

// Per-species marginal log-lik summed over cells, for species s, given its
// arm coefficient triples (mu+b on each arm), the per-species field loadings
// L_s (length K) on occupancy and Lpos_s on cover (or nullptr), the field W
// (n_sites x K), and the dispersion. Mirrors .occu_cover_site_ll with the field
// offset injected on psi (and ep when a cover factor is present).
inline double ms_ocs_species_ll(const MsOcsData& d, int s,
                                const double* th_occ, const double* th_p,
                                const double* th_pos, const double* L_s,
                                const double* Lpos_s, double log_disp) {
    const double disp = std::exp(log_disp);
    const NumericMatrix& y  = d.y[s];
    const NumericMatrix& yp = d.y_pos[s];
    const NumericMatrix& vv = d.valid[s];
    double ll = 0.0;
    std::vector<double> eta_p_buf(d.max_visits);
    for (int i = 0; i < d.n_sites; ++i) {
        // occupancy predictor + shared field offset
        double eta_occ = row_dot(d.X_occ, i, th_occ, d.P_occ);
        for (int k = 0; k < d.K; ++k) eta_occ += d.W(i, k) * L_s[k];
        const double psi = sigmoid_(clamp30(eta_occ));

        // detection / cover predictors (cell-level, constant across visits)
        const double eta_p   = row_dot(d.X_p,   i, th_p,   d.P_p);
        double eta_pos       = row_dot(d.X_pos, i, th_pos, d.P_pos);
        if (d.cover_factor)
            for (int k = 0; k < d.K; ++k) eta_pos += d.W(i, k) * Lpos_s[k];
        const double p_i = sigmoid_(clamp30(eta_p));

        bool any_det = false;
        for (int v = 0; v < d.max_visits; ++v)
            if (vv(i, v) > 0.5 && y(i, v) > 0.5) { any_det = true; break; }

        if (any_det) {
            ll += log_safe_(psi);
            for (int v = 0; v < d.max_visits; ++v) {
                if (vv(i, v) < 0.5) continue;
                ll += (y(i, v) > 0.5) ? log_safe_(p_i) : log_safe_(1.0 - p_i);
                if (y(i, v) > 0.5) {
                    ll += d.is_beta
                        ? BetaPositive::log_density(yp(i, v), eta_pos, disp)
                        : LognormalPositive::log_density(yp(i, v), eta_pos, disp);
                }
            }
        } else {
            // no-detection occupancy mixture L = psi P0 + (1 - psi)
            int nv = 0;
            for (int v = 0; v < d.max_visits; ++v)
                if (vv(i, v) > 0.5) eta_p_buf[nv++] = eta_p;
            double g_w = 0.0, nh_w = 0.0;
            std::vector<double> g_p(nv > 0 ? nv : 1), nh_p(nv > 0 ? nv : 1);
            ll += nodet_mixture_block(psi, eta_p_buf.data(), nv, false, false,
                                      g_w, nh_w, g_p.data(), nh_p.data(),
                                      nullptr, nullptr);
        }
    }
    return ll;
}

} // namespace tulpaObs

// Total marginal log-likelihood (data term only, no priors) for the spatial
// community occu_cover model at the packed inner latent theta_inner =
// c(mu[P], b[S*P] species-major, vec(L)[S*K], [vec(Lpos)[S*K]], vec(W)[N*K],
// log_disp). Unconstrained loadings (constrain = FALSE). Cross-check entry for
// the R oracle sum_s sum_cells .occu_cover_site_ll.
// [[Rcpp::export]]
double cpp_ms_ocs_marginal_ll(Rcpp::List spec, Rcpp::NumericVector theta_inner) {
    tulpaObs::MsOcsData d = tulpaObs::ms_ocs_build_data(spec);
    const int P = d.P_occ + d.P_p + d.P_pos;
    const int N = d.n_sites, S = d.S, K = d.K;
    const double* th = theta_inner.begin();

    int off = 0;
    const double* mu = th + off; off += P;
    const double* b  = th + off; off += S * P;          // species-major, length P each
    const double* Lv = th + off; off += S * K;          // vec(L), column-major (S x K)
    const double* Lpv = nullptr;
    if (d.cover_factor) { Lpv = th + off; off += S * K; }
    // fill W (n_sites x K) from vec(W)
    d.W = Rcpp::NumericMatrix(N, K);
    for (int k = 0; k < K; ++k)
        for (int i = 0; i < N; ++i) d.W(i, k) = th[off + k * N + i];
    off += N * K;
    const double log_disp = th[off];

    std::vector<double> th_occ(d.P_occ), th_p(d.P_p), th_pos(d.P_pos);
    std::vector<double> L_s(K), Lpos_s(K);
    double total = 0.0;
    for (int s = 0; s < S; ++s) {
        const double* b_s = b + s * P;
        for (int j = 0; j < d.P_occ; ++j) th_occ[j] = mu[j] + b_s[j];
        for (int j = 0; j < d.P_p; ++j)   th_p[j]   = mu[d.P_occ + j] + b_s[d.P_occ + j];
        for (int j = 0; j < d.P_pos; ++j) th_pos[j] = mu[d.P_occ + d.P_p + j] +
                                                      b_s[d.P_occ + d.P_p + j];
        for (int k = 0; k < K; ++k) L_s[k] = Lv[k * S + s];   // L(s,k), col-major
        if (d.cover_factor)
            for (int k = 0; k < K; ++k) Lpos_s[k] = Lpv[k * S + s];
        total += tulpaObs::ms_ocs_species_ll(
            d, s, th_occ.data(), th_p.data(), th_pos.data(),
            L_s.data(), d.cover_factor ? Lpos_s.data() : nullptr, log_disp);
    }
    return total;
}
