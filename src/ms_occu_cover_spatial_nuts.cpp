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
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/likelihood.h>
#include <tulpa/nuts_api.h>
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

    // Areal field structure: 0 icar (Q), 1 car_proper (D - h A), 2 bym2
    // (V diag(1/((1-h)+h s)) V'). has_hyper -> a per-factor logit_h coordinate.
    int field_type = 0;
    bool has_hyper = false;
    NumericMatrix Q, A, V;                         // icar Q; car A; bym2 V
    std::vector<double> deg, gamma, s;             // car deg/gamma; bym2 s
    double logdetD = 0.0;                          // car log|D|
    double field_rank = 0.0;                       // N-1 (icar) or N (proper)
};

// Shared field W is passed per call as a column-major n_sites x K buffer
// (W[k * n_sites + i]) rather than held on the data struct, so the gradient is
// re-entrant (the FullGradFn may run concurrently across NUTS chains).
inline double Wik(const double* W, int N, int i, int k) {
    return W[(std::size_t) k * N + i];
}

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
    // Field structure
    if (spec.containsElementNamed("Q")) d.Q = as<NumericMatrix>(spec["Q"]);
    if (spec.containsElementNamed("field_rank"))
        d.field_rank = as<double>(spec["field_rank"]);
    d.has_hyper = spec.containsElementNamed("has_hyper")
                  && as<bool>(spec["has_hyper"]);
    std::string ft = spec.containsElementNamed("field_type")
                     ? as<std::string>(spec["field_type"]) : "icar";
    if (ft == "car_proper") {
        d.field_type = 1;
        d.A = as<NumericMatrix>(spec["A"]);
        d.deg = as<std::vector<double>>(spec["deg"]);
        d.gamma = as<std::vector<double>>(spec["gamma"]);
        d.logdetD = as<double>(spec["logdetD"]);
    } else if (ft == "bym2") {
        d.field_type = 2;
        d.V = as<NumericMatrix>(spec["V"]);
        d.s = as<std::vector<double>>(spec["s"]);
    } else {
        d.field_type = 0;
    }
    return d;
}

// ---- Areal field structure R(h): products, quadratic, log-determinant, and
// hyperparameter derivatives. h is ignored for icar. ----

// R(h) W_k -> out (length N).
inline void field_RWk(const MsOcsData& d, double h, const double* Wk, double* out) {
    const int N = d.n_sites;
    if (d.field_type == 1) {                       // car_proper: D - h A
        for (int i = 0; i < N; ++i) {
            double aw = 0.0;
            for (int j = 0; j < N; ++j) aw += d.A(i, j) * Wk[j];
            out[i] = d.deg[i] * Wk[i] - h * aw;
        }
    } else if (d.field_type == 2) {                // bym2: V diag(1/m) V'
        std::vector<double> t(N, 0.0);
        for (int j = 0; j < N; ++j) {
            double acc = 0.0;
            for (int i = 0; i < N; ++i) acc += d.V(i, j) * Wk[i];
            t[j] = acc / ((1.0 - h) + h * d.s[j]);
        }
        for (int i = 0; i < N; ++i) {
            double acc = 0.0;
            for (int j = 0; j < N; ++j) acc += d.V(i, j) * t[j];
            out[i] = acc;
        }
    } else {                                        // icar: Q
        for (int i = 0; i < N; ++i) {
            double acc = 0.0;
            for (int j = 0; j < N; ++j) acc += d.Q(i, j) * Wk[j];
            out[i] = acc;
        }
    }
}

// log|R(h)| for the proper fields (constant for icar -> returns 0, dropped).
inline double field_logdetR(const MsOcsData& d, double h) {
    const int N = d.n_sites;
    if (d.field_type == 1) {
        double acc = d.logdetD;
        for (int j = 0; j < N; ++j) acc += std::log(1.0 - h * d.gamma[j]);
        return acc;
    }
    if (d.field_type == 2) {
        double acc = 0.0;
        for (int j = 0; j < N; ++j) acc -= std::log((1.0 - h) + h * d.s[j]);
        return acc;
    }
    return 0.0;
}

// dquad = W_k'(dR/dh)W_k and dlogdet = d log|R(h)|/dh (proper fields only).
inline void field_hyper_terms(const MsOcsData& d, double h, const double* Wk,
                              double& dquad, double& dlogdet) {
    const int N = d.n_sites;
    if (d.field_type == 1) {                       // dR/dh = -A
        double aw_sum = 0.0;
        for (int i = 0; i < N; ++i) {
            double aw = 0.0;
            for (int j = 0; j < N; ++j) aw += d.A(i, j) * Wk[j];
            aw_sum += Wk[i] * aw;
        }
        dquad = -aw_sum;
        double dl = 0.0;
        for (int j = 0; j < N; ++j) dl -= d.gamma[j] / (1.0 - h * d.gamma[j]);
        dlogdet = dl;
    } else if (d.field_type == 2) {
        double dq = 0.0, dl = 0.0;
        for (int j = 0; j < N; ++j) {
            double t = 0.0;
            for (int i = 0; i < N; ++i) t += d.V(i, j) * Wk[i];
            const double m = (1.0 - h) + h * d.s[j];
            dq -= (t * t) * (d.s[j] - 1.0) / (m * m);
            dl -= (d.s[j] - 1.0) / m;
        }
        dquad = dq; dlogdet = dl;
    } else { dquad = 0.0; dlogdet = 0.0; }
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
                                const double* Lpos_s, double log_disp,
                                const double* W) {
    const double disp = std::exp(log_disp);
    const int N = d.n_sites;
    const NumericMatrix& y  = d.y[s];
    const NumericMatrix& yp = d.y_pos[s];
    const NumericMatrix& vv = d.valid[s];
    double ll = 0.0;
    std::vector<double> eta_p_buf(d.max_visits);
    for (int i = 0; i < d.n_sites; ++i) {
        // occupancy predictor + shared field offset
        double eta_occ = row_dot(d.X_occ, i, th_occ, d.P_occ);
        for (int k = 0; k < d.K; ++k) eta_occ += Wik(W, N, i, k) * L_s[k];
        const double psi = sigmoid_(clamp30(eta_occ));

        // detection / cover predictors (cell-level, constant across visits)
        const double eta_p   = row_dot(d.X_p,   i, th_p,   d.P_p);
        double eta_pos       = row_dot(d.X_pos, i, th_pos, d.P_pos);
        if (d.cover_factor)
            for (int k = 0; k < d.K; ++k) eta_pos += Wik(W, N, i, k) * Lpos_s[k];
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

// Accumulate species s's data-log-lik gradient into the per-arm coefficient
// score (g_occ_s / g_p_rowsum-chained / g_pos_rowsum-chained), the per-species
// loading scores (g_L_s, g_Lpos_s, length K), the shared-field score (g_W,
// n_sites x K, accumulated across species), and the dispersion score. Mirrors
// .occu_cover_eta_grad + the chain in .ms_ocs_penll_grad. `g_occ_s/g_p_s/g_pos_s`
// are length P_occ/P_p/P_pos outputs for THIS species (mu and b_s share them).
inline double ms_ocs_species_grad(const MsOcsData& d, int s,
                                  const double* th_occ, const double* th_p,
                                  const double* th_pos, const double* L_s,
                                  const double* Lpos_s, double log_disp,
                                  double* g_occ_s, double* g_p_s, double* g_pos_s,
                                  double* g_L_s, double* g_Lpos_s,
                                  double* g_W, double& g_ld, const double* W) {
    const double disp   = std::exp(log_disp);
    const double sigma  = disp;
    const double inv_s2 = (sigma > 0.0) ? 1.0 / (sigma * sigma) : 0.0;
    const NumericMatrix& y  = d.y[s];
    const NumericMatrix& yp = d.y_pos[s];
    const NumericMatrix& vv = d.valid[s];
    const int N = d.n_sites, J = d.max_visits, K = d.K;

    for (int j = 0; j < d.P_occ; ++j) g_occ_s[j] = 0.0;
    for (int j = 0; j < d.P_p;   ++j) g_p_s[j]   = 0.0;
    for (int j = 0; j < d.P_pos; ++j) g_pos_s[j] = 0.0;
    for (int k = 0; k < K; ++k) { g_L_s[k] = 0.0; if (g_Lpos_s) g_Lpos_s[k] = 0.0; }
    double ll = 0.0;

    for (int i = 0; i < N; ++i) {
        double eta_occ = row_dot(d.X_occ, i, th_occ, d.P_occ);
        for (int k = 0; k < K; ++k) eta_occ += Wik(W, N, i, k) * L_s[k];
        const double psi = sigmoid_(clamp30(eta_occ));
        const double eta_p = row_dot(d.X_p, i, th_p, d.P_p);
        double eta_pos     = row_dot(d.X_pos, i, th_pos, d.P_pos);
        if (d.cover_factor)
            for (int k = 0; k < K; ++k) eta_pos += Wik(W, N, i, k) * Lpos_s[k];
        const double p_i = sigmoid_(clamp30(eta_p));

        bool any_det = false; int n_valid = 0;
        for (int v = 0; v < J; ++v) {
            if (vv(i, v) < 0.5) continue;
            ++n_valid;
            if (y(i, v) > 0.5) any_det = true;
        }

        double g_psi = 0.0, gp_rowsum = 0.0, gpos_rowsum = 0.0;
        if (any_det) {
            g_psi = 1.0 - psi;
            ll += log_safe_(psi);
            for (int v = 0; v < J; ++v) {
                if (vv(i, v) < 0.5) continue;
                gp_rowsum += (y(i, v) > 0.5) ? (1.0 - p_i) : (-p_i);
                ll += (y(i, v) > 0.5) ? log_safe_(p_i) : log_safe_(1.0 - p_i);
                if (y(i, v) > 0.5) {
                    ll += d.is_beta
                        ? BetaPositive::log_density(yp(i, v), eta_pos, disp)
                        : LognormalPositive::log_density(yp(i, v), eta_pos, disp);
                    if (d.is_beta) {
                        gpos_rowsum += BetaPositive::grad_eta(yp(i, v), eta_pos, disp);
                        // dispersion score, beta (per detected visit)
                        const double mu = sigmoid_(clamp30(eta_pos));
                        const double a = mu * disp, b = (1.0 - mu) * disp;
                        const double ly = log_safe_(yp(i, v));
                        const double l1my = log_safe_(1.0 - yp(i, v));
                        g_ld += disp * (tulpa::math::portable_digamma(disp)
                                        - tulpa::math::portable_digamma(a) * mu
                                        - tulpa::math::portable_digamma(b) * (1.0 - mu)
                                        + mu * ly + (1.0 - mu) * l1my);
                    } else {
                        const double r = (log_safe_(yp(i, v)) - eta_pos) / sigma;
                        gpos_rowsum += r / sigma;       // (log y - eta) / sigma^2
                        g_ld += r * r - 1.0;
                    }
                }
            }
        } else {
            double log_P0 = (double) n_valid * log_safe_(1.0 - p_i);
            const double P0 = std::exp(log_P0);
            const double A = psi * P0, L = A + (1.0 - psi);
            const double invL = (L > 0.0) ? 1.0 / L : 0.0;
            g_psi = psi * (1.0 - psi) * (P0 - 1.0) * invL;
            gp_rowsum = -(A * invL) * p_i * (double) n_valid;
            ll += log_safe_(L);
        }

        // chain to coefficients / loadings / field
        for (int j = 0; j < d.P_occ; ++j) g_occ_s[j] += d.X_occ(i, j) * g_psi;
        for (int j = 0; j < d.P_p;   ++j) g_p_s[j]   += d.X_p(i, j)   * gp_rowsum;
        for (int j = 0; j < d.P_pos; ++j) g_pos_s[j] += d.X_pos(i, j) * gpos_rowsum;
        for (int k = 0; k < K; ++k) {
            g_L_s[k] += Wik(W, N, i, k) * g_psi;
            g_W[k * N + i] += g_psi * L_s[k];
            if (d.cover_factor) {
                g_Lpos_s[k] += Wik(W, N, i, k) * gpos_rowsum;
                g_W[k * N + i] += gpos_rowsum * Lpos_s[k];
            }
        }
    }
    (void) inv_s2;
    return ll;
}

// ---------------------------------------------------------------------------
// Small dense linear algebra on a community covariance Cholesky factor (P_arm is
// tiny: the number of coefficients on one arm). Matrices are row-major P*P.
// ---------------------------------------------------------------------------

// Packed column-major lower-triangle (diagonal on the log scale) -> lower
// Cholesky factor C (row-major), matching .ms_ocs_chol_unpack.
inline void chol_unpack_cpp(const double* vec, int P, std::vector<double>& C) {
    C.assign((std::size_t) P * P, 0.0);
    int pos = 0;
    for (int j = 0; j < P; ++j) {
        C[(std::size_t) j * P + j] = std::exp(vec[pos++]);
        for (int i = j + 1; i < P; ++i) C[(std::size_t) i * P + j] = vec[pos++];
    }
}

// Inverse of a lower-triangular matrix (row-major), by forward substitution.
inline void lower_tri_inv(const std::vector<double>& C, int P,
                          std::vector<double>& M) {
    M.assign((std::size_t) P * P, 0.0);
    for (int j = 0; j < P; ++j) {
        M[(std::size_t) j * P + j] = 1.0 / C[(std::size_t) j * P + j];
        for (int i = j + 1; i < P; ++i) {
            double s = 0.0;
            for (int k = j; k < i; ++k)
                s += C[(std::size_t) i * P + k] * M[(std::size_t) k * P + j];
            M[(std::size_t) i * P + j] = -s / C[(std::size_t) i * P + i];
        }
    }
}

// Sigma^{-1} = Cinv' Cinv (Cinv lower-tri, row-major) -> row-major P*P.
inline void sinv_from_cinv(const std::vector<double>& Mi, int P,
                           std::vector<double>& Si) {
    Si.assign((std::size_t) P * P, 0.0);
    for (int a = 0; a < P; ++a)
        for (int b = 0; b < P; ++b) {
            double s = 0.0;
            for (int k = 0; k < P; ++k)
                s += Mi[(std::size_t) k * P + a] * Mi[(std::size_t) k * P + b];
            Si[(std::size_t) a * P + b] = s;
        }
}

// Cholesky-coordinate gradient of T = -0.5 tr(Sigma^{-1} M) - 0.5 S log|Sigma|
// plus the coordinate hyperprior, written into `out` (packed column-major lower
// triangle). Mirrors .ms_ocs_chol_block_grad: G = 0.5 Si M Si - 0.5 S Si,
// dC = 2 G C, with the log-link chain on the diagonal.
inline void chol_block_grad_cpp(const std::vector<double>& C,
                                const std::vector<double>& Si,
                                const std::vector<double>& M, int P, double S,
                                const double* vec, double logdiag_mean,
                                double logdiag_sd, double offdiag_sd,
                                double* out) {
    std::vector<double> SM((std::size_t) P * P, 0.0), G((std::size_t) P * P, 0.0),
                        dC((std::size_t) P * P, 0.0);
    for (int a = 0; a < P; ++a)
        for (int b = 0; b < P; ++b) {
            double s = 0.0;
            for (int k = 0; k < P; ++k) s += Si[a * P + k] * M[k * P + b];
            SM[a * P + b] = s;
        }
    for (int a = 0; a < P; ++a)
        for (int b = 0; b < P; ++b) {
            double s = 0.0;
            for (int k = 0; k < P; ++k) s += SM[a * P + k] * Si[k * P + b];
            G[a * P + b] = 0.5 * s - 0.5 * S * Si[a * P + b];
        }
    for (int a = 0; a < P; ++a)
        for (int b = 0; b < P; ++b) {
            double s = 0.0;
            for (int k = 0; k < P; ++k) s += G[a * P + k] * C[k * P + b];
            dC[a * P + b] = 2.0 * s;
        }
    int pos = 0;
    for (int j = 0; j < P; ++j) {
        const double cjj = C[(std::size_t) j * P + j];
        out[pos] = dC[(std::size_t) j * P + j] * cjj
                 - (vec[pos] - logdiag_mean) / (logdiag_sd * logdiag_sd);
        ++pos;
        for (int i = j + 1; i < P; ++i) {
            out[pos] = dC[(std::size_t) i * P + j]
                     - vec[pos] / (offdiag_sd * offdiag_sd);
            ++pos;
        }
    }
}

} // namespace tulpaObs

// Data-log-lik gradient (no priors) over the packed inner latent, same layout as
// cpp_ms_ocs_marginal_ll's theta_inner. Cross-check entry for the R oracle
// .ms_ocs_penll_grad with the prior precisions zeroed.
// [[Rcpp::export]]
Rcpp::NumericVector cpp_ms_ocs_marginal_grad(Rcpp::List spec,
                                             Rcpp::NumericVector theta_inner) {
    tulpaObs::MsOcsData d = tulpaObs::ms_ocs_build_data(spec);
    const int P = d.P_occ + d.P_p + d.P_pos;
    const int N = d.n_sites, S = d.S, K = d.K;
    const double* th = theta_inner.begin();

    int off = 0;
    const double* mu = th + off; off += P;
    const double* b  = th + off; off += S * P;
    const double* Lv = th + off; off += S * K;
    const double* Lpv = nullptr;
    if (d.cover_factor) { Lpv = th + off; off += S * K; }
    std::vector<double> Wbuf((std::size_t) N * K);
    for (int t = 0; t < N * K; ++t) Wbuf[t] = th[off + t];
    const int w_off = off; off += N * K;
    const double log_disp = th[off];

    const int Lw = d.cover_factor ? 2 * S * K : S * K;
    Rcpp::NumericVector grad(P + S * P + Lw + N * K + 1);
    double* g = grad.begin();
    double* g_mu = g;
    double* g_b  = g + P;
    double* g_L  = g + P + S * P;            // vec(g_L) (S x K col-major)
    double* g_Lpos = d.cover_factor ? (g_L + S * K) : nullptr;
    double* g_W  = g + w_off;                // reuse packed W offset for field score
    double& g_ld = g[P + S * P + Lw + N * K];

    std::vector<double> th_occ(d.P_occ), th_p(d.P_p), th_pos(d.P_pos);
    std::vector<double> L_s(K), Lpos_s(K);
    std::vector<double> g_occ_s(d.P_occ), g_p_s(d.P_p), g_pos_s(d.P_pos);
    std::vector<double> g_L_s(K), g_Lpos_s(K);
    for (int s = 0; s < S; ++s) {
        const double* b_s = b + s * P;
        for (int j = 0; j < d.P_occ; ++j) th_occ[j] = mu[j] + b_s[j];
        for (int j = 0; j < d.P_p; ++j)   th_p[j]   = mu[d.P_occ + j] + b_s[d.P_occ + j];
        for (int j = 0; j < d.P_pos; ++j) th_pos[j] = mu[d.P_occ + d.P_p + j] +
                                                      b_s[d.P_occ + d.P_p + j];
        for (int k = 0; k < K; ++k) L_s[k] = Lv[k * S + s];
        if (d.cover_factor) for (int k = 0; k < K; ++k) Lpos_s[k] = Lpv[k * S + s];

        tulpaObs::ms_ocs_species_grad(
            d, s, th_occ.data(), th_p.data(), th_pos.data(),
            L_s.data(), d.cover_factor ? Lpos_s.data() : nullptr, log_disp,
            g_occ_s.data(), g_p_s.data(), g_pos_s.data(),
            g_L_s.data(), d.cover_factor ? g_Lpos_s.data() : nullptr,
            g_W, g_ld, Wbuf.data());

        // mu and b_s share the per-arm coefficient score.
        double* gb = g_b + s * P;
        for (int j = 0; j < d.P_occ; ++j) { g_mu[j] += g_occ_s[j]; gb[j] = g_occ_s[j]; }
        for (int j = 0; j < d.P_p; ++j) {
            g_mu[d.P_occ + j] += g_p_s[j]; gb[d.P_occ + j] = g_p_s[j];
        }
        for (int j = 0; j < d.P_pos; ++j) {
            g_mu[d.P_occ + d.P_p + j] += g_pos_s[j];
            gb[d.P_occ + d.P_p + j] = g_pos_s[j];
        }
        for (int k = 0; k < K; ++k) {
            g_L[k * S + s] = g_L_s[k];
            if (d.cover_factor) g_Lpos[k * S + s] = g_Lpos_s[k];
        }
    }
    return grad;
}

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
    std::vector<double> Wbuf((std::size_t) N * K);
    for (int t = 0; t < N * K; ++t) Wbuf[t] = th[off + t];
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
            L_s.data(), d.cover_factor ? Lpos_s.data() : nullptr, log_disp,
            Wbuf.data());
    }
    return total;
}


// Full-vector joint log-posterior + gradient for the spatial-factor community
// occu_cover NUTS target -- the C++ mirror of .ms_ocs_joint_logpost (the FullGradFn
// core). `theta` packs c(par_inner, chol_occ, chol_p, chol_pos, log_tau). icar
// fields, unconstrained loadings (the proper-CAR / BYM2 logit_h block and the
// constrained loadings are later increments). `spec` additionally carries `Q`
// (the ICAR structure, n_sites x n_sites) and `field_rank` (N - 1); `pri` carries
// the hyperprior scalars. Returns list(lp, grad).
namespace tulpaObs {

struct PriScalars {
    double logdiag_mean, logdiag_sd, offdiag_sd,
           log_tau_mean, log_tau_sd, log_disp_mean, log_disp_sd,
           logit_h_mean, logit_h_sd;
};

// Free triangular-loading length: K*S - K(K-1)/2 (matches .ms_ocs_lfree_dim).
inline int ms_ocs_lfree_dim(int S, int K) { return K * S - K * (K - 1) / 2; }

// Free triangular vector -> full S x K loading matrix (col-major), exp() on the
// diagonal, zeros above it. Mirrors .ms_ocs_lfree_to_L.
inline void lfree_to_L_cpp(const double* lf, int S, int K, double* L) {
    for (int t = 0; t < S * K; ++t) L[t] = 0.0;
    int pos = 0;
    for (int k = 0; k < K; ++k) {
        L[(std::size_t) k * S + k] = std::exp(lf[pos++]);
        for (int i = k + 1; i < S; ++i) L[(std::size_t) k * S + i] = lf[pos++];
    }
}

// Gradient wrt full L (S x K col-major) -> gradient wrt the free triangular
// vector: drop the structural zeros, log-link chain on the diagonal. Mirrors
// .ms_ocs_gL_to_glfree.
inline void gL_to_glfree_cpp(const double* gL, const double* L, int S, int K,
                             double* glf) {
    int pos = 0;
    for (int k = 0; k < K; ++k) {
        glf[pos++] = gL[(std::size_t) k * S + k] * L[(std::size_t) k * S + k];
        for (int i = k + 1; i < S; ++i) glf[pos++] = gL[(std::size_t) k * S + i];
    }
}

// Length of the full NUTS coordinate vector. `constrain` swaps the S*K
// unconstrained loading block for the K*S - K(K-1)/2 triangular free block.
inline int ms_ocs_nuts_total(const MsOcsData& d, bool constrain = false) {
    const int P = d.P_occ + d.P_p + d.P_pos;
    const int Lbase = constrain ? ms_ocs_lfree_dim(d.S, d.K) : d.S * d.K;
    const int Lw = Lbase + (d.cover_factor ? d.S * d.K : 0);
    const int chol = d.P_occ * (d.P_occ + 1) / 2 + d.P_p * (d.P_p + 1) / 2
                   + d.P_pos * (d.P_pos + 1) / 2;
    const int hyper = d.has_hyper ? d.K : 0;       // logit_h block
    return P + d.S * P + Lw + d.n_sites * d.K + 1 + chol + d.K + hyper;
}

// Full-vector joint log-posterior + gradient core (icar, unconstrained). `g` is a
// pre-zeroed length-ms_ocs_nuts_total(d) buffer; returns the log-posterior.
// Shared by the Rcpp cross-check entry and the NUTS FullGradFn.
inline double ms_ocs_joint_eval(const MsOcsData& d, const Rcpp::NumericMatrix& Q,
                                double field_rank, const PriScalars& pr,
                                double sigma_beta, double sd_L,
                                const double* th, double* g) {
    using std::vector;
    const int P = d.P_occ + d.P_p + d.P_pos;
    const int N = d.n_sites, S = d.S, K = d.K;
    const double inv_sb2 = 1.0 / (sigma_beta * sigma_beta);
    const double inv_sdL2 = 1.0 / (sd_L * sd_L);

    const double logdiag_mean = pr.logdiag_mean, logdiag_sd = pr.logdiag_sd;
    const double offdiag_sd   = pr.offdiag_sd;
    const double log_tau_mean = pr.log_tau_mean, log_tau_sd = pr.log_tau_sd;
    const double log_disp_mean = pr.log_disp_mean, log_disp_sd = pr.log_disp_sd;

    int off = 0;
    const double* mu = th + off; off += P;
    const double* b  = th + off; off += S * P;
    const double* Lv = th + off; off += S * K;
    const double* Lpv = nullptr;
    if (d.cover_factor) { Lpv = th + off; off += S * K; }
    std::vector<double> Wbuf((std::size_t) N * K);
    for (int t = 0; t < N * K; ++t) Wbuf[t] = th[off + t];
    const int w_off = off; off += N * K;
    const double log_disp = th[off]; const int ld_idx = off; off += 1;
    // hyperparameter blocks
    const int P_arm[3] = {d.P_occ, d.P_p, d.P_pos};
    int q_off[3]; int q_dim[3];
    for (int a = 0; a < 3; ++a) {
        q_dim[a] = P_arm[a] * (P_arm[a] + 1) / 2;
        q_off[a] = off; off += q_dim[a];
    }
    const int tau_off = off; off += K;
    int logit_h_off = -1;
    if (d.has_hyper) { logit_h_off = off; off += K; }
    (void) off;

    double* g_W = g + w_off;
    double& g_ld = g[ld_idx];

    // ---- data log-lik + inner gradient ----
    vector<double> th_occ(d.P_occ), th_p(d.P_p), th_pos(d.P_pos);
    vector<double> L_s(K), Lpos_s(K);
    vector<double> g_occ_s(d.P_occ), g_p_s(d.P_p), g_pos_s(d.P_pos);
    vector<double> g_L_s(K), g_Lpos_s(K);
    double lp = 0.0;
    const int arm_start[3] = {0, d.P_occ, d.P_occ + d.P_p};
    double* g_mu = g; double* g_b = g + P;
    double* g_L = g + P + S * P;
    double* g_Lpos = d.cover_factor ? (g_L + S * K) : nullptr;
    for (int s = 0; s < S; ++s) {
        const double* b_s = b + s * P;
        for (int j = 0; j < d.P_occ; ++j) th_occ[j] = mu[j] + b_s[j];
        for (int j = 0; j < d.P_p; ++j)   th_p[j]   = mu[d.P_occ + j] + b_s[d.P_occ + j];
        for (int j = 0; j < d.P_pos; ++j) th_pos[j] = mu[arm_start[2] + j] + b_s[arm_start[2] + j];
        for (int k = 0; k < K; ++k) L_s[k] = Lv[k * S + s];
        if (d.cover_factor) for (int k = 0; k < K; ++k) Lpos_s[k] = Lpv[k * S + s];

        lp += tulpaObs::ms_ocs_species_grad(
            d, s, th_occ.data(), th_p.data(), th_pos.data(),
            L_s.data(), d.cover_factor ? Lpos_s.data() : nullptr, log_disp,
            g_occ_s.data(), g_p_s.data(), g_pos_s.data(),
            g_L_s.data(), d.cover_factor ? g_Lpos_s.data() : nullptr, g_W, g_ld,
            Wbuf.data());

        double* gb = g_b + s * P;
        for (int j = 0; j < d.P_occ; ++j) { g_mu[j] += g_occ_s[j]; gb[j] = g_occ_s[j]; }
        for (int j = 0; j < d.P_p; ++j) {
            g_mu[d.P_occ + j] += g_p_s[j]; gb[d.P_occ + j] = g_p_s[j];
        }
        for (int j = 0; j < d.P_pos; ++j) {
            g_mu[arm_start[2] + j] += g_pos_s[j]; gb[arm_start[2] + j] = g_pos_s[j];
        }
        for (int k = 0; k < K; ++k) {
            g_L[k * S + s] = g_L_s[k];
            if (d.cover_factor) g_Lpos[k * S + s] = g_Lpos_s[k];
        }
    }

    // ---- community covariance: per-arm b-quadratic + log-det normaliser + chol
    //      block gradient; accumulate the b-prior into g_b. ----
    for (int a = 0; a < 3; ++a) {
        const int Pa = P_arm[a]; if (Pa == 0) continue;
        vector<double> C, Cinv, Si;
        tulpaObs::chol_unpack_cpp(th + q_off[a], Pa, C);
        tulpaObs::lower_tri_inv(C, Pa, Cinv);
        tulpaObs::sinv_from_cinv(Cinv, Pa, Si);
        double logdet = 0.0;
        for (int j = 0; j < Pa; ++j) logdet += 2.0 * std::log(C[(std::size_t) j * Pa + j]);

        // M_arm = sum_s b_{s,arm} b_{s,arm}', and the b-quadratic + its gradient.
        vector<double> M((std::size_t) Pa * Pa, 0.0);
        double quad_sum = 0.0;
        for (int s = 0; s < S; ++s) {
            const double* bsa = b + s * P + arm_start[a];
            // Si b_{s,arm}
            for (int u = 0; u < Pa; ++u) {
                double sib = 0.0;
                for (int v = 0; v < Pa; ++v) sib += Si[(std::size_t) u * Pa + v] * bsa[v];
                g_b[s * P + arm_start[a] + u] -= sib;        // b-prior gradient
                quad_sum += bsa[u] * sib;
                for (int v = 0; v < Pa; ++v) M[(std::size_t) u * Pa + v] += bsa[u] * bsa[v];
            }
        }
        lp += -0.5 * quad_sum - 0.5 * S * logdet;

        tulpaObs::chol_block_grad_cpp(C, Si, M, Pa, (double) S, th + q_off[a],
                                      logdiag_mean, logdiag_sd, offdiag_sd,
                                      g + q_off[a]);
        // chol coordinate hyperprior contribution to lp
        int pos = q_off[a];
        for (int j = 0; j < Pa; ++j) {
            const double vd = th[pos++];
            lp += -0.5 * ((vd - logdiag_mean) / logdiag_sd) * ((vd - logdiag_mean) / logdiag_sd);
            for (int i = j + 1; i < Pa; ++i) {
                const double vo = th[pos++];
                lp += -0.5 * (vo / offdiag_sd) * (vo / offdiag_sd);
            }
        }
    }

    // ---- mu / loading Gaussian priors (gradients + lp) ----
    for (int j = 0; j < P; ++j) { g_mu[j] -= inv_sb2 * mu[j]; lp += -0.5 * inv_sb2 * mu[j] * mu[j]; }
    {
        double sumL2 = 0.0;
        for (int t = 0; t < S * K; ++t) { g_L[t] -= inv_sdL2 * Lv[t]; sumL2 += Lv[t] * Lv[t]; }
        lp += -0.5 * inv_sdL2 * sumL2;
        if (d.cover_factor) {
            double sumLp2 = 0.0;
            for (int t = 0; t < S * K; ++t) { g_Lpos[t] -= inv_sdL2 * Lpv[t]; sumLp2 += Lpv[t] * Lpv[t]; }
            lp += -0.5 * inv_sdL2 * sumLp2;
        }
    }

    // ---- field GMRF: per-factor tau quadratic + rank normaliser + (proper
    //      fields) the h-dependent log|R(h)| and the logit_h gradient. The field
    //      structure R(h) is icar Q, proper-CAR D - h A, or bym2 V diag(1/m) V'.
    (void) Q;
    const double logit_h_mean = pr.logit_h_mean, logit_h_sd = pr.logit_h_sd;
    std::vector<double> RWk(N);
    for (int k = 0; k < K; ++k) {
        const double log_tau = th[tau_off + k];
        const double tau = std::exp(log_tau);
        const double h = d.has_hyper
                       ? 1.0 / (1.0 + std::exp(-th[logit_h_off + k])) : 0.0;
        field_RWk(d, h, Wbuf.data() + (std::size_t) k * N, RWk.data());
        double quad = 0.0;
        for (int i = 0; i < N; ++i) {
            quad += Wbuf[(std::size_t) k * N + i] * RWk[i];
            g_W[(std::size_t) k * N + i] -= tau * RWk[i];    // field-prior gradient
        }
        lp += -0.5 * tau * quad + 0.5 * field_rank * log_tau;
        lp += -0.5 * ((log_tau - log_tau_mean) / log_tau_sd) * ((log_tau - log_tau_mean) / log_tau_sd);
        g[tau_off + k] = -0.5 * tau * quad + 0.5 * field_rank
                       - (log_tau - log_tau_mean) / log_tau_sd / log_tau_sd;
        if (d.has_hyper) {
            lp += 0.5 * field_logdetR(d, h);
            double dquad = 0.0, dlogdet = 0.0;
            field_hyper_terms(d, h, Wbuf.data() + (std::size_t) k * N, dquad, dlogdet);
            const double dlp_dh = -0.5 * tau * dquad + 0.5 * dlogdet;
            const double dh = h * (1.0 - h);                 // plogis derivative
            const double zc = th[logit_h_off + k];
            lp += -0.5 * ((zc - logit_h_mean) / logit_h_sd) * ((zc - logit_h_mean) / logit_h_sd);
            g[logit_h_off + k] = dlp_dh * dh
                               - (zc - logit_h_mean) / logit_h_sd / logit_h_sd;
        }
    }

    // ---- log_disp prior ----
    lp += -0.5 * ((log_disp - log_disp_mean) / log_disp_sd) * ((log_disp - log_disp_mean) / log_disp_sd);
    g_ld -= (log_disp - log_disp_mean) / (log_disp_sd * log_disp_sd);

    return lp;
}

// Constrained-loadings dispatch: when `constrain`, the L block is the triangular
// free vector. Expand it to a full loading matrix, defer to the (validated)
// unconstrained core, then map the L-block gradient back to the free vector by
// the chain rule -- the reparameterisation adapter, mirror of .ms_ocs_penll_grad_c.
inline double ms_ocs_joint_eval_c(const MsOcsData& d, const Rcpp::NumericMatrix& Q,
                                  double field_rank, const PriScalars& pr,
                                  double sigma_beta, double sd_L, bool constrain,
                                  const double* th, double* g) {
    if (!constrain)
        return ms_ocs_joint_eval(d, Q, field_rank, pr, sigma_beta, sd_L, th, g);
    const int P = d.P_occ + d.P_p + d.P_pos;
    const int S = d.S, K = d.K, N = d.n_sites;
    const int head = P + S * P;
    const int nL = ms_ocs_lfree_dim(S, K);
    const int SK = S * K;
    const int chol = d.P_occ*(d.P_occ+1)/2 + d.P_p*(d.P_p+1)/2 + d.P_pos*(d.P_pos+1)/2;
    const int tail = (d.cover_factor ? SK : 0) + N * K + 1 + chol + K
                   + (d.has_hyper ? K : 0);

    std::vector<double> Lfull((std::size_t) SK);
    lfree_to_L_cpp(th + head, S, K, Lfull.data());
    std::vector<double> th_full((std::size_t) head + SK + tail);
    std::copy(th, th + head, th_full.begin());
    std::copy(Lfull.begin(), Lfull.end(), th_full.begin() + head);
    std::copy(th + head + nL, th + head + nL + tail, th_full.begin() + head + SK);

    std::vector<double> g_full((std::size_t) head + SK + tail, 0.0);
    const double lp = ms_ocs_joint_eval(d, Q, field_rank, pr, sigma_beta, sd_L,
                                        th_full.data(), g_full.data());
    std::copy(g_full.begin(), g_full.begin() + head, g);
    gL_to_glfree_cpp(g_full.data() + head, Lfull.data(), S, K, g + head);
    std::copy(g_full.begin() + head + SK, g_full.end(), g + head + nL);
    return lp;
}

// NUTS model carrying the marshalled data + hyperparameters; the FullGradFn
// reaches it through ModelData.model_response_data.
struct MsOcsNutsModel {
    MsOcsData d;
    Rcpp::NumericMatrix Q;
    double field_rank = 0.0, sigma_beta = 5.0, sd_L = 1.0;
    PriScalars pr;
    bool constrain = false;
    int total = 0;
};

// FullGradFn: log-posterior + gradient over the entire parameter vector. NUTS
// maximises the value, so this returns the log-posterior (no negation).
inline void ms_ocs_full_grad(const std::vector<double>& params,
                             const tulpa::ModelData& data,
                             const tulpa::ParamLayout& /*layout*/,
                             std::vector<double>& grad, double* log_post_out) {
    const MsOcsNutsModel* m =
        static_cast<const MsOcsNutsModel*>(data.model_response_data);
    grad.assign((std::size_t) m->total, 0.0);
    const double lp = ms_ocs_joint_eval_c(m->d, m->Q, m->field_rank, m->pr,
                                          m->sigma_beta, m->sd_L, m->constrain,
                                          params.data(), grad.data());
    if (log_post_out) *log_post_out = lp;
}

inline PriScalars ms_ocs_pri_from_list(const Rcpp::List& pri) {
    return PriScalars{ Rcpp::as<double>(pri["chol_logdiag_mean"]),
                       Rcpp::as<double>(pri["chol_logdiag_sd"]),
                       Rcpp::as<double>(pri["chol_offdiag_sd"]),
                       Rcpp::as<double>(pri["log_tau_mean"]),
                       Rcpp::as<double>(pri["log_tau_sd"]),
                       Rcpp::as<double>(pri["log_disp_mean"]),
                       Rcpp::as<double>(pri["log_disp_sd"]),
                       Rcpp::as<double>(pri["logit_h_mean"]),
                       Rcpp::as<double>(pri["logit_h_sd"]) };
}

} // namespace tulpaObs

// Full-vector joint log-posterior + gradient, the Rcpp cross-check entry for the
// R oracle .ms_ocs_joint_logpost. icar, unconstrained loadings.
// [[Rcpp::export]]
Rcpp::List cpp_ms_ocs_joint_logpost(Rcpp::List spec, Rcpp::NumericVector theta,
                                    Rcpp::List pri, double sigma_beta,
                                    double sd_L, bool constrain = false) {
    tulpaObs::MsOcsData d = tulpaObs::ms_ocs_build_data(spec);
    Rcpp::NumericMatrix Q = spec["Q"];
    const double field_rank = Rcpp::as<double>(spec["field_rank"]);
    tulpaObs::PriScalars pr = tulpaObs::ms_ocs_pri_from_list(pri);
    const int total = tulpaObs::ms_ocs_nuts_total(d, constrain);
    Rcpp::NumericVector grad(total);
    const double lp = tulpaObs::ms_ocs_joint_eval_c(
        d, Q, field_rank, pr, sigma_beta, sd_L, constrain,
        theta.begin(), grad.begin());
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// Run NUTS on the spatial-factor community occu_cover target via tulpa's engine
// and the FullGradFn (gradient mode "H"). icar, unconstrained loadings. `theta0`
// is the warm-start (the Laplace mode); `inv_metric` an optional length-n_params
// inverse-mass diagonal (the Laplace curvature). Returns draws + diagnostics.
// [[Rcpp::export]]
Rcpp::List cpp_ms_ocs_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                           Rcpp::List pri, double sigma_beta, double sd_L,
                           Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                           int n_iter, int n_warmup, int max_treedepth,
                           double adapt_delta, int seed, bool verbose,
                           bool constrain = false) {
    tulpaObs::MsOcsNutsModel m;
    m.d = tulpaObs::ms_ocs_build_data(spec);
    m.Q = Rcpp::as<Rcpp::NumericMatrix>(spec["Q"]);
    m.field_rank = Rcpp::as<double>(spec["field_rank"]);
    m.sigma_beta = sigma_beta; m.sd_L = sd_L;
    m.pr = tulpaObs::ms_ocs_pri_from_list(pri);
    m.constrain = constrain;
    m.total = tulpaObs::ms_ocs_nuts_total(m.d, constrain);
    if ((int) theta0.size() != m.total)
        Rcpp::stop("theta0 length %d != expected %d", (int) theta0.size(), m.total);

    tulpa::LikelihoodSpec lspec;
    lspec.name = "ms_occu_cover_spatial";
    lspec.n_processes = 1;
    lspec.gradient_fn = &tulpaObs::ms_ocs_full_grad;

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
    layout.total_params = m.total;

    tulpa::set_gradient_mode_str("H");

    std::vector<double> init(theta0.begin(), theta0.end());
    std::vector<double> imv;
    const double* im = nullptr;
    if (inv_metric.isNotNull()) {
        Rcpp::NumericVector v(inv_metric);
        imv.assign(v.begin(), v.end()); im = imv.data();
    }

    tulpa::NUTSFn run_nuts = tulpa::get_nuts_fn();
    tulpa::NUTSResult result = {};
    run_nuts(&data, &layout, init.data(), m.total, n_iter, n_warmup,
             max_treedepth, adapt_delta, static_cast<unsigned int>(seed),
             verbose ? 1 : 0, im, &result);

    const int n_samples = result.n_sample, np = m.total;
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

