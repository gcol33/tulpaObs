// nuts_field_hyper.h
// Non-centered areal field block for every NUTS target that carries one. The
// field marginal SD (sigma), the mixing / spatial-correlation parameter (rho)
// and the cross-arm copy amplitude (alpha) may each be a coordinate of the
// sampled vector, so the sampler integrates the outer hyperparameter layer
// itself instead of conditioning on the deterministic backend's summary of it.
//
// Pinning every hyper is the degenerate configuration of this same block, not a
// second implementation: sigma = 1 with a constant column scaling and no iid
// block reduces the forward to z = B1 raw and the backward to
// grad += B1' grad_z - raw, lp -= 0.5 ||raw||^2, term for term and in the same
// summation order. That is what the count / observation families use, with the
// loading carrying their nested-Laplace tau Q(rho) already baked into its
// columns.
//
// The field is
//
//   z = sigma * ( B1 %*% (s1(rho) * raw1) + s2(rho) * raw2 ),
//   raw1 ~ N(0, I_m1), raw2 ~ N(0, I_n)          (raw2 present for bym2 only)
//
// with B1 a FIXED n x m1 basis and s1 / s2 scalar-or-per-column scalings. Every
// areal kind this package fixes hypers for factors that way, so no leapfrog step
// re-decomposes anything:
//
//   icar        B1 = U_+ diag(1 / sqrt(lambda_+))  (sum-to-zero eigen-loading of
//               the intrinsic precision Q), s1 = 1, no raw2.
//               Cov(z) = sigma^2 Q^+ on the sum-to-zero subspace.
//   bym2        B1 = the same centred ICAR basis, s1 = sqrt(rho / scale_factor),
//               raw2 = the unstructured iid block with s2 = sqrt(1 - rho)
//               (Riebler 2016), so rho re-weights the two blocks in place.
//   car_proper  Q(rho) = D - rho W = D^{1/2} (I - rho Lambda) D^{1/2} in the
//               eigenbasis of the symmetrically normalised adjacency
//               D^{-1/2} W D^{-1/2} = U Lambda U'. Hence
//               (tau Q(rho))^{-1} = sigma^2 B1 diag(1 / (1 - rho lambda_j)) B1'
//               with B1 = D^{-1/2} U fixed, i.e. s1_j = (1 - rho lambda_j)^{-1/2}.
//               rho rescales eigenvalues in a fixed basis -- no per-step
//               Cholesky, contrary to the O(n^3) reading of Q(rho) = D - rho W.
//
// Each sampled hyper rides an unconstrained coordinate u through a doubly
// bounded transform, so the sampler sees no wall:
//
//   t = t_lo + (t_hi - t_lo) * expit(u),   value = inv_link(t)
//
// with `t` the axis's own coordinate (log for sigma / alpha, logit for rho) and
// [t_lo, t_hi] its bounds. The prior is FLAT in t over those bounds -- the
// measure the nested-Laplace outer grid integrates against, which weights its
// log-spaced sigma / alpha and logit-spaced rho nodes equally. Flat prior plus
// change-of-variables leaves a normalised log-density of log(e) + log(1 - e),
// e = expit(u), per sampled hyper.
//
// A block may carry a per-site design WEIGHT, which is what makes it a
// spatially-varying coefficient rather than a second intercept field: site i
// loads w_i * z[unit(i)] on the state arm and alpha * w_i * z[unit(i)] on the
// copied arm. An absent weight vector is the unit weight, so an intercept field
// carries none. Several blocks stack -- each with its own basis, site -> unit
// map, weight and (sigma, rho, alpha) coordinates -- and the flat vector lays
// them out back to back, so a fit that carries one block has the layout it had
// before the second was possible.

#ifndef TULPAOBS_NUTS_FIELD_HYPER_H
#define TULPAOBS_NUTS_FIELD_HYPER_H

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <cstddef>

namespace tulpaObs {

// Per-column scaling of a basis block as a function of rho.
enum HyperFieldScale {
    HF_SCALE_CONST    = 0,   // s = 1                       (icar / pinned loading)
    HF_SCALE_BYM2_STR = 1,   // s = sqrt(rho / scale_factor)
    HF_SCALE_BYM2_IID = 2,   // s = sqrt(1 - rho)
    HF_SCALE_CAR      = 3    // s_j = (1 - rho lambda_j)^{-1/2}
};

// One bounded, sampled (or pinned) hyperparameter. `link` is 0 for log (sigma,
// alpha) and 1 for logit (rho): the coordinate the outer grid spaces its nodes
// in, hence the coordinate the flat prior is flat in. `coord < 0` marks a pinned
// hyper, whose value is `fixed` and which contributes nothing to the gradient.
struct HyperCoord {
    int coord = -1;          // index in theta, or < 0 when pinned
    int link = 0;            // 0 = log, 1 = logit
    double fixed = 0.0;      // value when pinned
    double t_lo = 0.0, t_hi = 0.0;

    bool sampled() const { return coord >= 0; }
};

// Value of a hyper at `theta`, plus the pieces the chain rule needs:
//   value       the natural-scale hyper
//   dvalue_dt   d value / d t             (value for log, value(1-value) for logit)
//   dt_du       d t / d u                 ((t_hi - t_lo) e (1 - e))
//   e           expit(u), for the Jacobian gradient (1 - 2e)
struct HyperValue {
    double value = 0.0, dvalue_dt = 0.0, dt_du = 0.0, e = 0.0;
};

inline HyperValue hyper_coord_value(const HyperCoord& h, const double* theta) {
    HyperValue hv;
    if (!h.sampled()) { hv.value = h.fixed; return hv; }
    const double u = theta[h.coord];
    const double e = 1.0 / (1.0 + std::exp(-u));
    const double t = h.t_lo + (h.t_hi - h.t_lo) * e;
    if (h.link == 1) {                       // logit
        const double v = 1.0 / (1.0 + std::exp(-t));
        hv.value = v; hv.dvalue_dt = v * (1.0 - v);
    } else {                                 // log
        const double v = std::exp(t);
        hv.value = v; hv.dvalue_dt = v;
    }
    hv.dt_du = (h.t_hi - h.t_lo) * e * (1.0 - e);
    hv.e = e;
    return hv;
}

// Add one sampled hyper's contribution to `grad` from d log p / d value, and
// return its (normalised) flat-prior-plus-Jacobian log-density. Pinned: no-op.
inline double hyper_coord_backward(const HyperCoord& h, const HyperValue& hv,
                                   double dlp_dvalue, double* grad) {
    if (!h.sampled()) return 0.0;
    grad[h.coord] += dlp_dvalue * hv.dvalue_dt * hv.dt_du + (1.0 - 2.0 * hv.e);
    const double e = hv.e;
    return std::log(e) + std::log(1.0 - e);
}

// Marshalled field block. n_units == 0 => inactive (every helper is a no-op and
// the flat vector carries no field coordinates).
struct HyperFieldBlock {
    int n_units = 0;                  // field nodes (rows of B1)
    int m1 = 0;                       // columns of B1
    int scale1 = HF_SCALE_CONST;      // rho-dependence of the B1 columns
    bool has_iid = false;             // bym2 unstructured block (n_units coords)
    int n_raw = 0;                    // m1 + (has_iid ? n_units : 0)
    int o_raw = 0;                    // first whitened coordinate in theta
    double sf = 1.0;                  // bym2 scale factor
    std::vector<double> B1;           // row-major n_units x m1
    std::vector<double> lambda;       // car_proper: normalised-adjacency eigenvalues
    std::vector<int> field_map;       // 0-based unit per site (length n_sites)
    std::vector<double> weight;       // per-site design weight (empty = unit)
    HyperCoord sigma, rho, alpha;

    bool active() const { return n_units > 0; }
    double site_weight(int i) const { return weight.empty() ? 1.0 : weight[i]; }
};

// Working state of one forward pass: the hyper values and the field itself.
struct HyperFieldState {
    HyperValue sigma, rho, alpha;
    std::vector<double> s1;           // per-column scaling of B1 (length m1)
    std::vector<double> ds1;          // d s1 / d rho
    double s2 = 0.0, ds2 = 0.0;       // iid block scaling and its rho-derivative
    std::vector<double> z;            // field value per unit
};

// Total flat coordinates the block contributes (whitened field + sampled hypers).
inline int hyper_field_size(const HyperFieldBlock& fb) {
    if (!fb.active()) return 0;
    int n = fb.n_raw;
    if (fb.sigma.sampled()) ++n;
    if (fb.rho.sampled())   ++n;
    if (fb.alpha.sampled()) ++n;
    return n;
}

// Read the block from a NUTS spec. `base` is the first free flat coordinate; the
// whitened field takes n_raw of them and each sampled hyper one more, in the
// order (sigma, rho, alpha). A spec carrying only the legacy `field_load` /
// `field_Linv` + `field_alpha` entries marshals as the pinned block -- sigma
// pinned at 1 with a constant scaling, which reproduces the fixed-hyper loading
// exactly (that loading already has sigma and rho baked into its columns).
inline HyperFieldBlock hyper_field_build(const Rcpp::List& spec, int base,
                                         int n_sites) {
    HyperFieldBlock fb;
    if (!spec.containsElementNamed("n_field_units")) return fb;
    fb.n_units = Rcpp::as<int>(spec["n_field_units"]);
    if (fb.n_units <= 0) { fb.n_units = 0; return fb; }

    Rcpp::IntegerVector fm = spec["field_map"];        // 1-based site -> unit
    if ((int) fm.size() != n_sites)
        Rcpp::stop("field_map must have length n_sites");
    fb.field_map.resize(n_sites);
    for (int i = 0; i < n_sites; ++i) {
        const int u = fm[i] - 1;
        if (u < 0 || u >= fb.n_units)
            Rcpp::stop("field_map values must lie in [1, n_field_units]");
        fb.field_map[i] = u;
    }

    Rcpp::NumericMatrix B = spec.containsElementNamed("field_load")
        ? Rcpp::as<Rcpp::NumericMatrix>(spec["field_load"])
        : Rcpp::as<Rcpp::NumericMatrix>(spec["field_Linv"]);
    if (B.nrow() != fb.n_units)
        Rcpp::stop("field loading must have n_field_units rows");
    fb.m1 = B.ncol();
    fb.B1.resize((std::size_t) fb.n_units * fb.m1);
    for (int u = 0; u < fb.n_units; ++u)
        for (int j = 0; j < fb.m1; ++j)
            fb.B1[(std::size_t) u * fb.m1 + j] = B(u, j);

    if (spec.containsElementNamed("field_weight")) {
        Rcpp::NumericVector w = spec["field_weight"];
        if ((int) w.size() != n_sites)
            Rcpp::stop("field_weight must have length n_sites");
        fb.weight.assign(w.begin(), w.end());
    }

    if (spec.containsElementNamed("field_scale1"))
        fb.scale1 = Rcpp::as<int>(spec["field_scale1"]);
    if (spec.containsElementNamed("field_has_iid"))
        fb.has_iid = Rcpp::as<int>(spec["field_has_iid"]) != 0;
    if (spec.containsElementNamed("field_sf"))
        fb.sf = Rcpp::as<double>(spec["field_sf"]);
    if (spec.containsElementNamed("field_lambda")) {
        Rcpp::NumericVector lam = spec["field_lambda"];
        fb.lambda.assign(lam.begin(), lam.end());
    }
    if (fb.scale1 == HF_SCALE_CAR && (int) fb.lambda.size() != fb.m1)
        Rcpp::stop("field_lambda must have one eigenvalue per basis column");
    fb.n_raw = fb.m1 + (fb.has_iid ? fb.n_units : 0);
    fb.o_raw = base;
    int k = base + fb.n_raw;

    // Hyper coordinates. A `<name>_lo` / `<name>_hi` pair marks the hyper as
    // sampled between those bounds in its own link coordinate; otherwise it is
    // pinned at `<name>_fixed`.
    struct Slot { const char* fixed; const char* lo; const char* hi; int link;
                  HyperCoord* h; double dflt; };
    const Slot slots[3] = {
        {"field_sigma_fixed", "field_sigma_lo", "field_sigma_hi", 0, &fb.sigma, 1.0},
        {"field_rho_fixed",   "field_rho_lo",   "field_rho_hi",   1, &fb.rho,   1.0},
        {"field_alpha_fixed", "field_alpha_lo", "field_alpha_hi", 0, &fb.alpha, 0.0}
    };
    for (int s = 0; s < 3; ++s) {
        const Slot& sl = slots[s];
        sl.h->link = sl.link;
        sl.h->fixed = spec.containsElementNamed(sl.fixed)
            ? Rcpp::as<double>(spec[sl.fixed]) : sl.dflt;
        if (spec.containsElementNamed(sl.lo) && spec.containsElementNamed(sl.hi)) {
            sl.h->t_lo = Rcpp::as<double>(spec[sl.lo]);
            sl.h->t_hi = Rcpp::as<double>(spec[sl.hi]);
            sl.h->coord = k++;
        }
    }
    // The legacy pinned spec carries the copy amplitude under `field_alpha`.
    if (!fb.alpha.sampled() && !spec.containsElementNamed("field_alpha_fixed") &&
        spec.containsElementNamed("field_alpha"))
        fb.alpha.fixed = Rcpp::as<double>(spec["field_alpha"]);
    return fb;
}

// Forward pass: hyper values, block scalings and the field z. No-op (clears the
// state) when inactive.
inline void hyper_field_forward(const HyperFieldBlock& fb, const double* theta,
                                HyperFieldState& st) {
    if (!fb.active()) { st.z.clear(); st.alpha.value = 0.0; return; }
    st.sigma = hyper_coord_value(fb.sigma, theta);
    st.rho   = hyper_coord_value(fb.rho,   theta);
    st.alpha = hyper_coord_value(fb.alpha, theta);

    const double rho = st.rho.value;
    st.s1.assign(fb.m1, 1.0);
    st.ds1.assign(fb.m1, 0.0);
    switch (fb.scale1) {
        case HF_SCALE_BYM2_STR: {
            const double s = std::sqrt(rho / fb.sf);
            const double d = 1.0 / (2.0 * std::sqrt(rho * fb.sf));
            for (int j = 0; j < fb.m1; ++j) { st.s1[j] = s; st.ds1[j] = d; }
            break;
        }
        case HF_SCALE_CAR: {
            for (int j = 0; j < fb.m1; ++j) {
                const double a = 1.0 - rho * fb.lambda[j];
                const double r = 1.0 / std::sqrt(a);
                st.s1[j] = r;
                st.ds1[j] = 0.5 * fb.lambda[j] * r / a;
            }
            break;
        }
        default: break;                       // HF_SCALE_CONST
    }
    if (fb.has_iid) {
        st.s2  = std::sqrt(1.0 - rho);
        st.ds2 = -1.0 / (2.0 * std::sqrt(1.0 - rho));
    } else { st.s2 = 0.0; st.ds2 = 0.0; }

    const double sigma = st.sigma.value;
    const double* raw1 = theta + fb.o_raw;
    const double* raw2 = theta + fb.o_raw + fb.m1;
    st.z.assign(fb.n_units, 0.0);
    for (int u = 0; u < fb.n_units; ++u) {
        const double* Bu = &fb.B1[(std::size_t) u * fb.m1];
        double zz = 0.0;
        for (int j = 0; j < fb.m1; ++j) zz += Bu[j] * st.s1[j] * raw1[j];
        if (fb.has_iid) zz += st.s2 * raw2[u];
        st.z[u] = sigma * zz;
    }
}

// Backward pass. `grad_z[u]` is d log L / d z[u] summed over the sites mapped to
// unit u (the caller folds the copy amplitude into it: psi-arm score plus alpha
// times the copied arm's score). `g_alpha_data` is d log L / d alpha from the
// data term (sum over sites of z[unit] times the copied arm's eta-score), which
// only the caller can form. Adds every field gradient into `grad` and returns
// the whitened-field prior plus the sampled hypers' log-densities.
inline double hyper_field_backward(const HyperFieldBlock& fb, const double* theta,
                                   const HyperFieldState& st,
                                   const std::vector<double>& grad_z,
                                   double g_alpha_data, double* grad) {
    if (!fb.active()) return 0.0;
    const double sigma = st.sigma.value;
    const double* raw1 = theta + fb.o_raw;
    const double* raw2 = theta + fb.o_raw + fb.m1;
    double lp = 0.0;

    // B1' grad_z, shared by the raw1 gradient and the rho gradient.
    std::vector<double> BtG(fb.m1, 0.0);
    for (int u = 0; u < fb.n_units; ++u) {
        const double gz = grad_z[u];
        if (gz == 0.0) continue;
        const double* Bu = &fb.B1[(std::size_t) u * fb.m1];
        for (int j = 0; j < fb.m1; ++j) BtG[j] += Bu[j] * gz;
    }

    double g_sigma_t = 0.0;      // d log p / d t_sigma = sum_u grad_z[u] z[u]
    double g_rho = 0.0;          // d log p / d rho
    for (int u = 0; u < fb.n_units; ++u) g_sigma_t += grad_z[u] * st.z[u];

    for (int j = 0; j < fb.m1; ++j) {
        const double rv = raw1[j];
        grad[fb.o_raw + j] += sigma * st.s1[j] * BtG[j] - rv;
        lp -= 0.5 * rv * rv;
        g_rho += sigma * st.ds1[j] * rv * BtG[j];
    }
    if (fb.has_iid) {
        for (int u = 0; u < fb.n_units; ++u) {
            const double rv = raw2[u];
            grad[fb.o_raw + fb.m1 + u] += sigma * st.s2 * grad_z[u] - rv;
            lp -= 0.5 * rv * rv;
            g_rho += sigma * st.ds2 * rv * grad_z[u];
        }
    }

    // d log p / d sigma = g_sigma_t / sigma, and d sigma / d t_sigma = sigma, so
    // the log-link chain hands the hyper helper g_sigma_t / sigma directly.
    lp += hyper_coord_backward(fb.sigma, st.sigma,
                               sigma > 0.0 ? g_sigma_t / sigma : 0.0, grad);
    lp += hyper_coord_backward(fb.rho,   st.rho,   g_rho,        grad);
    lp += hyper_coord_backward(fb.alpha, st.alpha, g_alpha_data, grad);
    return lp;
}

// --------------------------------------------------------------------------
// Several field blocks on one arm pair
//
// Block b contributes w_b(i) * z_b[unit_b(i)] to the state arm and
// alpha_b * w_b(i) * z_b[unit_b(i)] to the copied arm. The three helpers below
// are the ONLY places that loading is written, so the target's eta assembly and
// its score cannot express it differently, and a one-block fit runs the same
// arithmetic it ran when the block count was fixed at one.
// --------------------------------------------------------------------------

// Site i's loading of block b on the state arm.
inline double hyper_field_site_value(const HyperFieldBlock& fb,
                                     const HyperFieldState& st, int i) {
    if (!fb.active()) return 0.0;
    return fb.site_weight(i) * st.z[fb.field_map[i]];
}

// Site i's total offset on the state arm and on the copied arm.
inline void hyper_field_site_offsets(const std::vector<HyperFieldBlock>& fbs,
                                     const std::vector<HyperFieldState>& sts,
                                     int i, double& state_off, double& copy_off) {
    state_off = 0.0;
    copy_off  = 0.0;
    for (std::size_t b = 0; b < fbs.size(); ++b) {
        const double v = hyper_field_site_value(fbs[b], sts[b], i);
        state_off += v;
        copy_off  += sts[b].alpha.value * v;
    }
}

// Scatter site i's arm scores onto each block's per-unit field score and its
// copy amplitude's data score. `g_state` is d log L / d eta_state at this site;
// `g_copy` is the copied arm's eta-score summed over the site's rows.
inline void hyper_field_site_score(const std::vector<HyperFieldBlock>& fbs,
                                   const std::vector<HyperFieldState>& sts,
                                   int i, double g_state, double g_copy,
                                   std::vector<std::vector<double> >& g_z,
                                   std::vector<double>& g_alpha_data) {
    for (std::size_t b = 0; b < fbs.size(); ++b) {
        if (!fbs[b].active()) continue;
        const double w = fbs[b].site_weight(i);
        g_z[b][fbs[b].field_map[i]] += w * (g_state + sts[b].alpha.value * g_copy);
        g_alpha_data[b] += hyper_field_site_value(fbs[b], sts[b], i) * g_copy;
    }
}

// Read every field block a spec declares. `field_blocks` is the plural spelling
// (one list per block, each carrying the same entries a single block does); a
// spec with the entries at top level is the one-block case.
inline std::vector<HyperFieldBlock> hyper_field_build_list(const Rcpp::List& spec,
                                                           int base, int n_sites) {
    std::vector<HyperFieldBlock> out;
    if (spec.containsElementNamed("field_blocks")) {
        Rcpp::List blocks = spec["field_blocks"];
        for (int b = 0; b < blocks.size(); ++b) {
            HyperFieldBlock fb = hyper_field_build(
                Rcpp::as<Rcpp::List>(blocks[b]), base, n_sites);
            if (!fb.active())
                Rcpp::stop("field_blocks[[%d]] declares no field units", b + 1);
            base += hyper_field_size(fb);
            out.push_back(fb);
        }
        return out;
    }
    HyperFieldBlock fb = hyper_field_build(spec, base, n_sites);
    if (fb.active()) out.push_back(fb);
    return out;
}

inline int hyper_field_list_size(const std::vector<HyperFieldBlock>& fbs) {
    int n = 0;
    for (std::size_t b = 0; b < fbs.size(); ++b) n += hyper_field_size(fbs[b]);
    return n;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_NUTS_FIELD_HYPER_H
