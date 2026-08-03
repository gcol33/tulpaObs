# Latent-N truncation policy for the binned distance-sampling family, the
# per-site headroom cap ported from R/nmix_site_marginal.R (gcol33/tulpaObs#168).
# The tolerances and the widening step are family-agnostic (.NMIX_SCORE_TOL,
# .NMIX_BOUNDARY_TOL, .nmix_widen_headroom(), R/nmix_site_marginal.R) and are
# reused here rather than re-declared, so distance and the N-mixture stay on
# ONE guard policy.
#
# Distance differs from the N-mixture in two ways a shared helper has to
# absorb:
#
#  1. K_lo == R_i (the site's own detected total) EXACTLY for every distance
#     site (see src/distance_kernel.h), so the #167 comb_table -- already
#     indexed by the offset k = N - K_lo, never by K_lo itself -- stays valid
#     unchanged under a per-site cap; only the loop's own upper bound shrinks
#     (compute_distance_site()'s K_hi). No separate per-site cache struct is
#     needed the way nmix_precompute_site() needs one.
#
#  2. The default K_max is MULTIPLICATIVE (3*R_max + 100, not R_max + 100):
#     the truncation bounds the latent TOTAL abundance in the covered region,
#     which exceeds the detected total by a factor of 1/p (undetected
#     individuals), so a "detected total + h" buffer is too tight whenever
#     detection is low. The default headroom is therefore derived so the
#     binding site (R_i = R_max, the one the shared ceiling already sizes
#     itself around) gets EXACTLY the states the old shared-K_max default gave
#     it: K_max = R_max + headroom = 3*R_max + 100 => headroom = 2*R_max + 100.
#     Every other site gets that SAME absolute headroom above its own total --
#     the nmix mechanism, re-derived for distance's multiplicative default
#     rather than nmix's fixed additive constant (.NMIX_HEADROOM).

.dist_truncation <- function(K_max, site_tot) {
  R_max <- if (length(site_tot)) max(site_tot) else 0L
  if (is.null(K_max)) {
    headroom <- as.integer(2L * R_max + 100L)
    return(list(K_max = as.integer(3L * R_max + 100L), headroom = headroom))
  }
  list(K_max = as.integer(K_max), headroom = -1L)
}

# Coefficient-space score gap between a capped-headroom and an uncapped
# (shared-ceiling) evaluation, from their eta-space gradients sandwiched by the
# design matrices -- the same crossprod nmix_laplace()'s inline guard forms,
# factored out here since every distance fitter (Laplace, spatial, NUTS,
# community) needs the identical comparison over its own (X_lambda, X_sigma).
# Largest absolute disagreement is what the guard escalates on, matching
# .NMIX_SCORE_TOL's meaning: an absolute tolerance on the log-likelihood
# gradient in coefficient units.
.dist_score_gap <- function(grad_lambda_h, grad_lambda_u, X_lambda,
                            grad_sigma_h, grad_sigma_u, X_sigma) {
  max(abs(c(
    as.numeric(crossprod(X_lambda, grad_lambda_h - grad_lambda_u)),
    as.numeric(crossprod(X_sigma,  grad_sigma_h  - grad_sigma_u))
  )))
}

# Largest disagreement between the capped and uncapped score across every
# species in a community distance fit, at the community EM's current
# per-species coefficients -- the distance counterpart of
# .nmix_community_score_gap() (R/nmix_site_marginal.R), built directly over
# .tobs_ms_distance_engine()'s exposed `y_s` / `quad_xptr` / `key_code` rather
# than a fresh per-species marginal object, since the shared comb_table (#167)
# makes constructing one unnecessary here. `site_off` (length n_sites) is any
# shared per-site abundance offset (a spatial field); `fac_off`
# (n_sites x n_species) is any per-species latent-factor offset. Both are
# folded into eta_lambda exactly as the community EM's own `eta_of()` does, so
# the guard checks the truncation at the predictor the fit actually ran under.
.dist_community_score_gap <- function(eng, em, X_lam, X_sig, lam_idx, sig_idx,
                                      hazard, S, site_off = NULL,
                                      fac_off = NULL) {
  if (eng$headroom < 0L) return(0)
  Ns <- nrow(X_lam)
  if (is.null(site_off)) site_off <- numeric(Ns)
  if (is.null(fac_off))  fac_off  <- matrix(0, Ns, S)
  eta_b <- if (hazard) em$global[1L] else 0
  worst <- 0
  for (s in seq_len(S)) {
    theta <- em$mu + em$b_list[[s]]
    el <- as.numeric(X_lam %*% theta[lam_idx]) + site_off + fac_off[, s]
    es <- as.numeric(X_sig %*% theta[sig_idx])
    sw_h <- cpp_distance_site_sweep(eng$y_s[[s]], el, es, eng$quad_xptr, eng$K_max,
                                    nb = FALSE, r = Inf, key = eng$key_code,
                                    eta_b = eta_b, headroom = eng$headroom)
    sw_u <- cpp_distance_site_sweep(eng$y_s[[s]], el, es, eng$quad_xptr, eng$K_max,
                                    nb = FALSE, r = Inf, key = eng$key_code,
                                    eta_b = eta_b, headroom = -1L)
    d <- .dist_score_gap(sw_h$grad_lam, sw_u$grad_lam, X_lam,
                         sw_h$grad_sig, sw_u$grad_sig, X_sig)
    if (is.finite(d) && d > worst) worst <- d
  }
  worst
}
