# =============================================================================
# Exact-marginal refinement for dynamic (multi-season HMM) occupancy
# (gcol33/tulpaObs#86).
#
# Like the single-season occupancy marginal (R/occu_marginal.R), the dynamic
# occupancy marginal is exact and cheap: the latent occupancy sequence
# integrates out by the HMM forward recursion -- `.tobs_ploglik_dynamic()`, the
# same kernel logLik() / WAIC / LOO score. The default deterministic engine fits
# this with a forward-backward E-step and a per-block pseudo-binomial Laplace
# M-step (the colonization / extinction blocks average over M = 1000 pseudo
# transitions). That encoding leaves a small EM-discretisation residual below
# the exact marginal MLE and reads the wrong curvature for the marginal, so the
# block standard errors are mis-scaled.
#
# Per tulpa's "nested approximation + debias" design, this is the debias step: a
# Newton refinement of the EM mode on the exact HMM-forward marginal
# log-posterior, with standard errors read from its Hessian. It moves the fixed
# effects onto the exact marginal MLE (matching unmarked::colext) and restores
# calibrated, near-nominal-coverage SEs. The refinement shares the
# optimise / vcov / draw-update / log-likelihood-refresh scaffold with the
# single-season path (`.tobs_marginal_refine_apply`); only the marginal
# objective differs (HMM forward vs the closed-form two-state mixture).
#
# Spatial (nested-Laplace) dynamic fits keep the EM field -- they carry a latent
# field, not a closed-form coefficient marginal -- so the refine is gated to the
# plain (non-spatial) Laplace path in `.tobs_laplace()`. Dynamic detection is
# site-level (the dyn_occu builder has no visit-level detection design), matching
# the HMM-forward kernel.
# =============================================================================


# Exact-marginal negative log-posterior for dynamic occupancy at the packed
# fixed-effect vector par = c(beta_psi1, beta_p, beta_gamma, beta_epsilon), in
# the process-block order `.tobs_ploglik_dynamic()` indexes. `pmean` / `pprec`
# are the Gaussian-prior mean and precision (1/sd^2, 0 = flat) aligned with par,
# so the refine maximises the same penalised marginal the EM M-steps target.
.tobs_dyn_occu_marginal_nlp <- function(par, model, pmean, pprec) {
  draw <- matrix(par, nrow = 1L)
  ll   <- .tobs_ploglik_dynamic(model, draw)   # [1 x n_sites]
  val  <- sum(ll[1L, ])
  if (!is.finite(val)) return(1e10)
  -val + 0.5 * sum(pprec * (par - pmean)^2)
}


# Refine a dynamic-occupancy `tobs_fit` on the exact HMM-forward marginal and
# overwrite its fixed-effect estimates + SEs + pseudo-draws with the debiased,
# Hessian-calibrated values. Falls back to the unmodified fit on any failure
# (never worse than the EM result).
.tobs_dyn_occu_marginal_refine <- function(fit, model, prior_spec = NULL) {
  tryCatch({
    pi_list <- model$process_info
    if (length(pi_list) < 4L) return(fit)

    # Fixed-effect names in process-block order (psi1, p, gamma, epsilon), the
    # layout fit$means / the draw matrix already carry.
    all_nm <- unlist(lapply(pi_list, function(p)
      paste0(p$name, "_", p$coef_names)), use.names = FALSE)

    # Gaussian-prior mean / precision per process block, aligned with par. The
    # bucket keying (psi1 / p / gamma / epsilon) matches `.attach_priors_to_blocks`
    # so the refine penalises exactly as the EM did; flat when prior_spec is NULL
    # (priors = FALSE), giving the exact marginal MLE.
    pr    <- lapply(pi_list, function(p)
      .prior_for_submodel(prior_spec, p$name, p$coef_names))
    pmean <- unlist(lapply(pr, `[[`, "mean"), use.names = FALSE)
    psd   <- unlist(lapply(pr, `[[`, "sd"),   use.names = FALSE)
    pprec <- ifelse(is.finite(psd) & psd > 0, 1 / psd^2, 0)

    nlp <- function(par) .tobs_dyn_occu_marginal_nlp(par, model, pmean, pprec)
    .tobs_marginal_refine_apply(fit, model, all_nm, nlp)
  }, error = function(e) fit)
}
