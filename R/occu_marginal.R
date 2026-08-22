# =============================================================================
# Closed-form marginal-likelihood refinement for single-season
# occupancy.
#
# For a single-season occupancy model the latent occupancy z integrates out in
# closed form, so the marginal likelihood is exact and cheap:
#
#   site with a detection : L_i = psi_i * prod_{valid j} Bern(y_ij; p_ij)
#   site with no detection: L_i = psi_i * prod_{valid j} (1 - p_ij) + (1 - psi_i)
#
# The default deterministic engine fits this via an EM whose occupancy M-step
# uses an M-inflated pseudo-binomial Laplace step; that encoding attenuates the
# detection coefficients at small J (an EM/M-step encoding artifact, not PQL) and
# -- because the M-step Hessian and the Louis observed info are the wrong
# curvature for the *marginal* likelihood -- under-disperses the detection-slope
# SE by ~3x
# regardless of N or J. Per tulpa's "nested approximation + debias" design, this
# is the debias step: a Newton refinement of the EM mode on the exact marginal
# log-posterior, with standard errors read from its Hessian. It restores
# unbiased point estimates and calibrated (near-nominal-coverage) SEs on the
# single-season path (verified against the analytic marginal MLE and the
# Monte-Carlo sd of the estimator). Spatial / RE / multi-season fits keep the
# EM path (no closed-form marginal); this refinement is single-season only.
# =============================================================================

# Exact marginal negative log-posterior for single-season occupancy at a packed
# parameter vector par = c(beta_occ, beta_det_site, beta_det_visit). `valid` is
# the n_sites x max_visits logical visit mask, `Y` the 0/1 detection matrix
# (invalid cells zeroed), `any_det` the per-site detection indicator. `pmean` /
# `pprec` are the Gaussian-prior mean and precision (1/sd^2, 0 = flat) aligned
# with par.
.tobs_occu_marginal_nlp <- function(par, X_occ, X_det, X_det_visit,
                                     valid, Y, any_det, n_sites, max_visits,
                                     p_occ, p_det, p_visit, pmean, pprec) {
  cl <- .tobs_clamp_eta
  bo <- par[seq_len(p_occ)]
  bd_site <- par[p_occ + seq_len(p_det)]
  psi <- plogis(cl(as.numeric(X_occ %*% bo)))
  eta_site <- as.numeric(X_det %*% bd_site)
  if (p_visit > 0L) {
    ev <- as.numeric(X_det_visit %*% par[p_occ + p_det + seq_len(p_visit)])
    logit_p <- matrix(eta_site, n_sites, max_visits) +
               matrix(ev, n_sites, max_visits, byrow = TRUE)
  } else {
    logit_p <- matrix(eta_site, n_sites, max_visits)
  }
  p <- plogis(cl(logit_p))
  logp   <- ifelse(valid, log(p),     0)
  log1mp <- ifelse(valid, log(1 - p), 0)
  det_ll   <- log(psi) + rowSums(Y * logp + (1 - Y) * log1mp)
  nodet_ll <- log(psi * exp(rowSums(log1mp)) + (1 - psi))
  ll <- sum(ifelse(any_det, det_ll, nodet_ll))
  penalty <- 0.5 * sum(pprec * (par - pmean)^2)
  -ll + penalty
}

# Shared scaffold for an exact-marginal Newton refinement of an EM mode. Given
# the fixed-effect names `all_nm` (the coefficients to refine, in fit$means /
# draw-column order) and an exact-marginal negative-log-posterior closure `nlp`
# of the packed parameter vector, optimise from the EM mode, read calibrated SEs
# from the Hessian, and overwrite the fit's fixed-effect estimates, SEs, and the
# matching pseudo-draws (drawn from the full joint covariance so derived
# quantities carry the coefficient correlation). The marginal log-likelihood
# logLik() reports is refreshed from the refined means through the same family
# pointwise kernel WAIC / LOO use. `refresh`, when supplied, is called as
# refresh(fit, par, V) to recompute any family-specific stored quantities (e.g.
# the per-site occupancy weights) at the refined mode. Falls back to the
# unmodified fit on any failure -- the refined fit is never worse than the EM
# result.
.tobs_marginal_refine_apply <- function(fit, model, all_nm, nlp, refresh = NULL) {
  tryCatch({
    if (!all(all_nm %in% names(fit$means))) return(fit)
    start <- as.numeric(fit$means[all_nm])
    if (any(!is.finite(start))) return(fit)

    # A tighter reltol than optim()'s default (relative to the function
    # VALUE, ~1e-8) matters when the starting point is farther from the
    # optimum along a shallow direction -- the default can pass a point whose
    # NLP value has stopped moving while a coordinate is still ~0.01 off:
    # int_occu()'s un-refined EM mode is farther from the exact-marginal
    # optimum than occu()'s own, and needed this to reach the same point
    # occu() already reached under the looser default.
    opt <- optim(start, nlp, method = "BFGS", hessian = TRUE,
                control = list(reltol = 1e-12, maxit = 500L))
    V <- tryCatch(solve(opt$hessian), error = function(e) NULL)
    if (is.null(V)) return(fit)
    se <- sqrt(pmax(diag(V), 0))
    if (any(!is.finite(c(opt$par, se)))) return(fit)

    fit$means[all_nm] <- opt$par
    fit$sds[all_nm]   <- se
    if (!is.null(fit$draws) && all(all_nm %in% colnames(fit$draws))) {
      Vsym <- (V + t(V)) / 2
      fit$draws[, all_nm] <- .rmvn(nrow(fit$draws), opt$par, Vsym)
    }
    if (is.function(refresh)) fit <- refresh(fit, opt$par, V)

    ml <- .tobs_laplace_marginal_loglik(model, fit$means)
    fit$log_lik  <- ml$loglik
    fit$log_prob <- rep(ml$loglik, length(fit$log_prob))
    fit
  }, error = function(e) fit)
}

# Refine a single-season occupancy `tobs_fit` on the exact marginal likelihood
# and overwrite its fixed-effect estimates + SEs (and the per-site occupancy
# weights / pseudo-draws) with the debiased, Hessian-calibrated values. Falls
# back to the unmodified fit on any failure (never worse than the EM result).
.tobs_occu_marginal_refine <- function(fit, model, prior_spec = NULL) {
  tryCatch({
    X_occ <- model$X_processes[[1]]
    X_det <- model$X_processes[[2]]
    X_det_visit <- model$X_det_visit
    y <- model$y
    n_sites <- nrow(X_occ); max_visits <- ncol(y)
    p_occ <- ncol(X_occ); p_det <- ncol(X_det)
    p_visit <- if (is.null(X_det_visit)) 0L else ncol(X_det_visit)

    valid <- y >= 0
    if (!any(valid)) return(fit)
    Y <- y; Y[!valid] <- 0
    any_det <- rowSums(Y * valid) > 0

    pi_list <- model$process_info
    occ_nm <- paste0(pi_list[[1]]$name, "_", pi_list[[1]]$coef_names)
    det_nm <- paste0(pi_list[[2]]$name, "_", pi_list[[2]]$coef_names)
    vis_nm <- if (p_visit > 0L) paste0("p_visit_", model$det_visit_names) else character(0)
    all_nm <- c(occ_nm, det_nm, vis_nm)

    # Gaussian-prior mean / precision aligned with par (flat when prior_spec
    # is NULL, i.e. priors = FALSE).
    pr_o <- .prior_for_submodel(prior_spec, "psi", pi_list[[1]]$coef_names)
    pr_d <- .prior_for_submodel(prior_spec, "p",
                                c(pi_list[[2]]$coef_names, model$det_visit_names))
    pmean <- c(pr_o$mean, pr_d$mean)
    psd   <- c(pr_o$sd,   pr_d$sd)
    pprec <- ifelse(is.finite(psd) & psd > 0, 1 / psd^2, 0)

    nlp <- function(par) .tobs_occu_marginal_nlp(
      par, X_occ = X_occ, X_det = X_det, X_det_visit = X_det_visit,
      valid = valid, Y = Y, any_det = any_det, n_sites = n_sites,
      max_visits = max_visits, p_occ = p_occ, p_det = p_det,
      p_visit = p_visit, pmean = pmean, pprec = pprec)

    # Refresh per-site occupancy weights P(z = 1 | y) at the refined estimate so
    # fitted / residuals stay consistent with the reported coefficients.
    refresh <- function(fit, par, V) {
      bo <- par[seq_len(p_occ)]
      bd <- par[p_occ + seq_len(p_det)]
      psi <- plogis(.tobs_clamp_eta(as.numeric(X_occ %*% bo)))
      eta_site <- as.numeric(X_det %*% bd)
      if (p_visit > 0L) {
        ev <- as.numeric(X_det_visit %*% par[p_occ + p_det + seq_len(p_visit)])
        logit_p <- matrix(eta_site, n_sites, max_visits) +
                   matrix(ev, n_sites, max_visits, byrow = TRUE)
      } else {
        logit_p <- matrix(eta_site, n_sites, max_visits)
      }
      p <- plogis(.tobs_clamp_eta(logit_p))
      log1mp <- ifelse(valid, log(1 - p), 0)
      prod0 <- exp(rowSums(log1mp))
      fit$weights <- ifelse(any_det, 1, psi * prod0 / (psi * prod0 + (1 - psi)))
      fit
    }

    .tobs_marginal_refine_apply(fit, model, all_nm, nlp, refresh)
  }, error = function(e) fit)
}
