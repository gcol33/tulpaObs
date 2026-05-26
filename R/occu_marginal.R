# =============================================================================
# Closed-form marginal-likelihood refinement for single-season occupancy
# (gcol33/tulpaObs#7).
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
  cl <- function(e) pmin(pmax(e, -30), 30)
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

# Refine a single-season occupancy `tobs_fit` on the exact marginal likelihood
# and overwrite its fixed-effect estimates + SEs (and the per-site occupancy
# weights / pseudo-draws) with the debiased, Hessian-calibrated values. Falls
# back to the unmodified fit on any failure (never worse than the EM result).
.tobs_occu_marginal_refine <- function(fit, model, prior_spec = NULL) {
  out <- tryCatch({
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
    if (!all(all_nm %in% names(fit$means))) return(fit)

    start <- as.numeric(fit$means[all_nm])
    if (any(!is.finite(start))) return(fit)

    # Gaussian-prior mean / precision aligned with par (flat when prior_spec
    # is NULL, i.e. priors = FALSE).
    pr_o <- .prior_for_submodel(prior_spec, "psi", pi_list[[1]]$coef_names)
    pr_d <- .prior_for_submodel(prior_spec, "p",
                                c(pi_list[[2]]$coef_names, model$det_visit_names))
    pmean <- c(pr_o$mean, pr_d$mean)
    psd   <- c(pr_o$sd,   pr_d$sd)
    pprec <- ifelse(is.finite(psd) & psd > 0, 1 / psd^2, 0)

    opt <- optim(start, .tobs_occu_marginal_nlp, method = "BFGS", hessian = TRUE,
                 X_occ = X_occ, X_det = X_det, X_det_visit = X_det_visit,
                 valid = valid, Y = Y, any_det = any_det, n_sites = n_sites,
                 max_visits = max_visits, p_occ = p_occ, p_det = p_det,
                 p_visit = p_visit, pmean = pmean, pprec = pprec)
    V <- tryCatch(solve(opt$hessian), error = function(e) NULL)
    if (is.null(V)) return(fit)
    se <- sqrt(pmax(diag(V), 0))
    if (any(!is.finite(c(opt$par, se)))) return(fit)

    # Overwrite the fixed-effect estimates, SEs, and matching pseudo-draws.
    fit$means[all_nm] <- opt$par
    fit$sds[all_nm]   <- se
    if (!is.null(fit$draws) && all(all_nm %in% colnames(fit$draws))) {
      n_pseudo <- nrow(fit$draws)
      for (j in seq_along(all_nm)) {
        fit$draws[, all_nm[j]] <- stats::rnorm(n_pseudo, opt$par[j], se[j])
      }
    }
    # Refresh per-site occupancy weights P(z = 1 | y) at the refined estimate so
    # fitted / residuals stay consistent with the reported coefficients.
    bo <- opt$par[seq_len(p_occ)]
    bd <- opt$par[p_occ + seq_len(p_det)]
    psi <- plogis(pmin(pmax(as.numeric(X_occ %*% bo), -30), 30))
    eta_site <- as.numeric(X_det %*% bd)
    if (p_visit > 0L) {
      ev <- as.numeric(X_det_visit %*% opt$par[p_occ + p_det + seq_len(p_visit)])
      logit_p <- matrix(eta_site, n_sites, max_visits) +
                 matrix(ev, n_sites, max_visits, byrow = TRUE)
    } else {
      logit_p <- matrix(eta_site, n_sites, max_visits)
    }
    p <- plogis(pmin(pmax(logit_p, -30), 30))
    log1mp <- ifelse(valid, log(1 - p), 0)
    prod0 <- exp(rowSums(log1mp))
    w <- ifelse(any_det, 1, psi * prod0 / (psi * prod0 + (1 - psi)))
    fit$weights <- w
    fit
  }, error = function(e) fit)
  out
}
