# =============================================================================
# Closed-form marginal-likelihood refinement for integrated (multi-source)
# occupancy (gcol33/tulpaObs#225).
#
# The latent occupancy z integrates out in closed form exactly as it does for
# single-season occu() (gcol33/tulpaObs#7) -- one shared z per site, observed
# through S independent detection sources:
#
#   site with ANY source detecting : L_i = psi_i * prod_s prod_{valid j} Bern(y_sij; p_sj)
#   site with no detection anywhere: L_i = psi_i * prod_s prod_{valid j} (1 - p_sj) + (1 - psi_i)
#
# .tobs_occu_marginal_refine() already debiases this EM/M-step artifact for
# model_type == "single" (and .tobs_dyn_occu_marginal_refine() for
# "dynamic"); model_type == "integrated" was missing this refinement
# entirely, so its reported fixed effects and SEs were the EM's raw
# (attenuated, mis-scaled) output. Diagnosed via #225: int_occu() with a
# single source did not numerically reduce to occu() on the same data even
# at full EM convergence with priors off, while the E-step weights and
# M-step block encodings matched between the two paths to numerical
# precision -- ruling out the EM machinery and pointing at this exact
# refinement occu() gets and int_occu() did not.
# =============================================================================

# Exact marginal negative log-posterior for integrated occupancy at a packed
# parameter vector par = c(beta_occ, beta_det_1, ..., beta_det_S). `site_maps`
# is the list of 0-based row -> site maps per source (site_maps[[s]][r] + 1 =
# the site row r of y_sources[[s]] belongs to); `y_sources` the per-source
# 0/1/-1 detection matrices. `pmean` / `pprec` are the Gaussian-prior mean and
# precision (1/sd^2, 0 = flat) aligned with par.
.tobs_integrated_marginal_nlp <- function(par, X_occ, X_det_list, y_sources,
                                          site_maps, n_sites, p_occ, p_det,
                                          pmean, pprec) {
  cl <- .tobs_clamp_eta
  bo <- par[seq_len(p_occ)]
  psi <- plogis(cl(as.numeric(X_occ %*% bo)))

  any_det <- logical(n_sites)
  det_ll_contrib   <- numeric(n_sites)   # sum_s sum_j y log p + (1-y) log(1-p)
  log1mp_sum       <- numeric(n_sites)   # sum_s sum_valid_j log(1 - p_sj)
  off <- p_occ
  n_src <- length(X_det_list)
  for (s in seq_len(n_src)) {
    ps <- p_det[s]
    bd <- par[off + seq_len(ps)]
    off <- off + ps
    eta <- as.numeric(X_det_list[[s]] %*% bd)
    p_s <- plogis(cl(eta))
    ys  <- y_sources[[s]]
    rows <- seq_len(nrow(ys))
    smap <- site_maps[[s]] + 1L
    for (r in rows) {
      i <- smap[r]
      yr <- ys[r, ]
      valid <- yr >= 0
      if (!any(valid)) next
      if (any(yr[valid] == 1)) any_det[i] <- TRUE
      pr <- p_s[r]
      pr_v <- rep(pr, sum(valid))
      yv <- yr[valid]
      det_ll_contrib[i] <- det_ll_contrib[i] +
        sum(yv * log(pr_v) + (1 - yv) * log(1 - pr_v))
      log1mp_sum[i] <- log1mp_sum[i] + sum(log(1 - pr_v))
    }
  }
  det_ll   <- log(psi) + det_ll_contrib
  nodet_ll <- log(psi * exp(log1mp_sum) + (1 - psi))
  ll <- sum(ifelse(any_det, det_ll, nodet_ll))
  penalty <- 0.5 * sum(pprec * (par - pmean)^2)
  -ll + penalty
}

# Refine an integrated-occupancy `tobs_fit` on the exact marginal likelihood
# and overwrite its fixed-effect estimates + SEs (and the per-site occupancy
# weights / pseudo-draws) with the debiased, Hessian-calibrated values. Falls
# back to the unmodified fit on any failure (never worse than the EM result).
.tobs_integrated_marginal_refine <- function(fit, model, prior_spec = NULL) {
  tryCatch({
    X_occ <- model$X_processes[[1]]
    n_sites <- model$n_sites
    n_src <- model$n_sources
    p_occ <- ncol(X_occ)

    X_det_list <- lapply(seq_len(n_src), function(s) model$X_processes[[1 + s]])
    p_det <- vapply(X_det_list, ncol, integer(1))

    pi_list <- model$process_info
    occ_nm <- paste0(pi_list[[1]]$name, "_", pi_list[[1]]$coef_names)
    det_nm <- unlist(lapply(seq_len(n_src), function(s)
      paste0(pi_list[[1 + s]]$name, "_", pi_list[[1 + s]]$coef_names)))
    all_nm <- c(occ_nm, det_nm)
    if (!all(all_nm %in% names(fit$means))) return(fit)

    pr_o <- .prior_for_submodel(prior_spec, "psi", pi_list[[1]]$coef_names)
    pr_list <- lapply(seq_len(n_src), function(s)
      .prior_for_submodel(prior_spec, paste0("det", s), pi_list[[1 + s]]$coef_names))
    pmean <- c(pr_o$mean, unlist(lapply(pr_list, `[[`, "mean")))
    psd   <- c(pr_o$sd,   unlist(lapply(pr_list, `[[`, "sd")))
    pprec <- ifelse(is.finite(psd) & psd > 0, 1 / psd^2, 0)

    nlp <- function(par) .tobs_integrated_marginal_nlp(
      par, X_occ = X_occ, X_det_list = X_det_list, y_sources = model$y_sources,
      site_maps = model$site_maps, n_sites = n_sites, p_occ = p_occ,
      p_det = p_det, pmean = pmean, pprec = pprec)

    # Refresh per-site occupancy weights P(z = 1 | y) at the refined estimate
    # so fitted / residuals stay consistent with the reported coefficients.
    refresh <- function(fit, par, V) {
      bo <- par[seq_len(p_occ)]
      psi <- plogis(.tobs_clamp_eta(as.numeric(X_occ %*% bo)))
      any_det <- logical(n_sites)
      log1mp_sum <- numeric(n_sites)
      off <- p_occ
      for (s in seq_len(n_src)) {
        ps <- p_det[s]
        bd <- par[off + seq_len(ps)]
        off <- off + ps
        p_s <- plogis(.tobs_clamp_eta(as.numeric(X_det_list[[s]] %*% bd)))
        ys <- model$y_sources[[s]]
        smap <- model$site_maps[[s]] + 1L
        for (r in seq_len(nrow(ys))) {
          i <- smap[r]
          yr <- ys[r, ]
          valid <- yr >= 0
          if (!any(valid)) next
          if (any(yr[valid] == 1)) any_det[i] <- TRUE
          log1mp_sum[i] <- log1mp_sum[i] + sum(log(1 - rep(p_s[r], sum(valid))))
        }
      }
      prod0 <- exp(log1mp_sum)
      fit$weights <- ifelse(any_det, 1, psi * prod0 / (psi * prod0 + (1 - psi)))
      fit
    }

    .tobs_marginal_refine_apply(fit, model, all_nm, nlp, refresh)
  }, error = function(e) fit)
}
