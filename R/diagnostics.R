# ============================================================================
# Occupancy-specific diagnostics
# Generic diagnostics (moran_i, durbin_watson, variogram, compare_models,
# modelAverage) are in tulpa — inherited via tulpa_fit class.
# ============================================================================

#' Model criteria for occupancy / cover models
#'
#' `tobs_waic()`, `tobs_dic()`, and `tobs_cpo()` build the family pointwise
#' log-likelihood matrix once and hand it to the engine's single criteria layer
#' [tulpa::tulpa_criteria()], which derives WAIC / DIC / CPO / LPML / PSIS-LOO.
#' DIC additionally evaluates the deviance at the posterior mean of the
#' parameters, supplied by the family-specific `.tobs_loglik_at_mean()`.
#'
#' @param object A `tobs_fit` object.
#' @param n.draws Posterior draws used to build the pointwise log-likelihood
#'   (the cover / occu_cover paths sample this many; the draw-matrix families use
#'   the first `n.draws` stored draws). Default 1000.
#' @param ... Forwarded to [tulpa::tulpa_criteria()] (e.g. `chunk_size`).
#' @return A `tulpa_criteria` object. `tobs_waic()` also carries `$elpd`
#'   (an alias for `elpd_waic`) for back-compatibility.
#' @seealso [tulpa::tulpa_criteria()]
#' @export
tobs_waic <- function(object, n.draws = 1000L, ...) {
  ll_mat <- .tobs_pointwise_loglik(object, n.draws = n.draws)
  cr <- tulpa::tulpa_criteria(ll_mat, criteria = "waic", ...)
  cr$elpd <- cr$elpd_waic
  cr
}

#' @rdname tobs_waic
#' @export
tobs_dic <- function(object, n.draws = 1000L, ...) {
  ll_mat <- .tobs_pointwise_loglik(object, n.draws = n.draws)
  lam <- .tobs_loglik_at_mean(object, n.draws = n.draws)
  tulpa::tulpa_criteria(ll_mat, criteria = "dic", loglik_at_mean = lam, ...)
}

#' @rdname tobs_waic
#' @export
tobs_cpo <- function(object, n.draws = 1000L, ...) {
  ll_mat <- .tobs_pointwise_loglik(object, n.draws = n.draws)
  tulpa::tulpa_criteria(ll_mat, criteria = c("loo", "cpo", "lpml"),
                        pointwise = TRUE, ...)
}

# Pointwise log-likelihood matrix [n_draws x n_obs], marginalized over the
# latent state, for every family. This is the input WAIC, PSIS-LOO, and LOO
# stacking all consume, so it lives in one place (tobs_stack reuses it). The
# per-family marginal likelihoods mirror the C++ kernels in src/*_likelihood.h.
#
# Fidelity note: the linear predictors are evaluated from the *process fixed-
# effect* coefficient draws (the leading columns of `draws`), exactly as
# `fitted()` and the historical WAIC do. Structured-term contributions
# (spatial / temporal / random-effect fields) are not added to the predictor,
# and dynamic visit-level detection covariates are not folded in. For models
# with those components the score is conditional on the fixed-effect predictor.
.tobs_pointwise_loglik <- function(object, n.draws = NULL) {
  nd <- n.draws %||% 1000L
  if (inherits(object, "cover_fit")) return(.tobs_ploglik_cover(object, nd))
  if (identical(object$model$model_type %||% "NULL", "occu_cover")) {
    return(.tobs_ploglik_occu_cover(object, nd))
  }
  # Three-level occupancy + cover: the WAIC / LOO unit is the cell (the top-level
  # marginalised observation); scored over the draw matrix (calibrated under the
  # NUTS path, pseudo-draws under Laplace).
  if (identical(object$model$model_type %||% "NULL", "occu_multiscale_cover")) {
    return(.tobs_ploglik_occu_multiscale_cover(object, nd))
  }
  # Spatial-factor community occu_cover: the per-cell likelihood needs the latent
  # field, so it is scored over the NUTS draws (calibrated WAIC / LOO -- the point
  # of the NUTS path).
  if (identical(object$model$model_type %||% "NULL", "ms_occu_cover_spatial")) {
    return(.tobs_ploglik_ms_occu_cover_spatial(object, nd))
  }
  # Community N-mixture (ms_abun): the per-(species, site) likelihood needs the
  # per-species deviations, so it is scored over the NUTS draws (the calibrated
  # WAIC / LOO the NUTS path enables; the Laplace community-mean draws omit them).
  if (identical(object$model$model_type %||% "NULL", "ms_nmix")) {
    return(.tobs_ploglik_ms_nmix(object, nd))
  }

  draws <- object$draws
  if (is.null(draws) || !is.matrix(draws)) {
    stop("Pointwise log-likelihood needs a posterior draw matrix; ",
         "`object$draws` is missing or not a matrix.", call. = FALSE)
  }
  if (!is.null(n.draws) && n.draws < nrow(draws)) {
    draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
  }
  .tobs_ploglik_from_draws(object$model, draws)
}

# Per-family pointwise log-likelihood given an explicit [n_draws x p] draw
# matrix. Split out from the dispatcher so the posterior-mean evaluation
# (.tobs_loglik_at_mean) drives the same per-family kernels with a one-row mean.
.tobs_ploglik_from_draws <- function(model, draws) {
  mt <- model$model_type %||% "NULL"
  switch(
    mt,
    single     = .tobs_ploglik_replicated(model, draws),
    dynamic    = .tobs_ploglik_dynamic(model, draws),
    integrated = .tobs_ploglik_integrated(model, draws),
    jsdm       = .tobs_ploglik_jsdm(model, draws),
    nmix       = .tobs_ploglik_nmix(model, draws),
    removal    = .tobs_ploglik_removal(model, draws),
    distance   = .tobs_ploglik_distance(model, draws),
    fp_occu    = .tobs_ploglik_fp_occu(model, draws),
    dyn_abun   = .tobs_ploglik_dyn_abun(model, draws),
    occu_multiscale_cover = {
      idx <- .tobs_occu_ms_cover_nuts_layout(model)
      t(apply(draws, 1L, function(par)
        .occu_ms_cover_nonspatial_ll(par, model, idx, per_cell = TRUE)))
    },
    stop("Pointwise log-likelihood is not implemented for model_type = '",
         mt, "'.", call. = FALSE)
  )
}

# Pointwise log-likelihood at the posterior mean of the parameters, the plug-in
# DIC needs (length n_obs). The draw-matrix families evaluate their per-family
# kernel at the column-mean draw; the cover / occu_cover families plug in the
# posterior-mean linear predictors (and mean dispersion) via family helpers.
.tobs_loglik_at_mean <- function(object, n.draws = 1000L) {
  if (inherits(object, "cover_fit")) {
    return(.tobs_cover_loglik_at_mean(object, n.draws))
  }
  if (identical(object$model$model_type %||% "NULL", "occu_cover")) {
    return(.tobs_occu_cover_loglik_at_mean(object, n.draws))
  }
  draws <- object$draws
  if (is.null(draws) || !is.matrix(draws)) {
    stop("DIC needs a posterior draw matrix to evaluate the deviance at the ",
         "posterior mean; `object$draws` is missing or not a matrix.",
         call. = FALSE)
  }
  mean_draw <- matrix(colMeans(draws), nrow = 1L)
  as.numeric(.tobs_ploglik_from_draws(object$model, mean_draw))
}

# Total marginal log-likelihood at a fixed-effect point estimate -- the value
# logLik() / AIC() / BIC() / glance() report for an EM+Laplace fit (single,
# dynamic, integrated, jsdm). Reuses the family pointwise kernel
# (.tobs_ploglik_from_draws) on a one-row draw matrix, so the reported marginal
# is the same likelihood the WAIC / LOO scoring uses -- one source of truth. The
# unit summed over is the model's marginal observation (the site, with its
# visits / seasons pooled into the closed-form or HMM-forward marginal), so
# `nobs` is the number of those units. `par` is the named fixed-effect vector in
# the process-block layout the draw matrix uses (betas concatenated in process
# order; any trailing visit / random-effect columns are ignored by the kernel,
# matching WAIC). Returns NA when the model_type has no pointwise kernel or the
# evaluation errors, so a fit never fails to assemble over a logLik plumbing
# detail (gcol33/tulpaObs#87).
.tobs_laplace_marginal_loglik <- function(model, par) {
  tryCatch({
    draw <- matrix(as.numeric(par), nrow = 1L)
    colnames(draw) <- names(par)
    ll  <- .tobs_ploglik_from_draws(model, draw)
    val <- sum(ll[1L, ])
    list(loglik = if (is.finite(val)) val else NA_real_,
         nobs   = ncol(ll))
  }, error = function(e) list(loglik = NA_real_, nobs = NA_integer_))
}

# --- shared numerics --------------------------------------------------------

# log(inv_logit(eta)) = log(p) and log(1 - inv_logit(eta)) = log(1 - p),
# computed stably via plogis(log.p) -- the R analogues of the C++
# log_inv_logit / log1m_inv_logit.
.tobs_log_p   <- function(eta) stats::plogis(eta,  log.p = TRUE)
.tobs_log_1mp <- function(eta) stats::plogis(-eta, log.p = TRUE)

# Numerically stable elementwise log(exp(a) + exp(b)).
.tobs_logaddexp <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  both_inf <- is.infinite(m) & (a == b)
  out[both_inf] <- m[both_inf]
  out
}

# Cumulative column offset of process k's beta block within the draw matrix
# (processes are laid out in order, fixed-effect betas first).
.tobs_proc_offset <- function(model, k) {
  if (k == 1L) return(0L)
  sum(vapply(model$process_info[seq_len(k - 1L)],
             function(p) p$p, integer(1)))
}

# Linear-predictor draws for process k: [n_draws x nrow(X_k)].
.tobs_eta_draws <- function(model, draws, k) {
  p_k    <- model$process_info[[k]]$p
  off    <- .tobs_proc_offset(model, k)
  beta_k <- draws[, off + seq_len(p_k), drop = FALSE]
  beta_k %*% t(model$X_processes[[k]])
}

# --- per-family marginal likelihoods ---------------------------------------

# Single-season occupancy: per replicate row, marginalized over z.
.tobs_ploglik_replicated <- function(model, draws) {
  eta_psi <- .tobs_eta_draws(model, draws, 1L)   # [S x n_obs]
  eta_p   <- .tobs_eta_draws(model, draws, 2L)
  y       <- model$y                             # [n_obs x max_visits], <0 = NA
  S <- nrow(eta_psi); n_obs <- nrow(y)
  ll <- matrix(0, S, n_obs)
  for (i in seq_len(n_obs)) {
    yi <- y[i, ]; valid <- yi >= 0; nv <- sum(valid)
    if (nv == 0L) next                           # no data -> 0 contribution
    log_p   <- .tobs_log_p(eta_p[, i]);   log_1mp <- .tobs_log_1mp(eta_p[, i])
    log_psi <- .tobs_log_p(eta_psi[, i]); log1m_psi <- .tobs_log_1mp(eta_psi[, i])
    k1 <- sum(yi[valid] == 1); k0 <- nv - k1
    if (k1 > 0L) {
      ll[, i] <- log_psi + k1 * log_p + k0 * log_1mp
    } else {
      ll[, i] <- .tobs_logaddexp(log_psi + nv * log_1mp, log1m_psi)
    }
  }
  ll
}

# JSDM: per site x species, y ~ Bernoulli(psi). No detection, no marginalization.
.tobs_ploglik_jsdm <- function(model, draws) {
  eta <- .tobs_eta_draws(model, draws, 1L)       # [S x N]
  y   <- model$y_jsdm
  Y   <- matrix(y, nrow(eta), length(y), byrow = TRUE)
  Y * .tobs_log_p(eta) + (1 - Y) * .tobs_log_1mp(eta)
}

# N-mixture: per site, the latent abundance N integrated out in closed form (the
# Royle 2004 marginal). The observation unit is the site (the per-site marginal
# pools that site's visits), so the pointwise log-likelihood is [n_draws x
# n_sites]. NB is detected by the trailing log_r draw column; the per-draw size is
# r = exp(log_r). Reuses the same nmix_site_marginal() kernel the fit used, so the
# WAIC / LOO scoring is on one source of truth. (For NUTS fits the draws are the
# exact posterior; the laplace path's Gaussian draws also score, the N-mixture
# coefficient marginal being well-behaved.)
.tobs_ploglik_nmix <- function(model, draws) {
  X_lambda <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
  p_lam <- ncol(X_lambda); p_p <- ncol(X_p)
  is_nb <- ("log_r" %in% colnames(draws)) || (ncol(draws) > p_lam + p_p)
  K_max <- as.integer(max(model$y_long) + 100L)
  marg  <- nmix_site_marginal(y = model$y_long, site_idx = model$site_idx,
                              X_lambda = X_lambda, X_p = X_p,
                              mixture = if (is_nb) "NB" else "P", K_max = K_max)
  S <- nrow(draws); n_sites <- model$n_sites
  ll <- matrix(0, S, n_sites)
  for (s in seq_len(S)) {
    bl <- draws[s, seq_len(p_lam)]
    bp <- draws[s, p_lam + seq_len(p_p)]
    r  <- if (is_nb) exp(draws[s, p_lam + p_p + 1L]) else Inf
    ll[s, ] <- marg$eval_beta(bl, bp, r = r)$log_lik_site
  }
  ll
}

# Removal sampling: per site, the latent abundance N integrated out in closed
# form over the depleting-binomial removal likelihood (the same marginal the fit
# used). The observation unit is the site (its passes are pooled), so the
# pointwise log-likelihood is [n_draws x n_sites]. NB is detected by the trailing
# log_r draw column.
.tobs_ploglik_removal <- function(model, draws) {
  X_lambda <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
  p_lam <- ncol(X_lambda); p_p <- ncol(X_p)
  is_nb <- ("log_r" %in% colnames(draws)) || (ncol(draws) > p_lam + p_p)
  marg  <- .tobs_removal_nuts_marginal(model, mixture = if (is_nb) "NB" else "P")
  S <- nrow(draws); n_sites <- model$n_sites
  ll <- matrix(0, S, n_sites)
  for (s in seq_len(S)) {
    bl <- draws[s, seq_len(p_lam)]
    bp <- draws[s, p_lam + seq_len(p_p)]
    r  <- if (is_nb) exp(draws[s, p_lam + p_p + 1L]) else Inf
    ll[s, ] <- marg$eval_beta(bl, bp, r = r)$log_lik_site
  }
  ll
}

# Distance sampling: per site, the latent abundance N integrated out in closed
# form over the binned multinomial-over-N detection likelihood (the same marginal
# the fit used). The observation unit is the site (its bins are pooled), so the
# pointwise log-likelihood is [n_draws x n_sites]. The hazard-rate shape is read
# from the model key; NB is detected by the trailing log_r draw column.
.tobs_ploglik_distance <- function(model, draws) {
  p_lam <- model$process_info[[1]]$p; p_sig <- model$process_info[[2]]$p
  hazard <- identical(model$key, "hazard")
  is_nb  <- "log_r" %in% colnames(draws)
  marg   <- .tobs_distance_nuts_marginal(model, mixture = if (is_nb) "NB" else "P")
  S <- nrow(draws); n_sites <- model$n_sites
  off <- p_lam + p_sig
  ll <- matrix(0, S, n_sites)
  for (s in seq_len(S)) {
    bl <- draws[s, seq_len(p_lam)]
    bs <- draws[s, p_lam + seq_len(p_sig)]
    eb <- if (hazard) draws[s, off + 1L] else 0
    r  <- if (is_nb) exp(draws[s, off + (if (hazard) 2L else 1L)]) else Inf
    ll[s, ] <- marg$eval_beta(bl, bs, eta_b = eb, r = r)$log_lik_site
  }
  ll
}

# False-positive occupancy: per site, the latent occupancy z integrated out in
# closed form over the Miller et al. (2011) multistate marginal. The observation
# unit is the site (its visits are pooled), so the pointwise log-likelihood is
# [n_draws x n_sites]. Coefficient layout: (psi, p11, p10, b) site-level arms.
.tobs_ploglik_fp_occu <- function(model, draws) {
  lay <- .tobs_fp_occu_nuts_layout(model$process_info[[1]]$p,
                                   model$process_info[[2]]$p,
                                   model$process_info[[3]]$p,
                                   model$process_info[[4]]$p)
  marg <- .tobs_fp_occu_nuts_marginal(model)
  S <- nrow(draws); n_sites <- model$n_sites
  ll <- matrix(0, S, n_sites)
  for (s in seq_len(S)) {
    ev <- marg$eval_beta(draws[s, lay$psi], draws[s, lay$p11],
                         draws[s, lay$p10], draws[s, lay$b])
    ll[s, ] <- ev$log_lik_site
  }
  ll
}

# Open N-mixture (dyn_abun): per site, the latent abundance sequence integrated
# out by the HMM forward recursion (the same marginal the fit used). The
# observation unit is the site (its seasons / visits are pooled), so the pointwise
# log-likelihood is [n_draws x n_sites]. Layout: (lambda, p, omega, gamma) arms.
.tobs_ploglik_dyn_abun <- function(model, draws) {
  lay <- .tobs_dyn_abun_nuts_layout(model$process_info[[1]]$p,
                                    model$process_info[[2]]$p,
                                    model$process_info[[3]]$p,
                                    model$process_info[[4]]$p)
  marg <- .tobs_dyn_abun_nuts_marginal(model)
  S <- nrow(draws); n_sites <- model$n_sites
  ll <- matrix(0, S, n_sites)
  for (s in seq_len(S)) {
    ev <- marg$eval_beta(draws[s, lay$lambda], draws[s, lay$p],
                         draws[s, lay$omega], draws[s, lay$gamma])
    ll[s, ] <- ev$log_lik_site
  }
  ll
}

# Community N-mixture (ms_abun): per (species, site), the latent abundance N
# integrated out in closed form, scored over the NUTS draws of the FULL parameter
# vector (community means + per-species deviations + community covariances). The
# pointwise unit is a (species, site): the per-species marginal pools that
# species-site's visits, so the matrix is [n_draws x (n_species * n_sites)] with
# the per-species blocks laid contiguously. NB is read from the layout. Needs a
# NUTS fit (object$nuts): the Laplace community-mean draws omit the per-species
# deviations, so the per-(species, site) likelihood is not identified from them
# (the same constraint as the spatial-factor community occu_cover path).
.tobs_ploglik_ms_nmix <- function(object, n.draws = 1000L) {
  nd <- object$nuts
  if (is.null(nd) || is.null(nd$draws)) {
    stop("WAIC / LOO for the community N-mixture ms_abun() needs a NUTS fit ",
         "(method = \"nuts\"): the Laplace community-mean draws omit the ",
         "per-species deviations, so the per-(species, site) likelihood is not ",
         "identified from them.", call. = FALSE)
  }
  model <- object$model
  lay   <- nd$layout
  is_nb <- isTRUE(lay$is_nb)
  lf    <- .tobs_ms_nmix_longform(model)
  K_max <- nd$K_max %||% as.integer(max(lf$y) + 100L)
  margs <- .tobs_ms_abun_nuts_marginals(lf, model$X_processes[[1]],
                                        model$n_sites,
                                        if (is_nb) "NB" else "P", K_max)
  draws <- nd$draws
  M <- nrow(draws)
  if (!is.null(n.draws) && as.integer(n.draws) < M) {
    idx <- unique(round(seq(1, M, length.out = as.integer(n.draws))))
    draws <- draws[idx, , drop = FALSE]; M <- nrow(draws)
  }
  S <- lay$n_species; n_sites <- model$n_sites
  out <- matrix(0, M, S * n_sites)
  for (m in seq_len(M)) {
    mu <- draws[m, lay$mu]
    # Non-centered draws store the whitened z; reconstruct b = C z per draw.
    B <- .tobs_ms_abun_nuts_b_from_z(draws[m, ], lay)
    for (s in seq_len(S)) {
      b_s <- B[s, ]
      r   <- if (is_nb) exp(mu[lay$logr] + b_s[lay$logr]) else Inf
      ev  <- margs[[s]]$eval_beta(mu[lay$lambda] + b_s[lay$lambda],
                                  mu[lay$p]      + b_s[lay$p], r = r)
      out[m, (s - 1L) * n_sites + seq_len(n_sites)] <- ev$log_lik_site
    }
  }
  out
}

# Integrated multi-source: per site, shared psi, detection summed over the
# sources that observed it (src/integrated_occ_likelihood.h).
.tobs_ploglik_integrated <- function(model, draws) {
  eta_psi <- .tobs_eta_draws(model, draws, 1L)   # [S x n_sites]
  S <- nrow(eta_psi); n_sites <- model$n_sites
  n_sources <- model$n_sources
  log_psi   <- .tobs_log_p(eta_psi); log1m_psi <- .tobs_log_1mp(eta_psi)
  # per-source log(p) / log(1-p) at every site (valid only where observed)
  log_p_src   <- vector("list", n_sources)
  log_1mp_src <- vector("list", n_sources)
  for (s in seq_len(n_sources)) {
    eta_s <- .tobs_eta_draws(model, draws, 1L + s)
    log_p_src[[s]]   <- .tobs_log_p(eta_s)
    log_1mp_src[[s]] <- .tobs_log_1mp(eta_s)
  }

  ll <- matrix(0, S, n_sites)
  for (i in seq_len(n_sites)) {
    log_det_occ <- numeric(S)   # sum_s log P(y_is | occupied)
    any_det <- FALSE
    for (s in seq_len(n_sources)) {
      local <- which(model$site_maps[[s]] + 1L == i)
      if (!length(local)) next
      yvec <- model$y_sources[[s]][local[1L], ]
      valid <- yvec >= 0; nv <- sum(valid)
      if (nv == 0L) next
      k1 <- sum(yvec[valid] == 1); k0 <- nv - k1
      log_det_occ <- log_det_occ + k1 * log_p_src[[s]][, i] +
                                   k0 * log_1mp_src[[s]][, i]
      if (k1 > 0L) any_det <- TRUE
    }
    if (any_det) {
      ll[, i] <- log_psi[, i] + log_det_occ
    } else {
      # log_det_occ here is sum_s nv_s * log(1-p_s) (all-zero across sources)
      ll[, i] <- .tobs_logaddexp(log_psi[, i] + log_det_occ, log1m_psi[, i])
    }
  }
  ll
}

# Dynamic (multi-season HMM): per-site forward recursion in log space,
# mirroring src/dyn_occ_likelihood.h. Site-level detection only.
.tobs_ploglik_dynamic <- function(model, draws) {
  eta_psi1 <- .tobs_eta_draws(model, draws, 1L)
  eta_p    <- .tobs_eta_draws(model, draws, 2L)
  eta_gam  <- .tobs_eta_draws(model, draws, 3L)
  eta_eps  <- .tobs_eta_draws(model, draws, 4L)
  S <- nrow(eta_psi1); n_sites <- model$n_sites; Tn <- model$n_seasons
  n_visits <- model$n_visits; any_det <- model$any_detected
  NEG <- -1e10

  out <- matrix(0, S, n_sites)
  for (i in seq_len(n_sites)) {
    lp    <- .tobs_log_p(eta_p[, i]);    l1mp   <- .tobs_log_1mp(eta_p[, i])
    lgam  <- .tobs_log_p(eta_gam[, i]);  l1mgam <- .tobs_log_1mp(eta_gam[, i])
    leps  <- .tobs_log_p(eta_eps[, i]);  l1meps <- .tobs_log_1mp(eta_eps[, i])
    a_occ <- .tobs_log_p(eta_psi1[, i]); a_un   <- .tobs_log_1mp(eta_psi1[, i])
    site_ll <- numeric(S)
    for (t in seq_len(Tn)) {
      idx <- (i - 1L) * Tn + (t - 1L) + 1L
      nv  <- n_visits[idx]
      if (nv > 0L) {
        yvec <- model$y[i, , t]; yvec[is.na(yvec)] <- -1L
        valid <- yvec >= 0
        k1 <- sum(yvec[valid] == 1); k0 <- sum(yvec[valid] == 0)
        if (any_det[idx]) {
          site_ll <- site_ll + a_occ + (k1 * lp + k0 * l1mp)
          a_occ <- numeric(S); a_un <- rep(NEG, S)   # z_t = 1 known
        } else {
          term1 <- a_occ + nv * l1mp                 # occupied, all non-detections
          term2 <- a_un                              # unoccupied
          lnorm <- .tobs_logaddexp(term1, term2)
          site_ll <- site_ll + lnorm
          a_occ <- term1 - lnorm
          a_un  <- term2 - lnorm
        }
      }
      if (t < Tn) {
        new_occ <- .tobs_logaddexp(a_occ + l1meps, a_un + lgam)
        new_un  <- .tobs_logaddexp(a_occ + leps,   a_un + l1mgam)
        a_occ <- new_occ; a_un <- new_un
      }
    }
    out[, i] <- site_ll
  }
  out
}

#' Posterior predictive check
#' @param object A `tobs_fit` object.
#' @param fit.stat `"freeman-tukey"` (default) or `"chi-squared"`.
#' @param n.samples Number of posterior samples (default 500).
#' @return A list with `fit.y`, `fit.y.rep`, and `bayesian.p`.
#' @export
tobs_ppc <- function(object, fit.stat = c("freeman-tukey", "chi-squared"),
                     n.samples = 500) {
  fit.stat <- match.arg(fit.stat)
  if (inherits(object, "cover_fit")) {
    return(.tobs_ppc_cover(object, fit.stat, n.samples))
  }
  if (identical(object$model$model_type %||% "NULL", "occu_cover")) {
    return(.tobs_ppc_occu_cover(object, fit.stat, n.samples))
  }
  model <- object$model
  if (model$model_type != "single") {
    stop("tobs_ppc supports single-season occupancy, cover(), and occu_cover() ",
         "fits.", call. = FALSE)
  }

  draws <- object$draws; pi_list <- model$process_info
  X_occ <- model$X_processes[[1]]; X_det <- model$X_processes[[2]]
  y <- model$y; n_sites <- model$n_sites; max_visits <- model$max_visits
  p_occ <- pi_list[[1]]$p; p_det <- pi_list[[2]]$p
  n.samples <- min(n.samples, nrow(draws))
  draw_idx <- sample.int(nrow(draws), n.samples)

  stat_fn <- if (fit.stat == "freeman-tukey") {
    function(obs, exp) sum((sqrt(obs) - sqrt(exp))^2, na.rm = TRUE)
  } else {
    function(obs, exp) sum((obs - exp)^2 / (exp + 1e-10), na.rm = TRUE)
  }

  # Per-site valid mask, visit count, and whether the species was ever
  # detected. The latent z is sampled from its full conditional given this
  # detection history (the spOccupancy ppcOcc construction), not from the
  # prior psi: a site with any detection is occupied with probability 1, an
  # all-zero history occupied with probability
  #   psi (1-p)^J / [psi (1-p)^J + (1-psi)].
  valid_mat <- y >= 0
  n_valid   <- rowSums(valid_mat)
  any_det   <- rowSums(y * valid_mat) > 0

  fit_y <- fit_y_rep <- numeric(n.samples)
  for (s in seq_len(n.samples)) {
    idx <- draw_idx[s]
    psi <- plogis(as.vector(X_occ %*% draws[idx, seq_len(p_occ)]))
    p   <- plogis(as.vector(X_det %*% draws[idx, p_occ + seq_len(p_det)]))
    z_prob <- ifelse(
      any_det, 1,
      psi * (1 - p)^n_valid / (psi * (1 - p)^n_valid + (1 - psi))
    )
    z_prob[n_valid == 0L] <- psi[n_valid == 0L]   # no data -> prior
    z <- rbinom(n_sites, 1, z_prob)
    expected <- y_rep <- matrix(NA, n_sites, max_visits)
    for (i in seq_len(n_sites)) for (j in seq_len(max_visits)) if (y[i,j] >= 0) {
      expected[i,j] <- z[i] * p[i]; y_rep[i,j] <- rbinom(1, 1, z[i] * p[i])
    }
    valid <- !is.na(expected)
    fit_y[s] <- stat_fn(y[valid], expected[valid])
    fit_y_rep[s] <- stat_fn(y_rep[valid], expected[valid])
  }
  list(fit.y = fit_y, fit.y.rep = fit_y_rep, bayesian.p = mean(fit_y_rep > fit_y))
}

#' PIT residuals
#' @param object A `tobs_fit` object.
#' @param n.samples Number of posterior samples (default 250).
#' @return Numeric vector of PIT residuals.
#' @export
tobs_pit_residuals <- function(object, n.samples = 250) {
  if (inherits(object, "cover_fit")) {
    return(.tobs_pit_cover(object, n.samples))
  }
  if (identical(object$model$model_type %||% "NULL", "occu_cover")) {
    return(.tobs_pit_occu_cover(object, n.samples))
  }
  model <- object$model
  if (model$model_type != "single") {
    stop("tobs_pit_residuals supports single-season occupancy, cover(), and ",
         "occu_cover() fits.", call. = FALSE)
  }
  draws <- object$draws; pi_list <- model$process_info
  X_occ <- model$X_processes[[1]]; X_det <- model$X_processes[[2]]
  y <- model$y; n_sites <- model$n_sites
  p_occ <- pi_list[[1]]$p; p_det <- pi_list[[2]]$p
  n_draws <- min(n.samples, nrow(draws))
  draw_idx <- sample.int(nrow(draws), n_draws)

  pit <- numeric(n_sites)
  for (i in seq_len(n_sites)) {
    yi <- y[i, ]; valid <- yi >= 0; n_det <- sum(yi[valid] == 1); n_valid <- sum(valid)
    if (n_valid == 0) { pit[i] <- runif(1); next }
    cdf_vals <- numeric(n_draws)
    for (s in seq_len(n_draws)) {
      psi <- plogis(sum(X_occ[i, ] * draws[draw_idx[s], seq_len(p_occ)]))
      p <- plogis(sum(X_det[i, ] * draws[draw_idx[s], p_occ + seq_len(p_det)]))
      cdf_vals[s] <- if (n_det > 0) 1 else psi * (1 - p)^n_valid + (1 - psi)
    }
    pit[i] <- min(1, max(0, mean(cdf_vals) + runif(1, 0, 1/n_draws)))
  }
  pit
}

#' Goodness-of-fit tests for a tobs fit
#'
#' Posterior-predictive goodness-of-fit checks for an occupancy / N-mixture
#' `tobs_fit`: dispersion, zero-inflation, and outlier counts compared against
#' the fitted model's simulated replicates, plus a Kolmogorov-Smirnov
#' uniformity test on PIT residuals.
#'
#' @param pit Numeric vector of PIT residuals (`tobs_test_uniformity`).
#' @param object A fitted `tobs_fit`.
#' @param n.samples Number of posterior-predictive replicates to simulate.
#' @return A list with the observed statistic, its posterior-predictive
#'   expectation, and a tail p-value; `tobs_test_uniformity` returns the
#'   `htest` object from [stats::ks.test()].
#' @name tobs_gof_tests
NULL

#' @rdname tobs_gof_tests
#' @export
tobs_test_uniformity <- function(pit) {
  # PIT residuals can have ties when the response is discrete (counts /
  # detections). The asymptotic KS p-value remains valid; only the
  # ties-warning is noisy, so suppress just that warning class.
  withCallingHandlers(
    ks.test(pit, "punif"),
    warning = function(w) {
      if (grepl("ties should not be present", conditionMessage(w),
                fixed = TRUE)) invokeRestart("muffleWarning")
    }
  )
}

#' @rdname tobs_gof_tests
#' @export
tobs_test_dispersion <- function(object, n.samples = 250) {
  sims <- simulate(object, nsim = n.samples); y_obs <- object$model$y
  obs_var <- var(rowSums(y_obs * (y_obs >= 0), na.rm = TRUE))
  sim_vars <- vapply(sims, function(ys) var(rowSums(ys * (ys >= 0), na.rm = TRUE)), double(1))
  list(observed = obs_var, expected = mean(sim_vars),
       ratio = obs_var / mean(sim_vars), p.value = mean(sim_vars >= obs_var))
}

#' @rdname tobs_gof_tests
#' @export
tobs_test_zero_inflation <- function(object, n.samples = 250) {
  sims <- simulate(object, nsim = n.samples); y_obs <- object$model$y
  count_zeros <- function(y) sum(apply(y, 1, function(r) { v <- r >= 0; all(r[v] == 0) }))
  obs <- count_zeros(y_obs); sim <- vapply(sims, count_zeros, integer(1))
  list(observed = obs, expected = mean(sim), ratio = obs / max(mean(sim), 1),
       p.value = mean(sim >= obs))
}

#' @rdname tobs_gof_tests
#' @export
tobs_test_outliers <- function(object, n.samples = 250) {
  sims <- simulate(object, nsim = n.samples); y_obs <- object$model$y
  n_sites <- nrow(y_obs)
  obs_det <- apply(y_obs, 1, function(r) sum(r[r >= 0] == 1))
  sim_det <- vapply(sims, function(ys) apply(ys, 1, function(r) sum(r[r >= 0] == 1)), double(n_sites))
  lower <- apply(sim_det, 1, quantile, 0.025); upper <- apply(sim_det, 1, quantile, 0.975)
  n_out <- sum(obs_det < lower | obs_det > upper)
  sim_out <- vapply(seq_len(n.samples), function(s) sum(sim_det[,s] < lower | sim_det[,s] > upper), integer(1))
  list(n_outliers = n_out, expected = mean(sim_out), p.value = mean(sim_out >= n_out))
}

#' Comprehensive model checking
#' @param object A `tobs_fit` object.
#' @param coords Optional coordinates for spatial diagnostics.
#' @param n.samples Posterior samples for simulation tests.
#' @return Invisibly, diagnostic results.
#' @export
tobs_check <- function(object, coords = NULL, n.samples = 250) {
  cat("=== tobs Model Diagnostics ===\n\n")
  if (identical(object$method, "nuts")) {
    cat(sprintf("Sampler: %d samples, %d divergent, mean accept = %.3f\n",
                object$n_samples, sum(object$divergent), mean(object$accept_prob)))
  } else {
    cat(sprintf("Fit: %s, %d posterior draws (no NUTS sampler diagnostics)\n",
                object$method %||% "laplace", object$n_samples))
  }

  w <- tryCatch(tobs_waic(object), error = function(e) NULL)
  if (!is.null(w)) cat(sprintf("\nWAIC: %.1f (p_waic = %.1f)\n", w$waic, w$p_waic))

  dic <- tryCatch(tobs_dic(object), error = function(e) NULL)
  if (!is.null(dic) && is.finite(dic$dic))
    cat(sprintf("DIC:  %.1f (p_DIC = %.1f)\n", dic$dic, dic$p_dic))

  cpo <- tryCatch(tobs_cpo(object), error = function(e) NULL)
  if (!is.null(cpo)) cat(sprintf("LPML: %.1f (elpd_loo = %.1f)\n",
                                 cpo$lpml, cpo$elpd_loo))

  ppc <- tryCatch(tobs_ppc(object, n.samples = n.samples), error = function(e) NULL)
  if (!is.null(ppc)) {
    cat(sprintf("\nPPC: Bayesian p = %.3f\n", ppc$bayesian.p))
    if (ppc$bayesian.p < 0.05 || ppc$bayesian.p > 0.95) cat("  WARNING: poor fit\n")
  }

  zi <- tryCatch(tobs_test_zero_inflation(object, n.samples = n.samples), error = function(e) NULL)
  if (!is.null(zi)) cat(sprintf("\nZero-inflation: obs=%d, exp=%.1f, p=%.3f\n", zi$observed, zi$expected, zi$p.value))

  if (!is.null(coords)) {
    mi <- tryCatch(tulpa::moran_i(residuals(object)$occ, coords), error = function(e) NULL)
    if (!is.null(mi)) {
      cat(sprintf("\nMoran's I: %.3f (p = %.3f)\n", unname(mi$statistic), mi$p.value))
      if (mi$p.value < 0.05) cat("  WARNING: spatial autocorrelation\n")
    }
  }
  cat("\n")
  invisible(list(waic = w, dic = dic, cpo = cpo, ppc = ppc,
                 zero_inflation = zi))
}
