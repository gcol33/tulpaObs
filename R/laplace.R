# ============================================================================
# Occupancy-specific callbacks for tulpa's EM+Laplace engine
# ============================================================================

#' Fit a tobs model with Laplace approximation (internal)
#'
#' Uses tulpa's generic EM+Laplace engine ([tulpa::tulpa_em_laplace()]) with
#' occupancy-specific E-step and M-step encoding, augmented with a
#' weakly-informative quadratic prior on the fixed-effect coefficients (see
#' [occu_priors()]). The prior is attached to each M-step block as a per-block
#' `beta_prior` (see `.attach_priors_to_blocks()` in `R/occu_priors.R`), which
#' tulpa threads through every phase -- the EM iterations and the MI / Gibbs
#' correction refits alike. Supports all built-in model
#' types. Called from `.tobs_fit_model()`; not user-facing.
#'
#' @param model A `tobs_model` from `.tobs_build_model()`.
#' @param spatial Optional `tobs_spatial` spec (NULL for non-spatial).
#' @param priors Optional prior spec from [occu_priors()]. `NULL` -> use
#'   the package defaults. Pass `FALSE` (or `"none"`) to disable the prior
#'   and recover the historical unpenalised MAP behavior — the penalised
#'   objective is `Q(beta) = -log L(beta) + sum (beta_j - mu_j)^2 / (2 sd_j^2)`,
#'   so `sd_j = Inf` yields a zero penalty term.
#' @param sigma_beta Reserved for future use (NUTS-side beta prior); ignored
#'   by the EM-Laplace path.
#' @param max_iter,tol,damping EM controls.
#' @param correction Post-EM correction (`"none"`, `"mi"`, `"gibbs"`). MI /
#'   Gibbs run tulpa's post-EM Rubin-pooled correction; the fixed-effect prior
#'   (when active) threads into the correction refits, so the corrected fit is
#'   penalised the same way as the EM point estimate.
#' @param n_imputations Number of MI draws when `correction = "mi"`.
#' @param verbose Print per-iteration progress.
#' @keywords internal
.tobs_laplace <- function(model, spatial = NULL, re = NULL,
                          priors = NULL,
                          sigma_beta = 10,
                          max_iter = 50L, tol = 1e-4, damping = 0.3,
                          correction = c("auto", "mi", "gibbs", "none"),
                          n_imputations = 20L, n_gibbs = 10L, seed = NULL,
                          approx = c("gaussian_laplace", "simplified_laplace"),
                          re_aghq = TRUE, n_quad = 9L, lkj_eta = 1.5,
                          latent_prior = NULL, heldout_state = NULL,
                          verbose = TRUE) {
  correction <- match.arg(correction)
  approx <- match.arg(approx)
  if (!inherits(model, "tobs_model")) stop("model must be a tobs_model object")

  # Nested-Laplace path: a multi-block latent prior is attached to the state
  # ("occ") M-step block so tulpa's EM dispatches it through
  # tulpa_nested_laplace() with hyperparameter integration (see
  # `.tobs_em_nested_laplace()`). The latent prior is mutually exclusive with an
  # SPDE `spatial` term (that is the single-Laplace continuous-spatial route)
  # and with the variance-component RE EM (RE blocks are folded into the multi-
  # block prior upstream), and it has no MI/Gibbs correction.
  if (!is.null(latent_prior)) {
    if (!is.null(spatial)) {
      stop("latent_prior and an SPDE spatial term cannot be combined on the ",
           "Laplace path.", call. = FALSE)
    }
    if (!is.null(re)) {
      stop("latent_prior is built from the spatial/temporal/re terms upstream; ",
           "pass re = NULL when latent_prior is set.", call. = FALSE)
    }
    if (!correction %in% c("auto", "none")) {
      stop("Nested Laplace (latent_prior) has no MI/Gibbs correction.",
           call. = FALSE)
    }
    return(.tobs_laplace_nested(model, latent_prior = latent_prior,
                                heldout_state = heldout_state,
                                priors = priors, max_iter = max_iter,
                                tol = tol, damping = damping, approx = approx,
                                verbose = verbose))
  }

  .validate_spatial_laplace(spatial, model$model_type)

  # Formula random effects on the deterministic path. Supported forms (iid
  # intercept, uncorrelated slopes, and correlated slopes on the occupancy
  # predictor of a single-season model) are fit via the variance-component EM
  # in R/em_laplace_re.R; everything else errors with a pointer to NUTS rather
  # than being silently dropped (gcol33/tulpaObs#11).
  if (!is.null(re)) {
    .validate_re_laplace(re, model, spatial, approx)
    em_result <- .tobs_em_laplace_re(model, re, priors = priors,
                                     max_iter = max_iter, tol = tol,
                                     damping = damping, aghq = re_aghq,
                                     n_quad = n_quad, lkj_eta = lkj_eta,
                                     verbose = verbose)
    re_block <- .tobs_re_param_block(em_result$re_post)
    fit <- build_laplace_fit(em_result, model, spatial,
                             c(occ = ncol(model$X_processes[[1]]),
                               det = ncol(model$X_processes[[2]])),
                             prior_spec = NULL, approx = "gaussian_laplace",
                             re_block = re_block)
    fit$re <- if (inherits(re, "tobs_re")) list(re) else re
    return(fit)
  }

  callbacks <- switch(model$model_type,
    single     = build_single_callbacks(model, spatial),
    dynamic    = build_dynamic_callbacks(model, spatial),
    integrated = build_integrated_callbacks(model, spatial),
    jsdm       = build_jsdm_callbacks(model, spatial),
    stop(sprintf("Laplace not supported for model_type '%s'", model$model_type))
  )

  prior_spec <- .resolve_occu_priors(priors)

  # Single engine for every Laplace fit: tulpa's generic EM+Laplace. The
  # fixed-effect prior is attached per M-step block as a `beta_prior` (see
  # .attach_priors_to_blocks); tulpa's block fitter applies it in every phase,
  # so a prior-aware MI/Gibbs correction comes for free (gcol33/tulpa#27).
  # Spatial fits are left unpenalised here -- the SPDE/NNGP solver carries its
  # own fixed-effect prior and tulpa_laplace() rejects `beta_prior` on the
  # spatial path (tulpaObs#5) -- so the prior is attached only when there is no
  # spatial term.
  m_step_encode <- if (is.null(spatial)) {
    function(weights, ...) {
      .attach_priors_to_blocks(callbacks$m_step_encode(weights, ...),
                               model, prior_spec)
    }
  } else {
    callbacks$m_step_encode
  }

  # MI / Gibbs draw hard z with R's RNG; seed it so the corrected fit
  # reproduces.
  if (correction %in% c("mi", "gibbs") && !is.null(seed)) {
    set.seed(as.integer(seed))
  }
  em_result <- tulpa::tulpa_em_laplace(
    e_step        = callbacks$e_step,
    m_step_encode = m_step_encode,
    draw_z        = callbacks$z_draw,
    max_iter      = max_iter,
    tol           = tol,
    damping       = damping,
    correction    = correction,
    n_imputations = n_imputations,
    n_gibbs       = n_gibbs,
    verbose       = verbose
  )

  # tulpa_em_laplace returns flat convergence fields; synthesize the nested
  # `convergence` list that build_laplace_fit() / summary() expect.
  if (is.null(em_result$convergence)) {
    em_result$convergence <- list(converged = em_result$converged,
                                  n_iter = em_result$n_iter,
                                  history = em_result$history)
  }

  fit <- build_laplace_fit(em_result, model, spatial, callbacks$p_per_submodel,
                           prior_spec = prior_spec,
                           approx = approx)
  fit$priors <- prior_spec

  # Debias step (gcol33/tulpaObs#7): for single-season occupancy the marginal
  # likelihood is closed-form, so refine the EM mode with an exact-marginal
  # Newton step and read calibrated SEs from its Hessian. The EM's M-inflated
  # pseudo-binomial Laplace attenuates the detection coefficients and
  # under-disperses their SEs at small J; the refinement restores unbiased,
  # near-nominal-coverage fixed effects. Spatial / multi-season fits have no
  # closed-form marginal and keep the EM result.
  if (identical(model$model_type, "single") && is.null(spatial) &&
      identical(approx, "gaussian_laplace") &&
      correction %in% c("auto", "none")) {
    fit <- .tobs_occu_marginal_refine(fit, model, prior_spec)
  }

  # Record the seed used for a stochastic correction so the run reproduces.
  if (correction %in% c("mi", "gibbs") && !is.null(seed)) {
    fit$seed <- as.integer(seed)
  }
  fit
}

# ============================================================================
# Nested-Laplace driver
#
# Shares the per-model-type E-step / M-step callbacks with `.tobs_laplace`
# (single source of truth: there is no `build_*_callbacks_nested` family). The
# only difference from the single-Laplace path is that the multi-block latent
# `prior` is attached to the state ("occ") M-step block, which makes tulpa's
# generic EM engine dispatch that block through `tulpa::tulpa_nested_laplace()`
# (hyperparameter integration) instead of `tulpa::tulpa_laplace()`. The
# detection block(s) stay single-Laplace (nested-Laplace does not accept the
# per-row observation weights they carry).
#
# Held-out state units (`heldout_state`, integer indices into the state block's
# rows) are encoded with `n_trials = 0`: they drop out of the likelihood but
# stay in the design so their latent value is pulled from the prior (spatial
# neighbours / shared field) and `X beta` -- the INLA NA-response mechanism that
# `.tobs_predict_heldout()` reads back marginalised over the grid.
# ============================================================================
.tobs_laplace_nested <- function(model, latent_prior, heldout_state = NULL,
                                 priors = NULL, max_iter = 25L, tol = 1e-3,
                                 damping = 0.3, approx = "gaussian_laplace",
                                 verbose = TRUE) {
  if (!is.null(priors) && !isFALSE(priors)) {
    message(".tobs_laplace_nested(): fixed-effect priors are not applied on ",
            "the nested-Laplace path; the latent block carries its own prior.")
  }

  callbacks <- switch(model$model_type,
    single     = build_single_callbacks(model, spatial = NULL,
                                        latent_prior = latent_prior),
    dynamic    = build_dynamic_callbacks(model, spatial = NULL),
    integrated = build_integrated_callbacks(model, spatial = NULL),
    jsdm       = build_jsdm_callbacks(model, spatial = NULL),
    stop(sprintf("nested Laplace not supported for model_type '%s'",
                 model$model_type), call. = FALSE)
  )

  m_step_encode <- function(weights, ...) {
    blocks <- callbacks$m_step_encode(weights, ...)
    blocks$occ$prior <- latent_prior
    if (!is.null(heldout_state) && length(heldout_state) > 0L) {
      blocks$occ$n_trials[heldout_state] <- 0L
      blocks$occ$y[heldout_state]        <- 0L
    }
    blocks
  }

  em_result <- tulpa::tulpa_em_laplace(
    e_step        = callbacks$e_step,
    m_step_encode = m_step_encode,
    max_iter      = max_iter,
    tol           = tol,
    damping       = damping,
    correction    = "none",
    verbose       = verbose
  )
  if (is.null(em_result$convergence)) {
    em_result$convergence <- list(converged = em_result$converged,
                                  n_iter = em_result$n_iter,
                                  history = em_result$history)
  }

  fit <- build_laplace_fit(em_result, model, spatial = NULL,
                           callbacks$p_per_submodel,
                           prior_spec = NULL, approx = approx,
                           latent_prior = latent_prior)
  fit$method <- "nested_laplace"

  # State-field posterior. For single-season occupancy, refine the EM field with
  # one exact-marginal occupancy-family pass (calibrated mode, fitted_eta_var and
  # grid weights with no M-inflation; see .tobs_occu_state_marginal_fit). Other
  # model types keep the EM occ fit (their NA-response mapping is not yet wired,
  # and the occupancy-family reduction is single-season specific).
  state_fit <- if (identical(model$model_type, "single")) {
    .tobs_occu_state_marginal_fit(model, em_result, latent_prior,
                                  max_iter = max_iter, tol = tol)
  } else {
    em_result$fits$occ
  }
  fit$nested_laplace <- list(multi_prior = latent_prior,
                             occ_fit     = state_fit,
                             heldout     = heldout_state)
  # Marginalised state-level psi posterior (per state row, integrated over the
  # hyperparameter grid). Computed here on the fitting-scale design so the
  # per-cell betas and the field share one space; psi is scale-invariant and
  # survives the per-process unscaling downstream. Held-out rows are the
  # NA-response prediction targets.
  fit$state_posterior <- .tobs_nested_state_posterior(
    model, state_fit, latent_prior, heldout = heldout_state)

  # Areal field summary: each areal block's field SD (sigma) integrated over the
  # outer grid, and the per-cell field(s) read off the calibrated state fit. A
  # single intercept field reports `sigma` + `spatial_field`; a varying-
  # coefficient structure (intercept + trend field(s)) adds `sigma_trend` and
  # `trend_field` so the occu() output mirrors what occu_cover() exposes.
  fit <- .tobs_nested_attach_field_summary(fit, model, state_fit, latent_prior)
  fit
}


# Surface the areal hyperparameters + per-cell field(s) of a nested-Laplace occu
# fit. `latent_prior` is the list of latent blocks; the areal (icar / car_proper)
# blocks each carry one field, in formula order (intercept field first, then any
# varying-coefficient trend fields). `state_fit` is the calibrated occupancy-
# marginal nested fit (theta_grid, weights, modes). Each areal block's SD is the
# derived sigma = 1/sqrt(tau) marginalized over the outer grid (the marginalize-
# derived-quantities rule: a weighted summary of sigma per cell, not a plug-in of
# the tau posterior mean); each field is the grid-weighted posterior mean of the
# block's latent slice, demeaned to the sum-to-zero convention the ICAR prior
# sits under. Temporal / iid blocks are skipped (they surface elsewhere).
.tobs_nested_attach_field_summary <- function(fit, model, state_fit,
                                              latent_prior) {
  blocks <- if (!is.null(latent_prior$type)) list(latent_prior) else latent_prior
  is_areal <- vapply(blocks, function(b)
    isTRUE(b$type %in% c("icar", "car_proper")), logical(1))
  if (!any(is_areal)) return(fit)
  if (is.null(state_fit$theta_grid) || is.null(state_fit$weights) ||
      is.null(state_fit$modes)) {
    return(fit)
  }

  w  <- state_fit$weights
  w  <- w / sum(w)
  tg <- state_fit$theta_grid
  if (!is.matrix(tg)) tg <- matrix(tg, ncol = 1L,
                                   dimnames = list(NULL, state_fit$theta_names))
  tg_names <- colnames(tg)
  modes <- state_fit$modes
  p_occ <- ncol(model$X_processes[[1]])

  # 0-based latent offset of each block's field slice in the mode tail. The
  # areal blocks contribute n_spatial_units each; bym2 contributes 2x (phi +
  # theta) but only the areal SVC path (icar / car_proper, d_fac = 1) is wired
  # to a per-field summary here.
  blk_len <- vapply(blocks, .nl_block_field_len, integer(1))
  blk_off <- cumsum(c(0L, blk_len))[seq_along(blocks)]

  sigma_means <- numeric(0); sigma_sds <- numeric(0); sigma_nms <- character(0)
  field_list  <- list()
  areal_idx   <- which(is_areal)
  for (j in seq_along(areal_idx)) {
    b <- areal_idx[j]
    n_units <- as.integer(blocks[[b]]$n_spatial_units %||% blk_len[b])

    # sigma = 1/sqrt(tau) marginalized over the grid. The ICAR axis is the
    # precision tau (joint-grid column `b<b>.tau`); fall back to the single-axis
    # grid when there is one block.
    tau_col <- match(sprintf("b%d.tau", b), tg_names)
    if (is.na(tau_col)) tau_col <- match("tau", tg_names)
    if (!is.na(tau_col)) {
      tau_vals <- as.numeric(tg[, tau_col])
      sig_vals <- 1 / sqrt(pmax(tau_vals, 1e-12))
      s_mean <- sum(w * sig_vals)
      s_sd   <- sqrt(max(sum(w * sig_vals^2) - s_mean^2, 0))
      nm <- if (j == 1L) "sigma"
            else if (length(areal_idx) == 2L) "sigma_trend"
            else sprintf("sigma_trend%d", j - 1L)
      sigma_means[[nm]] <- s_mean
      sigma_sds  [[nm]] <- s_sd
      sigma_nms <- c(sigma_nms, nm)
    }

    # Grid-weighted posterior mean of this block's field slice, demeaned.
    cols <- p_occ + blk_off[b] + seq_len(n_units)
    if (max(cols) <= ncol(modes)) {
      fld <- as.numeric(crossprod(modes[, cols, drop = FALSE], w))
      fld <- fld - mean(fld)
      field_list[[j]] <- fld
    } else {
      field_list[[j]] <- rep(NA_real_, n_units)
    }
  }

  if (length(sigma_nms)) {
    keep <- sigma_nms %in% names(sigma_means)
    fit$means <- c(fit$means, unlist(sigma_means)[sigma_nms[keep]])
    fit$sds   <- c(fit$sds,   unlist(sigma_sds)[sigma_nms[keep]])
  }

  field_z_table <- function(fld) {
    n <- length(fld)
    data.frame(cell = seq_len(n), z_mean = fld,
               z_sd = rep(NA_real_, n),
               z_lower = rep(NA_real_, n), z_upper = rep(NA_real_, n))
  }
  if (length(field_list) >= 1L) {
    fit$spatial_field <- field_list[[1L]]
    fit$field_table   <- field_z_table(field_list[[1L]])
  }
  if (length(field_list) >= 2L) {
    trend_means <- field_list[-1L]
    fit$trend_field   <- trend_means[[1L]]
    fit$trend_fields  <- trend_means
    fit$trend_field_table  <- field_z_table(trend_means[[1L]])
    fit$trend_field_tables <- lapply(trend_means, field_z_table)
  }
  fit
}

# ============================================================================
# Single-season callbacks
# ============================================================================
build_single_callbacks <- function(model, spatial = NULL, latent_prior = NULL) {
  y <- model$y
  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  X_det_visit <- model$X_det_visit  # NULL when no visit-level covariates
  max_visits <- ncol(y)
  n_sites <- model$n_sites
  p_occ <- ncol(X_occ)
  p_det <- ncol(X_det)
  p_det_visit <- if (is.null(X_det_visit)) 0L else ncol(X_det_visit)
  p_det_total <- p_det + p_det_visit

  # An SPDE term may enter the state arm, the detection arm, or both (its
  # `$shared = c(occ, det)` membership). Resolve the per-arm field once: the
  # state field attaches to the occ block, the detection field to the det
  # block. A detection field with visit-level detection covariates is not yet
  # wired (the field is site-indexed; the det block is per (site, visit) and
  # would need a row-expanded mesh projection).
  spatial_occ <- .spatial_for_arm(spatial, 1L)
  spatial_det <- .spatial_for_arm(spatial, 2L)
  if (!is.null(spatial_det) && p_det_visit > 0L) {
    stop("SPDE on the detection process with visit-level detection covariates ",
         "is not yet plumbed in .tobs_laplace; use shared occupancy-arm SPDE ",
         "or method = 'nuts'.", call. = FALSE)
  }

  # A continuous Matern (SPDE) block on the nested-Laplace latent prior needs
  # the same modest pseudo-binomial inflation the single-Laplace SPDE path uses
  # (M = 4): at M = 1000 the data signal swamps the SPDE prior precision, the
  # mesh field over-fits, and the occupancy slope inflates (the field absorbs
  # the covariate signal). The areal icar/bym2/car_proper blocks are strongly
  # informative at the grid scale and tolerate the sharp M = 1000 encoding, so
  # the modest M is gated on a continuous-field block only.
  nested_has_spde <- .tobs_latent_prior_has_spde(latent_prior)

  n_valid <- integer(n_sites)
  n_det <- integer(n_sites)
  any_det <- logical(n_sites)
  valid_mat <- matrix(FALSE, n_sites, max_visits)
  for (i in seq_len(n_sites)) {
    v <- y[i, ] >= 0
    valid_mat[i, ] <- v
    n_valid[i] <- sum(v)
    n_det[i] <- sum(y[i, v] == 1)
    any_det[i] <- n_det[i] > 0
  }
  keep <- n_valid > 0

  # Per-visit indexing for the X_det_visit path. Row r of X_det_visit
  # corresponds to (site = (r-1) %/% max_visits + 1, visit = (r-1) %% max_visits + 1).
  if (p_det_visit > 0L) {
    site_idx_all <- rep(seq_len(n_sites), each = max_visits)
    visit_idx_all <- rep(seq_len(max_visits), times = n_sites)
  }

  e_step <- function(fits, ...) {
    beta_occ <- extract_beta(fits$occ, p_occ)
    eta_occ <- as.vector(X_occ %*% beta_occ)
    sp_off <- .spatial_eta_offset(spatial_occ, fits$occ, p_occ)
    if (length(sp_off) == n_sites) eta_occ <- eta_occ + sp_off
    # Nested-Laplace: make the E-step weight P(z_i = 1 | y_i) field-aware (the
    # field informs which undetected sites are occupied; without this the EM
    # converges to the fixed-effect-only fixed point and the field cannot track
    # the data). Prefer the engine's exact per-cell fitted eta, marginalised
    # over the hyperparameter grid -- correct for every prior including bym2.
    # Fall back to the grid-weighted mode field offset (exact for d_fac = 1
    # priors; skips bym2) when the engine did not return fitted_eta.
    if (!is.null(latent_prior)) {
      eta_marg <- .nested_eta_marginal(fits$occ, n_sites)
      if (!is.null(eta_marg)) {
        eta_occ <- eta_marg
      } else {
        lat_off <- .nested_eta_offset(latent_prior, fits$occ, p_occ, n_sites)
        if (length(lat_off) == n_sites) eta_occ <- eta_occ + lat_off
      }
    }
    psi <- plogis(eta_occ)

    if (p_det_visit == 0L) {
      beta_det <- extract_beta(fits$det, p_det)
      eta_det <- as.vector(X_det %*% beta_det)
      det_off <- .spatial_eta_offset(spatial_det, fits$det, p_det)
      if (length(det_off) == n_sites) eta_det <- eta_det + det_off
      p_site <- plogis(eta_det)
      return(list(weights = occ_weights(psi, p_site, n_sites,
                                        n_valid, n_det, any_det)))
    }

    # Visit-level path: logit(p_ij) = X_det[i,] beta_site + X_det_visit[(i-1)*J + j,] beta_visit
    beta_det <- extract_beta(fits$det, p_det_total)
    eta_site <- as.vector(X_det %*% beta_det[seq_len(p_det)])
    eta_visit_long <- as.vector(X_det_visit %*%
                                  beta_det[(p_det + 1L):p_det_total])
    # eta_visit_long is in site-major order: reshape so [i, j] = visit (i, j)
    eta_visit_mat <- matrix(eta_visit_long, n_sites, max_visits, byrow = TRUE)
    logit_p_ij <- matrix(eta_site, n_sites, max_visits) + eta_visit_mat
    logit_p_ij <- pmin(pmax(logit_p_ij, -30), 30)
    # log(1 - plogis(eta)) = -log1pexp(eta) computed stably as -pmax(eta,0) - log1p(exp(-|eta|))
    log_1mp <- -(pmax(logit_p_ij, 0) + log1p(exp(-abs(logit_p_ij))))
    log_1mp[!valid_mat] <- 0
    log_prod_1mp <- rowSums(log_1mp)
    weights <- numeric(n_sites)
    for (i in seq_len(n_sites)) {
      if (any_det[i]) {
        weights[i] <- 1
      } else if (n_valid[i] == 0L) {
        weights[i] <- psi[i]
      } else {
        num <- psi[i] * exp(log_prod_1mp[i])
        weights[i] <- num / (num + (1 - psi[i]))
      }
    }
    list(weights = weights)
  }

  m_step_encode <- function(weights, ...) {
    if (is.null(spatial_occ) && !nested_has_spde) {
      # Pseudo-binomial encoding: y = round(M*w), n_trials = M. The
      # M-inflation makes the M-step into a sharp binomial whose mode
      # equals the weighted mean, and is the historical encoding used
      # everywhere else in the package (areal nested-Laplace blocks too --
      # their intrinsic-field prior is strong at the grid scale).
      M <- 1000L
      y_occ <- ifelse(any_det, M, as.integer(round(weights * M)))
      y_occ <- pmin(pmax(y_occ, 0L), M)
      occ_block <- list(y = y_occ, n_trials = rep(M, n_sites), X = X_occ,
                        family = "binomial")
    } else if (is.null(spatial_occ) && nested_has_spde) {
      # Nested-Laplace continuous SPDE block: modest M (= 4) so the mesh-field
      # prior is not swamped, mirroring the single-Laplace SPDE encoding. The
      # block prior is attached to occ$prior upstream (.tobs_laplace_nested),
      # so no .attach_spatial_spde() here.
      M <- 4L
      y_occ <- ifelse(any_det, M, as.integer(round(weights * M)))
      y_occ <- pmin(pmax(y_occ, 0L), M)
      occ_block <- list(y = y_occ, n_trials = rep(M, n_sites), X = X_occ,
                        family = "binomial")
    } else {
      # Modest pseudo-binomial encoding for SPDE: M = 4 gives some
      # fractional resolution on the weights while keeping the per-site
      # effective sample size O(1), so the SPDE prior precision is not
      # swamped by the data signal as it would be at M = 1000.
      M <- 4L
      y_occ <- ifelse(any_det, M, as.integer(round(weights * M)))
      y_occ <- pmin(pmax(y_occ, 0L), M)
      occ_block <- list(y = y_occ, n_trials = rep(M, n_sites), X = X_occ,
                        family = "binomial")
      occ_block <- .attach_spatial_spde(occ_block, spatial_occ)
    }
    # Detection block: weight by w_i = P(z_i = 1 | y_i, theta). Sites that
    # the E-step thinks are likely empty (w_i ~ 0) must drop out of the
    # detection fit, otherwise they bias p_hat downward by feeding their
    # all-zero detection history as evidence about (1 - p)^J. Sites with
    # any detection have w_i = 1 (the E-step sets this).
    w_det <- weights
    w_det[any_det] <- 1

    if (!is.null(spatial_det)) {
      # SPDE detection field: the single-Laplace spatial solver consumes no
      # per-observation `weights`, so the occupancy weight is folded into the
      # binomial response by scaling both successes and trials by w_i
      # (y = round(w_i n_det_i), n = round(w_i n_valid_i)). This is the
      # frequency-weight-as-counts identity for a binomial mode/Hessian, and it
      # keeps ALL n_sites rows so the per-site rows stay aligned with the full
      # mesh projection A (n_sites x n_mesh). A near-empty site (w_i ~ 0)
      # collapses to a (0, 0) row that contributes nothing to the likelihood,
      # score, or Hessian -- the analogue of dropping it under the explicit
      # weight on the non-spatial path.
      y_det_w <- as.integer(round(w_det * n_det))
      n_det_w <- as.integer(round(w_det * n_valid))
      y_det_w <- pmin(pmax(y_det_w, 0L), n_det_w)
      det_block <- list(y = y_det_w, n_trials = n_det_w, X = X_det,
                        family = "binomial")
      det_block <- .attach_spatial_spde(det_block, spatial_det)
    } else if (p_det_visit == 0L) {
      keep_det <- keep & (w_det > 1e-6)
      det_block <- list(y = n_det[keep_det], n_trials = n_valid[keep_det],
                        X = X_det[keep_det, , drop = FALSE],
                        weights = w_det[keep_det], family = "binomial")
    } else {
      # Per-visit detection block: one Bernoulli row per (site, valid visit)
      # whose weight is the site's posterior occupancy w_i. Combined design
      # matrix stacks site-level X_det (replicated across visits) and
      # visit-level X_det_visit (already in site-major order).
      keep_visit <- valid_mat & (w_det >= 1e-6)
      site_kept <- site_idx_all[as.vector(t(keep_visit))]
      visit_kept <- visit_idx_all[as.vector(t(keep_visit))]
      n_kept <- length(site_kept)
      if (n_kept > 0L) {
        visit_row_idx <- (site_kept - 1L) * max_visits + visit_kept
        X_combined <- cbind(
          X_det[site_kept, , drop = FALSE],
          X_det_visit[visit_row_idx, , drop = FALSE]
        )
        y_kept <- y[cbind(site_kept, visit_kept)]
        det_block <- list(y = as.integer(y_kept),
                          n_trials = rep(1L, n_kept),
                          X = X_combined,
                          weights = w_det[site_kept],
                          family = "binomial")
      } else {
        det_block <- list(y = integer(0),
                          n_trials = integer(0),
                          X = matrix(0, 0, p_det_total),
                          weights = numeric(0),
                          family = "binomial")
      }
    }
    list(occ = occ_block, det = det_block)
  }

  z_draw <- function(weights, ...) {
    z <- as.integer(any_det)
    z[!any_det] <- rbinom(sum(!any_det), 1, clamp_w(weights[!any_det]))
    z
  }

  hard_encode <- function(z, ...) {
    occ_sites <- which(z == 1L)
    det_keep <- occ_sites[n_valid[occ_sites] > 0]
    occ_block <- list(y = z, n_trials = rep(1L, n_sites), X = X_occ,
                      family = "binomial")
    occ_block <- .attach_spatial_spde(occ_block, spatial_occ)
    if (!is.null(spatial_det)) {
      # Hard-z detection field: keep ALL n_sites rows aligned with the mesh
      # projection by zeroing the trials of sites that contribute no detection
      # evidence (z = 0 or no valid visits); a (0, 0) row drops out cleanly.
      keep_det_mask <- (z == 1L) & (n_valid > 0L)
      n_det_h <- ifelse(keep_det_mask, n_valid, 0L)
      y_det_h <- ifelse(keep_det_mask, n_det, 0L)
      det_block <- list(y = as.integer(y_det_h), n_trials = as.integer(n_det_h),
                        X = X_det, family = "binomial")
      det_block <- .attach_spatial_spde(det_block, spatial_det)
    } else if (p_det_visit == 0L) {
      det_block <- if (length(det_keep) > 0)
        list(y = n_det[det_keep], n_trials = n_valid[det_keep],
             X = X_det[det_keep, , drop = FALSE], family = "binomial")
      else NULL
    } else {
      keep_visit <- valid_mat & matrix(z == 1L, n_sites, max_visits)
      site_kept <- site_idx_all[as.vector(t(keep_visit))]
      visit_kept <- visit_idx_all[as.vector(t(keep_visit))]
      n_kept <- length(site_kept)
      det_block <- if (n_kept > 0L) {
        visit_row_idx <- (site_kept - 1L) * max_visits + visit_kept
        X_combined <- cbind(
          X_det[site_kept, , drop = FALSE],
          X_det_visit[visit_row_idx, , drop = FALSE]
        )
        y_kept <- y[cbind(site_kept, visit_kept)]
        list(y = as.integer(y_kept), n_trials = rep(1L, n_kept),
             X = X_combined, family = "binomial")
      } else NULL
    }
    list(occ = occ_block, det = det_block)
  }

  init <- glm_init(X_occ, X_det, any_det, n_det, n_valid, keep, p_occ, p_det)
  if (p_det_visit > 0L) {
    # Pad det init with zeros for visit-level cols so warm-start shapes
    # match the combined design when the penalized driver pulls beta_init
    # from the previous EM iteration.
    init$det$beta <- c(init$det$beta, rep(0, p_det_visit))
    init$det$se   <- c(init$det$se,   rep(1, p_det_visit))
  }

  list(e_step = e_step, m_step_encode = m_step_encode, z_draw = z_draw,
       hard_encode = hard_encode, init = init,
       p_per_submodel = c(occ = p_occ, det = p_det_total))
}

# ============================================================================
# Dynamic occupancy callbacks
# ============================================================================
build_dynamic_callbacks <- function(model, spatial = NULL) {
  y_flat <- model$y_flat
  n_sites <- model$n_sites
  n_seasons <- model$n_seasons
  max_visits <- model$max_visits
  X_occ <- model$X_processes[[1]]  # psi1
  X_det <- model$X_processes[[2]]  # p
  X_col <- model$X_processes[[3]]  # gamma
  X_ext <- model$X_processes[[4]]  # epsilon
  p_occ <- ncol(X_occ); p_det <- ncol(X_det)
  p_col <- ncol(X_col); p_ext <- ncol(X_ext)

  # The state field enters season-1 occupancy psi1 only (one psi1 row per
  # site, the identity map). The colonization (gamma) and extinction
  # (epsilon) transition predictors are separate latent processes whose own
  # mesh fields are not wired here; `.validate_spatial_laplace` only routes a
  # state-arm (shared[1]) SPDE term to the dynamic path, and the term
  # constructor maps it to the psi1 predictor.
  spatial_occ <- .spatial_for_arm(spatial, 1L)

  # Precompute per site-season
  nv <- model$n_visits
  ad <- model$any_detected

  e_step <- function(fits, ...) {
    beta_occ <- extract_beta(fits$occ, p_occ)
    beta_det <- extract_beta(fits$det, p_det)
    beta_col <- extract_beta(fits$col, p_col)
    beta_ext <- extract_beta(fits$ext, p_ext)

    eta_psi1 <- as.vector(X_occ %*% beta_occ)
    sp_off <- .spatial_eta_offset(spatial_occ, fits$occ, p_occ)
    if (length(sp_off) == n_sites) eta_psi1 <- eta_psi1 + sp_off
    psi1 <- plogis(eta_psi1)
    p <- plogis(as.vector(X_det %*% beta_det))
    gam <- plogis(as.vector(X_col %*% beta_col))
    eps <- plogis(as.vector(X_ext %*% beta_ext))

    # HMM forward pass to get P(z_it = 1 | y)
    w <- matrix(NA_real_, n_sites, n_seasons)
    for (i in seq_len(n_sites)) {
      alpha_occ <- psi1[i]; alpha_unocc <- 1 - psi1[i]
      for (t in seq_len(n_seasons)) {
        idx <- (i - 1) * n_seasons + t
        nv_it <- nv[idx]; det_it <- ad[idx]
        if (nv_it > 0) {
          prob_y_occ <- if (det_it) 1 else (1 - p[i])^nv_it
          prob_y_unocc <- if (det_it) 0 else 1
          post_occ <- alpha_occ * prob_y_occ
          post_unocc <- alpha_unocc * prob_y_unocc
          total <- post_occ + post_unocc
          w[i, t] <- post_occ / total
          alpha_occ <- post_occ / total
          alpha_unocc <- post_unocc / total
        } else {
          w[i, t] <- alpha_occ
        }
        if (t < n_seasons) {
          new_occ <- alpha_occ * (1 - eps[i]) + alpha_unocc * gam[i]
          new_unocc <- alpha_occ * eps[i] + alpha_unocc * (1 - gam[i])
          alpha_occ <- new_occ; alpha_unocc <- new_unocc
        }
      }
    }
    list(weights = w)
  }

  m_step_encode <- function(weights, ...) {
    w <- weights  # n_sites x n_seasons matrix
    # Occupancy: psi1 from season 1 weights. The pseudo-binomial inflation is
    # M = 1000 without a field; with the psi1 SPDE field it drops to M = 4 so
    # the field prior precision is not swamped by the data signal (the same
    # modest-M encoding the single-season state arm uses).
    M <- 1000L
    M_occ <- if (is.null(spatial_occ)) 1000L else 4L
    w1 <- w[, 1]
    y_occ <- ifelse(ad[seq(1, by = n_seasons, length.out = n_sites)], M_occ,
                    as.integer(round(w1 * M_occ)))
    y_occ <- pmin(pmax(y_occ, 0L), M_occ)

    # Colonization: from transitions where z_{t-1}=0, z_t=1
    # Extinction: from transitions where z_{t-1}=1, z_t=0
    # Approximate: site-level average
    col_y <- integer(n_sites); col_n <- integer(n_sites)
    ext_y <- integer(n_sites); ext_n <- integer(n_sites)
    for (i in seq_len(n_sites)) {
      for (t in 2:n_seasons) {
        p_prev_occ <- w[i, t - 1]
        p_curr_occ <- w[i, t]
        # P(colonization event) ~= (1-w_{t-1}) * w_t
        col_y[i] <- col_y[i] + as.integer(round((1 - p_prev_occ) * p_curr_occ * M))
        col_n[i] <- col_n[i] + as.integer(round((1 - p_prev_occ) * M))
        ext_y[i] <- ext_y[i] + as.integer(round(p_prev_occ * (1 - p_curr_occ) * M))
        ext_n[i] <- ext_n[i] + as.integer(round(p_prev_occ * M))
      }
    }
    col_n <- pmax(col_n, 1L); ext_n <- pmax(ext_n, 1L)
    col_y <- pmin(pmax(col_y, 0L), col_n)
    ext_y <- pmin(pmax(ext_y, 0L), ext_n)

    # Detection: per-(site, season) rows weighted by w[i, t] = P(z_it = 1 | y).
    # Replaces the legacy hard threshold (w > 0.5) which silently dropped
    # site-seasons in the boundary regime and double-counted detection
    # evidence for site-seasons in the high-confidence regime. X_det is
    # site-indexed in this model, so per-season rows just replicate the
    # site's covariates.
    rows_i <- integer(n_sites * n_seasons)
    det_count <- integer(n_sites * n_seasons)
    vis_count <- integer(n_sites * n_seasons)
    w_it <- numeric(n_sites * n_seasons)
    n_rows <- 0L
    for (i in seq_len(n_sites)) {
      for (t in seq_len(n_seasons)) {
        idx <- (i - 1) * n_seasons + t
        if (nv[idx] <= 0) next
        base <- (i - 1) * n_seasons * max_visits + (t - 1) * max_visits
        dc <- 0L; vc <- 0L
        for (j in seq_len(nv[idx])) {
          v <- y_flat[base + j]
          if (v >= 0) { vc <- vc + 1L; if (v == 1) dc <- dc + 1L }
        }
        if (vc == 0L) next
        w_eff <- if (dc > 0L) 1 else w[i, t]
        if (w_eff <= 1e-6) next
        n_rows <- n_rows + 1L
        rows_i[n_rows] <- i
        det_count[n_rows] <- dc
        vis_count[n_rows] <- vc
        w_it[n_rows] <- w_eff
      }
    }
    if (n_rows > 0L) {
      rows_i <- rows_i[seq_len(n_rows)]
      det_count <- det_count[seq_len(n_rows)]
      vis_count <- vis_count[seq_len(n_rows)]
      w_it <- w_it[seq_len(n_rows)]
      det_block <- list(y = det_count, n_trials = vis_count,
                        X = X_det[rows_i, , drop = FALSE],
                        weights = w_it, family = "binomial")
    } else {
      det_block <- list(y = integer(0), n_trials = integer(0),
                        X = X_det[integer(0), , drop = FALSE],
                        weights = numeric(0), family = "binomial")
    }

    occ_block <- list(y = y_occ, n_trials = rep(M_occ, n_sites), X = X_occ,
                      family = "binomial")
    occ_block <- .attach_spatial_spde(occ_block, spatial_occ)

    list(
      occ = occ_block,
      det = det_block,
      col = list(y = col_y, n_trials = col_n, X = X_col,
                 family = "binomial"),
      ext = list(y = ext_y, n_trials = ext_n, X = X_ext,
                 family = "binomial")
    )
  }

  z_draw <- function(weights, ...) {
    w <- weights
    z <- matrix(0L, n_sites, n_seasons)
    for (i in seq_len(n_sites)) {
      for (t in seq_len(n_seasons)) {
        idx <- (i - 1) * n_seasons + t
        if (ad[idx]) z[i, t] <- 1L
        else z[i, t] <- rbinom(1, 1, clamp_w(w[i, t]))
      }
    }
    z
  }

  hard_encode <- function(z, ...) {
    z1 <- z[, 1]
    # Colonization/extinction from hard transitions
    col_y <- integer(n_sites); col_n <- integer(n_sites)
    ext_y <- integer(n_sites); ext_n <- integer(n_sites)
    for (i in seq_len(n_sites)) {
      for (t in 2:n_seasons) {
        if (z[i, t - 1] == 0) { col_n[i] <- col_n[i] + 1L; if (z[i, t] == 1) col_y[i] <- col_y[i] + 1L }
        if (z[i, t - 1] == 1) { ext_n[i] <- ext_n[i] + 1L; if (z[i, t] == 0) ext_y[i] <- ext_y[i] + 1L }
      }
    }
    col_n <- pmax(col_n, 1L); ext_n <- pmax(ext_n, 1L)

    total_det <- integer(n_sites); total_vis <- integer(n_sites)
    for (i in seq_len(n_sites)) {
      for (t in seq_len(n_seasons)) {
        if (z[i, t] == 1 && nv[(i-1)*n_seasons+t] > 0) {
          base <- (i-1)*n_seasons*max_visits + (t-1)*max_visits
          for (j in seq_len(nv[(i-1)*n_seasons+t])) {
            v <- y_flat[base + j]
            if (v >= 0) { total_vis[i] <- total_vis[i]+1L; if (v==1) total_det[i] <- total_det[i]+1L }
          }
        }
      }
    }
    dk <- total_vis > 0

    occ_block <- list(y = z1, n_trials = rep(1L, n_sites), X = X_occ,
                      family = "binomial")
    occ_block <- .attach_spatial_spde(occ_block, spatial_occ)
    list(
      occ = occ_block,
      det = if (sum(dk) > 0)
        list(y = total_det[dk], n_trials = total_vis[dk],
             X = X_det[dk,,drop=FALSE], family = "binomial")
      else NULL,
      col = list(y = col_y, n_trials = col_n, X = X_col,
                 family = "binomial"),
      ext = list(y = ext_y, n_trials = ext_n, X = X_ext,
                 family = "binomial")
    )
  }

  init <- list(
    occ = list(beta = rep(0, p_occ), se = rep(1, p_occ)),
    det = list(beta = rep(0, p_det), se = rep(1, p_det)),
    col = list(beta = rep(0, p_col), se = rep(1, p_col)),
    ext = list(beta = rep(0, p_ext), se = rep(1, p_ext))
  )

  list(e_step = e_step, m_step_encode = m_step_encode, z_draw = z_draw,
       hard_encode = hard_encode, init = init,
       p_per_submodel = c(occ = p_occ, det = p_det, col = p_col, ext = p_ext))
}

# ============================================================================
# Integrated occupancy callbacks
# ============================================================================
build_integrated_callbacks <- function(model, spatial = NULL) {
  y_sources <- model$y_sources
  site_maps <- model$site_maps
  X_occ <- model$X_processes[[1]]
  n_sites <- model$n_sites
  n_sources <- model$n_sources
  p_occ <- ncol(X_occ)

  # The shared psi field enters the state arm (one state row per site, the
  # identity map -- the proven single-season path). A field on the detection
  # arm enters every source's per-source detection block, broadcast onto that
  # source's sites (`src_rows`) and folded into the response by count-scaling,
  # exactly as the single-season detection arm does. The field is fit
  # independently per source block (one realization per submodel block; a
  # genuinely shared realization across sources needs the copy() path).
  spatial_occ <- .spatial_for_arm(spatial, 1L)
  spatial_det <- .spatial_for_arm(spatial, 2L)

  # Per-source detection info
  src_info <- lapply(seq_len(n_sources), function(s) {
    ys <- y_sources[[s]]; ns <- nrow(ys); mv <- ncol(ys)
    nv <- integer(ns); nd <- integer(ns); ad <- logical(ns)
    for (i in seq_len(ns)) {
      valid <- ys[i, ] >= 0; nv[i] <- sum(valid)
      nd[i] <- sum(ys[i, valid] == 1); ad[i] <- nd[i] > 0
    }
    X_det <- model$X_processes[[1 + s]]
    src_rows <- site_maps[[s]] + 1L
    # Detection-arm field projection broadcast onto this source's sites.
    spatial_det_s <- if (!is.null(spatial_det))
      .tobs_spde_broadcast_spec(spatial_det, src_rows) else NULL
    list(nv = nv, nd = nd, ad = ad, X_det = X_det[src_rows, , drop = FALSE],
         p_det = ncol(X_det), keep = nv > 0, src_rows = src_rows,
         spatial_det = spatial_det_s)
  })

  # Global detection status per site
  any_det_global <- logical(n_sites)
  for (s in seq_len(n_sources)) {
    for (j in seq_along(src_info[[s]]$src_rows)) {
      if (src_info[[s]]$ad[j]) any_det_global[src_info[[s]]$src_rows[j]] <- TRUE
    }
  }

  e_step <- function(fits, ...) {
    beta_occ <- extract_beta(fits$occ, p_occ)
    eta_occ <- as.vector(X_occ %*% beta_occ)
    sp_off <- .spatial_eta_offset(spatial_occ, fits$occ, p_occ)
    if (length(sp_off) == n_sites) eta_occ <- eta_occ + sp_off
    psi <- plogis(eta_occ)
    weights <- psi  # Prior occupancy
    for (s in seq_len(n_sources)) {
      si <- src_info[[s]]
      beta_det <- extract_beta(fits[[paste0("det", s)]], si$p_det)
      eta_det <- as.vector(si$X_det %*% beta_det)
      det_off <- .spatial_eta_offset(si$spatial_det, fits[[paste0("det", s)]],
                                     si$p_det)
      if (length(det_off) == length(si$src_rows)) eta_det <- eta_det + det_off
      p_s <- plogis(eta_det)
      for (j in seq_along(si$src_rows)) {
        i <- si$src_rows[j]
        if (si$ad[j]) { weights[i] <- 1 }
        else if (si$nv[j] > 0) {
          prod_1mp <- (1 - p_s[j])^si$nv[j]
          num <- weights[i] * prod_1mp
          weights[i] <- num / (num + (1 - weights[i]) + 1e-10)
        }
      }
    }
    list(weights = weights)
  }

  m_step_encode <- function(weights, ...) {
    if (is.null(spatial_occ)) {
      M <- 1000L
      y_occ <- ifelse(any_det_global, M, as.integer(round(weights * M)))
      y_occ <- pmin(pmax(y_occ, 0L), M)
      occ_block <- list(y = y_occ, n_trials = rep(M, n_sites), X = X_occ,
                        family = "binomial")
    } else {
      M <- 4L
      y_occ <- ifelse(any_det_global, M, as.integer(round(weights * M)))
      y_occ <- pmin(pmax(y_occ, 0L), M)
      occ_block <- list(y = y_occ, n_trials = rep(M, n_sites), X = X_occ,
                        family = "binomial")
      occ_block <- .attach_spatial_spde(occ_block, spatial_occ)
    }
    specs <- list(occ = occ_block)
    # Per-source detection blocks: weight each row by w_i at the global
    # site mapped through src_rows. Sites where the E-step says "almost
    # certainly empty" drop out of every source's detection fit.
    for (s in seq_len(n_sources)) {
      si <- src_info[[s]]
      w_src <- weights[si$src_rows]
      w_src[si$ad] <- 1
      if (!is.null(si$spatial_det)) {
        # SPDE detection field on this source: fold the occupancy weight into
        # the binomial response by count-scaling (y = round(w nd), n =
        # round(w nv)) so every source row stays aligned with the broadcast
        # mesh projection A; a near-empty site collapses to a (0, 0) row.
        y_det_w <- as.integer(round(w_src * si$nd))
        n_det_w <- as.integer(round(w_src * si$nv))
        y_det_w <- pmin(pmax(y_det_w, 0L), n_det_w)
        det_block <- list(y = y_det_w, n_trials = n_det_w, X = si$X_det,
                          family = "binomial")
        specs[[paste0("det", s)]] <- .attach_spatial_spde(det_block,
                                                          si$spatial_det)
      } else {
        dk <- si$keep & (w_src > 1e-6)
        specs[[paste0("det", s)]] <- list(y = si$nd[dk], n_trials = si$nv[dk],
                                          X = si$X_det[dk,,drop=FALSE],
                                          weights = w_src[dk],
                                          family = "binomial")
      }
    }
    specs
  }

  z_draw <- function(weights, ...) {
    z <- as.integer(any_det_global)
    undet <- !any_det_global
    z[undet] <- rbinom(sum(undet), 1, clamp_w(weights[undet]))
    z
  }

  hard_encode <- function(z, ...) {
    occ_block <- list(y = z, n_trials = rep(1L, n_sites), X = X_occ,
                      family = "binomial")
    occ_block <- .attach_spatial_spde(occ_block, spatial_occ)
    specs <- list(occ = occ_block)
    for (s in seq_len(n_sources)) {
      si <- src_info[[s]]
      if (!is.null(si$spatial_det)) {
        keep_det_mask <- (z[si$src_rows] == 1L) & (si$nv > 0L)
        n_det_h <- ifelse(keep_det_mask, si$nv, 0L)
        y_det_h <- ifelse(keep_det_mask, si$nd, 0L)
        det_block <- list(y = as.integer(y_det_h),
                          n_trials = as.integer(n_det_h),
                          X = si$X_det, family = "binomial")
        specs[[paste0("det", s)]] <- .attach_spatial_spde(det_block,
                                                          si$spatial_det)
      } else {
        occ_local <- z[si$src_rows] == 1L & si$nv > 0
        if (any(occ_local)) {
          specs[[paste0("det", s)]] <- list(y = si$nd[occ_local],
                                            n_trials = si$nv[occ_local],
                                            X = si$X_det[occ_local,,drop=FALSE],
                                            family = "binomial")
        }
      }
    }
    specs
  }

  init <- list(occ = list(beta = rep(0, p_occ), se = rep(1, p_occ)))
  p_sub <- c(occ = p_occ)
  for (s in seq_len(n_sources)) {
    nm <- paste0("det", s)
    init[[nm]] <- list(beta = rep(0, src_info[[s]]$p_det), se = rep(1, src_info[[s]]$p_det))
    p_sub[nm] <- src_info[[s]]$p_det
  }

  list(e_step = e_step, m_step_encode = m_step_encode, z_draw = z_draw,
       hard_encode = hard_encode, init = init, p_per_submodel = p_sub)
}

# ============================================================================
# JSDM callbacks (no detection — simple Bernoulli)
# ============================================================================
build_jsdm_callbacks <- function(model, spatial = NULL) {
  y_jsdm <- model$y_jsdm
  X_occ <- model$X_processes[[1]]
  N <- model$N
  p_occ <- ncol(X_occ)

  # JSDM has no detection process: y is observed directly, so the state arm is
  # a plain Bernoulli on the observed presence/absence. A site-level field is
  # shared across the species at a site -- the state block carries
  # N = n_sites * n_species rows in site-major order, so the site-indexed mesh
  # projection A is broadcast onto those rows via `site_of_row`. No occupancy
  # weight to fold in (no latent state), so the response is the observed y at
  # unit trials with the broadcast field attached directly.
  spatial_occ <- .spatial_for_arm(spatial, 1L)
  if (!is.null(spatial_occ)) {
    site_of_row <- .tobs_state_block_dims(model)$site_of_row
    spatial_occ <- .tobs_spde_broadcast_spec(spatial_occ, site_of_row)
  }

  # No E-step needed — no latent variable (y is observed directly)
  # But we still use EM framework for consistency with species RE
  e_step <- function(fits, ...) {
    list(weights = as.numeric(y_jsdm))
  }

  m_step_encode <- function(weights, ...) {
    occ_block <- list(y = as.integer(y_jsdm), n_trials = rep(1L, N), X = X_occ,
                      family = "binomial")
    occ_block <- .attach_spatial_spde(occ_block, spatial_occ)
    list(occ = occ_block)
  }

  z_draw <- function(weights, ...) as.integer(y_jsdm)
  hard_encode <- function(z, ...) {
    occ_block <- list(y = z, n_trials = rep(1L, N), X = X_occ,
                      family = "binomial")
    occ_block <- .attach_spatial_spde(occ_block, spatial_occ)
    list(occ = occ_block)
  }

  init <- list(occ = list(beta = rep(0, p_occ), se = rep(1, p_occ)))

  list(e_step = e_step, m_step_encode = m_step_encode, z_draw = z_draw,
       hard_encode = hard_encode, init = init, p_per_submodel = c(occ = p_occ))
}

# ============================================================================
# Shared helpers
# ============================================================================

# SPDE (continuous mesh field) coverage on the .tobs_laplace path, by model
# type and arm (each cell is "wired" or an honest stop()):
#
#   model_type   state arm (shared[1])   detection arm (shared[2])
#   single       wired                   wired
#   integrated   wired                   wired (per source)
#   jsdm         wired                   n/a (no detection process)
#   community    wired                   stop()
#   dynamic      wired (psi1 only)       stop()
#
# The state field on jsdm / community broadcasts the site-indexed mesh
# projection A onto the N = n_sites * n_species state rows (one shared
# site-level field across the species at a site, via .tobs_spde_broadcast_spec
# and .tobs_state_block_dims). The dynamic state field enters season-1 psi1
# only; the colonization / extinction transition predictors are separate latent
# processes whose own mesh fields are not wired (a state-arm spde() term maps to
# psi1). A single field shared across both arms at once (shared = c(TRUE, TRUE))
# is a stop() everywhere: the single-Laplace block fitter fits one field
# realization per submodel block, so a genuinely shared realization needs the
# copy() path, not two independent blocks. The areal path (icar/bym2/car_proper
# via nested_laplace) is wider; this matrix is the continuous-mesh SPDE path
# only. The continuous gp()/spde() fields on the N-mixture arms (abun /
# em_nested / ms_abun) are tracked separately (gcol33/tulpaObs#21).

# Validate that `spatial` (a `tobs_spatial` or NULL) can be consumed by the
# Laplace path. Wired by arm and model type:
#   state arm (shared[1]): single, integrated, jsdm, community, dynamic (psi1)
#   detection arm (shared[2]): single, integrated (per source)
# jsdm has no detection process; community and dynamic do not yet carry a
# detection-arm field. A single realization shared across both arms at once
# needs the copy() path (the single-Laplace block fitter fits one realization
# per submodel block), so c(TRUE, TRUE) errors here rather than silently fitting
# two independent fields. Other combinations error explicitly.
.validate_spatial_laplace <- function(spatial, model_type) {
  if (is.null(spatial)) return(invisible())
  if (!inherits(spatial, "tobs_spatial")) {
    stop("spatial must be a tobs_spatial term (from a spde() formula term)",
         call. = FALSE)
  }
  if (!identical(spatial$type, "spde")) {
    stop(sprintf(
      ".tobs_laplace currently supports spatial$type == 'spde' only (got '%s'). Use method = 'nuts' for other spatial types.",
      spatial$type), call. = FALSE)
  }
  on_occ <- isTRUE(spatial$shared[1])
  on_det <- length(spatial$shared) >= 2 && isTRUE(spatial$shared[2])
  if (!on_occ && !on_det) {
    stop("SPDE must be attached to the occupancy/state or detection submodel.",
         call. = FALSE)
  }
  if (on_occ && on_det) {
    stop("A single SPDE field shared across the occupancy and detection arms is not plumbed in .tobs_laplace; attach the field to one arm, or use method = 'nuts'.",
         call. = FALSE)
  }
  if (on_det && model_type %in% c("jsdm", "community", "dynamic")) {
    stop(sprintf(
      "SPDE on the detection process is plumbed for single-season and integrated occupancy only in .tobs_laplace (got model_type = '%s'). Attach the field to the state arm, or use method = 'nuts'.",
      model_type), call. = FALSE)
  }
  invisible()
}

# Gate the deterministic random-effect path. The variance-component EM in
# R/em_laplace_re.R fits iid intercept, uncorrelated slopes, and correlated
# slopes (a full RE covariance) on EITHER the occupancy or the detection
# predictor of a single-season model (each arm carries its own RE block). Forms
# it cannot fit -- non-single families, RE + spatial, RE + visit-level
# detection, a single RE shared across both predictors -- error here with a
# pointer to `method = "nuts"` (which fits every RE form) rather than being
# silently dropped (gcol33/tulpaObs#11). The raw EM variance components (sigma,
# correlation) carry the Laplace small-cluster bias for binary data (the glmer
# nAGQ=1 regime, not Breslow-Clayton PQL); the default re.aghq = TRUE refines
# them on the exact-marginal adaptive Gauss-Hermite likelihood (R/re_aghq.R),
# removing the attenuation, with a default LKJ(re.lkj = 1.5) penalty
# regularizing a weakly-identified RE correlation off the +-1 boundary.
.validate_re_laplace <- function(re, model, spatial, approx) {
  re_list <- if (inherits(re, "tobs_re")) list(re) else re

  if (!identical(model$model_type, "single")) {
    stop(sprintf(
      "Random effects under method = 'laplace' are wired for single-season occupancy only (got model_type = '%s'). Use method = 'nuts' for random effects on this family.",
      model$model_type), call. = FALSE)
  }
  if (!is.null(spatial)) {
    stop("A random effect combined with a spatial term is not supported on the Laplace path. Use method = 'nuts'.",
         call. = FALSE)
  }
  if (!is.null(model$X_det_visit)) {
    stop("Random effects with visit-level detection covariates are not supported on the Laplace path. Use method = 'nuts'.",
         call. = FALSE)
  }
  for (r in re_list) {
    if (length(r$shared) >= 2L && isTRUE(r$shared[1]) && isTRUE(r$shared[2])) {
      stop("A single random effect shared across occupancy and detection is not supported on the Laplace path (each arm fits its own RE block). Use method = 'nuts'.",
           call. = FALSE)
    }
  }
  invisible()
}

# Linear-predictor offset induced by the SPDE mesh field at the current fit.
# After tulpa_laplace returns mode = c(beta, u_mesh), the spatial contribution
# to eta at the observed locations is A %*% u_mesh.
.spatial_eta_offset <- function(spatial, fits_sub, p_fixed) {
  if (is.null(spatial) || is.null(fits_sub) || is.null(fits_sub$mode)) {
    return(rep(0, 0))
  }
  if (!identical(spatial$type, "spde")) return(rep(0, 0))
  mode_vec <- fits_sub$mode
  if (length(mode_vec) <= p_fixed) return(rep(0, 0))
  u <- mode_vec[(p_fixed + 1L):length(mode_vec)]
  as.numeric(spatial$tulpa_spec$A %*% u)
}

# Attach the tulpa-side spatial spec to an M-step block. The block's `spatial`
# field is forwarded as-is by tulpa_em_laplace -> tulpa_laplace.
.attach_spatial_spde <- function(block, spatial) {
  if (is.null(spatial) || !identical(spatial$type, "spde")) return(block)
  block$spatial <- spatial$tulpa_spec
  block
}

# Broadcast a site-indexed SPDE field onto a state block whose rows are
# (site, species) -- community and jsdm carry N = n_sites * n_species rows
# ordered site-major, so a single site-level field is shared across the species
# at a site. The mesh / FEM matrices (C, G, n_mesh, nu, priors) describe the
# field on the sites and are unchanged; only the projection A (and its
# pre-extracted CSC slots A_x / A_i / A_p, which `laplace_spde_at()` hands to
# the SPDE solver) is re-rowed so row r of the block projects through the mesh
# basis at `site_of_row[r]`. The eta offset reader `.spatial_eta_offset()` uses
# the same broadcast A, so the field contribution lands on every species row of
# a site identically. `site_of_row` is 1-based into the n_sites mesh rows.
.tobs_spde_broadcast_spec <- function(spatial, site_of_row) {
  if (is.null(spatial) || !identical(spatial$type, "spde")) return(spatial)
  sp <- spatial$tulpa_spec
  A_b <- sp$A[site_of_row, , drop = FALSE]
  A_csc <- methods::as(A_b, "CsparseMatrix")
  sp$A   <- A_b
  sp$A_x <- A_csc@x
  sp$A_i <- A_csc@i
  sp$A_p <- A_csc@p
  spatial$tulpa_spec <- sp
  spatial
}

# Select the SPDE spec for a given arm (1 = occupancy/state, 2 = detection)
# from a `tobs_spatial` term carrying a `$shared = c(occ, det)` membership.
# Returns the spec when the arm carries the field, NULL otherwise. One mesh
# term is shared across the arms it enters, so each arm references the same
# `tulpa_spec` (the field realization is fit independently per arm here -- two
# separate mesh blocks, not a copied realization, which the single-Laplace path
# does not support across submodels).
.spatial_for_arm <- function(spatial, arm) {
  if (is.null(spatial) || !identical(spatial$type, "spde")) return(NULL)
  sh <- spatial$shared
  if (length(sh) >= arm && isTRUE(sh[arm])) spatial else NULL
}

# TRUE when a nested-Laplace multi-block latent prior carries a continuous
# Matern (SPDE) field block. Used to switch the occupancy M-step to the modest
# pseudo-binomial inflation that keeps the mesh-field prior from being swamped.
.tobs_latent_prior_has_spde <- function(latent_prior) {
  if (is.null(latent_prior)) return(FALSE)
  blocks <- if (!is.null(latent_prior$type)) list(latent_prior) else latent_prior
  any(vapply(blocks, function(b) identical(b$type, "spde"), logical(1)))
}

# NUTS sampler-health diagnostics for a non-sampled (Laplace / nested-Laplace)
# fit. No HMC trajectory exists, so acceptance, divergence, tree depth and the
# integrator step size are unavailable; they are NA rather than 0/1 so a user
# inspecting sampler health does not read "no sampler ran" as "sampler ran
# cleanly" (NA-on-unavailable, the same rule as .se_from_laplace_fit). Splice
# into a tobs_fit build with !!! / do.call so the named fields land directly.
.tobs_na_nuts_diagnostics <- function(n_draws) {
  list(
    accept_prob = rep(NA_real_, n_draws),
    divergent   = rep(NA_real_, n_draws),
    treedepth   = rep(NA_integer_, n_draws),
    epsilon     = NA_real_
  )
}

extract_beta <- function(sub, p) {
  if (is.null(sub)) return(rep(0, p))
  if (!is.null(sub$beta)) return(sub$beta)
  if (!is.null(sub$mean)) return(sub$mean)
  if (!is.null(sub$mode)) return(sub$mode[seq_len(p)])
  rep(0, p)
}

# SE for the fixed-effect block of a tulpa_laplace() fit. Reads the
# negative-log-posterior Hessian (`H_beta`, the precision matrix), inverts
# it, and returns sqrt(diag(.)) restricted to the first `p` fixed effects.
# When the inner fit had a spatial mesh field attached (`spde` / `gp`),
# tulpa_laplace skips H_beta — return NA so callers can flag the
# uncertainty as unavailable instead of carrying a placeholder.
.se_from_laplace_fit <- function(fi, p) {
  if (!is.null(fi$se)) {
    se <- as.numeric(fi$se)
    if (length(se) >= p) return(se[seq_len(p)])
  }
  H <- fi$H_beta
  if (is.null(H)) return(rep(NA_real_, p))
  cov <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(cov)) return(rep(NA_real_, p))
  d <- sqrt(pmax(diag(cov), 0))
  if (length(d) >= p) return(d[seq_len(p)])
  c(d, rep(NA_real_, p - length(d)))
}

# Louis-corrected observed Fisher info for the occupancy fixed-effect block
# of a single-season occu fit (tulpaObs#7).
#
# Why this is needed. The inner M-step encodes the soft-imputed P(z_i = 1 | y_i)
# as a pseudo-binomial likelihood with n_trials = M (M = 1000 non-spatial,
# M = 4 spatial). The resulting inner Hessian is
#
#   H_inner = M * X' diag(psi (1 - psi)) X + P_prior
#
# i.e. the complete-data Fisher info inflated by the M trick, plus the prior
# precision. This is the wrong object for SE reporting on two counts: the M
# factor is an artefact of the M-step encoding (not data information), and
# the complete-data info ignores the missing-z variance.
#
# Louis identity for the occupancy score s_i = x_i (z_i - psi_i) gives the
# observed Fisher info at the EM stationary point:
#
#   I_obs(beta_psi) = E[-d2 log f / dbeta2 | y] - Var(s_complete | y)
#                   = X' diag(psi (1 - psi)) X - X' diag(w (1 - w)) X
#                   = X' diag(psi (1 - psi) - w (1 - w)) X
#
# where w_i = P(z_i = 1 | y_i, theta_hat) is the converged E-step weight. The
# per-site `psi(1-psi) - w(1-w)` term can be negative (the marginal log-lik
# can be locally convex at a single site), but the aggregate X' D X is PSD at
# the MLE because it equals minus the marginal log-lik Hessian at its max.
.louis_info_psi_single <- function(X_occ, beta_psi, weights,
                                   spatial = NULL, spatial_fit = NULL,
                                   prior_spec = NULL,
                                   coef_names = NULL) {
  p_psi <- length(beta_psi)
  if (p_psi == 0L) return(NULL)
  if (is.null(X_occ) || nrow(X_occ) == 0L) return(NULL)
  if (is.null(weights) || length(weights) != nrow(X_occ)) return(NULL)

  eta <- as.numeric(X_occ %*% beta_psi)
  sp_off <- .spatial_eta_offset(spatial, spatial_fit, p_psi)
  if (length(sp_off) == nrow(X_occ)) eta <- eta + sp_off
  eta <- pmin(pmax(eta, -30), 30)
  psi <- plogis(eta)

  d <- psi * (1 - psi) - weights * (1 - weights)
  I_obs <- as.matrix(crossprod(X_occ, d * X_occ))

  if (!is.null(prior_spec)) {
    if (is.null(coef_names)) coef_names <- colnames(X_occ) %||% paste0("x", seq_len(p_psi))
    pr <- .prior_for_submodel(prior_spec, "psi", coef_names)
    if (!is.null(pr)) {
      pen_prec <- ifelse(is.finite(pr$sd), 1 / (pr$sd^2), 0)
      diag(I_obs) <- diag(I_obs) + pen_prec[seq_len(p_psi)]
    }
  }
  I_obs
}

# Marginal Louis observed Fisher info for the SITE-LEVEL detection block of a
# single-season fit (tulpaObs#7, detection arm). The detection M-step fits a
# weighted binomial whose returned H_beta = X_det' diag(w n_valid p(1-p)) X_det
# is the *complete-data* info: it treats the soft occupancy weight
# w_i = P(z_i = 1 | y) as known and so under-states the SE (the M-step Hessian
# is the wrong object for SEs, exactly as on the psi arm). The occupancy and
# detection estimating equations both depend on the latent z, so the detection
# SE must come from the JOINT (psi, det) Louis observed info, marginalized over
# psi -- the diagonal det block alone fixes the intercept but leaves the slope
# under-dispersed.
#
# Complete-data scores: s_psi,i = x_psi,i (z_i - psi_i),
# s_det,i = z_i (n_det_i - n_valid_i p_i) x_det,i. With z_i | y ~ Bern(w_i) the
# Louis identity (E[I_complete | y] - Var(s_complete | y)) gives the joint
# observed info in three blocks (the complete-data cross block is 0):
#
#   I_pp = X_psi' diag( psi(1-psi) - w(1-w) ) X_psi                 (+ psi prior)
#   I_dd = X_det' diag( w n_valid p(1-p) - (n_valid p)^2 w(1-w) ) X_det (+ p prior)
#   I_pd = - X_psi' diag( n_valid p w(1-w) ) X_det
#
# (a detected site has w_i = 1 so its w(1-w) terms vanish.) The marginal
# detection info is the Schur complement I_dd - I_pd' I_pp^{-1} I_pd, whose
# inverse is the (beta_det) block of the full joint covariance.
.louis_info_det_single <- function(X_occ, beta_psi, X_det, beta_det,
                                   weights, n_valid, prior_spec = NULL,
                                   occ_coef_names = NULL, det_coef_names = NULL,
                                   spatial = NULL, spatial_fit = NULL) {
  p_det <- length(beta_det)
  p_psi <- length(beta_psi)
  if (p_det == 0L) return(NULL)
  if (is.null(X_det) || nrow(X_det) == 0L) return(NULL)
  if (is.null(weights) || length(weights) != nrow(X_det)) return(NULL)
  if (is.null(n_valid) || length(n_valid) != nrow(X_det)) return(NULL)

  w  <- weights
  nv <- as.numeric(n_valid)
  p  <- plogis(pmin(pmax(as.numeric(X_det %*% beta_det), -30), 30))

  add_prior <- function(I, arm, p_k, coef_names, Xcols) {
    if (is.null(prior_spec)) return(I)
    if (is.null(coef_names)) coef_names <- Xcols %||% paste0("x", seq_len(p_k))
    pr <- .prior_for_submodel(prior_spec, arm, coef_names)
    if (!is.null(pr)) {
      pen <- ifelse(is.finite(pr$sd), 1 / (pr$sd^2), 0)
      diag(I) <- diag(I) + pen[seq_len(p_k)]
    }
    I
  }

  # Detection diagonal block.
  d_dd <- w * nv * p * (1 - p) - (nv * p)^2 * w * (1 - w)
  d_dd[nv <= 0] <- 0
  I_dd <- add_prior(as.matrix(crossprod(X_det, d_dd * X_det)),
                    "p", p_det, det_coef_names, colnames(X_det))

  # Couple with the occupancy block via the joint Louis cross term, then
  # marginalize psi out by Schur complement. Skip when the occupancy inputs are
  # unavailable (fall back to the diagonal block, which still fixes the level).
  if (p_psi > 0L && !is.null(X_occ) && nrow(X_occ) == nrow(X_det)) {
    eta_o <- as.numeric(X_occ %*% beta_psi)
    sp_off <- .spatial_eta_offset(spatial, spatial_fit, p_psi)
    if (length(sp_off) == nrow(X_occ)) eta_o <- eta_o + sp_off
    psi <- plogis(pmin(pmax(eta_o, -30), 30))

    d_pp <- psi * (1 - psi) - w * (1 - w)
    I_pp <- add_prior(as.matrix(crossprod(X_occ, d_pp * X_occ)),
                      "psi", p_psi, occ_coef_names, colnames(X_occ))
    d_pd <- -nv * p * w * (1 - w)
    d_pd[nv <= 0] <- 0
    I_pd <- as.matrix(crossprod(X_occ, d_pd * X_det))    # p_psi x p_det

    schur <- tryCatch(I_dd - crossprod(I_pd, solve(I_pp, I_pd)),
                      error = function(e) NULL)
    if (!is.null(schur)) I_dd <- schur
  }
  I_dd
}

# SE vector from an observed-info matrix; returns NA of length p on failure.
.se_from_info <- function(I, p) {
  if (is.null(I)) return(rep(NA_real_, p))
  cov <- tryCatch(solve(I), error = function(e) NULL)
  if (is.null(cov)) return(rep(NA_real_, p))
  d <- sqrt(pmax(diag(cov), 0))
  if (length(d) >= p) d[seq_len(p)] else c(d, rep(NA_real_, p - length(d)))
}

# Within-arm covariance block carrying the off-diagonal correlation of the
# observed-information inverse, rescaled so its diagonal matches the reported
# marginal SEs exactly (`sds`). `prec` is the precision / observed-information
# matrix the SEs were derived from (cov = solve(prec)); NULL or a non-invertible
# `prec` yields the diagonal block diag(sds^2), i.e. the previous behaviour. This
# keeps the marginal SEs byte-identical while restoring the joint correlation the
# diagonal pseudo-draws used to discard (gcol33/tulpaObs#44). NA / non-finite SEs
# map to a zero-variance coordinate (drawn as a point mass downstream).
.cor_scaled_cov <- function(prec, sds) {
  p <- length(sds)
  s <- ifelse(is.finite(sds), sds, 0)
  if (is.null(prec)) return(diag(s^2, nrow = p))
  cov <- tryCatch(solve(prec), error = function(e) NULL)
  if (is.null(cov)) return(diag(s^2, nrow = p))
  cov <- as.matrix(cov)[seq_len(p), seq_len(p), drop = FALSE]
  d <- sqrt(pmax(diag(cov), 0))
  ok <- is.finite(d) & d > 0
  R <- diag(p)
  if (sum(ok) > 1L) R[ok, ok] <- cov[ok, ok, drop = FALSE] / tcrossprod(d[ok])
  outer(s, s) * R
}

# Assemble a block-diagonal covariance from per-block matrices in append order,
# flooring zero / NA variances so the matrix is PD and chol-decomposable in
# .rmvn (mirrors the old `max(sd_j, 1e-4)` point-mass floor for NA-SE columns).
.assemble_block_diag <- function(blocks, n_params) {
  V <- matrix(0, n_params, n_params)
  off <- 0L
  for (B in blocks) {
    b <- nrow(B)
    if (b == 0L) next
    idx <- off + seq_len(b)
    V[idx, idx] <- B
    off <- off + b
  }
  dd <- diag(V)
  dd[!is.finite(dd) | dd <= 0] <- 1e-8
  diag(V) <- dd
  V
}

clamp_w <- function(w) pmin(pmax(w, 0.001), 0.999)

occ_weights <- function(psi, p, N, n_valid, n_det, any_det) {
  weights <- numeric(N)
  for (i in seq_len(N)) {
    if (any_det[i]) { weights[i] <- 1 }
    else if (n_valid[i] == 0) { weights[i] <- psi[i] }
    else {
      prod_1mp <- (1 - p[i])^n_valid[i]
      num <- psi[i] * prod_1mp
      weights[i] <- num / (num + (1 - psi[i]))
    }
  }
  weights
}

glm_init <- function(X_occ, X_det, any_det, n_det, n_valid, keep, p_occ, p_det) {
  tryCatch({
    occ_glm <- glm(any_det ~ X_occ[, -1, drop = FALSE] - 1 + X_occ[, 1], family = binomial)
    det_glm <- glm(cbind(n_det[keep], n_valid[keep] - n_det[keep]) ~
                      X_det[keep, -1, drop = FALSE] - 1 + X_det[keep, 1], family = binomial)
    list(occ = list(beta = unname(coef(occ_glm)), se = rep(1, p_occ)),
         det = list(beta = unname(coef(det_glm)), se = rep(1, p_det)))
  }, error = function(e) {
    list(occ = list(beta = rep(0, p_occ), se = rep(1, p_occ)),
         det = list(beta = rep(0, p_det), se = rep(1, p_det)))
  })
}

# Build tobs_fit from EM result
build_laplace_fit <- function(em_result, model, spatial, p_per_submodel,
                              prior_spec = NULL,
                              approx = "gaussian_laplace",
                              re_block = NULL, latent_prior = NULL) {
  pi_list <- model$process_info

  # Per-arm SPDE membership: the field may sit on the state arm, the detection
  # arm, or both. The fixed-effect SE machinery is arm-specific (the Louis
  # observed-info correction assumes a non-spatial M-step Hessian on that arm),
  # so route each arm with its own field.
  spatial_occ <- .spatial_for_arm(spatial, 1L)
  spatial_det <- .spatial_for_arm(spatial, 2L)

  # Collect betas from correction (if available) or EM fits. Each arm also
  # contributes its within-arm covariance block (`prec_k` = the precision the
  # arm's SEs came from), so the pseudo-draws carry the joint correlation rather
  # than a diagonal stand-in (gcol33/tulpaObs#44). Cross-arm covariance stays
  # zero, matching the EM's separate-arm M-step factorization.
  means <- numeric()
  sds <- numeric()
  nms <- character()
  louis_psi_se <- NULL
  cov_blocks <- list()

  for (k in seq_along(pi_list)) {
    pi <- pi_list[[k]]
    sub_name <- names(p_per_submodel)[k]
    if (is.null(sub_name)) sub_name <- names(p_per_submodel)[min(k, length(p_per_submodel))]
    prec_k <- NULL

    if (is.list(em_result$pooled) && !is.null(em_result$pooled[[sub_name]])) {
      # MI/Gibbs correction pool from rubins_pool().
      cr <- em_result$pooled[[sub_name]]
      means <- c(means, cr$mean)
      sds_k <- cr$se
      cov_k <- if (!is.null(cr$vcov) &&
                   all(dim(as.matrix(cr$vcov)) == length(sds_k)))
                 as.matrix(cr$vcov) else .cor_scaled_cov(NULL, sds_k)
    } else if (!is.null(em_result$fits[[sub_name]])) {
      fi <- em_result$fits[[sub_name]]
      beta <- extract_beta(fi, pi$p)
      means <- c(means, beta)

      # Louis-corrected SE on the psi block of a single-season fit. The inner
      # M-step Hessian is M * I_complete + P_prior (pseudo-binomial trick);
      # I_obs = X' diag(psi(1-psi) - w(1-w)) X + P_prior is the right object
      # for SEs. See `.louis_info_psi_single` and tulpaObs#7.
      # Louis-corrected occupancy SE applies to the fixed-effect-only fit. When
      # random effects are present the occupancy block's fixed-effect SE comes
      # from the GLMM marginal precision (`H_beta`, Schur over the RE block)
      # that tulpa_laplace returns, so skip Louis on the RE path.
      use_louis <- identical(model$model_type, "single") &&
                   identical(sub_name, "occ") &&
                   !is.null(em_result$weights) &&
                   is.null(re_block)
      if (!is.null(re_block) && identical(sub_name, "occ")) {
        # Occupancy fixed-effect SE on the RE path: natural-scale observed info
        # marginalised over the random-effect block (the M-step H_beta is
        # M-inflated). Computed in .tobs_re_occ_fixed_se().
        sds_k <- re_block$occ_se
      } else if (use_louis) {
        I_obs <- .louis_info_psi_single(
          X_occ       = model$X_processes[[1]],
          beta_psi    = beta,
          weights     = em_result$weights,
          spatial     = spatial_occ,
          spatial_fit = fi,
          prior_spec  = prior_spec,
          coef_names  = pi$coef_names
        )
        louis_psi_se <- .se_from_info(I_obs, pi$p)
        sds_k <- louis_psi_se
        prec_k <- I_obs
      } else if (identical(model$model_type, "single") &&
                 identical(sub_name, "det") &&
                 is.null(spatial_det) &&
                 is.null(re_block) &&
                 is.null(model$X_det_visit) &&
                 !is.null(em_result$weights) &&
                 nrow(model$X_processes[[2]]) == length(em_result$weights)) {
        # Site-level detection SE via the marginal Louis observed info
        # (tulpaObs#7, detection arm). The detection M-step's H_beta is the
        # complete-data info (soft occupancy weight treated as known), which
        # under-states the SE the same way the psi arm did; recompute the
        # observed info, marginalizing over the coupled occupancy block.
        beta_psi_fit <- extract_beta(em_result$fits[["occ"]],
                                     ncol(model$X_processes[[1]]))
        I_obs <- .louis_info_det_single(
          X_occ          = model$X_processes[[1]],
          beta_psi       = beta_psi_fit,
          X_det          = model$X_processes[[2]],
          beta_det       = beta,
          weights        = em_result$weights,
          n_valid        = rowSums(model$y >= 0),
          prior_spec     = prior_spec,
          occ_coef_names = pi_list[[1]]$coef_names,
          det_coef_names = pi$coef_names,
          spatial        = spatial_occ,
          spatial_fit    = em_result$fits[["occ"]]
        )
        se_det <- .se_from_info(I_obs, pi$p)
        if (any(!is.finite(se_det))) {
          sds_k <- .se_from_laplace_fit(fi, pi$p)
          prec_k <- fi$H_beta
        } else {
          sds_k <- se_det
          prec_k <- I_obs
        }
      } else {
        sds_k <- .se_from_laplace_fit(fi, pi$p)
        prec_k <- fi$H_beta
      }
      cov_k <- .cor_scaled_cov(prec_k, sds_k)
    } else {
      means <- c(means, rep(0, pi$p))
      sds_k <- rep(NA_real_, pi$p)
      cov_k <- .cor_scaled_cov(NULL, sds_k)
    }
    sds <- c(sds, sds_k)
    cov_blocks[[length(cov_blocks) + 1L]] <- cov_k
    nms <- c(nms, paste0(pi$name, "_", pi$coef_names))
  }

  # Append visit-level detection coefficients when X_det_visit is present.
  # The detection M-step block has X of width p_det + p_det_visit; the main
  # loop above extracts only the first p_det elements (the site-level
  # detection coefs). Pull the visit-level tail and label as `p_visit_<name>`
  # so the public output matches the NUTS engine's column layout.
  if (!is.null(model$det_visit_names) && length(model$det_visit_names) > 0L) {
    p_det_visit <- length(model$det_visit_names)
    pi_p <- pi_list[[2]]  # detection process metadata
    p_det <- pi_p$p
    p_det_total <- as.integer(p_per_submodel[["det"]] %||% (p_det + p_det_visit))
    visit_idx <- (p_det + 1L):p_det_total
    visit_nms <- paste0("p_visit_", model$det_visit_names)

    if (is.list(em_result$pooled) && !is.null(em_result$pooled[["det"]])) {
      cr <- em_result$pooled[["det"]]
      visit_means <- cr$mean[visit_idx]
      visit_sds   <- cr$se[visit_idx]
    } else if (!is.null(em_result$fits[["det"]])) {
      fi_det <- em_result$fits[["det"]]
      beta_full <- extract_beta(fi_det, p_det_total)
      se_full <- .se_from_laplace_fit(fi_det, p_det_total)
      visit_means <- beta_full[visit_idx]
      visit_sds   <- se_full[visit_idx]
    } else {
      visit_means <- rep(0, p_det_visit)
      visit_sds   <- rep(NA_real_, p_det_visit)
    }

    means <- c(means, visit_means)
    sds   <- c(sds, visit_sds)
    nms   <- c(nms, visit_nms)
    # Visit-level detection coefs are carried diagonal (their cross-covariance
    # with the site-level det block is not surfaced separately from the M-step).
    cov_blocks[[length(cov_blocks) + 1L]] <- .cor_scaled_cov(NULL, visit_sds)
  }

  # Append the deterministic random-effect block (sigma hyperparameters +
  # per-group BLUPs) so the public output matches the NUTS column layout and
  # ranef() / summary() can name them (gcol33/tulpaObs#11).
  if (!is.null(re_block)) {
    means <- c(means, re_block$means)
    sds   <- c(sds, re_block$sds)
    nms   <- c(nms, re_block$names)
    cov_blocks[[length(cov_blocks) + 1L]] <- .cor_scaled_cov(NULL, re_block$sds)
  }

  names(means) <- nms
  names(sds)   <- nms
  n_params <- length(means)

  # Pseudo-draws from the block-diagonal joint covariance: full within each
  # fixed-effect arm (so derived quantities like predicted psi = plogis(X beta)
  # propagate the coefficient correlation), zero across arms. Coordinates with
  # an unavailable SE (NA) floor to a near-constant point mass, the same as the
  # previous per-coefficient draw (gcol33/tulpaObs#44).
  n_pseudo <- 1000L
  V_draw <- .assemble_block_diag(cov_blocks, n_params)
  dimnames(V_draw) <- list(nms, nms)
  draws <- .rmvn(n_pseudo, means, V_draw)
  colnames(draws) <- nms

  # Simplified-Laplace skewness correction
  # Computes gamma_j at the original observation likelihood (NOT the M-step
  # pseudo-binomial encoding — see dev_notes/simplified_laplace_derivation.md
  # §3 and dev_notes/upstream_tulpa_sla_spec.md §3 for why).
  sla_gamma <- NULL
  sla_status <- "off"
  # Simplified-Laplace skewness correction is not wired for the random-effect
  # path (the gamma derivation assumes a fixed-effect-only M-step).
  if (identical(approx, "simplified_laplace") && is.null(re_block)) {
    sla_res <- switch(model$model_type,
      single     = .sla_compute_occu_single(model, em_result,
                                            spatial = spatial,
                                            prior_spec = prior_spec),
      dynamic    = .sla_compute_dyn_occu(model, em_result,
                                         spatial = spatial,
                                         prior_spec = prior_spec),
      integrated = .sla_compute_int_occu(model, em_result,
                                         spatial = spatial,
                                         prior_spec = prior_spec),
      list(gamma = NULL, valid = FALSE,
           reason = sprintf("simplified Laplace not yet supported for model_type '%s'",
                            model$model_type))
    )
    if (isTRUE(sla_res$valid)) {
      sla_gamma  <- sla_res$gamma
      sla_status <- "simplified_laplace"
      # Align gamma names with the joint parameter ordering used in `means`
      sla_gamma <- sla_gamma[intersect(names(means), names(sla_gamma))]
      gamma_full <- setNames(rep(0, n_params), nms)
      gamma_full[names(sla_gamma)] <- sla_gamma
      sla_gamma <- gamma_full
      draws <- .sla_replace_draws(draws, means, sds, sla_gamma)
    } else {
      sla_status <- paste0("fallback_gaussian (", sla_res$reason, ")")
    }
  }

  intercepts <- compute_intercepts(model, means)

  # When SPDE is attached to the occ submodel, the M-step mode is
  # c(beta_occ, u_mesh). Extract u_mesh so callers can inspect or project
  # the latent field to observation locations via A %*% u_mesh.
  spatial_field <- NULL
  if (!is.null(spatial_occ) && !is.null(em_result$fits$occ$mode)) {
    p_occ <- pi_list[[1]]$p
    mode_vec <- em_result$fits$occ$mode
    if (length(mode_vec) > p_occ) {
      spatial_field <- mode_vec[(p_occ + 1L):length(mode_vec)]
    }
  }
  # Nested-Laplace: the occ block mode is c(beta_occ, latent units...) where the
  # latent tail is the grid-weighted posterior mean of the multi-block field.
  if (is.null(spatial_field) && !is.null(latent_prior) &&
      !is.null(em_result$fits$occ$mode)) {
    p_occ <- pi_list[[1]]$p
    mode_vec <- em_result$fits$occ$mode
    if (length(mode_vec) > p_occ) {
      spatial_field <- mode_vec[(p_occ + 1L):length(mode_vec)]
    }
  }

  # When SPDE is attached to the detection submodel, the det M-step mode is
  # c(beta_det, u_mesh_det); extract the detection field tail. The field is
  # identified off its own proper Matern (range, sigma) PC prior the same way
  # the state field is -- no separate sum-to-zero constraint is imposed (the
  # detection intercept absorbs the field level under the mean-zero prior).
  spatial_field_det <- NULL
  if (!is.null(spatial_det) && !is.null(em_result$fits$det$mode)) {
    p_det <- pi_list[[2]]$p
    mode_det <- em_result$fits$det$mode
    if (length(mode_det) > p_det) {
      spatial_field_det <- mode_det[(p_det + 1L):length(mode_det)]
    }
  }

  structure(c(list(
    draws = draws, means = means, sds = sds,
    skew = sla_gamma, sla_status = sla_status,
    n_samples = n_pseudo, n_params = n_params,
    log_prob = rep(NA_real_, n_pseudo)),
    .tobs_na_nuts_diagnostics(n_pseudo),
    list(
    col_names = nms, param_names = nms,
    intercepts = intercepts,
    model = model, spatial = spatial,
    spatial_field = spatial_field,
    spatial_field_det = spatial_field_det,
    process_info = model$process_info,
    method = "laplace",
    re_effects = re_block$re_effects,
    aghq = em_result$aghq,
    convergence = em_result$convergence,
    correction = em_result$correction
  )), class = c("tobs_fit", "tulpa_fit"))
}
