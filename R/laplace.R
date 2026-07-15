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
    count      = build_count_callbacks(model, spatial),
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

  # Debias step: the occupancy marginal is exact (single-season closed-form
  # two-state mixture, gcol33/tulpaObs#7; dynamic HMM forward, gcol33/tulpaObs#86),
  # so refine the EM mode with an exact-marginal Newton step and read calibrated
  # SEs from its Hessian. The EM's pseudo-binomial Laplace M-steps leave a small
  # discretisation residual below the marginal MLE and mis-scale the block SEs;
  # the refinement restores unbiased, near-nominal-coverage fixed effects
  # (matching unmarked::occu / colext). Spatial (nested-Laplace) fits carry a
  # latent field, not a closed-form coefficient marginal, and keep the EM result.
  if (is.null(spatial) && identical(approx, "gaussian_laplace") &&
      correction %in% c("auto", "none")) {
    if (identical(model$model_type, "single")) {
      fit <- .tobs_occu_marginal_refine(fit, model, prior_spec)
    } else if (identical(model$model_type, "dynamic")) {
      fit <- .tobs_dyn_occu_marginal_refine(fit, model, prior_spec)
    }
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

