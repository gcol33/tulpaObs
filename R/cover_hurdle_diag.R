# ---------------------------------------------------------------------------
# Pointwise log-likelihood (for WAIC / PSIS-LOO / tobs_stack)
# ---------------------------------------------------------------------------

# Pointwise log-likelihood [n_draws x N] for a cover hurdle fit. The hurdle is
# fully observed (nothing to marginalize): per site, an occurrence Bernoulli
# plus, when y > 0, the positive-part density. Per-arm betas are sampled from
# the Gaussian-Laplace posterior in the optimizer's (scaled) parameterization
# and paired with the stored scaled designs, so eta matches the fitted model
# without needing the raw data frame. Dispersion (sigma_pos / phi_pos) is held
# at its fitted value; SLA skew is not applied to the LOO marginals.
#
# Only the separate-Laplace path (method = "laplace" / "laplace_sla") carries
# the per-arm mode + Hessian this needs; the nested-joint path errors.
# Shared hurdle log-likelihood kernel: given per-draw linear predictors
# `eta_occ` [S x N] (occurrence) and `eta_pos` [S x N_pos] (positive part),
# dispersion `disp` (scalar or length-S), and the response encoding, accumulate
# the [S x N] pointwise hurdle log-likelihood -- log(1 - p) at absent sites,
# log(p) + positive-part log-density at occupied sites. The single source of
# truth for both the separate-Laplace and nested-joint cover paths and for the
# posterior-mean plug-in (S = 1).
.tobs_cover_hurdle_ll <- function(eta_occ, eta_pos, disp, occur, y_pos, idx_pos,
                                  positive, bounds = NULL) {
  S <- nrow(eta_occ); N <- ncol(eta_occ)
  sd_disp <- if (length(disp) == 1L) rep(disp, S) else disp
  log_p   <- .tobs_log_p(eta_occ)
  log_1mp <- .tobs_log_1mp(eta_occ)

  ll <- matrix(0, S, N)
  absent <- occur == 0L
  if (any(absent)) ll[, absent] <- log_1mp[, absent, drop = FALSE]

  pos_col <- match(seq_len(N), idx_pos)   # eta_pos column per site (NA if absent)
  for (i in which(occur == 1L)) {
    j <- pos_col[i]
    if (positive == "lognormal") {
      # y_pos = log(y); density of natural-scale y is the Gaussian on log(y)
      # times the Jacobian 1/y, i.e. dnorm(log y, eta, sigma, log) - log(y).
      dens <- stats::dnorm(y_pos[j], mean = eta_pos[, j], sd = sd_disp,
                           log = TRUE) - y_pos[j]
    } else if (positive == "gaussian") {
      # Identity-Gaussian arm (gcol33/tulpaObs#112): y_pos is the raw response,
      # so the density is the plain Gaussian with no change-of-variable Jacobian.
      dens <- stats::dnorm(y_pos[j], mean = eta_pos[, j], sd = sd_disp,
                           log = TRUE)
    } else if (positive == "lognormal_trunc") {
      # Upper-truncated lognormal: the lognormal density divided by the retained
      # mass Phi((u - eta)/sigma), u the log-cover ceiling (log(1) = 0). Adds
      # -log Phi((u - eta)/sigma) to the lognormal log-density above.
      u <- bounds$trunc_upper[j]
      dens <- stats::dnorm(y_pos[j], mean = eta_pos[, j], sd = sd_disp,
                           log = TRUE) - y_pos[j] -
              stats::pnorm((u - eta_pos[, j]) / sd_disp, log.p = TRUE)
    } else if (positive == "ordinal") {
      # Interval-censored Gaussian: the observed class is a probability MASS,
      # P(latent log-cover in (lower, upper]) = Phi((upper - eta)/sigma) -
      # Phi((lower - eta)/sigma), with +/-Inf the open outer classes (pnorm
      # handles them). A genuine PMF over classes -- no change-of-variable
      # Jacobian, so the score is measure-invariant and comparable across arms.
      zl <- (bounds$lower[j] - eta_pos[, j]) / sd_disp
      zu <- (bounds$upper[j] - eta_pos[, j]) / sd_disp
      dens <- log(pmax(stats::pnorm(zu) - stats::pnorm(zl), 1e-300))
    } else {
      mu   <- stats::plogis(eta_pos[, j])
      dens <- stats::dbeta(y_pos[j], mu * sd_disp, (1 - mu) * sd_disp,
                           log = TRUE)
    }
    ll[, i] <- log_p[, i] + dens
  }
  ll
}

# The per-plot positive-arm bounds, or NULL for beta / lognormal -- the
# pointwise-loglik and PIT consumers pass this to .tobs_cover_hurdle_ll. Ordinal
# carries the interval (lower, upper]; lognormal_trunc carries the truncation
# ceiling on the log-cover scale.
.tobs_cover_bounds <- function(object) {
  positive <- object$positive %||% "lognormal"
  enc <- object$encoding
  if (identical(positive, "ordinal")) {
    return(list(lower = enc$pos_data$lower, upper = enc$pos_data$upper))
  }
  if (identical(positive, "lognormal_trunc")) {
    return(list(trunc_upper = enc$pos_data$trunc_upper))
  }
  NULL
}

# Per-draw cover linear predictors [S x N] / [S x N_pos]. Separate-Laplace path:
# sample each arm's Gaussian-Laplace posterior at its scaled design. Nested-joint
# path: sample the grid-integrated joint and project the shared field at each
# observation's spatial unit (`spi_full` / `spi_pos`); the dispersion is then the
# per-draw grid value. Returns list(eta_occ, eta_pos, disp).
.tobs_cover_eta_draws <- function(object, n.draws = 1000L) {
  enc      <- object$encoding
  positive <- object$positive %||% "lognormal"
  if (!is.null(.tobs_joint_fit(object))) {
    bundle  <- .tobs_joint_draws(object, n = n.draws)
    if (isTRUE(object$armspecific)) {
      # Arm-specific separate latents store no node map / weight on the bundle
      # block (gcol33/tulpaObs#95); the pointwise-loglik consumer runs over the
      # fit's observations, so it rebuilds the per-arm per-observation node map
      # and covariate-weight lookup from `armspec_blocks` and hands them to
      # .tobs_joint_arm_eta exactly as the shared-field path passes spi_* / wfun.
      u_occ  <- .tobs_armspec_obs_units(object, 1L, nrow(enc$occ_data$X))
      u_pos  <- .tobs_armspec_obs_units(object, 2L, nrow(enc$pos_data$X))
      wf_occ <- .tobs_armspec_obs_wfun(object, 1L)
      wf_pos <- .tobs_armspec_obs_wfun(object, 2L)
      eta_occ <- t(.tobs_joint_arm_eta(bundle, enc$occ_data$X, "occ", u_occ, wf_occ))
      eta_pos <- t(.tobs_joint_arm_eta(bundle, enc$pos_data$X, "pos", u_pos, wf_pos))
      return(list(eta_occ = eta_occ, eta_pos = eta_pos, disp = bundle$disp))
    }
    spi_full <- object$spi_full
    spi_pos  <- object$spi_pos
    if (is.null(spi_full) || is.null(spi_pos)) {
      stop("Pointwise log-likelihood for the nested-joint cover() fit needs the ",
           "per-observation spatial-unit index (`spi_full` / `spi_pos`); refit ",
           "with the current tulpaObs so they are stored on the fit.",
           call. = FALSE)
    }
    # A coupled trend field carries a per-observation weight; resolve it per arm
    # via wfun (the engine's svc_weight replayed at predict time). NULL when the
    # fit has no trend field, in which case .tobs_joint_arm_eta never calls it.
    w_occ_fun <- if (!is.null(object$trend_w_occ))
                   function(col) as.numeric(object$trend_w_occ) else NULL
    w_pos_fun <- if (!is.null(object$trend_w_pos))
                   function(col) as.numeric(object$trend_w_pos) else NULL
    eta_occ <- t(.tobs_joint_arm_eta(bundle, enc$occ_data$X, "occ", spi_full,
                                     wfun = w_occ_fun))
    eta_pos <- t(.tobs_joint_arm_eta(bundle, enc$pos_data$X, "pos", spi_pos,
                                     wfun = w_pos_fun))
    return(list(eta_occ = eta_occ, eta_pos = eta_pos, disp = bundle$disp))
  }
  if (!is.null(object$nuts) && is.matrix(object$draws)) {
    # NUTS path: project the exact posterior coefficient draws through the
    # natural-scale presence / positive designs. The dispersion is the per-draw
    # exp(log_disp) trailing column, so the score marginalizes the calibrated
    # NUTS posterior (the point of the sampler) rather than the Laplace mode.
    draws  <- object$draws
    if (!is.null(n.draws) && n.draws < nrow(draws)) {
      draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
    }
    p_pres <- ncol(enc$occ_data$X); p_pos <- ncol(enc$pos_data$X)
    B_occ <- draws[, seq_len(p_pres), drop = FALSE]
    B_pos <- draws[, p_pres + seq_len(p_pos), drop = FALSE]
    disp  <- exp(draws[, ncol(draws)])
    return(list(eta_occ = B_occ %*% t(enc$occ_data$X),
                eta_pos = B_pos %*% t(enc$pos_data$X), disp = disp))
  }
  if (is.null(enc) || is.null(object$occ$mode) || is.null(object$occ$H_beta) ||
      is.null(object$pos$mode) || is.null(object$pos$H_beta)) {
    stop("Pointwise log-likelihood for cover() is implemented for the ",
         "separate-Laplace path (method = 'laplace' / 'laplace_sla') and the ",
         "nested-joint shared-field path (method = 'nested_laplace'); this fit ",
         "carries neither a per-arm mode + Hessian nor a joint object.",
         call. = FALSE)
  }
  X_occ <- enc$occ_data$X
  X_pos <- enc$pos_data$X
  p_occ <- ncol(X_occ); p_pos <- ncol(X_pos)
  mode_occ   <- object$occ$mode[seq_len(p_occ)]
  mode_pos   <- object$pos$mode[seq_len(p_pos)]
  pos_vscale <- if (positive %in% c("lognormal", "gaussian"))
                  (object$sigma_pos %||% 1)^2 else 1
  V_occ <- tryCatch(solve(object$occ$H_beta), error = function(e) NULL)
  V_pos <- tryCatch(pos_vscale * solve(object$pos$H_beta), error = function(e) NULL)
  S <- as.integer(n.draws)
  B_occ <- .tobs_mvn_draws(mode_occ, V_occ, S)   # [S x p_occ]
  B_pos <- .tobs_mvn_draws(mode_pos, V_pos, S)   # [S x p_pos]
  disp  <- if (positive %in% c("lognormal", "gaussian")) object$sigma_pos
           else object$phi_pos
  list(eta_occ = B_occ %*% t(X_occ), eta_pos = B_pos %*% t(X_pos), disp = disp)
}

# Pointwise log-likelihood [n_draws x N] for a cover hurdle fit (separate-Laplace
# or nested-joint shared-field), assembled from the per-draw linear predictors
# through the shared hurdle kernel.
.tobs_ploglik_cover <- function(object, n.draws = 1000L, n.threads = 1L) {
  enc <- object$encoding
  e   <- .tobs_cover_eta_draws(object, n.draws)
  .cover_hurdle_ploglik_core(e$eta_occ, e$eta_pos, e$disp, enc$occ_data$y,
                             enc$pos_data$y, enc$idx_pos,
                             object$positive %||% "lognormal",
                             bounds = .tobs_cover_bounds(object),
                             n_threads = n.threads)
}

# Parallel C++ pointwise log-likelihood for the cover() hurdle, over draws
# (cpp_cover_hurdle_ploglik). Mirrors .tobs_cover_hurdle_ll (the R oracle, kept
# for the posterior-mean plug-in and the tests); every positive family routes
# through the one kernel. The absent/present split and the four family densities
# match the R kernel, so the two agree to libm rounding.
.cover_hurdle_ploglik_core <- function(eta_occ, eta_pos, disp, occur, y_pos,
                                       idx_pos, positive, bounds = NULL,
                                       n_threads = 1L) {
  S <- nrow(eta_occ); N <- ncol(eta_occ)
  disp_full <- if (length(disp) == 1L) rep(as.numeric(disp), S) else as.numeric(disp)
  pos_col <- match(seq_len(N), idx_pos)
  pos_col[is.na(pos_col)] <- 0L
  code <- switch(positive,
                 lognormal = 0L, lognormal_trunc = 1L, ordinal = 2L, beta = 3L,
                 gaussian = 4L,
                 stop("cover pointwise loglik: unknown positive family '",
                      positive, "'.", call. = FALSE))
  num <- function(x) if (is.null(x)) numeric(0) else as.numeric(x)
  cpp_cover_hurdle_ploglik(
    eta_occ, eta_pos, disp_full, as.integer(occur), as.numeric(y_pos),
    as.integer(pos_col), code,
    num(bounds$lower), num(bounds$upper), num(bounds$trunc_upper),
    max(1L, as.integer(n_threads)))
}

# Pointwise log-likelihood at the posterior mean of the parameters (length N):
# the hurdle kernel evaluated at the mean linear predictors and mean dispersion.
.tobs_cover_loglik_at_mean <- function(object, n.draws = 1000L) {
  enc <- object$encoding
  e   <- .tobs_cover_eta_draws(object, n.draws)
  mean_eta_occ <- matrix(colMeans(e$eta_occ), nrow = 1L)
  mean_eta_pos <- matrix(colMeans(e$eta_pos), nrow = 1L)
  as.numeric(.tobs_cover_hurdle_ll(
    mean_eta_occ, mean_eta_pos, mean(e$disp), enc$occ_data$y,
    enc$pos_data$y, enc$idx_pos, object$positive %||% "lognormal",
    bounds = .tobs_cover_bounds(object)
  ))
}


# ---------------------------------------------------------------------------
# Posterior predictive check + PIT for the cover hurdle (gcol33/tulpaObs#27)
# ---------------------------------------------------------------------------

# Randomized PIT for a cover() hurdle fit (length N). The predictive CDF mixes a
# point mass 1 - p at the structural zero with p * F_pos on the positive part
# (lognormal / beta CDF at the fitted per-draw predictor), projecting the shared
# field per observation for the nested-joint fit. Absent sites use the left /
# right limits [0, 1 - p] around the zero mass; occupied sites are continuous
# (F = 1 - p + p F_pos), so the engine's randomized PIT is degenerate there.
.tobs_pit_cover <- function(object, n.samples = 250) {
  enc      <- object$encoding
  positive <- object$positive %||% "lognormal"
  bounds   <- .tobs_cover_bounds(object)
  e <- .tobs_cover_eta_draws(object, n.draws = n.samples)
  S <- nrow(e$eta_occ); N <- ncol(e$eta_occ)
  sd_disp <- if (length(e$disp) == 1L) rep(e$disp, S) else e$disp
  pos_col <- match(seq_len(N), enc$idx_pos); pos_col[is.na(pos_col)] <- 0L
  code <- switch(positive, lognormal = 0L, lognormal_trunc = 1L, ordinal = 2L,
                 beta = 3L, gaussian = 4L,
                 stop("cover PIT: unknown positive family '", positive, "'.",
                      call. = FALSE))
  num <- function(x) if (is.null(x)) numeric(0) else as.numeric(x)
  # The per-observation predictive-CDF limits are deterministic; the former R
  # loop over occupied plots now runs in cpp_cover_pit_cdf.
  lim <- cpp_cover_pit_cdf(e$eta_occ, e$eta_pos, as.integer(enc$occ_data$y),
                           as.numeric(enc$pos_data$y), as.integer(pos_col),
                           sd_disp, code, num(bounds$lower), num(bounds$upper),
                           num(bounds$trunc_upper))
  tulpa::tulpa_pit(lim$cdf_upper, cdf_lower = lim$cdf_lower)
}

# Posterior predictive check for a cover() hurdle fit. Per draw, occurrence
# replicates from Bernoulli(p) and cover replicates from the fitted positive
# part; the discrepancy is a Freeman-Tukey (or chi-squared) sum on the
# occurrence arm (all N sites, against p) plus the positive arm (occupied
# subset, against the positive-part mean), returning a Bayesian p-value.
.tobs_ppc_cover <- function(object,
                            fit.stat = c("freeman-tukey", "chi-squared"),
                            n.samples = 500) {
  fit.stat <- match.arg(fit.stat)
  enc      <- object$encoding
  positive <- object$positive %||% "lognormal"
  if (identical(positive, "gaussian")) {
    stop("tobs_ppc() is not defined for cover(response = \"gaussian\"): the ",
         "Freeman-Tukey and chi-squared discrepancies assume a non-negative ",
         "response, but the identity-Gaussian arm is unbounded. Use tobs_waic() ",
         "/ tobs_loo() for model comparison.", call. = FALSE)
  }
  e <- .tobs_cover_eta_draws(object, n.draws = n.samples)
  S <- nrow(e$eta_occ)
  y_pos <- enc$pos_data$y
  sd_disp <- if (length(e$disp) == 1L) rep(e$disp, S) else e$disp
  # Ordinal / lognormal store per-plot log-cover; exp() recovers natural scale.
  y_pos_nat <- if (positive %in% c("lognormal", "lognormal_trunc", "ordinal"))
                 exp(y_pos) else y_pos
  trunc_u <- if (positive == "lognormal_trunc") enc$pos_data$trunc_upper else numeric(0)
  code <- switch(positive, lognormal = 0L, lognormal_trunc = 1L, ordinal = 2L,
                 beta = 3L, gaussian = 4L,
                 stop("cover PPC: unknown positive family '", positive, "'.",
                      call. = FALSE))
  # The occurrence + cover replicates draw from R's RNG stream in the C++ kernel
  # in the same order as the former R loop, so under a fixed seed the discrepancy
  # is byte-identical.
  r <- cpp_cover_ppc(e$eta_occ, e$eta_pos, as.integer(enc$occ_data$y),
                     as.numeric(y_pos_nat), sd_disp, as.numeric(trunc_u), code,
                     identical(fit.stat, "freeman-tukey"))
  list(fit.y = r$fit.y, fit.y.rep = r$fit.y.rep,
       bayesian.p = mean(r$fit.y.rep > r$fit.y))
}
# mode (point mass) if V is unavailable, and jitters a near-singular V.
.tobs_mvn_draws <- function(mu, V, S) {
  p <- length(mu)
  if (is.null(V)) return(matrix(mu, S, p, byrow = TRUE))
  R <- tryCatch(chol(V), error = function(e) {
    jit <- 1e-8 * (mean(diag(V)) + 1e-8)
    chol(V + diag(jit, p))
  })
  Z <- matrix(stats::rnorm(S * p), S, p)
  matrix(mu, S, p, byrow = TRUE) + Z %*% R
}


# ---------------------------------------------------------------------------
# Predict
# ---------------------------------------------------------------------------

#' Predict cover from a cover_fit
#'
#' Occurrence probability is always `p = plogis(X * beta_occ)`. The
#' conditional positive cover `mu` depends on the positive-part family:
#'
#' * `positive = "lognormal"`: `mu = exp(eta_pos + sigma_pos^2 / 2)`
#'   (lognormal back-transform on log-cover).
#' * `positive = "beta"`: `mu = plogis(eta_pos)` (mean of the beta on
#'   the natural cover scale with logit link).
#'
#' Expected cover is `E[y] = p * mu` under both positive parts.
#'
#' The separate-Laplace fit (`method = "laplace"`) returns a fixed-effects-only
#' numeric vector. The nested-Laplace shared-field fit
#' (`method = "nested_laplace"`) instead projects the shared occupancy-cover
#' field and returns a `tobs_prediction` of posterior draws -- the same tidy /
#' `change` contract as [predict.tobs_fit()] for `occu_cover()`: pass
#' `type = "change"` with `times = c(t1, t2)` and `time_col` for a per-cell
#' delta map. Each prediction unit is a row of `newdata` (or a `cell` column
#' indexing the field cells), and every quantity is marginalized per draw over
#' the grid-integrated joint posterior.
#'
#' @param object A `cover_fit`.
#' @param newdata A data frame of covariates matching the original formula. For
#'   the nested-Laplace fit, one row per spatial unit (or a `cell` column).
#' @param type Separate-Laplace fit: one of `"expected"`, `"occupancy"`,
#'   `"conditional"`. Nested-Laplace fit: `"occurrence"`, `"cover_cond"`,
#'   `"cover_exp"`, or `"change"` (the legacy aliases are accepted and mapped).
#' @param include_RE Ignored for the separate-Laplace fit (no spatial
#'   projection); the nested-Laplace fit always projects the shared field.
#' @param times,time_col,level,nsim,draws Nested-Laplace fit only: `times =
#'   c(t1, t2)` and `time_col` drive the `"change"` map; `level` is the credible
#'   level, `nsim` the draw count, `draws` whether to attach the draw matrices.
#' @param ... Unused.
#' @return Separate-Laplace fit: a numeric vector. Nested-Laplace fit: a
#'   `tobs_prediction`.
#' @export
predict.cover_fit <- function(object, newdata = NULL,
                                     type = NULL, include_RE = FALSE,
                                     times = NULL, time_col = NULL,
                                     level = 0.95, nsim = 1000L, draws = TRUE,
                                     ...) {
  # Nested-Laplace shared-field fit: route through the unified joint predict
  # substrate (gcol33/tulpaObs#23). Map the legacy fixed-effects type names onto
  # the joint vocabulary so old calls keep working.
  if (!is.null(.tobs_joint_fit(object))) {
    if (is.null(type)) type <- "occurrence"
    type <- switch(type,
                   expected    = "cover_exp",
                   occupancy   = "occurrence",
                   conditional = "cover_cond",
                   type)
    return(.tobs_predict_joint(object, newdata = newdata, type = type,
                               times = times, level = level, nsim = nsim,
                               draws = draws, time_col = time_col))
  }

  if (is.null(type)) type <- "expected"
  type <- match.arg(type, c("expected", "occupancy", "conditional"))
  if (missing(newdata) || is.null(newdata)) {
    stop("`newdata` is required.", call. = FALSE)
  }
  if (isTRUE(include_RE) && !is.null(object$encoding$spatial_spec)) {
    message("predict.cover_fit(): spatial RE projection at new ",
            "locations is not implemented for the separate-Laplace fit; ",
            "returning fixed-effects-only predictions.")
  }

  X <- stats::model.matrix(object$encoding$formula, newdata)
  if (ncol(X) != length(object$beta_occ)) {
    stop("Design-matrix column count (", ncol(X), ") does not match the ",
         "fitted model (", length(object$beta_occ), "). Check `newdata`.",
         call. = FALSE)
  }

  eta_occ <- as.numeric(X %*% object$beta_occ)
  eta_pos <- as.numeric(X %*% object$beta_pos)
  p  <- stats::plogis(eta_occ)
  positive <- object$positive %||% "lognormal"
  mu <- if (positive %in% c("beta", "beta_oi")) {
    stats::plogis(eta_pos)
  } else if (identical(positive, "gaussian")) {
    # Identity-Gaussian arm (gcol33/tulpaObs#112): mu = eta on the response scale.
    eta_pos
  } else {
    exp(eta_pos + object$sigma_pos^2 / 2)
  }
  # One-inflated Beta: conditional cover mixes the ceiling mass (cover = 1) with
  # the interior Beta mean, E[cover | y > 0] = pi + (1 - pi) * mu.
  pi1  <- if (identical(positive, "beta_oi")) (object$pi_one %||% 0) else 0
  cond <- pi1 + (1 - pi1) * mu

  switch(
    type,
    expected    = p * cond,
    occupancy   = p,
    conditional = cond
  )
}


# ---------------------------------------------------------------------------
# Print / summary
# ---------------------------------------------------------------------------

#' @export
print.cover_fit <- function(x, ...) {
  positive <- x$positive %||% "lognormal"
  cat(sprintf("<cover_fit (%s positive part)>\n", positive))
  cat(sprintf("  N total      : %d\n", x$n_total))
  cat(sprintf("  N positive   : %d (%.1f%%)\n",
              x$n_positive, 100 * x$n_positive / x$n_total))
  if (positive %in% c("lognormal", "lognormal_trunc", "ordinal", "gaussian")) {
    cat(sprintf("  sigma_pos    : %.4f\n", x$sigma_pos))
  } else {
    cat(sprintf("  phi_pos      : %.4f\n", x$phi_pos))
  }
  if (identical(positive, "beta_oi") && !is.null(x$pi_one) && is.finite(x$pi_one)) {
    cat(sprintf("  pi_one       : %.4f (SE %.4f; %d of %d positive plots at cover = 1)\n",
                x$pi_one, x$pi_one_sd %||% NA_real_, x$n_ceiling %||% NA_integer_,
                x$n_positive %||% NA_integer_))
  }
  cat(sprintf("  converged    : %s\n",
              if (isTRUE(x$converged)) "yes" else "no"))
  if (!is.null(x$sla_status) && !identical(x$sla_status, "off")) {
    cat(sprintf("  marginals    : %s\n", x$sla_status))
  }
  cat("\nPresence (binomial logit):\n")
  print(.coef_table(x$beta_occ, x$se_occ))
  cat("\n", .cover_pos_header(positive), "\n", sep = "")
  print(.coef_table(x$beta_pos, x$se_pos))
  invisible(x)
}

#' @export
summary.cover_fit <- function(object, ...) {
  # NUTS fit: return the per-parameter posterior table (mean / sd / quantiles
  # plus the cross-chain Rhat / ESS the convergence list carries), matching the
  # generic NUTS summary surface so the sampler diagnostics are visible.
  if (!is.null(object$nuts) && !is.null(object$draws)) {
    return(.tobs_cover_nuts_summary(object))
  }
  # Arm labels (gcol33/tulpaObs#61): the two hurdle arms are `presence`
  # (the y > 0 Bernoulli arm) and `positive` (the y | y > 0 arm). The `to =`
  # argument of a spatial() bar validates against these labels, so summary()
  # prints the same names (formula label == output label).
  out <- list(
    family       = object$family,
    positive     = object$positive %||% "lognormal",
    n_total      = object$n_total,
    n_positive   = object$n_positive,
    sigma_pos    = object$sigma_pos,
    phi_pos      = object$phi_pos,
    converged    = object$converged,
    presence     = .coef_table(object$beta_occ, object$se_occ),
    positive_arm = .coef_table(object$beta_pos, object$se_pos),
    log_marginal = object$log_marginal,
    hyperpar     = object$hyperpar
  )
  class(out) <- "summary.cover_fit"
  out
}

#' @export
print.summary.cover_fit <- function(x, ...) {
  cat("Cover hurdle fit summary\n")
  cat(sprintf("  positive part: %s\n", x$positive))
  cat(sprintf("  N total = %d, N positive = %d\n", x$n_total, x$n_positive))
  if (x$positive %in% c("lognormal", "lognormal_trunc", "ordinal", "gaussian")) {
    cat(sprintf("  sigma_pos = %.4f\n", x$sigma_pos))
  } else {
    cat(sprintf("  phi_pos   = %.4f\n", x$phi_pos))
  }
  cat(sprintf("  log marginal: occ = %.3f, pos = %.3f\n",
              x$log_marginal["occ"], x$log_marginal["pos"]))
  cat("\nPresence:\n"); print(x$presence)
  cat("\n", .cover_pos_header(x$positive), "\n", sep = "")
  print(x$positive_arm)
  invisible(x)
}

# Header naming the positive-arm density and the support it is fit on. Shared by
# print.cover_fit() and print.summary.cover_fit() so a new `positive` family is
# labelled once. The unmatched default is the lognormal arm (Gaussian on log y).
.cover_pos_header <- function(positive) {
  switch(positive,
    beta     = "Positive (beta, logit link, on y > 0):",
    beta_oi  = "Positive interior (one-inflated beta, logit link, on 0 < y < 1):",
    ordinal  = "Positive (ordinal interval-censored Gaussian on log y > 0):",
    gaussian = "Positive (identity-Gaussian on y != 0):",
    "Positive (Gaussian on log y > 0):")
}


