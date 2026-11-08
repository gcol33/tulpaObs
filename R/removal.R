# =============================================================================
# removal.R — removal-sampling abundance family (sequential depletion)
#
# Latent abundance N_i ~ Poisson(lambda_i) (or NegBin), with K ordered removal
# passes. Pass k removes y_{ik} of the A_{ik} = N_i - sum_{l<k} y_{il}
# individuals still present, each with probability p_{ik}. The depleting-binomial
# product equals the multinomial-removal likelihood (Royle 2004); the latent N is
# summed out in closed form (sum to K_max), so -- as for the static N-mixture --
# there is no EM: the package-internal `removal_laplace()` fits the marginal
# directly with analytic gradients and observed-Fisher curvature. This file owns
# the family interface; the per-site math is `src/removal_kernel.h` and the
# (shared) Laplace driver is `src/marginal_count_laplace.h`.
#
#   .tobs_build_removal()  data binder -> model_type = "removal"
#   .tobs_fit_removal()    dispatch to the removal Laplace fit
#   removal_laplace()      R wrapper over cpp_removal_laplace_fixed
# =============================================================================


# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# Bind a removal-sampling model. The abundance predictor is site-level
# (X_lambda, n_sites rows); the detection predictor is per-pass long form (one
# row per pass), so site-level detection covariates are replicated across a
# site's passes and pass-level covariates (X_det_visit, site-major) are stacked
# on. Unlike the N-mixture, removal needs COMPLETE pass sequences (depletion
# accumulates earlier passes), so an NA in `y` is an error rather than a dropped
# visit, and the passes are taken in column order (pass 1 = column 1).
.tobs_build_removal <- function(abund_formula, det_formula, data, y,
                                det_visit_formula = NULL, det_visit_data = NULL) {
  if (!is.matrix(y)) {
    stop("y must be a matrix (n_sites x n_passes) of integer removal counts.",
         call. = FALSE)
  }
  .tobs_check_site_count(nrow(y), nrow(data), "rows")
  if (anyNA(y)) {
    stop("removal() needs complete pass sequences: y must not contain NA. ",
         "Depletion accumulates the removals of earlier passes, so a missing ",
         "pass cannot be dropped the way an N-mixture visit can.", call. = FALSE)
  }
  y_int <- matrix(as.integer(round(y)), nrow(y), ncol(y))
  if (any(y_int < 0L)) {
    stop("y must contain nonnegative integer removal counts.", call. = FALSE)
  }
  # Each site's passes reach the kernel in depletion order with no gaps.
  .tobs_build_count_model("removal", abund_formula, det_formula, data, y_int,
                          det_visit_formula, det_visit_data, arm = "removal pass",
                          extra = list(n_passes = ncol(y_int)))
}


# ---------------------------------------------------------------------------
# Fitter
# ---------------------------------------------------------------------------

# Called from `.tobs_fit_model()` for `model$model_type == "removal"`, after the
# per-process design autoscaling. Wraps the removal Laplace fit into a tobs_fit
# (reusing build_nmix_fit -- the coefficient layout, joint vcov, draws, and NB
# dispersion handling are the (lambda, p[, log_r]) shape both families share).
.tobs_fit_removal <- function(model, mixture = "poisson", K_max = NULL,
                              max_iter = 100L, tol = 1e-6, verbose = TRUE) {
  mix_code <- switch(mixture,
                     poisson = "P", negbin = "NB",
                     stop(sprintf("Unknown mixture '%s' (use \"poisson\" or \"negbin\").",
                                  mixture), call. = FALSE))
  raw <- removal_laplace(
    y         = model$y_long,
    site_idx  = model$site_idx,
    X_lambda  = model$X_processes[[1]],
    X_p       = model$X_processes[[2]],
    mixture   = mix_code,
    K_max     = K_max,
    max_iter  = as.integer(max_iter),
    tol       = as.numeric(tol),
    verbose   = isTRUE(verbose)
  )
  build_nmix_fit(raw, model, spatial = NULL)
}


# No-RE removal Laplace warm start for the AGHQ RE path (the shape
# .tobs_fit_count_re's warm_fun expects: (model, mixture, K_max, max_iter, tol)).
# `mixture` is the "P"/"NB" code.
.tobs_removal_re_warm <- function(model, mixture, K_max, max_iter, tol) {
  removal_laplace(y = model$y_long, site_idx = model$site_idx,
                  X_lambda = model$X_processes[[1]], X_p = model$X_processes[[2]],
                  mixture = mixture, K_max = K_max,
                  max_iter = as.integer(max_iter), tol = as.numeric(tol),
                  verbose = FALSE)
}

# Fit a removal model with a site-level grouped random effect on the abundance
# (lambda) OR detection (p) arm. Routes through the shared
# count-model RE path: the per-site removal marginal carries the same
# (lambda, p[, log_r]) coefficient layout and REGroupOracle interface as the
# N-mixture, so the only family-specific pieces are the removal warm start and
# the removal AGHQ helper (its native RemovalGroupedOracle + depletion-aware
# K_max default).
.tobs_fit_removal_re <- function(model, re, mixture = "poisson", K_max = NULL,
                                 max_iter = 100L, tol = 1e-6, verbose = TRUE,
                                 n_quad = 1L, lkj_eta = 1.5,
                                 theta_prior_sd = 100) {
  mix_code <- switch(mixture,
                     poisson = "P", negbin = "NB", P = "P", NB = "NB",
                     stop(sprintf("Unknown mixture '%s' (use \"poisson\" or \"negbin\").",
                                  mixture), call. = FALSE))
  .tobs_fit_count_re(model, re,
                     warm_fun = .tobs_removal_re_warm,
                     aghq_fun = .tobs_removal_re_aghq,
                     family_label = "removal",
                     mixture = mix_code, K_max = K_max,
                     max_iter = max_iter, tol = tol, verbose = verbose,
                     n_quad = n_quad, lkj_eta = lkj_eta,
                     theta_prior_sd = theta_prior_sd)
}


# ---------------------------------------------------------------------------
# Laplace fit (R wrapper over cpp_removal_laplace_fixed)
# ---------------------------------------------------------------------------

#' Laplace fit of the removal-sampling abundance model
#'
#' @description
#' Maximum-likelihood fit (non-spatial, fixed effects only) of the removal /
#' sequential-depletion abundance model with a Poisson or negative-binomial
#' abundance mixing distribution. Pass `k` removes
#' `y_{ik} ~ Binomial(N_i - sum_{l<k} y_{il}, p_{ik})` of the individuals still
#' present, and the latent `N_i ~ Poisson(lambda_i)` (or `NegBin(lambda_i, r)`)
#' is summed out in closed form (truncation `K_max`). Optimisation is the same
#' analytic inner-Newton / observed-Fisher path the N-mixture uses; the per-site
#' detection arm sees the depleting available count rather than the full `N`.
#'
#' @param y Integer vector of per-pass removals (long form), passes in order.
#' @param site_idx Integer vector, same length as `y`, 1-based site index; a
#'   site's passes must appear in pass order (the depletion order).
#' @param X_lambda Numeric matrix `[n_sites x p_lambda]` of abundance covariates.
#' @param X_p Numeric matrix `[n_obs x p_p]` of per-pass detection covariates.
#' @param mixture `"P"` (Poisson, default) or `"NB"` (negative binomial).
#' @param beta_lambda_init,beta_p_init,log_r_init Optional warm starts.
#' @param r_max Upper bound on the NB size `r` (NB only, default `1e5`).
#' @param K_max Marginal-sum truncation; defaults to `max(site removal total) +
#'   100`. Must be `>=` the largest per-site removal total.
#' @param max_iter Newton iteration budget (default 100).
#' @param tol Gradient-norm convergence tolerance (default 1e-6).
#' @param verbose Print per-iteration progress.
#'
#' @return A list of class `nmix_fit` (the shared count-marginal fit shape):
#'   `beta_lambda`, `beta_p`, `mixture`, `log_r` / `r` (NB), `log_lik`, `vcov`,
#'   `H_obs`, `n_iter`, `converged`, `grad_norm`, per-site `mean_N` / `var_N` /
#'   `boundary_weight`.
#'
#' @references
#' Royle, J. A. (2004). N-mixture models for estimating population size from
#'   spatially replicated counts. *Biometrics* 60, 108-115.
#' Dorazio, R. M., Jelks, H. L., Jordan, F. (2005). Improving removal-based
#'   estimates of abundance. *Biometrics* 61, 1093-1101.
#'
#' @keywords internal
removal_laplace <- function(y, site_idx, X_lambda, X_p,
                            mixture = c("P", "NB"),
                            beta_lambda_init = NULL, beta_p_init = NULL,
                            log_r_init = NULL, r_max = 1e5, K_max = NULL,
                            max_iter = 100L, tol = 1e-6, verbose = FALSE) {
  mixture <- match.arg(mixture)
  nb <- identical(mixture, "NB")
  y        <- as.integer(y)
  site_idx <- as.integer(site_idx)
  if (!is.matrix(X_lambda)) stop("`X_lambda` must be a numeric matrix.", call. = FALSE)
  if (!is.matrix(X_p))      stop("`X_p` must be a numeric matrix.", call. = FALSE)
  if (length(y) != length(site_idx)) stop("length(y) must equal length(site_idx).", call. = FALSE)
  if (length(y) != nrow(X_p))        stop("length(y) must equal nrow(X_p).", call. = FALSE)
  if (any(y < 0L) || anyNA(y)) stop("`y` must be nonnegative integers with no NA.", call. = FALSE)
  n_sites <- nrow(X_lambda)
  if (min(site_idx) < 1L || max(site_idx) > n_sites) {
    stop("`site_idx` values must lie in [1, nrow(X_lambda)].", call. = FALSE)
  }

  p_lambda <- ncol(X_lambda); p_p <- ncol(X_p)
  site_tot <- tapply(y, factor(site_idx, levels = seq_len(n_sites)), sum)
  site_tot[is.na(site_tot)] <- 0L
  R_max <- max(as.integer(site_tot))
  if (is.null(beta_lambda_init)) {
    beta_lambda_init <- c(log(max(mean(as.numeric(site_tot)), 0.1) + 0.1),
                          rep(0, p_lambda - 1L))
  }
  if (is.null(beta_p_init)) beta_p_init <- rep(0, p_p)
  if (length(beta_lambda_init) != p_lambda) stop("length(beta_lambda_init) must equal ncol(X_lambda).", call. = FALSE)
  if (length(beta_p_init) != p_p)           stop("length(beta_p_init) must equal ncol(X_p).", call. = FALSE)

  if (is.null(K_max)) {
    K_max <- as.integer(R_max + 100L)
  } else {
    K_max <- as.integer(K_max)
    if (K_max < R_max) stop("`K.max` must be >= the largest per-site removal total.", call. = FALSE)
  }
  if (!is.numeric(r_max) || length(r_max) != 1L || r_max <= 0) {
    stop("`r_max` must be a positive scalar.", call. = FALSE)
  }
  if (is.null(log_r_init)) log_r_init <- 0
  else if (length(log_r_init) != 1L || !is.finite(log_r_init)) {
    stop("`log_r_init` must be a finite scalar.", call. = FALSE)
  }

  fit <- cpp_removal_laplace_fixed(
    y = y, site_idx = site_idx, X_lambda_R = X_lambda, X_p_R = X_p,
    beta_lambda_init = as.numeric(beta_lambda_init),
    beta_p_init = as.numeric(beta_p_init),
    K_max = K_max, max_iter = as.integer(max_iter), tol = as.numeric(tol),
    verbose = isTRUE(verbose), nb = nb,
    log_r_init = as.numeric(log_r_init), theta_max = log(r_max))

  nm_lam <- colnames(X_lambda); nm_p <- colnames(X_p)
  if (is.null(nm_lam)) nm_lam <- paste0("lam_", seq_len(p_lambda))
  if (is.null(nm_p))   nm_p   <- paste0("p_", seq_len(p_p))
  names(fit$beta_lambda) <- nm_lam
  names(fit$beta_p)      <- nm_p
  coef_names <- c(nm_lam, nm_p, if (nb) "log_r")
  rownames(fit$vcov)  <- colnames(fit$vcov)  <- coef_names
  rownames(fit$H_obs) <- colnames(fit$H_obs) <- coef_names

  fit$mixture <- mixture
  if (!nb) { fit$log_r <- NA_real_; fit$r <- NA_real_ }
  fit$K_max <- K_max
  fit$n_sites <- n_sites
  fit$n_obs <- length(y)
  fit$call <- match.call()
  if (nb && isTRUE(fit$dispersion_boundary)) {
    warning(sprintf(paste0("NB dispersion pinned at the boundary (r = r_max = %.3g); the data ",
            "are consistent with Poisson. Consider mixture = \"P\"."), r_max), call. = FALSE)
  }
  max_bw <- max(fit$boundary_weight, na.rm = TRUE)
  if (is.finite(max_bw) && max_bw > 1e-4) {
    warning(sprintf("Max posterior weight on N = K.max is %.2e at %d sites; raise K.max.",
                    max_bw, sum(fit$boundary_weight > 1e-4)), call. = FALSE)
  }
  class(fit) <- c("nmix_fit", "list")
  fit
}


# ---------------------------------------------------------------------------
# S3 helpers (routed from methods.R by model_type == "removal")
# ---------------------------------------------------------------------------

# Per-pass removal-cell probabilities pi_k = p_k prod_{l<k}(1 - p_l) for a vector
# of per-pass detection probabilities `p` in pass order.
.removal_pi <- function(p) {
  K <- length(p)
  if (K == 0L) return(numeric(0))
  surv <- c(1, cumprod(1 - p)[-K])     # prod_{l<k}(1 - p_l)
  p * surv
}

# The removal() generator, `cpp_simulate_removal()`: latent N per site, then
# each site's passes in turn, every pass removing a binomial share of the
# individuals not yet removed. theta = list(beta, off, r): the coefficient
# vector in the fit's draw layout, each arm's realised offset (or NULL), and
# the NB size (NA for Poisson).
.tobs_gen_removal <- function(model, theta) {
  p <- .tobs_arm_p(model, 2L)
  r <- theta$r %||% NA_real_
  cpp_simulate_removal(model$X_processes[[1L]], model$X_processes[[2L]],
    as.numeric(theta$beta), .tobs_off_arg(theta$off[[1L]]),
    .tobs_off_arg(theta$off[[2L]]),
    as.integer(model$site_idx), as.integer(model$visit_idx),
    model$n_sites, model$n_passes, p[1L], p[2L], is.finite(r), as.numeric(r))
}

.tobs_simulate_removal <- function(object, nsim = 1) {
  r_size <- object$nmix_dispersion$r
  .tobs_simulate_draw_gen(object, nsim, .tobs_gen_removal, 2L,
    list(r = if (!is.null(r_size) && is.finite(r_size)) as.numeric(r_size)
             else NA_real_))
}

# residuals() for removal. Under Poisson abundance the pass-k removal is
# marginally y_ik ~ Poisson(lambda_i * pi_ik), pi_ik the removal-cell prob; score
# the fit on that mean with Poisson deviance / Pearson residuals.
.tobs_residuals_removal <- function(object, type = c("deviance", "pearson",
                                                    "response")) {
  type  <- match.arg(type)
  model <- object$model
  fitv  <- .tobs_fitted_nmix(object)
  lambda <- fitv$lambda
  p_obs  <- fitv$p
  site_idx  <- model$site_idx
  visit_idx <- model$visit_idx
  y_long    <- model$y_long
  n_sites   <- model$n_sites
  n_pass    <- model$n_passes

  p_mat <- matrix(NA_real_, n_sites, n_pass)
  p_mat[cbind(site_idx, visit_idx)] <- p_obs
  pi_mat <- t(apply(p_mat, 1, .removal_pi))
  mu_mat <- lambda * pi_mat
  mu_long <- mu_mat[cbind(site_idx, visit_idx)]
  mu_long <- pmax(mu_long, 1e-10)

  # NB is closed under binomial thinning with the same size, so a negbin fit's
  # per-pass marginal is NB(r, lambda * pi_k) and is scored at its own variance.
  r_long <- .tobs_count_residual(y_long, mu_long, type,
                                 size = object$nmix_dispersion$r %||% Inf,
                                 eps = 1e-10)
  out <- matrix(NA_real_, n_sites, n_pass)
  out[cbind(site_idx, visit_idx)] <- r_long
  list(occ = NULL, det = out)
}


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

#' Simulate removal-sampling abundance data
#'
#' Latent abundance `N_i ~ Poisson(lambda_i)` (or `NegBin(lambda_i, size)`) with
#' `log lambda_i = X_lambda beta_lambda`, observed through `K` ordered removal
#' passes: pass `k` removes `Binomial(N_i - sum_{l<k} y_{il}, p_{ik})` of the
#' individuals still present, with `logit p_{ik} = X_p beta_p` (site-level
#' detection here). The returned `y` is an `N x K` integer matrix of per-pass
#' removals suitable for [tobs()] with [removal()].
#'
#' @param N Number of sites (default 100).
#' @param K Number of removal passes (default 4).
#' @param n.abund.covs Number of abundance covariates (default 2).
#' @param n.det.covs Number of detection covariates (default 1).
#' @param beta.lambda Abundance coefficients on the log scale. Default
#'   `c(log(6), runif(n.abund.covs, -0.5, 0.5))`.
#' @param beta.p Detection coefficients on the logit scale. Default
#'   `c(0.4, runif(n.det.covs, -0.5, 0.5))`.
#' @param mixture `"poisson"` (default) or `"negbin"`.
#' @param size Negative-binomial size `r` (`mixture = "negbin"` only, default 3).
#' @param seed Optional random seed.
#' @return A list with `y` (N x K removal matrix), `data` (covariate data frame),
#'   and `truth` (coefficients, per-site `lambda`, `p`, latent `N`, mixture/size).
#' @examples
#' sim <- simulate_removal(N = 50, K = 3, seed = 1)
#' dim(sim$y)
#' @export
simulate_removal <- function(N = 100, K = 4,
                             n.abund.covs = 2, n.det.covs = 1,
                             beta.lambda = NULL, beta.p = NULL,
                             mixture = c("poisson", "negbin"), size = 3,
                             seed = NULL) {
  mixture <- match.arg(mixture)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(beta.lambda)) beta.lambda <- c(log(6), stats::runif(n.abund.covs, -0.5, 0.5))
  if (is.null(beta.p))      beta.p      <- c(0.4, stats::runif(n.det.covs, -0.5, 0.5))

  abund_covs <- data.frame(matrix(stats::rnorm(N * n.abund.covs), N, n.abund.covs))
  names(abund_covs) <- paste0("abund_cov", seq_len(n.abund.covs))
  det_covs <- data.frame(matrix(stats::rnorm(N * n.det.covs), N, n.det.covs))
  names(det_covs) <- paste0("det_cov", seq_len(n.det.covs))
  data <- cbind(abund_covs, det_covs)

  X_lambda <- stats::model.matrix(~ ., abund_covs)
  X_det    <- stats::model.matrix(~ ., det_covs)
  lambda <- exp(as.vector(X_lambda %*% beta.lambda))
  p      <- plogis(as.vector(X_det %*% beta.p))
  model <- list(X_processes = list(X_lambda, X_det[rep(seq_len(N), each = K), ,
                                                   drop = FALSE]),
                process_info = list(list(p = ncol(X_lambda)), list(p = ncol(X_det))),
                site_idx = rep(seq_len(N), each = K), visit_idx = rep(seq_len(K), N),
                n_sites = N, n_passes = K)
  y <- .tobs_gen_removal(model, list(
    beta = c(beta.lambda, beta.p), off = list(NULL, NULL),
    r = if (identical(mixture, "negbin")) size else NA_real_))
  Nlat <- attr(y, "latent")
  y <- .tobs_gen_strip(y)

  list(
    y = y, data = data,
    truth = list(beta_lambda = beta.lambda, beta_p = beta.p,
                 lambda = lambda, p = p, N = Nlat, mixture = mixture,
                 size = if (identical(mixture, "negbin")) size else NA_real_)
  )
}
