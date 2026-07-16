# occu_ttd.R - Time-to-detection occupancy (Garrard et al. 2008; Bornand et al.
# 2014; unmarked occuTTD). A survey of a site records the TIME to first detection
# rather than a 0/1 detection; the time carries more information than the binary
# outcome. Exponential time-to-detection (constant hazard):
#
#   z_i               ~ Bernoulli(psi_i),          logit psi_i = X_psi_i . beta_psi
#   t_ij | z_i = 1    ~ Exponential(rate lambda_ij), log lambda_i = X_rate_i . beta_rate
#   observe t_ij if t_ij <= Tmax_ij (a detection), else censored at Tmax_ij.
#
# An unoccupied site (z = 0) never detects (all visits censored). The latent z
# integrates out in closed form (two states), so per site with detection rate
# lambda_i and survey lengths Tmax_ij:
#
#   occupied emission = k_i log lambda_i - lambda_i (S_det_i + Tcens_i)
#       k_i     = number of detected visits
#       S_det_i = sum of detection times over detected visits
#       Tcens_i = sum of Tmax over censored (surveyed, undetected) visits
#   L_i = psi_i exp(occ_emission_i)                       if any visit detected
#   L_i = psi_i exp(occ_emission_i) + (1 - psi_i)         if none detected
#
# Fit maximises the exact marginal (optim BFGS), observed-information vcov -- the
# same closed-form-marginal recipe as royle_nichols(). Rate is site-level;
# a Weibull shape and visit-varying rate are documented follow-ups.
#
#   .tobs_build_occu_ttd()   data binder -> model_type = "occu_ttd"
#   .tobs_fit_occu_ttd()     optim over the closed-form two-state marginal
#   .dispatch_occu_ttd()     tobs() entry (bind + fit + assemble)

# ---------------------------------------------------------------------------
# Marginal log-likelihood (per site), the single source of truth reused by the
# fitter, the pointwise log-likelihood, and simulate().
# ---------------------------------------------------------------------------

# Per-site marginal log-likelihood from psi / lambda vectors and the per-site
# sufficient statistics: k (detections), S_det (sum detection times), Tcens (sum
# survey length over censored visits). Returns length-n_sites.
.ttd_site_loglik <- function(psi, lambda, k, S_det, Tcens) {
  psi     <- pmin(pmax(psi, 1e-12), 1 - 1e-12)
  log_psi <- log(psi); log_1mpsi <- log1p(-psi)
  occ_emis <- k * log(lambda) - lambda * (S_det + Tcens)
  any_det  <- k > 0
  # No detection: psi * exp(occ_emis) + (1 - psi), logsumexp for stability.
  a  <- log_psi + occ_emis; b <- log_1mpsi
  mx <- pmax(a, b)
  ll_nodet <- mx + log(exp(a - mx) + exp(b - mx))
  ifelse(any_det, log_psi + occ_emis, ll_nodet)
}

# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# `y` is an n_sites x max_visits time-to-detection matrix: a value in
# (0, Tmax) is a detection time; a value >= Tmax (or exactly Tmax) is a
# non-detection censored at Tmax; NA is a visit not conducted. `surveyLength`
# (Tmax) is a scalar, a per-site vector, or an n_sites x max_visits matrix.
.tobs_build_occu_ttd <- function(state_formula, rate_formula, data, y,
                                 surveyLength = 1) {
  y  <- as.matrix(y)
  n_sites    <- nrow(y)
  max_visits <- ncol(y)
  .tobs_check_site_count(n_sites, nrow(data), "sites")

  Tmax <- surveyLength
  if (length(Tmax) == 1L) Tmax <- matrix(Tmax, n_sites, max_visits)
  else if (is.null(dim(Tmax)) && length(Tmax) == n_sites)
    Tmax <- matrix(Tmax, n_sites, max_visits)
  else Tmax <- as.matrix(Tmax)
  if (!all(dim(Tmax) == dim(y))) {
    stop("surveyLength must be a scalar, a length-n_sites vector, or an ",
         "n_sites x max_visits matrix matching `y`.", call. = FALSE)
  }
  if (any(Tmax[!is.na(y)] <= 0)) {
    stop("surveyLength (Tmax) must be positive at every surveyed visit.",
         call. = FALSE)
  }
  if (any(y[!is.na(y)] < 0)) {
    stop("occu_ttd() y must be non-negative detection times (>= Tmax = ",
         "non-detection, NA = visit not conducted).", call. = FALSE)
  }

  surveyed <- !is.na(y)
  detected <- surveyed & (y < Tmax)                       # 0 < t < Tmax
  k     <- rowSums(detected)                              # detections per site
  Sdet  <- rowSums(ifelse(detected, y, 0))               # sum detection times
  cens  <- surveyed & !detected                          # surveyed, undetected
  Tcens <- rowSums(ifelse(cens, Tmax, 0))                # sum Tmax over censored

  bind    <- .tobs_bind_formulas(list(psi = state_formula, rate = rate_formula),
                                 data)
  X_psi  <- stats::model.matrix(bind$fe$psi,  data)
  X_rate <- stats::model.matrix(bind$fe$rate, data)

  structure(list(
    model_type  = "occu_ttd",
    y           = y,
    Tmax        = Tmax,
    k_site      = as.integer(k),
    Sdet_site   = as.numeric(Sdet),
    Tcens_site  = as.numeric(Tcens),
    X_processes = list(X_psi, X_rate),
    formulas    = list(psi = bind$fe$psi, rate = bind$fe$rate),
    structured_terms = bind$terms,
    data        = data,
    n_sites     = n_sites,
    max_visits  = max_visits,
    process_info = list(
      list(name = "psi",  p = ncol(X_psi),
           coef_names = colnames(X_psi),  link = "logit"),
      list(name = "rate", p = ncol(X_rate),
           coef_names = colnames(X_rate), link = "log")
    )
  ), class = "tobs_model")
}

# ---------------------------------------------------------------------------
# Fitter: maximise the closed-form marginal, observed-information vcov
# ---------------------------------------------------------------------------

.tobs_fit_occu_ttd <- function(model, verbose = TRUE, ...) {
  X_psi  <- model$X_processes[[1L]]
  X_rate <- model$X_processes[[2L]]
  p_psi  <- ncol(X_psi); p_rate <- ncol(X_rate)
  k <- model$k_site; Sdet <- model$Sdet_site; Tcens <- model$Tcens_site

  nll <- function(theta) {
    psi    <- stats::plogis(as.vector(X_psi  %*% theta[seq_len(p_psi)]))
    lambda <- exp(as.vector(X_rate %*% theta[p_psi + seq_len(p_rate)]))
    lambda <- pmin(pmax(lambda, 1e-10), 1e10)
    ll     <- .ttd_site_loglik(psi, lambda, k, Sdet, Tcens)
    val    <- -sum(ll)
    if (is.finite(val)) val else 1e10
  }

  # Moment initialisation. Detected fraction seeds psi; the mean detection time
  # among detections seeds the rate (lambda ~ 1 / mean time).
  det_frac  <- mean(k > 0)
  psi0      <- min(max(det_frac, 0.05), 0.95)
  mean_tdet <- if (sum(k) > 0) sum(Sdet) / sum(k) else 1
  lam0      <- 1 / max(mean_tdet, 1e-3)
  init <- c(stats::qlogis(psi0), rep(0, p_psi - 1L),
            log(max(lam0, 1e-3)), rep(0, p_rate - 1L))

  opt <- stats::optim(init, nll, method = "BFGS", hessian = TRUE,
                      control = list(maxit = 500L))
  converged <- opt$convergence == 0L

  par_names <- c(paste0("psi_",  model$process_info[[1L]]$coef_names),
                 paste0("rate_", model$process_info[[2L]]$coef_names))
  means <- opt$par; names(means) <- par_names
  V <- tryCatch(solve(opt$hessian),
                error = function(e) diag(NA_real_, length(means)))
  V <- (V + t(V)) / 2
  dimnames(V) <- list(par_names, par_names)
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- par_names

  n_draws <- 1000L
  draws <- .occu_cover_rmvn(n_draws, means, V)
  colnames(draws) <- par_names

  intercepts <- list(
    psi  = stats::setNames(means[1L], par_names[1L]),
    rate = stats::setNames(means[p_psi + 1L], par_names[p_psi + 1L]))

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = length(means),
    log_prob     = rep(-opt$value, n_draws),
    log_lik      = -opt$value,
    N            = model$n_sites),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names    = par_names,
    param_names  = par_names,
    n_fixed      = length(means),
    fixed_names  = par_names,
    intercepts   = intercepts,
    process_info = model$process_info,
    model        = model,
    spatial      = NULL,
    method       = "laplace",
    convergence  = list(converged = converged, n_iter = opt$counts[[1L]])
  )), class = c("tobs_fit", "tulpa_fit"))
}

# ---------------------------------------------------------------------------
# tobs() dispatcher
# ---------------------------------------------------------------------------

.dispatch_occu_ttd <- function(formula, data, family, detection, y, visits,
                               engine, priors, control,
                               approx = "gaussian_laplace",
                               correction = "none", ...) {
  if (is.null(detection)) {
    stop("occu_ttd() requires a `detection` formula (the site-level ",
         "log detection-rate model).", call. = FALSE)
  }
  if (is.null(y)) {
    stop("occu_ttd() requires `y` (an N x J time-to-detection matrix; a value ",
         ">= surveyLength is a non-detection, NA a visit not conducted).",
         call. = FALSE)
  }
  if (!is.null(visits)) {
    stop("occu_ttd() detection is site-level; visit-level rate covariates ",
         "(`visits`) are not yet supported.", call. = FALSE)
  }
  if (!identical(.map_engine(engine, family = "occu_ttd"), "laplace")) {
    stop("occu_ttd() supports method = \"laplace\" only.", call. = FALSE)
  }
  model <- .tobs_build_occu_ttd(
    state_formula = formula, rate_formula = detection, data = data, y = y,
    surveyLength = family$params$surveyLength %||% 1)
  .tobs_fit_occu_ttd(model, verbose = isTRUE(control$verbose))
}

# ---------------------------------------------------------------------------
# fitted / predict / residuals / pointwise log-likelihood / simulate
# ---------------------------------------------------------------------------

# fitted(): per-site psi (occupancy), lambda (detection rate), and the marginal
# per-visit detection probability within a survey of length mean(Tmax),
# p = 1 - exp(-lambda * Tmax).
.tobs_fitted_occu_ttd <- function(object) {
  model  <- object$model
  p_psi  <- model$process_info[[1L]]$p
  beta   <- object$means
  psi    <- stats::plogis(as.vector(model$X_processes[[1L]] %*% beta[seq_len(p_psi)]))
  lambda <- exp(as.vector(model$X_processes[[2L]] %*%
                            beta[p_psi + seq_len(model$process_info[[2L]]$p)]))
  Tbar   <- rowMeans(model$Tmax)
  list(psi = psi, rate = lambda, p = 1 - exp(-lambda * Tbar))
}

# predict(): psi (occupancy, default) or rate (detection rate) at the fitted or
# new design, on the response scale.
.tobs_predict_occu_ttd <- function(object, newdata = NULL,
                                   type = c("state", "detection")) {
  type  <- match.arg(type)
  model <- object$model
  k     <- if (identical(type, "detection")) 2L else 1L
  X     <- if (is.null(newdata)) model$X_processes[[k]]
           else stats::model.matrix(model$formulas[[if (k == 1L) "psi" else "rate"]],
                                    newdata)
  off   <- if (k == 1L) 0L else model$process_info[[1L]]$p
  eta   <- as.vector(X %*% object$means[off + seq_len(model$process_info[[k]]$p)])
  if (identical(type, "detection")) exp(eta) else stats::plogis(eta)
}

# residuals(): per-site response / pearson / deviance on the observed
# any-detection indicator against the marginal detection probability
# psi * (1 - exp(-lambda * sum Tmax)).
.tobs_residuals_occu_ttd <- function(object, type) {
  model <- object$model
  p_psi <- model$process_info[[1L]]$p
  beta  <- object$means
  psi   <- stats::plogis(as.vector(model$X_processes[[1L]] %*% beta[seq_len(p_psi)]))
  lam   <- exp(as.vector(model$X_processes[[2L]] %*%
                           beta[p_psi + seq_len(model$process_info[[2L]]$p)]))
  Ttot  <- rowSums(model$Tmax * !is.na(model$y))
  pdet  <- psi * (1 - exp(-lam * Ttot))                  # P(>=1 detection)
  obs   <- as.numeric(model$k_site > 0)
  eps   <- 1e-10; p <- pmin(pmax(pdet, eps), 1 - eps)
  res <- switch(type,
    response = obs - p,
    pearson  = (obs - p) / sqrt(p * (1 - p) + eps),
    deviance = sign(obs - p) * sqrt(2 * abs(
      ifelse(obs > 0, obs * log(obs / p), 0) +
      ifelse(obs < 1, (1 - obs) * log((1 - obs) / (1 - p)), 0))))
  list(occ = res, det = NULL)
}

# Posterior replicate TTD histories: draw a coefficient vector, z ~ Bern(psi),
# and per surveyed visit an exponential time (censored at Tmax) when occupied.
.tobs_simulate_occu_ttd <- function(object, nsim = 1) {
  model <- object$model
  p_psi <- model$process_info[[1L]]$p; p_rate <- model$process_info[[2L]]$p
  Xp <- model$X_processes[[1L]]; Xr <- model$X_processes[[2L]]
  surveyed <- !is.na(model$y); Tmax <- model$Tmax
  draw_one <- function() {
    idx    <- sample.int(nrow(object$draws), 1L)
    beta   <- object$draws[idx, ]
    psi    <- stats::plogis(as.vector(Xp %*% beta[seq_len(p_psi)]))
    lambda <- exp(as.vector(Xr %*% beta[p_psi + seq_len(p_rate)]))
    z      <- stats::rbinom(model$n_sites, 1L, psi)
    yy     <- matrix(NA_real_, model$n_sites, model$max_visits)
    for (i in seq_len(model$n_sites)) {
      for (j in seq_len(model$max_visits)) {
        if (!surveyed[i, j]) next
        if (z[i] == 0L) { yy[i, j] <- Tmax[i, j]; next }   # unoccupied: censored
        t <- stats::rexp(1L, lambda[i])
        yy[i, j] <- if (t < Tmax[i, j]) t else Tmax[i, j]
      }
    }
    yy
  }
  if (nsim == 1L) return(draw_one())
  lapply(seq_len(nsim), function(s) draw_one())
}

# Pointwise log-likelihood [n_draws x n_sites] over the posterior draws.
.tobs_ploglik_occu_ttd <- function(object, n.draws = 1000L, n.threads = 1L) {
  model <- object$model
  draws <- object$draws
  if (!is.null(n.draws) && as.integer(n.draws) < nrow(draws)) {
    draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
  }
  p_psi <- model$process_info[[1L]]$p; p_rate <- model$process_info[[2L]]$p
  Xp <- model$X_processes[[1L]]; Xr <- model$X_processes[[2L]]
  k <- model$k_site; Sdet <- model$Sdet_site; Tcens <- model$Tcens_site
  t(vapply(seq_len(nrow(draws)), function(d) {
    psi    <- stats::plogis(as.vector(Xp %*% draws[d, seq_len(p_psi)]))
    lambda <- exp(as.vector(Xr %*% draws[d, p_psi + seq_len(p_rate)]))
    lambda <- pmin(pmax(lambda, 1e-10), 1e10)
    .ttd_site_loglik(psi, lambda, k, Sdet, Tcens)
  }, numeric(model$n_sites)))
}

# ---------------------------------------------------------------------------
# Simulator for recovery tests
# ---------------------------------------------------------------------------

#' Simulate a time-to-detection occupancy data set
#'
#' Draws from the [occu_ttd()] model: site occupancy `z ~ Bernoulli(psi)` and,
#' at occupied sites, exponential time-to-detection at rate `lambda` over a
#' survey of length `Tmax` (a non-detection is censored at `Tmax`).
#'
#' @param N Number of sites (default 200).
#' @param J Number of replicate surveys per site (default 4).
#' @param n_psi_covs,n_rate_covs Number of occupancy / rate covariates.
#' @param beta_psi Occupancy coefficients `c(intercept, slopes...)` (logit).
#'   Default `c(qlogis(0.6), runif(n_psi_covs, -0.5, 0.5))`.
#' @param beta_rate Log detection-rate coefficients `c(intercept, slopes...)`.
#'   Default `c(log(0.6), runif(n_rate_covs, -0.4, 0.4))`.
#' @param Tmax Survey length (scalar). Default 3.
#' @param seed Optional random seed.
#' @return A list with `y` (N x J TTD matrix), `data`, `Tmax`, and `truth`.
#' @export
simulate_occu_ttd <- function(N = 200, J = 4, n_psi_covs = 1, n_rate_covs = 1,
                              beta_psi = NULL, beta_rate = NULL, Tmax = 3,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(beta_psi))
    beta_psi  <- c(stats::qlogis(0.6), stats::runif(n_psi_covs, -0.5, 0.5))
  if (is.null(beta_rate))
    beta_rate <- c(log(0.6), stats::runif(n_rate_covs, -0.4, 0.4))

  psi_covs  <- data.frame(matrix(stats::rnorm(N * n_psi_covs), N, n_psi_covs))
  names(psi_covs) <- paste0("psi_cov", seq_len(n_psi_covs))
  rate_covs <- data.frame(matrix(stats::rnorm(N * n_rate_covs), N, n_rate_covs))
  names(rate_covs) <- paste0("rate_cov", seq_len(n_rate_covs))
  data <- cbind(psi_covs, rate_covs)

  X_psi  <- stats::model.matrix(~ ., psi_covs)
  X_rate <- stats::model.matrix(~ ., rate_covs)
  psi    <- plogis(as.vector(X_psi  %*% beta_psi))
  lambda <- exp(as.vector(X_rate %*% beta_rate))
  z      <- stats::rbinom(N, 1L, psi)

  y <- matrix(NA_real_, N, J)
  for (i in seq_len(N)) {
    for (j in seq_len(J)) {
      if (z[i] == 0L) { y[i, j] <- Tmax; next }
      t <- stats::rexp(1L, lambda[i])
      y[i, j] <- if (t < Tmax) t else Tmax
    }
  }

  list(y = y, data = data, Tmax = Tmax,
       truth = list(beta_psi = beta_psi, beta_rate = beta_rate,
                    psi = psi, rate = lambda, z = z, Tmax = Tmax))
}
