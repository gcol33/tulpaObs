# gdistremoval.R - joint distance + removal sampling (unmarked gdistremoval,
# Amundson et al. 2014). SINGLE-SEASON (not the open-population distsampOpen HMM).
#
# One site i, latent abundance N_i ~ Poisson(lambda_i). The detected birds are
# cross-classified two ways: by a distance band (distance sampling) AND by the
# removal period of first detection (removal sampling). With overall distance
# detection pdist_i, overall removal detection prem_i, and availability phi
# (fixed 1 in this single-primary-period model), the total detected count is a
# binomial thinning of N_i, and Poisson is closed under binomial thinning, so
#
#   ysum_i ~ Poisson(lambda_i * pdist_i * prem_i)                 (total detected)
#   yDist_i | ysum_i ~ Multinomial(ysum_i; cpd_i / pdist_i)       (band allocation)
#   yRem_i  | ysum_i ~ Multinomial(ysum_i; cpr_i / prem_i)        (period allocation)
#
# where the three pieces are conditionally independent, so the marginal is a
# closed-form product per site -- the double_observer() Poisson-multinomial
# pattern, with a distance multinomial + a depleting-removal multinomial. The
# distance band probs cpd_b = integral_band g(x; sigma) f(x) dx have a closed
# form under the half-normal key (line: f = 1/W; point: f = 2x/W^2); the removal
# probs are the depleting-catch multinomial cpr_k = r (1 - r)^(k-1).
#
#   log lambda_i = X_lambda_i . beta_lambda         (abundance)
#   log sigma_i  = X_sigma_i  . beta_sigma          (distance detection scale)
#   logit r_i    = X_r_i      . beta_r              (per-period removal capture)
#
#   .tobs_build_gdistremoval()   data binder -> model_type = "gdistremoval"
#   .tobs_fit_gdistremoval()     optim over the closed-form marginal
#   .dispatch_gdistremoval()     tobs() entry (bind + fit + assemble)
#
# Scope (v1): half-normal key, line / point transect, Poisson abundance, constant
# per-period removal capture, availability fixed at 1 (a single primary period
# does not identify phi separately from the detection probabilities). Hazard-rate
# key, NB / ZIP abundance, and a phi arm over multiple primary periods are
# follow-ups (each closed-form under the same thinning).

# ---------------------------------------------------------------------------
# Cell-probability helpers (vectorised over sites)
# ---------------------------------------------------------------------------

# Half-normal per-band detection-cell probabilities cpd[i, b] for per-site sigma,
# by the closed-form bin integral of g(x) f(x). Returns an [n_sites x n_bins]
# matrix; row sums are the overall distance detection probs pdist_i.
.gdr_dist_cp <- function(sigma, cutpoints, transect) {
  W <- cutpoints[length(cutpoints)]
  if (identical(transect, "point")) {
    # f(x) = 2x / W^2; integral exp(-x^2/2s^2) 2x/W^2 dx = 2 s^2 / W^2 [e(a) - e(b)]
    E   <- exp(-outer(1 / (2 * sigma^2), cutpoints^2))         # [n x (B+1)]
    cp  <- (2 * sigma^2 / W^2) * (E[, -ncol(E), drop = FALSE] -
                                  E[, -1, drop = FALSE])
  } else {
    # line: f(x) = 1/W; integral exp(-x^2/2s^2) dx = s sqrt(2pi) [Phi(b/s) - Phi(a/s)]
    Phi <- stats::pnorm(outer(1 / sigma, cutpoints))           # [n x (B+1)]
    cp  <- (sigma * sqrt(2 * pi) / W) * (Phi[, -1, drop = FALSE] -
                                         Phi[, -ncol(Phi), drop = FALSE])
  }
  pmax(cp, 0)
}

# Depleting-removal per-period cell probabilities cpr[i, k] = r_i (1 - r_i)^(k-1)
# for per-site capture r over Jrem periods. Row sums are prem_i = 1 - (1-r)^Jrem.
.gdr_rem_cp <- function(r, Jrem) {
  r <- pmin(pmax(r, 1e-10), 1 - 1e-10)
  k <- matrix(seq_len(Jrem) - 1L, nrow = length(r), ncol = Jrem, byrow = TRUE)
  r * (1 - r)^k
}

# Per-site marginal log-likelihood (vectorised). lambda / sigma / r length
# n_sites; yDist [n_sites x Jdist]; yRem [n_sites x Jrem].
.gdr_site_loglik <- function(lambda, sigma, r, yDist, yRem, cutpoints, transect) {
  cpd   <- .gdr_dist_cp(sigma, cutpoints, transect)
  cpr   <- .gdr_rem_cp(r, ncol(yRem))
  pdist <- rowSums(cpd); prem <- rowSums(cpr)
  ok    <- is.finite(pdist) & pdist > 1e-12 & prem > 1e-12
  ysum  <- rowSums(yDist)

  ll <- rep(-1e10, length(lambda))
  if (!any(ok)) return(ll)
  mu <- lambda * pdist * prem                    # E[ysum] with phi = 1
  llc <- stats::dpois(ysum, mu, log = TRUE)      # total-detected Poisson
  # conditional band / period multinomials (only where individuals were detected)
  pid <- cpd / pdist; pir <- cpr / prem
  mdist <- lgamma(ysum + 1) - rowSums(lgamma(yDist + 1)) +
           rowSums(yDist * log(pmax(pid, 1e-300)))
  mrem  <- lgamma(ysum + 1) - rowSums(lgamma(yRem + 1)) +
           rowSums(yRem * log(pmax(pir, 1e-300)))
  has   <- ysum > 0
  site  <- llc + ifelse(has, mdist + mrem, 0)
  ll[ok] <- site[ok]
  ll
}

# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# `y` = yDist (n_sites x Jdist distance-band counts); `y_rem` = yRem
# (n_sites x Jrem removal-period counts). Per-site row totals must match (the
# same detected birds cross-classified). `abund_formula` models log lambda,
# `det_formula` log sigma, `rem_formula` logit r.
.tobs_build_gdistremoval <- function(abund_formula, det_formula, rem_formula,
                                     data, y, y_rem, cutpoints, transect) {
  yD <- as.matrix(y); yR <- as.matrix(y_rem)
  if (any(yD < 0 | yD != round(yD), na.rm = TRUE) ||
      any(yR < 0 | yR != round(yR), na.rm = TRUE)) {
    stop("gdistremoval() y / y_rem must be non-negative integer counts.",
         call. = FALSE)
  }
  if (nrow(yD) != nrow(yR)) {
    stop("gdistremoval() y and y_rem must have the same number of sites (rows).",
         call. = FALSE)
  }
  # The same detected birds are cross-classified, so the per-site totals MUST
  # agree; a mismatch is a data error, surfaced up front (input-totals rule).
  if (!isTRUE(all.equal(rowSums(yD), rowSums(yR)))) {
    stop("gdistremoval() per-site totals of y (distance) and y_rem (removal) ",
         "must match: rowSums(y) == rowSums(y_rem) (the same detected birds ",
         "cross-classified by distance band and removal period).", call. = FALSE)
  }
  if (length(cutpoints) != ncol(yD) + 1L) {
    stop(sprintf(paste0("gdistremoval() cutpoints must have length ncol(y) + 1 ",
         "= %d (the distance-bin edges 0 = c_0 < ... < c_B)."), ncol(yD) + 1L),
         call. = FALSE)
  }
  cutpoints <- as.numeric(cutpoints)
  if (any(diff(cutpoints) <= 0) || cutpoints[1] < 0) {
    stop("gdistremoval() cutpoints must be strictly increasing and start >= 0.",
         call. = FALSE)
  }
  n_sites <- nrow(yD)
  .tobs_check_site_count(n_sites, nrow(data), "sites")
  storage.mode(yD) <- "integer"; storage.mode(yR) <- "integer"

  bind <- .tobs_bind_formulas(
    list(lambda = abund_formula, sigma = det_formula, r = rem_formula), data)
  X_lambda <- stats::model.matrix(bind$fe$lambda, data)
  X_sigma  <- stats::model.matrix(bind$fe$sigma, data)
  X_r      <- stats::model.matrix(bind$fe$r, data)

  structure(list(
    model_type  = "gdistremoval",
    y           = yD,
    y_rem       = yR,
    cutpoints   = cutpoints,
    transect    = transect,
    n_bins      = ncol(yD),
    n_periods   = ncol(yR),
    X_processes = list(X_lambda, X_sigma, X_r),
    formulas    = list(lambda = bind$fe$lambda, sigma = bind$fe$sigma,
                       r = bind$fe$r),
    structured_terms = bind$terms,
    data        = data,
    n_sites     = n_sites,
    process_info = list(
      list(name = "lambda", p = ncol(X_lambda),
           coef_names = colnames(X_lambda), link = "log"),
      list(name = "sigma", p = ncol(X_sigma),
           coef_names = colnames(X_sigma), link = "log"),
      list(name = "r", p = ncol(X_r), coef_names = colnames(X_r),
           link = "logit")
    )
  ), class = "tobs_model")
}

.gdr_unpack <- function(theta, model) {
  Xl <- model$X_processes[[1L]]; Xs <- model$X_processes[[2L]]
  Xr <- model$X_processes[[3L]]
  pl <- ncol(Xl); ps <- ncol(Xs); pr <- ncol(Xr)
  list(
    lambda = exp(as.vector(Xl %*% theta[seq_len(pl)])),
    sigma  = exp(as.vector(Xs %*% theta[pl + seq_len(ps)])),
    r      = stats::plogis(as.vector(Xr %*% theta[pl + ps + seq_len(pr)]))
  )
}

# ---------------------------------------------------------------------------
# Fitter (called from .dispatch_gdistremoval)
# ---------------------------------------------------------------------------

.tobs_fit_gdistremoval <- function(model, verbose = TRUE, ...) {
  Xl <- model$X_processes[[1L]]; Xs <- model$X_processes[[2L]]
  Xr <- model$X_processes[[3L]]
  pl <- ncol(Xl); ps <- ncol(Xs); pr <- ncol(Xr)
  yD <- model$y; yR <- model$y_rem
  cutpoints <- model$cutpoints; transect <- model$transect

  nll <- function(theta) {
    up  <- .gdr_unpack(theta, model)
    val <- -sum(.gdr_site_loglik(up$lambda, up$sigma, up$r, yD, yR,
                                 cutpoints, transect))
    if (is.finite(val)) val else 1e10
  }

  # Moment init: total detected seeds lambda given a rough overall detection; the
  # distance scale seeds at the median cutpoint; the depletion fraction (first
  # period share of the removal counts) seeds r.
  ntot   <- rowSums(yD)
  r_hat  <- {
    tot <- sum(yR)
    if (tot > 0) min(max(sum(yR[, 1L]) / tot, 0.05), 0.9) else 0.4
  }
  sig0   <- stats::median(cutpoints[-1])
  prem0  <- 1 - (1 - r_hat)^model$n_periods
  pdist0 <- sum(.gdr_dist_cp(sig0, cutpoints, transect))
  lam0   <- mean(ntot) / max(pdist0 * prem0, 0.05)
  init   <- c(log(max(lam0, 1e-2)), rep(0, pl - 1L),
              log(max(sig0, 1e-2)),  rep(0, ps - 1L),
              stats::qlogis(r_hat),  rep(0, pr - 1L))

  par_names <- unlist(lapply(model$process_info, function(pp)
    paste0(pp$name, "_", pp$coef_names)))

  .tobs_bfgs_marginal_fit(nll, init, par_names, model, N = model$n_sites)
}

# ---------------------------------------------------------------------------
# tobs() dispatcher
# ---------------------------------------------------------------------------

.dispatch_gdistremoval <- function(formula, data, family, detection, y, visits,
                                   engine, priors, control,
                                   approx = "gaussian_laplace",
                                   correction = "none", ...) {
  dots <- list(...)
  if (is.null(detection))
    stop("gdistremoval() requires a `detection` formula (the site-level ",
         "log-sigma distance-scale model).", call. = FALSE)
  if (is.null(y))
    stop("gdistremoval() requires `y` (an n_sites x n_bins integer matrix of ",
         "per-distance-band detected counts).", call. = FALSE)
  if (is.null(dots$y_rem))
    stop("gdistremoval() requires `y_rem` (an n_sites x n_periods integer ",
         "matrix of per-removal-period detected counts).", call. = FALSE)
  if (!is.null(visits))
    stop("gdistremoval() detection is site-level; visit-level covariates ",
         "(`visits`) are not yet supported.", call. = FALSE)
  cutpoints <- family$params$cutpoints
  if (is.null(cutpoints))
    stop("gdistremoval() requires `cutpoints` (the distance-bin edges); pass ",
         "gdistremoval(cutpoints = ...).", call. = FALSE)
  model <- .tobs_build_gdistremoval(
    abund_formula = formula, det_formula = detection,
    rem_formula = dots$removal %||% ~1, data = data, y = y, y_rem = dots$y_rem,
    cutpoints = cutpoints, transect = family$params$transect)
  .tobs_reject_unwired_structs(
    model, "gdistremoval()",
    hint = paste0("the joint distance-removal marginal is fitted on fixed ",
                  "effects only, so drop the term"))
  .tobs_fit_gdistremoval(model, verbose = isTRUE(control$verbose))
}

# ---------------------------------------------------------------------------
# fitted / predict / residuals / pointwise log-likelihood / simulate
# ---------------------------------------------------------------------------

# fitted(): per-site abundance lambda, distance scale sigma, removal capture r,
# and the overall detection probs pdist / prem and expected total count.
.tobs_fitted_gdistremoval <- function(object) {
  model <- object$model
  up <- .gdr_unpack(object$means, model)
  pdist <- rowSums(.gdr_dist_cp(up$sigma, model$cutpoints, model$transect))
  prem  <- rowSums(.gdr_rem_cp(up$r, model$n_periods))
  list(lambda = up$lambda, sigma = up$sigma, r = up$r,
       pdist = pdist, prem = prem, count = up$lambda * pdist * prem)
}

# predict(): abundance (default), distance scale, or removal capture.
.tobs_predict_gdistremoval <- function(object, newdata = NULL,
                                       type = c("abundance", "distance",
                                                "removal")) {
  type  <- match.arg(type)
  model <- object$model
  b <- object$means
  pl <- ncol(model$X_processes[[1L]]); ps <- ncol(model$X_processes[[2L]])
  pr <- ncol(model$X_processes[[3L]])
  Xnew <- function(f, Xdefault)
    if (is.null(newdata)) Xdefault else stats::model.matrix(f, newdata)
  switch(type,
    abundance = exp(as.vector(Xnew(model$formulas$lambda,
                 model$X_processes[[1L]]) %*% b[seq_len(pl)])),
    distance  = exp(as.vector(Xnew(model$formulas$sigma,
                 model$X_processes[[2L]]) %*% b[pl + seq_len(ps)])),
    removal   = stats::plogis(as.vector(Xnew(model$formulas$r,
                 model$X_processes[[3L]]) %*% b[pl + ps + seq_len(pr)])))
}

# residuals(): per-site Pearson / deviance on the total detected count against
# its expected value lambda * pdist * prem.
.tobs_residuals_gdistremoval <- function(object, type) {
  fv   <- .tobs_fitted_gdistremoval(object)
  mtot <- fv$count
  obs  <- rowSums(object$model$y)
  eps  <- 1e-10
  res <- switch(type,
    response = obs - mtot,
    pearson  = (obs - mtot) / sqrt(mtot + eps),
    deviance = sign(obs - mtot) * sqrt(2 * abs(
      ifelse(obs > 0, obs * log(obs / pmax(mtot, eps)), 0) - (obs - mtot))))
  list(occ = res, det = NULL)
}

# Pointwise log-likelihood [n_draws x n_sites] over the posterior draws.
.tobs_ploglik_gdistremoval <- function(object, n.draws = 1000L, n.threads = 1L) {
  model <- object$model
  draws <- object$draws
  if (!is.null(n.draws) && as.integer(n.draws) < nrow(draws)) {
    draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
  }
  yD <- model$y; yR <- model$y_rem
  t(vapply(seq_len(nrow(draws)), function(d) {
    up <- .gdr_unpack(draws[d, ], model)
    .gdr_site_loglik(up$lambda, up$sigma, up$r, yD, yR,
                     model$cutpoints, model$transect)
  }, numeric(model$n_sites)))
}

# Posterior replicate counts: draw a coefficient vector, then per site draw N,
# the detected total, and the distance-band + removal-period allocations.
.tobs_simulate_gdistremoval <- function(object, nsim = 1) {
  model <- object$model
  draw_one <- function() {
    idx <- sample.int(nrow(object$draws), 1L)
    up  <- .gdr_unpack(object$draws[idx, ], model)
    .gdr_draw(up$lambda, up$sigma, up$r, model$cutpoints, model$transect,
              model$n_periods)
  }
  if (nsim == 1L) return(draw_one())
  lapply(seq_len(nsim), function(s) draw_one())
}

# Core draw: per site N ~ Pois(lambda), detected ~ Binom(N, pdist*prem), then the
# detected allocated to distance bands ~ Mult(cpd/pdist) and removal periods ~
# Mult(cpr/prem) independently (the unmarked model's cross-classification).
.gdr_draw <- function(lambda, sigma, r, cutpoints, transect, n_periods) {
  n <- length(lambda)
  cpd <- .gdr_dist_cp(sigma, cutpoints, transect); pdist <- rowSums(cpd)
  cpr <- .gdr_rem_cp(r, n_periods);                prem  <- rowSums(cpr)
  N   <- stats::rpois(n, lambda)
  yD  <- matrix(0L, n, ncol(cpd)); yR <- matrix(0L, n, n_periods)
  for (i in seq_len(n)) {
    if (N[i] == 0L) next
    det <- stats::rbinom(1L, N[i], pdist[i] * prem[i])
    if (det > 0L) {
      yD[i, ] <- as.integer(stats::rmultinom(1L, det, cpd[i, ] / pdist[i]))
      yR[i, ] <- as.integer(stats::rmultinom(1L, det, cpr[i, ] / prem[i]))
    }
  }
  list(yDist = yD, yRem = yR)
}

# ---------------------------------------------------------------------------
# Simulator for recovery tests
# ---------------------------------------------------------------------------

#' Simulate a joint distance + removal data set
#'
#' Draws from the [gdistremoval()] model: site abundance `N ~ Poisson(lambda)`,
#' the detected birds cross-classified by a distance band (half-normal key, scale
#' `sigma`) and by a removal period (per-period capture `r`).
#'
#' @param N Number of sites (default 200).
#' @param cutpoints Distance-bin edges `0 = c_0 < ... < c_B`.
#' @param n_periods Number of removal periods.
#' @param transect `"line"` (default) or `"point"`.
#' @param n_abund_covs,n_det_covs,n_rem_covs Number of abundance / distance /
#'   removal covariates.
#' @param beta_lambda,beta_sigma,beta_r Coefficients on the log-abundance,
#'   log-scale, and logit-removal arms (`c(intercept, slopes...)`). Defaults give
#'   moderate abundance / detection.
#' @param seed Optional random seed.
#' @return A list with `y` (distance-band counts), `y_rem` (removal-period
#'   counts), `data`, and `truth`.
#' @export
simulate_gdistremoval <- function(N = 200, cutpoints = c(0, 10, 20, 30, 40),
                                  n_periods = 4L, transect = "line",
                                  n_abund_covs = 1, n_det_covs = 1, n_rem_covs = 1,
                                  beta_lambda = NULL, beta_sigma = NULL,
                                  beta_r = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(beta_lambda))
    beta_lambda <- c(log(25), stats::runif(n_abund_covs, -0.4, 0.4))
  if (is.null(beta_sigma))
    beta_sigma <- c(log(stats::median(cutpoints[-1])),
                    stats::runif(n_det_covs, -0.2, 0.2))
  if (is.null(beta_r))
    beta_r <- c(stats::qlogis(0.4), stats::runif(n_rem_covs, -0.3, 0.3))

  mk <- function(k, tag) {
    d <- data.frame(matrix(stats::rnorm(N * k), N, k))
    names(d) <- paste0(tag, seq_len(k)); d
  }
  ac <- mk(n_abund_covs, "abund_cov"); dc <- mk(n_det_covs, "det_cov")
  rc <- mk(n_rem_covs, "rem_cov")
  data <- cbind(ac, dc, rc)

  lambda <- exp(as.vector(stats::model.matrix(~ ., ac) %*% beta_lambda))
  sigma  <- exp(as.vector(stats::model.matrix(~ ., dc) %*% beta_sigma))
  r      <- plogis(as.vector(stats::model.matrix(~ ., rc) %*% beta_r))

  dr <- .gdr_draw(lambda, sigma, r, as.numeric(cutpoints), transect,
                  as.integer(n_periods))
  colnames(dr$yDist) <- paste0("band", seq_len(ncol(dr$yDist)))
  colnames(dr$yRem)  <- paste0("period", seq_len(n_periods))
  list(y = dr$yDist, y_rem = dr$yRem, data = data,
       truth = list(beta_lambda = beta_lambda, beta_sigma = beta_sigma,
                    beta_r = beta_r, lambda = lambda, sigma = sigma, r = r))
}
