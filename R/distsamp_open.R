# distsamp_open.R - open-population distance sampling (unmarked distsampOpen).
# A Dail-Madsen open N-mixture (as dyn_abun(), #37) with a distance-bin
# MULTINOMIAL emission at each primary period instead of the binomial. This is the
# open-population counterpart of the single-season gdistremoval().
#
#   N_{i,1}   ~ Poisson(lambda_i)                                   (initial abundance)
#   N_{i,t}   = Binomial(N_{i,t-1}, omega_i) + Poisson(gamma_i)     (Dail-Madsen "constant")
#   n_{i,t}   ~ Binomial(N_{i,t}, pdist_i)                          (detected total, distance)
#   yDist_{i,t,.} | n_{i,t} ~ Multinomial(cpd_i / pdist_i)         (distance-band allocation)
#
# The band allocation is CONDITIONAL on the period total, independent of the latent
# N, so it factors out of the HMM sum (the gdistremoval trick):
#
#   site_ll = sum_t dmultinom(yDist[i,t,.] ; cpd_i/pdist_i)                 (band allocations)
#           + HMM_forward( y = n_{i,.}, p_t = pdist_i, lambda, omega, gamma )
#
# and the HMM-forward over the latent N sequence with a per-period binomial
# emission dbinom(n_{i,t}; N_{i,t}, pdist_i) is EXACTLY the dyn_abun() marginal, so
# the existing (validated) cpp_dyn_abun_total_log_lik kernel is reused: feed it
# eta_p = logit(pdist_i) (the overall distance detection enters as the detection
# logit) and y = the per-period detected totals. No new HMM kernel; the only new
# pieces are pdist(sigma), the per-period band multinomials, and their assembly.
#
#   log lambda_i = X_lambda_i . beta_lambda    (initial abundance)
#   log sigma_i  = X_sigma_i  . beta_sigma     (distance detection scale)
#   logit omega_i= X_omega_i  . beta_omega     (apparent survival)
#   log gamma_i  = X_gamma_i  . beta_gamma     (recruitment)
#
#   .tobs_build_distsamp_open()   data binder -> model_type = "distsamp_open"
#   .tobs_fit_distsamp_open()     optim over the composed marginal
#   .dispatch_distsamp_open()     tobs() entry
#
# Scope (v1): half-normal key, line / point transect, Poisson initial abundance,
# constant Dail-Madsen dynamics, site-level arms. NB / ZIP initial abundance, other
# dynamics (autoreg / ricker / gompertz), and season-varying sigma are follow-ups.

# ---------------------------------------------------------------------------
# Composed marginal
# ---------------------------------------------------------------------------

# Per-site log-likelihood of the open-population distance model. `theta` packs
# c(beta_lambda, beta_sigma, beta_omega, beta_gamma); returns the summed
# log-likelihood (or a large penalty on an infeasible detection).
.dso_negll <- function(theta, model) {
  Xl <- model$X_processes[[1L]]; Xs <- model$X_processes[[2L]]
  Xo <- model$X_processes[[3L]]; Xg <- model$X_processes[[4L]]
  pl <- ncol(Xl); ps <- ncol(Xs); po <- ncol(Xo); pg <- ncol(Xg)
  bl <- theta[seq_len(pl)]
  bs <- theta[pl + seq_len(ps)]
  bo <- theta[pl + ps + seq_len(po)]
  bg <- theta[pl + ps + po + seq_len(pg)]

  sigma <- exp(as.vector(Xs %*% bs))
  cpd   <- .gdr_dist_cp(sigma, model$cutpoints, model$transect)   # [n x Jbin]
  pdist <- rowSums(cpd)
  if (any(!is.finite(pdist)) || any(pdist <= 1e-8) || any(pdist >= 1 - 1e-10))
    return(1e10)

  ev <- cpp_dyn_abun_total_log_lik(
    model$y_flat, model$n_sites, model$n_seasons, 1L, model$K_max,
    as.vector(Xl %*% bl),                 # eta_lambda (log link)
    stats::qlogis(pdist),                 # eta_p = logit(pdist)  (the reuse)
    as.vector(Xo %*% bo),                 # eta_omega (logit link)
    as.vector(Xg %*% bg),                 # eta_gamma (log link)
    use_nb = FALSE, eta_logr = 0)
  hmm <- sum(ev$log_lik)

  # Per-period distance-band multinomials (a function of sigma only).
  pid  <- cpd / pdist                                             # [n x Jbin]
  band <- 0
  for (t in seq_len(model$n_seasons)) {
    yb   <- model$y[, , t]                                        # [n x Jbin]
    band <- band + sum(lgamma(rowSums(yb) + 1) - rowSums(lgamma(yb + 1)) +
                       rowSums(yb * log(pmax(pid, 1e-300))))
  }
  val <- -(hmm + band)
  if (is.finite(val)) val else 1e10
}

# Analytic gradient of the LOG-likelihood wrt theta = c(beta_lambda, beta_sigma,
# beta_omega, beta_gamma). The dyn_abun kernel returns grad_eta_{lambda,p,omega,
# gamma} in the same call as the value, so lambda / omega / gamma chain straight
# through their designs. sigma enters only via pdist (the detection logit) and the
# band multinomials; the distance-integral derivative dcpd/dsigma is a cheap
# central finite difference on .gdr_dist_cp (NO extra HMM evaluation), so the whole
# gradient costs one kernel call. Returns a length-n_par numeric vector.
.dso_grad <- function(theta, model) {
  Xl <- model$X_processes[[1L]]; Xs <- model$X_processes[[2L]]
  Xo <- model$X_processes[[3L]]; Xg <- model$X_processes[[4L]]
  pl <- ncol(Xl); ps <- ncol(Xs); po <- ncol(Xo); pg <- ncol(Xg)
  bl <- theta[seq_len(pl)]
  bs <- theta[pl + seq_len(ps)]
  bo <- theta[pl + ps + seq_len(po)]
  bg <- theta[pl + ps + po + seq_len(pg)]

  sigma <- exp(as.vector(Xs %*% bs))
  cpd   <- .gdr_dist_cp(sigma, model$cutpoints, model$transect)
  pdist <- rowSums(cpd)
  pd    <- pmin(pmax(pdist, 1e-10), 1 - 1e-10)

  ev <- cpp_dyn_abun_total_log_lik(
    model$y_flat, model$n_sites, model$n_seasons, 1L, model$K_max,
    as.vector(Xl %*% bl), stats::qlogis(pd),
    as.vector(Xo %*% bo), as.vector(Xg %*% bg), use_nb = FALSE, eta_logr = 0)

  # Cheap central FD of the distance cell probs wrt sigma (distance integral only).
  h    <- 1e-5 * pmax(sigma, 1)
  cpp1 <- .gdr_dist_cp(sigma + h, model$cutpoints, model$transect)
  cpm1 <- .gdr_dist_cp(sigma - h, model$cutpoints, model$transect)
  dcpd <- (cpp1 - cpm1) / (2 * h)                    # [n x Jbin]
  dpd  <- rowSums(dcpd)                              # dpdist/dsigma

  # HMM contribution to dL/dsigma via eta_p = logit(pdist).
  dL_dsig <- ev$grad_eta_p / (pd * (1 - pd)) * dpd
  # Band contribution: sum_b Y_b d(log cpd_b)/dsigma - Ntot d(log pdist)/dsigma,
  # with Y_b the band totals across periods and Ntot the total detected per site.
  Yb   <- apply(model$y, c(1L, 2L), sum)             # [n x Jbin] band totals
  Ntot <- rowSums(Yb)
  dL_dsig <- dL_dsig +
    rowSums(Yb / pmax(cpd, 1e-300) * dcpd) - Ntot / pd * dpd

  g <- numeric(pl + ps + po + pg)
  g[seq_len(pl)]                <- as.numeric(crossprod(Xl, ev$grad_eta_lambda))
  g[pl + seq_len(ps)]           <- as.numeric(crossprod(Xs, dL_dsig * sigma))
  g[pl + ps + seq_len(po)]      <- as.numeric(crossprod(Xo, ev$grad_eta_omega))
  g[pl + ps + po + seq_len(pg)] <- as.numeric(crossprod(Xg, ev$grad_eta_gamma))
  g
}

# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# `y` is a 3D array [n_sites x n_bins x n_seasons] of per-distance-band counts at
# each primary period. abund_formula = log lambda, det_formula = log sigma,
# omega_formula = logit survival, gamma_formula = log recruitment.
.tobs_build_distsamp_open <- function(abund_formula, det_formula, omega_formula,
                                      gamma_formula, data, y, cutpoints, transect,
                                      K_max = NULL) {
  if (length(dim(y)) != 3L) {
    stop("distsamp_open() y must be a 3D array [n_sites x n_bins x n_seasons] of ",
         "per-distance-band counts.", call. = FALSE)
  }
  if (any(y < 0 | y != round(y), na.rm = TRUE)) {
    stop("distsamp_open() y must be non-negative integer counts.", call. = FALSE)
  }
  n_sites   <- dim(y)[1L]; n_bins <- dim(y)[2L]; n_seasons <- dim(y)[3L]
  if (n_seasons < 2L) {
    stop("distsamp_open() needs >= 2 primary periods (n_seasons = dim(y)[3]); ",
         "for a single period use distance().", call. = FALSE)
  }
  if (length(cutpoints) != n_bins + 1L) {
    stop(sprintf(paste0("distsamp_open() cutpoints must have length dim(y)[2] + 1 ",
         "= %d (the distance-bin edges)."), n_bins + 1L), call. = FALSE)
  }
  cutpoints <- as.numeric(cutpoints)
  if (any(diff(cutpoints) <= 0) || cutpoints[1] < 0) {
    stop("distsamp_open() cutpoints must be strictly increasing and start >= 0.",
         call. = FALSE)
  }
  .tobs_check_site_count(n_sites, nrow(data), "sites")
  storage.mode(y) <- "integer"

  # Per-period detected totals feed the HMM kernel (secondary occasions absorbed,
  # J = 1). Layout matches cpp_dyn_abun_total_log_lik: aperm(y[,1,,], c(2,3,1)).
  ntot   <- apply(y, c(1L, 3L), sum)                     # [n_sites x n_seasons]
  y_kern <- array(as.integer(ntot), c(n_sites, 1L, n_seasons))
  y_flat <- as.integer(aperm(y_kern, c(2L, 3L, 1L)))
  # The abundance-HMM forward is cubic in K_max, so the truncation is chosen from
  # the DETECTION-CORRECTED abundance scale, not a blunt multiple of the detected
  # total: max_i N_i ~ max(ntot) / pdist, plus a few Poisson SDs of headroom. A
  # rough pdist at the median-cutpoint sigma sets the scale (the fit re-integrates
  # the true pdist per iteration; the truncation only needs to bound N's tail).
  if (is.null(K_max)) {
    pd0  <- sum(.gdr_dist_cp(stats::median(cutpoints[-1]), cutpoints, transect))
    maxN <- max(ntot) / max(pd0, 0.1)
    K_max <- as.integer(ceiling(maxN + 4 * sqrt(maxN) + 10))
  }
  K_max <- as.integer(K_max)

  bind <- .tobs_bind_formulas(
    list(lambda = abund_formula, sigma = det_formula,
         omega = omega_formula, gamma = gamma_formula), data)
  X_lambda <- stats::model.matrix(bind$fe$lambda, data)
  X_sigma  <- stats::model.matrix(bind$fe$sigma, data)
  X_omega  <- stats::model.matrix(bind$fe$omega, data)
  X_gamma  <- stats::model.matrix(bind$fe$gamma, data)

  structure(list(
    model_type  = "distsamp_open",
    y           = y,
    y_flat      = y_flat,
    ntot        = ntot,
    cutpoints   = cutpoints,
    transect    = transect,
    n_bins      = n_bins,
    n_seasons   = n_seasons,
    K_max       = K_max,
    X_processes = list(X_lambda, X_sigma, X_omega, X_gamma),
    formulas    = list(lambda = bind$fe$lambda, sigma = bind$fe$sigma,
                       omega = bind$fe$omega, gamma = bind$fe$gamma),
    structured_terms = bind$terms,
    data        = data,
    n_sites     = n_sites,
    process_info = list(
      list(name = "lambda", p = ncol(X_lambda),
           coef_names = colnames(X_lambda), link = "log"),
      list(name = "sigma", p = ncol(X_sigma),
           coef_names = colnames(X_sigma), link = "log"),
      list(name = "omega", p = ncol(X_omega),
           coef_names = colnames(X_omega), link = "logit"),
      list(name = "gamma", p = ncol(X_gamma),
           coef_names = colnames(X_gamma), link = "log")
    )
  ), class = "tobs_model")
}

.dso_unpack <- function(theta, model) {
  Xl <- model$X_processes[[1L]]; Xs <- model$X_processes[[2L]]
  Xo <- model$X_processes[[3L]]; Xg <- model$X_processes[[4L]]
  pl <- ncol(Xl); ps <- ncol(Xs); po <- ncol(Xo); pg <- ncol(Xg)
  list(
    lambda = exp(as.vector(Xl %*% theta[seq_len(pl)])),
    sigma  = exp(as.vector(Xs %*% theta[pl + seq_len(ps)])),
    omega  = stats::plogis(as.vector(Xo %*% theta[pl + ps + seq_len(po)])),
    gamma  = exp(as.vector(Xg %*% theta[pl + ps + po + seq_len(pg)]))
  )
}

# ---------------------------------------------------------------------------
# Fitter
# ---------------------------------------------------------------------------

.tobs_fit_distsamp_open <- function(model, verbose = TRUE, ...) {
  ps <- vapply(model$process_info, function(pp) pp$p, integer(1))
  # Moment init: mean detected total seeds lambda given a rough detection; the
  # distance scale seeds at the median cutpoint; moderate survival / recruitment.
  ntot   <- model$ntot
  sig0   <- stats::median(model$cutpoints[-1])
  pdist0 <- sum(.gdr_dist_cp(sig0, model$cutpoints, model$transect))
  lam0   <- mean(ntot[, 1L]) / max(pdist0, 0.05)
  init <- c(log(max(lam0, 1e-2)), rep(0, ps[1] - 1L),
            log(max(sig0, 1e-2)),  rep(0, ps[2] - 1L),
            stats::qlogis(0.6),     rep(0, ps[3] - 1L),
            log(max(mean(ntot), 1)), rep(0, ps[4] - 1L))

  # Analytic gradient (one HMM kernel call per evaluation) drives BFGS; the
  # observed information is the FD-Jacobian of the negative gradient at the mode
  # (2*n_par cheap gradient calls), as in fp_occu() -- far cheaper than a numeric
  # Hessian of the value over the cubic-in-K forward recursion.
  fn  <- function(theta) .dso_negll(theta, model)
  ngr <- function(theta) -.dso_grad(theta, model)
  opt <- stats::optim(init, fn, gr = ngr, method = "BFGS",
                      control = list(maxit = 500L))
  converged <- opt$convergence == 0L

  par_names <- unlist(lapply(model$process_info, function(pp)
    paste0(pp$name, "_", pp$coef_names)))
  means <- opt$par; names(means) <- par_names
  Hobs <- .fp_fd_jacobian(ngr, opt$par)
  V <- tryCatch(solve(Hobs),
                error = function(e) diag(NA_real_, length(means)))
  V <- (V + t(V)) / 2
  dimnames(V) <- list(par_names, par_names)
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- par_names

  n_draws <- 1000L
  draws <- .occu_cover_rmvn(n_draws, means, V)
  colnames(draws) <- par_names

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

.dispatch_distsamp_open <- function(formula, data, family, detection, y, visits,
                                    engine, priors, control,
                                    approx = "gaussian_laplace",
                                    correction = "none", ...) {
  dots <- list(...)
  if (is.null(detection))
    stop("distsamp_open() requires a `detection` formula (the site-level ",
         "log-sigma distance-scale model).", call. = FALSE)
  if (is.null(y))
    stop("distsamp_open() requires `y` (a 3D array [n_sites x n_bins x ",
         "n_seasons] of per-distance-band counts).", call. = FALSE)
  if (!is.null(visits))
    stop("distsamp_open() detection is site-level; visit-level covariates ",
         "(`visits`) are not yet supported.", call. = FALSE)
  if (!identical(.map_engine(engine, family = "distsamp_open"), "laplace"))
    stop("distsamp_open() supports method = \"laplace\" only.", call. = FALSE)
  cutpoints <- family$params$cutpoints
  if (is.null(cutpoints))
    stop("distsamp_open() requires `cutpoints` (the distance-bin edges); pass ",
         "distsamp_open(cutpoints = ...).", call. = FALSE)
  model <- .tobs_build_distsamp_open(
    abund_formula = formula, det_formula = detection,
    omega_formula = dots$omega %||% ~1, gamma_formula = dots$gamma %||% ~1,
    data = data, y = y, cutpoints = cutpoints,
    transect = family$params$transect, K_max = family$params$K_max)
  .tobs_fit_distsamp_open(model, verbose = isTRUE(control$verbose))
}

# ---------------------------------------------------------------------------
# fitted / predict / residuals / pointwise log-likelihood / simulate
# ---------------------------------------------------------------------------

# fitted(): per-site initial abundance lambda, distance scale sigma, survival
# omega, recruitment gamma, and the overall distance detection pdist.
.tobs_fitted_distsamp_open <- function(object) {
  model <- object$model
  up <- .dso_unpack(object$means, model)
  pdist <- rowSums(.gdr_dist_cp(up$sigma, model$cutpoints, model$transect))
  list(lambda = up$lambda, sigma = up$sigma, omega = up$omega, gamma = up$gamma,
       pdist = pdist)
}

# predict(): abundance (default), distance scale, survival, or recruitment.
.tobs_predict_distsamp_open <- function(object, newdata = NULL,
                                        type = c("abundance", "distance",
                                                 "survival", "recruitment")) {
  type  <- match.arg(type)
  model <- object$model; b <- object$means
  p <- vapply(model$process_info, function(pp) pp$p, integer(1))
  off <- cumsum(c(0L, p))
  arm <- switch(type, abundance = 1L, distance = 2L, survival = 3L,
                recruitment = 4L)
  X <- if (is.null(newdata)) model$X_processes[[arm]]
       else stats::model.matrix(model$formulas[[arm]], newdata)
  eta <- as.vector(X %*% b[off[arm] + seq_len(p[arm])])
  if (identical(type, "survival")) stats::plogis(eta) else exp(eta)
}

# residuals(): per-site Pearson / deviance on the first-period detected total
# against its expected value lambda * pdist.
.tobs_residuals_distsamp_open <- function(object, type) {
  fv   <- .tobs_fitted_distsamp_open(object)
  m1   <- fv$lambda * fv$pdist
  obs  <- object$model$ntot[, 1L]
  eps  <- 1e-10
  res <- switch(type,
    response = obs - m1,
    pearson  = (obs - m1) / sqrt(m1 + eps),
    deviance = sign(obs - m1) * sqrt(2 * abs(
      ifelse(obs > 0, obs * log(obs / pmax(m1, eps)), 0) - (obs - m1))))
  list(occ = res, det = NULL)
}

# Pointwise (per-site) log-likelihood [n_draws x n_sites] over posterior draws.
.tobs_ploglik_distsamp_open <- function(object, n.draws = 1000L, n.threads = 1L) {
  model <- object$model
  draws <- object$draws
  if (!is.null(n.draws) && as.integer(n.draws) < nrow(draws)) {
    draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
  }
  t(vapply(seq_len(nrow(draws)), function(d)
    .dso_site_loglik(draws[d, ], model), numeric(model$n_sites)))
}

# Per-site log-likelihood vector at one coefficient draw (for WAIC / ploglik).
.dso_site_loglik <- function(theta, model) {
  up    <- .dso_unpack(theta, model)
  cpd   <- .gdr_dist_cp(up$sigma, model$cutpoints, model$transect)
  pdist <- rowSums(cpd)
  pdist <- pmin(pmax(pdist, 1e-10), 1 - 1e-10)
  pl <- length(up$lambda)
  ev <- cpp_dyn_abun_total_log_lik(
    model$y_flat, model$n_sites, model$n_seasons, 1L, model$K_max,
    log(up$lambda), stats::qlogis(pdist), stats::qlogis(up$omega),
    log(up$gamma), use_nb = FALSE, eta_logr = 0)
  pid  <- cpd / pdist
  band <- numeric(model$n_sites)
  for (t in seq_len(model$n_seasons)) {
    yb   <- model$y[, , t]
    band <- band + lgamma(rowSums(yb) + 1) - rowSums(lgamma(yb + 1)) +
            rowSums(yb * log(pmax(pid, 1e-300)))
  }
  as.numeric(ev$log_lik) + band
}

# Posterior replicate distance-bin arrays: draw a coefficient vector, then per
# site draw the open-population N sequence and the per-period distance obs.
.tobs_simulate_distsamp_open <- function(object, nsim = 1) {
  model <- object$model
  draw_one <- function() {
    idx <- sample.int(nrow(object$draws), 1L)
    up  <- .dso_unpack(object$draws[idx, ], model)
    .dso_draw(up$lambda, up$sigma, up$omega, up$gamma, model$cutpoints,
              model$transect, model$n_seasons)
  }
  if (nsim == 1L) return(draw_one())
  lapply(seq_len(nsim), function(s) draw_one())
}

# Core draw: per site the Dail-Madsen N sequence + per-period distance obs.
.dso_draw <- function(lambda, sigma, omega, gamma, cutpoints, transect, T) {
  n   <- length(lambda)
  cpd <- .gdr_dist_cp(sigma, cutpoints, transect); pdist <- rowSums(cpd)
  Jb  <- ncol(cpd)
  N   <- matrix(0L, n, T); N[, 1L] <- stats::rpois(n, lambda)
  for (t in 2:T)
    N[, t] <- stats::rbinom(n, N[, t - 1L], omega) + stats::rpois(n, gamma)
  y <- array(0L, c(n, Jb, T))
  for (i in seq_len(n)) for (t in seq_len(T)) {
    det <- stats::rbinom(1L, N[i, t], pdist[i])
    if (det > 0L)
      y[i, , t] <- as.integer(stats::rmultinom(1L, det, cpd[i, ] / pdist[i]))
  }
  y
}

# ---------------------------------------------------------------------------
# Simulator for recovery tests
# ---------------------------------------------------------------------------

#' Simulate an open-population distance-sampling data set
#'
#' Draws from the [distsamp_open()] model: an open metapopulation
#' (`N_1 ~ Poisson(lambda)`, `N_t = Binomial(N_{t-1}, omega) + Poisson(gamma)`)
#' observed by distance sampling (half-normal key, scale `sigma`) at each primary
#' period.
#'
#' @param N Number of sites (default 200).
#' @param cutpoints Distance-bin edges `0 = c_0 < ... < c_B`.
#' @param n_seasons Number of primary periods (default 4).
#' @param transect `"line"` (default) or `"point"`.
#' @param n_abund_covs,n_det_covs Number of abundance / distance covariates.
#' @param beta_lambda,beta_sigma Coefficients on the log-abundance and log-scale
#'   arms. Defaults give moderate abundance / detection.
#' @param omega,gamma Apparent survival probability and recruitment rate
#'   (intercept-only defaults 0.7 / 2.5).
#' @param seed Optional random seed.
#' @return A list with `y` (`[n_sites x n_bins x n_seasons]` distance-band
#'   counts), `data`, and `truth`.
#' @export
simulate_distsamp_open <- function(N = 200, cutpoints = c(0, 10, 20, 30, 40),
                                   n_seasons = 4L, transect = "line",
                                   n_abund_covs = 1, n_det_covs = 1,
                                   beta_lambda = NULL, beta_sigma = NULL,
                                   omega = 0.7, gamma = 2.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(beta_lambda))
    beta_lambda <- c(log(15), stats::runif(n_abund_covs, -0.3, 0.3))
  if (is.null(beta_sigma))
    beta_sigma <- c(log(stats::median(cutpoints[-1])),
                    stats::runif(n_det_covs, -0.2, 0.2))

  mk <- function(k, tag) {
    d <- data.frame(matrix(stats::rnorm(N * k), N, k))
    names(d) <- paste0(tag, seq_len(k)); d
  }
  ac <- mk(n_abund_covs, "abund_cov"); dc <- mk(n_det_covs, "det_cov")
  data <- cbind(ac, dc)
  lambda <- exp(as.vector(stats::model.matrix(~ ., ac) %*% beta_lambda))
  sigma  <- exp(as.vector(stats::model.matrix(~ ., dc) %*% beta_sigma))

  y <- .dso_draw(lambda, sigma, rep(omega, N), rep(gamma, N),
                 as.numeric(cutpoints), transect, as.integer(n_seasons))
  dimnames(y) <- list(NULL, paste0("band", seq_len(dim(y)[2])),
                      paste0("period", seq_len(n_seasons)))
  list(y = y, data = data,
       truth = list(beta_lambda = beta_lambda, beta_sigma = beta_sigma,
                    omega = omega, gamma = gamma, lambda = lambda, sigma = sigma))
}
