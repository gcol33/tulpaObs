# =============================================================================
# ms_count.R - community / multispecies relative-abundance GLMM (spAbundance
# msAbund). The community analogue of count(): per-species GLMM on an observed
# count / continuous response, no detection and no latent state.
#
#   y_{s,i}        ~ Poisson / NB / Gaussian(mu_{s,i})
#   g(mu_{s,i})    = X_i . (mu_beta + b_s)
#   b_s ~ N(0, Sigma_beta)                                   (community RE)
#
# The response is observed directly, so the per-species contribution is a plain
# GLMM log-likelihood (no marginalisation). The per-species coefficient
# deviations b_s are the random effects, fit by the shared community Laplace-EM
# (.tobs_community_em in R/community_em.R) -- exactly as ms_occu does, with a
# count sp_ll / sp_grad instead of the occupancy marginal. Negbin carries a
# per-species dispersion RE (log_r_s ~ N(mu_log_r, sigma_log_r)) as a second arm
# (matching ms_abun); Gaussian carries a per-species residual variance estimated
# in an outer loop. Non-spatial Laplace. y is [n_sites x n_species] (or a named
# list of n_species vectors), one value per site per species.
# =============================================================================


# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

.tobs_build_ms_count <- function(formula, data, y, species, response = "poisson",
                                 trials = NULL, structured_terms = list()) {
  # "bernoulli" is the jsdm() front door: an observed presence/absence community
  # GLMM (the spOccupancy lfJSDM / sfJSDM model class). It is the same community
  # model as msAbund -- per-species coefficients with a Gaussian community
  # covariance, no detection, no latent state -- with a logit link, so it shares
  # this binder, the community EM, the latent driver, and every S3 method rather
  # than carrying a parallel implementation (gcol33/tulpaObs#121). "binomial" is
  # the k-of-n generalization (community svcPGBinom, gcol33/tulpaObs#125): the
  # same logit-link path with a per-(site, species) trial count; "bernoulli" is
  # its trials = 1 special case (and stays the jsdm() alias).
  response <- match.arg(response,
                        c("poisson", "negbin", "gaussian", "bernoulli",
                          "binomial"))

  # y -> [n_sites x n_species] matrix. A named list of vectors becomes columns.
  to_mat <- function(z) {
    if (is.list(z) && !is.data.frame(z) && !is.array(z)) {
      m <- do.call(cbind, lapply(z, as.numeric))
      colnames(m) <- names(z)
      return(m)
    }
    z <- as.matrix(z)
    storage.mode(z) <- "double"
    z
  }
  y <- to_mat(y)
  n_sites   <- nrow(y)
  n_species <- ncol(y)

  species_names <- if (is.character(species)) species
                   else if (!is.null(colnames(y))) colnames(y)
                   else paste0("sp", seq_len(n_species))
  if (length(species_names) != n_species) {
    stop(sprintf("species has %d names but y has %d species (columns)",
                 length(species_names), n_species), call. = FALSE)
  }
  .tobs_check_site_count(n_sites, if (is.data.frame(data)) nrow(data) else n_sites,
                         "sites")

  X <- stats::model.matrix(formula, data)
  if (nrow(X) != n_sites) {
    stop(sprintf("design has %d rows but y has %d sites.", nrow(X), n_sites),
         call. = FALSE)
  }

  is_count_fam <- response %in% c("poisson", "negbin")
  is_binary    <- identical(response, "bernoulli")
  is_binom     <- identical(response, "binomial")
  valid <- !is.na(y)
  if (is_count_fam) {
    yv <- y[valid]
    if (any(yv < 0) || any(abs(yv - round(yv)) > 1e-8)) {
      stop(sprintf(paste0("ms_count(response = \"%s\"): y must be non-negative ",
                          "integer counts."), response), call. = FALSE)
    }
  }
  if (is_binary) {
    yv <- y[valid]
    if (any(yv != 0 & yv != 1)) {
      stop("jsdm(): y must contain only 0, 1, or NA (observed presence / ",
           "absence).", call. = FALSE)
    }
  }

  # Per-(site, species) trial count for the binomial response. A scalar recycles;
  # a length-n_sites vector is per-site (shared across species); a matrix is
  # per-(site, species). Bernoulli is trials = 1 (never carries a trial count).
  n_trials_mat <- NULL
  if (is_binom) {
    if (is.null(trials)) trials <- 1L
    if (is.matrix(trials)) {
      n_trials_mat <- trials
    } else if (length(trials) == 1L) {
      n_trials_mat <- matrix(as.numeric(trials), n_sites, n_species)
    } else if (length(trials) == n_sites) {
      n_trials_mat <- matrix(as.numeric(trials), n_sites, n_species)
    } else {
      stop("ms_count(response = \"binomial\"): `trials` must be a scalar, a ",
           "length-n_sites vector, or an n_sites x n_species matrix.",
           call. = FALSE)
    }
    if (nrow(n_trials_mat) != n_sites || ncol(n_trials_mat) != n_species) {
      stop("ms_count(response = \"binomial\"): the `trials` matrix must be ",
           "n_sites x n_species.", call. = FALSE)
    }
    storage.mode(n_trials_mat) <- "double"
    ntv <- n_trials_mat[valid]
    if (any(ntv < 1 | abs(ntv - round(ntv)) > 1e-8)) {
      stop("ms_count(response = \"binomial\"): `trials` must be positive ",
           "integers.", call. = FALSE)
    }
    yv <- y[valid]
    if (any(yv < 0) || any(abs(yv - round(yv)) > 1e-8)) {
      stop("ms_count(response = \"binomial\"): y must be non-negative integer ",
           "success counts.", call. = FALSE)
    }
    if (any(yv > ntv)) {
      stop("ms_count(response = \"binomial\"): every success count must be <= ",
           "its trial count (0 <= k <= n).", call. = FALSE)
    }
  }

  link <- switch(response, gaussian = "identity",
                 bernoulli = , binomial = "logit", "log")

  # Per-species valid rows + design + response, so each sp_ll skips NA sites. The
  # binomial arm also carries that species' per-site trial counts.
  summaries <- lapply(seq_len(n_species), function(s) {
    v  <- valid[, s]
    ys <- if (is_count_fam || is_binary || is_binom) as.numeric(round(y[v, s]))
          else as.numeric(y[v, s])
    su <- list(y = ys, X = X[v, , drop = FALSE], valid = v, n = sum(v))
    if (is_binom) su$n_trials <- as.numeric(round(n_trials_mat[v, s]))
    su
  })

  structure(list(
    model_type    = "ms_count",
    y             = y,
    n_trials      = n_trials_mat,
    valid         = valid,
    response      = response,
    link          = link,
    n_sites       = n_sites,
    n_species     = n_species,
    species_names = species_names,
    X             = X,
    summaries     = summaries,
    structured_terms = structured_terms,
    formulas      = list(mu = formula),
    data          = data,
    process_info  = list(
      list(name = "mu", p = ncol(X), coef_names = colnames(X), link = link)
    )
  ), class = "tobs_model")
}


# ---------------------------------------------------------------------------
# Per-species GLMM log-likelihood + analytic gradient kernels
# ---------------------------------------------------------------------------

# Poisson: theta = beta_s. log mu = X beta.
.ms_count_ll_pois <- function(su, beta) {
  eta <- as.numeric(su$X %*% beta)
  mu  <- exp(pmin(eta, 700))
  sum(stats::dpois(su$y, mu, log = TRUE))
}
.ms_count_grad_pois <- function(su, beta) {
  mu <- exp(pmin(as.numeric(su$X %*% beta), 700))
  as.numeric(crossprod(su$X, su$y - mu))
}

# Bernoulli (logit): theta = beta_s. logit psi = X beta. The jsdm() response --
# presence/absence observed directly, so this is a plain Bernoulli GLMM
# log-likelihood with no latent state to marginalise.
.ms_count_ll_bern <- function(su, beta) {
  eta <- as.numeric(su$X %*% beta)
  sum(ifelse(su$y > 0, stats::plogis(eta, log.p = TRUE),
                       stats::plogis(-eta, log.p = TRUE)))
}
.ms_count_grad_bern <- function(su, beta) {
  psi <- stats::plogis(as.numeric(su$X %*% beta))
  as.numeric(crossprod(su$X, su$y - psi))
}

# Binomial (logit): theta = beta_s, per-site trial count n_i in su$n_trials.
# y_i successes out of n_i; the Bernoulli kernel is n_i == 1. score = X'(y - n*p).
.ms_count_ll_binom <- function(su, beta) {
  p <- stats::plogis(as.numeric(su$X %*% beta))
  sum(stats::dbinom(su$y, size = su$n_trials,
                    prob = pmin(pmax(p, 1e-12), 1 - 1e-12), log = TRUE))
}
.ms_count_grad_binom <- function(su, beta) {
  p <- stats::plogis(as.numeric(su$X %*% beta))
  as.numeric(crossprod(su$X, su$y - su$n_trials * p))
}

# Gaussian (identity): theta = beta_s, per-species residual variance phi.
.ms_count_ll_gauss <- function(su, beta, phi) {
  mu <- as.numeric(su$X %*% beta)
  sum(stats::dnorm(su$y, mu, sqrt(max(phi, 1e-8)), log = TRUE))
}
.ms_count_grad_gauss <- function(su, beta, phi) {
  mu <- as.numeric(su$X %*% beta)
  as.numeric(crossprod(su$X, su$y - mu)) / max(phi, 1e-8)
}

# Negbin (NB2): theta = c(beta_s, log_r_s). Var = mu + mu^2 / r.
.ms_count_ll_nb <- function(su, beta, log_r) {
  r  <- exp(min(log_r, 30))
  mu <- exp(pmin(as.numeric(su$X %*% beta), 700))
  sum(stats::dnbinom(su$y, size = r, mu = pmax(mu, 1e-10), log = TRUE))
}
.ms_count_grad_nb <- function(su, beta, log_r) {
  r  <- exp(min(log_r, 30))
  mu <- exp(pmin(as.numeric(su$X %*% beta), 700))
  mu <- pmax(mu, 1e-10)
  # d/deta = r (y - mu) / (r + mu);  eta = log mu
  s_eta  <- r * (su$y - mu) / (r + mu)
  g_beta <- as.numeric(crossprod(su$X, s_eta))
  # d/dr = digamma(y+r) - digamma(r) + log(r/(r+mu)) + 1 - (y+r)/(r+mu)
  dLL_dr <- digamma(su$y + r) - digamma(r) + log(r / (r + mu)) + 1 -
            (su$y + r) / (r + mu)
  g_logr <- r * sum(dLL_dr)                     # chain rule log_r -> r
  c(g_beta, g_logr)
}


# ---------------------------------------------------------------------------
# Laplace-EM fitter (shared community engine)
# ---------------------------------------------------------------------------

# Run the shared community Laplace-EM for the response family and return the raw
# EM output alongside the arm layout / dispersion summary. The single source of
# the EM setup for both the Laplace front door (.tobs_fit_ms_count) and the NUTS
# warm start (.tobs_fit_ms_count_nuts): the NUTS path needs the raw fit (the
# per-species deviations, including the negbin log_r arm the Laplace summary
# drops) to pack its initial position.
.tobs_ms_count_run_em <- function(model, priors = NULL, max.iter = 200L,
                                  tol = 1e-4, sigma.beta = 5, verbose = TRUE,
                                  newton.max = 30L) {
  response <- model$response %||% "poisson"
  su       <- model$summaries
  S        <- model$n_species
  P_beta   <- model$process_info[[1L]]$p
  is_log   <- identical(model$link, "log")

  # Warm start: intercept at the pooled mean on the link scale, slopes 0. For the
  # binomial response the pooled mean is a PROPORTION (successes / trials), not a
  # raw count, so the logit intercept starts sensibly.
  yv      <- model$y[model$valid]
  mbar    <- if (identical(response, "binomial"))
               mean(yv / pmax(model$n_trials[model$valid], 1)) else mean(yv)
  mu0     <- numeric(P_beta)
  mu0[1L] <- switch(model$link %||% "log",
                    log      = log(max(mbar, 0.1)),
                    logit    = stats::qlogis(min(max(mbar, 1e-3), 1 - 1e-3)),
                    identity = mbar)

  run_em <- function(P, arm_idx, sp_ll, sp_grad, init_mu) {
    .tobs_community_em(
      S = S, P = P, arm_idx = arm_idx, sp_ll = sp_ll, sp_grad = sp_grad,
      init_mu = init_mu, init_global = numeric(0),
      penalize_global = FALSE, sigma_beta = sigma.beta, priors = priors,
      sigma_init = 0.3, max_iter = as.integer(max.iter), tol = as.numeric(tol),
      newton_max = as.integer(newton.max), verbose = isTRUE(verbose))
  }

  if (identical(response, "poisson")) {
    arm_idx <- list(mu = seq_len(P_beta))
    sp_ll   <- function(s, theta, global) .ms_count_ll_pois(su[[s]], theta)
    sp_grad <- function(s, theta, global) .ms_count_grad_pois(su[[s]], theta)
    fit <- run_em(P_beta, arm_idx, sp_ll, sp_grad, mu0)
    disp <- NULL

  } else if (identical(response, "bernoulli")) {
    # jsdm(): observed presence/absence, no dispersion parameter.
    arm_idx <- list(mu = seq_len(P_beta))
    sp_ll   <- function(s, theta, global) .ms_count_ll_bern(su[[s]], theta)
    sp_grad <- function(s, theta, global) .ms_count_grad_bern(su[[s]], theta)
    fit <- run_em(P_beta, arm_idx, sp_ll, sp_grad, mu0)
    disp <- NULL

  } else if (identical(response, "binomial")) {
    # k-of-n binomial community GLMM (community svcPGBinom): logit link, the
    # trial count pins the variance, so no dispersion parameter.
    arm_idx <- list(mu = seq_len(P_beta))
    sp_ll   <- function(s, theta, global) .ms_count_ll_binom(su[[s]], theta)
    sp_grad <- function(s, theta, global) .ms_count_grad_binom(su[[s]], theta)
    fit <- run_em(P_beta, arm_idx, sp_ll, sp_grad, mu0)
    disp <- NULL

  } else if (identical(response, "gaussian")) {
    arm_idx <- list(mu = seq_len(P_beta))
    # Per-species residual variance, refined in an outer loop around the EM (the
    # coefficients depend on phi only through the RE shrinkage, so a few passes
    # converge). Captured by reference so sp_ll sees the current phi.
    phi <- rep(max(stats::var(yv), 1e-4), S)
    fit <- NULL
    for (outer in seq_len(15L)) {
      sp_ll   <- function(s, theta, global) .ms_count_ll_gauss(su[[s]], theta, phi[s])
      sp_grad <- function(s, theta, global) .ms_count_grad_gauss(su[[s]], theta, phi[s])
      fit <- run_em(P_beta, arm_idx, sp_ll, sp_grad, mu0)
      phi_new <- vapply(seq_len(S), function(s) {
        mu <- as.numeric(su[[s]]$X %*% (fit$mu + fit$b_list[[s]]))
        max(mean((su[[s]]$y - mu)^2), 1e-8)
      }, numeric(1))
      if (max(abs(log(phi_new) - log(phi))) < 1e-4) { phi <- phi_new; break }
      phi <- phi_new
    }
    disp <- list(response = "gaussian", variance = phi)

  } else { # negbin: per-species log_r as a second community RE arm
    P <- P_beta + 1L
    arm_idx <- list(mu = seq_len(P_beta), disp = P)
    sp_ll   <- function(s, theta, global)
      .ms_count_ll_nb(su[[s]], theta[seq_len(P_beta)], theta[P])
    sp_grad <- function(s, theta, global)
      .ms_count_grad_nb(su[[s]], theta[seq_len(P_beta)], theta[P])
    init_mu <- c(mu0, log(2))                    # community mean log_r
    fit <- run_em(P, arm_idx, sp_ll, sp_grad, init_mu)
    disp <- list(response = "negbin",
                 mu_log_r = fit$mu[P], sigma_log_r = sqrt(fit$Sigma$disp[1, 1]),
                 r_s = exp(vapply(fit$b_list, function(b) fit$mu[P] + b[P],
                                  numeric(1))))
  }

  list(fit = fit, arm_idx = arm_idx, disp = disp, response = response,
       P_beta = P_beta)
}


.tobs_fit_ms_count <- function(model,
                               priors     = NULL,
                               max.iter   = 200L,
                               tol        = 1e-4,
                               sigma.beta = 5,
                               verbose    = TRUE,
                               ...) {
  dots <- list(...)
  em   <- .tobs_ms_count_run_em(
    model, priors = priors, max.iter = max.iter, tol = tol,
    sigma.beta = sigma.beta, verbose = verbose,
    newton.max = as.integer(dots$newton.max %||% 30L))
  build_ms_count_fit(model, em$fit, em$arm_idx, em$disp)
}


# ---------------------------------------------------------------------------
# Wrap the EM output into a tobs_fit
# ---------------------------------------------------------------------------

build_ms_count_fit <- function(model, fit, arm_idx, disp = NULL) {
  pi_list  <- model$process_info
  cn       <- pi_list[[1L]]$coef_names
  P_beta   <- pi_list[[1L]]$p
  beta_idx <- seq_len(P_beta)

  # Report only the mean-coefficient arm as the community-mean coefficients; the
  # negbin log_r community mean rides fit$ms_dispersion.
  par_names <- paste0("mu_", cn)
  means <- fit$mu[beta_idx]; names(means) <- par_names
  V <- fit$Vf[beta_idx, beta_idx, drop = FALSE]
  dimnames(V) <- list(par_names, par_names); V <- (V + t(V)) / 2
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- par_names

  n_draws <- 1000L
  draws <- .occu_cover_rmvn(n_draws, means, V)
  colnames(draws) <- par_names

  B    <- do.call(rbind, fit$b_list)            # S x P
  blup <- B[, beta_idx, drop = FALSE]
  coef <- sweep(blup, 2L, means, "+")
  rownames(blup) <- rownames(coef) <- model$species_names
  colnames(blup) <- colnames(coef) <- cn

  Sigma_mu <- fit$Sigma$mu
  dimnames(Sigma_mu) <- list(cn, cn)

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = length(means),
    log_prob     = rep(fit$logML, n_draws),
    log_lik      = fit$logML,
    N            = sum(model$valid)),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names    = par_names,
    param_names  = par_names,
    n_fixed      = length(means),
    fixed_names  = par_names,
    process_info = pi_list,
    model        = model,
    spatial      = NULL,
    method       = "laplace",
    ms_community = list(
      Sigma_mu = Sigma_mu,
      sd_mu    = sqrt(pmax(diag(Sigma_mu), 0)),
      coef_mu  = coef, blup_mu = blup
    ),
    ms_dispersion = disp,
    convergence  = list(converged = isTRUE(fit$converged), n_iter = fit$n_iter)
  )), class = c("tobs_fit", "tulpa_fit"))
}


# ---------------------------------------------------------------------------
# S3 helpers (routed from methods.R by model_type == "ms_count")
# ---------------------------------------------------------------------------

# Response-scale mean from the linear predictor, per the family link. logit is
# the jsdm() (bernoulli) response; log the count responses; identity Gaussian.
.ms_count_linkinv <- function(model, eta) {
  switch(model$link %||% "log",
         log      = exp(eta),
         logit    = stats::plogis(eta),
         identity = eta)
}

.tobs_ranef_ms_count <- function(object) {
  .tobs_ranef_ms_long(object$ms_community, c(mu = "blup_mu"))
}

# Per-species fitted mean mu [n_sites x n_species] on the response scale. A
# shared areal field (nested_laplace path) adds its per-site offset to every
# species' predictor.
.tobs_fitted_ms_count <- function(object) {
  model <- object$model
  cm    <- object$ms_community
  eta   <- model$X %*% t(cm$coef_mu)            # n_sites x n_species
  fld   <- object$spatial_field
  if (!is.null(fld) && length(fld) == nrow(eta)) {
    eta <- sweep(eta, 1L, as.numeric(fld), "+")
  }
  # Latent factors add a per-(species, site) residual offset (eta lambda').
  fo <- model$count_factor_offset
  if (!is.null(fo) && all(dim(fo) == dim(eta))) eta <- eta + fo
  mu    <- .ms_count_linkinv(model, eta)
  # Binomial: the fitted quantity on the y-scale is expected successes n * p.
  if (identical(model$response %||% "poisson", "binomial") &&
      !is.null(model$n_trials)) {
    mu <- mu * model$n_trials
  }
  dimnames(mu) <- list(NULL, model$species_names)
  list(mu = mu)
}

# predict(): per-species mean [n_sites x n_species] on the response scale --
# probabilities for the jsdm() (bernoulli) response, counts for the others. In
# sample this is fitted()$mu, which carries the shared field / latent factor
# offsets. For newdata the design is rebuilt at the per-species coefficients; a
# new site has no field node and no factor score, so those offsets are dropped
# (field interpolation is a separate concern, as on the single-species count).
.tobs_predict_ms_count <- function(object, newdata = NULL) {
  if (is.null(newdata)) return(fitted(object)$mu)
  model <- object$model
  X   <- stats::model.matrix(model$formulas$mu, newdata)
  eta <- X %*% t(object$ms_community$coef_mu)
  mu  <- .ms_count_linkinv(model, eta)
  dimnames(mu) <- list(NULL, model$species_names)
  mu
}

# residuals(): per-(site, species) residuals against the fitted mean. Deviance
# (the default) uses each family's saturated-model form; Pearson scales by the
# family SD (bernoulli mu(1-mu), Poisson mu, negbin mu + mu^2/r, gaussian the
# per-species residual variance); response is the raw y - mu. NA cells of y stay
# NA. Mirrors the single-species .tobs_residuals_count.
.tobs_residuals_ms_count <- function(object,
                                     type = c("deviance", "pearson",
                                              "response")) {
  type     <- match.arg(type)
  model    <- object$model
  response <- model$response %||% "poisson"
  disp     <- object$ms_dispersion
  mu  <- fitted(object)$mu
  y   <- matrix(as.numeric(model$y), model$n_sites, model$n_species)
  mup <- pmax(mu, 1e-8)
  rs  <- if (identical(response, "negbin"))
           matrix(disp$r_s, nrow(mu), ncol(mu), byrow = TRUE) else NULL
  nt  <- if (identical(response, "binomial") && !is.null(model$n_trials))
           pmax(model$n_trials, 1) else NULL

  r <- switch(type,
    response = y - mu,
    pearson  = (y - mu) / sqrt(switch(response,
                 bernoulli = pmax(mu * (1 - mu), 1e-8),
                 binomial  = pmax(mu * (1 - mu / nt), 1e-8),
                 poisson   = mup,
                 negbin    = pmax(mu + mu^2 / rs, 1e-8),
                 gaussian  = matrix(pmax(disp$variance, 1e-8), nrow(mu),
                                    ncol(mu), byrow = TRUE))),
    deviance = switch(response,
      bernoulli = {
        p <- pmin(pmax(mu, 1e-10), 1 - 1e-10)
        sign(y - mu) * sqrt(pmax(-2 * (y * log(p) + (1 - y) * log1p(-p)), 0))
      },
      binomial = {
        t1  <- ifelse(y > 0, y * log(y / mup), 0)
        fy  <- nt - y; fmu <- pmax(nt - mu, 1e-8)
        t2  <- ifelse(fy > 0, fy * log(fy / fmu), 0)
        sign(y - mu) * sqrt(pmax(2 * (t1 + t2), 0))
      },
      poisson  = {
        term <- ifelse(y > 0, y * log(y / mup), 0)
        sign(y - mu) * sqrt(pmax(2 * (term - (y - mu)), 0))
      },
      negbin   = {
        term <- ifelse(y > 0, y * log(y / mup), 0)
        d <- 2 * (term - (y + rs) * log((y + rs) / (mup + rs)))
        sign(y - mu) * sqrt(pmax(d, 0))
      },
      gaussian = y - mu))
  dimnames(r) <- list(NULL, model$species_names)
  list(mu = r)
}

# Draw community count data under the fitted per-species coefficients.
.tobs_simulate_ms_count <- function(object, nsim = 1) {
  model <- object$model
  cm    <- object$ms_community
  n_sites <- model$n_sites; n_species <- model$n_species
  response <- model$response %||% "poisson"
  disp <- object$ms_dispersion
  one <- function() {
    out <- matrix(NA_real_, n_sites, n_species,
                  dimnames = list(NULL, model$species_names))
    for (s in seq_len(n_species)) {
      eta <- as.numeric(model$X %*% cm$coef_mu[s, ])
      mu  <- .ms_count_linkinv(model, eta)
      nts <- if (identical(response, "binomial") && !is.null(model$n_trials))
               as.integer(round(model$n_trials[, s])) else NULL
      out[, s] <- switch(response,
        poisson   = stats::rpois(n_sites, mu),
        bernoulli = stats::rbinom(n_sites, 1L, mu),
        binomial  = stats::rbinom(n_sites, size = nts, prob = mu),
        negbin    = stats::rnbinom(n_sites, size = disp$r_s[s], mu = mu),
        gaussian  = stats::rnorm(n_sites, mu, sqrt(disp$variance[s])))
    }
    out
  }
  if (nsim == 1L) one() else lapply(seq_len(nsim), function(i) one())
}


# Pointwise log-likelihood [n_draws x (n_sites * n_species)] for WAIC / LOO: the
# per-species deviations held at their BLUPs, the community-mean coefficients
# varying over the draws. One observation unit per (species, valid site).
.tobs_ploglik_ms_count <- function(object, n.draws = 1000L) {
  model <- object$model
  cm    <- object$ms_community
  draws <- object$draws
  if (!is.null(n.draws) && n.draws < nrow(draws)) {
    draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
  }
  response <- model$response %||% "poisson"
  is_log   <- identical(model$link, "log")
  disp     <- object$ms_dispersion
  # Shared areal field: a per-site offset added to every species' predictor. The
  # field per-species valid rows follow su$valid (all TRUE on the spatial path).
  # Latent factors add a per-(species, site) offset (eta lambda').
  fld_full <- object$spatial_field
  fac_off  <- object$model$count_factor_offset
  cols <- lapply(seq_len(model$n_species), function(s) {
    su <- model$summaries[[s]]
    if (su$n == 0L) return(NULL)
    eta <- draws %*% t(su$X) + matrix(as.numeric(su$X %*% cm$blup_mu[s, ]),
                                      nrow(draws), su$n, byrow = TRUE)
    if (!is.null(fld_full) && length(fld_full) == length(su$valid)) {
      eta <- eta + matrix(as.numeric(fld_full[su$valid]),
                          nrow(draws), su$n, byrow = TRUE)
    }
    if (!is.null(fac_off) && nrow(fac_off) == length(su$valid)) {
      eta <- eta + matrix(as.numeric(fac_off[su$valid, s]),
                          nrow(draws), su$n, byrow = TRUE)
    }
    mu  <- if (is_log) pmin(pmax(exp(eta), 1e-300), 1e8)
           else .ms_count_linkinv(model, eta)
    Y   <- matrix(su$y, nrow(draws), su$n, byrow = TRUE)
    switch(response,
      poisson   = stats::dpois(Y, pmax(mu, 1e-300), log = TRUE),
      # Bernoulli scored on the log scale from eta directly (plogis(log.p) is
      # stable where mu saturates at 0 / 1).
      bernoulli = ifelse(Y > 0, stats::plogis(eta, log.p = TRUE),
                                stats::plogis(-eta, log.p = TRUE)),
      binomial  = {
        NT <- matrix(su$n_trials, nrow(draws), su$n, byrow = TRUE)
        stats::dbinom(Y, size = NT,
                      prob = pmin(pmax(mu, 1e-12), 1 - 1e-12), log = TRUE)
      },
      negbin    = stats::dnbinom(Y, size = disp$r_s[s], mu = pmax(mu, 1e-8),
                                 log = TRUE),
      gaussian  = stats::dnorm(Y, mu, sqrt(max(disp$variance[s], 1e-8)),
                               log = TRUE))
  })
  do.call(cbind, cols)
}


# ---------------------------------------------------------------------------
# Family constructor
# ---------------------------------------------------------------------------

#' Multispecies (community) relative-abundance GLMM family
#'
#' The community analogue of [count()]: a per-species GLMM on an observed count
#' or continuous response with Gaussian community hyperpriors on the per-species
#' coefficients (the spAbundance `msAbund` model). No detection and no latent
#' state -- the abundance counterpart of [ms_occu()] without the occupancy layer.
#' Poisson / negative-binomial (log link) or Gaussian (identity). `y` is an
#' `n_sites x n_species` matrix (or a named list of `n_species` count vectors),
#' one value per site per species; `species` names the columns.
#'
#' As with the other community Laplace-EM families ([ms_occu()], [ms_occu_cover()]),
#' the community means are recovered essentially unbiased for the Gaussian
#' (identity link, exact Laplace) and Poisson responses; the negative-binomial
#' slope carries a mild first-order-Laplace (PQL) attenuation of order the
#' community variance (a few percent, shrinking with the number of species and
#' observations per species). The community-mean Wald intervals are calibrated to
#' the package rubric (pooled coverage at least 0.85).
#'
#' `method = "nuts"` samples the exact joint posterior (community means,
#' per-species deviations, and the community covariance) for all three responses
#' -- the negative binomial carrying a per-species dispersion random effect, the
#' Gaussian a per-species free residual variance -- which removes the
#' Laplace-EM's negative-binomial attenuation and returns calibrated,
#' non-Gaussian community intervals. A shared areal field or latent factors are
#' available through `method = "nested_laplace"` / the `latent()` term; see the
#' package overview for the spatial and factor variants.
#'
#' @param response One of `"poisson"`, `"negbin"`, `"gaussian"`, or
#'   `"binomial"`. The binomial response is the community `k`-of-`n` GLMM
#'   (community `svcPGBinom`): supply the per-site (or per-`site x species`)
#'   trial count as `trials =` on [tobs()] (default 1, i.e. Bernoulli, which is
#'   the [jsdm()] response). With `trials > 1` the community-mean intercept
#'   carries a small first-order-Laplace bias of order `1 / n_species` (a few
#'   hundredths on the logit scale at 20 species, shrinking with more species;
#'   the slope and the `trials = 1` case are unbiased) -- the same character as
#'   the negative-binomial slope attenuation noted above.
#' @return A `tobs_family` object.
#' @seealso [count()] (single species), [ms_occu()], [ms_abun()]
#' @examples
#' \donttest{
#' sim <- simulate_ms_count(N = 120, n_species = 8, seed = 1)
#' fit <- tobs(~ x, data = sim$data, family = ms_count(), y = sim$y,
#'             species = colnames(sim$y), method = "laplace")
#' summary(fit)
#' }
#' @export
ms_count <- function(response = c("poisson", "negbin", "gaussian",
                                  "binomial")) {
  response <- match.arg(response)
  obs_family(
    name           = "ms_count",
    class_long     = "multispecies count / relative-abundance GLMM",
    latent         = "none",
    observation    = response,
    replicates     = "single",
    default_engine = "laplace",
    status         = "working",
    params         = list(response = response),
    response       = "matrix"
  )
}
