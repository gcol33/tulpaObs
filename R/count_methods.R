# count_methods.R - fitted() / predict() / residuals() for the count /
# relative-abundance GLMM. count() observes the response directly (no detection
# process, no latent state): the fitted quantity is the per-site mean mu on the
# response scale (exp(X beta) for the log-link Poisson / negbin, X beta for the
# identity-link Gaussian). The residual is the response residual y - mu.

# Per-site linear predictor eta at a design X [n x p]. `add_field = TRUE` adds
# the areal latent field (one node per site) for an in-sample nested-Laplace
# fit; the log-mean per site is then X beta + sum_k W[i,k] f_k[i]. The fit
# records that per-site sum as `count_field_offset` -- for a varying-coefficient
# (SVC) field the contribution is the weighted sum over the intercept and trend
# fields, not the intercept field alone. `spatial_field` is the fallback for a
# fit made before the offset was recorded (a plain intercept field, weight 1).
# The field is skipped for newdata (a new site has no field node -- field
# interpolation is a separate concern).
.tobs_count_eta <- function(object, X, add_field = FALSE) {
  beta <- object$means[seq_len(object$model$process_info[[1L]]$p)]
  eta  <- as.numeric(X %*% beta)
  if (isTRUE(add_field)) {
    fld <- object$model$count_field_offset %||% object$spatial_field
    if (!is.null(fld) && length(fld) == length(eta)) eta <- eta + as.numeric(fld)
  }
  eta
}

# Response-scale mean from eta, per the family link. For the log link this is the
# expected count exp(eta); for the identity link eta; for the binomial logit link
# it is the per-trial success PROBABILITY plogis(eta) (the fitted / residual code
# scales it by the trial count to expected successes where the trials are known).
.tobs_count_mu <- function(object, eta) {
  link <- object$model$link %||% "log"
  if (identical(link, "log")) exp(eta)
  else if (identical(link, "logit")) stats::plogis(eta)
  else eta
}

.tobs_fitted_count <- function(object) {
  eta <- .tobs_count_eta(object, object$model$X_occ, add_field = TRUE)
  mu  <- .tobs_count_mu(object, eta)
  # Binomial: the fitted quantity on the y-scale is the expected number of
  # successes n_i * p_i, so residuals compare the observed successes to it.
  if (identical(object$model$link %||% "log", "logit")) {
    nt <- as.numeric(object$model$n_trials %||% rep(1, length(mu)))
    mu <- mu * nt
  }
  list(mu = mu)
}

# predict(): in-sample returns fitted()$mu (field-aware); newdata recomputes mu
# from the design at the fixed effects only (no field node at a new site). For a
# binomial fit a new site has no trial count, so predict returns the per-trial
# success probability plogis(eta) (multiply by your own trials for expected
# successes); the in-sample fitted()$mu is expected successes.
.tobs_predict_count <- function(object, newdata = NULL) {
  if (is.null(newdata)) return(fitted(object)$mu)
  X <- stats::model.matrix(object$model$formulas$occ, newdata)
  .tobs_count_mu(object, .tobs_count_eta(object, X))
}

# Poisson / negative-binomial residual on a count response, given the fitted
# mean and the estimated NB size (`size = Inf` is Poisson). The single source of
# truth for the count families that score a thinned-abundance mean: count()
# below, and the latent-N families whose per-cell marginal is the abundance
# distribution thinned by a detection probability -- abun() (Binomial(N, p)),
# removal() (depleting binomial) and distance() (bin multinomial). The negative
# binomial is closed under binomial thinning with the SAME size, so a per-cell
# marginal of an NB(r, lambda) abundance is NB(r, lambda * pi) and the variance
# is mu + mu^2/r; scoring it at the Poisson variance inflates every residual by
# exactly the overdispersion the fit estimated. `eps` floors the mean inside the
# logs and the variance only; `mu` itself is used for y - mu and its sign, so a
# caller that already floored its mean passes the same value for both.
.tobs_count_residual <- function(y, mu, type, size = Inf, eps = 1e-8) {
  mup <- pmax(mu, eps)
  switch(type,
    response = y - mu,
    pearson  = (y - mu) / sqrt(if (is.finite(size)) mup + mup^2 / size else mup),
    deviance = {
      term <- ifelse(y > 0, y * log(y / mup), 0)
      d <- if (is.finite(size))
             2 * (term - (y + size) * log((y + size) / (mup + size)))
           else 2 * (term - (y - mu))
      sign(y - mu) * sqrt(pmax(d, 0))
    })
}

# simulate(): a replicate response drawn at one posterior beta draw (a random
# row of object$draws, the .tobs_simulate_distance() / .tobs_simulate_nmix()
# convention) under the family's own distribution -- count() has no detection
# process and no latent state, so this IS its whole posterior predictive.
# Reuses .tobs_count_eta() / .tobs_count_mu() (field offset + link) against a
# lightweight stand-in carrying that one draw as `means`, rather than
# R/sbc.R's .tobs_sbc_replicate_count() duplicating the mean/field/link
# computation inline (which it did without the field offset, and read a
# negbin dispersion column `object$draws` for count() has never carried).
# Drives test_dispersion() / test_zero_inflation() / test_outliers() via
# .tobs_count_gof_families (R/diagnostics_count_gof.R).
.tobs_simulate_count <- function(object, nsim = 1) {
  model    <- object$model
  response <- model$response %||% "poisson"
  draws    <- object$draws
  n_draws  <- nrow(draws)
  p        <- model$process_info[[1L]]$p
  n        <- nrow(model$X_occ)
  disp     <- object$count_dispersion
  one <- function() {
    beta <- draws[sample.int(n_draws, 1L), seq_len(p)]
    fake <- list(model = model, means = beta,
                spatial_field = object$spatial_field)
    eta  <- .tobs_count_eta(fake, model$X_occ, add_field = TRUE)
    mu   <- .tobs_count_mu(fake, eta)
    if (identical(model$link %||% "log", "logit")) {
      nt <- as.numeric(model$n_trials %||% rep(1, n)); mu <- mu * nt
    }
    switch(response,
      poisson  = stats::rpois(n, mu),
      negbin   = stats::rnbinom(n, mu = mu, size = disp$phi %||% Inf),
      gaussian = stats::rnorm(n, mu, sqrt(disp$phi %||% 1)),
      binomial = {
        nt <- as.numeric(model$n_trials %||% rep(1, n))
        stats::rbinom(n, size = round(nt), prob = mu / pmax(nt, 1))
      },
      stop("simulate() for count() is not implemented for response = '",
           response, "'.", call. = FALSE))
  }
  if (nsim == 1L) one() else lapply(seq_len(nsim), function(i) one())
}

# residuals(): deviance (default), Pearson, or response (raw y - mu), per the
# family. A count GLMM has one series, one value per response row: it fills the
# unit-level `occ` slot of the residuals() contract and leaves `det` NULL. Poisson / negbin deviance uses the standard GLM saturated-model form;
# the Gaussian deviance residual is the raw y - mu.
.tobs_residuals_count <- function(object, type = c("deviance", "pearson",
                                                   "response")) {
  type     <- match.arg(type)
  response <- object$model$response %||% "poisson"
  mu   <- fitted(object)$mu
  y    <- as.numeric(object$model$y_count)
  size <- object$count_dispersion$phi %||% Inf   # negbin size
  phi  <- object$count_dispersion$phi %||% 1      # gaussian variance
  mup  <- pmax(mu, 1e-8)
  # Binomial trial counts (expected successes mu = n * p, so p = mu / n).
  nt   <- as.numeric(object$model$n_trials %||% rep(1, length(mu)))
  ntp  <- pmax(nt, 1)

  # The Poisson / negbin branches are .tobs_count_residual(); gaussian and
  # binomial are count()'s own responses and stay here.
  pois_nb <- function() .tobs_count_residual(
    y, mu, type, size = if (identical(response, "negbin")) size else Inf)

  r <- switch(type,
    response = y - mu,
    pearson  = switch(response,
                 poisson  = pois_nb(),
                 negbin   = pois_nb(),
                 gaussian = (y - mu) / sqrt(rep(max(phi, 1e-8), length(mu))),
                 binomial = (y - mu) / sqrt(pmax(mu * (1 - mu / ntp), 1e-8))),
    deviance = switch(response,
      poisson  = pois_nb(),
      negbin   = pois_nb(),
      gaussian = y - mu,
      binomial = {
        # Standard binomial deviance residual with n - y failures; the
        # 0 * log(0/.) terms are taken as 0.
        t1 <- ifelse(y > 0, y * log(y / mup), 0)
        fy <- nt - y; fmu <- pmax(nt - mu, 1e-8)
        t2 <- ifelse(fy > 0, fy * log(fy / fmu), 0)
        sign(y - mu) * sqrt(pmax(2 * (t1 + t2), 0))
      }))
  list(occ = r, det = NULL)
}
