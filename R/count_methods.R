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

# residuals(): deviance (default), Pearson, or response (raw y - mu), per the
# family. Poisson / negbin deviance uses the standard GLM saturated-model form;
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

  r <- switch(type,
    response = y - mu,
    pearson  = (y - mu) / sqrt(switch(response,
                 poisson  = mup,
                 negbin   = mup + mup^2 / size,
                 gaussian = rep(max(phi, 1e-8), length(mu)),
                 binomial = pmax(mu * (1 - mu / ntp), 1e-8))),
    deviance = switch(response,
      poisson  = {
        term <- ifelse(y > 0, y * log(y / mup), 0)
        sign(y - mu) * sqrt(pmax(2 * (term - (y - mu)), 0))
      },
      negbin   = {
        term <- ifelse(y > 0, y * log(y / mup), 0)
        d <- 2 * (term - (y + size) * log((y + size) / (mup + size)))
        sign(y - mu) * sqrt(pmax(d, 0))
      },
      gaussian = y - mu,
      binomial = {
        # Standard binomial deviance residual with n - y failures; the
        # 0 * log(0/.) terms are taken as 0.
        t1 <- ifelse(y > 0, y * log(y / mup), 0)
        fy <- nt - y; fmu <- pmax(nt - mu, 1e-8)
        t2 <- ifelse(fy > 0, fy * log(fy / fmu), 0)
        sign(y - mu) * sqrt(pmax(2 * (t1 + t2), 0))
      }))
  list(mu = r, det = NULL)
}
