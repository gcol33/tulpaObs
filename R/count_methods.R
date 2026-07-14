# count_methods.R - fitted() / predict() / residuals() for the count /
# relative-abundance GLMM. count() observes the response directly (no detection
# process, no latent state): the fitted quantity is the per-site mean mu on the
# response scale (exp(X beta) for the log-link Poisson / negbin, X beta for the
# identity-link Gaussian). The residual is the response residual y - mu.

# Per-site linear predictor eta at a design X [n x p].
.tobs_count_eta <- function(object, X) {
  beta <- object$means[seq_len(object$model$process_info[[1L]]$p)]
  as.numeric(X %*% beta)
}

# Response-scale mean from eta, per the family link.
.tobs_count_mu <- function(object, eta) {
  if (identical(object$model$link %||% "log", "log")) exp(eta) else eta
}

.tobs_fitted_count <- function(object) {
  mu <- .tobs_count_mu(object, .tobs_count_eta(object, object$model$X_occ))
  list(mu = mu)
}

# predict(): in-sample returns fitted()$mu; newdata recomputes mu from the design.
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

  r <- switch(type,
    response = y - mu,
    pearson  = (y - mu) / sqrt(switch(response,
                 poisson  = mup,
                 negbin   = mup + mup^2 / size,
                 gaussian = rep(max(phi, 1e-8), length(mu)))),
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
      gaussian = y - mu))
  list(mu = r, det = NULL)
}
