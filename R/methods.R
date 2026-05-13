# ============================================================================
# Occupancy-specific S3 methods for TulpaObs_fit objects
# Generic S3 (coef, confint, vcov, logLik, summary, tidy, glance, ranef, plot)
# are inherited from tulpa::tulpa_fit via class = c("TulpaObs_fit", "tulpa_fit")
# ============================================================================

#' Number of observations
#' @param object A `TulpaObs_fit` object.
#' @param ... Ignored.
#' @return Integer count of non-NA detection history entries.
#' @export
nobs.TulpaObs_fit <- function(object, ...) {
  model <- object$model
  if (model$model_type == "single" || model$model_type == "community") {
    y <- model$y
    sum(y >= 0)
  } else if (model$model_type == "dynamic") {
    sum(model$y_flat >= 0)
  } else {
    NA_integer_
  }
}

#' Fitted values (occupancy and detection probabilities)
#' @param object A `TulpaObs_fit` object.
#' @param ... Ignored.
#' @return A list with `psi` (occupancy probabilities), `p` (detection probabilities),
#'   and `z` (posterior P(z=1)) at posterior mean.
#' @export
fitted.TulpaObs_fit <- function(object, ...) {
  model <- object$model
  means <- object$means
  pi_list <- model$process_info

  # Extract occupancy and detection linear predictors at posterior mean
  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  beta_occ <- means[seq_len(pi_list[[1]]$p)]
  beta_det <- means[pi_list[[1]]$p + seq_len(pi_list[[2]]$p)]

  eta_occ <- as.vector(X_occ %*% beta_occ)
  eta_det <- as.vector(X_det %*% beta_det)

  psi <- plogis(eta_occ)
  p <- plogis(eta_det)

  # Compute z posterior: P(z=1 | y)
  if (model$model_type %in% c("single", "community")) {
    y <- model$y
    n_obs <- nrow(y)
    max_visits <- ncol(y)
    z <- numeric(n_obs)
    for (i in seq_len(n_obs)) {
      yi <- y[i, ]
      valid <- yi >= 0
      if (any(yi[valid] == 1)) {
        z[i] <- 1  # Detected → occupied
      } else {
        # P(z=1|all zeros) = psi * prod(1-p) / [psi*prod(1-p) + (1-psi)]
        n_valid <- sum(valid)
        prod_1mp <- (1 - p[i])^n_valid
        z[i] <- psi[i] * prod_1mp / (psi[i] * prod_1mp + (1 - psi[i]))
      }
    }
  } else {
    z <- psi  # Approximate for dynamic models
  }

  list(psi = psi, p = p, z = z)
}

#' Residuals from occupancy model
#' @param object A `TulpaObs_fit` object.
#' @param type One of `"deviance"` (default), `"pearson"`, or `"response"`.
#' @param ... Ignored.
#' @return A list with `occ` (site-level) and `det` (visit-level) residuals.
#' @export
residuals.TulpaObs_fit <- function(object, type = c("deviance", "pearson", "response"), ...) {
  type <- match.arg(type)
  fit_vals <- fitted(object)
  model <- object$model

  # Occupancy residuals (site-level)
  z_obs <- if (model$model_type %in% c("single", "community")) {
    apply(model$y, 1, function(row) as.integer(any(row[row >= 0] == 1)))
  } else {
    rep(NA_real_, model$n_sites)
  }

  occ_resid <- switch(type,
    response = z_obs - fit_vals$z,
    pearson = (z_obs - fit_vals$z) / sqrt(fit_vals$z * (1 - fit_vals$z) + 1e-10),
    deviance = {
      sign(z_obs - fit_vals$z) * sqrt(2 * abs(
        ifelse(z_obs == 1,
               -log(fit_vals$z + 1e-10),
               -log(1 - fit_vals$z + 1e-10))
      ))
    }
  )

  # Detection residuals (visit-level, for single/community)
  det_resid <- NULL
  if (model$model_type %in% c("single", "community")) {
    y <- model$y
    p_hat <- fit_vals$p
    n_obs <- nrow(y)
    max_visits <- ncol(y)
    det_resid <- matrix(NA_real_, n_obs, max_visits)
    for (i in seq_len(n_obs)) {
      for (j in seq_len(max_visits)) {
        if (y[i, j] >= 0) {
          expected <- fit_vals$z[i] * p_hat[i]
          det_resid[i, j] <- switch(type,
            response = y[i, j] - expected,
            pearson = (y[i, j] - expected) / sqrt(expected * (1 - expected) + 1e-10),
            deviance = {
              sign(y[i, j] - expected) * sqrt(2 * abs(
                ifelse(y[i, j] == 1, -log(expected + 1e-10), -log(1 - expected + 1e-10))
              ))
            }
          )
        }
      }
    }
  }

  list(occ = occ_resid, det = det_resid)
}

#' Simulate replicate datasets from posterior
#' @param object A `TulpaObs_fit` object.
#' @param nsim Number of simulated datasets (default 1).
#' @param seed Optional random seed.
#' @param ... Ignored.
#' @return A list of simulated detection history matrices.
#' @export
simulate.TulpaObs_fit <- function(object, nsim = 1, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  model <- object$model
  draws <- object$draws
  n_samples <- nrow(draws)
  pi_list <- model$process_info

  if (model$model_type != "single") {
    stop("simulate() currently only supports single-season models")
  }

  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  n_sites <- model$n_sites
  max_visits <- model$max_visits

  result <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    # Sample a posterior draw
    draw_idx <- sample.int(n_samples, 1)
    beta_occ <- draws[draw_idx, seq_len(pi_list[[1]]$p)]
    beta_det <- draws[draw_idx, pi_list[[1]]$p + seq_len(pi_list[[2]]$p)]

    psi <- plogis(as.vector(X_occ %*% beta_occ))
    p <- plogis(as.vector(X_det %*% beta_det))
    z <- rbinom(n_sites, 1, psi)

    y_sim <- matrix(NA_integer_, n_sites, max_visits)
    for (i in seq_len(n_sites)) {
      # Respect original visit structure (NA pattern)
      for (j in seq_len(max_visits)) {
        if (model$y[i, j] >= 0) {
          y_sim[i, j] <- rbinom(1, 1, z[i] * p[i])
        }
      }
    }
    result[[s]] <- y_sim
  }
  if (nsim == 1) result[[1]] else result
}

#' Predict from occupancy model
#'
#' Three modes:
#' - **In-sample**: `predict(fit)` returns fitted values.
#' - **Design-matrix**: `predict(fit, X.0 = ...)` predicts at new covariate values.
#' - **Terms-based**: `predict(fit, terms = "elev")` varies one covariate, others at mean.
#'
#' @param object A `TulpaObs_fit` object.
#' @param X.0 Optional design matrix for occupancy prediction.
#' @param type `"occupancy"` (default), `"detection"`, or `"both"`.
#' @param quantiles Quantile levels for credible intervals.
#' @param terms Character vector of terms to vary (ggpredict-style).
#' @param n_points Number of prediction points per continuous term.
#' @param ... Ignored.
#' @return Depends on mode. In-sample: `fitted()` result. Design-matrix/terms:
#'   data.frame with estimate and CIs.
#' @export
predict.TulpaObs_fit <- function(object, X.0 = NULL,
                                 type = c("occupancy", "detection", "both"),
                                 quantiles = c(0.025, 0.5, 0.975),
                                 terms = NULL, n_points = 50L, ...) {
  type <- match.arg(type)

  # In-sample mode
  if (is.null(X.0) && is.null(terms)) return(fitted(object))

  model <- object$model
  draws <- object$draws
  pi_list <- model$process_info

  # Terms-based mode
  if (!is.null(terms)) {
    return(predict_terms(object, terms, type, quantiles, n_points))
  }

  # Design-matrix mode
  n_pred <- nrow(X.0)
  n_draws <- nrow(draws)
  p_occ <- pi_list[[1]]$p

  if (ncol(X.0) != p_occ) {
    stop(sprintf("X.0 has %d columns but model has %d occupancy coefficients",
                 ncol(X.0), p_occ))
  }

  # Compute predictions for each draw
  psi_draws <- matrix(NA_real_, n_draws, n_pred)
  for (s in seq_len(n_draws)) {
    beta <- draws[s, seq_len(p_occ)]
    psi_draws[s, ] <- plogis(as.vector(X.0 %*% beta))
  }

  data.frame(
    mean = colMeans(psi_draws),
    sd = apply(psi_draws, 2, sd),
    q2.5 = apply(psi_draws, 2, quantile, quantiles[1]),
    q50 = apply(psi_draws, 2, quantile, quantiles[2]),
    q97.5 = apply(psi_draws, 2, quantile, quantiles[3])
  )
}

# Terms-based prediction (ggpredict-style)
predict_terms <- function(object, terms, type, quantiles, n_points) {
  model <- object$model
  draws <- object$draws
  pi_list <- model$process_info

  proc_idx <- if (type == "detection") 2 else 1
  p_proc <- pi_list[[proc_idx]]$p
  X_orig <- model$X_processes[[proc_idx]]
  coef_names <- pi_list[[proc_idx]]$coef_names
  beta_offset <- if (proc_idx > 1) sum(vapply(pi_list[1:(proc_idx-1)], function(pi) pi$p, integer(1))) else 0

  # Parse the first term (simple: just a variable name for now)
  term_var <- terms[1]
  col_idx <- match(term_var, coef_names)
  if (is.na(col_idx)) {
    stop(sprintf("term '%s' not found in %s coefficients: %s",
                 term_var, pi_list[[proc_idx]]$name,
                 paste(coef_names, collapse = ", ")))
  }

  # Create prediction grid: vary term_var, hold others at mean
  x_range <- range(X_orig[, col_idx])
  x_grid <- seq(x_range[1], x_range[2], length.out = n_points)

  X_pred <- matrix(colMeans(X_orig), nrow = n_points, ncol = p_proc, byrow = TRUE)
  X_pred[, col_idx] <- x_grid

  # Predict from each draw
  n_draws <- nrow(draws)
  pred_draws <- matrix(NA_real_, n_draws, n_points)
  for (s in seq_len(n_draws)) {
    beta <- draws[s, beta_offset + seq_len(p_proc)]
    pred_draws[s, ] <- plogis(as.vector(X_pred %*% beta))
  }

  result <- data.frame(
    x = x_grid,
    estimate = colMeans(pred_draws),
    lower = apply(pred_draws, 2, quantile, quantiles[1]),
    upper = apply(pred_draws, 2, quantile, quantiles[3])
  )
  attr(result, "term") <- term_var
  attr(result, "process") <- pi_list[[proc_idx]]$name
  class(result) <- c("occu_prediction", "data.frame")
  result
}

#' @export
plot.occu_prediction <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    plot(x$x, x$estimate, type = "l", ylim = range(c(x$lower, x$upper)),
         xlab = attr(x, "term"), ylab = attr(x, "process"),
         main = sprintf("Effect of %s on %s", attr(x, "term"), attr(x, "process")))
    polygon(c(x$x, rev(x$x)), c(x$lower, rev(x$upper)),
            col = rgb(0, 0, 0, 0.1), border = NA)
    return(invisible(x))
  }
  p <- ggplot2::ggplot(x, ggplot2::aes(x = .data$x, y = .data$estimate)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), alpha = 0.2) +
    ggplot2::geom_line() +
    ggplot2::labs(x = attr(x, "term"), y = attr(x, "process"))
  print(p)
  invisible(x)
}

#' Compute marginal effect of a covariate
#' @param object A `TulpaObs_fit` object.
#' @param covariate Name of covariate.
#' @param process `"occupancy"` (default) or `"detection"`.
#' @param n_points Number of prediction points (default 100).
#' @return A data.frame with covariate values and predicted probabilities.
#' @export
marginal_effect <- function(object, covariate,
                            process = c("occupancy", "detection"),
                            n_points = 100L) {
  process <- match.arg(process)
  type <- if (process == "detection") "detection" else "occupancy"
  predict_terms(object, terms = covariate, type = type,
                quantiles = c(0.025, 0.5, 0.975), n_points = n_points)
}

#' Estimate species richness from community model
#' @param object A `TulpaObs_fit` object from a community model.
#' @return A data.frame with site-level richness estimates.
#' @export
richness <- function(object) {
  model <- object$model
  if (model$model_type != "community") {
    stop("richness() requires a community model")
  }

  draws <- object$draws
  pi_list <- model$process_info
  X_occ <- model$X_processes[[1]]
  p_occ <- pi_list[[1]]$p
  n_draws <- nrow(draws)
  n_sites <- model$n_sites
  n_species <- model$n_species
  N <- model$N

  # For each draw, compute psi per site-species, sum across species per site
  richness_draws <- matrix(0, n_draws, n_sites)
  for (s in seq_len(n_draws)) {
    beta_occ <- draws[s, seq_len(p_occ)]
    psi <- plogis(as.vector(X_occ %*% beta_occ))
    # psi is length N = n_sites * n_species
    # obs = (site-1)*n_species + species
    for (i in seq_len(n_sites)) {
      idx <- (i - 1) * n_species + seq_len(n_species)
      richness_draws[s, i] <- sum(psi[idx])
    }
  }

  data.frame(
    site = seq_len(n_sites),
    mean = colMeans(richness_draws),
    sd = apply(richness_draws, 2, sd),
    q2.5 = apply(richness_draws, 2, quantile, 0.025),
    q97.5 = apply(richness_draws, 2, quantile, 0.975)
  )
}

# tidy, glance, ranef inherited from tulpa::tulpa_fit

#' Update and refit an occupancy model
#' @param object A `TulpaObs_fit` object.
#' @param ... Named arguments to override in `occu_fit()`.
#' @param evaluate If TRUE (default), refit the model.
#' @return Updated `TulpaObs_fit` object (or call if `evaluate = FALSE`).
#' @export
update.TulpaObs_fit <- function(object, ..., evaluate = TRUE) {
  args <- list(model = object$model,
               spatial = object$spatial,
               temporal = object$temporal,
               re = object$re,
               svc = object$svc,
               latent = object$latent)
  dots <- list(...)
  for (nm in names(dots)) args[[nm]] <- dots[[nm]]

  if (!evaluate) return(args)
  do.call(occu_fit, args)
}

#' Check model identifiability
#'
#' Diagnostics for potential identifiability issues in occupancy models.
#' Checks for: confounded covariates, low detection rates, sparse data.
#'
#' @param model A `TulpaObs` model object (before fitting).
#' @param fit Optional `TulpaObs_fit` object (for post-fit diagnostics).
#' @return A list with diagnostic messages and flags.
#' @export
checkIdentifiability <- function(model, fit = NULL) {
  issues <- character()

  if (!inherits(model, "TulpaObs")) {
    stop("model must be a TulpaObs object")
  }

  # Pre-fit checks
  if (model$model_type %in% c("single", "community")) {
    y <- model$y
    n_sites <- nrow(y)
    max_visits <- ncol(y)

    # Naive occupancy
    n_detected <- sum(apply(y, 1, function(row) any(row[row >= 0] == 1)))
    naive_occ <- n_detected / n_sites
    if (naive_occ < 0.05) {
      issues <- c(issues, sprintf("Very low naive occupancy (%.1f%%). Model may struggle to estimate occupancy coefficients.", 100 * naive_occ))
    }
    if (naive_occ > 0.95) {
      issues <- c(issues, sprintf("Very high naive occupancy (%.1f%%). Little information to estimate occupancy effects.", 100 * naive_occ))
    }

    # Mean visits per site
    n_visits <- apply(y, 1, function(row) sum(row >= 0))
    if (mean(n_visits) < 2) {
      issues <- c(issues, sprintf("Very few visits per site (mean %.1f). Detection and occupancy may be confounded.", mean(n_visits)))
    }

    # Check for collinearity in occupancy covariates
    X_occ <- model$X_processes[[1]]
    if (ncol(X_occ) > 2) {
      cors <- cor(X_occ[, -1, drop = FALSE])
      high_cor <- which(abs(cors) > 0.8 & upper.tri(cors), arr.ind = TRUE)
      if (nrow(high_cor) > 0) {
        names_occ <- colnames(X_occ)[-1]
        for (k in seq_len(nrow(high_cor))) {
          issues <- c(issues, sprintf("High correlation (%.2f) between %s and %s.",
                                      cors[high_cor[k, 1], high_cor[k, 2]],
                                      names_occ[high_cor[k, 1]],
                                      names_occ[high_cor[k, 2]]))
        }
      }
    }
  }

  # Post-fit checks
  if (!is.null(fit) && inherits(fit, "TulpaObs_fit")) {
    if (sum(fit$divergent) > 0) {
      issues <- c(issues, sprintf("%d divergent transitions. Consider increasing adapt_delta or reparameterizing.", sum(fit$divergent)))
    }
    if (mean(fit$accept_prob) < 0.5) {
      issues <- c(issues, sprintf("Low mean acceptance probability (%.2f). Model may be poorly specified.", mean(fit$accept_prob)))
    }
  }

  result <- list(
    identifiable = length(issues) == 0,
    issues = issues
  )
  if (length(issues) > 0) {
    for (msg in issues) message("- ", msg)
  } else {
    message("No identifiability issues detected.")
  }
  invisible(result)
}

#' Specify prior distributions for occupancy models
#'
#' @param beta.normal Prior for occupancy fixed effects: `list(mean, sd)`.
#' @param alpha.normal Prior for detection fixed effects: `list(mean, sd)`.
#' @param sigma.sq.psi Prior for occupancy RE variance: `c(shape, rate)`.
#' @param sigma.sq.p Prior for detection RE variance: `c(shape, rate)`.
#' @return A `TulpaObs_priors` object.
#' @export
occu_priors <- function(beta.normal = list(mean = 0, sd = sqrt(2.72)),
                        alpha.normal = list(mean = 0, sd = sqrt(2.72)),
                        sigma.sq.psi = c(0.1, 0.1),
                        sigma.sq.p = c(0.1, 0.1)) {
  structure(list(
    beta_mean = beta.normal$mean,
    beta_sd = beta.normal$sd,
    alpha_mean = alpha.normal$mean,
    alpha_sd = alpha.normal$sd,
    sigma_sq_psi_shape = sigma.sq.psi[1],
    sigma_sq_psi_rate = sigma.sq.psi[2],
    sigma_sq_p_shape = sigma.sq.p[1],
    sigma_sq_p_rate = sigma.sq.p[2]
  ), class = "TulpaObs_priors")
}

#' @export
print.TulpaObs_priors <- function(x, ...) {
  cat("TulpaObs priors:\n")
  cat(sprintf("  beta ~ Normal(%.2f, %.2f)\n", x$beta_mean, x$beta_sd))
  cat(sprintf("  alpha ~ Normal(%.2f, %.2f)\n", x$alpha_mean, x$alpha_sd))
  invisible(x)
}


# ============================================================================
# spOccupancy $ compatibility accessor
# ============================================================================

#' Access spOccupancy-compatible fields from TulpaObs fits
#'
#' Allows accessing spOccupancy-style fields (e.g., `$beta.samples`,
#' `$psi.samples`) on TulpaObs_fit objects. Since TulpaObs stores actual
#' posterior draws, this is a thin remapping layer.
#'
#' @param x A `TulpaObs_fit` object.
#' @param name Field name to access.
#' @return The requested field value.
#' @export
`$.TulpaObs_fit` <- function(x, name) {
  # First check native fields
  val <- .subset2(x, name)
  if (!is.null(val)) return(val)

  model <- .subset2(x, "model")
  draws <- .subset2(x, "draws")
  pi_list <- if (!is.null(model)) model$process_info else NULL

  switch(name,
    # Occupancy fixed effect draws
    "beta.samples" = {
      if (is.null(pi_list)) return(NULL)
      p_occ <- pi_list[[1]]$p
      draws[, seq_len(p_occ), drop = FALSE]
    },

    # Detection fixed effect draws
    "alpha.samples" = {
      if (is.null(pi_list) || length(pi_list) < 2) return(NULL)
      p_occ <- pi_list[[1]]$p
      p_det <- pi_list[[2]]$p
      draws[, p_occ + seq_len(p_det), drop = FALSE]
    },

    # Occupancy probabilities (n_draws x n_sites)
    "psi.samples" = {
      fit_vals <- fitted(x)
      if (!is.null(fit_vals$psi)) {
        # Recompute from draws for proper uncertainty
        if (is.null(pi_list)) return(NULL)
        p_occ <- pi_list[[1]]$p
        X_occ <- model$X_processes[[1]]
        n_draws <- nrow(draws)
        psi_mat <- matrix(NA_real_, n_draws, nrow(X_occ))
        for (s in seq_len(n_draws)) {
          beta <- draws[s, seq_len(p_occ)]
          psi_mat[s, ] <- plogis(as.vector(X_occ %*% beta))
        }
        psi_mat
      }
    },

    # Latent occupancy state
    "z.samples" = {
      psi <- x$psi.samples
      if (!is.null(psi)) {
        matrix(rbinom(length(psi), 1, psi), nrow = nrow(psi))
      }
    },

    # Detection probabilities
    "p.samples" = {
      if (is.null(pi_list) || length(pi_list) < 2) return(NULL)
      p_occ <- pi_list[[1]]$p
      p_det <- pi_list[[2]]$p
      X_det <- model$X_processes[[2]]
      n_draws <- nrow(draws)
      p_mat <- matrix(NA_real_, n_draws, nrow(X_det))
      for (s in seq_len(n_draws)) {
        alpha <- draws[s, p_occ + seq_len(p_det)]
        p_mat[s, ] <- plogis(as.vector(X_det %*% alpha))
      }
      p_mat
    },

    # Computation time
    "run.time" = .subset2(x, "elapsed"),

    # Not a compat field
    NULL
  )
}


# ============================================================================
# Spatial prediction at new locations
# ============================================================================

#' Predict occupancy at new spatial locations
#'
#' Generates occupancy predictions at new coordinates, including the
#' spatial random effect interpolated from the fitted field.
#'
#' @param object A `TulpaObs_fit` object fitted with a spatial component.
#' @param newcoords Matrix of new coordinates (n_new x 2).
#' @param newocc.covs Optional data.frame of covariates at new locations.
#' @param quantiles Quantiles for credible intervals (default 0.025, 0.5, 0.975).
#' @return A data.frame with `mean`, `sd`, and quantile columns.
#' @export
predict_spatial <- function(object, newcoords, newocc.covs = NULL,
                            quantiles = c(0.025, 0.5, 0.975)) {
  if (is.null(object$spatial)) {
    stop("predict_spatial requires a model fitted with a spatial component", call. = FALSE)
  }

  # Build design matrix at new locations
  if (!is.null(newocc.covs)) {
    # Use model's formula to build X
    model <- object$model
    occ_formula <- model$occ_formula
    if (!is.null(occ_formula)) {
      X.0 <- model.matrix(occ_formula, data = newocc.covs)
    } else {
      X.0 <- as.matrix(cbind(1, newocc.covs))
    }
  } else {
    n_new <- nrow(newcoords)
    p_occ <- object$model$process_info[[1]]$p
    X.0 <- matrix(0, n_new, p_occ)
    X.0[, 1] <- 1  # Intercept only
  }

  # Fixed effect prediction
  draws <- object$draws
  p_occ <- object$model$process_info[[1]]$p
  n_draws <- nrow(draws)
  n_new <- nrow(newcoords)

  eta_draws <- matrix(NA_real_, n_draws, n_new)
  for (s in seq_len(n_draws)) {
    beta <- draws[s, seq_len(p_occ)]
    eta_draws[s, ] <- as.vector(X.0 %*% beta)
  }

  # Interpolate spatial field to new locations using nearest-neighbor
  sp_type <- object$spatial$type
  cn <- colnames(draws)
  sp_cols <- grep("^phi_spatial\\[|^w_gp\\[|^gp_w\\[", cn)

  if (length(sp_cols) > 0 && !is.null(object$spatial$coords)) {
    fit_coords <- object$spatial$coords
    n_fit <- nrow(fit_coords)

    # Compute distances from new points to fitted points
    for (s in seq_len(n_draws)) {
      sp_effects <- draws[s, sp_cols]
      # Nearest-neighbor interpolation
      for (i in seq_len(n_new)) {
        dists <- sqrt((fit_coords[, 1] - newcoords[i, 1])^2 +
                       (fit_coords[, 2] - newcoords[i, 2])^2)
        # IDW with k=5 nearest neighbors
        k <- min(5, n_fit)
        nn <- order(dists)[seq_len(k)]
        w <- 1 / (dists[nn] + 1e-10)
        w <- w / sum(w)
        eta_draws[s, i] <- eta_draws[s, i] + sum(w * sp_effects[nn])
      }
    }
  }

  psi_draws <- plogis(eta_draws)

  result <- data.frame(
    mean = colMeans(psi_draws),
    sd = apply(psi_draws, 2, sd)
  )
  for (q in quantiles) {
    qname <- paste0("q", gsub("\\.", "", format(q * 100, nsmall = 1)))
    result[[qname]] <- apply(psi_draws, 2, quantile, q)
  }
  result
}
