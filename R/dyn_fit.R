#' Fit a multi-season dynamic occupancy model
#'
#' @param model A `tulpaOcc_dynmodel` object from [dynOcc()]
#' @param sigma_beta Prior SD for regression coefficients (default 10)
#' @param iter Total NUTS iterations (default 2000)
#' @param warmup Warmup iterations (default 1000)
#' @param max_treedepth Maximum NUTS tree depth (default 10)
#' @param adapt_delta Target acceptance rate (default 0.8)
#' @param seed Random seed
#' @param verbose Print sampler progress (default TRUE)
#'
#' @return A `tulpaOcc_dynfit` object with posterior draws and diagnostics
#' @export
dynOcc_fit <- function(model, sigma_beta = 10, iter = 2000, warmup = 1000,
                       max_treedepth = 10, adapt_delta = 0.8, seed = 42,
                       verbose = TRUE) {

  if (!inherits(model, "tulpaOcc_dynmodel")) {
    stop("model must be a tulpaOcc_dynmodel object from dynOcc()")
  }

  fit <- cpp_dyn_occ_fit(
    y_flat_r = model$y_flat,
    n_visits_r = model$n_visits,
    any_detected_r = model$any_detected,
    X_occ_r = model$X_occ,
    X_det_r = model$X_det,
    X_col_r = model$X_col,
    X_ext_r = model$X_ext,
    n_sites = model$n_sites,
    n_seasons = model$n_seasons,
    max_visits = model$max_visits,
    sigma_beta = sigma_beta,
    n_iter = as.integer(iter),
    n_warmup = as.integer(warmup),
    max_treedepth = as.integer(max_treedepth),
    adapt_delta = adapt_delta,
    seed = as.integer(seed),
    verbose = verbose
  )

  # Add readable parameter names
  param_names <- c(
    paste0("psi1_", model$occ_names),
    paste0("p_", model$det_names),
    paste0("gamma_", model$col_names),
    paste0("epsilon_", model$ext_names)
  )
  fit$param_names <- param_names

  # Back-transform intercepts
  fit$mean_psi1 <- plogis(fit$means[1])
  fit$mean_p <- plogis(fit$means[model$p_occ + 1])
  fit$mean_gamma <- plogis(fit$means[model$p_occ + model$p_det + 1])
  fit$mean_epsilon <- plogis(fit$means[model$p_occ + model$p_det + model$p_col + 1])

  fit$model <- model
  class(fit) <- "tulpaOcc_dynfit"
  fit
}

#' @export
print.tulpaOcc_dynfit <- function(x, ...) {
  cat("tulpaOcc fit (dynamic occupancy, NUTS)\n")
  cat(sprintf("  Sites: %d, Seasons: %d, Max visits: %d\n",
              x$model$n_sites, x$model$n_seasons, x$model$max_visits))
  cat(sprintf("  Samples: %d, Step size: %.4f\n", x$n_samples, x$epsilon))
  n_div <- sum(x$divergent)
  if (n_div > 0) cat(sprintf("  WARNING: %d divergent transitions\n", n_div))
  cat("\n")

  cat(sprintf("Mean initial occupancy: %.3f\n", x$mean_psi1))
  cat(sprintf("Mean detection: %.3f\n", x$mean_p))
  cat(sprintf("Mean colonization: %.3f\n", x$mean_gamma))
  cat(sprintf("Mean extinction: %.3f\n", x$mean_epsilon))
  invisible(x)
}

#' @export
summary.tulpaOcc_dynfit <- function(object, ...) {
  draws <- object$draws
  n_params <- ncol(draws)
  n_named <- min(length(object$param_names), n_params)

  result <- data.frame(
    mean = colMeans(draws),
    sd = apply(draws, 2, sd),
    q2.5 = apply(draws, 2, quantile, 0.025),
    q50 = apply(draws, 2, quantile, 0.50),
    q97.5 = apply(draws, 2, quantile, 0.975)
  )
  if (n_named > 0) {
    rownames(result)[1:n_named] <- object$param_names[1:n_named]
  }
  result
}
