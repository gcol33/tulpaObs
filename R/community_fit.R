#' Fit a community (multi-species) occupancy model
#'
#' Species-level variation is modeled as random effects (community-level
#' hierarchical prior). Each species gets a random intercept on both
#' occupancy and detection.
#'
#' @param model A `tulpaOcc_community` object from [communityOcc()]
#' @param sigma_beta Prior SD for community-mean regression coefficients (default 10)
#' @param sigma_re_scale Prior scale for species RE standard deviations (default 1)
#' @param iter Total NUTS iterations (default 2000)
#' @param warmup Warmup iterations (default 1000)
#' @param max_treedepth Maximum NUTS tree depth (default 10)
#' @param adapt_delta Target acceptance rate (default 0.8)
#' @param seed Random seed
#' @param verbose Print sampler progress (default TRUE)
#'
#' @return A `tulpaOcc_communityfit` object
#' @export
communityOcc_fit <- function(model, sigma_beta = 10, sigma_re_scale = 1,
                             iter = 2000, warmup = 1000,
                             max_treedepth = 10, adapt_delta = 0.8,
                             seed = 42, verbose = TRUE) {

  if (!inherits(model, "tulpaOcc_community")) {
    stop("model must be a tulpaOcc_community object from communityOcc()")
  }

  fit <- cpp_community_occ_fit(
    y_r = model$y,
    X_occ_r = model$X_occ,
    X_det_r = model$X_det,
    species_group_r = model$species_group,
    n_species = model$n_species,
    sigma_beta = sigma_beta,
    sigma_re_scale = sigma_re_scale,
    n_iter = as.integer(iter),
    n_warmup = as.integer(warmup),
    max_treedepth = as.integer(max_treedepth),
    adapt_delta = adapt_delta,
    seed = as.integer(seed),
    verbose = verbose
  )

  # Add readable parameter names
  param_names <- c(
    paste0("psi_", model$occ_names),
    paste0("p_", model$det_names)
  )
  fit$param_names <- param_names

  fit$mean_psi <- plogis(fit$means[1])
  fit$mean_p <- plogis(fit$means[model$p_occ + 1])

  fit$model <- model
  class(fit) <- "tulpaOcc_communityfit"
  fit
}

#' @export
print.tulpaOcc_communityfit <- function(x, ...) {
  cat("tulpaOcc fit (community occupancy, NUTS)\n")
  cat(sprintf("  Sites: %d, Species: %d\n",
              x$model$n_sites, x$model$n_species))
  cat(sprintf("  Samples: %d, Step size: %.4f\n", x$n_samples, x$epsilon))
  n_div <- sum(x$divergent)
  if (n_div > 0) cat(sprintf("  WARNING: %d divergent transitions\n", n_div))
  cat("\n")

  cat(sprintf("Community-mean occupancy (intercept): %.3f\n", x$mean_psi))
  cat(sprintf("Community-mean detection (intercept): %.3f\n", x$mean_p))
  invisible(x)
}

#' @export
summary.tulpaOcc_communityfit <- function(object, ...) {
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
