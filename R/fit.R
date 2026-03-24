#' Fit a single-season occupancy model
#'
#' Uses tulpa's full NUTS sampler with dual averaging and mass matrix adaptation.
#'
#' @param model A `tulpaOcc_model` object from [occ()]
#' @param spatial Optional spatial specification from [occ_icar()], [occ_bym2()],
#'   [occ_gp()], or [occ_multiscale_gp()]
#' @param sigma_beta Prior SD for regression coefficients (default 10)
#' @param iter Total NUTS iterations (default 2000)
#' @param warmup Warmup iterations (default 1000)
#' @param max_treedepth Maximum NUTS tree depth (default 10)
#' @param adapt_delta Target acceptance rate (default 0.8)
#' @param seed Random seed
#' @param verbose Print sampler progress (default TRUE)
#'
#' @return A `tulpaOcc_fit` object with posterior draws and diagnostics
#' @export
occ_fit <- function(model, spatial = NULL, sigma_beta = 10,
                    iter = 2000, warmup = 1000,
                    max_treedepth = 10, adapt_delta = 0.8, seed = 42,
                    verbose = TRUE) {

  if (!inherits(model, "tulpaOcc_model")) {
    stop("model must be a tulpaOcc_model object from occ()")
  }

  if (!is.null(spatial) && !inherits(spatial, "tulpaOcc_spatial")) {
    stop("spatial must be a tulpaOcc_spatial object (occ_icar, occ_bym2, occ_gp, occ_multiscale_gp)")
  }

  # Build spatial params list for C++
  spatial_params <- build_spatial_params(spatial, model$n_sites)

  fit <- cpp_occ_fit(
    y_r = model$y,
    X_occ_r = model$X_occ,
    X_det_r = model$X_det,
    X_det_visit_r = model$X_det_visit,
    spatial_params_r = spatial_params,
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
    paste0("psi_", model$occ_names),
    paste0("p_", model$det_names)
  )
  if (model$p_det_visit > 0) {
    param_names <- c(param_names,
                     paste0("p_visit_", model$det_visit_names))
  }
  fit$param_names <- param_names

  # Back-transform occupancy and detection intercepts to probability scale
  occ_intercept_idx <- 1
  det_intercept_idx <- model$p_occ + 1
  fit$mean_psi <- plogis(fit$means[occ_intercept_idx])
  fit$mean_p <- plogis(fit$means[det_intercept_idx])

  fit$model <- model
  fit$spatial <- spatial
  class(fit) <- "tulpaOcc_fit"
  fit
}

# Build spatial params list for C++ from spatial spec (or NULL)
build_spatial_params <- function(spatial, n_sites) {
  if (is.null(spatial)) {
    return(list(type = "none"))
  }

  params <- list(type = spatial$type)

  if (spatial$type %in% c("icar", "bym2")) {
    if (spatial$n_units != n_sites) {
      stop(sprintf("spatial has %d units but model has %d sites",
                   spatial$n_units, n_sites))
    }
    params$n_units <- spatial$n_units
    params$adj_row_ptr <- spatial$adj_row_ptr
    params$adj_col_idx <- spatial$adj_col_idx
    params$n_neighbors <- spatial$n_neighbors
    params$spatial_shared_occ <- spatial$shared[1]
    params$spatial_shared_det <- spatial$shared[2]
    if (spatial$type == "bym2") {
      params$scale_factor <- spatial$scale_factor
    }
  } else if (spatial$type == "gp") {
    if (spatial$n_obs != n_sites) {
      stop(sprintf("spatial has %d locations but model has %d sites",
                   spatial$n_obs, n_sites))
    }
    params$n_obs <- spatial$n_obs
    params$nn <- spatial$nn
    params$coords <- spatial$coords
    params$nn_idx <- spatial$nn_idx
    params$nn_dist <- spatial$nn_dist
    params$nn_neighbor_dist <- spatial$nn_neighbor_dist
    params$nn_order <- spatial$nn_order
    params$nn_order_inv <- spatial$nn_order_inv
    params$cov_type <- spatial$cov_type
    params$nu <- spatial$nu
    params$spatial_shared_occ <- spatial$shared[1]
    params$spatial_shared_det <- spatial$shared[2]
    params$sigma2_prior_U <- spatial$sigma2_prior_U
    params$sigma2_prior_alpha <- spatial$sigma2_prior_alpha
    params$phi_prior_lower <- spatial$phi_prior_lower
    params$phi_prior_upper <- spatial$phi_prior_upper
  } else if (spatial$type == "multiscale_gp") {
    if (spatial$n_obs != n_sites) {
      stop(sprintf("spatial has %d locations but model has %d sites",
                   spatial$n_obs, n_sites))
    }
    params$n_obs <- spatial$n_obs
    params$coords <- spatial$coords
    params$nn_local <- spatial$nn_local
    params$nn_idx_local <- spatial$nn_idx_local
    params$nn_dist_local <- spatial$nn_dist_local
    params$nn_neighbor_dist_local <- spatial$nn_neighbor_dist_local
    params$nn_order_local <- spatial$nn_order_local
    params$nn_order_inv_local <- spatial$nn_order_inv_local
    params$nn_regional <- spatial$nn_regional
    params$nn_idx_regional <- spatial$nn_idx_regional
    params$nn_dist_regional <- spatial$nn_dist_regional
    params$nn_neighbor_dist_regional <- spatial$nn_neighbor_dist_regional
    params$nn_order_regional <- spatial$nn_order_regional
    params$nn_order_inv_regional <- spatial$nn_order_inv_regional
    params$cov_type <- spatial$cov_type
    params$nu <- spatial$nu
    params$spatial_shared_occ <- spatial$shared[1]
    params$spatial_shared_det <- spatial$shared[2]
    params$range_local_lower <- spatial$range_local_lower
    params$range_local_upper <- spatial$range_local_upper
    params$range_regional_lower <- spatial$range_regional_lower
    params$range_regional_upper <- spatial$range_regional_upper
    params$sigma2_local_prior_U <- spatial$sigma2_local_prior_U
    params$sigma2_local_prior_alpha <- spatial$sigma2_local_prior_alpha
    params$sigma2_regional_prior_U <- spatial$sigma2_regional_prior_U
    params$sigma2_regional_prior_alpha <- spatial$sigma2_regional_prior_alpha
  }

  params
}

#' @export
print.tulpaOcc_fit <- function(x, ...) {
  cat("tulpaOcc fit (single-season occupancy, NUTS)\n")
  cat(sprintf("  Sites: %d, Max visits: %d\n",
              x$model$n_sites, x$model$max_visits))
  cat(sprintf("  Samples: %d, Step size: %.4f\n", x$n_samples, x$epsilon))
  n_div <- sum(x$divergent)
  if (n_div > 0) cat(sprintf("  WARNING: %d divergent transitions\n", n_div))
  cat("\n")

  cat("Occupancy coefficients (logit scale):\n")
  occ_idx <- seq_len(x$model$p_occ)
  occ_means <- x$means[occ_idx]
  names(occ_means) <- x$model$occ_names
  print(round(occ_means, 4))

  cat(sprintf("\nMean occupancy (intercept-only): %.3f\n", x$mean_psi))

  cat("\nDetection coefficients (logit scale):\n")
  det_idx <- x$model$p_occ + seq_len(x$model$p_det)
  det_means <- x$means[det_idx]
  names(det_means) <- x$model$det_names
  print(round(det_means, 4))

  cat(sprintf("\nMean detection (intercept-only): %.3f\n", x$mean_p))

  if (x$model$p_det_visit > 0) {
    cat("\nVisit-level detection coefficients:\n")
    visit_idx <- x$model$p_occ + x$model$p_det + seq_len(x$model$p_det_visit)
    visit_means <- x$means[visit_idx]
    names(visit_means) <- x$model$det_visit_names
    print(round(visit_means, 4))
  }

  invisible(x)
}

#' @export
summary.tulpaOcc_fit <- function(object, ...) {
  draws <- object$draws
  n_params <- ncol(draws)

  result <- data.frame(
    mean = colMeans(draws),
    sd = apply(draws, 2, sd),
    q2.5 = apply(draws, 2, quantile, 0.025),
    q50 = apply(draws, 2, quantile, 0.50),
    q97.5 = apply(draws, 2, quantile, 0.975)
  )
  rownames(result) <- object$param_names
  result
}
