#' Internal engine entry point
#'
#' Dispatches to Laplace (default) or NUTS for a built `tobs_model`. Not
#' user-facing; called from `tobs()` via the per-family `.dispatch_*` helpers.
#'
#' @keywords internal
.tobs_fit_model <- function(model, spatial = NULL, temporal = NULL,
                            re = NULL, svc = NULL, latent = NULL,
                            method = c("laplace", "nuts"),
                            priors = NULL,
                            sigma_beta = 10, sigma_re_scale = 1,
                            max_iter = 100L, tol = 1e-4, damping = 0.7,
                            iter = 2000, warmup = 1000,
                            max_treedepth = 10, adapt_delta = 0.8, seed = 42,
                            verbose = TRUE) {

  method <- match.arg(method)

  if (!inherits(model, "tobs_model")) {
    stop("model must be a tobs_model object (from `.tobs_build_model()`)")
  }

  if (method == "laplace") {
    return(.tobs_laplace(model, spatial = spatial, sigma_beta = sigma_beta,
                         max_iter = max_iter, tol = tol, damping = damping,
                         verbose = verbose))
  }

  if (!is.null(spatial) && !inherits(spatial, "tobs_spatial")) {
    stop("spatial must be a tobs_spatial object")
  }
  if (!is.null(temporal) && !inherits(temporal, "tobs_temporal")) {
    stop("temporal must be a tobs_temporal object from tobs_temporal()")
  }
  if (!is.null(svc) && !inherits(svc, "tobs_svc")) {
    stop("svc must be a tobs_svc object from tobs_svc()")
  }
  if (!is.null(latent) && !inherits(latent, "tobs_latent")) {
    stop("latent must be a tobs_latent object from tobs_latent()")
  }

  model_type <- model$model_type

  # ---- Build spec list for C++ ----
  spec <- list(
    model_type = model_type,
    sigma_beta = sigma_beta,
    n_iter = as.integer(iter),
    n_warmup = as.integer(warmup),
    max_treedepth = as.integer(max_treedepth),
    adapt_delta = adapt_delta,
    seed = as.integer(seed),
    verbose = verbose
  )

  # ---- Process design matrices ----
  spec$X_processes <- model$X_processes

  # ---- Process names for column labels ----
  spec$process_names <- lapply(model$process_info, function(pi) {
    paste0(pi$name, "_", pi$coef_names)
  })

  # ---- Model-type-specific fields ----
  if (model_type == "single") {
    spec$y <- model$y
    if (!is.null(model$X_det_visit)) {
      spec$X_det_visit <- model$X_det_visit
      spec$extra_param_names <- paste0("p_visit_", model$det_visit_names)
    }

  } else if (model_type == "dynamic") {
    spec$y_flat <- model$y_flat
    spec$n_visits <- model$n_visits
    spec$any_detected <- model$any_detected
    spec$n_sites <- model$n_sites
    spec$n_seasons <- model$n_seasons
    spec$max_visits <- model$max_visits

  } else if (model_type == "community") {
    spec$y <- model$y
    spec$re_group <- model$species_group
    spec$n_re_groups <- model$n_species
    spec$sigma_re_scale <- sigma_re_scale
    spec$re_shared_occ <- TRUE
    spec$re_shared_det <- TRUE
    spec$n_sites_raw <- model$n_sites

  } else if (model_type == "integrated") {
    spec$y_sources <- model$y_sources
    spec$site_maps <- model$site_maps
    spec$n_sources <- model$n_sources
    spec$n_sites <- model$n_sites

  } else if (model_type == "jsdm") {
    spec$y_jsdm <- model$y_jsdm
    # Species RE for JSDM (like community)
    spec$re_group <- model$species_group
    spec$n_re_groups <- model$n_species
    spec$sigma_re_scale <- sigma_re_scale
    spec$re_shared_occ <- TRUE
  }

  # ---- Spatial ----
  if (!is.null(spatial)) {
    spatial_params <- build_spatial_params(spatial, model$n_sites)
    spec$spatial_params <- spatial_params

    # For community models, build spatial_group mapping site-species -> site
    if (model_type == "community") {
      n_species <- model$n_species
      n_sites <- model$n_sites
      # spatial_group[obs] = site index (1-based)
      # obs = (site-1)*n_species + species, so site = floor((obs-1)/n_species) + 1
      spec$spatial_group <- as.integer(rep(seq_len(n_sites), each = n_species))
    }
  }

  # ---- Temporal ----
  if (!is.null(temporal)) {
    temp_spec <- list(type = temporal$type, shared = temporal$shared,
                      cyclic = temporal$cyclic,
                      tau_shape = temporal$tau_shape,
                      tau_rate = temporal$tau_rate)
    # Resolve time index from variable name or direct vector
    if (is.character(temporal$time)) {
      temp_spec$time_idx <- as.integer(as.factor(model$data[[temporal$time]]))
      temp_spec$n_times <- length(unique(temp_spec$time_idx))
    } else {
      temp_spec$time_idx <- as.integer(temporal$time)
      temp_spec$n_times <- max(temp_spec$time_idx)
    }
    # Optional group index
    if (!is.null(temporal$group)) {
      if (is.character(temporal$group)) {
        temp_spec$group_idx <- as.integer(as.factor(model$data[[temporal$group]]))
        temp_spec$n_groups <- length(unique(temp_spec$group_idx))
      } else {
        temp_spec$group_idx <- as.integer(temporal$group)
        temp_spec$n_groups <- max(temp_spec$group_idx)
      }
    }
    spec$temporal_spec <- temp_spec
  }

  # ---- Random effects ----
  if (!is.null(re)) {
    # Accept single tobs_re or list of tobs_re
    if (inherits(re, "tobs_re")) re <- list(re)
    re_spec <- build_re_spec(re, model)
    spec$re_spec <- re_spec
  }

  # ---- SVC ----
  if (!is.null(svc)) {
    # Build X_svc from the design matrix columns
    X_occ <- model$X_processes[[1]]
    svc_indices_0based <- svc$indices - 1L  # C++ 0-based
    X_svc_flat <- numeric(nrow(X_occ) * svc$n_svc)
    for (j in seq_along(svc$indices)) {
      col <- svc$indices[j]
      for (i in seq_len(nrow(X_occ))) {
        X_svc_flat[(i - 1) * svc$n_svc + j] <- X_occ[i, col]
      }
    }
    spec$svc_spec <- list(
      n_obs = svc$n_obs, n_svc = svc$n_svc, nn = svc$nn,
      coords = svc$coords, svc_indices = svc_indices_0based,
      X_svc = X_svc_flat,
      nn_idx = svc$nn_idx, nn_dist = svc$nn_dist,
      nn_order = svc$nn_order, nn_order_inv = svc$nn_order_inv,
      cov_type = svc$cov_type, shared = svc$shared,
      sigma2_prior_scale = svc$sigma2_prior_scale,
      phi_prior_lower = svc$phi_prior_lower,
      phi_prior_upper = svc$phi_prior_upper
    )
  }

  # ---- Latent factors ----
  if (!is.null(latent)) {
    spec$latent_spec <- list(
      n_factors = latent$n_factors,
      shared = latent$shared,
      constraint = latent$constraint,
      sigma_prior_rate = latent$sigma_prior_rate
    )
  }

  # ---- Call unified C++ entry point ----
  fit <- cpp_occu_fit(spec)

  # ---- Build R parameter names ----
  param_names <- unlist(spec$process_names)
  if (!is.null(model$det_visit_names) && length(model$det_visit_names) > 0) {
    param_names <- c(param_names, paste0("p_visit_", model$det_visit_names))
  }
  fit$param_names <- param_names

  # ---- Compute probability-scale intercepts ----
  fit$intercepts <- compute_intercepts(model, fit$means)

  fit$model <- model
  fit$spatial <- spatial
  fit$temporal <- temporal
  fit$re <- re
  fit$svc <- svc
  fit$latent <- latent
  # Expose process_info at top level for tulpa generic S3 methods
  fit$process_info <- model$process_info
  class(fit) <- c("tobs_fit", "tulpa_fit")
  fit
}

# Build multi-term RE spec for C++ from list of occu_re objects
build_re_spec <- function(re_list, model) {
  n_terms <- length(re_list)
  N <- if (!is.null(model$N)) model$N else model$n_sites

  groups <- vector("list", n_terms)
  n_groups <- integer(n_terms)
  has_slopes <- FALSE
  n_coefs <- integer(n_terms)
  slope_matrices <- vector("list", n_terms)

  # Aggregate shared across all terms (union)
  max_proc <- length(model$process_info)
  shared <- rep(FALSE, max_proc)

  for (t in seq_len(n_terms)) {
    re <- re_list[[t]]

    # Resolve group assignment
    if (is.character(re$group)) {
      grp <- as.integer(as.factor(model$data[[re$group]]))
    } else {
      grp <- as.integer(re$group)
    }
    if (length(grp) != N) {
      stop(sprintf("RE group vector has %d elements but model has %d observations",
                   length(grp), N))
    }
    groups[[t]] <- grp
    n_groups[t] <- max(grp)

    # Sharing
    for (k in seq_along(re$shared)) {
      if (re$shared[k]) shared[k] <- TRUE
    }

    # Slopes
    if (re$type == "slope" && !is.null(re$covariate)) {
      has_slopes <- TRUE
      n_coefs[t] <- 2L  # Intercept + 1 slope
      # Build slope design matrix from model data
      X_slope <- model.matrix(as.formula(paste("~", re$covariate)), model$data)
      # Only the slope column (not intercept)
      slope_matrices[[t]] <- X_slope[, -1, drop = FALSE]
    } else {
      n_coefs[t] <- 1L  # Intercept only
    }
  }

  spec <- list(
    n_terms = n_terms,
    groups = groups,
    n_groups = n_groups,
    shared = shared,
    sigma_re_scale = re_list[[1]]$sigma_scale
  )

  if (has_slopes) {
    spec$has_slopes <- TRUE
    spec$n_coefs <- n_coefs
    spec$slope_matrices <- slope_matrices
    spec$correlated <- any(vapply(re_list, function(r) r$correlated, logical(1)))
  }

  spec
}

# Compute back-transformed intercepts on probability scale
compute_intercepts <- function(model, means) {
  result <- list()
  offset <- 0
  for (pi in model$process_info) {
    result[[pi$name]] <- plogis(means[offset + 1])
    offset <- offset + pi$p
  }
  result
}

# Build spatial params list for C++ from spatial spec (or NULL)
build_spatial_params <- function(spatial, n_sites) {
  if (is.null(spatial)) return(list(type = "none"))

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
    if (spatial$type == "bym2") params$scale_factor <- spatial$scale_factor

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
