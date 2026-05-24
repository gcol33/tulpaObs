#' Internal engine entry point
#'
#' Dispatches to Laplace (default) or NUTS for a built `tobs_model`. Not
#' user-facing; called from `tobs()` via the per-family `.dispatch_*` helpers.
#' Spatial / temporal / random-effect / SVC / latent structure is read from
#' the structured terms the formula carried (`model$structured_terms`), not
#' from arguments — there is a single user-facing specification path.
#'
#' @keywords internal
.tobs_fit_model <- function(model,
                            method = c("laplace", "nested_laplace", "nuts"),
                            priors = NULL,
                            sigma_beta = 10, sigma_re_scale = 1,
                            max_iter = 100L, tol = 1e-4, damping = 0.7,
                            iter = 2000, warmup = 1000,
                            max_treedepth = 10, adapt_delta = 0.8, seed = 42,
                            approx = c("gaussian_laplace", "simplified_laplace"),
                            verbose = TRUE) {

  method <- match.arg(method)
  approx <- match.arg(approx)

  if (!inherits(model, "tobs_model")) {
    stop("model must be a tobs_model object (from `.tobs_build_model()`)")
  }

  # Engine-shaped structure specs derived from the formula's structured terms.
  structs  <- .tobs_structures_from_model(model)
  spatial  <- structs$spatial
  temporal <- structs$temporal
  re       <- structs$re
  svc      <- structs$svc
  latent   <- structs$latent

  # Autoscale every per-process design matrix before the engine sees it
  # (gcol33/tulpaObs#9). The engine optimizes on the centered+scaled
  # design; per-process betas / SEs / draws are transformed back to the
  # user-facing natural scale below. `model` (natural-scale) is restored
  # on the returned fit so `fitted()`, `residuals()`, `predict()`, and
  # diagnostics see the same X they would have without this hook.
  scale_info   <- .autoscale_model_X(model)
  fit_model    <- scale_info$model
  scales       <- scale_info$scales
  process_info <- model$process_info

  if (method == "laplace") {
    fit <- .tobs_laplace(fit_model, spatial = spatial, re = re,
                         priors = priors,
                         sigma_beta = sigma_beta,
                         max_iter = max_iter, tol = tol, damping = damping,
                         approx = approx,
                         verbose = verbose)
    fit <- .unscale_fit_per_process(fit, scales, process_info)
    fit$model      <- model
    fit$intercepts <- compute_intercepts(model, fit$means)
    return(fit)
  }

  if (method == "nested_laplace") {
    # Nested-Laplace path: single-season occupancy only. The driver
    # builds a multi-block latent prior from spatial + temporal + re
    # and routes the occupancy M-step block through
    # tulpa::tulpa_nested_laplace() via the per-block dispatcher in
    # tulpa::tulpa_em_laplace().
    nl_max_iter <- min(as.integer(max_iter), 25L)
    fit <- .tobs_em_nested_laplace(
      model    = fit_model,
      spatial  = spatial,
      temporal = temporal,
      re       = re,
      priors   = priors,
      sigma_beta = sigma_beta,
      max_iter = nl_max_iter,
      tol      = tol,
      damping  = damping,
      verbose  = verbose
    )
    fit <- .unscale_fit_per_process(fit, scales, process_info)
    fit$model      <- model
    fit$intercepts <- compute_intercepts(model, fit$means)
    return(fit)
  }

  # spatial / temporal / re / svc / latent are produced by
  # .tobs_structures_from_model(), which guarantees their classes; no
  # user-input validation is needed here.
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

  # ---- Process design matrices (autoscaled; gcol33/tulpaObs#9) ----
  spec$X_processes <- fit_model$X_processes

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
  # Index codes were resolved when the temporal() term was constructed.
  if (!is.null(temporal)) {
    temp_spec <- list(type = temporal$type, shared = temporal$shared,
                      cyclic = temporal$cyclic,
                      tau_shape = temporal$tau_shape,
                      tau_rate = temporal$tau_rate,
                      time_idx = temporal$time_idx,
                      n_times  = temporal$n_times)
    if (!is.null(temporal$group_idx)) {
      temp_spec$group_idx <- temporal$group_idx
      temp_spec$n_groups  <- temporal$n_groups
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
    # Build X_svc from the design matrix columns. Pull from the
    # (autoscaled) `fit_model` so the SVC base column values match the
    # global beta's parameterization seen by the optimizer
    # (gcol33/tulpaObs#9). For svc on the intercept the values are 1.0 in
    # both natural and scaled spaces; for svc on a non-intercept numeric
    # column the per-location offsets land on the scaled-column scale.
    X_occ <- fit_model$X_processes[[1]]
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

  # Unscale per-process beta slices in means / sds / draws (the engine
  # optimized on the centered+scaled design; gcol33/tulpaObs#9).
  fit <- .unscale_fit_per_process(fit, scales, process_info)

  # ---- Build R parameter names ----
  param_names <- unlist(spec$process_names)
  if (!is.null(model$det_visit_names) && length(model$det_visit_names) > 0) {
    param_names <- c(param_names, paste0("p_visit_", model$det_visit_names))
  }

  # Name the random-effect block (log_sigma / chol / z, type-blocked per
  # tulpa's layout) and reconstruct per-group BLUPs into `re_effects` so
  # summary() / ranef() label them instead of showing param[i]
  # (gcol33/tulpaObs#11). Counts and positions are unchanged.
  if (!is.null(re)) {
    re_design <- .tobs_re_design(if (inherits(re, "tobs_re")) list(re) else re,
                                 model)
    n_lead <- length(param_names)
    re_nms <- .tobs_re_nuts_param_names(re_design)
    if (n_lead + length(re_nms) <= length(fit$means)) {
      param_names <- c(param_names, re_nms)
      names(fit$means)[seq_along(param_names)] <- param_names
      if (!is.null(fit$draws) && ncol(fit$draws) >= length(param_names)) {
        colnames(fit$draws)[seq_along(param_names)] <- param_names
      }
      fit$re_effects <- tryCatch(
        .tobs_re_nuts_effects(fit$draws, re_design, n_lead),
        error = function(e) NULL)
    }
  }
  fit$param_names <- param_names

  # ---- Compute probability-scale intercepts (on natural-scale means) ----
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

# Translate the structured terms a formula carried (`model$structured_terms`)
# into the engine-shaped specs the fitter consumes. A term's process
# membership (`$processes`, set by the formula parser) becomes the length-2
# `shared = c(occ, det)` vector the C++ engine and Laplace paths read.
#
# The engine carries two sharing booleans (occupancy/state + detection), so a
# structured effect on a later process (colonization / extinction / extra
# integrated sources) is rejected rather than silently dropped. At most one
# spatial / temporal / svc / latent term is allowed; random effects may be
# multiple.
.tobs_structures_from_model <- function(model) {
  out <- list(spatial = NULL, temporal = NULL, re = NULL,
              svc = NULL, latent = NULL)
  terms <- model$structured_terms
  if (is.null(terms) || length(terms) == 0L) return(out)

  proc_names <- vapply(model$process_info, function(pi) pi$name, character(1))
  re_list <- list()

  one <- function(slot, label) {
    if (!is.null(out[[slot]])) {
      stop(sprintf("Only one %s term is supported per model.", label),
           call. = FALSE)
    }
  }

  for (t in terms) {
    spec  <- t$spec
    procs <- t$processes
    if (any(procs > 2L)) {
      late <- proc_names[procs[procs > 2L]]
      stop(sprintf(
        "Structured term `%s` enters the '%s' predictor; structured effects ",
        spec$label %||% class(spec)[1], paste(late, collapse = "', '")),
        "are supported on the occupancy/state and detection predictors only.",
        call. = FALSE)
    }
    shared <- c(1L %in% procs, 2L %in% procs)

    if (inherits(spec, "tobs_spatial")) {
      one("spatial", "spatial"); spec$shared <- shared; out$spatial <- spec
    } else if (inherits(spec, "tobs_temporal")) {
      one("temporal", "temporal"); spec$shared <- shared; out$temporal <- spec
    } else if (inherits(spec, "tobs_svc")) {
      one("svc", "svc"); spec$shared <- shared; out$svc <- spec
    } else if (inherits(spec, "tobs_latent")) {
      one("latent", "latent"); spec$shared <- any(shared); out$latent <- spec
    } else if (inherits(spec, "tobs_re")) {
      spec$shared <- shared
      re_list[[length(re_list) + 1L]] <- spec
    }
  }
  if (length(re_list)) out$re <- re_list
  out
}

# Build multi-term RE spec for C++ from list of tobs_re objects
build_re_spec <- function(re_list, model) {
  n_terms <- length(re_list)
  N <- if (!is.null(model$N)) model$N else model$n_sites

  groups <- vector("list", n_terms)
  n_groups <- integer(n_terms)
  has_slopes <- FALSE
  n_coefs <- integer(n_terms)
  has_intercept <- rep(TRUE, n_terms)
  slope_matrices <- vector("list", n_terms)

  # Aggregate shared across all terms (union)
  max_proc <- length(model$process_info)
  shared <- rep(FALSE, max_proc)

  for (t in seq_len(n_terms)) {
    re <- re_list[[t]]

    # Group codes were resolved when the re() term was constructed.
    grp <- as.integer(re$group_idx)
    if (length(grp) != N) {
      stop(sprintf("RE group vector has %d elements but model has %d observations",
                   length(grp), N))
    }
    groups[[t]] <- grp
    n_groups[t] <- if (!is.null(re$n_groups)) re$n_groups else max(grp)

    # Sharing
    for (k in seq_along(re$shared)) {
      if (re$shared[k]) shared[k] <- TRUE
    }

    # Slopes. The covariate was resolved at construction: a numeric column or
    # a multi-column matrix (a bare symbol / cbind() in the formula), or column
    # names for a direct re() call. Each column is one slope; n_coefs is the
    # slope count plus the implicit intercept unless the block is slope-only.
    if (re$type == "slope" && !is.null(re$covariate)) {
      has_slopes <- TRUE
      Xs <- .tobs_re_slope_matrix(re$covariate, model$data)
      has_intercept[t] <- isTRUE(re$intercept)
      if (ncol(Xs) == 0L) {
        stop("re(): random slope resolved to zero covariate columns.",
             call. = FALSE)
      }
      n_coefs[t] <- (if (has_intercept[t]) 1L else 0L) + ncol(Xs)
      slope_matrices[[t]] <- Xs
    } else {
      n_coefs[t] <- 1L  # Intercept only
    }
  }

  spec <- list(
    n_terms = n_terms,
    groups = groups,
    n_groups = n_groups,
    shared = shared,
    re_has_intercept = as.integer(has_intercept),
    sigma_re_scale = re_list[[1]]$sigma_scale
  )

  if (has_slopes) {
    spec$has_slopes <- TRUE
    spec$n_coefs <- n_coefs
    spec$slope_matrices <- slope_matrices
    # Per-term correlation flag (0/1): a term is correlated only if it asked
    # for it and has more than one coefficient. The engine reads this per term
    # (re_correlated[t] / re_n_chol[t]), so mixed `|` / `||` blocks are honoured.
    spec$correlated <- as.integer(
      vapply(re_list, function(r) isTRUE(r$correlated), logical(1)) &
      (n_coefs > 1L))
  }

  spec
}

# Resolve an re() slope covariate to a numeric [N x n_slopes] design matrix.
# Accepts column name(s) (resolved via model.matrix, intercept dropped), a
# numeric matrix (e.g. cbind(x, z) from bar desugaring), or a single numeric
# vector. Column names are carried through for ranef() labelling.
.tobs_re_slope_matrix <- function(cov, data) {
  if (is.character(cov)) {
    X <- stats::model.matrix(stats::reformulate(cov), data)
    icpt <- match("(Intercept)", colnames(X))
    if (!is.na(icpt)) X <- X[, -icpt, drop = FALSE]
    return(X)
  }
  if (is.matrix(cov)) {
    storage.mode(cov) <- "double"
    if (is.null(colnames(cov))) {
      colnames(cov) <- paste0("slope", seq_len(ncol(cov)))
    }
    return(cov)
  }
  m <- matrix(as.numeric(cov), ncol = 1L)
  colnames(m) <- "slope1"
  m
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
