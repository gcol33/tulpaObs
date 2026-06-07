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
                            sigma.beta = 10, sigma.re.scale = 1,
                            max.iter = 100L, tol = 1e-4, damping = 0.7,
                            n.iter = 2000, n.warmup = 1000, n.thin = 1L,
                            n.chains = 1L, n.threads = 1L,
                            max.treedepth = 10, adapt.delta = 0.8, seed = 42,
                            approx = c("gaussian_laplace", "simplified_laplace"),
                            correction = "none",
                            n.gibbs = 10L, n.imputations = 20L,
                            re.aghq = TRUE, n.quad = 9L, re.lkj = 1.5,
                            K.max = NULL, mixture = "poisson",
                            verbose = TRUE, ...) {

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

  # N-mixture abundance: a closed-form marginal Laplace fit (tulpa owns the
  # likelihood). No EM, no NUTS yet. `method` is "laplace" (non-spatial) or
  # "nested_laplace" (areal spatial offset). Reuses the per-process autoscaler;
  # the coefficient covariance is transformed back to natural scale alongside
  # the means / draws.
  if (identical(model$model_type, "nmix")) {
    # NUTS: sample the exact coefficient posterior of the non-spatial N-mixture
    # via the in-tree C++ FullGradFn over the closed-form marginal (R/abun_nuts.R).
    # Spatial / RE / temporal terms are not yet wired on the sampler (#51).
    if (identical(method, "nuts")) {
      if (!is.null(spatial)) {
        stop("method = \"nuts\" for abun() is the non-spatial N-mixture sampler; ",
             "a spatial term (icar()/bym2()/car_proper()) on the abundance arm ",
             "fits under method = \"nested_laplace\".", call. = FALSE)
      }
      if (!is.null(temporal)) {
        stop("method = \"nuts\" for abun() does not yet support temporal terms ",
             "(#51); use method = \"laplace\".", call. = FALSE)
      }
      # Random effects (tulpaObs#51): a single intercept RE on one arm samples
      # under NUTS (non-centered per-site offset + log_sigma hyperparameter).
      # Slopes / multi-term / both-arm RE stay on the AGHQ Laplace path.
      fit <- .tobs_fit_abun_nuts(
        fit_model, mixture = mixture, K_max = K.max, sigma.beta = sigma.beta,
        re = re,
        n.iter = n.iter, n.warmup = n.warmup, n.chains = n.chains,
        max.treedepth = max.treedepth, adapt.delta = adapt.delta,
        seed = seed, verbose = verbose)
      fit <- .unscale_fit_per_process(fit, scales, process_info)
      fit$vcov   <- .unscale_vcov(fit$vcov, scales, process_info)
      fit$model  <- model
      fit$intercepts <- compute_intercepts(model, fit$means)
      return(fit)
    }
    nmix_method <- if (is.null(spatial)) "laplace" else "nested_laplace"
    fit <- .tobs_fit_nmix(fit_model, method = nmix_method, spatial = spatial,
                          temporal = temporal, re = re, priors = priors,
                          mixture = mixture, K_max = K.max,
                          max_iter = max.iter, tol = tol,
                          n_quad = n.quad, lkj_eta = re.lkj,
                          sigma_beta = sigma.beta,
                          verbose = verbose)
    fit <- .unscale_fit_per_process(fit, scales, process_info)
    fit$vcov   <- .unscale_vcov(fit$vcov, scales, process_info)
    fit$model  <- model
    fit$intercepts <- compute_intercepts(model, fit$means)
    return(fit)
  }

  # Removal sampling: the sequential-depletion abundance marginal (its latent N
  # summed out in closed form, like the N-mixture). Non-spatial fixed effects
  # only this round (gcol33/tulpaObs#51); "laplace" or "nuts".
  if (identical(model$model_type, "removal")) {
    if (!is.null(spatial) || !is.null(temporal)) {
      stop("removal() currently supports non-spatial fixed effects only; a ",
           "spatial / temporal term is not yet wired. (#51)", call. = FALSE)
    }
    if (!is.null(re) && !identical(method, "nuts")) {
      stop("removal() random effects fit under method = \"nuts\" (a single ",
           "intercept RE on one arm, tulpaObs#51); the Laplace path is ",
           "non-spatial fixed effects only.", call. = FALSE)
    }
    if (identical(method, "nuts")) {
      fit <- .tobs_fit_removal_nuts(
        fit_model, mixture = mixture, K_max = K.max, sigma.beta = sigma.beta,
        re = re,
        n.iter = n.iter, n.warmup = n.warmup, n.chains = n.chains,
        max.treedepth = max.treedepth, adapt.delta = adapt.delta,
        seed = seed, verbose = verbose)
    } else {
      fit <- .tobs_fit_removal(fit_model, mixture = mixture, K_max = K.max,
                               max_iter = max.iter, tol = tol, verbose = verbose)
    }
    fit <- .unscale_fit_per_process(fit, scales, process_info)
    fit$vcov   <- .unscale_vcov(fit$vcov, scales, process_info)
    fit$model  <- model
    fit$intercepts <- compute_intercepts(model, fit$means)
    return(fit)
  }

  # Distance sampling: the binned multinomial-over-N marginal (its latent N
  # summed out in closed form, like the N-mixture). Non-spatial fixed effects
  # only this round (gcol33/tulpaObs#51); "laplace" or "nuts".
  if (identical(model$model_type, "distance")) {
    if (!is.null(spatial) || !is.null(re) || !is.null(temporal)) {
      stop("distance() currently supports non-spatial fixed effects only; a ",
           "spatial / random-effect / temporal term is not yet wired. (#51)",
           call. = FALSE)
    }
    if (identical(method, "nuts")) {
      fit <- .tobs_fit_distance_nuts(
        fit_model, mixture = mixture, K_max = K.max, sigma.beta = sigma.beta,
        n.iter = n.iter, n.warmup = n.warmup, n.chains = n.chains,
        max.treedepth = max.treedepth, adapt.delta = adapt.delta,
        seed = seed, verbose = verbose)
    } else {
      fit <- .tobs_fit_distance(fit_model, mixture = mixture, K_max = K.max,
                                max_iter = max.iter, tol = tol, verbose = verbose)
    }
    fit <- .unscale_fit_per_process(fit, scales, process_info)
    fit$vcov   <- .unscale_vcov(fit$vcov, scales, process_info)
    fit$model  <- model
    fit$intercepts <- compute_intercepts(model, fit$means)
    return(fit)
  }

  # Open-population (Dail-Madsen) N-mixture: the latent abundance sequence summed
  # out by an exact HMM forward recursion (not closed form). Non-spatial fixed
  # effects only this round (gcol33/tulpaObs#51); "laplace" or "nuts".
  if (identical(model$model_type, "dyn_abun")) {
    if (!is.null(spatial) || !is.null(re) || !is.null(temporal)) {
      stop("dyn_abun() currently supports non-spatial fixed effects only; a ",
           "spatial / random-effect / temporal term is not yet wired. (#51)",
           call. = FALSE)
    }
    if (identical(method, "nuts")) {
      fit <- .tobs_fit_dyn_abun_nuts(
        fit_model, sigma.beta = sigma.beta,
        n.iter = n.iter, n.warmup = n.warmup, n.chains = n.chains,
        max.treedepth = max.treedepth, adapt.delta = adapt.delta,
        seed = seed, verbose = verbose)
    } else {
      fit <- .tobs_fit_dyn_abun(fit_model, max_iter = 300L, tol = 1e-8,
                                verbose = verbose)
    }
    fit <- .unscale_fit_per_process(fit, scales, process_info)
    fit$vcov   <- .unscale_vcov(fit$vcov, scales, process_info)
    fit$model  <- model
    fit$intercepts <- compute_intercepts(model, fit$means)
    return(fit)
  }

  # False-positive occupancy: the Miller et al. (2011) multistate marginal (its
  # latent occupancy z summed out in closed form). Non-spatial fixed effects only
  # this round (gcol33/tulpaObs#51); "laplace" (analytic-gradient BFGS over the
  # exact marginal) or "nuts".
  if (identical(model$model_type, "fp_occu")) {
    if (!is.null(spatial) || !is.null(re) || !is.null(temporal)) {
      stop("fp_occu() currently supports non-spatial fixed effects only; a ",
           "spatial / random-effect / temporal term is not yet wired. (#51)",
           call. = FALSE)
    }
    if (identical(method, "nuts")) {
      fit <- .tobs_fit_fp_occu_nuts(
        fit_model, sigma.beta = sigma.beta,
        n.iter = n.iter, n.warmup = n.warmup, n.chains = n.chains,
        max.treedepth = max.treedepth, adapt.delta = adapt.delta,
        seed = seed, verbose = verbose)
    } else {
      fit <- .tobs_fit_fp_occu(fit_model, max_iter = 500L, tol = 1e-8,
                               sigma.beta = NULL, verbose = verbose)
    }
    fit <- .unscale_fit_per_process(fit, scales, process_info)
    fit$vcov   <- .unscale_vcov(fit$vcov, scales, process_info)
    fit$model  <- model
    fit$intercepts <- compute_intercepts(model, fit$means)
    return(fit)
  }

  if (method == "laplace") {
    fit <- .tobs_laplace(fit_model, spatial = spatial, re = re,
                         priors = priors,
                         sigma_beta = sigma.beta,
                         max_iter = max.iter, tol = tol, damping = damping,
                         approx = approx,
                         correction = correction,
                         n_imputations = n.imputations, n_gibbs = n.gibbs,
                         seed = seed,
                         re_aghq = re.aghq, n_quad = n.quad, lkj_eta = re.lkj,
                         verbose = verbose)
    fit <- .unscale_fit_per_process(fit, scales, process_info)
    fit$model      <- model
    fit$intercepts <- compute_intercepts(model, fit$means)
    return(fit)
  }

  if (method == "nested_laplace") {
    # Nested-Laplace path: single-season, integrated, or dynamic
    # occupancy. The driver builds a multi-block latent prior from spatial +
    # temporal + re and attaches it to the state ("occ") M-step block, which
    # tulpa::tulpa_em_laplace() routes through tulpa::tulpa_nested_laplace().
    nl_max_iter <- min(as.integer(max.iter), 25L)
    # INLA NA-response prediction: single-season sites with an all-missing
    # detection history are held out of the likelihood (n_trials = 0) but stay
    # in the design so the latent field interpolates their occupancy from
    # neighbours / shared structure (`.tobs_heldout_sites()`).
    heldout_state <- .tobs_heldout_sites(model)
    fit <- .tobs_em_nested_laplace(
      model    = fit_model,
      spatial  = spatial,
      temporal = temporal,
      re       = re,
      priors   = priors,
      sigma_beta = sigma.beta,
      max_iter = nl_max_iter,
      tol      = tol,
      damping  = damping,
      heldout_state = heldout_state,
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
    sigma_beta = sigma.beta,
    n_iter = as.integer(n.iter),               # total iterations (incl. warmup)
    n_warmup = as.integer(n.warmup),
    max_treedepth = as.integer(max.treedepth),
    adapt_delta = adapt.delta,
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
    spec$sigma_re_scale <- sigma.re.scale
    spec$re_shared_occ <- TRUE
  }

  # ---- Spatial ----
  if (!is.null(spatial)) {
    spatial_params <- build_spatial_params(spatial, model$n_sites)
    spec$spatial_params <- spatial_params
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

  # ---- Run NUTS: one or more chains, pooled in R (see nuts_chains.R) ----
  fit <- .tobs_run_chains(spec, n.chains = n.chains, n.thin = n.thin,
                          n.threads = n.threads, verbose = verbose)

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
  # Resolved per-chain seeds (chain c used seed + c - 1) for reproducibility.
  fit$seeds <- as.integer(seed) + seq_len(as.integer(n.chains)) - 1L
  # Expose process_info at top level for tulpa generic S3 methods
  fit$process_info <- model$process_info
  # Cross-chain convergence diagnostics (Rhat / bulk + tail ESS). tulpa owns
  # the generic estimator (mcmc_diagnostics); it reads fit$draws + fit$chain_id.
  # Computed on the named, natural-scale draws (Rhat / ESS are scale-invariant).
  fit$convergence <- tryCatch(tulpa::mcmc_diagnostics(fit),
                              error = function(e) NULL)
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
      .tobs_reject_weighted_spatial(spec, "occupancy/abundance spatial")
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

# Compute back-transformed intercepts on the response scale. The link is read
# per process (`pi$link`, default "logit"): logit-link processes (occupancy /
# detection probabilities) back-transform with plogis(); the log-link
# abundance process (lambda) with exp(), giving the mean per-site abundance.
compute_intercepts <- function(model, means) {
  result <- list()
  offset <- 0
  for (pi in model$process_info) {
    b0 <- means[offset + 1]
    link <- pi$link %||% "logit"
    result[[pi$name]] <- if (identical(link, "log")) exp(b0) else plogis(b0)
    offset <- offset + pi$p
  }
  result
}

# Transform a coefficient covariance from the autoscaled design back to the
# natural scale. The per-process scaling is a block-diagonal linear map T (each
# block's `.scale_transform()`), so vcov_natural = T vcov_scaled T'. Cross-arm
# blocks are transformed exactly (block-diagonal T preserves them).
.unscale_vcov <- function(vcov, scales, process_info) {
  if (is.null(vcov) || is.null(scales) || is.null(process_info)) return(vcov)
  p_tot <- nrow(vcov)
  Tfull <- diag(p_tot)
  off <- 0L
  for (k in seq_along(process_info)) {
    p_k <- as.integer(process_info[[k]]$p)
    if (p_k == 0L) next
    sc <- scales[[k]]
    if (!is.null(sc) && length(sc$cols) > 0L) {
      idx <- off + seq_len(p_k)
      Tfull[idx, idx] <- .scale_transform(sc)
    }
    off <- off + p_k
  }
  out <- Tfull %*% vcov %*% t(Tfull)
  dimnames(out) <- dimnames(vcov)
  out
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
