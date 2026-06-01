# =============================================================================
# abun.R — N-mixture abundance family (Royle 2004)
#
# Latent abundance N_i ~ Poisson(lambda_i), counts y_ij | N_i ~ Binomial(N_i,
# p_ij). The marginal likelihood integrates N out exactly (closed-form sum to
# K_max), so there is no EM: the package-internal `nmix_laplace()` family fits
# the marginal directly with analytical gradients and observed Fisher curvature.
# This file owns the family interface — data binding, the thin call into the
# (spatial) N-mixture Laplace fitter, and the `tobs_fit` wrapper.
#
#   .tobs_build_abun()   data binder -> model_type = "nmix"
#   .tobs_fit_nmix()     dispatch to the package's (spatial) N-mixture Laplace
#   build_nmix_fit()     wrap a nmix_fit into a tobs_fit
# =============================================================================


# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# Bind an N-mixture model. The abundance predictor is site-level
# (X_lambda, n_sites rows); the detection predictor is long-form (one row per
# observed visit), so site-level detection covariates are replicated across a
# site's visits and visit-level covariates (X_det_visit, site-major) are
# stacked on. Missing visits (NA in `y`) drop out of the long form, which is
# exactly the ragged-visit handling tulpa's marginal-sum likelihood expects.
.tobs_build_abun <- function(abund_formula, det_formula, data, y,
                             det_visit_formula = NULL, det_visit_data = NULL) {
  if (!is.matrix(y)) {
    stop("y must be a matrix (n_sites x max_visits) of integer counts.",
         call. = FALSE)
  }
  if (nrow(y) != nrow(data)) {
    stop(sprintf("y has %d rows but data has %d rows", nrow(y), nrow(data)),
         call. = FALSE)
  }
  y_int <- matrix(as.integer(round(y)), nrow(y), ncol(y))
  if (any(y_int < 0L, na.rm = TRUE)) {
    stop("y must contain nonnegative integer counts (or NA for unobserved ",
         "visits).", call. = FALSE)
  }
  n_sites    <- nrow(y_int)
  max_visits <- ncol(y_int)

  # Strip structured terms (spatial / re / temporal) into `bind$terms`; the
  # fixed-effect formulas build the design matrices. Mirrors .tobs_build_single
  # so icar()/bym2() on the abundance formula are picked up by
  # .tobs_structures_from_model() exactly as for occupancy.
  bind       <- .tobs_bind_formulas(list(lambda = abund_formula, p = det_formula),
                                    data)
  X_lambda   <- model.matrix(bind$fe$lambda, data)
  X_det_site <- model.matrix(bind$fe$p, data)

  X_det_visit <- NULL
  if (!is.null(det_visit_formula) && !is.null(det_visit_data)) {
    mf <- stats::model.frame(det_visit_formula, det_visit_data,
                             na.action = stats::na.pass)
    X_det_visit <- stats::model.matrix(det_visit_formula, mf)
    X_det_visit[is.na(X_det_visit)] <- 0
    expected_rows <- n_sites * max_visits
    if (nrow(X_det_visit) != expected_rows) {
      stop(sprintf("det_visit_data must have %d rows (n_sites * max_visits), got %d",
                   expected_rows, nrow(X_det_visit)), call. = FALSE)
    }
  }

  # Long form in site-major order (site varies slowest), dropping NA visits.
  valid_t    <- as.vector(t(!is.na(y_int)))                       # site-major
  site_mat   <- matrix(seq_len(n_sites), n_sites, max_visits)
  visit_mat  <- matrix(seq_len(max_visits), n_sites, max_visits, byrow = TRUE)
  y_long     <- as.vector(t(y_int))[valid_t]
  site_idx   <- as.vector(t(site_mat))[valid_t]
  visit_idx  <- as.vector(t(visit_mat))[valid_t]
  if (length(y_long) == 0L) {
    stop("y contains no observed visits (all NA).", call. = FALSE)
  }

  X_p <- X_det_site[site_idx, , drop = FALSE]
  if (!is.null(X_det_visit)) {
    visit_row <- (site_idx - 1L) * max_visits + visit_idx
    X_p <- cbind(X_p, X_det_visit[visit_row, , drop = FALSE])
  }

  det_coef_names <- colnames(X_det_site)
  if (!is.null(X_det_visit)) {
    det_coef_names <- c(det_coef_names, colnames(X_det_visit))
  }

  # X_processes carries (X_lambda, X_p) so the per-process autoscaler and
  # unscaler in occu_fit.R operate on both arms (the scaling is link-agnostic,
  # a linear reparameterization of the design columns).
  structure(list(
    model_type = "nmix",
    y          = y_int,
    y_long     = as.integer(y_long),
    site_idx   = as.integer(site_idx),
    visit_idx  = as.integer(visit_idx),
    X_processes = list(X_lambda, X_p),
    formulas   = list(lambda = bind$fe$lambda, det = bind$fe$p),
    structured_terms = bind$terms,
    data       = data,
    n_sites    = n_sites,
    max_visits = max_visits,
    process_info = list(
      list(name = "lambda", p = ncol(X_lambda), coef_names = colnames(X_lambda),
           link = "log"),
      list(name = "p",      p = ncol(X_p),      coef_names = det_coef_names,
           link = "logit")
    ),
    mean_count = mean(y_long),
    max_count  = if (length(y_long)) max(y_long) else 0L
  ), class = "tobs_model")
}


# ---------------------------------------------------------------------------
# Fitter: dispatch to tulpa's marginal N-mixture Laplace engine
# ---------------------------------------------------------------------------

# Called from `.tobs_fit_model()` for `model$model_type == "nmix"`, after the
# per-process design autoscaling. `model` here is the autoscaled model;
# unscaling back to the natural coefficient scale happens in the caller.
.tobs_fit_nmix <- function(model, method = c("laplace", "nested_laplace"),
                           spatial = NULL, temporal = NULL, re = NULL,
                           priors = NULL, mixture = "poisson",
                           K_max = NULL, max_iter = 100L, tol = 1e-6,
                           n_quad = 9L, lkj_eta = 1.5, sigma_beta = 10,
                           verbose = TRUE) {
  method <- match.arg(method)

  # tulpaObs vocabulary ("poisson" / "negbin") -> tulpa's mixing-distribution
  # code ("P" / "NB"). The NB marginal sum and its analytic dispersion score
  # live in tulpa's N-mixture kernel; both the non-spatial and the areal
  # nested-Laplace paths accept it.
  mix_code <- switch(mixture,
                     poisson = "P", negbin = "NB",
                     stop(sprintf("Unknown mixture '%s' (use \"poisson\" or \"negbin\").",
                                  mixture), call. = FALSE))

  # Capability gates. tulpa's N-mixture engine carries fixed effects plus an
  # optional areal spatial offset and (since gcol33/tulpaObs#13) site-level
  # random effects on one arm; temporal is the open upstream extension. Error
  # with a pointer rather than silently dropping the requested structure.
  if (!is.null(temporal)) {
    stop("A temporal term on N-mixture abundance is not yet supported.",
         call. = FALSE)
  }
  if (!is.null(re) && !is.null(spatial)) {
    stop("Random effects together with an areal spatial term on N-mixture ",
         "abundance are not yet supported. Use one or the other.",
         call. = FALSE)
  }
  if (!is.null(re) && !is.null(model$X_det_visit)) {
    stop("Random effects with visit-level detection covariates ",
         "(det_visit_formula) are not yet supported; the site-level make_site ",
         "path assumes a per-site detection design.", call. = FALSE)
  }
  if (!is.null(priors) && !isFALSE(priors)) {
    message(".tobs_fit_nmix(): fixed-effect priors are not applied on the ",
            "N-mixture path; tulpa's marginal Laplace fit is unpenalised.")
  }

  X_lambda <- model$X_processes[[1]]
  X_p      <- model$X_processes[[2]]
  y_long   <- model$y_long
  site_idx <- model$site_idx

  # Random-effect path: warm-start the betas with the no-RE Laplace fit, then
  # refine on the exact-marginal AGHQ likelihood with a variance-component
  # update. n_quad = 1 is the joint Laplace (the small-cluster sigma
  # attenuation regime); n_quad > 1 debiases it.
  if (!is.null(re)) {
    return(.tobs_fit_nmix_re(model, re = re,
                             mixture = mix_code, K_max = K_max,
                             max_iter = max_iter, tol = tol,
                             n_quad = n_quad, lkj_eta = lkj_eta,
                             theta_prior_sd = sigma_beta,
                             verbose = verbose))
  }

  if (is.null(spatial)) {
    raw <- nmix_laplace(
      y         = y_long,
      site_idx  = site_idx,
      X_lambda  = X_lambda,
      X_p       = X_p,
      mixture   = mix_code,
      K_max     = K_max,
      max_iter  = as.integer(max_iter),
      tol       = as.numeric(tol),
      verbose   = isTRUE(verbose)
    )
    return(build_nmix_fit(raw, model, spatial = NULL))
  }

  # Areal spatial offset on the abundance arm via tulpa's nested-Laplace
  # N-mixture fitters (icar / bym2 / car_proper).
  raw <- .tobs_fit_nmix_spatial(model, spatial, X_lambda, X_p, y_long, site_idx,
                                mixture = mix_code, K_max = K_max,
                                max_iter = max_iter, tol = tol,
                                verbose = verbose)
  build_nmix_fit(raw, model, spatial = spatial)
}


# ---------------------------------------------------------------------------
# Random-effect path: warm-start + AGHQ refinement
# ---------------------------------------------------------------------------

# Abundance(lambda)/detection(p) arm split (N-mixture AGHQ path). A shared term
# is rejected -- the make_site path integrates one arm at a time; a cross-arm RE
# would need a joint two-arm oracle. Shares .tobs_split_re_arms (R/em_laplace_re.R).
.tobs_nmix_re_split_arms <- function(re_list, model) {
  .tobs_split_re_arms(
    re_list, model, "lambda", "p",
    paste0("A random effect shared across the abundance (lambda) and ",
           "detection (p) arms is not supported on the N-mixture AGHQ path."))
}


# Fit an N-mixture with grouped random effects. Warm-starts the betas with the
# no-RE Laplace fit (so the AGHQ joint optimization sits near the MLE rather
# than starting from origin), seeds Sigma at a diagonal 0.25 per term, and
# routes through .tobs_nmix_re_aghq() for the variance-component fit.
.tobs_fit_nmix_re <- function(model, re, mixture = "P", K_max = NULL,
                              max_iter = 100L, tol = 1e-6, verbose = TRUE,
                              n_quad = 1L, lkj_eta = 1.5,
                              theta_prior_sd = 100) {
  if (inherits(re, "tobs_re")) re <- list(re)
  arms <- .tobs_nmix_re_split_arms(re, model)
  # v1: single arm. Both populated -> reject (cross-arm AGHQ needs make_group
  # with a joint two-arm oracle, a separate engine path).
  if (length(arms$lambda) && length(arms$p)) {
    stop("Random effects on BOTH the abundance and detection arms in one ",
         "N-mixture fit are not yet supported; the AGHQ path integrates one ",
         "arm at a time. Put the RE on lambda OR p, not both.",
         call. = FALSE)
  }
  design <- if (length(arms$lambda)) arms$lambda else arms$p
  arm    <- if (length(arms$lambda)) "lambda" else "p"

  X_lambda <- model$X_processes[[1]]
  X_p      <- model$X_processes[[2]]
  y_long   <- model$y_long
  site_idx <- model$site_idx

  # Warm-start the betas (and log_r under NB) from the no-RE fit. A well-placed
  # start matters more here than for the closed-form non-spatial fit: the AGHQ
  # objective wraps n_groups inner Newton mode-finds at every theta candidate,
  # so a far start can spend an extra 5-10 outer iterations on betas the no-RE
  # fit already locates cheaply.
  warm <- tryCatch(
    nmix_laplace(y = y_long, site_idx = site_idx,
                 X_lambda = X_lambda, X_p = X_p,
                 mixture = mixture, K_max = K_max,
                 max_iter = as.integer(max_iter), tol = as.numeric(tol),
                 verbose = FALSE),
    error = function(e) NULL)
  beta_lambda_init <- if (!is.null(warm)) warm$beta_lambda
                      else c(log(max(mean(y_long), 0.1)), rep(0, ncol(X_lambda) - 1L))
  beta_p_init      <- if (!is.null(warm)) warm$beta_p else rep(0, ncol(X_p))
  r_init <- if (!is.null(warm) && identical(mixture, "NB") &&
                is.finite(warm$r %||% NA_real_)) as.numeric(warm$r) else 10

  Sigma_init <- lapply(design, function(d) diag(0.25, d$n_coefs))
  b_init <- numeric(sum(vapply(design,
                               function(d) as.integer(d$n_groups * d$n_coefs),
                               integer(1))))

  ref <- .tobs_nmix_re_aghq(model, design,
                            beta_lambda = beta_lambda_init,
                            beta_p      = beta_p_init,
                            Sigma_list  = Sigma_init,
                            b           = b_init,
                            mixture     = mixture,
                            r_init      = r_init,
                            K_max       = K_max,
                            n_quad      = as.integer(n_quad),
                            lkj_eta     = lkj_eta,
                            theta_prior_sd = theta_prior_sd,
                            max_iter    = as.integer(max_iter),
                            verbose     = isTRUE(verbose))
  if (is.null(ref) || !isTRUE(ref$ok)) {
    stop("N-mixture AGHQ random-effect refinement did not produce a usable ",
         "fit (singular marginal Hessian or non-finite optimum). Try a ",
         "different K_max or simplify the RE structure.", call. = FALSE)
  }

  # Assemble a `raw`-shaped object build_nmix_fit consumes. The AGHQ marginal
  # Hessian over theta is the joint (lambda coefs, p coefs[, log_r]) covariance
  # the engine returned -- the same coordinate set the non-spatial fit's
  # observed-Fisher inverse covers.
  raw <- list(
    mixture     = mixture,
    beta_lambda = ref$beta_lambda,
    beta_p      = ref$beta_p,
    log_r       = ref$log_r,
    r           = ref$r,
    vcov        = ref$vcov,
    log_lik     = ref$log_marginal,
    converged   = ref$converged,
    K_max       = K_max
  )
  re_post <- list(arm = arm, design = design,
                  Sigma_list = ref$Sigma_list,
                  b = ref$b, b_var = ref$b_var,
                  n_quad = ref$n_quad, lkj_eta = ref$lkj_eta)
  build_nmix_fit(raw, model, spatial = NULL, re_post = re_post)
}


# Areal spatial N-mixture: one spatial unit per site (identity map). Routes to
# the tulpa fitter matching the spatial term type.
.tobs_fit_nmix_spatial <- function(model, spatial, X_lambda, X_p, y_long,
                                   site_idx, mixture = "P", K_max, max_iter,
                                   tol, verbose) {
  .tobs_reject_weighted_spatial(spatial, "N-mixture abundance spatial")
  # Continuous Matern field: the SPDE FEM precision Q(range, sigma) and the
  # mesh projection A are threaded through the same outer-grid integrator the
  # areal path uses, with (range, sigma) as the grid axes.
  if (identical(spatial$type, "spde")) {
    return(nmix_laplace_spde(
      y = y_long, site_idx = site_idx, X_lambda = X_lambda, X_p = X_p,
      spatial = spatial, mixture = mixture, K_max = K_max,
      max_iter = max_iter, tol = tol, verbose = verbose))
  }
  if (identical(spatial$type, "gp") || identical(spatial$type, "multiscale_gp")) {
    stop(sprintf(
      "N-mixture abundance supports the continuous field via spde(); the dense ",
      "GP term '%s' is not wired. Use spde() for a mesh-based continuous Matern ",
      "field, or icar() / bym2() / car_proper() for an areal field.",
      spatial$type), call. = FALSE)
  }
  if (!spatial$type %in% c("icar", "bym2", "car_proper")) {
    stop(sprintf(
      "N-mixture abundance supports the areal spatial terms icar() / bym2() / %s",
      "car_proper() and the continuous mesh field spde() under method = ",
      "\"nested_laplace\"; got '%s'. (car() is the improper non-intrinsic CAR; ",
      "use icar() for the intrinsic field.)"),
      spatial$type, call. = FALSE)
  }
  n_sites <- model$n_sites
  if (spatial$n_units != n_sites) {
    stop(sprintf("spatial term has %d units but the model has %d sites; one ",
                 "spatial unit per site is required for N-mixture.",
                 spatial$n_units, n_sites), call. = FALSE)
  }
  # icar/bym2/car_proper precompute the CSR adjacency; fall back to deriving
  # it from the graph for any spatial term that carries only the graph.
  csr <- if (!is.null(spatial$adj_row_ptr)) {
    list(row_ptr = spatial$adj_row_ptr, col_idx = spatial$adj_col_idx,
         n_neighbors = spatial$n_neighbors)
  } else {
    adjacency_to_csr(spatial$graph)
  }
  common <- list(
    y = y_long, site_idx = site_idx, map_site_to_unit = seq_len(n_sites),
    X_lambda = X_lambda, X_p = X_p,
    adj_row_ptr = csr$row_ptr, adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors, n_spatial = spatial$n_units,
    mixture = mixture,
    K_max = K_max, max_iter = as.integer(max_iter), tol = as.numeric(tol),
    verbose = isTRUE(verbose)
  )
  scale_factor <- spatial$scale_factor %||% compute_bym2_scale(spatial$graph)
  switch(
    spatial$type,
    icar = do.call(nmix_laplace_icar, common),
    bym2 = do.call(nmix_laplace_bym2,
                   c(common, list(scale_factor = scale_factor))),
    car_proper = do.call(nmix_laplace_car_proper, common)
  )
}


# ---------------------------------------------------------------------------
# Wrap a tulpa N-mixture fit into a tobs_fit
# ---------------------------------------------------------------------------

# Build the public fit object. Parameter layout is (lambda coefs, p coefs) to
# match `process_info`. Pseudo-draws are drawn from the JOINT (lambda, p)
# covariance so derived quantities (e.g. lambda * detection-corrected counts,
# total abundance) propagate the cross-arm covariance rather than treating the
# arms as independent.
build_nmix_fit <- function(raw, model, spatial = NULL, re_post = NULL) {
  pi_list <- model$process_info
  p_lam   <- pi_list[[1]]$p
  p_p     <- pi_list[[2]]$p
  is_nb   <- identical(raw$mixture, "NB")
  nms <- c(paste0("lambda_", pi_list[[1]]$coef_names),
           paste0("p_",      pi_list[[2]]$coef_names))

  beta_lambda <- if (!is.null(raw$beta_lambda)) raw$beta_lambda else raw$beta_lambda_mean
  beta_p      <- if (!is.null(raw$beta_p))      raw$beta_p      else raw$beta_p_mean
  means <- c(as.numeric(beta_lambda), as.numeric(beta_p))

  # NB dispersion. The non-spatial fit estimates log_r jointly with the betas
  # and returns it as the last vcov coordinate, so it is carried as a model
  # coefficient (reported by coef() / vcov() / confint() with an SE). The
  # spatial fit integrates the NB size r over the outer hyperparameter grid
  # instead (alongside tau / rho / sigma), so it carries no log_r coordinate
  # and is summarized as a grid hyperparameter (r_mean / r_sd) below.
  has_logr <- is_nb && is.null(spatial) && is.finite(raw$log_r %||% NA_real_)
  if (has_logr) {
    nms   <- c(nms, "log_r")
    means <- c(means, as.numeric(raw$log_r))
  }
  names(means) <- nms

  # Joint covariance over (lambda, p[, log_r]). The non-spatial fit returns the
  # full marginal observed-Fisher inverse; the spatial fits return the
  # grid-integrated coefficient covariance over the full (beta_lambda, beta_p)
  # block (the constrained top block of each grid cell's joint Hessian inverse,
  # with the shared field folded out, then law-of-total-covariance integrated --
  # `.nmix_grid_vcov()`). Both carry the cross-arm (lambda, p) covariance, so
  # derived quantities that combine the two arms (expected counts lambda * p,
  # detection-corrected abundance) propagate the correlation.
  vcov <- .nmix_vcov(raw, p_lam, p_p, p_extra = if (has_logr) 1L else 0L)
  rownames(vcov) <- colnames(vcov) <- nms
  sds <- sqrt(pmax(diag(vcov), 0))
  names(sds) <- nms

  # The fixed-effect block is the leading (lambda, p[, log_r]) coordinates; the
  # RE variance-component / BLUP columns appended below are not fixed effects.
  # Recording the count + names lets the generic coefficient table
  # (.fit_fixed_table -> .fixed_draws_mat) restrict coef() / confint() / vcov()
  # to the fixed block instead of defaulting to every draw column.
  n_fixed     <- length(nms)
  fixed_names <- nms

  n_pseudo <- 1000L
  draws <- .rmvn(n_pseudo, means, vcov)
  colnames(draws) <- nms

  # Random-effect block (gcol33/tulpaObs#13): append the variance components
  # (sigma_g_*, cor_g_*_* for a correlated block) and per-group BLUPs to the
  # public parameter layout the same way the occupancy RE path does. The
  # per-group BLUPs carry their AGHQ marginal posterior SD, so their pseudo-draws
  # are Gaussian at the BLUP mean with that SD. The variance components are AGHQ
  # marginal-optimum point estimates whose joint posterior covariance with the
  # fixed effects the marginal-Hessian path does not surface; their SD is NA, so
  # their draws are left NA rather than a fabricated near-degenerate column that
  # would read as "known almost exactly" (NA-on-unavailable, the same rule as the
  # NUTS-only sampler diagnostics, tulpaObs#17). ranef() / summary() read the
  # point estimate and BLUP SE by name; no consumer of the fixed-effect draws
  # touches these trailing columns.
  re_block <- NULL
  if (!is.null(re_post) && length(re_post$design)) {
    re_block <- .tobs_re_param_block(list(design = re_post$design,
                                          b      = re_post$b,
                                          b_var  = re_post$b_var,
                                          Sigma  = re_post$Sigma_list))
    re_means <- re_block$means; re_sds <- re_block$sds; re_names <- re_block$names
    means <- c(means, re_means); sds <- c(sds, re_sds)
    nms   <- c(nms, re_names)
    names(means) <- nms; names(sds) <- nms
    n_re <- length(re_means)
    re_draws <- matrix(NA_real_, n_pseudo, n_re)
    for (j in seq_len(n_re)) {
      if (is.finite(re_sds[j]))
        re_draws[, j] <- stats::rnorm(n_pseudo, re_means[j], re_sds[j])
    }
    colnames(re_draws) <- re_names
    draws <- cbind(draws, re_draws)
  }

  spatial_field <- raw$phi_mean %||% raw$z_mean %||% raw$u_mean %||% raw$v_mean
  hyper <- .nmix_hyper(raw)

  # NB dispersion summary on the natural (r) scale, for simulate() and print().
  # Non-spatial: r = exp(log_r), SE by the delta method r * SE(log_r). Spatial:
  # r_mean / r_sd straight from the grid weights. Poisson -> NULL.
  dispersion <- NULL
  if (is_nb) {
    if (has_logr) {
      se_logr <- sqrt(pmax(vcov["log_r", "log_r"], 0))
      dispersion <- list(r = as.numeric(raw$r %||% exp(raw$log_r)),
                         log_r = as.numeric(raw$log_r),
                         r_sd = as.numeric(exp(raw$log_r) * se_logr))
    } else if (is.finite(raw$r_mean %||% NA_real_)) {
      dispersion <- list(r = as.numeric(raw$r_mean), log_r = NA_real_,
                         r_sd = as.numeric(raw$r_sd))
    }
  }

  # The non-spatial engine returns a scalar marginal log-likelihood at the mode
  # (N summed out); the spatial fitters return a per-grid vector, reduced here to
  # the grid-weighted data log-likelihood. Replicate into `log_prob` so the
  # generic logLik.tulpa_fit() (mean(log_prob)) surfaces it -- the deterministic
  # fit has no per-draw likelihood variation. Likewise reduce the per-grid
  # convergence flags to a single summary.
  ll <- raw$log_lik
  if (length(ll) > 1L) {
    w <- raw$weights %||% rep(1 / length(ll), length(ll))
    w[!is.finite(w)] <- 0
    ll <- if (sum(w) > 0) sum(w * ll, na.rm = TRUE) else NA_real_
  }
  ll <- if (length(ll) == 1L) ll else NA_real_

  conv <- raw$converged %||% TRUE
  if (length(conv) > 1L) conv <- all(conv, na.rm = TRUE)
  n_it <- raw$n_iter %||% NA_integer_
  if (length(n_it) > 1L) n_it <- max(n_it, na.rm = TRUE)

  structure(c(list(
    draws = draws, means = means, sds = sds, vcov = vcov,
    n_samples = n_pseudo, n_params = length(means),
    log_prob = rep(ll, n_pseudo),
    N = length(model$y_long)),
    .tobs_na_nuts_diagnostics(n_pseudo),
    list(
    col_names = nms, param_names = nms,
    n_fixed = n_fixed, fixed_names = fixed_names,
    process_info = pi_list,
    model = model, spatial = spatial,
    spatial_field = spatial_field,
    method = if (is.null(spatial)) "laplace" else "nested_laplace",
    log_lik = raw$log_lik %||% NA_real_,
    K_max = raw$K_max,
    mean_N = raw$mean_N, var_N = raw$var_N,
    boundary_weight = raw$boundary_weight,
    nmix_hyper = hyper,
    mixture = if (is_nb) "negbin" else "poisson",
    nmix_dispersion = dispersion,
    re_effects = re_block$re_effects,
    nmix_re = if (!is.null(re_post))
      list(arm = re_post$arm, n_quad = re_post$n_quad,
           lkj_eta = re_post$lkj_eta,
           Sigma_list = re_post$Sigma_list)
      else NULL,
    convergence = list(converged = raw$converged %||% TRUE,
                       n_iter = raw$n_iter %||% NA_integer_)
  )), class = c("tobs_fit", "tulpa_fit"))
}

# Assemble the coefficient covariance from a tulpa N-mixture fit. Non-spatial:
# the engine returns the full joint `vcov` over (lambda, p) and, under NB, the
# trailing `log_r` coordinate (`p_extra = 1`). Spatial: the grid integrator
# returns the full joint (lambda, p) `vcov` (the field-folded cross-arm block,
# from .nmix_grid_vcov over each cell's constrained joint Hessian inverse),
# which is used directly. The per-arm `vcov_lambda` / `vcov_p` block-diagonal
# stitch (and the diagonal from per-arm sds) is only a fallback for an older
# engine that does not surface the full joint block.
.nmix_vcov <- function(raw, p_lam, p_p, p_extra = 0L) {
  p_tot <- p_lam + p_p + p_extra
  if (!is.null(raw$vcov) && all(dim(as.matrix(raw$vcov)) == p_tot)) {
    return(unname(as.matrix(raw$vcov)))
  }
  V <- matrix(0, p_tot, p_tot)
  lam_idx <- seq_len(p_lam)
  p_idx   <- p_lam + seq_len(p_p)
  if (!is.null(raw$vcov_lambda) && all(dim(raw$vcov_lambda) == p_lam)) {
    V[lam_idx, lam_idx] <- as.matrix(raw$vcov_lambda)
  } else if (!is.null(raw$beta_lambda_sd)) {
    diag(V)[lam_idx] <- raw$beta_lambda_sd^2
  }
  if (!is.null(raw$vcov_p) && all(dim(raw$vcov_p) == p_p)) {
    V[p_idx, p_idx] <- as.matrix(raw$vcov_p)
  } else if (!is.null(raw$beta_p_sd)) {
    diag(V)[p_idx] <- raw$beta_p_sd^2
  }
  V
}

# Spatial hyperparameter posterior summary (NULL for the non-spatial fit).
.nmix_hyper <- function(raw) {
  out <- list()
  if (!is.null(raw$sigma_mean)) out$sigma <- c(mean = raw$sigma_mean, sd = raw$sigma_sd)
  if (!is.null(raw$rho_mean))   out$rho   <- c(mean = raw$rho_mean,   sd = raw$rho_sd)
  if (!is.null(raw$tau_mean))   out$tau   <- c(mean = raw$tau_mean,   sd = raw$tau_sd)
  if (!is.null(raw$range_mean)) out$range <- c(mean = raw$range_mean, sd = raw$range_sd)
  # Spatial NB: the size r is integrated over the outer grid alongside the
  # spatial hyperparameters, so it is reported here rather than in vcov.
  if (is.finite(raw$r_mean %||% NA_real_)) out$r <- c(mean = raw$r_mean, sd = raw$r_sd)
  if (length(out) == 0L) NULL else out
}

# Draw from a multivariate normal via the Cholesky of `sigma`; falls back to
# independent normals (diagonal) when `sigma` is not PD (e.g. a non-converged
# observed-info Hessian).
.rmvn <- function(n, mu, sigma) {
  p <- length(mu)
  L <- tryCatch(chol(sigma), error = function(e) NULL)
  z <- matrix(stats::rnorm(n * p), n, p)
  if (is.null(L)) {
    sds <- sqrt(pmax(diag(sigma), 1e-8))
    return(sweep(z * rep(sds, each = n), 2, mu, "+"))
  }
  sweep(z %*% L, 2, mu, "+")
}


# ---------------------------------------------------------------------------
# N-mixture S3 helpers (routed to from methods.R by model_type == "nmix")
# ---------------------------------------------------------------------------

# Per-arm linear predictors at the posterior mean: site-level expected
# abundance `lambda`, per-visit detection `p`, and the posterior mean abundance
# `N` (E[N_i | y], from tulpa's marginal-sum fit).
.tobs_fitted_nmix <- function(object) {
  model <- object$model
  means <- object$means
  p_lam <- model$process_info[[1]]$p
  p_p   <- model$process_info[[2]]$p
  beta_lambda <- means[seq_len(p_lam)]
  beta_p      <- means[p_lam + seq_len(p_p)]
  lambda <- exp(as.vector(model$X_processes[[1]] %*% beta_lambda))
  p_obs  <- plogis(as.vector(model$X_processes[[2]] %*% beta_p))
  list(lambda = lambda, p = p_obs, N = object$mean_N)
}

# Posterior draws of the per-arm linear predictor at a design matrix, returned
# on the response scale ([n_draws x nrow(X.0)]): exp() for the abundance arm,
# plogis() for detection.
.tobs_nmix_response_draws <- function(draws, X.0, beta_offset, p_proc, link) {
  beta <- draws[, beta_offset + seq_len(p_proc), drop = FALSE]
  eta  <- beta %*% t(X.0)
  if (identical(link, "log")) exp(eta) else plogis(eta)
}

# predict() for N-mixture. In-sample -> fitted(); design-matrix / terms modes
# act on the abundance arm by default (`type = "abundance"`) or the detection
# arm (`type = "detection"`).
.tobs_predict_nmix <- function(object, X.0 = NULL, type = "abundance",
                               quantiles = c(0.025, 0.5, 0.975),
                               terms = NULL, n_points = 50L) {
  if (!type %in% c("abundance", "detection")) {
    stop("predict(type=) for N-mixture must be \"abundance\" or \"detection\".",
         call. = FALSE)
  }
  if (is.null(X.0) && is.null(terms)) return(fitted(object))

  model   <- object$model
  draws   <- object$draws
  pi_list <- model$process_info
  proc_idx <- if (type == "detection") 2L else 1L
  link     <- pi_list[[proc_idx]]$link %||% "logit"
  p_proc   <- pi_list[[proc_idx]]$p
  beta_off <- if (proc_idx > 1L) pi_list[[1]]$p else 0L

  if (!is.null(terms)) {
    coef_names <- pi_list[[proc_idx]]$coef_names
    X_orig     <- model$X_processes[[proc_idx]]
    col_idx <- match(terms[1], coef_names)
    if (is.na(col_idx)) {
      stop(sprintf("term '%s' not found in %s coefficients: %s", terms[1],
                   pi_list[[proc_idx]]$name,
                   paste(coef_names, collapse = ", ")), call. = FALSE)
    }
    x_grid <- seq(min(X_orig[, col_idx]), max(X_orig[, col_idx]),
                  length.out = n_points)
    X_pred <- matrix(colMeans(X_orig), n_points, p_proc, byrow = TRUE)
    X_pred[, col_idx] <- x_grid
    pred <- .tobs_nmix_response_draws(draws, X_pred, beta_off, p_proc, link)
    result <- data.frame(
      x = x_grid,
      estimate = colMeans(pred),
      lower = apply(pred, 2, quantile, quantiles[1]),
      upper = apply(pred, 2, quantile, quantiles[length(quantiles)])
    )
    attr(result, "term") <- terms[1]
    attr(result, "process") <- pi_list[[proc_idx]]$name
    class(result) <- c("tobs_prediction", "data.frame")
    return(result)
  }

  if (ncol(X.0) != p_proc) {
    stop(sprintf("X.0 has %d columns but the %s arm has %d coefficients",
                 ncol(X.0), pi_list[[proc_idx]]$name, p_proc), call. = FALSE)
  }
  pred <- .tobs_nmix_response_draws(draws, X.0, beta_off, p_proc, link)
  data.frame(
    mean  = colMeans(pred),
    sd    = apply(pred, 2, sd),
    q2.5  = apply(pred, 2, quantile, quantiles[1]),
    q50   = apply(pred, 2, quantile, quantiles[min(2L, length(quantiles))]),
    q97.5 = apply(pred, 2, quantile, quantiles[length(quantiles)])
  )
}

# simulate() for N-mixture: draw N_i from the abundance mixing distribution
# (Poisson, or NegBin(mu = lambda_i, size = r) under mixture = "negbin"), then
# y_ij ~ Binomial(N_i, p_ij) at the observed visits, respecting the NA pattern.
.tobs_simulate_nmix <- function(object, nsim = 1) {
  model   <- object$model
  draws   <- object$draws
  n_draws <- nrow(draws)
  p_lam   <- model$process_info[[1]]$p
  p_p     <- model$process_info[[2]]$p
  X_lambda <- model$X_processes[[1]]
  X_p      <- model$X_processes[[2]]
  n_sites    <- model$n_sites
  max_visits <- model$max_visits
  site_idx   <- model$site_idx
  visit_idx  <- model$visit_idx
  r_size     <- object$nmix_dispersion$r   # NULL under Poisson

  result <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    di <- sample.int(n_draws, 1L)
    beta_lambda <- draws[di, seq_len(p_lam)]
    beta_p      <- draws[di, p_lam + seq_len(p_p)]
    lambda <- exp(as.vector(X_lambda %*% beta_lambda))
    p_obs  <- plogis(as.vector(X_p %*% beta_p))
    N <- if (!is.null(r_size) && is.finite(r_size)) {
      stats::rnbinom(n_sites, size = r_size, mu = lambda)
    } else {
      stats::rpois(n_sites, lambda)
    }
    y_sim <- matrix(NA_integer_, n_sites, max_visits)
    for (k in seq_along(site_idx)) {
      i <- site_idx[k]; j <- visit_idx[k]
      y_sim[i, j] <- stats::rbinom(1L, N[i], p_obs[k])
    }
    result[[s]] <- y_sim
  }
  if (nsim == 1L) result[[1]] else result
}

# residuals() for N-mixture. A single visit of a Poisson-thinned count is
# marginally y_ij ~ Poisson(lambda_i * p_ij), so the fit is scored on that mean
# with Poisson deviance / Pearson residuals.
.tobs_residuals_nmix <- function(object, type = c("deviance", "pearson",
                                                  "response")) {
  type  <- match.arg(type)
  model <- object$model
  fitv  <- .tobs_fitted_nmix(object)
  lambda <- fitv$lambda
  p_obs  <- fitv$p
  site_idx  <- model$site_idx
  visit_idx <- model$visit_idx
  y_long    <- model$y_long

  mu <- lambda[site_idx] * p_obs
  mu <- pmax(mu, 1e-10)
  r_long <- switch(type,
    response = y_long - mu,
    pearson  = (y_long - mu) / sqrt(mu),
    deviance = {
      d <- 2 * (ifelse(y_long > 0, y_long * log(y_long / mu), 0) - (y_long - mu))
      sign(y_long - mu) * sqrt(pmax(d, 0))
    }
  )
  out <- matrix(NA_real_, model$n_sites, model$max_visits)
  out[cbind(site_idx, visit_idx)] <- r_long
  out
}


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

#' Simulate Royle (2004) N-mixture abundance data
#'
#' Latent abundance `N_i ~ Poisson(lambda_i)` (or, under
#' `mixture = "negbin"`, `N_i ~ NegBin(mean = lambda_i, size = size)` with
#' variance `lambda + lambda^2 / size`) with `log lambda_i = X_lambda
#' beta_lambda`, and replicate counts `y_ij ~ Binomial(N_i, p_i)` with
#' `logit p_i = X_p beta_p` (site-level detection). The returned `y` is an
#' `N x J` integer count matrix suitable for [tobs()] with [abun()].
#'
#' @param N Number of sites (default 100).
#' @param J Number of replicate visits (default 4).
#' @param n_abund_covs Number of abundance covariates (default 2).
#' @param n_det_covs Number of detection covariates (default 1).
#' @param beta_lambda Abundance coefficients `c(intercept, slopes...)` on the
#'   log scale. Default `c(log(3), runif(n_abund_covs, -0.5, 0.5))`.
#' @param beta_p Detection coefficients `c(intercept, slopes...)` on the logit
#'   scale. Default `c(0, runif(n_det_covs, -0.5, 0.5))` (intercept 0 = p 0.5).
#' @param mixture Abundance mixing distribution: `"poisson"` (default) or
#'   `"negbin"` (negative binomial, overdispersed).
#' @param size Negative-binomial size `r` (the `mixture = "negbin"` dispersion;
#'   variance `lambda + lambda^2 / r`). Smaller `r` means more overdispersion;
#'   as `r` grows large the model approaches Poisson. Default 2. Ignored under
#'   Poisson.
#' @param seed Optional random seed.
#' @return A list with `y` (N x J count matrix), `data` (covariate data frame),
#'   and `truth` (the coefficients, per-site `lambda`, `p`, latent `N`, and the
#'   `mixture` / `size` used).
#' @export
simulate_abun <- function(N = 100, J = 4,
                          n_abund_covs = 2, n_det_covs = 1,
                          beta_lambda = NULL, beta_p = NULL,
                          mixture = c("poisson", "negbin"), size = 2,
                          seed = NULL) {
  mixture <- match.arg(mixture)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(beta_lambda)) beta_lambda <- c(log(3), stats::runif(n_abund_covs, -0.5, 0.5))
  if (is.null(beta_p))      beta_p      <- c(0, stats::runif(n_det_covs, -0.5, 0.5))

  abund_covs <- data.frame(matrix(stats::rnorm(N * n_abund_covs), N, n_abund_covs))
  names(abund_covs) <- paste0("abund_cov", seq_len(n_abund_covs))
  det_covs <- data.frame(matrix(stats::rnorm(N * n_det_covs), N, n_det_covs))
  names(det_covs) <- paste0("det_cov", seq_len(n_det_covs))
  data <- cbind(abund_covs, det_covs)

  X_lambda <- stats::model.matrix(~ ., abund_covs)
  X_det    <- stats::model.matrix(~ ., det_covs)

  lambda <- exp(as.vector(X_lambda %*% beta_lambda))
  p      <- plogis(as.vector(X_det %*% beta_p))
  Nlat   <- if (identical(mixture, "negbin")) {
    stats::rnbinom(N, size = size, mu = lambda)
  } else {
    stats::rpois(N, lambda)
  }

  y <- matrix(NA_integer_, N, J)
  for (i in seq_len(N)) {
    y[i, ] <- stats::rbinom(J, Nlat[i], p[i])
  }

  list(
    y = y,
    data = data,
    truth = list(beta_lambda = beta_lambda, beta_p = beta_p,
                 lambda = lambda, p = p, N = Nlat,
                 mixture = mixture,
                 size = if (identical(mixture, "negbin")) size else NA_real_)
  )
}
