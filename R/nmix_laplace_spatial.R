# =============================================================================
# Areal count-marginal abundance via nested Laplace: one engine for the Royle
# (2004) N-mixture and the removal-sampling marginal, with an ICAR, proper-CAR
# or BYM2 field on the abundance arm.
#
#   N_i ~ Poisson(lambda_i) (or NB with size r),
#   log lambda_i = X_lambda_i beta_lambda + phi_u(i),
#   y_i | N_i ~ the family's observation kernel (binomial visits for the
#                N-mixture; depleting passes for removal), logit p from X_p.
#
# Field kinds (Besag, York and Mollie 1991; Cressie 1993; Riebler et al. 2016):
#   icar        phi = z, z ~ ICAR(tau): sum-to-zero per connected component,
#               the outer grid over tau.
#   car_proper  phi = z, z ~ N(0, [tau (D - rho W)]^-1): full rank for rho < 1,
#               no centring; the outer grid over (tau, rho), log|Q(rho)|
#               precomputed per rho.
#   bym2        phi = sigma (sqrt(rho / s) v + sqrt(1 - rho) w), v ICAR
#               (per-component scaled, sum-to-zero), w iid N(0, I), s the
#               Riebler scaling factor (geometric mean of diag(Q^+)); the
#               outer grid over (sigma, rho).
# Under NB the size r is one more outer axis. At each grid point the inner
# Newton finds the joint mode of the coefficients and the field (marginal
# observed information with the Var[N | y] cross-arm correction, falling back to
# the complete-data Fisher block where that is not PSD) and the Laplace
# log-marginal is accumulated; the grid weights integrate the hyperparameters
# (Rue, Martino and Chopin 2009).
#
# The two families share the C++ inner Newton / Laplace driver
# (cpp_nested_laplace_<family>_<kind>) and the grid summaries; they differ in
# the C++ entry and the latent-N floor: the N-mixture count must clear the
# site's largest visit count, the removal count the site's total removal.
# Default axes are one table (.count_spatial_default_grid()), so both families
# integrate the same grids. `K_max` defaults to the family floor + 100; the fit
# warns when the grid-weighted posterior mass on K_max exceeds 1e-4.
# =============================================================================

# The C++ nested-Laplace entry of a family and field kind.
.count_spatial_cpp <- function(family, kind) {
  switch(paste(family, kind),
    "nmix icar"          = cpp_nested_laplace_nmix_icar,
    "nmix car_proper"    = cpp_nested_laplace_nmix_car_proper,
    "nmix bym2"          = cpp_nested_laplace_nmix_bym2,
    "removal icar"       = cpp_nested_laplace_removal_icar,
    "removal car_proper" = cpp_nested_laplace_removal_car_proper,
    "removal bym2"       = cpp_nested_laplace_removal_bym2,
    stop(sprintf("No areal count fitter for family '%s', field '%s'.", family, kind),
         call. = FALSE))
}

# The outer hyperparameter axes of a field kind: each given axis checked, each
# missing one taken from the default table.
.count_spatial_grids <- function(kind, tau_grid = NULL, rho_grid = NULL,
                                 sigma_grid = NULL) {
  dflt <- function(x, axis) if (is.null(x)) .count_spatial_default_grid(axis) else x
  switch(kind,
    icar = list(tau_grid = .count_spatial_check_grid(dflt(tau_grid, "tau_icar"),
                                                     "tau_grid", 0, Inf)),
    car_proper = list(
      tau_grid = .count_spatial_check_grid(dflt(tau_grid, "tau_car"), "tau_grid", 0, Inf),
      rho_grid = .count_spatial_check_grid(
        dflt(rho_grid, "rho_car"), "rho_grid", 0, 1, open = TRUE,
        hint = "Pass explicit eigenvalue bounds via spatial_car_proper().")),
    bym2 = list(
      sigma_grid = .count_spatial_check_grid(dflt(sigma_grid, "sigma_bym2"),
                                             "sigma_grid", 0, Inf),
      rho_grid   = .count_spatial_check_grid(dflt(rho_grid, "rho_bym2"),
                                             "rho_grid", 0, 1, open = FALSE)))
}

# Fit one areal count model. `family` is "nmix" or "removal", `kind` the field.
# Data arguments: long-form counts `y` with their 1-based `site_idx`, the
# site -> spatial unit map, the abundance design (one row per site) and the
# detection design (one row per visit / pass), and the CSR adjacency of the
# `n_spatial` units (`adj_col_idx` 0-based). `tau_grid` / `rho_grid` /
# `sigma_grid` override the kind's default axes, `r_grid` the NB size axis;
# `scale_factor` / `node_prec` the BYM2 scaling (computed from the graph when
# NULL); `beta_*_init`, `z_init` (icar / car_proper) and `v_init` / `w_init`
# (bym2) are warm starts. Returns a list of class `nmix_spatial_fit`: the C++
# grid result (theta_grid, log_marginal, modes, per-grid diagnostics) plus the
# grid-weighted summaries (weights, hyperparameter moments, coefficient means
# and the law-of-total-covariance vcov, the posterior-mean field, the NB size
# when integrated), `n_sites`, `n_obs`, `prior_type` and `call`.
.count_spatial_fit <- function(family = c("nmix", "removal"),
                               kind = c("icar", "car_proper", "bym2"),
                               y, site_idx, map_site_to_unit, X_lambda, X_p,
                               adj_row_ptr, adj_col_idx, n_neighbors, n_spatial,
                               tau_grid = NULL, rho_grid = NULL, sigma_grid = NULL,
                               mixture = c("P", "NB"), r_grid = NULL,
                               scale_factor = NULL, node_prec = NULL,
                               beta_lambda_init = NULL, beta_p_init = NULL,
                               z_init = NULL, v_init = NULL, w_init = NULL,
                               K_max = NULL, max_iter = 100L, tol = 1e-6,
                               verbose = FALSE) {
  family  <- match.arg(family)
  kind    <- match.arg(kind)
  mixture <- match.arg(mixture)
  y                <- as.integer(y)
  site_idx         <- as.integer(site_idx)
  map_site_to_unit <- as.integer(map_site_to_unit)
  bym2 <- identical(kind, "bym2")
  pp <- .count_spatial_prep(
    y, site_idx, X_lambda, X_p, mixture, beta_lambda_init, beta_p_init, K_max,
    r_grid,
    K_floor = if (identical(family, "removal")) .count_K_floor_site_total
              else .count_K_floor_max_y,
    map_site_to_unit = map_site_to_unit, n_spatial = n_spatial,
    adj_row_ptr = adj_row_ptr, n_neighbors = n_neighbors,
    latent = if (bym2) list(v_init = v_init, w_init = w_init) else list(z_init = z_init),
    n_latent = c(n_spatial = n_spatial))
  grids <- .count_spatial_grids(kind, tau_grid, rho_grid, sigma_grid)
  num_or_null <- function(x) if (is.null(x)) NULL else as.numeric(x)

  args <- list(
    y = y, site_idx = site_idx, map_site_to_unit_R = map_site_to_unit,
    X_lambda_R = X_lambda, X_p_R = X_p,
    adj_row_ptr = as.integer(adj_row_ptr), adj_col_idx = as.integer(adj_col_idx),
    n_neighbors = as.integer(n_neighbors), n_spatial = as.integer(n_spatial))
  if (!identical(kind, "car_proper")) {
    args$node_comp <- .tobs_node_comp_csr(adj_row_ptr, adj_col_idx, n_spatial)
  }
  if (bym2) {
    sc <- .bym2_resolve_scaling(scale_factor, node_prec, adj_row_ptr, adj_col_idx,
                                n_spatial)
    scale_factor <- sc$scale
    args$node_prec <- as.numeric(sc$node_prec)
  }
  args <- c(args, lapply(grids, as.numeric), list(r_grid = as.numeric(pp$r_grid)))
  if (bym2) args$scale_factor <- as.numeric(scale_factor)
  args <- c(args, list(
    beta_lambda_init = as.numeric(pp$beta_lambda_init),
    beta_p_init      = as.numeric(pp$beta_p_init)))
  args <- c(args, if (bym2) list(v_init = num_or_null(v_init), w_init = num_or_null(w_init))
                  else list(z_init = num_or_null(z_init)))
  args <- c(args, list(K_max = pp$K_max, max_iter = as.integer(max_iter),
                       tol = as.numeric(tol), verbose = isTRUE(verbose)))
  fit <- do.call(.cpp_nmix_progress, c(list(.count_spatial_cpp(family, kind)), args))

  packed <- if (bym2) {
    .count_spatial_pack_bym2_common(fit, pp$p_lam, pp$p_p, n_spatial, X_lambda,
                                    X_p, mixture, scale_factor)
  } else {
    common <- .count_spatial_pack_common(fit, pp$p_lam, pp$p_p, n_spatial,
                                         X_lambda, X_p, mixture)
    if (identical(kind, "car_proper")) {
      rho <- .tobs_weighted_moment(common$weights, fit$theta_grid[, "rho"])
      common <- c(common, list(rho_mean = unname(rho["mean"]),
                               rho_sd = unname(rho["sd"])))
    }
    common
  }
  out <- c(fit, packed,
           list(n_sites = pp$n_sites, n_obs = pp$n_obs, prior_type = kind,
                call = match.call()))
  .count_spatial_warn_boundary(out)
  class(out) <- c("nmix_spatial_fit", "list")
  out
}

#' @exportS3Method print nmix_spatial_fit
print.nmix_spatial_fit <- function(x, ...) {
  ptype <- x$prior_type %||% "icar"
  label <- switch(ptype,
                  icar       = "ICAR",
                  car_proper = "CAR(rho)",
                  bym2       = "BYM2",
                  toupper(ptype))
  mix <- x$mixture %||% "P"
  cat(sprintf("tulpa spatial N-mixture (%s, mixture = %s) nested-Laplace fit\n",
              label, mix))
  cat(sprintf("  n_sites = %d   n_obs = %d   n_spatial = %d   K.max = %d\n",
              x$n_sites, x$n_obs, x$n_spatial, x$K_max))
  if (identical(mix, "NB") && is.finite(x$r_mean %||% NA_real_)) {
    cat(sprintf("  NB size r: mean = %.3g   sd = %.3g (grid-integrated)\n",
                x$r_mean, x$r_sd))
  }
  if (ptype == "icar") {
    cat(sprintf("  tau grid = [%.3g, %.3g] over %d points\n",
                min(x$tau_grid), max(x$tau_grid), x$n_grid))
    cat(sprintf("  tau_mean = %.3g   tau_sd = %.3g\n",
                x$tau_mean, x$tau_sd))
  } else if (ptype == "bym2") {
    cat(sprintf("  sigma grid = [%.3g, %.3g] x %d   rho grid = [%.3g, %.3g] x %d\n",
                min(x$sigma_grid), max(x$sigma_grid), length(x$sigma_grid),
                min(x$rho_grid), max(x$rho_grid), length(x$rho_grid)))
    cat(sprintf("  scale_factor = %.4g\n", x$scale_factor))
    cat(sprintf("  sigma_mean = %.3g   sigma_sd = %.3g\n",
                x$sigma_mean, x$sigma_sd))
    cat(sprintf("  rho_mean   = %.3g   rho_sd   = %.3g\n",
                x$rho_mean, x$rho_sd))
  } else {
    cat(sprintf("  tau grid = [%.3g, %.3g] x %d   rho grid = [%.3g, %.3g] x %d\n",
                min(x$tau_grid), max(x$tau_grid), length(x$tau_grid),
                min(x$rho_grid), max(x$rho_grid), length(x$rho_grid)))
    cat(sprintf("  tau_mean = %.3g   tau_sd = %.3g\n",
                x$tau_mean, x$tau_sd))
    if (!is.null(x$rho_mean)) {
      cat(sprintf("  rho_mean = %.3g   rho_sd = %.3g\n",
                  x$rho_mean, x$rho_sd))
    }
  }
  cat("\nabundance (log lambda):\n")
  print(x$beta_lambda_mean)
  cat("\ndetection (logit p):\n")
  print(x$beta_p_mean)
  invisible(x)
}

# Resolve the NB-dispersion grid for the spatial fits. Poisson -> a single
# r = Inf node (the kernel takes its Poisson branch). NB -> a user grid or a
# default log-spaced grid of NB sizes spanning strong overdispersion to nearly
# Poisson. The outer nested-Laplace integration treats r like tau/rho/sigma.
.nmix_resolve_r_grid <- function(mixture, r_grid) {
  if (identical(mixture, "P")) return(Inf)
  if (is.null(r_grid)) {
    return(exp(seq(log(0.5), log(40), length.out = 6L)))
  }
  if (!is.numeric(r_grid) || length(r_grid) < 1L || any(r_grid <= 0) || anyNA(r_grid)) {
    stop("`r_grid` must be a vector of positive NB sizes.", call. = FALSE)
  }
  r_grid
}

# ---------------------------------------------------------------------------
# Shared preamble for the count-marginal spatial fitters
#
# The areal N-mixture wrappers (icar / car_proper / bym2), the SPDE mesh wrapper
# and the removal wrappers fit the same count marginal over the same designs,
# and so validate the same arguments. Everything below is the part they share;
# what differs is passed in, so a default or a check is a property of the axis
# or of the family rather than of which wrapper the caller happened to reach.
# ---------------------------------------------------------------------------

# Default outer-grid axes for the areal count fitters.
#
# The tau axis carries 9 nodes on the 1D ICAR grid and 7 on the 2D proper-CAR
# grid: the proper-CAR grid is the product of both axes, so 7 x 5 = 35 inner
# Newton solves already costs four times the ICAR grid's 9. Both count families
# read these, so the axis is what decides, not the family.
.count_spatial_default_grid <- function(axis) {
  switch(
    axis,
    tau_icar   = exp(seq(log(0.3), log(30), length.out = 9L)),
    tau_car    = exp(seq(log(0.3), log(30), length.out = 7L)),
    rho_car    = c(0.1, 0.3, 0.5, 0.75, 0.95),
    sigma_bym2 = exp(seq(log(0.2), log(3), length.out = 5L)),
    rho_bym2   = c(0.05, 0.3, 0.5, 0.7, 0.95),
    stop(sprintf("Unknown count-spatial grid axis '%s'.", axis), call. = FALSE)
  )
}

# Range check for a resolved outer-grid axis. `open` = the proper-CAR rho
# interval, which excludes both ends (rho = 1 is the intrinsic limit, where
# Q(rho) drops rank and the dense Cholesky behind log|Q(rho)| fails); the closed
# form is the BYM2 mixing weight, where both ends are meaningful (all-iid and
# all-structured). Applied at every door, so a grid valid for one areal count
# family is valid for the other.
.count_spatial_check_grid <- function(x, name, lo, hi, open = TRUE, hint = "") {
  if (!is.numeric(x) || !length(x) || anyNA(x)) {
    stop(sprintf("`%s` must be a non-empty numeric vector.", name), call. = FALSE)
  }
  bad <- if (open) any(x <= lo) || any(x >= hi) else any(x < lo) || any(x > hi)
  if (bad) {
    stop(sprintf("`%s` values must lie %s.%s", name,
                 if (open) sprintf("strictly in (%g, %g)", lo, hi)
                 else sprintf("in [%g, %g]", lo, hi),
                 if (nzchar(hint)) paste0(" ", hint) else ""),
         call. = FALSE)
  }
  as.numeric(x)
}

# Latent-N truncation floor: the largest count at any single site, which the
# N-mixture marginal must be able to represent.
.count_K_floor_max_y <- function(y, site_idx, n_sites) {
  list(value = max(y), label = "max(y)")
}

# Removal's floor: the depleting passes at a site sum, so the latent N clears
# that site's TOTAL removal, not its largest single pass.
.count_K_floor_site_total <- function(y, site_idx, n_sites) {
  tot <- tapply(as.integer(y),
                factor(as.integer(site_idx), levels = seq_len(n_sites)), sum)
  tot[is.na(tot)] <- 0L
  list(value = max(as.integer(tot)),
       label = "the largest per-site removal total")
}

# Resolve the latent-N truncation. An explicit K_max below the family's floor
# cannot represent the data; the default clears the floor by the same 100
# headroom at every door.
.count_spatial_K_max <- function(K_max, floor) {
  if (is.null(K_max)) return(as.integer(floor$value + 100L))
  K_max <- as.integer(K_max)
  if (K_max < floor$value) {
    stop(sprintf("K.max must be >= %s.", floor$label), call. = FALSE)
  }
  K_max
}

# Resolve and validate everything the count-marginal spatial fitters share:
# the two designs against the declared dimensions, the site-to-unit map and the
# CSR graph against `n_spatial`, the coefficient warm starts, the NB-size grid
# and the latent-N truncation, plus any latent-field warm start.
#
# `map_site_to_unit` / `adj_row_ptr` are NULL on the SPDE door (a mesh carries a
# projector, not a graph), and the checks they gate are skipped there.
# `K_floor` names the family's truncation rule; `latent` is the named list of
# field warm starts, each checked against `n_latent` -- whose NAME is what the
# message calls that length (`n_spatial` areal, `n_mesh` on the mesh).
.count_spatial_prep <- function(y, site_idx, X_lambda, X_p, mixture,
                                beta_lambda_init, beta_p_init, K_max, r_grid,
                                K_floor = .count_K_floor_max_y,
                                map_site_to_unit = NULL, n_spatial = NULL,
                                adj_row_ptr = NULL, n_neighbors = NULL,
                                latent = list(), n_latent = NULL) {
  if (!is.matrix(X_lambda)) stop("`X_lambda` must be a numeric matrix.", call. = FALSE)
  if (!is.matrix(X_p))      stop("`X_p` must be a numeric matrix.", call. = FALSE)
  n_sites <- nrow(X_lambda); n_obs <- nrow(X_p)
  p_lam   <- ncol(X_lambda); p_p   <- ncol(X_p)
  if (length(y) != n_obs) stop("length(y) must equal nrow(X_p).", call. = FALSE)
  if (length(site_idx) != n_obs) {
    stop("length(site_idx) must equal nrow(X_p).", call. = FALSE)
  }
  if (!is.null(map_site_to_unit)) {
    if (length(map_site_to_unit) != n_sites) {
      stop("length(map_site_to_unit) must equal nrow(X_lambda).", call. = FALSE)
    }
    if (any(map_site_to_unit < 1L) || any(map_site_to_unit > n_spatial)) {
      stop("map_site_to_unit values must lie in [1, n_spatial].", call. = FALSE)
    }
  }
  if (!is.null(adj_row_ptr)) {
    if (length(adj_row_ptr) != n_spatial + 1L) {
      stop("length(adj_row_ptr) must equal n_spatial + 1.", call. = FALSE)
    }
    if (length(n_neighbors) != n_spatial) {
      stop("length(n_neighbors) must equal n_spatial.", call. = FALSE)
    }
  }
  if (is.null(beta_lambda_init)) {
    beta_lambda_init <- c(log(mean(y) + 0.1), rep(0, p_lam - 1L))
  }
  if (is.null(beta_p_init)) beta_p_init <- rep(0, p_p)
  if (length(beta_lambda_init) != p_lam) {
    stop("length(beta_lambda_init) must equal ncol(X_lambda).", call. = FALSE)
  }
  if (length(beta_p_init) != p_p) {
    stop("length(beta_p_init) must equal ncol(X_p).", call. = FALSE)
  }
  for (nm in names(latent)) {
    v <- latent[[nm]]
    if (!is.null(v) && length(v) != n_latent) {
      stop(sprintf("length(%s) must equal %s.", nm, names(n_latent)),
           call. = FALSE)
    }
  }
  list(n_sites = n_sites, n_obs = n_obs, p_lam = p_lam, p_p = p_p,
       beta_lambda_init = beta_lambda_init, beta_p_init = beta_p_init,
       K_max = .count_spatial_K_max(K_max, K_floor(y, site_idx, n_sites)),
       r_grid = .nmix_resolve_r_grid(mixture, r_grid))
}

# The boundary diagnostic every count-marginal spatial fitter reports: posterior
# mass sitting on N = K_max means the truncation clipped the latent-N sum, so the
# fit is answering a question the window could not hold. Every reported summary
# is a grid-weighted mixture over the outer nodes, so the truncation error it
# carries is each node's boundary mass weighted by that node's posterior grid
# weight; a node the posterior gives no weight contributes nothing however far
# its field reaches into the latent-N tail. One wording, one threshold, one call.
.count_spatial_warn_boundary <- function(out) {
  b <- out$boundary_max
  w <- out$weights
  if (length(w) != length(b)) {
    stop("boundary check needs one grid weight per boundary_max entry.",
         call. = FALSE)
  }
  ok <- is.finite(b) & is.finite(w)
  contrib <- w[ok] * b[ok]
  mass <- sum(contrib)
  if (mass > 1e-4) {
    k <- which(ok)[which.max(contrib)]
    warning(sprintf(paste0("Posterior-weighted mass on N = K.max is %.2e ",
                           "(largest contribution: grid point %d, weight %.2e, ",
                           "boundary mass %.2e); raise K.max."),
                    mass, k, w[k], b[k]), call. = FALSE)
  }
  invisible(out)
}

# Shared grid-summarisation for the single-field areal count fits (ICAR /
# proper-CAR, abundance-arm field on one spatial unit per site). The N-mixture
# and removal nested-Laplace wrappers differ only in the C++ entry point and the
# prior-specific hyperparameter summaries (tau here; rho added by the proper-CAR
# wrapper); the weights, weighted-mean coefficients, weighted field, grid-
# integrated coefficient covariance, and NB-size summary are identical. Single
# source of truth. Returns the common block the wrappers splice into `out`.
.count_spatial_pack_common <- function(fit, p_lam, p_p, n_spatial,
                                       X_lambda, X_p, mixture) {
  weights <- .tobs_grid_weights(fit, "tau_grid / data")
  tau  <- .tobs_weighted_moment(weights, fit$theta_grid[, "tau"])
  disp <- .nmix_dispersion_summary(mixture, fit$theta_grid, weights)

  modes <- fit$modes
  beta_lambda_mean <- as.numeric(crossprod(weights, modes[, seq_len(p_lam), drop = FALSE]))
  beta_p_mean      <- as.numeric(crossprod(
    weights, modes[, p_lam + seq_len(p_p), drop = FALSE]))
  z_mean <- as.numeric(crossprod(
    weights, modes[, p_lam + p_p + seq_len(n_spatial), drop = FALSE]))

  nm_lam <- colnames(X_lambda); nm_p <- colnames(X_p)
  if (is.null(nm_lam)) nm_lam <- paste0("lam_", seq_len(p_lam))
  if (is.null(nm_p))   nm_p   <- paste0("p_", seq_len(p_p))
  names(beta_lambda_mean) <- nm_lam
  names(beta_p_mean)      <- nm_p

  list(
    weights          = weights,
    tau_mean         = unname(tau["mean"]),
    tau_sd           = unname(tau["sd"]),
    mixture          = mixture,
    r_mean           = disp$r_mean,
    r_sd             = disp$r_sd,
    beta_lambda_mean = beta_lambda_mean,
    beta_p_mean      = beta_p_mean,
    vcov             = .nmix_grid_vcov(fit$cov_blocks, modes, weights,
                                       p_lam, p_p, c(nm_lam, nm_p)),
    grid_mixture     = .tobs_grid_mixture(weights,
                                          modes[, seq_len(p_lam + p_p), drop = FALSE],
                                          fit$cov_blocks),
    z_mean           = z_mean
  )
}

# Shared grid-summarisation for the BYM2 areal count fits (two-field v / w on the
# abundance arm). Like .count_spatial_pack_common but with the BYM2-specific
# sigma / rho summaries, the per-component (v, w) weighted means, and the total
# offset phi = sigma (sqrt(rho/scale) v + sqrt(1-rho) w) integrated over the grid.
# Single source of truth for the N-mixture and removal BYM2 wrappers.
.count_spatial_pack_bym2_common <- function(fit, p_lam, p_p, n_spatial,
                                            X_lambda, X_p, mixture, scale_factor) {
  weights <- .tobs_grid_weights(fit)
  sigma_vec <- fit$theta_grid[, "sigma"]; rho_vec <- fit$theta_grid[, "rho"]
  sigma <- .tobs_weighted_moment(weights, sigma_vec)
  rho   <- .tobs_weighted_moment(weights, rho_vec)
  disp  <- .nmix_dispersion_summary(mixture, fit$theta_grid, weights)

  modes <- fit$modes
  beta_lambda_mean <- as.numeric(crossprod(weights, modes[, seq_len(p_lam), drop = FALSE]))
  beta_p_mean      <- as.numeric(crossprod(
    weights, modes[, p_lam + seq_len(p_p), drop = FALSE]))
  v_idx <- p_lam + p_p + seq_len(n_spatial)
  w_idx <- p_lam + p_p + n_spatial + seq_len(n_spatial)
  v_mean <- as.numeric(crossprod(weights, modes[, v_idx, drop = FALSE]))
  w_mean <- as.numeric(crossprod(weights, modes[, w_idx, drop = FALSE]))

  phi_grid <- matrix(0, nrow = nrow(modes), ncol = n_spatial)
  for (k in seq_len(nrow(modes))) {
    sg <- sigma_vec[k]; rg <- rho_vec[k]
    phi_grid[k, ] <- (sg * sqrt(rg / scale_factor)) * modes[k, v_idx] +
                     (sg * sqrt(1 - rg)) * modes[k, w_idx]
  }
  phi_mean <- as.numeric(crossprod(weights, phi_grid))

  nm_lam <- colnames(X_lambda); nm_p <- colnames(X_p)
  if (is.null(nm_lam)) nm_lam <- paste0("lam_", seq_len(p_lam))
  if (is.null(nm_p))   nm_p   <- paste0("p_", seq_len(p_p))
  names(beta_lambda_mean) <- nm_lam; names(beta_p_mean) <- nm_p

  list(
    weights = weights,
    sigma_mean = unname(sigma["mean"]), sigma_sd = unname(sigma["sd"]),
    rho_mean = unname(rho["mean"]), rho_sd = unname(rho["sd"]),
    mixture = mixture,
    r_mean = disp$r_mean, r_sd = disp$r_sd,
    beta_lambda_mean = beta_lambda_mean, beta_p_mean = beta_p_mean,
    vcov = .nmix_grid_vcov(fit$cov_blocks, modes, weights, p_lam, p_p,
                           c(nm_lam, nm_p)),
    grid_mixture = .tobs_grid_mixture(
      weights, modes[, seq_len(p_lam + p_p), drop = FALSE], fit$cov_blocks),
    v_mean = v_mean, w_mean = w_mean, phi_mean = phi_mean)
}

# Posterior mean / sd of the NB size r from the outer grid weights (the "r"
# column of theta_grid). Poisson -> NA. Weighted moments mirror tau/rho/sigma.
.nmix_dispersion_summary <- function(mixture, theta_grid, weights) {
  if (identical(mixture, "P") || !("r" %in% colnames(theta_grid))) {
    return(list(r_mean = NA_real_, r_sd = NA_real_))
  }
  r <- .tobs_weighted_moment(weights, theta_grid[, "r"])
  list(r_mean = unname(r["mean"]), r_sd = unname(r["sd"]))
}

# Grid-integrated covariance of beta = (beta_lambda, beta_p) over the outer
# hyperparameter grid. `cov_blocks[[k]]` is the within-grid Laplace covariance
# at the k-th mode (the beta block of the joint H^{-1} the C++ kernel returns);
# .tobs_grid_vcov() adds the between-grid mode spread. A grid with no usable
# weight yields an all-NA covariance.
.nmix_grid_vcov <- function(cov_blocks, modes, weights, p_lam, p_p, nm) {
  p_beta     <- p_lam + p_p
  beta_modes <- modes[, seq_len(p_beta), drop = FALSE]
  ok <- is.finite(weights) & weights > 0
  V  <- if (any(ok)) {
    w <- weights; w[!ok] <- 0
    .tobs_grid_vcov(beta_modes, w, cov_blocks)
  } else {
    matrix(NA_real_, p_beta, p_beta)
  }
  dimnames(V) <- list(nm, nm)
  V
}
