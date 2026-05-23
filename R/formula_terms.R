# =============================================================================
# formula_terms.R — Structured terms for tobs() formulas
#
# tulpaObs specifies spatial / temporal / random-effect / latent structure
# *inside* the process formulas, the way lme4 (`(1|g)`), mgcv (`s()`), brms
# (`gp()`), and INLA (`f()`) do. Each special term is a constructor that
# returns a classed `tobs_*` spec. The parser (`.tobs_parse_formula`) detects
# these by name against the registry below and evaluates them with the model
# `data` columns in scope, so bare symbols (`gp(lon, lat)`, `re(observer)`)
# resolve to columns exactly like `s()` inside `gam()`.
#
# The term constructors are NOT exported. They are only meaningful inside a
# tobs() formula and are resolved through the registry, never called by the
# user directly. Process membership (which linear predictor an effect enters)
# is determined by *which process formula* the term appears in — there is no
# `shared` argument. A single realization is shared across processes with
# `copy("id")`.
# =============================================================================


# ---------------------------------------------------------------------------
# Internal: stamp class + id onto a spec list
# ---------------------------------------------------------------------------
.tobs_term <- function(spec, class, id = NULL, label = NULL) {
  spec$id    <- id
  spec$label <- label
  structure(spec, class = class)
}

# Collect a coords matrix from either an explicit `coords =` matrix or two or
# more bare coordinate columns passed positionally (resolved to vectors).
.tobs_collect_coords <- function(dots, coords, term) {
  if (!is.null(coords)) {
    if (!is.matrix(coords) || ncol(coords) != 2L) {
      stop(sprintf("%s(): `coords` must be a matrix with 2 columns.", term),
           call. = FALSE)
    }
    return(coords)
  }
  if (length(dots) < 2L) {
    stop(sprintf(
      "%s(): supply two coordinate columns, e.g. %s(lon, lat), or coords = M.",
      term, term), call. = FALSE)
  }
  m <- tryCatch(do.call(cbind, dots[1:2]),
                error = function(e) stop(sprintf(
                  "%s(): could not bind coordinate columns: %s",
                  term, conditionMessage(e)), call. = FALSE))
  if (ncol(m) != 2L) {
    stop(sprintf("%s(): coordinates must resolve to 2 columns.", term),
         call. = FALSE)
  }
  m
}

# Resolve a grouping/time vector to 1-based integer codes plus a level count.
.tobs_index_codes <- function(x, term, arg) {
  if (is.factor(x) || is.character(x)) {
    f <- as.factor(x)
    return(list(idx = as.integer(f), n = nlevels(f)))
  }
  if (is.numeric(x)) {
    f <- as.factor(x)
    return(list(idx = as.integer(f), n = nlevels(f)))
  }
  stop(sprintf("%s(): `%s` must be a factor, character, or numeric column.",
               term, arg), call. = FALSE)
}


# ---------------------------------------------------------------------------
# Spatial terms (areal)
# ---------------------------------------------------------------------------

# icar(graph)              — intrinsic CAR over an adjacency graph
.tobs_term_icar <- function(graph, id = NULL) {
  if (!is.matrix(graph)) stop("icar(): `graph` must be an adjacency matrix.", call. = FALSE)
  if (!isSymmetric(graph)) stop("icar(): `graph` must be symmetric.", call. = FALSE)
  csr <- adjacency_to_csr(graph)
  .tobs_term(list(
    type = "icar", n_units = nrow(graph),
    adj_row_ptr = csr$row_ptr, adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors
  ), class = "tobs_spatial", id = id, label = "icar")
}

# bym2(graph, scale_factor) — BYM2 reparameterization of ICAR + IID
.tobs_term_bym2 <- function(graph, scale_factor = NULL, id = NULL) {
  if (!is.matrix(graph)) stop("bym2(): `graph` must be an adjacency matrix.", call. = FALSE)
  if (!isSymmetric(graph)) stop("bym2(): `graph` must be symmetric.", call. = FALSE)
  csr <- adjacency_to_csr(graph)
  if (is.null(scale_factor)) scale_factor <- compute_bym2_scale(graph)
  .tobs_term(list(
    type = "bym2", n_units = nrow(graph),
    adj_row_ptr = csr$row_ptr, adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors, scale_factor = scale_factor
  ), class = "tobs_spatial", id = id, label = "bym2")
}


# ---------------------------------------------------------------------------
# Spatial terms (continuous, NNGP / SPDE)
# ---------------------------------------------------------------------------

# gp(lon, lat, ...)        — NNGP-approximated Gaussian process
.tobs_term_gp <- function(..., coords = NULL, cov = "exponential", nu = 1.5,
                          nn = 15, id = NULL,
                          sigma2_prior_U = 1.0, sigma2_prior_alpha = 0.01,
                          phi_prior_lower = 0.01, phi_prior_upper = 10.0) {
  coords <- .tobs_collect_coords(list(...), coords, "gp")
  cov <- match.arg(cov, c("exponential", "matern", "gaussian", "spherical"))
  n <- nrow(coords); nn <- min(nn, n - 1L)
  nngp <- compute_nngp_neighbors(coords, nn)
  .tobs_term(list(
    type = "gp", n_obs = n, nn = nn,
    coords = as.vector(t(coords)),
    nn_idx = as.vector(t(nngp$nn_idx)),
    nn_dist = as.vector(t(nngp$nn_dist)),
    nn_neighbor_dist = as.vector(nngp$nn_neighbor_dist),
    nn_order = nngp$nn_order, nn_order_inv = nngp$nn_order_inv,
    cov_type = cov, nu = nu,
    sigma2_prior_U = sigma2_prior_U, sigma2_prior_alpha = sigma2_prior_alpha,
    phi_prior_lower = phi_prior_lower, phi_prior_upper = phi_prior_upper
  ), class = "tobs_spatial", id = id, label = "gp")
}

# multiscale_gp(lon, lat, ...) — two-scale (local + regional) NNGP
.tobs_term_multiscale_gp <- function(..., coords = NULL, cov = "exponential",
                                     nu = 1.5, nn_local = 15, nn_regional = 15,
                                     id = NULL,
                                     range_local_lower = 0.01, range_local_upper = 10.0,
                                     range_regional_lower = 0.01, range_regional_upper = 100.0,
                                     sigma2_local_prior_U = 1.0, sigma2_local_prior_alpha = 0.01,
                                     sigma2_regional_prior_U = 1.0, sigma2_regional_prior_alpha = 0.01) {
  coords <- .tobs_collect_coords(list(...), coords, "multiscale_gp")
  cov <- match.arg(cov, c("exponential", "matern", "gaussian", "spherical"))
  n <- nrow(coords)
  nn_local <- min(nn_local, n - 1L); nn_regional <- min(nn_regional, n - 1L)
  loc <- compute_nngp_neighbors(coords, nn_local)
  reg <- compute_nngp_neighbors(coords, nn_regional)
  .tobs_term(list(
    type = "multiscale_gp", n_obs = n, coords = as.vector(t(coords)),
    nn_local = nn_local,
    nn_idx_local = as.vector(t(loc$nn_idx)),
    nn_dist_local = as.vector(t(loc$nn_dist)),
    nn_neighbor_dist_local = as.vector(loc$nn_neighbor_dist),
    nn_order_local = loc$nn_order, nn_order_inv_local = loc$nn_order_inv,
    nn_regional = nn_regional,
    nn_idx_regional = as.vector(t(reg$nn_idx)),
    nn_dist_regional = as.vector(t(reg$nn_dist)),
    nn_neighbor_dist_regional = as.vector(reg$nn_neighbor_dist),
    nn_order_regional = reg$nn_order, nn_order_inv_regional = reg$nn_order_inv,
    cov_type = cov, nu = nu,
    range_local_lower = range_local_lower, range_local_upper = range_local_upper,
    range_regional_lower = range_regional_lower, range_regional_upper = range_regional_upper,
    sigma2_local_prior_U = sigma2_local_prior_U, sigma2_local_prior_alpha = sigma2_local_prior_alpha,
    sigma2_regional_prior_U = sigma2_regional_prior_U, sigma2_regional_prior_alpha = sigma2_regional_prior_alpha
  ), class = "tobs_spatial", id = id, label = "multiscale_gp")
}

# spde(lon, lat, ...)      — continuous Matern field via triangular mesh
.tobs_term_spde <- function(..., coords = NULL, mesh = NULL, max_edge = NULL,
                            cutoff = 0, nu = 1, id = NULL,
                            prior_range = c(0.5, 0.5), prior_sigma = c(1, 0.5)) {
  if (is.null(coords) && is.null(mesh)) {
    coords <- .tobs_collect_coords(list(...), coords, "spde")
  }
  tulpa_spec <- tulpa::spatial_spde(
    coords = coords, data = NULL, mesh = mesh,
    max_edge = max_edge, cutoff = cutoff, nu = nu,
    prior_range = prior_range, prior_sigma = prior_sigma
  )
  .tobs_term(list(
    type = "spde", tulpa_spec = tulpa_spec, n_units = tulpa_spec$n_mesh,
    nu = nu, prior_range = prior_range, prior_sigma = prior_sigma
  ), class = "tobs_spatial", id = id, label = "spde")
}


# ---------------------------------------------------------------------------
# Random effects, temporal structure, SVC, latent factors
# ---------------------------------------------------------------------------

# re(group, ...)           — grouped random effect (intercept / slope / iid)
.tobs_term_re <- function(group, type = c("intercept", "slope", "iid"),
                          covariate = NULL, model = "iid",
                          correlated = TRUE, sigma_scale = 1, id = NULL) {
  type  <- match.arg(type)
  model <- match.arg(model, c("iid", "ar1", "rw1", "rw2"))
  if (type == "slope" && is.null(covariate)) {
    stop("re(): `covariate` must be given for a random slope.", call. = FALSE)
  }
  codes <- .tobs_index_codes(group, "re", "group")
  .tobs_term(list(
    group_idx = codes$idx, n_groups = codes$n,
    type = type, covariate = covariate, model = model,
    correlated = correlated, sigma_scale = sigma_scale
  ), class = "tobs_re", id = id, label = "re")
}

# temporal(time, ...)      — AR1 / RW1 / RW2 / IID temporal field
.tobs_term_temporal <- function(time, type = c("ar1", "rw1", "rw2", "iid"),
                                group = NULL, cyclic = FALSE,
                                tau_shape = 1, tau_rate = 0.01, id = NULL) {
  type <- match.arg(type)
  tcodes <- .tobs_index_codes(time, "temporal", "time")
  spec <- list(
    type = type, time_idx = tcodes$idx, n_times = tcodes$n,
    cyclic = cyclic, tau_shape = tau_shape, tau_rate = tau_rate
  )
  if (!is.null(group)) {
    gcodes <- .tobs_index_codes(group, "temporal", "group")
    spec$group_idx <- gcodes$idx
    spec$n_groups  <- gcodes$n
  }
  .tobs_term(spec, class = "tobs_temporal", id = id, label = "temporal")
}

# svc(lon, lat, indices)   — spatially varying coefficients on design columns
.tobs_term_svc <- function(..., indices, coords = NULL, cov = "exponential",
                           nn = 15, id = NULL, sigma2_prior_scale = 1,
                           phi_prior_lower = 0.01, phi_prior_upper = 10) {
  coords <- .tobs_collect_coords(list(...), coords, "svc")
  cov <- match.arg(cov, c("exponential", "matern", "gaussian"))
  n <- nrow(coords); nn <- min(nn, n - 1L)
  nngp <- compute_nngp_neighbors(coords, nn)
  .tobs_term(list(
    indices = as.integer(indices), n_svc = length(indices),
    n_obs = n, nn = nn, coords = as.vector(t(coords)),
    nn_idx = as.vector(t(nngp$nn_idx)), nn_dist = as.vector(t(nngp$nn_dist)),
    nn_order = nngp$nn_order, nn_order_inv = nngp$nn_order_inv,
    cov_type = cov, sigma2_prior_scale = sigma2_prior_scale,
    phi_prior_lower = phi_prior_lower, phi_prior_upper = phi_prior_upper
  ), class = "tobs_svc", id = id, label = "svc")
}

# latent(n_factors, ...)   — latent factors for community models
.tobs_term_latent <- function(n_factors, constraint = 0,
                              sigma_prior_rate = 1, id = NULL) {
  .tobs_term(list(
    n_factors = as.integer(n_factors), constraint = as.integer(constraint),
    sigma_prior_rate = sigma_prior_rate
  ), class = "tobs_latent", id = id, label = "latent")
}

# copy(id)                 — share one realization of a named term across
#                            another process's linear predictor
.tobs_term_copy <- function(id, scale = NULL) {
  if (!is.character(id) || length(id) != 1L) {
    stop("copy(): `id` must be a single string naming a term's `id =`.",
         call. = FALSE)
  }
  .tobs_term(list(ref = id, scale = scale),
             class = "tobs_copy", id = id, label = "copy")
}


# ---------------------------------------------------------------------------
# Registry — name -> constructor. Adding a component is one entry here.
# ---------------------------------------------------------------------------
.tobs_terms <- list(
  icar          = .tobs_term_icar,
  bym2          = .tobs_term_bym2,
  gp            = .tobs_term_gp,
  multiscale_gp = .tobs_term_multiscale_gp,
  spde          = .tobs_term_spde,
  re            = .tobs_term_re,
  temporal      = .tobs_term_temporal,
  svc           = .tobs_term_svc,
  latent        = .tobs_term_latent,
  copy          = .tobs_term_copy
)

# Names of the registered special terms (used by the parser to detect them).
.tobs_term_names <- function() names(.tobs_terms)


# ---------------------------------------------------------------------------
# Print methods for the specs (used when echoing a parsed model)
# ---------------------------------------------------------------------------

#' @export
print.tobs_spatial <- function(x, ...) {
  n <- if (!is.null(x$n_units)) x$n_units else x$n_obs
  cat(sprintf("tobs spatial term: %s (%d units)%s\n", x$type, n,
              if (!is.null(x$id)) sprintf(" [id: %s]", x$id) else ""))
  if (identical(x$type, "spde")) {
    cat(sprintf("  Matern nu=%g, mesh=%d nodes\n", x$nu, x$n_units))
  }
  invisible(x)
}

#' @export
print.tobs_re <- function(x, ...) {
  cat(sprintf("tobs re term: %s (%s, %d groups)%s\n", x$type, x$model,
              x$n_groups, if (!is.null(x$id)) sprintf(" [id: %s]", x$id) else ""))
  if (!is.null(x$covariate)) cat(sprintf("  Covariate: %s\n", x$covariate))
  invisible(x)
}

#' @export
print.tobs_temporal <- function(x, ...) {
  cat(sprintf("tobs temporal term: %s (%d times)%s\n", x$type, x$n_times,
              if (!is.null(x$id)) sprintf(" [id: %s]", x$id) else ""))
  if (isTRUE(x$cyclic)) cat("  Cyclic: yes\n")
  invisible(x)
}

#' @export
print.tobs_copy <- function(x, ...) {
  cat(sprintf("tobs copy term: -> %s%s\n", x$ref,
              if (!is.null(x$scale)) sprintf(" (scaled)") else ""))
  invisible(x)
}
