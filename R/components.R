#' Specify random effects for tobs models
#'
#' @param group Character name of grouping variable in data, or integer vector
#'   of group assignments (1-based).
#' @param type One of `"intercept"` (default), `"slope"`, or `"iid"`.
#' @param covariate For random slopes, the covariate name.
#' @param model Temporal structure: `"iid"` (default), `"ar1"`, `"rw1"`, `"rw2"`.
#' @param shared Logical vector: which processes get this RE.
#'   Default `c(TRUE, FALSE)` = occupancy only.
#' @param correlated Logical; for slopes, estimate correlations? Default TRUE.
#' @param sigma_scale Prior scale for RE standard deviation (default 1).
#'
#' @return A `tobs_re` object
#' @export
tobs_re <- function(group, type = c("intercept", "slope", "iid"),
                    covariate = NULL, model = "iid",
                    shared = c(TRUE, FALSE), correlated = TRUE,
                    sigma_scale = 1) {
  type <- match.arg(type)
  model <- match.arg(model, c("iid", "ar1", "rw1", "rw2"))

  if (type == "slope" && is.null(covariate)) {
    stop("covariate must be specified for random slopes")
  }

  structure(list(
    group = group,
    type = type,
    covariate = covariate,
    model = model,
    shared = shared,
    correlated = correlated,
    sigma_scale = sigma_scale
  ), class = "tobs_re")
}

#' Specify temporal structure for tobs models
#'
#' @param type One of `"ar1"`, `"rw1"`, `"rw2"`, `"iid"`.
#' @param time Character name of time variable in data, or integer vector
#'   of time indices (1-based).
#' @param group Optional character name of panel/group variable.
#' @param shared Logical vector: which processes get the temporal effect.
#'   Default `c(TRUE, FALSE)` = occupancy only.
#' @param cyclic Logical; cyclic RW1 for seasonal patterns? Default FALSE.
#' @param tau_shape Shape parameter for Gamma prior on temporal precision.
#' @param tau_rate Rate parameter for Gamma prior on temporal precision.
#'
#' @return A `tobs_temporal` object
#' @export
tobs_temporal <- function(type = c("ar1", "rw1", "rw2", "iid"),
                          time, group = NULL,
                          shared = c(TRUE, FALSE),
                          cyclic = FALSE,
                          tau_shape = 1, tau_rate = 0.01) {
  type <- match.arg(type)

  structure(list(
    type = type,
    time = time,
    group = group,
    shared = shared,
    cyclic = cyclic,
    tau_shape = tau_shape,
    tau_rate = tau_rate
  ), class = "tobs_temporal")
}

#' Specify spatially-varying coefficients
#'
#' @param indices Integer vector of design matrix column indices to vary
#'   spatially. 1 = intercept.
#' @param coords Matrix of coordinates (n_sites x 2).
#' @param cov Covariance function: `"exponential"` (default), `"matern"`, `"gaussian"`.
#' @param nn Number of nearest neighbors for NNGP (default 15).
#' @param shared Logical vector: which processes get SVC.
#'   Default `c(TRUE, FALSE)` = occupancy only.
#' @param sigma2_prior_scale Prior scale for SVC variance (default 1).
#' @param phi_prior_lower Lower bound for range prior (default 0.01).
#' @param phi_prior_upper Upper bound for range prior (default 10).
#'
#' @return A `tobs_svc` object
#' @export
tobs_svc <- function(indices, coords, cov = "exponential", nn = 15,
                     shared = c(TRUE, FALSE),
                     sigma2_prior_scale = 1,
                     phi_prior_lower = 0.01, phi_prior_upper = 10) {
  if (!is.matrix(coords) || ncol(coords) != 2) {
    stop("coords must be a matrix with 2 columns")
  }
  cov <- match.arg(cov, c("exponential", "matern", "gaussian"))
  n <- nrow(coords)
  nn <- min(nn, n - 1)

  nngp <- compute_nngp_neighbors(coords, nn)

  structure(list(
    indices = as.integer(indices),
    n_svc = length(indices),
    n_obs = n,
    nn = nn,
    coords = as.vector(t(coords)),
    nn_idx = as.vector(t(nngp$nn_idx)),
    nn_dist = as.vector(t(nngp$nn_dist)),
    nn_order = nngp$nn_order,
    nn_order_inv = nngp$nn_order_inv,
    cov_type = cov,
    shared = shared,
    sigma2_prior_scale = sigma2_prior_scale,
    phi_prior_lower = phi_prior_lower,
    phi_prior_upper = phi_prior_upper
  ), class = "tobs_svc")
}

#' Specify latent factors for community models
#'
#' @param n_factors Number of latent factors.
#' @param shared Logical; enter both processes? Default TRUE.
#' @param constraint `0` = sum-to-zero (default), `1` = first-zero.
#' @param sigma_prior_rate Rate for Gamma prior on factor precision (default 1).
#'
#' @return A `tobs_latent` object
#' @export
tobs_latent <- function(n_factors, shared = TRUE, constraint = 0,
                        sigma_prior_rate = 1) {
  structure(list(
    n_factors = as.integer(n_factors),
    shared = shared,
    constraint = as.integer(constraint),
    sigma_prior_rate = sigma_prior_rate
  ), class = "tobs_latent")
}

#' @export
print.tobs_re <- function(x, ...) {
  cat(sprintf("tobs RE: %s (%s model)\n", x$type, x$model))
  cat(sprintf("  Group: %s\n", if (is.character(x$group)) x$group else "custom"))
  if (!is.null(x$covariate)) cat(sprintf("  Covariate: %s\n", x$covariate))
  invisible(x)
}

#' @export
print.tobs_temporal <- function(x, ...) {
  cat(sprintf("tobs temporal: %s\n", x$type))
  cat(sprintf("  Time: %s\n", if (is.character(x$time)) x$time else "custom"))
  if (x$cyclic) cat("  Cyclic: yes\n")
  invisible(x)
}


#' Specify community-level random effects
#'
#' Adds a species-level random effect for community occupancy models. A
#' convenience wrapper around [tobs_re()] that sets the group to the species
#' identifier.
#'
#' @param type `"intercept"` (default) or `"slope"`.
#' @param covariate For random slopes, the covariate name.
#' @param shared Logical vector: which processes get this RE.
#'   Default `c(TRUE, TRUE)` = both occupancy and detection.
#' @param sigma_scale Prior scale for RE standard deviation (default 1).
#'
#' @return A `tobs_re` object with group set to `"species"`.
#' @export
tobs_community_re <- function(type = c("intercept", "slope"),
                              covariate = NULL,
                              shared = c(TRUE, TRUE),
                              sigma_scale = 1) {
  type <- match.arg(type)
  tobs_re(group = "species", type = type, covariate = covariate,
          shared = shared, sigma_scale = sigma_scale)
}


#' Specify areal spatial structure
#'
#' Convenience wrapper for BYM2 or ICAR spatial models using an adjacency
#' matrix. Equivalent to [tobs_bym2()] or [tobs_icar()].
#'
#' @param adj Adjacency matrix (n_sites x n_sites, symmetric, 0/1).
#' @param model `"bym2"` (default) or `"icar"`.
#' @param ... Additional arguments passed to [tobs_bym2()] or [tobs_icar()].
#'
#' @return A `tobs_spatial` object.
#' @export
tobs_areal <- function(adj, model = c("bym2", "icar"), ...) {
  model <- match.arg(model)
  if (model == "bym2") {
    tobs_bym2(adjacency = adj, ...)
  } else {
    tobs_icar(adjacency = adj, ...)
  }
}
