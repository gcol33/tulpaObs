# =============================================================================
# sim_cover_hurdle.R — simulator for the cover-hurdle family (Phase 1a)
#
# Generates vegetation-cover style data with a Bernoulli occurrence layer and
# a lognormal positive-cover layer, optionally with a shared exponential
# spatial field across both layers.
# =============================================================================


#' Simulate cover-hurdle data (lognormal positive part)
#'
#' Generates synthetic data matching the `cover(positive = "lognormal")`
#' generative model used in Phase 1a:
#'
#' \describe{
#'   \item{eta_occ}{= X %*% beta_occ + spatial_occ}
#'   \item{eta_pos}{= X %*% beta_pos + spatial_pos}
#'   \item{p}{= plogis(eta_occ)}
#'   \item{occur}{~ Bernoulli(p)}
#'   \item{log_cover}{~ Normal(eta_pos, sigma_pos^2) | occur = 1}
#'   \item{cover}{= exp(log_cover) when occur = 1, else 0}
#' }
#'
#' Both layers share the same design matrix by default (one continuous
#' covariate plus an intercept). When `spatial_range` is supplied, a simple
#' exponential-kernel Gaussian field on a random unit-square coordinate set
#' is added to both linear predictors; the two layers see independently
#' simulated draws of that field (matching the Phase 1a separate-fit story;
#' the shared-field model is Phase 1c).
#'
#' @param N Number of sites (default 200).
#' @param beta_occ Length-2 occurrence coefficients (intercept, slope on x).
#' @param beta_pos Length-2 log-cover coefficients (intercept, slope on x).
#' @param sigma_pos Lognormal residual standard deviation (default 0.4).
#' @param spatial_range Length scale of the exponential spatial field
#'   (in unit-square distance). `NULL` disables the spatial layer.
#' @param spatial_var Marginal variance of the spatial field (default 1).
#' @param seed Optional integer seed.
#' @return A list with:
#' \describe{
#'   \item{data}{A data frame with `cover`, covariate `x`, and `lon`, `lat`.}
#'   \item{y}{`data$cover` (length-N numeric vector, for passing to `tobs()`).}
#'   \item{coords}{An N x 2 matrix of coordinates.}
#'   \item{truth}{A list with `beta_occ`, `beta_pos`, `sigma_pos`,
#'     `p`, `mu`, `occur`, `spatial_occ`, `spatial_pos`.}
#' }
#' @export
#' @examples
#' sim <- simulate_cover(N = 200, seed = 1)
#' head(sim$data)
simulate_cover <- function(N             = 200L,
                                  beta_occ      = c(-0.5, 0.8),
                                  beta_pos      = c(-1.0, 0.3),
                                  sigma_pos     = 0.4,
                                  spatial_range = NULL,
                                  spatial_var   = 1,
                                  seed          = NULL) {
  if (!is.null(seed)) set.seed(seed)
  N <- as.integer(N)
  if (length(beta_occ) != 2L || length(beta_pos) != 2L) {
    stop("`beta_occ` and `beta_pos` must each be length-2 ",
         "(intercept + slope on x).", call. = FALSE)
  }

  x      <- stats::rnorm(N)
  coords <- cbind(lon = stats::runif(N), lat = stats::runif(N))
  X      <- cbind(`(Intercept)` = 1, x = x)

  eta_occ <- as.numeric(X %*% beta_occ)
  eta_pos <- as.numeric(X %*% beta_pos)

  spatial_occ <- rep(0, N)
  spatial_pos <- rep(0, N)
  if (!is.null(spatial_range)) {
    spatial_occ <- .draw_exp_field(coords, spatial_range, spatial_var)
    spatial_pos <- .draw_exp_field(coords, spatial_range, spatial_var)
    eta_occ <- eta_occ + spatial_occ
    eta_pos <- eta_pos + spatial_pos
  }

  p     <- stats::plogis(eta_occ)
  occur <- stats::rbinom(N, 1L, p)

  log_cover <- stats::rnorm(N, eta_pos, sigma_pos)
  cover     <- ifelse(occur == 1L, exp(log_cover), 0)
  cover     <- pmin(cover, 1)

  data <- data.frame(
    cover = cover,
    x     = x,
    lon   = coords[, 1L],
    lat   = coords[, 2L]
  )

  list(
    data   = data,
    y      = cover,
    coords = coords,
    truth  = list(
      beta_occ    = beta_occ,
      beta_pos    = beta_pos,
      sigma_pos   = sigma_pos,
      p           = p,
      mu          = exp(eta_pos + sigma_pos^2 / 2),
      occur       = occur,
      spatial_occ = spatial_occ,
      spatial_pos = spatial_pos
    )
  )
}


.draw_exp_field <- function(coords, range, var) {
  N <- nrow(coords)
  D <- as.matrix(stats::dist(coords))
  Sigma <- var * exp(-D / range)
  L <- tryCatch(t(chol(Sigma + diag(1e-8, N))),
                error = function(e) {
                  ev <- eigen(Sigma, symmetric = TRUE)
                  ev$vectors %*% diag(sqrt(pmax(ev$values, 0)))
                })
  as.numeric(L %*% stats::rnorm(N))
}
