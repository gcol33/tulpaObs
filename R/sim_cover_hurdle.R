# =============================================================================
# sim_cover_hurdle.R — simulator for the cover-hurdle family (Phase 1a)
#
# Generates vegetation-cover style data with a Bernoulli occurrence layer and
# a lognormal positive-cover layer, optionally with a shared exponential
# spatial field across both layers.
# =============================================================================


#' Simulate cover-hurdle data (lognormal positive part)
#'
#' Generates synthetic data matching the `cover(response = "lognormal")`
#' generative model used in Phase 1a:
#'
#' \describe{
#'   \item{eta_occ}{`= X %*% beta_occ + spatial_occ`}
#'   \item{eta_pos}{`= X %*% beta_pos + spatial_pos`}
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
                                  response      = c("lognormal", "gaussian"),
                                  spatial_range = NULL,
                                  spatial_var   = 1,
                                  seed          = NULL) {
  response <- match.arg(response)
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

  if (identical(response, "gaussian")) {
    # Delta-normal hurdle (gcol33/tulpaObs#112): the positive magnitude is a plain
    # Gaussian on the raw response (no log, no [0, 1] clamp); absence is the 0
    # sentinel. mu on the response scale is eta_pos.
    mag   <- stats::rnorm(N, eta_pos, sigma_pos)
    cover <- ifelse(occur == 1L, mag, 0)
    mu_truth <- eta_pos
  } else {
    log_cover <- stats::rnorm(N, eta_pos, sigma_pos)
    cover     <- ifelse(occur == 1L, exp(log_cover), 0)
    cover     <- pmin(cover, 1)
    mu_truth  <- exp(eta_pos + sigma_pos^2 / 2)
  }

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
      response    = response,
      p           = p,
      mu          = mu_truth,
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


#' Simulate cover-hurdle data with a shared BYM2 spatial field
#'
#' Generates synthetic data matching the joint nested-Laplace cover-hurdle
#' parameterisation (`sigma_occ * z` on the occurrence arm, `alpha *
#' sigma * z` on the cover arm) where `z` is a BYM2 latent field on the
#' supplied region adjacency. The cover-arm field amplitude
#' `sigma_pos = alpha * sigma` is the derived ratio surfaced as the
#' `"alpha"` column of `fit$joint$theta_grid`.
#'
#' Both BYM2 sub-blocks (`phi`, `theta`) are drawn IID Normal **and
#' demeaned to mean zero** before scaling, matching the sum-to-zero
#' constraint applied to both sub-blocks by `tulpa::spatial_bym2()` inside
#' the nested-Laplace engine (see `.joint_inner_var()`). Without the
#' demean each seed carries `mean(w_s) ~ N(0, sigma^2 / n_s)` and the
#' constrained intercept identified by the engine targets
#' `beta_pos_0_truth + alpha * mean(w_s)` rather than `beta_pos_0_truth`;
#' coverage of the *population* truth then collapses with alpha while
#' the engine-reported posterior is correctly calibrated for the
#' constrained parameter. See `example/validation/SUMMARY.md` in
#' INLAabun (Demo 3) for the diagnosis.
#'
#' @param N Number of sites (default 300).
#' @param adj `n_s x n_s` integer adjacency matrix for the BYM2 field.
#' @param beta_occ Length-2 occurrence coefficients (intercept + slope).
#' @param beta_pos Length-2 cover-arm coefficients (intercept + slope).
#' @param sigma Marginal spatial-field amplitude (default 0.6); the
#'   occurrence-arm linear predictor adds `sigma * z[region]`.
#' @param rho BYM2 mixing parameter, in `[0, 1]` (default 0.7).
#'   `rho = 1` is pure ICAR; `rho = 0` is pure IID.
#' @param alpha Cover-arm scaling: cover-arm linear predictor adds
#'   `alpha * sigma * z[region]` (so `sigma_pos = alpha * sigma`).
#' @param response Likelihood for the positive arm: `"beta"` or
#'   `"lognormal"`.
#' @param phi Beta precision when `positive = "beta"` (default 30).
#' @param sigma_pos_resid Lognormal residual SD when `positive =
#'   "lognormal"` (default 0.4); independent of the spatial `sigma`.
#' @param seed Optional integer seed.
#' @return A list with:
#' \describe{
#'   \item{data}{Data frame with `x` and `region` (factor).}
#'   \item{y}{Length-N numeric vector of cover values (zeros where
#'     `occur = 0`).}
#'   \item{adj}{The adjacency matrix passed in (for downstream
#'     `tulpa::spatial_bym2()`).}
#'   \item{truth}{Named list with `beta_occ`, `beta_pos`, `sigma`,
#'     `rho`, `alpha`, `sigma_pos = alpha * sigma`, `positive`, plus the
#'     simulated `phi_f`, `theta_f`, `w_s`, `region`.}
#' }
#' @export
#' @examples
#' adj <- matrix(0L, 10, 10)
#' for (s in 1:10) for (j in setdiff(c(s - 1L, s + 1L), c(0L, 11L)))
#'   adj[s, j] <- 1L
#' sim <- simulate_cover_joint(N = 100, adj = adj, alpha = 1.0, seed = 1)
#' head(sim$data)
simulate_cover_joint <- function(N               = 300L,
                                 adj,
                                 beta_occ        = c(-0.3, 0.7),
                                 beta_pos        = c( 0.4, -0.5),
                                 sigma           = 0.6,
                                 rho             = 0.7,
                                 alpha           = 1.0,
                                 positive        = c("beta", "lognormal",
                                                     "gaussian"),
                                 phi             = 30,
                                 sigma_pos_resid = 0.4,
                                 seed            = NULL) {
  positive <- match.arg(positive)
  if (!is.null(seed)) set.seed(seed)
  N <- as.integer(N)
  if (!is.matrix(adj) || nrow(adj) != ncol(adj)) {
    stop("`adj` must be a square integer adjacency matrix.", call. = FALSE)
  }
  n_s <- nrow(adj)
  if (length(beta_occ) != 2L || length(beta_pos) != 2L) {
    stop("`beta_occ` and `beta_pos` must each be length-2 ",
         "(intercept + slope on x).", call. = FALSE)
  }
  if (!is.numeric(rho) || rho < 0 || rho > 1) {
    stop("`rho` must be in [0, 1].", call. = FALSE)
  }

  region <- sample.int(n_s, N, replace = TRUE)

  # Demean each BYM2 sub-block before scaling — matches the hard sum-to-zero
  # constraint applied to both `phi` and `theta` by the joint nested-Laplace
  # engine. See `.joint_inner_var()` for the constraint correction and the
  # function docstring above for the diagnosis.
  phi_f   <- stats::rnorm(n_s, 0, 1); phi_f   <- phi_f   - mean(phi_f)
  theta_f <- stats::rnorm(n_s, 0, 1); theta_f <- theta_f - mean(theta_f)
  w_s     <- sigma * (sqrt(rho) * phi_f + sqrt(1 - rho) * theta_f)

  x       <- stats::rnorm(N)
  eta_occ <- beta_occ[1] + beta_occ[2] * x + w_s[region]
  occur   <- stats::rbinom(N, 1L, stats::plogis(eta_occ))
  eta_pos <- beta_pos[1] + beta_pos[2] * x + alpha * w_s[region]
  is_pos  <- occur == 1L

  y <- numeric(N)
  if (positive == "beta") {
    mu_pos <- stats::plogis(eta_pos)
    y[is_pos] <- stats::rbeta(sum(is_pos),
                              mu_pos[is_pos] * phi,
                              (1 - mu_pos[is_pos]) * phi)
    y <- pmin(pmax(y, 0), 1 - 1e-6)
  } else if (positive == "gaussian") {
    # Identity-Gaussian arm (gcol33/tulpaObs#112): raw magnitude, no log, no clamp.
    mag       <- stats::rnorm(N, eta_pos, sigma_pos_resid)
    y[is_pos] <- mag[is_pos]
  } else {
    log_y     <- stats::rnorm(N, eta_pos, sigma_pos_resid)
    y[is_pos] <- exp(log_y[is_pos])
    y         <- pmin(y, 1 - 1e-6)
  }

  list(
    data  = data.frame(x = x, region = factor(region)),
    y     = y,
    adj   = adj,
    truth = list(
      beta_occ        = beta_occ,
      beta_pos        = beta_pos,
      sigma           = sigma,
      rho             = rho,
      alpha           = alpha,
      sigma_pos       = alpha * sigma,
      positive        = positive,
      phi             = if (positive == "beta") phi else NA_real_,
      sigma_pos_resid = if (positive %in% c("lognormal", "gaussian"))
                        sigma_pos_resid else NA_real_,
      phi_f           = phi_f,
      theta_f         = theta_f,
      w_s             = w_s,
      region          = region
    )
  )
}
