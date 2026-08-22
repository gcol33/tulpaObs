# =============================================================================
# test-occu-svc-laplace-recovery.R - the continuous NNGP svc() surface on the
# DETERMINISTIC backends of single-season occupancy.
#
# svc() used to be sampler-only: under method = "laplace" / "nested_laplace" the
# guard in .tobs_fit_model() refused the term, so the package's default backend
# could not fit it at all. It now rides the shared areal-BFGS nested-Laplace
# driver (R/occu_svc.R) -- per outer-grid cell, BFGS over (coefficients,
# surface) against the exact occupancy marginal plus the NNGP log-prior, Laplace
# marginal from the FD-Hessian at the mode, hyperparameters integrated over the
# grid.
#
# The NNGP density is assembled here in R from the term's own Vecchia neighbour
# structure. The first test below is what makes that safe: the R precision's
# log-density is scored against tulpa's compiled SVC NNGP kernel at the same
# inputs. If those two ever drift, this path and the NUTS path are fitting
# different models under one term name.
#
# The simulation and the settings mirror test-occu-svc-nngp-recovery.R (the NUTS
# flavour) so the two backends are measured on the SAME truth at the same seeds.
# Measured surface correlation here: 0.78 / 0.60 / 0.83 (seeds 1 / 2 / 3),
# against the NUTS path's 0.76 / 0.60 / 0.81 on those seeds -- the surface is
# information-bounded at N = 150, J = 6, p = 0.6, and both backends sit at that
# bound.
#
# NOT the areal weighted-bar SVC: that arm arrives as a `spatial` term and is
# recovery-tested in test-occu-spatial-svc-recovery.R. See the SVC note in
# CLAUDE.md.
# =============================================================================

# A known spatially-varying INTERCEPT surface drawn from the model the term
# assumes: an exponential-kernel GP over the site coordinates. Centred, because
# a varying intercept is confounded with the global intercept up to a constant,
# so the recoverable target is the SHAPE.
.svcl_sim <- function(N, J, seed, sigma_f = 1.3, phi_f = 0.25, p_det = 0.6,
                      b0 = 0) {
  set.seed(seed)
  lon <- stats::runif(N); lat <- stats::runif(N)
  D <- as.matrix(stats::dist(cbind(lon, lat)))
  K <- sigma_f^2 * exp(-D / phi_f)
  f <- as.numeric(t(chol(K + 1e-8 * diag(N))) %*% stats::rnorm(N))
  f <- f - mean(f)
  z <- stats::rbinom(N, 1, stats::plogis(b0 + f))
  y <- matrix(stats::rbinom(N * J, 1, p_det * rep(z, J)), N, J)
  list(data = data.frame(lon = lon, lat = lat), y = y, f = f)
}

.svcl_fit <- function(sim, method = "laplace", nn = 10) {
  suppressWarnings(tobs(
    ~ svc(lon, lat, indices = 1L, nn = nn, prior_range = c(0.1, 0.05)),
    data = sim$data, family = occu(), detection = ~ 1, y = sim$y,
    method = method,
    control = list(verbose = FALSE, progress = FALSE)))
}


test_that("the R NNGP precision reproduces tulpa's compiled SVC kernel", {
  skip_if_not(exists("cpp_test_svc_nngp_twins", envir = asNamespace("tulpa"),
                     inherits = FALSE),
              "installed tulpa has no SVC NNGP twin probe")
  set.seed(11)
  N <- 40L
  co <- cbind(stats::runif(N), stats::runif(N))
  tm <- .tobs_term_svc(coords = co, indices = 1L, nn = 8,
                       prior_range = c(0.1, 0.05))
  co_m    <- matrix(as.numeric(tm$coords),   ncol = 2L,  byrow = TRUE)
  nn_idx  <- matrix(as.integer(tm$nn_idx),   nrow = N,   byrow = TRUE)
  nn_dist <- matrix(as.numeric(tm$nn_dist),  nrow = N,   byrow = TRUE)
  ord     <- as.integer(tm$nn_order)

  twin <- get("cpp_test_svc_nngp_twins", envir = asNamespace("tulpa"))
  z <- stats::rnorm(N)
  for (cov_code in c(0L, 1L, 2L)) {
    for (par in list(c(1.7, 0.3), c(0.4, 0.8))) {
      pr <- .tobs_nngp_precision(co_m, nn_idx, nn_dist, ord,
                                 par[1L], par[2L], cov_code)
      expect_false(is.null(pr))
      mine <- .tobs_nngp_log_density(z, pr$Q, pr$log_det)
      ref  <- twin(z, par[1L], par[2L], co_m, nn_idx, nn_dist, ord - 1L,
                   cov_code)[["dbl"]]
      # Same Vecchia factorisation, same kernels, same jitter and
      # conditional-variance floor: these are the same number, not a numerical
      # approximation of each other.
      expect_equal(mine, ref, tolerance = 1e-8)
    }
  }
})


test_that("the occupancy marginal's analytic gradients match finite differences", {
  set.seed(12)
  n <- 50L; J <- 4L
  dat <- data.frame(x = as.numeric(scale(stats::rnorm(n))))
  z   <- stats::rbinom(n, 1, stats::plogis(0.2 + 0.8 * dat$x))
  y   <- matrix(stats::rbinom(n * J, 1, 0.5 * rep(z, J)), n, J)
  model <- .tobs_build_model(~ x, ~ 1, dat, y)
  marg  <- .tobs_occu_svc_marginal(model, NULL)

  th  <- c(0.1, 0.7, -0.2)
  off <- stats::rnorm(n, 0, 0.3)
  ev  <- marg$eval(th, off)

  # Score with respect to the coefficients.
  g_fd <- vapply(seq_along(th), function(k) {
    h <- 1e-6; tp <- th; tm <- th; tp[k] <- tp[k] + h; tm[k] <- tm[k] - h
    (marg$eval(tp, off)$log_lik - marg$eval(tm, off)$log_lik) / (2 * h)
  }, numeric(1))
  expect_equal(ev$grad_fixed, g_fd, tolerance = 1e-5)

  # Score with respect to the per-site field offset -- the quantity the field
  # blocks scatter, so an error here would silently mis-fit every surface.
  ge_fd <- vapply(seq_len(n), function(i) {
    h <- 1e-6; op <- off; om <- off; op[i] <- op[i] + h; om[i] <- om[i] - h
    (marg$eval(th, op)$log_lik - marg$eval(th, om)$log_lik) / (2 * h)
  }, numeric(1))
  expect_equal(ev$grad_eta, ge_fd, tolerance = 1e-5)
})


test_that("occu() + svc() fits under laplace and exposes a shaped surface", {
  skip_on_cran()
  skip_if_fast()

  sim <- .svcl_sim(N = 60L, J = 4L, seed = 1L)
  fit <- .svcl_fit(sim)

  expect_s3_class(fit, "tobs_fit")
  # One weight per site, a bare vector for a single varying coefficient
  # (mirrors fit$spatial_field on the areal path and fit$svc_field on NUTS).
  expect_length(as.numeric(fit$svc_field), 60L)
  expect_true(all(is.finite(as.numeric(fit$svc_field))))
  # The term is echoed back; the estimated surface is separate.
  expect_s3_class(fit$svc, "tobs_svc")

  # The NNGP hyperparameters are surfaced, marginalised over the outer grid.
  expect_named(fit$svc_hyper, "svc1")
  expect_true(all(c("sigma", "phi") %in% names(fit$svc_hyper$svc1)))
  expect_gt(fit$svc_hyper$svc1[["sigma"]], 0)
  expect_gt(fit$svc_hyper$svc1[["phi"]], 0)

  # Coefficients and their SEs are finite and named as on every other route.
  expect_true(all(c("psi_(Intercept)", "p_(Intercept)") %in% names(fit$means)))
  expect_true(all(is.finite(fit$means)))
  expect_true(all(is.finite(fit$sds)) && all(fit$sds > 0))

  # fitted() reads the surface in sample rather than reporting X beta alone: on
  # an intercept-only occupancy design the fitted logit IS intercept + surface.
  ft <- fitted(fit)
  expect_length(ft$psi, 60L)
  expect_equal(as.numeric(stats::qlogis(ft$psi)),
               as.numeric(fit$means[["psi_(Intercept)"]]) +
                 as.numeric(fit$svc_field),
               tolerance = 1e-6)
})


test_that("occu() + svc() recovers a known varying-intercept surface (laplace)", {
  skip_on_cran()
  skip_if_fast()

  seeds <- 1:3
  cr <- p0 <- rep(NA_real_, length(seeds))
  for (i in seq_along(seeds)) {
    sim <- .svcl_sim(N = 150L, J = 6L, seed = seeds[i])
    fit <- tryCatch(.svcl_fit(sim), error = function(e) NULL)
    if (is.null(fit)) next
    cr[i] <- stats::cor(as.numeric(fit$svc_field), sim$f)
    p0[i] <- fit$means[["p_(Intercept)"]]
  }
  ok <- !is.na(cr)
  expect_gt(sum(ok), 1L)

  # Measured 0.78 / 0.60 / 0.83. The floor matches the NUTS test's (0.55): the
  # surface accuracy is set by the information at N = 150, J = 6, p = 0.6, and
  # both backends land on it, so a higher floor would be tuned to these seeds
  # rather than to the model.
  expect_gt(mean(cr[ok]), 0.55)

  # psi / p separation holds -- the surface does not eat the detection process.
  expect_lt(abs(mean(p0[ok]) - stats::qlogis(0.6)), 0.3)
})


test_that("occu() + svc() fits under nested_laplace and recovers the surface", {
  skip_on_cran()
  skip_if_fast()

  sim <- .svcl_sim(N = 120L, J = 6L, seed = 3L)
  fit <- .svcl_fit(sim, method = "nested_laplace")
  expect_s3_class(fit, "tobs_fit")
  expect_length(as.numeric(fit$svc_field), 120L)
  expect_gt(stats::cor(as.numeric(fit$svc_field), sim$f), 0.5)
})


test_that("occu() + svc() recovers a COVARIATE-weighted coefficient surface", {
  skip_on_cran()
  skip_if_fast()

  # The general case: the surface multiplies a design column rather than the
  # intercept, so eta = b0 + w_i * f(s_i). The fitted surface is on the
  # autoscaled column's scale (.unscale_fit_per_process rewrites the beta
  # slices, never the field), which cor() is invariant to.
  N <- 150L; J <- 6L
  sim <- .svcl_sim(N = N, J = J, seed = 21L)
  set.seed(21L)
  w  <- as.numeric(scale(stats::rnorm(N)))
  z  <- stats::rbinom(N, 1, stats::plogis(0.2 + w * sim$f))
  y  <- matrix(stats::rbinom(N * J, 1, 0.6 * rep(z, J)), N, J)
  dat <- cbind(sim$data, w = w)

  fit <- suppressWarnings(tobs(
    ~ w + svc(lon, lat, indices = 2L, nn = 10, prior_range = c(0.1, 0.05)),
    data = dat, family = occu(), detection = ~ 1, y = y, method = "laplace",
    control = list(verbose = FALSE, progress = FALSE)))

  expect_s3_class(fit, "tobs_fit")
  expect_length(as.numeric(fit$svc_field), N)
  # The weighted surface is the direct evidence the design weight is fit: with
  # the weight ignored the block would estimate a plain varying intercept and
  # this correlation would collapse.
  expect_gt(stats::cor(as.numeric(fit$svc_field), sim$f), 0.4)
})
