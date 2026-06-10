# =============================================================================
# test-ms-occu-spatial.R - areal-spatial community single-season occupancy
# (ms_occu() + shared field; the occupancy analogue of sfMsNMix, tulpaObs#75).
#
# A per-species two-state occupancy model with Gaussian community hyperpriors on
# the per-species coefficients AND one shared ICAR / BYM2 / proper-CAR field on
# the occupancy arm, fit by the in-tree community-spatial nested Laplace-EM
# (R/ms_occu_spatial.R, src/ms_occu_spatial.cpp).
#
# Coverage: (1) the per-site occupancy cell score + observed -Hessian byte-exact
# vs the closed-form marginal (FD), (2) community-mean + shared-field recovery
# (ICAR), (3) the proper-CAR + BYM2 field kinds run + recover the field, (4) S3
# (fitted carries the field), (5) gates (field on detection / NUTS + field).
# =============================================================================


# --- shared fixtures -------------------------------------------------------

.msocs_grid_graph <- function(side) {
  N <- side * side
  A <- matrix(0L, N, N)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    i <- idx(r, c)
    if (r < side) { j <- idx(r + 1L, c); A[i, j] <- 1L; A[j, i] <- 1L }
    if (c < side) { j <- idx(r, c + 1L); A[i, j] <- 1L; A[j, i] <- 1L }
  }
  A
}

# Community occupancy with a smooth shared field on the occupancy arm.
.msocs_sim <- function(side = 8L, J = 4L, n_species = 14L,
                       mu_psi = c(0, 0.5), mu_p = 0.2,
                       sd_psi = c(0.5, 0.3), sd_p = 0.4,
                       field_sd = 0.8, seed = 1L) {
  set.seed(seed)
  A <- .msocs_grid_graph(side); N <- nrow(A)
  coord <- expand.grid(r = seq_len(side), c = seq_len(side))
  f <- field_sd * scale(sin(coord$r / side * pi) + cos(coord$c / side * pi))[, 1]
  f <- f - mean(f)
  data <- data.frame(x = stats::rnorm(N))
  X_occ <- stats::model.matrix(~ x, data); X_det <- stats::model.matrix(~ 1, data)
  beta_psi <- cbind(stats::rnorm(n_species, mu_psi[1], sd_psi[1]),
                    stats::rnorm(n_species, mu_psi[2], sd_psi[2]))
  beta_p   <- matrix(stats::rnorm(n_species, mu_p, sd_p), n_species, 1)
  y <- array(NA_integer_, dim = c(N, J, n_species))
  for (s in seq_len(n_species)) {
    psi <- stats::plogis(as.numeric(X_occ %*% beta_psi[s, ]) + f)
    p   <- stats::plogis(beta_p[s, 1])
    z <- stats::rbinom(N, 1, psi)
    for (i in seq_len(N)) y[i, , s] <- stats::rbinom(J, 1, z[i] * p)
  }
  list(y = y, data = data, graph = A, field = f,
       beta_psi = beta_psi, species = paste0("sp", seq_len(n_species)))
}


# --- (1) per-site cell score + observed -Hessian vs FD ---------------------

test_that("ms_occu spatial cell score + curvature match finite differences", {
  skip_on_cran()
  ll_site <- function(eta_psi, eta_p, nv, nd) {
    clp <- function(x) pmin(pmax(x, 1e-12), 1 - 1e-12)
    psi <- clp(stats::plogis(eta_psi)); p <- clp(stats::plogis(eta_p))
    if (nd > 0) log(psi) + nd * log(p) + (nv - nd) * log1p(-p)
    else log(psi * (1 - p)^nv + (1 - psi))
  }
  max_g <- 0; max_h <- 0
  grid <- expand.grid(eta_psi = c(-1.2, 0, 0.8), eta_p = c(-0.5, 0.3, 1.0),
                      nv = c(2L, 4L), nd = c(0L, 1L, 2L))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]; if (g$nd > g$nv) next
    cell <- cpp_ms_occu_site_cell(g$eta_psi, g$eta_p, g$nv, g$nd, observed = TRUE)
    h <- 1e-6
    gx <- (ll_site(g$eta_psi + h, g$eta_p, g$nv, g$nd) -
           ll_site(g$eta_psi - h, g$eta_p, g$nv, g$nd)) / (2 * h)
    gp <- (ll_site(g$eta_psi, g$eta_p + h, g$nv, g$nd) -
           ll_site(g$eta_psi, g$eta_p - h, g$nv, g$nd)) / (2 * h)
    max_g <- max(max_g, abs(cell$grad[1] - gx), abs(cell$grad[2] - gp))
    hh <- 1e-4
    f <- function(a, b) ll_site(a, b, g$nv, g$nd)
    hxx <- (f(g$eta_psi + hh, g$eta_p) - 2 * f(g$eta_psi, g$eta_p) +
            f(g$eta_psi - hh, g$eta_p)) / hh^2
    hpp <- (f(g$eta_psi, g$eta_p + hh) - 2 * f(g$eta_psi, g$eta_p) +
            f(g$eta_psi, g$eta_p - hh)) / hh^2
    hxp <- (f(g$eta_psi + hh, g$eta_p + hh) - f(g$eta_psi + hh, g$eta_p - hh) -
            f(g$eta_psi - hh, g$eta_p + hh) + f(g$eta_psi - hh, g$eta_p - hh)) / (4 * hh^2)
    nh <- -matrix(c(hxx, hxp, hxp, hpp), 2, 2)
    max_h <- max(max_h, max(abs(cell$neg_hess - nh)))
  }
  expect_lt(max_g, 1e-6)
  expect_lt(max_h, 1e-4)
})


# --- (2) community-mean + shared-field recovery (ICAR) ---------------------

test_that("ms_occu + icar() recovers community means + the shared field", {
  skip_on_cran()
  skip_if_fast()
  sim <- .msocs_sim(side = 8L, J = 4L, n_species = 14L, seed = 11L)
  fit <- tobs(~ x + icar(graph = sim$graph), data = sim$data, family = ms_occu(),
              detection = ~ 1, y = sim$y, species = sim$species,
              method = "nested_laplace", control = list(verbose = FALSE))
  expect_equal(fit$method, "nested_laplace")
  expect_false(is.null(fit$spatial_field))

  truth <- c("psi_(Intercept)" = 0, "psi_x" = 0.5, "p_(Intercept)" = 0.2)
  m <- fit$means[names(truth)]; s <- fit$sds[names(truth)]
  expect_true(all(abs(m - truth) / s < 2.5))

  # The shared occupancy field is recovered (up to the sum-to-zero constraint).
  expect_gt(cor(fit$spatial_field, sim$field), 0.80)

  cm <- fit$ms_community
  expect_gt(cor(cm$coef_psi[, 1], sim$beta_psi[, 1]), 0.45)
})


# --- (3) proper-CAR + BYM2 field kinds run + recover the field -------------

test_that("ms_occu + car_proper()/bym2() recover the shared field", {
  skip_on_cran()
  skip_if_fast()
  sim <- .msocs_sim(side = 7L, J = 4L, n_species = 12L, seed = 3L)
  for (term in c("car_proper", "bym2")) {
    f <- stats::as.formula(sprintf("~ x + %s(graph = sim$graph)", term))
    fit <- tobs(f, data = sim$data, family = ms_occu(), detection = ~ 1,
                y = sim$y, species = sim$species, method = "nested_laplace",
                control = list(verbose = FALSE))
    expect_equal(fit$method, "nested_laplace")
    expect_gt(cor(fit$spatial_field, sim$field), 0.70)
  }
})


# --- (4) S3: fitted carries the field --------------------------------------

test_that("ms_occu spatial S3 works, fitted adds the field offset", {
  skip_on_cran()
  skip_if_fast()
  sim <- .msocs_sim(side = 6L, J = 3L, n_species = 8L, seed = 7L)
  fit <- tobs(~ x + icar(graph = sim$graph), data = sim$data, family = ms_occu(),
              detection = ~ 1, y = sim$y, species = sim$species,
              method = "nested_laplace", control = list(verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_no_error(print(fit))
  fv <- fitted(fit)
  expect_equal(dim(fv$psi), c(36L, 8L))
  expect_true(all(fv$psi > 0 & fv$psi < 1))
  re <- ranef(fit)
  expect_s3_class(re, "data.frame")
  cf <- coef(fit)
  expect_setequal(names(cf), c("psi", "p"))
})


# --- (5) gates -------------------------------------------------------------

test_that("ms_occu spatial gates: field on detection / NUTS + field", {
  skip_on_cran()
  sim <- .msocs_sim(side = 5L, J = 3L, n_species = 4L, seed = 1L)
  # field on the detection formula is rejected.
  expect_error(
    tobs(~ x, data = sim$data, family = ms_occu(),
         detection = ~ icar(graph = sim$graph), y = sim$y, species = sim$species,
         method = "nested_laplace", control = list(verbose = FALSE)),
    "occupancy arm")
  # NUTS + a field points to nested_laplace.
  expect_error(
    tobs(~ x + icar(graph = sim$graph), data = sim$data, family = ms_occu(),
         detection = ~ 1, y = sim$y, species = sim$species, method = "nuts"),
    "nested_laplace")
})
