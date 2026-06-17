# occu_cover() shared spatial field + per-group random intercept on the
# occupancy arm (gcol33/tulpaObs#56, the tulpaObs consumer of tulpa#86's field +
# per-group RE composition in the joint cell-coupling engine). The RE joins the
# joint fit as one `iid` prior block whose variance integrates on the outer grid
# alongside the field sigma / alpha.

.ocfr_grid_adj <- function(side) {
  ng <- side * side
  co <- expand.grid(x = seq_len(side), y = seq_len(side))
  adj <- matrix(0L, ng, ng)
  for (i in seq_len(ng)) for (j in seq_len(ng))
    if (i != j && abs(co$x[i]-co$x[j]) + abs(co$y[i]-co$y[j]) == 1L) adj[i, j] <- 1L
  adj
}

# Simulate occu_cover with a shared smoothed field on psi (+ alpha coupling onto
# cover) and a per-group random intercept on occupancy.
.ocfr_sim <- function(seed, side = 8L, J = 5L, n_g = 10L,
                      psi0 = 0.0, psi_x = 0.7, p0 = 0.3, pos0 = log(25),
                      sigma_u = 0.7, alpha = 0.5, sigma_cov = 0.4) {
  set.seed(seed)
  adj <- .ocfr_grid_adj(side); n <- nrow(adj)
  phi <- as.numeric(scale(rnorm(n)))
  for (r in 1:4) { pn <- phi; for (i in seq_len(n)) {
    nb <- which(adj[i, ] == 1L); pn[i] <- 0.5*phi[i] + 0.5*mean(phi[nb]) }; phi <- pn }
  phi <- 0.8 * as.numeric(scale(phi)); phi <- phi - mean(phi)
  grp <- factor(rep(seq_len(n_g), length.out = n))
  u <- rnorm(n_g, 0, sigma_u); u <- u - mean(u)
  x <- rnorm(n)
  psi <- plogis(psi0 + psi_x*x + phi + u[as.integer(grp)])
  z <- rbinom(n, 1, psi); p <- plogis(p0)
  y <- matrix(0L, n, J)
  for (i in seq_len(n)) y[i, ] <- if (z[i] == 1L) rbinom(J, 1, p) else 0L
  mu_cov <- pos0 + alpha*phi
  ypos <- matrix(0, n, J)
  for (i in seq_len(n)) for (j in seq_len(J))
    if (y[i, j] == 1L) ypos[i, j] <- exp(rnorm(1, mu_cov[i], sigma_cov))
  list(adj = adj, data = data.frame(x = x, g = grp), y = y, y_pos = ypos, phi = phi)
}

test_that("occu_cover() shared field + per-group RE: fit runs and reports the RE", {
  skip_on_cran()
  sim <- .ocfr_sim(1L, side = 6L, n_g = 6L)
  fit <- tobs(~ x + icar(graph = sim$adj) + re(g), data = sim$data,
              family = occu_cover("lognormal"), detection = ~ 1, positive = ~ 1,
              y = sim$y, y_pos = sim$y_pos, method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "nested_laplace")
  # The RE block is reported: sigma_re hyperparameter + per-group BLUPs.
  expect_true("sigma_re" %in% names(fit$means))
  expect_true(is.finite(fit$means[["sigma_re"]]))
  expect_gt(fit$means[["sigma_re"]], 0)
  expect_false(is.null(fit$re))
  expect_identical(fit$re$arm, "psi")
  expect_length(fit$re$blup, 6L)
  # ranef() surfaces the BLUP table.
  rf <- ranef(fit)
  expect_s3_class(rf, "data.frame")
  expect_equal(nrow(rf), 6L)
  expect_true(all(c("group", "blup", "blup_sd") %in% names(rf)))
})

test_that("occu_cover() spatial + RE recovers the means, field, and RE variance", {
  skip_on_cran()
  skip_if_fast()
  seeds <- 1:6
  truth <- c(psi_x = 0.7, p0 = 0.3, pos0 = log(25))
  res <- t(vapply(seeds, function(s) {
    sim <- .ocfr_sim(s, side = 8L, J = 5L, n_g = 10L)
    fit <- tobs(~ x + icar(graph = sim$adj) + re(g), data = sim$data,
                family = occu_cover("lognormal"), detection = ~ 1,
                positive = ~ 1 + copy(spatial()),
                y = sim$y, y_pos = sim$y_pos, method = "nested_laplace",
                control = list(verbose = FALSE, progress = FALSE))
    c(psi_x = unname(fit$means[["psi_x"]]),
      p0    = unname(fit$means[["p_(Intercept)"]]),
      pos0  = unname(fit$means[["pos_(Intercept)"]]),
      fcor  = stats::cor(fit$spatial_field, sim$phi),
      sd_re = unname(fit$means[["sigma_re"]]))
  }, numeric(5)))
  mn <- colMeans(res)
  # Community means recover (the field absorbs the spatial structure, the RE the
  # grouped structure); tolerances accommodate small-sample binary-occupancy noise.
  expect_lt(abs(mn[["psi_x"]] - truth[["psi_x"]]), 0.25)
  expect_lt(abs(mn[["p0"]]    - truth[["p0"]]),    0.30)
  expect_lt(abs(mn[["pos0"]]  - truth[["pos0"]]),  0.15)
  # The shared field is recovered (positive correlation with the truth).
  expect_gt(mn[["fcor"]], 0.6)
  # The RE variance is detected (positive); it carries the known binary-RE
  # Laplace small-cluster attenuation, so it is a lower bound, not exact.
  expect_gt(mn[["sd_re"]], 0.15)
})

test_that("occu_cover() spatial + RE gates the unsupported configurations", {
  skip_on_cran()
  sim <- .ocfr_sim(2L, side = 5L, n_g = 5L)
  # A random slope has no scalar-per-group iid form on the shared-field engine.
  expect_error(
    tobs(~ x + icar(graph = sim$adj) + re(g, type = "slope", covariate = "x"),
         data = sim$data, family = occu_cover("lognormal"),
         detection = ~ 1, positive = ~ 1, y = sim$y, y_pos = sim$y_pos,
         method = "nested_laplace", control = list(progress = FALSE)),
    "random INTERCEPT only")
  # The v3 escape hatch has no RE block.
  expect_error(
    tobs(~ x + icar(graph = sim$adj) + re(g), data = sim$data,
         family = occu_cover("lognormal"), detection = ~ 1, positive = ~ 1,
         y = sim$y, y_pos = sim$y_pos, method = "nested_laplace",
         control = list(engine = "v3_nested", progress = FALSE)),
    "joint_coupled engine")
})
