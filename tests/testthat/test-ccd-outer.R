# Mode-centred CCD for the outer field-hyperparameter integration. The helper
# reuses tulpa's exported CCD primitives; these tests pin its moment-recovery
# and decline behaviour (fast, no model fitting), plus an end-to-end
# engage/decline check on an areal family (slow).

test_that(".tobs_ccd_outer_grid reproduces the Gaussian outer moments", {
  # Known outer posterior in u = (log tau, rho) ~ N(mu, S): the CCD design must
  # recover the mode and the analytic integral moments exactly.
  mu <- c(log(5), 0.6)
  S  <- matrix(c(0.25, 0.05, 0.05, 0.02), 2, 2)
  Q  <- solve(S)
  eval_logm <- function(theta) {
    u <- c(log(theta[1]), theta[2]); d <- u - mu
    as.numeric(-0.5 * t(d) %*% Q %*% d)
  }
  axes <- list(
    tulpaObs:::.tobs_ccd_axis("tau", "log",      lower = 0.3, upper = 30,   start = 3),
    tulpaObs:::.tobs_ccd_axis("rho", "identity", lower = 0.1, upper = 0.95, start = 0.5))

  set.seed(1)
  cc <- tulpaObs:::.tobs_ccd_outer_grid(eval_logm, axes, k_samples = 200L)
  expect_false(is.null(cc))
  expect_equal(unname(cc$u_hat), mu, tolerance = 1e-3)        # mode

  lm <- apply(cc$nodes, 1L, eval_logm)
  w  <- cc$dnode * exp(lm - max(lm)); w <- w / sum(w)
  expect_equal(sum(w * log(cc$nodes[, 1])), mu[1], tolerance = 1e-2)   # E[log tau]
  expect_equal(sum(w * cc$nodes[, 2]),      mu[2], tolerance = 1e-2)   # E[rho]
  # E[tau] = exp(mu1 + 0.5 S11) for a log-normal margin.
  expect_equal(sum(w * cc$nodes[, 1]), exp(mu[1] + 0.5 * S[1, 1]), tolerance = 0.2)
  expect_true(is.finite(cc$pareto_k))                         # diagnosed
})

test_that(".tobs_ccd_outer_grid declines on a flat or single-axis problem", {
  axes2 <- list(
    tulpaObs:::.tobs_ccd_axis("tau", "log",      lower = 0.3, upper = 30,   start = 3),
    tulpaObs:::.tobs_ccd_axis("rho", "identity", lower = 0.1, upper = 0.95, start = 0.5))
  # Uninformative (flat) outer log-marginal -> no usable curvature -> decline.
  expect_null(tulpaObs:::.tobs_ccd_outer_grid(function(theta) 0, axes2,
                                              diagnose_k = FALSE))
  # A single positive axis is left to the 1D grid (CCD declines by design).
  axes1 <- list(tulpaObs:::.tobs_ccd_axis("tau", "log", lower = 0.3, upper = 30, start = 3))
  expect_null(tulpaObs:::.tobs_ccd_outer_grid(function(theta) -(log(theta[1]))^2,
                                              axes1, diagnose_k = FALSE))
})

test_that("areal fit: default grid and opt-in CCD both recover the field + slope", {
  skip_if_fast()
  skip_on_cran()
  skip_if_not_installed("MASS")
  .adj <- function(side) {
    ng <- side*side; co <- expand.grid(x = seq_len(side), y = seq_len(side))
    a <- matrix(0L, ng, ng)
    for (i in seq_len(ng)) for (j in seq_len(ng))
      if (i != j && abs(co$x[i]-co$x[j]) + abs(co$y[i]-co$y[j]) == 1L) a[i,j] <- 1L
    a
  }
  .sim <- function(adj, cuts, seed, sd_phi = 0.6) {
    set.seed(seed); ng <- nrow(adj); nb <- length(cuts) - 1L
    phi <- as.numeric(scale(rnorm(ng)))
    for (rep in 1:3) { pn <- phi; for (i in seq_len(ng)) { nb_i <- which(adj[i,]==1L); pn[i] <- 0.5*phi[i]+0.5*mean(phi[nb_i]) }; phi <- pn }
    phi <- sd_phi*as.numeric(scale(phi)); phi <- phi-mean(phi)
    x_ab <- rnorm(ng); x_sg <- rnorm(ng)
    lam <- exp(log(50) + 0.5*x_ab + phi); sg <- exp(log(0.4) + 0.2*x_sg)
    y <- matrix(0L, ng, nb)
    for (i in seq_len(ng)) {
      pib <- tulpaObs:::.distance_pi(sg[i], cuts, "halfnorm", "line")
      probs <- c(pib, max(1-sum(pib),0)); N <- rpois(1, lam[i])
      if (N>0) { cc <- rmultinom(1,N,probs); y[i,] <- cc[seq_len(nb)] }
    }
    list(y = y, data = data.frame(abund_cov1 = x_ab, sigma_cov1 = x_sg), phi = phi)
  }
  cuts <- seq(0, 1, length.out = 6); adj <- .adj(7L)
  fam  <- distance(key = "halfnorm", transect = "line", cutpoints = cuts)
  sim  <- .sim(adj, cuts, seed = 402)
  fm   <- ~ abund_cov1 + car_proper(graph = adj)

  fg <- tobs(formula = fm, data = sim$data, family = fam, detection = ~ sigma_cov1,
             y = sim$y, method = "nested_laplace",
             control = list(progress = FALSE, verbose = FALSE))           # default grid
  fc <- tobs(formula = fm, data = sim$data, family = fam, detection = ~ sigma_cov1,
             y = sim$y, method = "nested_laplace",
             control = list(progress = FALSE, verbose = FALSE, integration = "ccd"))

  # Default path is the fixed grid; the CCD path either engages or declines to it.
  expect_identical(fg$spatial_integration, "grid")
  expect_true(fc$spatial_integration %in% c("ccd", "grid"))
  # Both recover the abundance slope (truth 0.5) and the field.
  for (f in list(fg, fc)) {
    expect_lt(abs(f$means[["lambda_abund_cov1"]] - 0.5), 0.25)
    expect_gt(cor(f$spatial_field, sim$phi), 0.5)
  }
  # CCD must not move the estimate materially versus the grid.
  expect_lt(abs(fg$means[["lambda_abund_cov1"]] - fc$means[["lambda_abund_cov1"]]), 0.1)
})
