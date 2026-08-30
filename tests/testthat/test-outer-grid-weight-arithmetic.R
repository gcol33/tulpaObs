# The reported hyperparameter posterior is a weighted sum over the outer grid,
# and the weights carry the declared prior as well as the likelihood. This file
# recomputes both from the raw nodes rather than from the code that produced
# them, so a consumer that rebuilds the weights its own way -- and drops the
# prior doing it -- fails here rather than reporting a posterior against a
# measure the fit never integrated.

.ogw_adj <- function(g) {
  N <- g * g
  A <- matrix(0L, N, N)
  idx <- function(r, c) (r - 1L) * g + c
  for (r in seq_len(g)) for (c in seq_len(g)) {
    s <- idx(r, c)
    if (r > 1L) A[s, idx(r - 1L, c)] <- 1L
    if (r < g) A[s, idx(r + 1L, c)] <- 1L
    if (c > 1L) A[s, idx(r, c - 1L)] <- 1L
    if (c < g) A[s, idx(r, c + 1L)] <- 1L
  }
  A
}

.ogw_fit <- function(seed = 4242L) {
  adj <- .ogw_adj(6L)
  N   <- nrow(adj); J <- 3L
  sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal", adj = adj,
                             sigma = 0.7, alpha = 1, seed = seed)
  long <- data.frame(site_id = rep(seq_len(N), each = J),
                     visit   = rep(seq_len(J), times = N),
                     y = as.vector(t(sim$y)),
                     det_cov1 = sim$visit_data$det_cov1,
                     pos_cov1 = sim$visit_data$pos_cov1)
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
  fit <- tobs(
    ~ occ_cov1 + icar(graph = adj),
    data = cbind(data.frame(site_id = seq_len(N)), sim$data),
    family = occu_cover("lognormal"), detection = ~ det_cov1,
    positive = ~ pos_cov1 + copy(spatial(), alpha = grid(
      c(0, exp(seq(log(0.1), log(3), length.out = 5))))),
    y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "nested_laplace",
    control = list(engine = "joint", verbose = FALSE, sigma.beta = 5,
                   sigma.grid = exp(seq(log(0.2), log(1.6), length.out = 5)),
                   adaptive.grid = FALSE, var.of.means.consistency = FALSE))
  list(fit = fit, adj = adj)
}

test_that("outer-grid weights are the declared measure, node by node", {
  skip_on_cran()
  skip_if_fast()
  f <- .ogw_fit()
  j <- f$fit$joint_fit
  expect_true(!is.null(j$log_quad))
  expect_length(j$log_quad, length(j$log_marginal))

  ok <- is.finite(j$log_marginal)
  expect_true(any(ok))

  # The engine's own weights, rebuilt here from the two pieces it stores.
  lw  <- j$log_marginal[ok] + j$log_quad[ok]
  ext <- exp(lw - max(lw)); ext <- ext / sum(ext)
  eng <- j$weights[ok]; eng <- eng / sum(eng)
  expect_equal(eng, ext, tolerance = 1e-10)

  # And they are NOT the likelihood alone: an evenly spaced grid would make the
  # two agree, so the fixture's copy axis (a point mass at zero plus a spaced
  # continuum) is what gives this check teeth.
  flat <- exp(j$log_marginal[ok] - max(j$log_marginal[ok]))
  flat <- flat / sum(flat)
  expect_gt(max(abs(flat - ext)), 1e-6)
})

test_that("the reported copy-scale posterior is that weighted sum", {
  skip_on_cran()
  skip_if_fast()
  f <- .ogw_fit()
  j <- f$fit$joint_fit
  ok <- which(is.finite(j$log_marginal))
  w  <- j$weights[ok]; w[!is.finite(w) | w < 0] <- 0; w <- w / sum(w)
  a  <- as.numeric(j$theta_grid[ok, "alpha"])

  m  <- sum(w * a)
  sd <- sqrt(max(sum(w * a^2) - m^2, 0))

  expect_equal(f$fit$spatial$alpha_mean, m,  tolerance = 1e-10)
  expect_equal(unname(f$fit$means[["alpha"]]), m,  tolerance = 1e-10)
  expect_equal(unname(f$fit$sds[["alpha"]]),   sd, tolerance = 1e-10)

  # Recomputed against the likelihood alone the mean moves, so a consumer that
  # renormalises `log_marginal` instead of `weights` cannot pass both checks.
  wf <- exp(j$log_marginal[ok] - max(j$log_marginal[ok])); wf <- wf / sum(wf)
  expect_gt(abs(sum(wf * a) - m), 1e-6)

  # The declared prior the weights carry, read back off the fit.
  expect_identical(j$copy_atom$prior_mass, tulpa:::.TULPA_COPY_ATOM_MASS)
  expect_equal(j$copy_atom$posterior_mass, sum(w[a == 0]), tolerance = 1e-10)
})
