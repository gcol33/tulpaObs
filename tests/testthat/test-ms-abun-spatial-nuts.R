# Community / multispecies N-mixture NUTS + a shared areal field on the abundance
# arm (ms_abun() + icar/bym2/car_proper under method = "nuts"; gcol33/tulpaObs#73,
# #113). The #73 sampler carries a FIXED-HYPER non-centered shared field: the
# field precision (tau, rho) is fixed at the sfMsNMix nested-Laplace estimate and
# the whitened raw ~ N(0, I) is sampled jointly with the community means, per-
# species deviations, and community covariances. car_proper uses a square inverse-
# Cholesky loading; the intrinsic icar / bym2 fields use the #71 sum-to-zero
# eigen-loading that drops the constant precision null direction (n_raw < n_units),
# so they sample with the same well-conditioned geometry (0 divergences) rather
# than routing to nested_laplace. This path is wired at R/tobs_dispatch.R but was
# previously untested; the recovery invariant is that NUTS reproduces the
# integrated shared field (cor high), recovers the community abundance slope, and
# samples cleanly.
#
# Each fit warm-starts from the sfMsNMix nested Laplace-EM (S species x an EM x an
# outer grid) then runs the sampler, so every block is skip_on_cran() + skip_if_fast().

# Dense rook (4-neighbour) adjacency for a g x g grid.
.man_rook_adj <- function(g) {
  n <- g * g
  A <- matrix(0L, n, n)
  idx <- function(r, c) (r - 1L) * g + c
  for (r in seq_len(g)) for (c in seq_len(g)) {
    i <- idx(r, c)
    if (r > 1) A[i, idx(r - 1L, c)] <- 1L
    if (r < g) A[i, idx(r + 1L, c)] <- 1L
    if (c > 1) A[i, idx(r, c - 1L)] <- 1L
    if (c < g) A[i, idx(r, c + 1L)] <- 1L
  }
  A
}

# Fit a shared-field ms_abun() under NUTS via the front door; returns the fit.
.man_fit_nuts <- function(sim, adj, field_kind, seed = 1L,
                          n.iter = 500L, n.warmup = 500L) {
  form <- switch(field_kind,
    icar       = ~ abund_cov1 + icar(graph = adj),
    bym2       = ~ abund_cov1 + bym2(graph = adj),
    car_proper = ~ abund_cov1 + car_proper(graph = adj))
  tobs(form, detection = ~ det_cov1, family = ms_abun(),
       data = sim$data, y = sim$y, species = sim$species,
       method = "nuts",
       control = list(verbose = FALSE, progress = FALSE,
                      n.iter = n.iter, n.warmup = n.warmup, seed = seed))
}


test_that("ms_abun() NUTS + car_proper shared field recovers means + field (#73)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .man_rook_adj(6L)
  sim <- simulate_ms_abun(n_species = 12, J = 5, n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = c(log(4), 0.5), mu_p = c(0.3, -0.3),
                          sd_lambda = 0.5, sd_p = 0.4,
                          graph = adj, sigma.field = 0.6, seed = 7)
  # Reference nested-Laplace (sfMsNMix) fit -- NUTS should reproduce its field.
  nl <- tobs(~ abund_cov1 + car_proper(graph = adj), detection = ~ det_cov1,
             family = ms_abun(), data = sim$data, y = sim$y, species = sim$species,
             method = "nested_laplace", control = list(verbose = FALSE, progress = FALSE))
  nu <- .man_fit_nuts(sim, adj, "car_proper", seed = 3L)

  expect_identical(nu$method, "nuts")
  expect_equal(mean(nu$nuts$divergent), 0)                        # well-conditioned
  expect_length(nu$spatial_field, nrow(adj))
  expect_gt(cor(nu$spatial_field, nl$spatial_field), 0.9)         # reproduces NL field
  expect_gt(cor(nu$spatial_field, sim$truth$field), 0.7)          # tracks truth

  # Community abundance slope (the strongly-identified coordinate) recovered.
  truth <- c(sim$truth$mu_lambda, sim$truth$mu_p)
  names(truth) <- names(nu$means)
  islope <- 2L   # mu_lambda covariate slope
  expect_lt(abs(nu$means[islope] - truth[islope]) / nu$sds[islope], 3)
})


test_that("ms_abun() NUTS + icar shared field samples clean + centred (#113)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .man_rook_adj(6L)
  sim <- simulate_ms_abun(n_species = 12, J = 5, n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = c(log(4), 0.5), mu_p = c(0.3, -0.3),
                          sd_lambda = 0.5, sd_p = 0.4,
                          graph = adj, sigma.field = 0.6, seed = 7)
  nl <- tobs(~ abund_cov1 + icar(graph = adj), detection = ~ det_cov1,
             family = ms_abun(), data = sim$data, y = sim$y, species = sim$species,
             method = "nested_laplace", control = list(verbose = FALSE, progress = FALSE))
  nu <- .man_fit_nuts(sim, adj, "icar", seed = 3L)

  expect_identical(nu$method, "nuts")
  expect_equal(mean(nu$nuts$divergent), 0)                        # sum-to-zero geometry
  expect_length(nu$spatial_field, nrow(adj))
  expect_lt(abs(mean(nu$spatial_field)), 1e-6)                    # sum z = 0 (intrinsic)
  expect_gt(cor(nu$spatial_field, nl$spatial_field), 0.85)       # reproduces NL field
  expect_gt(cor(nu$spatial_field, sim$truth$field), 0.7)

  truth <- c(sim$truth$mu_lambda, sim$truth$mu_p)
  expect_lt(abs(nu$means[2L] - truth[2L]) / nu$sds[2L], 3)
})


test_that("ms_abun() NUTS + bym2 shared field samples clean + recovers field (#113)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .man_rook_adj(6L)
  sim <- simulate_ms_abun(n_species = 12, J = 5, n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = c(log(4), 0.5), mu_p = c(0.3, -0.3),
                          sd_lambda = 0.5, sd_p = 0.4,
                          graph = adj, sigma.field = 0.6, seed = 9)
  nl <- tobs(~ abund_cov1 + bym2(graph = adj), detection = ~ det_cov1,
             family = ms_abun(), data = sim$data, y = sim$y, species = sim$species,
             method = "nested_laplace", control = list(verbose = FALSE, progress = FALSE))
  nu <- .man_fit_nuts(sim, adj, "bym2", seed = 4L)

  expect_identical(nu$method, "nuts")
  expect_equal(mean(nu$nuts$divergent), 0)
  expect_length(nu$spatial_field, nrow(adj))
  # bym2 = structured phi + iid theta; the reconstructed unit field is not
  # sum-to-zero, so assert on shape recovery, not centring.
  expect_gt(cor(nu$spatial_field, nl$spatial_field), 0.8)
  expect_gt(cor(nu$spatial_field, sim$truth$field), 0.6)

  truth <- c(sim$truth$mu_lambda, sim$truth$mu_p)
  expect_lt(abs(nu$means[2L] - truth[2L]) / nu$sds[2L], 3)
})
