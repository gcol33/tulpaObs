# Opt-in exact-Newton inner solver for the areal shared-field community N-mixture
# (control$inner.solver = "newton"). Same model as the default Laplace-EM path --
# the field hyperparameter is outer-grid integrated either way -- but the inner
# step alternates a tulpa AGHQ community solve with an exact-Newton shared-field
# solve. The two solvers must produce the same fit (community means, shared field,
# fit shape); only the numbers differ to solver tolerance. The Newton path runs an
# FD-gradient AGHQ profile loop per grid node, so it is much slower than EM: the
# equivalence block is skip_if_fast() and drives the internal driver with a coarse
# grid to stay affordable in CI.

test_that("inner.solver = \"newton\" is a recognized control with regime guards", {
  skip_on_cran()
  skip_if_fast()
  adj <- rook_adj(3L)
  sim <- simulate_ms_abun(n_species = 4, J = 2, graph = adj,
                          sigma.field = 0.3, seed = 5)
  # The control key is admitted (not an "unknown control option"); an unknown
  # VALUE hits match.arg -- both fire before any heavy fitting.
  expect_error(
    tobs(~ abund_cov1 + icar(graph = adj), detection = ~ det_cov1,
         family = ms_abun(), data = sim$data, y = sim$y, species = sim$species,
         method = "nested_laplace", control = list(inner.solver = "nope")),
    "should be one of|'arg'")
  # Newton is Poisson-only: negbin + newton errors clearly.
  expect_error(
    tobs(~ abund_cov1 + icar(graph = adj), detection = ~ det_cov1,
         family = ms_abun(mixture = "negbin"), data = sim$data, y = sim$y,
         species = sim$species, method = "nested_laplace",
         control = list(verbose = FALSE, inner.solver = "newton")),
    "Poisson-only")
})

test_that("the Newton solver agrees with EM on the shared-field community fit", {
  skip_on_cran()
  skip_if_fast()
  adj <- rook_adj(4L)                       # 16 sites
  sim <- simulate_ms_abun(n_species = 8, J = 4, n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = c(log(4), 0.5), mu_p = c(0.3, -0.3),
                          sd_lambda = 0.5, sd_p = 0.4,
                          graph = adj, sigma.field = 0.6, seed = 7)
  truth <- c(sim$truth$mu_lambda, sim$truth$mu_p)

  # Build the model once; drive both inner solvers on the SAME data. The Newton
  # driver is called directly with a coarse (tau) grid so the FD-AGHQ profile
  # loop stays affordable (the front-door default grid is minutes per fit).
  model   <- .tobs_build_ms_abun(abund_formula = ~ abund_cov1 + icar(graph = adj),
                                 det_formula = ~ det_cov1, data = sim$data,
                                 y = sim$y, species = sim$species)
  spatial <- .tobs_structures_from_model(model)$spatial

  fit_em <- .tobs_fit_ms_nmix_spatial(model, spatial, inner_solver = "em",
                                      verbose = FALSE)
  fit_nt <- .tobs_fit_ms_nmix_spatial_newton(model, spatial, tau_grid = c(1, 5),
                                             inner_iter = 1L, inner_maxit = 10L,
                                             verbose = FALSE)

  # Identical fit shape; the inner solver is recorded on each.
  expect_identical(sort(names(fit_em)), sort(names(fit_nt)))
  expect_identical(fit_em$ms_community$optimizer, "em")
  expect_identical(fit_nt$ms_community$optimizer, "newton")
  expect_identical(fit_nt$method, "nested_laplace")
  expect_length(fit_nt$spatial_field, nrow(adj))
  expect_true(all(c("tau", "sigma") %in% names(fit_nt$ms_hyper)))

  # Same model, two inner solvers -> community means + shared field agree. The
  # coarse two-node Newton grid integrates the field hyperparameter loosely
  # (the full grids agree to <0.05; here ~0.3), so the means bound is a gross-
  # breakage backstop; the field correlation and the recovery check below are
  # the sharp guards (a wrong field scaling / offset would crater both).
  expect_lt(max(abs(fit_em$means - fit_nt$means)), 0.4)
  expect_gt(cor(fit_em$spatial_field, fit_nt$spatial_field), 0.9)
  # Newton recovers the community means about as well as EM.
  expect_true(all(abs(fit_nt$means - truth) / fit_nt$sds < 3.5))
})
