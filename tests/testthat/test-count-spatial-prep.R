# The shared preamble behind the count-marginal spatial fitters
# (gcol33/tulpaObs#229). The areal N-mixture wrappers, the SPDE wrapper and the
# removal wrappers used to inline the same coercions, dimension checks, beta
# inits and truncation rule, and the copies had drifted on the outer-grid
# defaults and on which range checks applied. These pin the settled values and
# the shared doors, so a fifth count family cannot reopen the drift.

cs_design <- function(n_sites = 6L, n_visits = 3L, seed = 1L) {
  set.seed(seed)
  X_lambda <- cbind(`(Intercept)` = rep(1, n_sites), x = rnorm(n_sites))
  site_idx <- rep(seq_len(n_sites), each = n_visits)
  X_p <- cbind(`(Intercept)` = rep(1, n_sites * n_visits))
  list(X_lambda = X_lambda, X_p = X_p, site_idx = site_idx,
       y = rpois(n_sites * n_visits, 3))
}

cs_prep <- function(d, ...) {
  .count_spatial_prep(d$y, d$site_idx, d$X_lambda, d$X_p, "P",
                      NULL, NULL, NULL, NULL, ...)
}

test_that("the outer-grid defaults are declared once per axis", {
  expect_equal(.count_spatial_default_grid("tau_icar"),
               exp(seq(log(0.3), log(30), length.out = 9L)))
  expect_equal(.count_spatial_default_grid("tau_car"),
               exp(seq(log(0.3), log(30), length.out = 7L)))
  expect_equal(.count_spatial_default_grid("rho_car"),
               c(0.1, 0.3, 0.5, 0.75, 0.95))
  expect_equal(.count_spatial_default_grid("sigma_bym2"),
               exp(seq(log(0.2), log(3), length.out = 5L)))
  expect_equal(.count_spatial_default_grid("rho_bym2"),
               c(0.05, 0.3, 0.5, 0.7, 0.95))
  expect_error(.count_spatial_default_grid("tau"), "Unknown count-spatial grid axis")
})

test_that("the tau axis is coarser on the 2D grid than on the 1D grid", {
  # The reason 9 and 7 are both defensible, and so the reason the choice has to
  # be a property of the axis rather than of which wrapper the caller reached:
  # the proper-CAR grid is a product over both axes.
  n_icar <- length(.count_spatial_default_grid("tau_icar"))
  n_car  <- length(.count_spatial_default_grid("tau_car")) *
            length(.count_spatial_default_grid("rho_car"))
  expect_lt(length(.count_spatial_default_grid("tau_car")), n_icar)
  expect_gt(n_car, n_icar)
})

test_that("both count families read the same default axes", {
  # removal_laplace_car_proper() used to default to 9 tau nodes and a 6-point
  # even rho grid where nmix_laplace_car_proper() used 7 and 5; neither door
  # recorded which value was chosen and which was inherited. Settled on the
  # nmix spelling for both, so the axis decides.
  src <- function(f) paste(deparse(body(f)), collapse = " ")
  for (pair in list(c("tau_car", "rho_car"), c("sigma_bym2", "rho_bym2"))) {
    for (axis in pair) {
      hits <- vapply(list(nmix_laplace_car_proper, removal_laplace_car_proper,
                          nmix_laplace_bym2, removal_laplace_bym2),
                     function(f) grepl(axis, src(f), fixed = TRUE), logical(1))
      expect_equal(sum(hits), 2L, info = axis)
    }
  }
  # No fitter carries a literal grid of its own any more.
  for (f in list(nmix_laplace_icar, nmix_laplace_car_proper, nmix_laplace_bym2,
                 removal_laplace_icar, removal_laplace_car_proper,
                 removal_laplace_bym2)) {
    expect_false(grepl("seq(log(0.3), log(30)", src(f), fixed = TRUE))
    expect_false(grepl("seq(log(0.2), log(3)", src(f), fixed = TRUE))
  }
})

test_that("the grid range check reports the open and closed intervals", {
  expect_identical(.count_spatial_check_grid(c(1, 2), "tau_grid", 0, Inf),
                   c(1, 2))
  expect_error(.count_spatial_check_grid(c(0.5, 1), "rho_grid", 0, 1),
               "strictly in (0, 1)", fixed = TRUE)
  expect_silent(.count_spatial_check_grid(c(0, 1), "rho_grid", 0, 1, open = FALSE))
  expect_error(.count_spatial_check_grid(c(0, 1.5), "rho_grid", 0, 1, open = FALSE),
               "in [0, 1]", fixed = TRUE)
  expect_error(.count_spatial_check_grid(c(0.5, 1), "rho_grid", 0, 1,
                                         hint = "Use spatial_car_proper()."),
               "Use spatial_car_proper().", fixed = TRUE)
  for (bad in list(numeric(0), NULL, c(1, NA), "1")) {
    expect_error(.count_spatial_check_grid(bad, "tau_grid", 0, Inf),
                 "non-empty numeric vector")
  }
})

test_that("the proper-CAR rho range check now applies at both doors", {
  # nmix_laplace_car_proper() rejected rho outside (0, 1); the removal copy did
  # not validate it at all, so rho = 1 reached the dense Cholesky behind
  # log|Q(rho)| at a rank-deficient Q.
  d <- cs_design()
  csr <- adjacency_to_csr(matrix(0, 6L, 6L) + (abs(outer(1:6, 1:6, "-")) == 1))
  for (f in list(nmix_laplace_car_proper, removal_laplace_car_proper)) {
    expect_error(
      f(d$y, d$site_idx, seq_len(6L), d$X_lambda, d$X_p,
        csr$row_ptr, csr$col_idx, csr$n_neighbors, 6L, rho_grid = c(0.5, 1)),
      "strictly in (0, 1)", fixed = TRUE)
  }
})

test_that("the BYM2 range checks now apply at both doors", {
  d <- cs_design()
  csr <- adjacency_to_csr(matrix(0, 6L, 6L) + (abs(outer(1:6, 1:6, "-")) == 1))
  for (f in list(nmix_laplace_bym2, removal_laplace_bym2)) {
    expect_error(
      f(d$y, d$site_idx, seq_len(6L), d$X_lambda, d$X_p,
        csr$row_ptr, csr$col_idx, csr$n_neighbors, 6L, sigma_grid = c(0, 1)),
      "strictly in (0, Inf)", fixed = TRUE)
    expect_error(
      f(d$y, d$site_idx, seq_len(6L), d$X_lambda, d$X_p,
        csr$row_ptr, csr$col_idx, csr$n_neighbors, 6L, rho_grid = c(0.5, 1.2)),
      "in [0, 1]", fixed = TRUE)
  }
})

test_that("the truncation floor is the family's own rule", {
  d <- cs_design()
  # N-mixture: the largest count at a single site.
  expect_identical(.count_K_floor_max_y(c(2L, 7L, 3L), c(1L, 1L, 2L), 2L),
                   list(value = 7L, label = "max(y)"))
  # Removal: the site's passes deplete, so the latent N clears their TOTAL.
  tot <- .count_K_floor_site_total(c(2L, 7L, 3L), c(1L, 1L, 2L), 2L)
  expect_identical(tot$value, 9L)
  expect_match(tot$label, "per-site removal total")
  # A site absent from site_idx contributes a zero total, not an NA.
  expect_identical(.count_K_floor_site_total(c(2L, 3L), c(1L, 1L), 3L)$value, 5L)
})

test_that("the default truncation clears the floor by the same headroom", {
  for (floor in list(.count_K_floor_max_y(c(1L, 4L), c(1L, 2L), 2L),
                     .count_K_floor_site_total(c(1L, 4L), c(1L, 1L), 1L))) {
    expect_identical(.count_spatial_K_max(NULL, floor),
                     as.integer(floor$value + 100L))
    expect_identical(.count_spatial_K_max(floor$value, floor),
                     as.integer(floor$value))
    expect_error(.count_spatial_K_max(floor$value - 1L, floor),
                 floor$label, fixed = TRUE)
  }
})

test_that("the preamble validates the designs, the map and the graph", {
  d <- cs_design()
  expect_error(cs_prep(list(y = d$y, site_idx = d$site_idx,
                            X_lambda = as.data.frame(d$X_lambda), X_p = d$X_p)),
               "`X_lambda` must be a numeric matrix")
  bad <- d; bad$y <- d$y[-1L]
  expect_error(cs_prep(bad), "length(y) must equal nrow(X_p)", fixed = TRUE)
  bad <- d; bad$site_idx <- d$site_idx[-1L]
  expect_error(cs_prep(bad), "length(site_idx) must equal nrow(X_p)", fixed = TRUE)
  expect_error(cs_prep(d, map_site_to_unit = seq_len(5L), n_spatial = 6L),
               "must equal nrow(X_lambda)", fixed = TRUE)
  expect_error(cs_prep(d, map_site_to_unit = c(1:5, 9L), n_spatial = 6L),
               "must lie in [1, n_spatial]", fixed = TRUE)
  expect_error(cs_prep(d, n_spatial = 6L, adj_row_ptr = integer(6L),
                       n_neighbors = integer(6L)),
               "must equal n_spatial + 1", fixed = TRUE)
  expect_error(cs_prep(d, n_spatial = 6L, adj_row_ptr = integer(7L),
                       n_neighbors = integer(5L)),
               "must equal n_spatial", fixed = TRUE)
})

test_that("the preamble resolves the coefficient warm starts", {
  d <- cs_design()
  pp <- cs_prep(d)
  expect_identical(pp$n_sites, 6L)
  expect_identical(pp$p_lam, 2L)
  expect_equal(pp$beta_lambda_init, c(log(mean(d$y) + 0.1), 0))
  expect_equal(pp$beta_p_init, 0)
  expect_identical(pp$r_grid, Inf)
  expect_error(
    .count_spatial_prep(d$y, d$site_idx, d$X_lambda, d$X_p, "P",
                        c(0, 0, 0), NULL, NULL, NULL),
    "must equal ncol(X_lambda)", fixed = TRUE)
  expect_error(
    .count_spatial_prep(d$y, d$site_idx, d$X_lambda, d$X_p, "P",
                        NULL, c(0, 0), NULL, NULL),
    "must equal ncol(X_p)", fixed = TRUE)
})

test_that("a latent warm start is checked against the length its door names", {
  d <- cs_design()
  expect_error(cs_prep(d, latent = list(z_init = rep(0, 5)),
                       n_latent = c(n_spatial = 6L)),
               "length(z_init) must equal n_spatial", fixed = TRUE)
  expect_error(cs_prep(d, latent = list(u_init = rep(0, 5)),
                       n_latent = c(n_mesh = 9L)),
               "length(u_init) must equal n_mesh", fixed = TRUE)
  # NULL is the absent warm start, not a length-0 one.
  expect_silent(cs_prep(d, latent = list(v_init = NULL, w_init = rep(0, 6)),
                        n_latent = c(n_spatial = 6L)))
})

test_that("the K_max boundary warning is one wording at every door", {
  src <- function(f) paste(deparse(body(f)), collapse = " ")
  fitters <- list(nmix_laplace_icar, nmix_laplace_car_proper, nmix_laplace_bym2,
                  nmix_laplace_spde, removal_laplace_icar,
                  removal_laplace_car_proper, removal_laplace_bym2)
  for (f in fitters) {
    expect_true(grepl(".count_spatial_warn_boundary", src(f), fixed = TRUE))
    expect_false(grepl("Max posterior weight on N = K_max", src(f), fixed = TRUE))
  }
  expect_warning(.count_spatial_warn_boundary(list(boundary_max = c(1e-3, 0))),
                 "raise K_max")
  expect_silent(.count_spatial_warn_boundary(list(boundary_max = c(1e-9, NA))))
})
