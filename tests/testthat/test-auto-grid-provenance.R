# Declare a defaulted outer-grid axis as OURS with tulpa::auto_grid().
#
# The engine's auto-recenter decides axis PROVENANCE, not field presence: it may
# move an axis that is absent, marked, or exactly equal to its own default, and
# never moves anything else. tulpaObs writes a grid on every joint fit -- it
# derives a second axis from the first and hands the same vector to several
# blocks -- so without the mark a package default is indistinguishable from a
# user pin and the rescue goes inert. Several tulpaObs defaults (the cover
# arm-specific 0.2-2.5 axis, the RE 0.05-2 axis, the EM-path per-component
# axes) never matched the engine default, so they read as pins.
#
# The marker is an attribute, dropped by sort() / [ / c() / as.numeric() /
# expand.grid(), which is the failure mode these tests exist to catch: a site
# that reshapes a defaulted grid has to re-apply it on the reshaped vector.

skip_if_no_auto_grid <- function() {
  skip_if_not(is.function(getExportedValue("tulpa", "auto_grid")),
              "tulpa::auto_grid() not available")
}

test_that("the shared default-grid helpers declare themselves", {
  skip_if_no_auto_grid()
  expect_true(tulpa::is_auto_grid(tulpaObs:::.tobs_default_sigma_grid()))
  expect_true(tulpa::is_auto_grid(tulpaObs:::.tobs_default_alpha_grid()))
  expect_true(tulpa::is_auto_grid(tulpaObs:::.tobs_default_bym2_rho_grid()))
})

test_that("the alpha and sigma defaults read the engine's axes (#209)", {
  skip_if_no_auto_grid()
  # Both are READ from the engine, so they cannot fall out of step with the
  # nodes `.nl_axis_matches_default()` compares against. The assertions read the
  # engine too -- a literal here would drift exactly as the helpers' own copies
  # would have.
  expect_equal(as.numeric(tulpaObs:::.tobs_default_alpha_grid()),
               as.numeric(tulpa:::.nl_grid_axis("copy_alpha")))
  expect_equal(as.numeric(tulpaObs:::.tobs_default_sigma_grid()),
               as.numeric(tulpa:::.nl_grid_axis("field_sd")))
  # The copy axis carries an exact 0, so the uncoupled model is ON the grid
  # rather than a limit of it. That is the engine's `prepend`, and reading the
  # axis is what keeps it.
  expect_true(any(as.numeric(tulpaObs:::.tobs_default_alpha_grid()) == 0))

  # What the read buys, asserted as the property rather than as a value: the
  # engine recognises both axes as ITS OWN default by value on the paths these
  # grids reach, which is the belt that still holds when a reshape has dropped
  # the auto_grid() marker.
  matches <- tulpa:::.nl_axis_matches_default
  expect_true(matches(tulpaObs:::.tobs_default_sigma_grid(),
                      "sigma_grid", ".joint_areal"))
  expect_true(matches(tulpaObs:::.tobs_default_sigma_grid(),
                      "sigma_grid", ".copy"))
  expect_true(matches(tulpaObs:::.tobs_default_alpha_grid(),
                      "alpha_grid", ".copy"))
})

test_that("the proper-CAR rho default is tulpaObs's own, not an engine read (#209)", {
  skip_if_no_auto_grid()
  # Deliberately OURS. The engine's `joint_car_rho` holds the same four nodes
  # where it exists, but it is absent at this package's declared Imports floor
  # (tulpa 0.0.136), and the engine binds neither it nor the `rho_car_grid`
  # field into `.NL_FAMILY_AXES` -- so there is no value-recognition belt to
  # keep in step and no floor-safe axis to read. A literal is the right
  # assertion here for the same reason the value is a literal in the source.
  expect_equal(as.numeric(tulpaObs:::.tobs_default_rho_car_grid()),
               c(0.5, 0.8, 0.95, 0.99))
  # No engine binding, so provenance rests entirely on the marker.
  expect_true(tulpa::is_auto_grid(tulpaObs:::.tobs_default_rho_car_grid()))
  expect_false(tulpa:::.nl_axis_matches_default(
    tulpaObs:::.tobs_default_rho_car_grid(), "rho_car_grid", "car_proper"))
})

test_that(".tobs_mark_auto marks only when the caller defaulted", {
  skip_if_no_auto_grid()
  expect_true(tulpa::is_auto_grid(tulpaObs:::.tobs_mark_auto(c(0.5, 1, 2), TRUE)))
  expect_false(tulpa::is_auto_grid(tulpaObs:::.tobs_mark_auto(c(0.5, 1, 2), FALSE)))
  # Values pass through untouched either way.
  expect_equal(as.numeric(tulpaObs:::.tobs_mark_auto(c(0.5, 1, 2), TRUE)),
               c(0.5, 1, 2))
})

test_that(".tobs_num_auto coerces and carries the source's own provenance", {
  skip_if_no_auto_grid()
  marked <- tulpa::auto_grid(c(0.5, 1, 2))
  expect_true(tulpa::is_auto_grid(tulpaObs:::.tobs_num_auto(marked)))
  expect_false(tulpa::is_auto_grid(tulpaObs:::.tobs_num_auto(c(0.5, 1, 2))))
  expect_equal(as.numeric(tulpaObs:::.tobs_num_auto(marked)), c(0.5, 1, 2))
  # An integer grid still reaches the block as a double.
  expect_type(tulpaObs:::.tobs_num_auto(1:3), "double")
})

test_that("the per-front-door default axes declare themselves", {
  skip_if_no_auto_grid()
  for (f in c(".tobs_default_rho_car_grid", ".tobs_default_temporal_tau_grid",
              ".tobs_default_temporal_rho_grid", ".tobs_default_temporal_sigma_grid",
              ".tobs_default_re_sigma_grid", ".tobs_default_occu_joint_sigma_grid")) {
    expect_true(tulpa::is_auto_grid(getFromNamespace(f, "tulpaObs")()),
                info = f)
  }
  # Values are the ones the front doors carried before they were single-sourced.
  expect_equal(as.numeric(tulpaObs:::.tobs_default_rho_car_grid()),
               c(0.5, 0.8, 0.95, 0.99))
  expect_equal(as.numeric(tulpaObs:::.tobs_default_temporal_tau_grid()), c(1, 4, 16))
  expect_equal(as.numeric(tulpaObs:::.tobs_default_temporal_rho_grid()), c(0.3, 0.7))
  expect_equal(as.numeric(tulpaObs:::.tobs_default_temporal_sigma_grid()),
               exp(seq(log(0.1), log(1), length.out = 3)))
  expect_equal(as.numeric(tulpaObs:::.tobs_default_re_sigma_grid()),
               exp(seq(log(0.1), log(1.5), length.out = 3)))
  expect_equal(as.numeric(tulpaObs:::.tobs_default_occu_joint_sigma_grid()),
               exp(seq(log(0.15), log(3), length.out = 4)))
})

# --------------------------------------------------------------------------- #
# cover() multi-block: the non-spatial blocks and the copy axis                 #
# #
# --------------------------------------------------------------------------- #

.agp_temporal <- function(type, n_times = 5L, n = 10L) {
  structure(list(type = type, time_idx = rep_len(seq_len(n_times), n),
                 n_times = n_times, shared = c(TRUE, FALSE)),
            class = "tobs_temporal")
}

test_that("cover() multi-block temporal blocks declare every defaulted axis", {
  skip_if_no_auto_grid()
  idx_pos <- 1:5

  ar1 <- tulpaObs:::.cover_temporal_block(.agp_temporal("ar1"), idx_pos, list())
  expect_true(tulpa::is_auto_grid(ar1$tau_grid))
  expect_true(tulpa::is_auto_grid(ar1$rho_grid))
  expect_equal(ar1$tau_grid, c(1, 4, 16), ignore_attr = TRUE)

  # Provenance is per axis: pinning the precision leaves the correlation ours.
  mixed <- tulpaObs:::.cover_temporal_block(
    .agp_temporal("ar1"), idx_pos, list(tau.temporal.grid = c(2, 8)))
  expect_false(tulpa::is_auto_grid(mixed$tau_grid))
  expect_true(tulpa::is_auto_grid(mixed$rho_grid))

  iid <- tulpaObs:::.cover_temporal_block(.agp_temporal("iid"), idx_pos, list())
  expect_true(tulpa::is_auto_grid(iid$sigma_grid))
  expect_false(tulpa::is_auto_grid(tulpaObs:::.cover_temporal_block(
    .agp_temporal("iid"), idx_pos, list(sigma.temporal.grid = c(0.2, 0.4)))$sigma_grid))

  for (ty in c("rw1", "rw2")) {
    rw <- tulpaObs:::.cover_temporal_block(.agp_temporal(ty), idx_pos, list())
    expect_true(tulpa::is_auto_grid(rw$tau_grid), info = ty)
    expect_false(tulpa::is_auto_grid(tulpaObs:::.cover_temporal_block(
      .agp_temporal(ty), idx_pos, list(tau.temporal.grid = c(2, 8)))$tau_grid),
      info = ty)
  }
})

test_that("cover() multi-block RE block declares its defaulted sigma axis", {
  skip_if_no_auto_grid()
  re <- structure(list(group = "obs", type = "intercept", model = "iid",
                       group_idx = rep_len(1:3, 10L), n_groups = 3L),
                  class = "tobs_re")
  blk <- tulpaObs:::.cover_re_block(re, 1:5, list())
  expect_true(tulpa::is_auto_grid(blk$sigma_grid))
  expect_equal(blk$sigma_grid, exp(seq(log(0.1), log(1.5), length.out = 3)),
               ignore_attr = TRUE)
  expect_false(tulpa::is_auto_grid(
    tulpaObs:::.cover_re_block(re, 1:5, list(sigma.re.grid = c(0.2, 0.4)))$sigma_grid))
})

test_that("cover() multi-block carries the copy axis's mark through as.numeric()", {
  skip_if_no_auto_grid()
  # `.cover_build_multi_prior()` coerced the vector its caller had just
  # marked, so the copy axis reached the engine looking like a pin. Both
  # provenances round-trip. The axis is `alpha_grid` since -- the copy
  # coefficient rides the one field tulpa's `.resolve_one_copy_spec()` reads a
  # grid off, and the retired `sigma_pos_grid` was never read on any cover()
  # route.
  g <- matrix(0L, 4L, 4L); g[1, 2] <- g[2, 1] <- g[2, 3] <- g[3, 2] <-
    g[3, 4] <- g[4, 3] <- 1L
  sp <- list(type = "icar", n_spatial_units = 4L,
             adj_row_ptr = c(0L, 1L, 3L, 5L, 6L),
             adj_col_idx = c(1L, 0L, 2L, 1L, 3L, 2L),
             n_neighbors = c(1L, 2L, 2L, 1L))
  build <- function(alpha_grid) {
    tulpaObs:::.cover_build_multi_prior(
      prior_spatial = sp, spi_full = rep_len(1:4, 10L),
      spi_pos = rep_len(1:4, 5L), idx_pos = 1:5,
      temporal = NULL, re = NULL, control = list(),
      alpha_grid = alpha_grid)
  }
  auto <- build(tulpaObs:::.tobs_default_alpha_grid())
  expect_true(tulpa::is_auto_grid(auto$copy$alpha_grid))
  expect_false(tulpa::is_auto_grid(build(c(0.4, 0.8, 1.2))$copy$alpha_grid))
  # The spatial block's own axis is defaulted (and declared) in the same call.
  expect_true(tulpa::is_auto_grid(auto$prior[[1L]]$sigma_grid))
})

# --------------------------------------------------------------------------- #
# cover() arm-specific block: the 0.2-2.5 axis that never matched the engine    #
# --------------------------------------------------------------------------- #

.agp_graph <- function(n = 6L) {
  g <- matrix(0L, n, n)
  for (s in seq_len(n)) {
    if (s > 1L) g[s, s - 1L] <- 1L
    if (s < n)  g[s, s + 1L] <- 1L
  }
  g
}

.agp_block <- function(type, control) {
  g <- .agp_graph()
  tulpaObs:::.cover_armspecific_block(
    type = type, graph = g, slot = 2L, idx_active = rep_len(1:6, 12L),
    n_occ = 12L, n_pos = 12L, svc_weight = NULL, control = control,
    block_label = "test")
}

test_that("arm-specific icar tau_grid is declared when defaulted", {
  skip_if_no_auto_grid()
  # tau = 1 / sigma^2 then sort() -- both drop the marker, so this asserts the
  # re-application on the translated vector, not on the source.
  b <- .agp_block("icar", list())
  expect_true(tulpa::is_auto_grid(b$tau_grid))
  expect_false(tulpa::is_auto_grid(.agp_block("icar", list(sigma.grid = c(0.4, 1, 2)))$tau_grid))
})

test_that("arm-specific car_proper declares both its axes independently", {
  skip_if_no_auto_grid()
  b <- .agp_block("car_proper", list())
  expect_true(tulpa::is_auto_grid(b$tau_grid))
  expect_true(tulpa::is_auto_grid(b$rho_car_grid))
  # Pinning one axis leaves the other ours: provenance is per axis, not per fit.
  mixed <- .agp_block("car_proper", list(rho.car.grid = c(0.5, 0.9)))
  expect_true(tulpa::is_auto_grid(mixed$tau_grid))
  expect_false(tulpa::is_auto_grid(mixed$rho_car_grid))
})

test_that("arm-specific bym2 declares the paired (sigma, rho) axes", {
  skip_if_no_auto_grid()
  # expand.grid() drops the marker off both columns, so this is the pairing
  # equivalent of the tau case above.
  b <- .agp_block("bym2", list())
  expect_true(tulpa::is_auto_grid(b$sigma_grid))
  expect_true(tulpa::is_auto_grid(b$rho_grid))
  pinned <- .agp_block("bym2", list(sigma.grid = c(0.4, 1), rho.grid = c(0.3, 0.9)))
  expect_false(tulpa::is_auto_grid(pinned$sigma_grid))
  expect_false(tulpa::is_auto_grid(pinned$rho_grid))
})

test_that("arm-specific bym2 defaults rho to the engine's axis (#206)", {
  skip_if_no_auto_grid()
  # The nodes are READ from the engine at fit time, so this block integrates the
  # same mixing-weight axis a bym2 block reaching the registry any other way
  # does. The assertion reads the engine too: a literal here would drift out of
  # step exactly as the block's own copy did when extended the axis to cover a
  # rho near 1.
  engine_rho <- tulpa:::.nl_grid_axis("bym2_rho")
  b <- .agp_block("bym2", list())
  expect_equal(as.numeric(unique(b$rho_grid)), as.numeric(engine_rho))
  # Paired, not two separate axes: the block carries one (sigma, rho) cell per
  # combination, which is what the registry consumes.
  expect_length(b$rho_grid,
                length(unique(as.numeric(b$sigma_grid))) * length(engine_rho))
  # A pin still wins, and still reads as a pin.
  pinned <- .agp_block("bym2", list(rho.grid = c(0.3, 0.9)))
  expect_equal(as.numeric(unique(pinned$rho_grid)), c(0.3, 0.9))
})

# --------------------------------------------------------------------------- #
# EM nested-Laplace path: per-component defaults built through expand.grid()    #
# --------------------------------------------------------------------------- #

test_that("the EM-path iid RE block declares its own sigma axis", {
  skip_if_no_auto_grid()
  # Group codes are resolved per site when the re() term is built, so the spec
  # is handed in already resolved rather than round-tripped through a formula.
  re <- structure(list(group = "g", type = "intercept", model = "iid",
                       group_idx = rep_len(1:4, 12L), n_groups = 4L),
                  class = "tobs_re")
  blk <- tulpaObs:::.tobs_block_from_re(re, model = list(), n_sites = 12L,
                                        site_of_row = seq_len(12L))
  expect_true(tulpa::is_auto_grid(blk$sigma_grid))

  # A user-supplied axis on the same block stays a pin.
  re$sigma_grid <- c(0.3, 1, 3)
  blk2 <- tulpaObs:::.tobs_block_from_re(re, model = list(), n_sites = 12L,
                                         site_of_row = seq_len(12L))
  expect_false(tulpa::is_auto_grid(blk2$sigma_grid))
})

# --------------------------------------------------------------------------- #
# End to end: a defaulted axis stays recentre-able, a pinned one does not       #
# --------------------------------------------------------------------------- #

test_that("occu_cover() declines to recenter a PINNED axis but not a defaulted one", {
  skip_if_fast()
  skip_on_cran()
  skip_if_no_auto_grid()

  N <- 40L; J <- 5L
  adj <- .agp_graph(N)
  sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                             adj = adj, sigma = 1, alpha = 1, seed = 186L)
  long <- data.frame(site_id = rep(seq_len(N), each = J),
                     visit = rep(seq_len(J), times = N),
                     y = as.vector(t(sim$y)),
                     det_cov1 = sim$visit_data$det_cov1,
                     pos_cov1 = sim$visit_data$pos_cov1)
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  run <- function(ctrl) suppressWarnings(tobs(
    occurrence = ~ occ_cov1 + icar(graph = adj), data = cell_dat,
    family = occu_cover("lognormal"), detection = ~ det_cov1,
    positive = ~ pos_cov1, y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "nested_laplace", control = c(list(verbose = FALSE), ctrl)))

  # A grid the user named is a pin: the engine says so and leaves it alone.
  expect_identical(run(list(sigma.grid = c(0.3, 0.9, 2.0)))$outer_grid_recenter_declined,
                   "axis_pinned")
  # A grid we defaulted is not. Whatever the engine then decides, "the user
  # pinned it" is not among the reasons -- which is the whole point of #186.
  auto <- run(list())
  expect_false(identical(auto$outer_grid_recenter_declined, "axis_pinned"))
})
