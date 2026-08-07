# Declare a defaulted outer-grid axis as OURS with tulpa::auto_grid()
# (gcol33/tulpaObs#186, consumer of gcol33/tulpa#293).
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
  # Marking must not change the values -- clause 3 of the engine's provenance
  # test still recognises the sigma axis by equality with its own default.
  expect_equal(as.numeric(tulpaObs:::.tobs_default_sigma_grid()),
               exp(seq(log(0.1), log(3), length.out = 5)))
})

test_that(".tobs_mark_auto marks only when the caller defaulted", {
  skip_if_no_auto_grid()
  expect_true(tulpa::is_auto_grid(tulpaObs:::.tobs_mark_auto(c(0.5, 1, 2), TRUE)))
  expect_false(tulpa::is_auto_grid(tulpaObs:::.tobs_mark_auto(c(0.5, 1, 2), FALSE)))
  # Values pass through untouched either way.
  expect_equal(as.numeric(tulpaObs:::.tobs_mark_auto(c(0.5, 1, 2), TRUE)),
               c(0.5, 1, 2))
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
