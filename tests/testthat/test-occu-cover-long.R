# Single occu_cover() fit from a long / plot-level frame.
#
# tobs() now accepts the same long-frame contract for ONE species that the by=
# batch path accepts for many: pass site / visit / response / y_pos / det.covs
# with a long `data`, and the paired occurrence / cover arms + the site-level
# design are built internally (via .occu_cover_arms_from_long(), the single-fit
# counterpart of the by= builder). The invariant these tests assert: that path
# produces a fit byte-identical to the hand-built route (two tobs_data() calls +
# y / y_pos / visits / data passed explicitly), because both feed the engine the
# same `responses`. occu_cover_inputs() exposes the same builder for inspection.

# A small line-graph grid: nc cells, n_int intervals, a few visits per site.
# Returns the long plot-level frame (one row per cell-interval-plot).
.mk_oc_long_eq <- function(nc, n_int, visits_per_site, seed = 1) {
  set.seed(seed)
  rows <- list()
  for (cell in seq_len(nc)) for (ti in seq_len(n_int)) {
    J   <- visits_per_site(cell, ti)
    tsc <- ti - (n_int + 1) / 2
    psi <- plogis(-0.3 + 0.25 * tsc + 0.15 * (cell - nc / 2))
    occ_state <- rbinom(1, 1, psi)
    x1  <- rnorm(J)
    hab <- factor(sample(c("A", "B", "C"), J, replace = TRUE),
                  levels = c("A", "B", "C"))
    p   <- plogis(0.2 + 0.4 * x1 + ifelse(hab == "B", 0.3, ifelse(hab == "C", -0.2, 0)))
    occur <- if (occ_state == 1) rbinom(J, 1, p) else rep(0L, J)
    mu  <- plogis(-0.5 + ifelse(hab == "B", 0.4, 0))
    cover <- ifelse(occur == 1,
                    pmin(pmax(rbeta(J, mu * 6, (1 - mu) * 6), 1e-4), 1 - 1e-4), 0)
    rows[[length(rows) + 1L]] <- data.frame(
      cell_idx = cell, yr.interval = ti, site_key = paste(cell, ti, sep = "_"),
      visit = seq_len(J), time.sc = tsc, x1 = x1, hab = hab,
      occur = as.integer(occur), cover.flat = cover, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

.line_graph_eq <- function(nc) {
  adj <- matrix(0L, nc, nc)
  for (i in seq_len(nc - 1L)) { adj[i, i + 1L] <- 1L; adj[i + 1L, i] <- 1L }
  adj
}


# ---- fast tier: the builder's shape + alignment, no fit ---------------------

test_that("occu_cover_inputs() builds aligned arms + a per-site design", {
  dd <- .mk_oc_long_eq(6L, 2L, function(cell, ti) 4L, seed = 2)
  nc <- length(unique(dd$cell_idx))
  n_sites <- length(unique(dd$site_key))

  res <- occu_cover_inputs(dd, site = "site_key", visit = "visit",
                           response = "occur", y_pos = "cover.flat",
                           det.covs = c("x1", "hab"), compact = TRUE)

  # Compact arms: ragged carriers, aligned occurrence / cover row order.
  expect_s3_class(res$y, "tobs_ragged")
  expect_s3_class(res$y_pos, "tobs_ragged")
  expect_identical(res$y$site, res$y_pos$site)
  expect_identical(res$y$visit, res$y_pos$visit)
  expect_equal(res$n_visits, nrow(dd))
  expect_equal(res$n_sites, n_sites)

  # Visit-level design: one length-V vector per det.cov, factor level set kept.
  expect_named(res$visits, c("x1", "hab"))
  expect_length(res$visits$x1, nrow(dd))
  expect_identical(attr(res$visits$hab, "tobs_levels"), c("A", "B", "C"))

  # Site-level design: first row per site (default = all columns).
  expect_equal(nrow(res$site_data), n_sites)
  expect_true(all(c("cell_idx", "time.sc") %in% names(res$site_data)))
  expect_true(isTRUE(res$compact))

  # occ.covs subsets the site-level frame; dense build matches the same dims.
  sub <- occu_cover_inputs(dd, site = "site_key", visit = "visit",
                           response = "occur", y_pos = "cover.flat",
                           occ.covs = c("cell_idx", "time.sc"),
                           det.covs = c("x1", "hab"), compact = FALSE)
  expect_identical(names(sub$site_data), c("cell_idx", "time.sc"))
  expect_equal(dim(sub$y), c(n_sites, max(table(dd$site_key))))
})


test_that("occu_cover_inputs() / long-frame tobs() reject bad keys", {
  dd <- .mk_oc_long_eq(4L, 2L, function(cell, ti) 3L, seed = 5)

  expect_error(
    occu_cover_inputs(dd, site = "site_key", visit = "visit",
                      response = "nope", y_pos = "cover.flat"),
    "not found")
  expect_error(
    occu_cover_inputs(dd, site = "site_key", visit = "visit",
                      response = "occur", y_pos = NULL),
    "single column name")

  # response= signals long-frame mode; supplying y too is contradictory.
  adj <- .line_graph_eq(4L)
  expect_error(
    tobs(occurrence = ~ time.sc + spatial(~ 1 || cell_idx, graph = adj),
         data = dd, family = occu_cover(response = "beta", cover_aggregate = "none"),
         detection = ~ x1, response = "occur", y_pos = "cover.flat",
         site = "site_key", visit = "visit", y = matrix(0, 1, 1)),
    "OR a pre-built")
})


# ---- equivalence gate: long-frame == hand-built (od/ocv) route --------------

test_that("single long-frame occu_cover fit == the hand-built tobs_data route", {
  skip_on_cran()
  nc  <- 12L
  adj <- .line_graph_eq(nc)
  dd  <- .mk_oc_long_eq(nc, 3L, function(cell, ti) sample(2:6, 1), seed = 3)

  ctrl <- list(engine = "joint", n.threads = 1L, n.threads.outer = 1L,
               max.iter = 200L, sigma.grid = c(0.5, 1.5), phi.grid.pos = c(2, 10),
               integration = "grid", adaptive.grid = FALSE, verbose = FALSE,
               progress = FALSE)

  # (1) hand-built mi: the construction users currently write by hand.
  od  <- tobs_data(dd, y = "occur", site = "site_key", visit = "visit",
                   type = "occurrence", occ.covs = c("cell_idx", "time.sc"),
                   det.covs = c("x1", "hab"), compact = TRUE)
  ocv <- suppressMessages(tobs_data(dd, y = "cover.flat", site = "site_key",
                   visit = "visit", type = "cover", compact = TRUE))
  fit_mi <- tobs(
    occurrence = ~ time.sc + spatial(~ 1 + time.sc || cell_idx, graph = adj),
    data = od$occ.covs,
    family = occu_cover(response = "beta", cover_aggregate = "none"),
    detection = ~ x1 + hab, positive = ~ hab,
    y = od$y, y_pos = ocv$y, visits = od$det.covs,
    method = "nested_laplace", control = ctrl)

  # (2) the new single-fit long-frame path (compact defaults on for nested_laplace).
  fit_long <- tobs(
    occurrence = ~ time.sc + spatial(~ 1 + time.sc || cell_idx, graph = adj),
    data = dd,
    family = occu_cover(response = "beta", cover_aggregate = "none"),
    detection = ~ x1 + hab, positive = ~ hab,
    site = "site_key", visit = "visit", response = "occur", y_pos = "cover.flat",
    det.covs = c("x1", "hab"),
    method = "nested_laplace", control = ctrl)

  # Same responses -> identical posterior, not merely close.
  expect_equal(fit_long$means, fit_mi$means, tolerance = 1e-8)
  expect_equal(fit_long$sds,   fit_mi$sds,   tolerance = 1e-8)
  expect_identical(fit_long$N, fit_mi$N)

  # Downstream predictions match (seed the draws so the comparison is the algebra).
  nd <- data.frame(cell = seq_len(nc), time.sc = 0,
                   hab = factor("A", levels = c("A", "B", "C")))
  set.seed(7); pl <- as.data.frame(predict(fit_long, newdata = nd,
                   type = "occurrence", nsim = 100, draws = FALSE))
  set.seed(7); pm <- as.data.frame(predict(fit_mi, newdata = nd,
                   type = "occurrence", nsim = 100, draws = FALSE))
  expect_equal(pl$mean, pm$mean, tolerance = 1e-8)
})


# ---- positive (unbounded) cover arm: lognormal / gamma ----------------------
#
# the cover arm's tobs_data() storage type follows the positive
# distribution -- beta is a [0, 1] proportion ("cover"), lognormal / gamma
# are positive reals ("positive"). The long-frame builder picks it from
# `positive = `, so a lognormal cover that exceeds 1 round-trips instead of
# being rejected by the [0, 1] check.

test_that("occu_cover_inputs() picks the cover type from `positive`", {
  dd <- .mk_oc_long_eq(6L, 2L, function(cell, ti) 4L, seed = 7)
  # Inflate cover above 1 where present (a lognormal cover is unbounded).
  dd$cover.flat[dd$occur == 1L] <- dd$cover.flat[dd$occur == 1L] * 5 + 1.2

  # Beta (default) is a proportion -> the > 1 values are rejected.
  expect_error(
    occu_cover_inputs(dd, site = "site_key", visit = "visit",
                      response = "occur", y_pos = "cover.flat"),
    "\\[0, 1\\]")

  # positive = "lognormal" stores the same column as a positive real.
  res <- occu_cover_inputs(dd, site = "site_key", visit = "visit",
                           response = "occur", y_pos = "cover.flat",
                           det.covs = c("x1", "hab"),
                           positive = "lognormal", compact = TRUE)
  expect_s3_class(res$y_pos, "tobs_ragged")
  expect_true(any(res$y_pos$values > 1, na.rm = TRUE))
})


test_that("single long-frame lognormal fit == the hand-built positive route", {
  skip_on_cran()
  N <- 24L; J <- 4L
  sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                             sigma = 0, alpha = 0, seed = 5L)
  long <- data.frame(
    site     = rep(seq_len(N), each = J),
    visit    = rep(seq_len(J), times = N),
    occur    = as.vector(t(sim$y)),
    cover    = as.vector(t(ifelse(is.na(sim$y_pos), 0, sim$y_pos))),
    det_cov1 = sim$visit_data$det_cov1,
    occ_cov1 = rep(sim$data$occ_cov1, each = J))
  # Force the cover arm above 1, the regime type = "cover" would reject.
  long$cover[long$occur == 1L] <- long$cover[long$occur == 1L] + 1.5
  expect_gt(max(long$cover), 1)

  ctrl <- list(verbose = FALSE, max.iter = 200L)

  # Hand-built positive arm (the construction users write today).
  od <- tobs_data(long, y = "occur", site = "site", visit = "visit",
                  type = "occurrence", occ.covs = "occ_cov1", det.covs = "det_cov1")
  op <- suppressMessages(tobs_data(long, y = "cover", site = "site",
                  visit = "visit", type = "positive"))
  yp <- op$y; yp[is.na(yp)] <- 0
  fit_mi <- tobs(
    occurrence = ~ occ_cov1, data = od$occ.covs, family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ 1,
    y = od$y, y_pos = yp, visits = od$det.covs,
    method = "laplace", control = ctrl)

  # The single-fit long-frame path, which picks type = "positive" from the family.
  fit_long <- tobs(
    occurrence = ~ occ_cov1, data = long, family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ 1,
    site = "site", visit = "visit", response = "occur", y_pos = "cover",
    det.covs = "det_cov1", method = "laplace", control = ctrl)

  expect_equal(fit_long$means, fit_mi$means, tolerance = 1e-8)
  expect_equal(fit_long$sds,   fit_mi$sds,   tolerance = 1e-8)
})
