# =============================================================================
# test-occu-cover-by.R - the `by = "<species_col>"` ergonomic batch path.
#
# `tobs(data, ..., by = "sp")` takes a long / plot-level frame, splits by the
# species column, builds each species' response onto one shared site x visit
# grid (via tobs_data()), and routes the per-species responses through the
# batched-independent driver. The contract: identical to (a) hand-building the
# per-species responses and passing the multi-response `y` list, and (b) fitting
# each species separately with single-species tobs(). The split only reorganises
# the input; it does not change the statistics.
# =============================================================================


# ---- fast-tier: argument validation + the long -> response build -------------

test_that("by = errors on an unsupported family / bad column / stray visits", {
  N <- 4L; J <- 2L
  df <- data.frame(site = rep(seq_len(N), each = J), visit = rep(seq_len(J), N),
                   sp = "a", occur = 0L, cov = 1)

  # occu() is not a per-plot-response family the species column splits.
  expect_error(
    tobs(~ cov, data = df, family = occu(), detection = ~ 1, by = "sp",
         site = "site", visit = "visit", response = "occur"),
    "wired for occu_cover")

  # `by` must name one real column.
  expect_error(
    tobs(~ cov, data = df, family = cover(), by = "nope",
         site = "site", response = "occur"),
    "must name a single column")

  # occu_cover() needs the pivot keys.
  expect_error(
    tobs(~ cov, data = df, family = occu_cover("lognormal"), detection = ~ 1,
         by = "sp", site = "site"),
    "supply")
})


test_that("tobs_data reuses one shared grid when sites / visits are supplied", {
  # The by = path pivots each species onto an external (sites, visits) set so the
  # rows align across species. A species absent at a site leaves that row NA;
  # a site value outside the set errors.
  df <- data.frame(site = c(1L, 1L, 3L), visit = c(1L, 2L, 1L),
                   y = c(1L, 0L, 1L))
  od <- tobs_data(df, y = "y", site = "site", visit = "visit",
                  sites = c(1L, 2L, 3L), visits = c(1L, 2L))
  expect_equal(dim(od$y), c(3L, 2L))
  expect_true(all(is.na(od$y[2L, ])))        # site 2 unobserved -> NA row
  expect_identical(od$y[1L, ], c(1L, 0L))
  expect_identical(od$y[3L, ], c(1L, NA_integer_))

  expect_error(
    tobs_data(df, y = "y", site = "site", visit = "visit",
              sites = c(1L, 2L)),
    "not in the supplied")
})


# ---- the equivalence gate: by == multi-response list == independent fits -----

test_that("by = (occu_cover) equals the multi-response list batch and independent fits", {
  skip_on_cran()

  N <- 20L; J <- 4L
  sim1 <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                              sigma = 0, alpha = 0, seed = 11L)
  sim2 <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                              sigma = 0, alpha = 0, seed = 22L)

  # Cell-level covariate (shared design) + per-species long records on a common
  # site x visit grid. Each species reports every plot (the visit covariates are
  # plot attributes), so the shared grid is unambiguous.
  cell_cov <- sim1$data$occ_cov1
  mk_long <- function(sim, lab) {
    data.frame(
      site  = rep(seq_len(N), each = J),
      visit = rep(seq_len(J), times = N),
      sp    = lab,
      occur = as.vector(t(sim$y)),
      cover = as.vector(t(ifelse(is.na(sim$y_pos), 0, sim$y_pos))),
      det_cov1 = sim$visit_data$det_cov1,
      occ_cov1 = rep(cell_cov, each = J)
    )
  }
  long <- rbind(mk_long(sim1, "a"), mk_long(sim2, "b"))

  ctrl <- list(verbose = FALSE, max.iter = 200L)

  # (1) the by = path.
  fit_by <- tobs(
    formula = ~ occ_cov1, data = long, family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ 1, method = "laplace",
    control = ctrl, by = "sp",
    site = "site", visit = "visit", response = "occur", y_pos = "cover",
    det.covs = "det_cov1")
  # The by = path builds visits internally; it must not also be handed `visits`.

  expect_s3_class(fit_by, "tobs_batch")
  expect_identical(fit_by$species, c("a", "b"))
  expect_identical(fit_by$n_species, 2L)

  # (2) hand-built multi-response `y` list on the same shared design.
  cell_dat <- data.frame(occ_cov1 = cell_cov)
  od1 <- tobs_data(mk_long(sim1, "a"), y = "occur", site = "site",
                   visit = "visit", det.covs = "det_cov1")
  y1 <- od1$y;  yp1 <- ifelse(is.na(sim1$y_pos), 0, sim1$y_pos)
  od2 <- tobs_data(mk_long(sim2, "b"), y = "occur", site = "site",
                   visit = "visit", det.covs = "det_cov1")
  y2 <- od2$y;  yp2 <- ifelse(is.na(sim2$y_pos), 0, sim2$y_pos)

  fit_list <- tobs(
    formula = ~ occ_cov1, data = cell_dat, family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ 1, method = "laplace", control = ctrl,
    y = list(a = y1, b = y2), y_pos = list(yp1, yp2), visits = od1$det.covs)

  expect_s3_class(fit_list, "tobs_batch")

  # (3) independent single-species fits.
  fit_one <- function(yy, ypp) tobs(
    formula = ~ occ_cov1, data = cell_dat, family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ 1, method = "laplace", control = ctrl,
    y = yy, y_pos = ypp, visits = od1$det.covs)
  ind_a <- fit_one(y1, yp1)
  ind_b <- fit_one(y2, yp2)

  # The looped backend replays the single-species pipeline, so each species'
  # point estimate / Hessian-derived SE is bit-identical across all three routes.
  expect_equal(fit_by$fits[["a"]]$means, fit_list$fits[["a"]]$means,
               tolerance = 1e-10)
  expect_equal(fit_by$fits[["b"]]$means, fit_list$fits[["b"]]$means,
               tolerance = 1e-10)
  expect_equal(fit_by$fits[["a"]]$means, ind_a$means, tolerance = 1e-10)
  expect_equal(fit_by$fits[["b"]]$means, ind_b$means, tolerance = 1e-10)
  expect_equal(fit_by$fits[["a"]]$sds,   ind_a$sds,   tolerance = 1e-10)
  expect_equal(fit_by$fits[["b"]]$sds,   ind_b$sds,   tolerance = 1e-10)

  # Genuinely two different fits.
  expect_false(isTRUE(all.equal(fit_by$fits[["a"]]$means,
                                fit_by$fits[["b"]]$means)))

  # Extractor.
  expect_identical(tobs_get(fit_by, "a"), fit_by$fits[["a"]])
})


test_that("by = (cover) equals independent single-species cover() fits", {
  skip_on_cran()

  N <- 40L
  set.seed(7)
  elev <- rnorm(N)
  mk_cover <- function(seed) {
    set.seed(seed)
    z <- rbinom(N, 1L, plogis(0.2 + 0.5 * elev))
    y <- ifelse(z == 1L, plogis(rnorm(N, -0.5 + 0.3 * elev, 0.4)), 0)
    pmin(pmax(y, 0), 0.999)
  }
  ya <- mk_cover(1L); yb <- mk_cover(2L)

  long <- rbind(
    data.frame(site = seq_len(N), sp = "a", cover = ya, elev = elev),
    data.frame(site = seq_len(N), sp = "b", cover = yb, elev = elev))

  ctrl <- list(verbose = FALSE)

  fit_by <- tobs(
    formula = ~ elev, data = long, family = cover("beta"),
    method = "laplace", control = ctrl, by = "sp",
    site = "site", response = "cover")

  expect_s3_class(fit_by, "tobs_batch")
  expect_identical(fit_by$species, c("a", "b"))

  cell_dat <- data.frame(elev = elev)
  ind_a <- tobs(~ elev, data = cell_dat, family = cover("beta"),
                y = ya, method = "laplace", control = ctrl)
  ind_b <- tobs(~ elev, data = cell_dat, family = cover("beta"),
                y = yb, method = "laplace", control = ctrl)

  # Cover fits carry per-arm coefficient vectors (beta_occ / beta_pos) + SEs.
  expect_equal(fit_by$fits[["a"]]$beta_occ, ind_a$beta_occ, tolerance = 1e-10)
  expect_equal(fit_by$fits[["a"]]$beta_pos, ind_a$beta_pos, tolerance = 1e-10)
  expect_equal(fit_by$fits[["a"]]$se_occ,   ind_a$se_occ,   tolerance = 1e-10)
  expect_equal(fit_by$fits[["b"]]$beta_occ, ind_b$beta_occ, tolerance = 1e-10)
  expect_equal(fit_by$fits[["b"]]$beta_pos, ind_b$beta_pos, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(fit_by$fits[["a"]]$beta_occ,
                                fit_by$fits[["b"]]$beta_occ)))
})


# ---- the by= looped backend is COMPACT on the nested-Laplace route ----------
#
# gcol33/tulpaObs#107: the default looped by= backend builds each species' arms
# (and the shared visit grid) compactly when the route is nested_laplace, so it
# no longer allocates B dense [n_sites x max_visits] response matrices plus a
# dense visit grid -- the padded-grid memory the single-fit compact path already
# avoids. The fit must still equal an independent dense single-species fit to the
# compact-vs-dense tolerance (the same invariant as test-occu-cover-compact.R).

test_that("by = (occu_cover, nested_laplace) is compact + equals dense single fits", {
  skip_on_cran()
  skip_if_fast()

  N <- 20L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }

  # Beta cover (a [0, 1] proportion -- the type = "cover" arm), the production
  # family. Each site is one graph cell (site_id == cell == graph node).
  sim1 <- simulate_occu_cover(N = N, J = J, positive = "beta", adj = adj,
                              sigma = 0.8, seed = 101L)
  sim2 <- simulate_occu_cover(N = N, J = J, positive = "beta", adj = adj,
                              sigma = 0.8, seed = 202L)

  # The by= contract: visit covariates are plot attributes shared across species
  # (here from sim1), so the species differ only in occur / cover. Long records
  # on the common site x visit grid.
  base <- data.frame(
    site_id  = rep(seq_len(N), each = J),
    visit    = rep(seq_len(J), times = N),
    occ_cov1 = rep(sim1$data$occ_cov1, each = J),
    det_cov1 = sim1$visit_data$det_cov1)
  mk <- function(sim, lab) {
    d <- base; d$sp <- lab
    d$occur <- as.vector(t(sim$y))
    d$cover <- as.vector(t(ifelse(is.na(sim$y_pos), 0, sim$y_pos)))
    d
  }
  long <- rbind(mk(sim1, "a"), mk(sim2, "b"))

  ctrl <- list(engine = "joint", n.threads = 1L, n.threads.outer = 1L,
               max.iter = 200L, sigma.grid = c(0.5, 1.5), phi.grid.pos = c(2, 10),
               integration = "grid", adaptive.grid = FALSE, verbose = FALSE,
               progress = FALSE)

  fit_by <- suppressWarnings(suppressMessages(tobs(
    occurrence = ~ occ_cov1 + spatial(~ 1 || site_id, graph = adj), data = long,
    family = occu_cover(positive = "beta", cover_aggregate = "none"),
    detection = ~ det_cov1, positive = ~ 1,
    method = "nested_laplace", control = ctrl, by = "sp",
    site = "site_id", visit = "visit", response = "occur", y_pos = "cover",
    det.covs = "det_cov1")))

  expect_s3_class(fit_by, "tobs_batch")
  expect_identical(fit_by$backend, "looped")
  # Compact build: the per-species observation count is the full plot set
  # (n_sites x J), the same N the single-fit compact path reports.
  expect_equal(fit_by$fits[["a"]]$N, N * J)

  # Independent dense single-species fits at the same shared design = the oracle.
  od1 <- tobs_data(mk(sim1, "a"), y = "occur", site = "site_id", visit = "visit",
                   occ.covs = "occ_cov1", det.covs = "det_cov1")
  cell_dat <- data.frame(site_id = seq_len(N), occ_cov1 = sim1$data$occ_cov1)
  yp1 <- ifelse(is.na(sim1$y_pos), 0, sim1$y_pos)
  yp2 <- ifelse(is.na(sim2$y_pos), 0, sim2$y_pos)
  fit_one <- function(yy, ypp) suppressWarnings(suppressMessages(tobs(
    occurrence = ~ occ_cov1 + spatial(~ 1 || site_id, graph = adj), data = cell_dat,
    family = occu_cover(positive = "beta", cover_aggregate = "none"),
    detection = ~ det_cov1, positive = ~ 1,
    y = yy, y_pos = ypp, visits = od1$det.covs,
    method = "nested_laplace", control = ctrl)))
  ind_a <- fit_one(sim1$y, yp1)
  ind_b <- fit_one(sim2$y, yp2)

  # Compact looped by= == dense independent fit, per species (float-summation
  # tolerance, the compact-vs-dense invariant).
  expect_equal(fit_by$fits[["a"]]$means, ind_a$means, tolerance = 1e-7)
  expect_equal(fit_by$fits[["b"]]$means, ind_b$means, tolerance = 1e-7)
  expect_equal(fit_by$fits[["a"]]$sds,   ind_a$sds,   tolerance = 1e-7)
  expect_false(isTRUE(all.equal(fit_by$fits[["a"]]$means,
                                fit_by$fits[["b"]]$means)))
})


# ---- by= with a positive (unbounded) cover arm ------------------------------
#
# gcol33/tulpaObs#107: the by= cover arm picks its tobs_data() storage type from
# the family's positive distribution, so a lognormal cover that exceeds 1 routes
# through type = "positive" instead of being rejected by the [0, 1] check that
# type = "cover" enforces (the regime the old by= path could not handle).

test_that("by = (occu_cover, lognormal) accepts cover > 1 via the positive type", {
  skip_on_cran()

  N <- 20L; J <- 4L
  sim1 <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                              sigma = 0, alpha = 0, seed = 11L)
  sim2 <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                              sigma = 0, alpha = 0, seed = 22L)
  cell_cov <- sim1$data$occ_cov1
  mk_long <- function(sim, lab) data.frame(
    site = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    sp = lab, occur = as.vector(t(sim$y)),
    cover = as.vector(t(ifelse(is.na(sim$y_pos), 0, sim$y_pos))),
    det_cov1 = sim$visit_data$det_cov1, occ_cov1 = rep(cell_cov, each = J))
  long <- rbind(mk_long(sim1, "a"), mk_long(sim2, "b"))
  long$cover[long$occur == 1L] <- long$cover[long$occur == 1L] + 1.5
  expect_gt(max(long$cover), 1)   # type = "cover" would reject this

  fit_by <- tobs(
    formula = ~ occ_cov1, data = long, family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ 1, method = "laplace",
    control = list(verbose = FALSE, max.iter = 200L), by = "sp",
    site = "site", visit = "visit", response = "occur", y_pos = "cover",
    det.covs = "det_cov1")

  expect_s3_class(fit_by, "tobs_batch")
  expect_identical(fit_by$species, c("a", "b"))
  expect_true(all(is.finite(fit_by$fits[["a"]]$means)))
  expect_true(all(is.finite(fit_by$fits[["b"]]$means)))
})
