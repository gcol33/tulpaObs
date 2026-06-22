# Compact (ragged) occu_cover input: one row per valid visit instead of a padded
# [n_sites x max_visits] grid, so there is no per-site visit cap. The joint
# nested-Laplace engine reads one valid visit at a time, so the compact build
# must produce a fit byte-identical to the dense build on the SAME data -- that
# equivalence is the invariant these tests assert (a divergence is a bug, not an
# "improvement", because both paths feed the engine the same `responses` list).
# The cap-removal test then fits a site whose visit count would force a wide
# dense grid, confirming the compact path needs only O(observations) memory.

# A small grid: nc cells on a line graph, n_int intervals, a handful of visits
# per site (well under any cap). Returns the long plot-level frame.
.mk_occu_cover_long <- function(nc, n_int, visits_per_site, seed = 1) {
  set.seed(seed)
  rows <- list()
  for (cell in seq_len(nc)) for (ti in seq_len(n_int)) {
    J   <- visits_per_site(cell, ti)
    tsc <- ti - (n_int + 1) / 2
    psi <- plogis(-0.3 + 0.25 * tsc + 0.15 * (cell - nc / 2))
    occ_state <- rbinom(1, 1, psi)
    x1  <- rnorm(J)
    hab <- factor(sample(c("A", "B", "C"), J, replace = TRUE), levels = c("A", "B", "C"))
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

.line_graph <- function(nc) {
  adj <- matrix(0L, nc, nc)
  for (i in seq_len(nc - 1L)) { adj[i, i + 1L] <- 1L; adj[i + 1L, i] <- 1L }
  adj
}

.fit_occu_cover <- function(dd, adj, compact) {
  od <- tobs_data(dd, y = "occur", site = "site_key", visit = "visit",
                  type = "occurrence", occ.covs = c("cell_idx", "time.sc"),
                  det.covs = c("x1", "hab"), compact = compact)
  ocv <- suppressMessages(tobs_data(dd, y = "cover.flat", site = "site_key",
                  visit = "visit", type = "cover", compact = compact))
  ctrl <- list(engine = "joint", n.threads = 1L, n.threads.outer = 1L,
               max.iter = 200L, sigma.grid = c(0.5, 1.5), phi.grid.pos = c(2, 10),
               integration = "grid", adaptive.grid = FALSE, verbose = FALSE,
               progress = FALSE)
  tobs(occurrence = ~ time.sc + spatial(~ 1 + time.sc || cell_idx, graph = adj),
       data = od$occ.covs,
       family = occu_cover(positive = "beta", cover_aggregate = "none"),
       detection = ~ x1 + hab, positive = ~ hab,
       y = od$y, y_pos = ocv$y, visits = od$det.covs,
       method = "nested_laplace", control = ctrl)
}


test_that("tobs_data(compact = TRUE) builds an aligned ragged carrier", {
  dd <- .mk_occu_cover_long(8L, 3L, function(cell, ti) 4L, seed = 2)
  od  <- tobs_data(dd, y = "occur", site = "site_key", visit = "visit",
                   type = "occurrence", occ.covs = c("cell_idx", "time.sc"),
                   det.covs = c("x1", "hab"), compact = TRUE)
  ocv <- suppressMessages(tobs_data(dd, y = "cover.flat", site = "site_key",
                   visit = "visit", type = "cover", compact = TRUE))
  expect_s3_class(od$y, "tobs_ragged")
  expect_equal(od$y$n_visits, nrow(dd))
  expect_equal(od$y$n_sites, length(unique(dd$site_key)))
  # occurrence and cover ragged carriers align row-for-row.
  expect_identical(od$y$site, ocv$y$site)
  expect_identical(od$y$visit, ocv$y$visit)
  # factor detection covariate keeps its level set for the downstream design.
  expect_true(isTRUE(attr(od$det.covs$hab, "tobs_factor")))
  expect_identical(attr(od$det.covs$hab, "tobs_levels"), c("A", "B", "C"))
  expect_length(od$det.covs$x1, nrow(dd))
})


test_that("compact == dense: byte-identical fit on uncapped data", {
  skip_on_cran()
  nc <- 12L
  adj <- .line_graph(nc)
  dd  <- .mk_occu_cover_long(nc, 3L, function(cell, ti) sample(2:6, 1), seed = 3)

  fit_d <- .fit_occu_cover(dd, adj, compact = FALSE)
  fit_c <- .fit_occu_cover(dd, adj, compact = TRUE)

  # The engine sees the same responses list either way, so the posterior is
  # identical, not merely close.
  expect_equal(fit_c$means, fit_d$means, tolerance = 1e-8)
  expect_equal(fit_c$sds,   fit_d$sds,   tolerance = 1e-8)
  expect_identical(fit_c$N, fit_d$N)

  # Downstream predictions and WAIC are the same (seed the draws so the
  # comparison is on the algebra, not Monte-Carlo noise).
  nd <- data.frame(cell = seq_len(nc), time.sc = 0,
                   hab = factor("A", levels = c("A", "B", "C")))
  set.seed(7); po_d <- as.data.frame(predict(fit_d, newdata = nd,
                   type = "occurrence", nsim = 100, draws = FALSE))
  set.seed(7); po_c <- as.data.frame(predict(fit_c, newdata = nd,
                   type = "occurrence", nsim = 100, draws = FALSE))
  expect_equal(po_c$mean, po_d$mean, tolerance = 1e-8)

  set.seed(11); wd <- tobs_waic(fit_d)
  set.seed(11); wc <- tobs_waic(fit_c)
  expect_equal(wc$waic,    wd$waic,    tolerance = 1e-8)
  expect_equal(wc$se_waic, wd$se_waic, tolerance = 1e-8)
  expect_equal(wc$n_obs,   wd$n_obs)   # = the site count, shared across species
})


test_that("compact removes the per-site visit cap", {
  skip_on_cran()
  nc <- 10L
  adj <- .line_graph(nc)
  # One site holds 4000 visits: the dense grid would be n_sites x 4000, the
  # compact build only carries the ~4000 + few-hundred actual observations.
  dd <- .mk_occu_cover_long(nc, 3L,
          function(cell, ti) if (cell == 1L && ti == 1L) 4000L else sample(2:5, 1),
          seed = 4)
  expect_gt(max(table(dd$site_key)), 1000)

  fit <- .fit_occu_cover(dd, adj, compact = TRUE)
  expect_equal(fit$N, nrow(dd))
  expect_true(all(is.finite(fit$means)))
})


test_that("compact == dense with an observation-arm random effect (reHab)", {
  skip_on_cran()
  nc  <- 12L
  adj <- .line_graph(nc)
  dd  <- .mk_occu_cover_long(nc, 3L, function(cell, ti) sample(2:6, 1), seed = 8)

  fit_re <- function(compact) {
    od <- tobs_data(dd, y = "occur", site = "site_key", visit = "visit",
                    type = "occurrence", occ.covs = c("cell_idx", "time.sc"),
                    det.covs = c("x1", "hab"), compact = compact)
    ocv <- suppressMessages(tobs_data(dd, y = "cover.flat", site = "site_key",
                    visit = "visit", type = "cover", compact = compact))
    ctrl <- list(engine = "joint", n.threads = 1L, n.threads.outer = 1L,
                 max.iter = 200L, sigma.grid = c(0.5, 1.5), phi.grid.pos = c(2, 10),
                 re.sigma.grid.p = c(0.2, 0.6, 1.5),
                 integration = "grid", adaptive.grid = FALSE, verbose = FALSE,
                 progress = FALSE)
    # habitat as a random detection intercept (the reHab variant); cover arm keeps
    # the fixed habitat term.
    tobs(occurrence = ~ time.sc + spatial(~ 1 + time.sc || cell_idx, graph = adj),
         data = od$occ.covs,
         family = occu_cover(positive = "beta", cover_aggregate = "none"),
         detection = ~ x1 + (1 | hab), positive = ~ hab,
         y = od$y, y_pos = ocv$y, visits = od$det.covs,
         method = "nested_laplace", control = ctrl)
  }

  fit_d <- fit_re(FALSE)
  fit_c <- fit_re(TRUE)
  # The detection RE rides the compacted visit rows exactly as it rides the dense
  # grid's valid cells, so the fit is the same.
  expect_equal(fit_c$means, fit_d$means, tolerance = 1e-8)
  expect_equal(fit_c$sds,   fit_d$sds,   tolerance = 1e-8)
  rd <- as.data.frame(ranef(fit_d)); rc <- as.data.frame(ranef(fit_c))
  expect_equal(dim(rc), dim(rd))
  num <- vapply(rd, is.numeric, logical(1))
  expect_equal(unname(as.matrix(rc[num])), unname(as.matrix(rd[num])), tolerance = 1e-8)
})


test_that("WAIC draw-chunking is exact (chunk size does not change the result)", {
  skip_on_cran()
  nc  <- 10L
  adj <- .line_graph(nc)
  dd  <- .mk_occu_cover_long(nc, 3L, function(cell, ti) sample(2:6, 1), seed = 9)
  fit <- .fit_occu_cover(dd, adj, compact = TRUE)

  set.seed(3)
  c0 <- .tobs_occu_cover_components(fit, 200L)
  core <- function(ch) .occu_cover_ploglik_core(
    fit$model, c0$b_occ, c0$b_det, c0$b_pos, c0$disp,
    c0$field_occ, c0$field_pos, chunk = ch)
  ll1 <- core(1L)                       # one draw per block
  llh <- core(37L)                      # an awkward block size that does not divide S
  llS <- core(nrow(c0$b_occ))           # a single block (the old, unchunked path)
  llA <- core(NULL)                     # the memory-adaptive default
  expect_identical(ll1, llS)
  expect_identical(llh, llS)
  expect_identical(llA, llS)

  # And the criterion built on it is unchanged.
  set.seed(11); wa <- tobs_waic(fit)
  expect_true(is.finite(wa$waic) && is.finite(wa$se_waic))
})


test_that("compact input is gated to the joint nested-Laplace path", {
  dd <- .mk_occu_cover_long(6L, 2L, function(cell, ti) 4L, seed = 5)
  od  <- tobs_data(dd, y = "occur", site = "site_key", visit = "visit",
                   type = "occurrence", occ.covs = c("cell_idx", "time.sc"),
                   det.covs = c("x1", "hab"), compact = TRUE)
  ocv <- suppressMessages(tobs_data(dd, y = "cover.flat", site = "site_key",
                   visit = "visit", type = "cover", compact = TRUE))
  # No spatial term on the occurrence formula -> the non-spatial laplace path,
  # which reads the dense grid. Compact must error clearly, not silently rebuild.
  expect_error(
    tobs(occurrence = ~ time.sc, data = od$occ.covs,
         family = occu_cover(positive = "beta", cover_aggregate = "none"),
         detection = ~ x1 + hab, positive = ~ hab,
         y = od$y, y_pos = ocv$y, visits = od$det.covs,
         method = "laplace", control = list(verbose = FALSE)),
    "compact")
})
