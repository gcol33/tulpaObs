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
       family = occu_cover(response = "beta", cover_aggregate = "none"),
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
         family = occu_cover(response = "beta", cover_aggregate = "none"),
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


# The per-visit diagnostics -- pointwise log-likelihood, PPC, PIT / LOO-PIT --
# all read .occu_cover_visit_view(), which is the compact fit's stored visit rows
# and a flattening of the dense fit's padded grid. Everything downstream of it
# therefore has to agree between the two builds; that is what these three assert
# (gcol33/tulpaObs#185, where tobs_cpo() / tobs_ppc() instead read model$y /
# model$valid and errored on a compact fit with "'x' must be an array of at least
# two dimensions" while tobs_waic() worked).

test_that(".occu_cover_visit_view is identical for a dense and a compact fit", {
  skip_on_cran()
  nc  <- 8L
  adj <- .line_graph(nc)
  dd  <- .mk_occu_cover_long(nc, 3L, function(cell, ti) sample(2:6, 1), seed = 21)
  vd <- tulpaObs:::.occu_cover_visit_view(.fit_occu_cover(dd, adj, FALSE)$model)
  vc <- tulpaObs:::.occu_cover_visit_view(.fit_occu_cover(dd, adj, TRUE)$model)

  expect_identical(vc$V, vd$V)
  expect_identical(vc$V, nrow(dd))
  expect_identical(vc$site_of_visit, vd$site_of_visit)
  expect_identical(vc$y_det_visit,   vd$y_det_visit)
  expect_identical(vc$y_pos_visit,   vd$y_pos_visit)
  expect_equal(unname(vc$X_det_visit), unname(vd$X_det_visit), tolerance = 0)
  expect_equal(unname(vc$X_pos_visit), unname(vd$X_pos_visit), tolerance = 0)
  # The per-site detection summaries the PPC and the PIT read off the view are
  # the ones the dense path derived from `rowSums(y * valid)`.
  expect_identical(vc$n_valid, vd$n_valid)
  expect_identical(vc$any_det, vd$any_det)
  expect_identical(sum(vd$n_valid), nrow(dd))
  expect_identical(vd$any_det,
                   as.integer(tapply(dd$occur, dd$site_key, max)[
                     unique(dd$site_key)] > 0L))
})


test_that("compact == dense for tobs_cpo() / tobs_ppc() / PIT (#185)", {
  skip_on_cran()
  nc  <- 10L
  adj <- .line_graph(nc)
  dd  <- .mk_occu_cover_long(nc, 3L, function(cell, ti) sample(2:6, 1), seed = 22)
  fit_d <- .fit_occu_cover(dd, adj, compact = FALSE)
  fit_c <- .fit_occu_cover(dd, adj, compact = TRUE)

  set.seed(31); cd <- tobs_cpo(fit_d, n.draws = 200L)
  set.seed(31); cc <- tobs_cpo(fit_c, n.draws = 200L)
  expect_equal(cc$elpd_loo, cd$elpd_loo, tolerance = 1e-10)
  expect_equal(cc$p_loo,    cd$p_loo,    tolerance = 1e-10)
  # The LOO-PIT is the piece the compact path could not build at all.
  expect_length(cc$pit, fit_c$model$n_sites)
  expect_true(all(cc$pit >= 0 & cc$pit <= 1))
  expect_equal(cc$pit, cd$pit, tolerance = 1e-10)

  set.seed(32); pd <- tobs_ppc(fit_d, n.samples = 150L)
  set.seed(32); pc <- tobs_ppc(fit_c, n.samples = 150L)
  expect_true(all(is.finite(pc$fit.y)) && all(is.finite(pc$fit.y.rep)))
  expect_equal(pc$fit.y,     pd$fit.y,     tolerance = 1e-10)
  expect_equal(pc$fit.y.rep, pd$fit.y.rep, tolerance = 1e-10)
  expect_identical(pc$bayesian.p, pd$bayesian.p)

  # Both discrepancies, so the chi-squared branch is covered too.
  set.seed(33); qd <- tobs_ppc(fit_d, fit.stat = "chi-squared", n.samples = 100L)
  set.seed(33); qc <- tobs_ppc(fit_c, fit.stat = "chi-squared", n.samples = 100L)
  expect_equal(qc$fit.y, qd$fit.y, tolerance = 1e-10)

  set.seed(34); rd <- tobs_pit_residuals(fit_d, n.samples = 150L)
  set.seed(34); rc <- tobs_pit_residuals(fit_c, n.samples = 150L)
  expect_equal(rc, rd, tolerance = 1e-10)
})


test_that("compact == dense for the per-visit RE offsets (#211)", {
  skip_on_cran()
  nc  <- 10L
  adj <- .line_graph(nc)
  dd  <- .mk_occu_cover_long(nc, 3L, function(cell, ti) sample(2:6, 1), seed = 8)

  fit_re <- function(compact) {
    od <- tobs_data(dd, y = "occur", site = "site_key", visit = "visit",
                    type = "occurrence", occ.covs = c("cell_idx", "time.sc"),
                    det.covs = c("x1", "hab"), compact = compact)
    ocv <- suppressMessages(tobs_data(dd, y = "cover.flat", site = "site_key",
                    visit = "visit", type = "cover", compact = compact))
    ctrl <- list(engine = "joint", n.threads = 1L, n.threads.outer = 1L,
                 max.iter = 60L, sigma.grid = c(0.5, 1.5), phi.grid.pos = c(2, 10),
                 re.sigma.grid.p = c(0.2, 0.6), integration = "grid",
                 adaptive.grid = FALSE, verbose = FALSE, progress = FALSE)
    tobs(occurrence = ~ time.sc + spatial(~ 1 + time.sc || cell_idx, graph = adj),
         data = od$occ.covs,
         family = occu_cover(response = "beta", cover_aggregate = "none"),
         detection = ~ x1 + (1 | hab), positive = ~ hab,
         y = od$y, y_pos = ocv$y, visits = od$det.covs,
         method = "nested_laplace", control = ctrl)
  }
  vd <- tulpaObs:::.occu_cover_visit_view(fit_re(FALSE)$model)
  vc <- tulpaObs:::.occu_cover_visit_view(fit_re(TRUE)$model)

  # The random effect's per-visit group codes are read off the same view the
  # per-visit diagnostics read, so the compact fit's stored visit rows and the
  # dense grid's valid cells carry the same codes in the same order.
  expect_length(vc$re$p, 1L)
  expect_identical(vc$re$p[[1L]]$codes, vd$re$p[[1L]]$codes)
  expect_identical(vc$re$p[[1L]]$n_groups, vd$re$p[[1L]]$n_groups)
  expect_length(vc$re$p[[1L]]$codes, nrow(dd))

  # And the [V x S] offsets built from one set of group draws therefore agree
  # to the bit between the two layouts.
  set.seed(7)
  D  <- matrix(stats::rnorm(20L * vd$re$p[[1L]]$n_groups), 20L)
  rd <- list(list(arm = "p", var = vd$re$p[[1L]]$var, draws = D))
  od_ <- tulpaObs:::.occu_cover_re_visit_offsets(vd, rd)
  oc_ <- tulpaObs:::.occu_cover_re_visit_offsets(vc, rd)
  expect_identical(dim(od_$p), c(nrow(dd), 20L))
  expect_identical(oc_$p, od_$p)
  expect_null(oc_$pos)
})


test_that("a detected visit with a missing cover carries no cover term (#185)", {
  skip_on_cran()
  nc  <- 8L
  adj <- .line_graph(nc)
  dd  <- .mk_occu_cover_long(nc, 3L, function(cell, ti) sample(2:6, 1), seed = 23)
  # Drop the cover value (keeping the detection) at a few detected visits. The
  # likelihood gates the cover density on `detected AND finite`, so the PPC has
  # to skip those cells; scoring the NA instead poisons the whole discrepancy.
  det <- which(dd$occur == 1L)
  expect_gt(length(det), 5L)
  dd$cover.flat[det[seq_len(5L)]] <- NA_real_

  fit_d <- .fit_occu_cover(dd, adj, compact = FALSE)
  fit_c <- .fit_occu_cover(dd, adj, compact = TRUE)
  vw <- tulpaObs:::.occu_cover_visit_view(fit_c$model)
  expect_equal(sum(vw$y_det_visit == 1L & !is.finite(vw$y_pos_visit)), 5L)

  set.seed(41); pd <- tobs_ppc(fit_d, n.samples = 150L)
  set.seed(41); pc <- tobs_ppc(fit_c, n.samples = 150L)
  expect_true(all(is.finite(pd$fit.y)) && all(is.finite(pd$fit.y.rep)))
  expect_equal(pc$fit.y,     pd$fit.y,     tolerance = 1e-10)
  expect_equal(pc$fit.y.rep, pd$fit.y.rep, tolerance = 1e-10)

  set.seed(42); cd <- tobs_cpo(fit_d, n.draws = 200L)
  set.seed(42); cc <- tobs_cpo(fit_c, n.draws = 200L)
  expect_true(is.finite(cc$elpd_loo))
  expect_equal(cc$elpd_loo, cd$elpd_loo, tolerance = 1e-10)
})


test_that("a compact `by =` batch scores tobs_cpo() / tobs_ppc() per species (#185)", {
  skip_on_cran()
  skip_if_fast()
  # `tobs(by = )` on the nested-Laplace route builds COMPACT per-species arms by
  # default, so every per-species fit in a batch hit the dense-grid read. This is
  # the shape the per-species `_loo.csv` is written from.
  nc  <- 8L
  adj <- .line_graph(nc)
  dd  <- .mk_occu_cover_long(nc, 3L, function(cell, ti) sample(2:5, 1), seed = 51)
  d2  <- rbind(transform(dd, sp = "spA"), transform(dd, sp = "spB"))
  set.seed(9)
  is_b <- d2$sp == "spB"
  d2$occur[is_b] <- rbinom(nrow(dd), 1L, 0.3)
  d2$cover.flat[is_b] <- ifelse(d2$occur[is_b] == 1L, 0.4, 0)

  fit_batch <- function(compact) {
    ctrl <- list(engine = "joint", n.threads = 1L, n.threads.outer = 1L,
                 max.iter = 200L, sigma.grid = c(0.5, 1.5),
                 phi.grid.pos = c(2, 10), integration = "grid",
                 adaptive.grid = FALSE, verbose = FALSE, progress = FALSE,
                 compact = compact)
    suppressMessages(suppressWarnings(tobs(
      occurrence = ~ time.sc + spatial(~ 1 + time.sc || cell_idx, graph = adj),
      data = d2, family = occu_cover(response = "beta", cover_aggregate = "none"),
      detection = ~ x1 + hab, positive = ~ hab,
      site = "site_key", visit = "visit", response = "occur",
      y_pos = "cover.flat", occ.covs = c("cell_idx", "time.sc"),
      det.covs = c("x1", "hab"), by = "sp",
      method = "nested_laplace", control = ctrl)))
  }
  b_c <- fit_batch(TRUE)
  b_d <- fit_batch(FALSE)
  expect_s3_class(b_c, "tobs_batch")

  for (sp in c("spA", "spB")) {
    fc <- tobs_get(b_c, sp); fd <- tobs_get(b_d, sp)
    expect_true(isTRUE(fc$model$ragged))
    expect_null(fc$model$y)                  # no padded grid to read
    expect_false(isTRUE(fd$model$ragged))

    set.seed(61); cc <- tobs_cpo(fc, n.draws = 200L)
    set.seed(61); cd <- tobs_cpo(fd, n.draws = 200L)
    expect_true(is.finite(cc$elpd_loo) && all(is.finite(cc$pit)))
    expect_equal(cc$elpd_loo, cd$elpd_loo, tolerance = 1e-10)
    expect_equal(cc$pit,      cd$pit,      tolerance = 1e-10)

    set.seed(62); pc <- tobs_ppc(fc, n.samples = 150L)
    set.seed(62); pd <- tobs_ppc(fd, n.samples = 150L)
    expect_equal(pc$fit.y, pd$fit.y, tolerance = 1e-10)
    expect_identical(pc$bayesian.p, pd$bayesian.p)
  }
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
         family = occu_cover(response = "beta", cover_aggregate = "none"),
         detection = ~ x1 + hab, positive = ~ hab,
         y = od$y, y_pos = ocv$y, visits = od$det.covs,
         method = "laplace", control = list(verbose = FALSE)),
    "compact")
})


# Detection-pattern compression: within a site, all-undetected visits that share
# a detection design row enter the occupancy mixture only through
# prod (1 - p)^weight, so they collapse to one weighted row (carried in the p
# arm's n_trials). Detected visits stay individual (each has its own cover). The
# reduction is exact sufficient statistics, so a compressed fit must match the
# uncompressed fit to floating-point (option tulpaObs.compress_nodet toggles it).

test_that(".occu_cover_compress_nodet_visits groups exchangeable nodet visits", {
  site <- c(1, 1, 1, 1, 1, 2, 2, 2)
  Xp   <- rbind(c(0, 0), c(0, 0), c(1, 0), c(0, 0), c(1, 0),   # site 1
                c(2, 1), c(2, 1), c(3, 0))                     # site 2
  ydet <- c(0, 0, 0, 1, 0, 0, 0, 1)                            # detected: rows 4, 8
  cmp  <- tulpaObs:::.occu_cover_compress_nodet_visits(site, Xp, ydet)
  expect_equal(sum(cmp$weight), length(site))          # weights partition all visits
  expect_true(all(c(4L, 8L) %in% cmp$sel[cmp$weight == 1L]))  # detected stay individual
  expect_lt(length(cmp$sel), length(site))             # some nodet rows collapsed
  # every compressed row's design equals its group members' (exact grouping)
  expect_false(any(duplicated(cmp$sel)))
})


test_that("detection-pattern compression == uncompressed fit (exact)", {
  skip_on_cran()
  nc  <- 10L
  adj <- .line_graph(nc)
  # Low-cardinality detection design (rounded x1 + 3-level hab) so many
  # all-undetected visits within a site share a detection row and compress.
  dd  <- .mk_occu_cover_long(nc, 3L, function(cell, ti) sample(8:16, 1), seed = 3)
  dd$x1 <- round(dd$x1)

  old <- options(tulpaObs.compress_nodet = TRUE)
  fit_on  <- .fit_occu_cover(dd, adj, compact = TRUE)
  options(tulpaObs.compress_nodet = FALSE)
  fit_off <- .fit_occu_cover(dd, adj, compact = TRUE)
  options(old)

  # Exact sufficient statistics: identical up to floating-point reassociation.
  expect_equal(fit_on$means, fit_off$means, tolerance = 1e-8)
  expect_equal(fit_on$sds,   fit_off$sds,   tolerance = 1e-8)
  expect_identical(fit_on$N, fit_off$N)

  # Predictions and WAIC also match (seed the draws so it is algebra, not MC).
  nd <- data.frame(cell = seq_len(nc), time.sc = 0,
                   hab = factor("A", levels = c("A", "B", "C")))
  set.seed(7); po_on  <- as.data.frame(predict(fit_on,  newdata = nd,
                   type = "occurrence", nsim = 100, draws = FALSE))
  set.seed(7); po_off <- as.data.frame(predict(fit_off, newdata = nd,
                   type = "occurrence", nsim = 100, draws = FALSE))
  expect_equal(po_on$mean, po_off$mean, tolerance = 1e-8)
})
