# Work-stealing parallel coupled-cell scatter ( flatten).
#
# When the outer nested-Laplace grid is UNDER-SATURATED (fewer active grid cells
# than team threads -- the tail of any grid, small grids, or a many-core server)
# the coupled-cell scatter (the occu_cover hot loop) splits its per-cell loop
# across the idle team threads as stealable tasks, accumulating into per-chunk
# partial Hessians reduced in a fixed chunk order. That reduce must be:
#   (1) run-to-run REPRODUCIBLE despite nondeterministic task stealing, and
#   (2) numerically equal to the serial scatter up to summation-order (~1e-6).
# A small outer grid solved on more outer threads than grid cells forces the
# parallel path on every cell, so this actually exercises it (not just the tail).

.mk_pc_long <- function(nc, n_int, seed = 5) {
  set.seed(seed); rows <- list()
  for (cell in seq_len(nc)) for (ti in seq_len(n_int)) {
    J <- sample(3:6, 1); tsc <- ti - (n_int + 1) / 2
    psi <- plogis(-0.3 + 0.25 * tsc + 0.15 * (cell - nc / 2))
    occ_state <- rbinom(1, 1, psi); x1 <- rnorm(J)
    hab <- factor(sample(c("A", "B", "C"), J, TRUE), levels = c("A", "B", "C"))
    p <- plogis(0.2 + 0.4 * x1 + ifelse(hab == "B", 0.3, ifelse(hab == "C", -0.2, 0)))
    occur <- if (occ_state == 1) rbinom(J, 1, p) else rep(0L, J)
    mu <- plogis(-0.5 + ifelse(hab == "B", 0.4, 0))
    cover <- ifelse(occur == 1, pmin(pmax(rbeta(J, mu * 6, (1 - mu) * 6), 1e-4), 1 - 1e-4), 0)
    rows[[length(rows) + 1L]] <- data.frame(
      cell_idx = cell, yr.interval = ti, site_key = paste(cell, ti, sep = "_"),
      visit = seq_len(J), time.sc = tsc, x1 = x1, hab = hab,
      occur = as.integer(occur), cover.flat = cover, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

.fit_pc <- function(dd, adj, n_out) {
  od  <- tobs_data(dd, y = "occur", site = "site_key", visit = "visit",
                   type = "occurrence", occ.covs = c("cell_idx", "time.sc"),
                   det.covs = c("x1", "hab"), compact = TRUE)
  ocv <- suppressMessages(tobs_data(dd, y = "cover.flat", site = "site_key",
                   visit = "visit", type = "cover", compact = TRUE))
  # A 4-cell outer grid (2 sigma x 2 phi) solved on n_out > 4 outer threads is
  # under-saturated on every cell, so the coupled-cell scatter runs parallel.
  ctrl <- list(engine = "joint", n.threads = n_out, n.threads.outer = n_out,
               max.iter = 200L, sigma.grid = c(0.6, 1.6), phi.grid.pos = c(3, 9),
               integration = "grid", adaptive.grid = FALSE, verbose = FALSE,
               progress = FALSE)
  tobs(occurrence = ~ time.sc + spatial(~ 1 + time.sc || cell_idx, graph = adj),
       data = od$occ.covs,
       family = occu_cover(response = "beta", cover_aggregate = "none"),
       detection = ~ x1 + hab, positive = ~ hab,
       y = od$y, y_pos = ocv$y, visits = od$det.covs,
       method = "nested_laplace", control = ctrl)
}

test_that("parallel coupled-cell scatter is reproducible and matches serial", {
  skip_on_cran()
  nc <- 120L
  dd  <- .mk_pc_long(nc, 3L, seed = 5)
  adj <- chain_adj(nc)
  n_out <- max(2L, min(8L, parallel::detectCores()))

  fit_serial <- .fit_pc(dd, adj, 1L)      # serial reference (no coupling chunking)
  fit_par_a  <- .fit_pc(dd, adj, n_out)   # under-saturated -> parallel coupling
  fit_par_b  <- .fit_pc(dd, adj, n_out)   # again, for reproducibility

  # (1) Deterministic per-chunk reduce -> bit-identical across runs despite
  #     nondeterministic task stealing.
  expect_identical(fit_par_a$means, fit_par_b$means)
  expect_identical(fit_par_a$sds,   fit_par_b$sds)

  # (2) The parallel scatter equals the serial scatter up to summation order.
  expect_equal(fit_par_a$means, fit_serial$means, tolerance = 1e-6)
  expect_equal(fit_par_a$sds,   fit_serial$sds,   tolerance = 1e-4)
  expect_identical(fit_par_a$N, fit_serial$N)
})
