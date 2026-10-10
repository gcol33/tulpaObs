# tulpa's non-copy `iid` prior block is non-centred: a random intercept enters
# the predictor as `sigma_k * u` with `u ~ N(0, 1)` and `sigma_k` the cell's
# `b<k>.sigma` node. The BLUPs, their SDs and the predict / WAIC draws must read
# the effect `sigma_k * u`, with the within-cell posterior variance included.
# Read as `u`, a detection random intercept on habitat inflated p_waic to
# 2-4e5 on the full MOT data and returned ranef() SDs of ~1e-5.

.re_scale_data <- function(seed = 415L, n_cells = 30L, J = 30L) {
  set.seed(seed)
  hab_eff <- c(-2.5, -1.2, -0.4, 0, 0.6, 1.1, 1.6, 2.1)
  lev <- paste0("H", seq_along(hab_eff))
  rows <- lapply(seq_len(n_cells), function(cell) {
    z <- rbinom(1, 1, plogis(0.6 + 0.1 * (cell - n_cells / 2)))
    x1 <- rnorm(J); h <- sample.int(length(lev), J, replace = TRUE)
    occur <- if (z) rbinom(J, 1, plogis(0.2 + 0.4 * x1 + hab_eff[h])) else rep(0L, J)
    data.frame(cell_idx = cell, site_key = as.character(cell), visit = seq_len(J),
               x1 = x1, habitat = factor(lev[h], levels = lev),
               occur = as.integer(occur),
               cover.flat = ifelse(occur == 1,
                                   pmin(pmax(rbeta(J, 2, 4), 1e-4), 1 - 1e-4), 0))
  })
  dd <- do.call(rbind, rows)
  adj <- matrix(0L, n_cells, n_cells)
  for (i in seq_len(n_cells - 1L)) adj[i, i + 1L] <- adj[i + 1L, i] <- 1L
  list(data = dd, adj = adj, hab_eff = hab_eff)
}

.re_scale_fit <- function(d, detection) {
  adj <- d$adj
  tobs(occurrence = ~ 1 + spatial(~ 1 || cell_idx, graph = adj), data = d$data,
       family = occu_cover(response = "beta", cover.aggregate = "none"),
       detection = detection, positive = ~ 1, response = "occur",
       y.pos = "cover.flat", site = "site_key", visit = "visit",
       det.covs = c("x1", "habitat"), occ.covs = "cell_idx",
       method = "nested_laplace",
       control = list(engine = "joint", max.iter = 200L, verbose = FALSE,
                      progress = FALSE))
}

test_that("a detection random intercept reads the scaled effect", {
  skip_on_cran()
  d <- .re_scale_data()
  fx <- suppressWarnings(.re_scale_fit(d, ~ x1 + habitat))
  re <- suppressWarnings(.re_scale_fit(d, ~ x1 + (1 | habitat)))

  t <- re$re[["p"]]
  expect_false(is.null(t))
  expect_identical(t$scale_axes, grep("^b[0-9]+[.]sigma$", t$scale_axes,
                                      value = TRUE))
  # Level means: the fixed fit's intercept + contrast against the RE fit's
  # intercept + centred BLUP. Unscaled, the BLUPs were the unit-scale latent.
  b_fx <- fx$means[grepl("^p_", names(fx$means))]
  lev_fx <- b_fx[["p_(Intercept)"]] +
    c(0, b_fx[paste0("p_habitat", levels(d$data$habitat)[-1])])
  lev_re <- re$means[["p_(Intercept)"]] + t$blup
  expect_gt(stats::cor(lev_fx, lev_re), 0.98)
  expect_lt(max(abs((lev_fx - mean(lev_fx)) - (lev_re - mean(lev_re)))), 0.5)

  # Posterior SDs of the BLUPs: the scale of a habitat contrast's SE, not 1e-5.
  expect_true(all(t$blup_sd > 0.05 & t$blup_sd < 1.5))

  # WAIC: shrinking eight habitat levels costs at most a few effective
  # parameters against estimating them as fixed effects.
  w_fx <- suppressWarnings(waic(fx)); w_re <- suppressWarnings(waic(re))
  p_fx <- w_fx$estimates["p_waic", "Estimate"]
  p_re <- w_re$estimates["p_waic", "Estimate"]
  expect_lt(abs(p_re - p_fx), 5)
  expect_lt(abs(w_re$estimates["elpd_waic", "Estimate"] -
                w_fx$estimates["elpd_waic", "Estimate"]), 10)
})
