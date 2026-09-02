# =============================================================================
# test-sbc.R
# -- posterior simulation-based calibration on tobs families, through
# the engine's `tulpa::sbc()` front door.
#
# Three tiers, deliberately:
#
#   CONTRACT (always on). The six callbacks plus `group_ids` exist and do what
#   the engine's contract says: the pooled data set carries BOTH data sets, the
#   replicate's cells are FRESH, the fitting call rebuilds with the pooled
#   graph. These cost one small fit or none.
#
#   GENERATOR (skip_on_cran). The replicate generator draws from the same law
#   the likelihood is written in. This is the one failure that would make SBC
#   report a false alarm against the engine, so it is checked against the
#   model's own eta assembly rather than assumed.
#
#   ACCEPTANCE (skip_if_fast + skip_on_cran). The calibration measurement
#   itself: the reported posterior uniform, a deliberately mis-scaled control
#   not. A band nothing can fail is not evidence, so both halves are asserted.
#   ~13 min at n.sim = 100.
# =============================================================================

# The dispersion grid the fixture puts phi_pos on. Without it the joint engine
# holds the cover dispersion at a value it sets per data set, so the replicate
# would be generated at one value and refitted under another.
.SBC_PHI_GRID <- exp(seq(log(0.20), log(0.90), length.out = 17L))

# The field-SD grid, PINNED. A defaulted axis is auto-recentred per fit, so the
# observed fit and the augmented refits would integrate `sigma` over different
# supports and the truth would be scored against a predictive on a support it is
# not from. A user grid is not recentred, which is the same reason tulpa's own
# SBC fixture fixes its hyperparameter support.
#
# BOTH AXES ARE RESOLVED, NOT MERELY PINNED: 21 sigma nodes and 17 phi_pos ones,
# where the fixture used to carry 9 and 7. An outer-grid axis IS the whole
# support of the quantity it carries, so its resolution is the resolution of the
# predictive the truth is ranked against, and a rank read against a continuous
# uniform scores a coarse axis as a departure whatever the fit does: at 9 sigma
# nodes 1000 draws hold 6 distinct values and the rank ECDF is a 6-step function.
# Measured on the SAME generator and the same seed, n.sim = 100, `sigma` p_unif
# 4.9e-04 -> 0.056 -> 0.714 and `disp` 8.7e-09 -> 1.9e-14 -> 0.819 across (9, 7)
# -> (21, 7) -> (21, 17) nodes, with every mean PIT already within noise of 0.5
# at the coarse grids -- a step-function ECDF, not a location shift. That is what
# lets the acceptance set below gate every scored quantity. The cost is a
# 13-minute acceptance run instead of a 3.6-minute one, and the whole point of
# the run is the read.
.SBC_SIGMA_GRID <- exp(seq(log(0.15), log(2.0), length.out = 21L))

# A small coupled occu_cover fixture: shared ICAR field on the occurrence arm,
# copied onto the cover arm, dispersion on the outer grid.
.sbc_fixture <- function(N = 30L, J = 4L, seed = 707L, sigma = 0.8,
                         alpha = 1.0, phi.grid = .SBC_PHI_GRID) {
  adj <- chain_adj(N)
  sim <- simulate_occu_cover(
    N = N, J = J, positive = "lognormal", adj = adj,
    beta_occ = c(0.2, 0.6), beta_p = c(0.4, -0.5),
    beta_pos = c(log(0.25), 0.3),
    sigma = sigma, alpha = alpha, sigma_pos = 0.4, seed = seed)
  long <- data.frame(site_id = rep(seq_len(N), each = J),
                     visit = rep(seq_len(J), times = N),
                     y = as.vector(t(sim$y)),
                     det_cov1 = sim$visit_data$det_cov1,
                     pos_cov1 = sim$visit_data$pos_cov1)
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
  ctl <- list(verbose = FALSE, engine = "joint", progress = FALSE)
  if (!is.null(phi.grid)) {
    ctl$phi.grid.pos <- phi.grid
    ctl$sigma.grid <- .SBC_SIGMA_GRID
  }
  fit <- suppressWarnings(tobs(
    formula = ~ occ_cov1 + icar(graph = adj),
    data = cbind(data.frame(site_id = seq_len(N)), sim$data),
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1 + share(spatial()),
    y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "nested_laplace", control = ctl))
  list(fit = fit, sim = sim, adj = adj, phi.grid = phi.grid)
}

# The refit control: the same pinned hyperparameter support the observed fit
# integrated, so both stages read the same grid.
.sbc_fit_control <- function(fx) {
  if (is.null(fx$phi.grid)) return(list())
  list(phi.grid.pos = fx$phi.grid, sigma.grid = .SBC_SIGMA_GRID)
}


# ---------------------------------------------------------------------------
# Contract
# ---------------------------------------------------------------------------

test_that("sbc() refuses a non-fit and an unregistered family", {
  expect_error(sbc(list()), "No sbc method is registered")

  # double_observer() graduated to registered; a synthetic fake family keeps
  # this test's premise -- a fit whose family is genuinely absent from the
  # registry -- true regardless of what graduates next, matching
  # test-sbc-registry.R's own "roster names what is registered" test.
  fake <- structure(list(model = list()), class = "tobs_fit",
                    tobs_family = list(name = "not_a_family"))
  expect_error(sbc(fake), "not registered for family")
  expect_error(sbc(fake), "occu_cover")
})


test_that("the callback list satisfies the engine's posterior contract", {
  skip_on_cran()
  fx <- .sbc_fixture(N = 20L, J = 3L)
  m <- sbc(fx$fit, model.only = TRUE,
                fit.control = .sbc_fit_control(fx))

  expect_true(all(c("data_obs", "fit", "draw_theta", "simulate", "pool",
                    "arms", "group_ids") %in% names(m)))
  expect_true(all(vapply(m[c("fit", "draw_theta", "simulate", "pool", "arms",
                             "group_ids")], is.function, logical(1))))
  expect_false(is.function(m$data_obs))

  # Both field scales and the copy are scored, and nothing is silently fixed
  # once the dispersion is on the grid.
  q <- attr(m, "quantities")
  expect_true(all(c("sigma", "sigma_pos_field", "alpha", "disp") %in% q))
  expect_true(all(c("psi_(Intercept)", "psi_occ_cov1", "p_(Intercept)",
                    "pos_(Intercept)") %in% q))
  expect_length(attr(m, "fixed"), 0L)

  # `draw_theta` is reproducible from its seed alone and moves with it, which
  # is what the driver's RNG split relies on.
  t1 <- m$draw_theta(fx$fit, 11L)
  expect_identical(t1, m$draw_theta(fx$fit, 11L))
  expect_false(isTRUE(all.equal(unname(t1), unname(m$draw_theta(fx$fit, 12L)))))
  expect_true(all(q %in% names(t1)))
})


test_that("pooling keeps both data sets and the replicate's cells are fresh", {
  skip_on_cran()
  fx <- .sbc_fixture(N = 20L, J = 3L)
  m <- sbc(fx$fit, model.only = TRUE,
                fit.control = .sbc_fit_control(fx))

  obs <- m$data_obs
  rep <- m$simulate(m$draw_theta(fx$fit, 5L), 660005L)
  pl  <- m$pool(obs, rep)

  n_o <- length(m$group_ids(obs)); n_r <- length(m$group_ids(rep))
  expect_identical(length(unique(m$group_ids(pl))),
                   length(unique(m$group_ids(obs))) +
                     length(unique(m$group_ids(rep))))
  expect_identical(nrow(pl$y), n_o + n_r)
  expect_identical(nrow(pl$cells), n_o + n_r)

  # The observed rows survive the pooling: fitting the replicate alone would be
  # ordinary SBC under a hand-made prior, not the posterior experiment.
  expect_equal(pl$y[seq_len(n_o), ], obs$y)
  expect_equal(pl$y[n_o + seq_len(n_r), ], rep$y)

  # Fresh cells means a BLOCK-DIAGONAL graph: the two field blocks are a-priori
  # independent, which is what makes the replicate conditionally independent of
  # the observed data given theta.
  expect_identical(dim(pl$graph), c(n_o + n_r, n_o + n_r))
  expect_true(all(pl$graph[seq_len(n_o), n_o + seq_len(n_r)] == 0L))
  expect_true(all(pl$graph[n_o + seq_len(n_r), seq_len(n_o)] == 0L))
  expect_equal(pl$graph[seq_len(n_o), seq_len(n_o)], unname(as.matrix(obs$graph)))
})


test_that("the refit rebuilds the same call on the pooled graph", {
  skip_on_cran()
  fx <- .sbc_fixture(N = 20L, J = 3L)
  m <- sbc(fx$fit, model.only = TRUE,
                fit.control = .sbc_fit_control(fx))

  # Refitting the observed data alone must reproduce the observed fit: the
  # rebuilt formulas, family, method and control are the ones it came from.
  again <- m$fit(m$data_obs)
  expect_equal(again$means, fx$fit$means, tolerance = 1e-8)

  rep <- m$simulate(m$draw_theta(fx$fit, 3L), 660003L)
  pooled <- m$fit(m$pool(m$data_obs, rep))
  expect_s3_class(pooled, "tobs_fit")
  expect_identical(pooled$model$n_sites, 40L)
  expect_true(all(is.finite(pooled$means)))
  # The cover arm still sees the field: without the copy `alpha` is pinned and
  # the coupled model has quietly become an uncoupled one.
  expect_true("alpha" %in% names(pooled$means))
  expect_gt(pooled$means[["alpha"]], 0)
})


test_that("a dispersion the engine holds fixed is reported, not scored", {
  skip_on_cran()
  fx <- .sbc_fixture(N = 20L, J = 3L, phi.grid = NULL)
  expect_warning(m <- sbc(fx$fit, model.only = TRUE),
                 "holds the cover dispersion fixed")
  expect_true("disp" %in% attr(m, "fixed"))
  expect_false("disp" %in% attr(m, "quantities"))
})


test_that("a visit design that cannot be inverted is refused, not guessed", {
  skip_on_cran()
  fx <- .sbc_fixture(N = 20L, J = 3L)
  bad <- fx$fit
  colnames(bad$model$X_det_visit) <- "poly(det_cov1, 2)1"
  expect_error(tulpaObs:::.tobs_sbc_visit_matrices(bad$model),
               "plain numeric")
})


# ---------------------------------------------------------------------------
# Generator
# ---------------------------------------------------------------------------

test_that("the replicate generator draws from the law the likelihood scores", {
  skip_on_cran()
  fx <- .sbc_fixture(N = 30L, J = 4L)
  m <- sbc(fx$fit, model.only = TRUE,
                fit.control = .sbc_fit_control(fx))
  spec <- environment(m$simulate)$spec

  # Field off, so every moment below has a closed form in theta alone.
  th <- m$draw_theta(fx$fit, 1L)
  th[["sigma"]] <- 0; th[["alpha"]] <- 0; th[["sigma_pos_field"]] <- 0
  th[["disp"]] <- 0.4
  bl <- attr(th, "blocks")
  e <- tulpaObs:::.occu_cover_eta_from_par(spec$model, th[bl$occ], th[bl$det],
                                           th[bl$pos])

  # EVERY assertion below is ABSOLUTE and sized against the estimator's own
  # Monte-Carlo error, not against a relative tolerance. Measured over 4000
  # replicates the detection rate converges on its target as 1/sqrt(R) (0.4901
  # at 120, 0.4776 at 1000, 0.4748 at 4000, target 0.4733), and the per-
  # replicate SD puts the standard error at 600 replicates near 0.0022 -- so a
  # 0.01 band is about four standard errors, while the relative 3% this used to
  # assert was under three and failed on a high draw.
  reps <- lapply(seq_len(600L), function(i)
    tulpaObs:::.tobs_sbc_sim_occu_cover(spec, th, 90000L + i))
  det_rate <- mean(vapply(reps, function(r) mean(r$y), numeric(1)))
  # E[y] = psi * p, averaged over the design.
  expect_lt(abs(det_rate - mean(e$psi * e$p_mat)), 0.01)

  # Cover at a detected visit is lognormal on the cover predictor, so the mean
  # of log-cover is that predictor and its SD is the dispersion.
  lc <- unlist(lapply(reps, function(r) log(r$y_pos[r$y == 1L])))
  ep <- unlist(lapply(reps, function(r) e$ep_mat[r$y == 1L]))
  expect_lt(abs(mean(lc - ep)), 0.02)
  expect_lt(abs(stats::sd(lc - ep) - 0.4), 0.02)
})


test_that("the field enters the replicate at the scale the engine's block carries", {
  skip_on_cran()
  fx <- .sbc_fixture(N = 40L, J = 4L)
  m <- sbc(fx$fit, model.only = TRUE,
                fit.control = .sbc_fit_control(fx))
  spec <- environment(m$simulate)$spec

  # The joint engine's ICAR block is the RAW graph precision Q = D - W at tau
  # = 1 with the amplitude in the arm scale, so the field `sigma` multiplies
  # is x ~ N(0, Q^+): geo-mean marginal variance scale_q, NOT 1.
  # `.occu_cover_draw_icar_field()` returns the Sorbye-Rue NORMALISED draw --
  # the convention `simulate_occu_cover()` states its `sigma` in -- so the
  # generator carries the constant that maps one to the other. Generating at
  # the normalised scale instead refits every replicate under a field
  # sqrt(scale_q) wider and piles both arm SDs at the top of their rank
  # support.
  #
  # scale_q is recomputed here from the graph, independently of the package
  # helper the generator uses, so the two cannot agree by sharing a bug.
  Q <- diag(rowSums(fx$adj)) - fx$adj
  eig <- eigen(Q, symmetric = TRUE)
  pos <- eig$values > 1e-10
  Qplus_diag <- rowSums(sweep(eig$vectors[, pos, drop = FALSE]^2, 2L,
                              eig$values[pos], "/"))
  scale_q <- exp(mean(log(Qplus_diag)))
  expect_equal(spec$field_scale, sqrt(scale_q), tolerance = 1e-8)

  # On this graph the two conventions are far apart, so a generator that lost
  # the constant could not pass by rounding.
  expect_gt(spec$field_scale, 2)

  # The field as GENERATED carries the engine's marginal width.
  set.seed(9L)
  F <- spec$field_scale * tulpaObs:::.occu_cover_draw_icar_field(fx$adj, 400L)
  expect_equal(exp(mean(log(apply(F, 1L, stats::var)))), scale_q,
               tolerance = 0.2)
  expect_equal(mean(colMeans(F)), 0, tolerance = 1e-8)
})


# ---------------------------------------------------------------------------
# Acceptance
# ---------------------------------------------------------------------------

test_that("occu_cover posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  fx <- .sbc_fixture(N = 50L, J = 6L, seed = 707L)
  res <- sbc(fx$fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = c("wide", "narrow"),
                  fit.control = .sbc_fit_control(fx), seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  pu <- function(arm, qs) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  # WHAT IS GATED, and why it is not the in-band indicator. Each band holds at
  # 0.95 SIMULTANEOUSLY within one quantity, so requiring six of them to hold
  # at once fails about a quarter of the time on a perfectly calibrated
  # algorithm. The gate is the smallest uniformity p-value over the set, at a
  # threshold well below any per-quantity level, which is a multiplicity-aware
  # regression guard rather than a coin flip.
  #
  # EVERY scored quantity is in the set: the six arm coefficients, the DERIVED
  # copy scale `alpha`, both field-scale reads (the occurrence-arm SD the outer
  # grid integrates and the cover-arm SD derived from it), and the cover
  # dispersion. The two field SDs used to pile at the top of their support
  # (81/100 and 87/100 in the top decile) because the replicate's field was
  # drawn one Sorbye-Rue constant narrower than the law the refit inverts -- the
  # engine's ICAR block is the raw Q = D - W, so its `sigma` multiplies a field
  # of geo-mean marginal SD sqrt(scale_q), not 1. Held at THIS grid with the
  # generator's constant put back to 1, they return to 42/100 and 52/100 in the
  # top decile, so what moved them is the generator and not the resolution.
  qs <- c("psi_(Intercept)", "psi_occ_cov1", "p_(Intercept)", "p_det_cov1",
          "pos_(Intercept)", "pos_pos_cov1", "alpha", "sigma",
          "sigma_pos_field", "disp")
  ok <- pu("posterior", qs)
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-3)

  # A deliberately mis-scaled posterior has to fail the same read, or the band
  # is not measuring anything. Same quantities, same simulations, same fits --
  # only the reported width differs.
  bad <- pu("narrow", qs)
  expect_lt(min(bad), 1e-3)
  # The copy scale carries it on its own: calibrated as reported, far outside
  # once its width is distorted.
  expect_gt(ok[["alpha"]], 0.05)
  expect_lt(bad[["alpha"]], 1e-6)
  expect_false(rp$inside[rp$arm == "narrow" & rp$quantity == "alpha"])
})
