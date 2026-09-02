# 95% CI coverage for mu_log_r under ms_abun(negbin) with a per-species
# dispersion RE (log_r_s ~ N(mu_log_r, sigma_log_r)). Registry/gate tests for
# this fixture are in test-ms-abun-nb-rs.R; point recovery is in
# test-ms-abun-nb-rs-recovery.R.

test_that("ms_abun(negbin) mu_log_r 95% CI covers at the nominal rate", {
  skip_on_cran()
  skip_if_fast()
  S <- 18L
  SIGMA_LOGR <- 0.5
  n_seed <- 20L
  covered <- logical(0)
  est <- numeric(0); se <- numeric(0); truth <- numeric(0); real <- numeric(0)
  kept <- integer(0)
  refused <- integer(0)
  unconverged <- integer(0)
  for (s in seq_len(n_seed)) {
    seed <- 500L + s
    sim <- simulate_ms_abun(n_species = S, N = 100, J = 5,
                            n_abund_covs = 1, n_det_covs = 1,
                            mu_lambda = c(log(5), 0.4), mu_p = c(0.3, -0.3),
                            sd_lambda = 0.4, sd_p = 0.35,
                            mixture = "negbin", size = 5,
                            sigma_logr = SIGMA_LOGR, seed = seed)
    fit <- tryCatch(
      suppressWarnings(
        tobs(~ abund_cov1, data = sim$data, y = sim$y,
             family = ms_abun(mixture = "negbin"),
             detection = ~ det_cov1,
             species = sim$species, method = "laplace",
             control = list(verbose = FALSE, n.quad = 3L))),
      error = function(e) NULL)
    # A seed whose per-species posterior solve fails delivers no interval at
    # all: the engine refuses an optimum carrying its failure sentinel and this
    # fitter errors, rather than returning the community dispersion block at
    # the values it STARTED from -- `mu_log_r` a hair below `log(r_init)`,
    # `sigma_log_r` exactly the initial 0.5, one `NA` in `r_s`,
    # `converged = TRUE` and the tightest `log_r` SE of the twenty -- which is
    # what let one such seed carry 74% of `sum(z^2)` into a coverage loop that
    # scored it as data.
    #
    # Neither of the two obvious gates catches such a fit: `converged` was
    # `TRUE`, and a `sigma_log_r`-collapse detector keyed on sigma being small
    # misses it because its `sigma_log_r` is exactly the initial value (#250
    # item 3). What catches it is the engine declining, and the per-species
    # solve status the fit now carries -- so both are read here, and both
    # counted, so the denominator is visible rather than assumed.
    if (is.null(fit)) { refused <- c(refused, seed); next }
    if (!converged(fit)) { unconverged <- c(unconverged, seed); next }
    lr   <- fit$means[["log_r"]]
    selr <- fit$sds[["log_r"]]
    lo <- lr - 1.96 * selr; hi <- lr + 1.96 * selr
    covered <- c(covered, sim$truth$mu_log_r >= lo && sim$truth$mu_log_r <= hi)
    est <- c(est, lr); se <- c(se, selr); kept <- c(kept, seed)
    truth <- c(truth, sim$truth$mu_log_r)
    # The mean of the 18 log-dispersions this seed actually drew. The interval
    # targets the population constant, so `truth` is what coverage is scored
    # against; `real` is what separates the estimator from the draw, below.
    real <- c(real, sim$truth$mu_log_r_real)
  }
  z <- (est - truth) / se

  # A seed can clear BOTH gates above -- the engine does not refuse it and
  # `converged()` is TRUE -- and still report `se` as exactly 0, when the
  # community dispersion component collapses onto its boundary. That makes `z`
  # infinite for that seed, and `shapiro.test()` returns a NaN p-value for a
  # sample containing an infinity rather than erroring, so `NaN <= 0.05` renders
  # a degenerate interval as "the errors are not normal" -- a true failure
  # pointing at the wrong quantity. NA and NaN do not do this: `shapiro.test()`
  # drops them. Only an infinity does, which is what a zero SE produces.
  #
  # Checked before the assertions that consume `z`, and it names the seed, so
  # the finding is the interval rather than the normality.
  bad <- which(!is.finite(z))
  expect_true(length(bad) == 0L,
              info = if (length(bad))
                sprintf("non-finite z from seed(s) %s: se = %s, est - truth = %s",
                        paste(kept[bad], collapse = ", "),
                        paste(format(se[bad]), collapse = ", "),
                        paste(format((est - truth)[bad]), collapse = ", ")))
  # A sweep that drops most of its seeds is not a coverage measurement any
  # more, whatever rate the survivors show.
  expect_gte(length(covered), n_seed - 3L)

  # WHY THIS BLOCK DOES NOT ASSERT NOMINAL COVERAGE, and why the assertion it
  # used to carry was not the tripwire it advertised itself as (#285).
  #
  # mu_log_r is a POPULATION mean. Each seed draws 18 log-dispersions around it
  # and the fit sees only those, so the error splits
  #
  #     est - mu  =  (est - mu_real)  +  (mu_real - mu)
  #
  # with `mu_real` the seed's own realized species mean. The second term is
  # N(0, sigma_logr^2 / S) = N(0, 0.118^2) here and belongs inside the interval
  # -- the SE does carry it, as `se^2 = (sigma_hat^2 + c) / S` with one c across
  # a factor of 4.5 in S. But it is measured over only 19 seeds, and it carries
  # 65% of the spread, so an ordinary fluctuation in the DRAW reads as a
  # miscalibrated SE. On seeds 501-520 that is what happened: they drew their
  # species means 18.5% wider than sigma_logr / sqrt(18) (chi^2_18 = 25.3,
  # p = 0.12), which is most of the k = 1.277 that #280 and #285 report.
  #
  # Measured on these 19 fits (dev_notes/_i285_*.R over the 97-fit sweep behind
  # #285; NOTES_measurements.md):
  #
  #   bias  mean(est) - mu                +0.017  (0.42 SE of the mean)
  #   sd(est)                              0.1733
  #   mean reported se                     0.1356
  #   k_raw = sd(est) / mean(se)           1.277
  #   coverage                             15/19 = 0.789
  #   realized draw sd(mu_real - mu)       0.1397  vs nominal 0.1179
  #   k with the draw at its expectation   1.188   (null band [0.81, 1.19])
  #   the same, SE rebuilt at sigma = 0.5  1.130
  #
  # Pooled over 97 fits at 8 / 18 / 36 species, putting the draw at its
  # expectation AND rebuilding the SE at the simulated sigma gives 0.990 /
  # 1.035 / 0.963 -- calibrated at every group count, and the S = 18 spike
  # #285 filed is the seed block's draw, not the estimator.
  #
  # What IS left is that `sigma_hat` is attenuated (0.458 here against a
  # simulated 0.5) and the SE is built on it. That is monotone in S
  # (0.418 / 0.448 / 0.487 at 8 / 18 / 36) and is documented on `?ms_abun`.
  #
  # So a coverage floor here cannot be a defect tripwire: rebuilding the SE at
  # the true sigma -- a strictly better interval than any estimator can deliver
  # -- still covers only 16/19 = 0.842 on these seeds, because their draw
  # excess is untouched by it. The old `expect_lt(mean(covered), 0.90)` was
  # advertised as going green once the SE was fixed and would not have.
  #
  # The four assertions below are therefore: unbiased, normal errors, a gross
  # floor on the raw rate, and -- the estimator property -- the interval scale
  # once the draw is put at its expectation.
  expect_lt(abs(mean(est) - mean(truth)) / (sd(est) / sqrt(length(est))), 2.5)
  expect_gt(stats::shapiro.test(z)$p.value, 0.05)
  expect_gt(mean(covered), 0.70)

  # Interval scale with the realized draw replaced by its expectation. This is
  # the one statistic here that is a property of the estimator rather than of
  # the seed block, so it is the one worth a band. Wide on purpose: at 19 seeds
  # its own null is [0.81, 1.19] at +-2.5 SE, so this is a gross-regression
  # guard (observed 1.19), not a calibration gate. Tightening it needs more
  # seeds, not a smaller number -- #285 measured that separating 1.28 from
  # noise here takes n >= 40.
  k_corr <- sqrt(stats::var(est - real) + SIGMA_LOGR^2 / S) / mean(se)
  expect_gt(k_corr, 0.70)
  expect_lt(k_corr, 1.45)
})
