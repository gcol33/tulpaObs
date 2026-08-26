# Shared community-mean recovery assertion.
#
# A community simulator draws its per-species coefficients around a POPULATION
# mean, so the mean a seed actually realizes sits `beta_sd / sqrt(S)` from that
# constant -- 0.10 at S = 16. Scoring the fit against the constant therefore
# spends most of the budget on draw noise before it can measure the estimator:
# on the ms_count factor fixture, seed 4 sits 0.142 from the nominal intercept
# and 0.002 from the realized one, so `tolerance = 0.2` left ~0.06 for estimator
# error. Score against `colMeans(bs)` -- the seed's own realized mean -- and the
# budget becomes a pure estimator budget.
#
# Compare ABSOLUTELY. testthat's numeric `tolerance` is relative only while the
# target exceeds the tolerance and absolute below it (all.equal.numeric), so one
# `tolerance = 0.35` silently meant +-0.28 on a 0.8 slope and +-0.35 on a 0.2
# intercept.
#
# `tol` is per coefficient and is NOT transferable between families: it has to
# come from that family's own multi-seed deviation spread, recorded at the call
# site.
#
# Where the realized mean comes from: a simulator that reports it names it after
# the constant it belongs to, so `mu_lambda` is paired with `mu_lambda_real`
# (`simulate_ms_abun()` reports `mu_lambda_real`, `mu_p_real`, `mu_log_r_real`,
# `mu_omega_real`). The single-arm fixtures built inside test files keep the
# older `beta_real` spelling for the one mean they carry.
#
# The same split governs INTERVALS, in the opposite direction. An interval for a
# community mean targets the population constant and its width already carries
# the `sd^2 / n_species` draw term, so coverage is scored against the constant.
# Interval SCALE is not: `sd(est - constant) / mean(se)` charges the estimator
# for a draw measured on however many seeds the sweep ran, which at ~20 seeds
# has a 95% range of about +-22% on its own. #280 and #285 read an 18-species
# block that drew 18-21% wide as a 1.28x-narrow SE; it calibrates against
# `mu_log_r_real` (NOTES_measurements.md).
expect_community_mean <- function(fit, real, tol) {
  est <- unname(fit$means[seq_along(real)])
  stopifnot(length(tol) == length(real))
  for (k in seq_along(real)) {
    testthat::expect_lt(abs(est[k] - real[k]), tol[k])
  }
  invisible(est)
}
