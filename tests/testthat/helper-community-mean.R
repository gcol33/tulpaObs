# Shared community-mean recovery assertion (gcol33/tulpaObs#155).
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
# site. Each family's simulator supplies the realized mean as `beta_real`.
expect_community_mean <- function(fit, real, tol) {
  est <- unname(fit$means[seq_along(real)])
  stopifnot(length(tol) == length(real))
  for (k in seq_along(real)) {
    testthat::expect_lt(abs(est[k] - real[k]), tol[k])
  }
  invisible(est)
}
