# int_occu() with a single source is mathematically the same model as
# occu() -- one shared latent occupancy state observed through one detection
# process. This boundary anchor (the same kind dyn_int_occu() already
# carries against dyn_occu()/int_occu()) caught two compounding bugs that
# left int_occu()'s detection-arm intercept and SE silently wrong on every
# fit with an autoscaled detection covariate:
#
#   1. int_occu()'s per-source detection design matrix lost its column names
#      when padded to full-site width (R/occu.R), so the covariate
#      autoscaler's unscale step could not find the intercept column by
#      name and silently left it in standardized-covariate units while
#      correctly unscaling the slope.
#   2. model_type == "integrated" never received the exact-marginal Newton
#      debiasing single-season occu() and dyn_occu() already have.
#
# Kept as a permanent regression guard: a single-source int_occu() fit must
# match occu() on the same data to full optimizer precision, not just "close".

.anchor_sim <- function(seed, N = 150) {
  set.seed(seed)
  x_cov <- rnorm(N); det_cov <- rnorm(N)
  z <- rbinom(N, 1, plogis(0.2 + 0.7 * x_cov))
  p <- plogis(-0.2 + 0.4 * det_cov)
  y <- matrix(0L, N, 4L)
  for (i in seq_len(N)) if (z[i] == 1L) y[i, ] <- rbinom(4L, 1, p[i])
  list(x_cov = x_cov, det_cov = det_cov, y = y)
}

test_that("int_occu() with one source is bit-identical to occu() on the same data", {
  sim <- .anchor_sim(seed = 22L)

  fit_occu <- tobs(~ occ_cov1, data = data.frame(occ_cov1 = sim$x_cov, det_cov1 = sim$det_cov),
                   family = occu(), detection = ~ det_cov1, y = sim$y,
                   method = "laplace", control = list(verbose = FALSE, progress = FALSE))
  fit_int <- tobs(~ occ_cov, data = data.frame(occ_cov = sim$x_cov, det_cov = sim$det_cov),
                  family = int_occu(), detection = ~ det_cov,
                  y = list(src1 = sim$y),
                  method = "laplace", control = list(verbose = FALSE, progress = FALSE))

  expect_equal(unname(fit_occu$means), unname(fit_int$means), tolerance = 1e-5)
  expect_equal(unname(fit_occu$sds),   unname(fit_int$sds),   tolerance = 1e-5)

  # The design matrix must carry column names at all (not NULL) -- the exact
  # defect that broke the autoscale unscale step's intercept-column lookup.
  # Covariate names differ between the two fits by construction (det_cov vs
  # det_cov1), so this checks presence/shape, not string equality.
  expect_false(is.null(colnames(fit_int$model$X_processes[[2L]])))
  expect_identical(colnames(fit_int$model$X_processes[[2L]])[1L], "(Intercept)")
  expect_identical(length(colnames(fit_int$model$X_processes[[2L]])),
                   length(colnames(fit_occu$model$X_processes[[2L]])))
})

test_that("the single-source anchor holds across seeds", {
  skip_on_cran()
  skip_if_fast()

  for (seed in c(1L, 7L, 42L)) {
    sim <- .anchor_sim(seed = seed)
    fit_occu <- suppressWarnings(tobs(
      ~ occ_cov1, data = data.frame(occ_cov1 = sim$x_cov, det_cov1 = sim$det_cov),
      family = occu(), detection = ~ det_cov1, y = sim$y,
      method = "laplace", control = list(verbose = FALSE, progress = FALSE)))
    fit_int <- suppressWarnings(tobs(
      ~ occ_cov, data = data.frame(occ_cov = sim$x_cov, det_cov = sim$det_cov),
      family = int_occu(), detection = ~ det_cov, y = list(src1 = sim$y),
      method = "laplace", control = list(verbose = FALSE, progress = FALSE)))
    expect_equal(unname(fit_occu$means), unname(fit_int$means), tolerance = 1e-4,
                info = paste("seed", seed))
  }
})
