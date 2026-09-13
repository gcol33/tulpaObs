# =============================================================================
# test-occu-field-offset-fitted.R - an areal field on occu() reaches the
# post-fit readers (#318).
#
#   - fitted()$psi is plogis(intercept + field), exactly, and tracks the truth
#   - the pointwise log-likelihood behind waic() scores the field: the field fit
#     beats an intercept-only fit on elpd_waic
# =============================================================================

.ofo_data <- function() {
  g <- 6L
  adj <- rook_adj(g)
  set.seed(11)
  xy <- expand.grid(c = seq_len(g), r = seq_len(g))
  f_true <- as.numeric(scale(sin(xy$r / 2) + cos(xy$c / 2))) * 1.2
  psi_true <- plogis(0.2 + f_true)
  n <- g * g
  z <- rbinom(n, 1, psi_true)
  y <- matrix(rbinom(n * 4L, 1, rep(z, 4L) * 0.45), n, 4L)
  list(adj = adj, psi_true = psi_true, y = y,
       d = data.frame(site = seq_len(n), cell = seq_len(n)))
}

test_that("fitted()$psi carries the areal field on an occu() fit", {
  s <- .ofo_data()
  adj <- s$adj
  fit <- suppressMessages(tobs(~ 1 + icar(graph = adj), data = s$d,
    family = occu(), detection = ~ 1, y = s$y, method = "nested_laplace",
    control = list(progress = FALSE, verbose = FALSE)))
  psi <- fitted(fit)$psi
  expect_equal(psi, plogis(as.numeric(fit$means[1]) + fit$spatial_field),
               tolerance = 1e-12)
  expect_gt(stats::sd(psi), 0)
  expect_gt(stats::cor(psi, s$psi_true), 0.6)
})

test_that("waic() scores the areal field on an occu() fit", {
  s <- .ofo_data()
  adj <- s$adj
  ctl <- list(progress = FALSE, verbose = FALSE)
  fit_field <- suppressMessages(tobs(~ 1 + icar(graph = adj), data = s$d,
    family = occu(), detection = ~ 1, y = s$y, method = "nested_laplace",
    control = ctl))
  fit_flat <- suppressMessages(tobs(~ 1, data = s$d, family = occu(),
    detection = ~ 1, y = s$y, method = "laplace", control = ctl))
  expect_gt(suppressWarnings(waic(fit_field))$elpd_waic,
            suppressWarnings(waic(fit_flat))$elpd_waic)
})
