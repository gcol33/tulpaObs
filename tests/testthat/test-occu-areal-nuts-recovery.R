# =============================================================================
# Areal-field exposure + recovery for single-season occupancy NUTS
# (gcol33/tulpaObs#142). The occu() NUTS path samples an icar/bym2 areal field
# but used to leave the field params unnamed ("param[k]") and fit$spatial_field
# NULL, so the field could only be smoke-tested. The engine already exports the
# block offsets on ParamLayout; occu_fit.cpp now emits them as `spatial_layout`
# and names the columns, and `.tobs_areal_field` slices `fit$spatial_field`.
#
# For icar/car_proper the field node enters the logit predictor directly, so the
# centred posterior-mean node is the per-cell surface; occupancy fields are
# weakly identified (one binary obs/node), so the recovery bar is a modest field
# correlation, not a tight point match. All heavy -> skip_if_fast()-gated.
# =============================================================================

.chain_adj <- function(N) {
  a <- matrix(0L, N, N)
  for (s in 1:N) { if (s > 1) a[s, s - 1] <- 1L; if (s < N) a[s, s + 1] <- 1L }
  a
}
.smooth_field <- function(N, sd, ph) {
  f <- sin(2 * pi * seq_len(N) / N + ph); f <- f - mean(f); f / sd(f) * sd
}

test_that("occu() + icar NUTS exposes and recovers the areal field", {
  skip_on_cran(); skip_if_fast()
  cors <- numeric(0)
  for (s in seq_len(6L)) {
    set.seed(100 + s)
    n <- 50L; adj <- .chain_adj(n); f <- .smooth_field(n, 1.0, 0.7)
    J <- 6L; psi <- plogis(0.2 + f); z <- rbinom(n, 1, psi)
    y <- matrix(0L, n, J)
    for (i in 1:n) y[i, ] <- if (z[i]) rbinom(J, 1, 0.55) else 0L
    fit <- tobs(~ icar(graph = adj), data = data.frame(idx = 1:n), family = occu(),
                y = y, detection = ~ 1, method = "nuts",
                control = list(n.iter = 600L, n.warmup = 300L, n.chains = 1L,
                               seed = 1L, verbose = FALSE))
    if (s == 1L) {
      # Structural: the field is exposed, named, and the right length.
      expect_false(is.null(fit$spatial_layout))
      expect_identical(fit$spatial_layout$type, "icar")
      expect_length(fit$spatial_field, n)
      expect_true(any(grepl("spatial_field\\[", names(fit$means))))
      expect_false(any(grepl("^param\\[", names(fit$means))))  # nothing unnamed
    }
    cors <- c(cors, suppressWarnings(cor(fit$spatial_field, f)))
  }
  # Measured mean ~0.81 (seeds vary 0.73-0.91); assert a conservative floor.
  expect_gt(mean(cors), 0.6)
})

test_that("occu() + bym2 NUTS names every field and hyperparameter column", {
  skip_on_cran(); skip_if_fast()
  set.seed(101)
  n <- 40L; adj <- .chain_adj(n); f <- .smooth_field(n, 1.0, 0.7)
  J <- 6L; psi <- plogis(0.2 + f); z <- rbinom(n, 1, psi)
  y <- matrix(0L, n, J)
  for (i in 1:n) y[i, ] <- if (z[i]) rbinom(J, 1, 0.55) else 0L
  fit <- tobs(~ bym2(graph = adj), data = data.frame(idx = 1:n), family = occu(),
              y = y, detection = ~ 1, method = "nuts",
              control = list(n.iter = 400L, n.warmup = 200L, n.chains = 1L,
                             seed = 1L, verbose = FALSE))
  nm <- names(fit$means)
  expect_identical(fit$spatial_layout$type, "bym2")
  expect_true(any(grepl("spatial_field\\[", nm)))   # structured (phi) block
  expect_true(any(grepl("spatial_theta\\[", nm)))   # unstructured block
  expect_true("log_sigma_spatial" %in% nm)
  expect_true("logit_rho_spatial" %in% nm)
  expect_false(any(grepl("^param\\[", nm)))          # nothing falls through
})
