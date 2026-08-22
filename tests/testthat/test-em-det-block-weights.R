# The M-step detection block of the EM-Laplace path.
#
# The E-step weight w_i = P(z_i = 1 | y_i, theta) is continuous in (0, 1) and
# belongs on the block as a per-observation `weights` vector: that is what makes
# the M-step maximise the expected complete-data log-likelihood. The single
# alternative -- folding w_i into the response by count-scaling -- rounds a
# continuous weight onto integer counts, which erases a low-weight row's TRIALS
# while keeping its (zero) successes and so biases p upward.
#
# The row set is the part a mesh field constrains: the block's projection A
# carries one row per input row, so a field pins the row set and a near-empty
# row must come through at weight ~0 rather than being dropped. Without a field
# the row is dropped, which is the same statement.
#
# The whole encoding rests on the engine consuming `weights` on the spatial
# Laplace route (tulpa >= 0.0.184); the last block asserts that contract
# directly so an engine regression fails here rather than silently reverting the
# M-step to an unweighted fit.

.det_block_data <- function(n = 80, J = 3, seed = 4) {
  set.seed(seed)
  coords <- cbind(runif(n), runif(n))
  x <- rnorm(n); dcov <- rnorm(n)
  z <- rbinom(n, 1, plogis(0.1 + 0.7 * x))
  p <- plogis(-0.3 + 0.4 * dcov)
  y <- matrix(0L, n, J)
  for (i in seq_len(n)) if (z[i] == 1L) y[i, ] <- rbinom(J, 1, p[i])
  y[1, ] <- NA_integer_                       # a site with no valid visit
  list(data = data.frame(occ_cov = x, det_cov = dcov,
                         lon = coords[, 1], lat = coords[, 2]),
       y = y, coords = coords, n = n)
}

.det_block_spde <- function(coords, arm = c(FALSE, TRUE)) {
  sp <- tulpaObs:::.tobs_term_spde(coords = coords, max_edge = c(0.4, 0.8),
                                   nu = 1, prior_range = c(0.3, 0.5),
                                   prior_sigma = c(0.8, 0.5))
  sp$shared <- arm
  sp
}

# Weights that put most rows at a genuinely fractional value, which is where a
# rounded encoding and a weighted one part company. A detection proves
# occupancy, so the M-step pins those rows at 1 whatever the E-step returned;
# `eff` is the weight the block actually sees.
.det_block_weights <- function(y, seed = 11) {
  set.seed(seed)
  w <- rbeta(nrow(y), 0.7, 1.6)
  w[seq_len(5)] <- c(1, 1, 1e-9, 0.02, 0.5)
  any_det <- rowSums(y == 1L, na.rm = TRUE) > 0
  eff <- w; eff[any_det] <- 1
  list(w = w, eff = eff, n_valid = rowSums(!is.na(y)))
}

test_that("the single-season detection block carries the E-step weight, field or not", {
  skip_on_cran()
  fx <- .det_block_data()
  model <- tulpaObs:::.tobs_build_single(~ occ_cov, ~ det_cov, fx$data, fx$y)
  wt <- .det_block_weights(fx$y)

  plain <- tulpaObs:::build_single_callbacks(model, NULL)$m_step_encode(wt$w)$det
  expect_false(is.null(plain$weights))
  expect_equal(length(plain$weights), length(plain$y))
  # y / n_trials stay the raw detection counts; the weight is not folded in.
  expect_true(all(plain$y <= plain$n_trials))
  expect_true(all(plain$n_trials > 0L))
  expect_equal(plain$weights, wt$eff[wt$eff > 1e-6 & wt$n_valid > 0L])

  skip_if_no_tulpamesh()
  sp <- .det_block_spde(fx$coords)
  fielded <- tulpaObs:::build_single_callbacks(model, sp)$m_step_encode(wt$w)$det
  expect_false(is.null(fielded$weights))
  expect_false(is.null(fielded$spatial))
  expect_equal(fielded$weights, wt$eff)
})

test_that("a mesh field pins the detection block to the projection's row set", {
  skip_on_cran()
  skip_if_no_tulpamesh()
  fx <- .det_block_data()
  model <- tulpaObs:::.tobs_build_single(~ occ_cov, ~ det_cov, fx$data, fx$y)
  sp <- .det_block_spde(fx$coords)
  wt <- .det_block_weights(fx$y)

  blk <- tulpaObs:::build_single_callbacks(model, sp)$m_step_encode(wt$w)$det
  n_proj <- nrow(sp$tulpa_spec$A)
  expect_equal(length(blk$y), n_proj)
  expect_equal(nrow(blk$X), n_proj)
  expect_equal(length(blk$weights), n_proj)

  # Rows a count-scaled encoding would erase (round(w * n_valid) == 0) keep
  # their real visit counts and enter at their real weight instead.
  erased <- which(round(wt$eff * wt$n_valid) == 0L)
  expect_gt(length(erased), 0L)
  expect_equal(blk$n_trials[erased], as.integer(wt$n_valid[erased]))
  expect_equal(blk$weights[erased], wt$eff[erased])

  # The site with no valid visit is the (0, 0) row the field forces us to keep.
  expect_equal(blk$n_trials[1L], 0L)
  expect_equal(blk$y[1L], 0L)

  # Without the field the same rows drop out, which is what weight ~0 expresses.
  plain <- tulpaObs:::build_single_callbacks(model, NULL)$m_step_encode(wt$w)$det
  expect_lt(length(plain$y), n_proj)
  expect_equal(sum(wt$eff > 1e-6 & wt$n_valid > 0L), length(plain$y))
})

test_that("the integrated per-source detection blocks encode the same way", {
  skip_on_cran()
  fx <- .det_block_data(n = 80, J = 2)
  y_list <- list(src1 = fx$y, src2 = fx$y[, 1, drop = FALSE])
  model <- tulpaObs:::.tobs_build_integrated(~ occ_cov, ~ det_cov, fx$data, y_list)
  wt <- .det_block_weights(fx$y)
  # A source's own detections pin its rows at 1; pool across sources the way
  # the E-step does.
  det_any <- lapply(y_list, function(m) rowSums(m == 1L, na.rm = TRUE) > 0)
  eff <- lapply(det_any, function(d) { e <- wt$w; e[d] <- 1; e })

  plain <- tulpaObs:::build_integrated_callbacks(model, NULL)$m_step_encode(wt$w)
  for (k in seq_along(y_list)) {
    blk <- plain[[paste0("det", k)]]
    expect_false(is.null(blk$weights))
    expect_equal(length(blk$weights), length(blk$y))
    expect_true(all(blk$n_trials > 0L))
  }

  skip_if_no_tulpamesh()
  sp <- .det_block_spde(fx$coords)
  fielded <- tulpaObs:::build_integrated_callbacks(model, sp)$m_step_encode(wt$w)
  for (k in seq_along(y_list)) {
    blk <- fielded[[paste0("det", k)]]
    expect_false(is.null(blk$weights))
    # One row per source row, matching the broadcast projection.
    expect_equal(length(blk$y), nrow(blk$spatial$A))
    expect_equal(blk$weights, eff[[k]])
    expect_equal(blk$n_trials, as.integer(rowSums(!is.na(y_list[[k]]))))
  }
})

test_that("the spatial Laplace kernel consumes the weights the block hands it", {
  skip_on_cran()
  skip_if_no_tulpamesh()
  # An integer weight is exactly row replication; the same fit through the
  # replicated design must land on the same mode. If the engine ever drops
  # `weights` on the spatial route again, this separates.
  fx <- .det_block_data(n = 80)
  sp <- .det_block_spde(fx$coords, arm = c(TRUE, FALSE))
  set.seed(3)
  X <- cbind(1, fx$data$det_cov)
  n_trials <- rep(5L, fx$n)
  y <- rbinom(fx$n, n_trials, plogis(-0.3 + 0.4 * fx$data$det_cov))
  wt <- sample(1:3, fx$n, replace = TRUE)

  weighted <- tulpa::tulpa_laplace(y = y, n_trials = n_trials, X = X,
                                   spatial = sp$tulpa_spec, family = "binomial",
                                   weights = as.numeric(wt))
  idx <- rep(seq_len(fx$n), times = wt)
  sp_rep <- tulpaObs:::.tobs_spde_broadcast_spec(sp, idx)
  replicated <- tulpa::tulpa_laplace(y = y[idx], n_trials = n_trials[idx],
                                     X = X[idx, , drop = FALSE],
                                     spatial = sp_rep$tulpa_spec,
                                     family = "binomial")
  expect_equal(weighted$mode[1:2], replicated$mode[1:2], tolerance = 1e-8)

  # ... and a fractional weight is not silently ignored.
  unweighted <- tulpa::tulpa_laplace(y = y, n_trials = n_trials, X = X,
                                     spatial = sp$tulpa_spec, family = "binomial")
  set.seed(5)
  fractional <- tulpa::tulpa_laplace(y = y, n_trials = n_trials, X = X,
                                     spatial = sp$tulpa_spec, family = "binomial",
                                     weights = runif(fx$n, 0.05, 1))
  expect_gt(max(abs(fractional$mode[1:2] - unweighted$mode[1:2])), 1e-3)
})
