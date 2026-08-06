# Batched per-cell pointwise log-likelihood for the three-level multiscale
# occupancy + cover family (occu_multiscale_cover). The R reference
# .occu_ms_cover_nonspatial_ll (per_cell = TRUE) is the oracle; the C++ kernel
# (via .occu_ms_cover_ploglik_core) mirrors it draw for draw over the cell ->
# plot -> visit marginal plus the cover density, and parallelises over draws.
# Byte-close (~1e-14) and thread-count invariant, both cover families and both
# visit-block layouts.

.mk_ms_cover <- function(positive, has_pv, n_cells = 15L, n_plots = 50L, J = 3L,
                         seed = 91) {
  set.seed(seed)
  plot_cell <- sort(sample(seq_len(n_cells), n_plots, replace = TRUE))
  valid <- matrix(runif(n_plots * J) < 0.85, n_plots, J)
  y <- matrix(0L, n_plots, J); ypos <- matrix(0, n_plots, J)
  for (i in seq_len(n_plots)) for (j in seq_len(J)) if (valid[i, j]) {
    y[i, j] <- rbinom(1L, 1L, 0.3)
    if (y[i, j] == 1L)
      ypos[i, j] <- if (positive == "beta") runif(1, 1e-3, 1 - 1e-3) else rlnorm(1)
  }
  X_p_visit   <- if (has_pv) cbind(rnorm(n_plots * J)) else NULL
  X_pos_visit <- if (has_pv) cbind(rnorm(n_plots * J)) else NULL
  p_p   <- 2L + (if (has_pv) 1L else 0L)
  p_pos <- 2L + (if (has_pv) 1L else 0L)
  model <- list(model_type = "occu_multiscale_cover", positive = positive,
                n_cells = n_cells, n_plots = n_plots, max_visits = J,
                plot_cell = plot_cell, valid = valid, y = y, y_pos = ypos,
                X_psi = cbind(1, rnorm(n_cells)), X_theta = cbind(1, rnorm(n_plots)),
                X_p_site = cbind(1, rnorm(n_plots)), X_p_visit = X_p_visit,
                X_pos_site = cbind(1, rnorm(n_plots)), X_pos_visit = X_pos_visit,
                process_info = list(list(p = 2L), list(p = 2L),
                                    list(p = p_p), list(p = p_pos)))
  list(model = model, total = 2L + 2L + p_p + p_pos + 1L)
}

test_that("multiscale occupancy+cover per-cell ploglik: C++ == R oracle", {
  for (positive in c("beta", "lognormal")) for (has_pv in c(FALSE, TRUE)) {
    m <- .mk_ms_cover(positive, has_pv); S <- 30L
    draws <- matrix(rnorm(S * m$total, 0, 0.5), S, m$total)
    idx <- .tobs_occu_ms_cover_nuts_layout(m$model)
    R <- t(apply(draws, 1L, function(par)
      .occu_ms_cover_nonspatial_ll(par, m$model, idx, per_cell = TRUE)))
    C1 <- .occu_ms_cover_ploglik_core(m$model, draws, 1L)
    expect_equal(C1, R, tolerance = 1e-8,
                 info = sprintf("%s pv=%s", positive, has_pv))
    expect_identical(.occu_ms_cover_ploglik_core(m$model, draws, 4L), C1)
  }
})
