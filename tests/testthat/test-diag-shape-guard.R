# Argument-shape guards on the diagnostic / pointwise-likelihood
# kernels.
#
# These kernels read their arguments through .begin() pointers and through
# Rcpp's operator(), neither of which bounds-checks, so a mis-shaped argument
# used to read heap bytes into a fitted statistic (a wrong Bayesian p-value or
# PIT residual, no error) or segfault. Every block below builds a consistent
# argument set, checks the kernel accepts it, then perturbs ONE shape and pins
# that the kernel errors naming that argument and the shape it wanted.
#
# The checks live in src/tobs_shape.h and run before any pointer is taken and
# before any OpenMP region, so a mismatch throws from the calling thread.
#
# No model is fit here, so no speed gate applies.

# One argument set with named overrides applied, for calling a kernel with a
# single shape perturbed.
.shape_call <- function(fn, defaults, overrides) {
  if (length(overrides)) defaults[names(overrides)] <- overrides
  do.call(fn, defaults)
}

test_that("cpp_occu_single_ploglik rejects mis-shaped eta_p / y", {
  S <- 3L; N <- 5L; mv <- 2L
  d <- list(eta_psi = matrix(0, S, N), eta_p = matrix(0, S, N),
            y = matrix(0L, N, mv), n_threads = 1L)
  run <- function(...) .shape_call(tulpaObs:::cpp_occu_single_ploglik, d, list(...))

  expect_equal(dim(run()), c(S, N))

  expect_error(run(eta_p = d$eta_p[, -1, drop = FALSE]),
               "eta_p must be [3 x 5]; got [3 x 4].", fixed = TRUE)
  expect_error(run(y = d$y[-1, , drop = FALSE]),
               "y must have 5 rows; got 4.", fixed = TRUE)
})

test_that("cpp_occu_dynamic_ploglik rejects mis-shaped arms / flat arrays", {
  S <- 3L; n_sites <- 4L; mv <- 2L; Tn <- 3L
  E <- matrix(0, S, n_sites)
  d <- list(eta_psi1 = E, eta_p = E, eta_gam = E, eta_eps = E,
            y = integer(n_sites * mv * Tn),
            n_visits = rep(mv, n_sites * Tn),
            any_detected = integer(n_sites * Tn),
            n_sites = n_sites, max_visits = mv, n_seasons = Tn, n_threads = 1L)
  run <- function(...) .shape_call(tulpaObs:::cpp_occu_dynamic_ploglik, d, list(...))

  expect_equal(dim(run()), c(S, n_sites))

  # A detection / colonization / extinction arm with the wrong number of sites.
  expect_error(run(eta_p = E[, -1, drop = FALSE]),
               "eta_p must be [3 x 4]; got [3 x 3].", fixed = TRUE)
  expect_error(run(eta_gam = E[-1, , drop = FALSE]),
               "eta_gam must be [3 x 4]; got [2 x 4].", fixed = TRUE)
  expect_error(run(eta_eps = E[, -1, drop = FALSE]),
               "eta_eps must be [3 x 4]; got [3 x 3].", fixed = TRUE)

  # n_sites is a plain argument, unrelated to the arm matrices until checked.
  expect_error(run(n_sites = n_sites + 1L),
               "eta_psi1 must be [3 x 5]; got [3 x 4].", fixed = TRUE)

  # The flat [n_sites x max_visits x n_seasons] detection array and the two
  # per-(site, season) summaries.
  expect_error(run(y = d$y[-1]), "y must have length 24; got 23.", fixed = TRUE)
  expect_error(run(n_visits = d$n_visits[-1]),
               "n_visits must have length 12; got 11.", fixed = TRUE)
  expect_error(run(any_detected = d$any_detected[-1]),
               "any_detected must have length 12; got 11.", fixed = TRUE)
})

test_that("cpp_occu_integrated_ploglik rejects a mis-sized source slab / counts", {
  S <- 3L; N <- 4L; D <- 2L
  d <- list(eta_psi = matrix(0, S, N), eta_src = numeric(S * N * D),
            K1 = matrix(0L, N, D), K0 = matrix(1L, N, D),
            n_sources = D, n_threads = 1L)
  run <- function(...) .shape_call(tulpaObs:::cpp_occu_integrated_ploglik, d, list(...))

  expect_equal(dim(run()), c(S, N))

  expect_error(run(eta_src = d$eta_src[-1]),
               "eta_src must have length 24; got 23.", fixed = TRUE)
  expect_error(run(n_sources = D + 1L),
               "eta_src must have length 36; got 24.", fixed = TRUE)
  expect_error(run(K1 = matrix(0L, N, D + 1L)),
               "K1 must be [4 x 2]; got [4 x 3].", fixed = TRUE)
  expect_error(run(K0 = d$K0[-1, , drop = FALSE]),
               "K0 must be [4 x 2]; got [3 x 2].", fixed = TRUE)
})

test_that("cpp_cover_pit_cdf bounds pos_col and checks the positive arm", {
  S <- 3L; N <- 5L; n_pos <- 2L
  d <- list(eta_occ = matrix(0, S, N), eta_pos = matrix(0, S, n_pos),
            occur = c(1L, 0L, 1L, 0L, 0L), y_pos = c(0.2, 0.4),
            pos_col = c(1L, 0L, 2L, 0L, 0L), disp = rep(1, S), positive = 0L,
            lower = numeric(0), upper = numeric(0), trunc_upper = numeric(0))
  run <- function(...) .shape_call(tulpaObs:::cpp_cover_pit_cdf, d, list(...))

  expect_equal(dim(run()$cdf_lower), c(S, N))

  # pos_col = 0 at an occupied plot is the j = -1 read: it must name the plot
  # rather than index the column before the first.
  expect_error(run(pos_col = c(0L, 0L, 2L, 0L, 0L)),
               "pos_col[1] = 0 at a plot with occurrence 1", fixed = TRUE)
  expect_error(run(pos_col = c(1L, 0L, 3L, 0L, 0L)),
               "pos_col[3] = 3 at a plot with occurrence 1", fixed = TRUE)

  expect_error(run(eta_pos = d$eta_pos[-1, , drop = FALSE]),
               "eta_pos must have 3 rows; got 2.", fixed = TRUE)
  expect_error(run(disp = d$disp[-1]),
               "disp must have length 3; got 2.", fixed = TRUE)
  expect_error(run(y_pos = d$y_pos[-1]),
               "y_pos must have length 2; got 1.", fixed = TRUE)
  expect_error(run(occur = d$occur[-1]),
               "occur must have length 5; got 4.", fixed = TRUE)

  # Ordinal reads a per-plot class interval; truncated lognormal a ceiling.
  expect_error(run(positive = 2L),
               "lower must have length 2; got 0.", fixed = TRUE)
  expect_error(run(positive = 1L),
               "trunc_upper must have length 2; got 0.", fixed = TRUE)
})

test_that("cpp_cover_ppc rejects a mis-shaped positive arm / dispersion", {
  set.seed(2371)
  S <- 3L; N <- 5L; n_pos <- 2L
  d <- list(eta_occ = matrix(0, S, N), eta_pos = matrix(0, S, n_pos),
            occur = c(1L, 0L, 1L, 0L, 0L), y_pos_nat = c(0.2, 0.4),
            disp = rep(1, S), trunc_upper = numeric(0), positive = 0L,
            freeman = TRUE)
  run <- function(...) .shape_call(tulpaObs:::cpp_cover_ppc, d, list(...))

  expect_length(run()$fit.y, S)

  expect_error(run(eta_pos = matrix(0, S, n_pos + 1L)),
               "eta_pos must be [3 x 2]; got [3 x 3].", fixed = TRUE)
  expect_error(run(eta_pos = d$eta_pos[-1, , drop = FALSE]),
               "eta_pos must be [3 x 2]; got [2 x 2].", fixed = TRUE)
  expect_error(run(disp = d$disp[-1]),
               "disp must have length 3; got 2.", fixed = TRUE)
  expect_error(run(occur = d$occur[-1]),
               "occur must have length 5; got 4.", fixed = TRUE)
  expect_error(run(positive = 1L),
               "trunc_upper must have length 2; got 0.", fixed = TRUE)
})

test_that("cpp_single_ppc checks the site block, the draw matrix and draw_idx", {
  set.seed(2372)
  n_sites <- 6L; mv <- 3L; ndr <- 10L
  y <- matrix(rbinom(n_sites * mv, 1L, 0.4), n_sites, mv)
  storage.mode(y) <- "integer"
  d <- list(X_occ = cbind(1, rnorm(n_sites)), X_det = cbind(1, rnorm(n_sites)),
            draws = matrix(rnorm(ndr * 4L, 0, 0.3), ndr, 4L), draw_idx = 1:5,
            y = y, n_valid = as.integer(rep(mv, n_sites)),
            any_det = as.integer(rowSums(y) > 0), freeman = TRUE)
  run <- function(...) .shape_call(tulpaObs:::cpp_single_ppc, d, list(...))

  expect_length(run()$fit.y, 5L)

  expect_error(run(X_det = d$X_det[-1, , drop = FALSE]),
               "X_det must have 6 rows; got 5.", fixed = TRUE)
  expect_error(run(y = d$y[-1, , drop = FALSE]),
               "y must have 6 rows; got 5.", fixed = TRUE)
  expect_error(run(draws = d$draws[, 1:3, drop = FALSE]),
               "draws must have at least 4 columns; got 3.", fixed = TRUE)
  expect_error(run(draw_idx = c(1L, 2L, ndr + 1L)),
               "draw_idx[3] = 11 is outside [1, 10].", fixed = TRUE)
  expect_error(run(draw_idx = c(0L, 2L)),
               "draw_idx[1] = 0 is outside [1, 10].", fixed = TRUE)
  expect_error(run(n_valid = d$n_valid[-1]),
               "n_valid must have length 6; got 5.", fixed = TRUE)
  expect_error(run(any_det = d$any_det[-1]),
               "any_det must have length 6; got 5.", fixed = TRUE)
})

test_that("cpp_single_pit_cdf checks the site block, the draw matrix and draw_idx", {
  set.seed(2373)
  n_sites <- 6L; mv <- 3L; ndr <- 10L
  y <- matrix(rbinom(n_sites * mv, 1L, 0.4), n_sites, mv)
  storage.mode(y) <- "integer"
  d <- list(X_occ = cbind(1, rnorm(n_sites)), X_det = cbind(1, rnorm(n_sites)),
            draws = matrix(rnorm(ndr * 4L, 0, 0.3), ndr, 4L), draw_idx = 1:5,
            y = y)
  run <- function(...) .shape_call(tulpaObs:::cpp_single_pit_cdf, d, list(...))

  expect_length(run()$cdf_lower, n_sites)

  expect_error(run(X_det = d$X_det[-1, , drop = FALSE]),
               "X_det must have 6 rows; got 5.", fixed = TRUE)
  expect_error(run(y = d$y[-1, , drop = FALSE]),
               "y must have 6 rows; got 5.", fixed = TRUE)
  expect_error(run(draws = d$draws[, 1:3, drop = FALSE]),
               "draws must have at least 4 columns; got 3.", fixed = TRUE)
  expect_error(run(draw_idx = c(2L, ndr + 5L)),
               "draw_idx[2] = 15 is outside [1, 10].", fixed = TRUE)
  # The limits are a posterior mean over the selection, so an empty one has no
  # value to divide by.
  expect_error(run(draw_idx = integer(0)),
               "draw_idx must select at least one draw", fixed = TRUE)
})

test_that("cpp_occu_cover_ppc_agg checks the slabs, pos_site and the CSR units", {
  set.seed(2374)
  n_sites <- 4L; max_v <- 2L; S <- 3L
  y <- matrix(c(1L, 0L), n_sites, max_v)
  d <- list(psi_all = matrix(0.5, n_sites, S),
            p_all = matrix(0.3, n_sites, S * max_v),
            ep_all = matrix(0.1, n_sites, S * max_v),
            y = y, valid = matrix(1L, n_sites, max_v),
            any_det = as.integer(rowSums(y) > 0),
            n_valid = as.integer(rep(max_v, n_sites)), disp = rep(1, S),
            mode_code = 1L, pos_site = c(0L, 2L), yv = c(0.2, 0.4),
            vals_flat = numeric(0), unit_off = 0L, disp2 = 0, positive = 0L,
            freeman = TRUE)
  run <- function(...) .shape_call(tulpaObs:::cpp_occu_cover_ppc_agg, d, list(...))

  expect_length(run()$fit.y, S)

  # pos_site is a raw offset into the per-draw cover-predictor slab.
  expect_error(run(pos_site = c(0L, n_sites)),
               "pos_site[2] = 4 is outside [0, 4).", fixed = TRUE)
  expect_error(run(pos_site = c(-1L, 2L)),
               "pos_site[1] = -1 is outside [0, 4).", fixed = TRUE)

  expect_error(run(p_all = d$p_all[, -1, drop = FALSE]),
               "p_all must be [4 x 6]; got [4 x 5].", fixed = TRUE)
  expect_error(run(ep_all = d$ep_all[-1, , drop = FALSE]),
               "ep_all must be [4 x 6]; got [3 x 6].", fixed = TRUE)
  expect_error(run(valid = d$valid[-1, , drop = FALSE]),
               "valid must have 4 rows; got 3.", fixed = TRUE)
  expect_error(run(y = d$y[, -1, drop = FALSE]),
               "y must be [4 x 2]; got [4 x 1].", fixed = TRUE)
  expect_error(run(disp = d$disp[-1]),
               "disp must have length 3; got 2.", fixed = TRUE)
  expect_error(run(any_det = d$any_det[-1]),
               "any_det must have length 4; got 3.", fixed = TRUE)
  expect_error(run(n_valid = d$n_valid[-1]),
               "n_valid must have length 4; got 3.", fixed = TRUE)
  expect_error(run(yv = d$yv[-1]),
               "yv must have length 2; got 1.", fixed = TRUE)

  # Latent mode reads the detected covers through CSR offsets.
  lat <- list(mode_code = 2L, yv = numeric(0), vals_flat = c(0.2, 0.3, 0.4),
              unit_off = c(0L, 1L, 3L), disp2 = 0.5)
  latent <- function(...) {
    ov <- list(...)
    a <- lat
    if (length(ov)) a[names(ov)] <- ov
    do.call(run, a)
  }
  expect_length(latent()$fit.y, S)

  expect_error(latent(unit_off = c(0L, 1L)),
               "unit_off must have length 3; got 2.", fixed = TRUE)
  expect_error(latent(unit_off = c(1L, 2L, 3L)),
               "unit_off[1] must be 0; got 1.", fixed = TRUE)
  expect_error(latent(unit_off = c(0L, 3L, 1L)),
               "unit_off must be non-decreasing", fixed = TRUE)
  expect_error(latent(vals_flat = c(0.2, 0.3)),
               "vals_flat must have length 3; got 2.", fixed = TRUE)
})

test_that("cpp_occu_cover_ppc checks the per-visit view and the dispersion", {
  set.seed(2375)
  n_sites <- 4L; V <- 6L; S <- 3L
  d <- list(
    X_occ = cbind(1, rnorm(n_sites)),
    X_det_site = cbind(1, rnorm(n_sites)),
    X_pos_site = cbind(1, rnorm(n_sites)),
    X_det_visit = matrix(numeric(0), V, 0L),
    X_pos_visit = matrix(numeric(0), V, 0L),
    site_of_visit = as.integer(c(1, 1, 2, 3, 4, 4)),
    y_det_visit = as.integer(c(1, 0, 1, 0, 1, 1)),
    y_pos_visit = c(0.2, NA, 0.3, NA, 0.4, 0.5),
    b_occ = matrix(rnorm(S * 2L, 0, 0.3), S, 2L),
    b_det = matrix(rnorm(S * 2L, 0, 0.3), S, 2L),
    b_pos = matrix(rnorm(S * 2L, 0, 0.3), S, 2L),
    disp = rep(1, S),
    field_occ = matrix(0, n_sites, S),
    field_pos = matrix(0, n_sites, S),
    off_det = matrix(numeric(0), V, 0L),
    off_pos = matrix(numeric(0), V, 0L),
    any_det = as.integer(c(1, 1, 0, 1)),
    n_valid = as.integer(c(2, 1, 1, 2)),
    positive = 0L, eta_bound = 30, freeman = TRUE)
  run <- function(...) .shape_call(tulpaObs:::cpp_occu_cover_ppc, d, list(...))

  expect_length(run()$fit.y, S)

  expect_error(run(y_det_visit = d$y_det_visit[-1]),
               "y_det_visit must have length 6; got 5.", fixed = TRUE)
  expect_error(run(y_pos_visit = d$y_pos_visit[-1]),
               "y_pos_visit must have length 6; got 5.", fixed = TRUE)
  expect_error(run(any_det = d$any_det[-1]),
               "any_det must have length 4; got 3.", fixed = TRUE)
  expect_error(run(n_valid = d$n_valid[-1]),
               "n_valid must have length 4; got 3.", fixed = TRUE)
  # disp is read once per draw and was the one length-S argument with no check.
  expect_error(run(disp = d$disp[-1]),
               "disp must have length 3; got 2.", fixed = TRUE)
})

# --- kernels reached by #263 ------------------------------------------------
# The blocks above cover what #236 / #237 swept. These cover the entry points
# those sweeps did not reach: a buffer sized from one argument and indexed with
# values taken from another, with nothing relating the two.

test_that("cpp_cover_hurdle_ploglik bounds pos_col like its PIT sibling", {
  S <- 3L; N <- 5L; n_pos <- 2L
  d <- list(eta_occ = matrix(0, S, N), eta_pos = matrix(0, S, n_pos),
            disp = rep(1, S), occur = c(1L, 0L, 1L, 0L, 0L),
            y_pos = c(0.2, 0.4), pos_col = c(1L, 0L, 2L, 0L, 0L),
            positive = 0L, lower = numeric(0), upper = numeric(0),
            trunc_upper = numeric(0), n_threads = 1L)
  run <- function(...) .shape_call(tulpaObs:::cpp_cover_hurdle_ploglik, d, list(...))

  expect_equal(dim(run()), c(S, N))

  # The WAIC path read the column before the first where the PIT path errored.
  expect_error(run(pos_col = c(0L, 0L, 2L, 0L, 0L)),
               "pos_col[1] = 0 at a plot with occurrence 1", fixed = TRUE)
  expect_error(run(pos_col = c(1L, 0L, 3L, 0L, 0L)),
               "pos_col[3] = 3 at a plot with occurrence 1", fixed = TRUE)

  expect_error(run(eta_pos = d$eta_pos[-1, , drop = FALSE]),
               "eta_pos must have 3 rows; got 2.", fixed = TRUE)
  expect_error(run(disp = d$disp[-1]),
               "disp must have length 3; got 2.", fixed = TRUE)
  expect_error(run(y_pos = d$y_pos[-1]),
               "y_pos must have length 2; got 1.", fixed = TRUE)
  expect_error(run(occur = d$occur[-1]),
               "occur must have length 5; got 4.", fixed = TRUE)
  expect_error(run(positive = 2L),
               "lower must have length 2; got 0.", fixed = TRUE)
  expect_error(run(positive = 1L),
               "trunc_upper must have length 2; got 0.", fixed = TRUE)
})

test_that("cpp_occu_cover_ploglik_ragged bounds site_of_visit and the draws", {
  S <- 3L; n_sites <- 4L; V <- 6L
  d <- list(X_occ = matrix(1, n_sites, 1L),
            X_det_site = matrix(1, n_sites, 1L),
            X_pos_site = matrix(1, n_sites, 1L),
            X_det_visit = matrix(0, V, 0L),
            X_pos_visit = matrix(0, V, 0L),
            site_of_visit = c(1L, 1L, 2L, 3L, 4L, 4L),
            y_det_visit = c(1L, 0L, 0L, 1L, 0L, 0L),
            y_pos_visit = c(0.3, 0, 0, 0.5, 0, 0),
            b_occ = matrix(0, S, 1L), b_det = matrix(0, S, 1L),
            b_pos = matrix(0, S, 1L), disp = rep(1, S),
            field_occ = matrix(0, n_sites, S),
            field_pos = matrix(0, n_sites, S),
            off_det = matrix(0, V, 0L), off_pos = matrix(0, V, 0L),
            positive = 0L, eta_bound = 30, n_threads = 1L)
  run <- function(...) .shape_call(tulpaObs:::cpp_occu_cover_ploglik_ragged, d,
                                   list(...))

  expect_equal(dim(run()), c(S, n_sites))

  # site() feeds the per-site accumulators of all three ragged kernels.
  expect_error(run(site_of_visit = c(1L, 1L, 2L, 3L, 4L, 5L)),
               "site_of_visit[6] = 5 is outside [1, 4].", fixed = TRUE)
  expect_error(run(site_of_visit = c(0L, 1L, 2L, 3L, 4L, 4L)),
               "site_of_visit[1] = 0 is outside [1, 4].", fixed = TRUE)
  expect_error(run(site_of_visit = c(NA_integer_, 1L, 2L, 3L, 4L, 4L)),
               "site_of_visit[1] is NA.", fixed = TRUE)

  # b_occ fixes the draw count; the other two arms have to agree with it.
  expect_error(run(b_det = matrix(0, S + 1L, 1L)),
               "b_det must carry the same 3 draws as b_occ; got 4.", fixed = TRUE)
  expect_error(run(b_pos = matrix(0, S - 1L, 1L)),
               "b_pos must carry the same 3 draws as b_occ; got 2.", fixed = TRUE)
  expect_error(run(disp = d$disp[-1]),
               "disp must have length 3; got 2.", fixed = TRUE)
})

test_that("cpp_occu_mscale_cover_ploglik checks its blocks, dims and plot_cell", {
  S <- 3L; n_cells <- 2L; n_plots <- 4L; J <- 2L
  total <- 6L
  d <- list(draws = matrix(0, S, total),
            X_psi = matrix(1, n_cells, 1L),
            X_theta = matrix(1, n_plots, 1L),
            X_p_site = matrix(1, n_plots, 1L),
            X_p_visit = matrix(0, n_plots * J, 0L),
            X_pos_site = matrix(1, n_plots, 1L),
            X_pos_visit = matrix(0, n_plots * J, 0L),
            y = matrix(0L, n_plots, J),
            y_pos = matrix(0, n_plots, J),
            valid = matrix(TRUE, n_plots, J),
            plot_cell = c(1L, 1L, 2L, 2L),
            positive = 0L,
            idx_psi = 0L, p_psi = 1L, idx_theta = 1L, p_theta = 1L,
            idx_p_site = 2L, p_p_site = 1L, idx_p_visit = 0L, p_p_visit = 0L,
            idx_pos_site = 3L, p_pos_site = 1L,
            idx_pos_visit = 0L, p_pos_visit = 0L,
            idx_disp = 4L, n_threads = 1L)
  run <- function(...) .shape_call(tulpaObs:::cpp_occu_mscale_cover_ploglik, d,
                                   list(...))

  expect_equal(dim(run()), c(S, n_cells))

  # A block is an offset and a width packed in R against a third argument.
  expect_error(run(idx_disp = total),
               "disp block [6, 7) does not fit a draw row of 6 columns.",
               fixed = TRUE)
  expect_error(run(p_psi = 7L),
               "psi block [0, 7) does not fit a draw row of 6 columns.",
               fixed = TRUE)
  expect_error(run(idx_theta = -1L),
               "theta block [-1, 0) does not fit a draw row of 6 columns.",
               fixed = TRUE)

  # plot_cell addresses the per-cell accumulator.
  expect_error(run(plot_cell = c(1L, 1L, 2L, 3L)),
               "plot_cell[4] = 3 is outside [1, 2].", fixed = TRUE)
  expect_error(run(plot_cell = c(1L, 1L, 2L, 2L, 2L)),
               "plot_cell must have length 4; got 5.", fixed = TRUE)

  expect_error(run(y = matrix(0L, n_plots + 1L, J)),
               "y must be [4 x 2]; got [5 x 2].", fixed = TRUE)
  expect_error(run(valid = matrix(TRUE, n_plots, J + 1L)),
               "valid must be [4 x 2]; got [4 x 3].", fixed = TRUE)
  # n_cells is X_psi's own row count, so the row half is a tautology; the
  # column half is not -- p_psi is packed separately against the draw row.
  expect_error(run(X_psi = matrix(1, n_cells, 2L)),
               "X_psi must be [2 x 1]; got [2 x 2].", fixed = TRUE)
  expect_error(run(X_theta = matrix(1, n_plots, 2L)),
               "X_theta must be [4 x 1]; got [4 x 2].", fixed = TRUE)
})

test_that("the grouped-RE oracle factories relate their counts to their designs", {
  n_sites <- 4L; n_groups <- 2L; J <- 2L
  n_obs <- n_sites * J
  d <- list(arm = 0L,
            y = rep(1L, n_obs),
            site_idx = rep(seq_len(n_sites), each = J),
            X_lambda = matrix(1, n_sites, 1L),
            X_p = matrix(1, n_obs, 1L),
            Z_site = matrix(1, n_sites, 1L),
            site_group = c(1L, 1L, 2L, 2L),
            n_sites = n_sites, n_groups = n_groups, K_max = 20L)
  run <- function(...) .shape_call(tulpaObs:::cpp_nmix_grouped_oracle, d, list(...))

  expect_true(inherits(run(), "externalptr"))

  # n_sites / n_groups are forwarded from further up, never derived from the
  # designs, so build_common is where the two meet.
  expect_error(run(X_lambda = matrix(1, n_sites + 1L, 1L)),
               "X_lambda must have 4 rows; got 5.", fixed = TRUE)
  expect_error(run(Z_site = matrix(1, n_sites - 1L, 1L)),
               "Z_site must have 4 rows; got 3.", fixed = TRUE)
  expect_error(run(site_group = c(1L, 1L, 2L, 3L)),
               "site_group[4] = 3 is outside [1, 2].", fixed = TRUE)
  # A site whose grouping value is NA used to be dropped from every group,
  # biasing the variance component with no error.
  expect_error(run(site_group = c(1L, 1L, 2L, NA_integer_)),
               "site_group[4] is NA.", fixed = TRUE)
  bad_site <- d$site_idx; bad_site[n_obs] <- 5L
  expect_error(run(site_idx = bad_site),
               "site_idx[8] = 5 is outside [1, 4].", fixed = TRUE)
  expect_error(run(X_p = matrix(1, n_obs - 1L, 1L)),
               "X_p must have 8 rows; got 7.", fixed = TRUE)

  # The removal factory shares the same builder.
  expect_error(do.call(tulpaObs:::cpp_removal_grouped_oracle,
                       utils::modifyList(d, list(site_group = c(1L, 1L, 2L, 3L)))),
               "site_group[4] = 3 is outside [1, 2].", fixed = TRUE)
})

test_that("cpp_nmix_community_oracle bounds its species / site codes", {
  n_sites <- 3L; n_species <- 2L; J <- 2L
  n_obs <- n_sites * n_species * J
  d <- list(y = rep(1L, n_obs),
            site_idx = rep(rep(seq_len(n_sites), each = J), times = n_species),
            species_idx = rep(seq_len(n_species), each = n_sites * J),
            X_lambda = matrix(1, n_sites, 1L),
            X_p = matrix(1, n_obs, 1L),
            n_sites = n_sites, n_species = n_species, K_max = 20L)
  run <- function(...) .shape_call(tulpaObs:::cpp_nmix_community_oracle, d, list(...))

  expect_true(inherits(run(), "externalptr"))

  # rows[species_idx[r] - 1][site_idx[r] - 1].push_back(r) is a heap WRITE
  # through both codes, and neither was checked.
  bad_sp <- d$species_idx; bad_sp[1] <- 3L
  expect_error(run(species_idx = bad_sp),
               "species_idx[1] = 3 is outside [1, 2].", fixed = TRUE)
  bad_site <- d$site_idx; bad_site[2] <- 0L
  expect_error(run(site_idx = bad_site),
               "site_idx[2] = 0 is outside [1, 3].", fixed = TRUE)
  expect_error(run(X_lambda = matrix(1, n_sites + 2L, 1L)),
               "X_lambda must have 3 rows; got 5.", fixed = TRUE)
  expect_error(run(X_p = matrix(1, n_obs + 1L, 1L)),
               "X_p must have 12 rows; got 13.", fixed = TRUE)
})
