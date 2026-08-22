# Recovery / diagnostic tests for the adaptive hyperparameter grid wrapped
# around `tulpa_nested_laplace_joint` on the cover() single-block route.
#
# The joint engine integrates the cross-arm coupling as the copy coefficient
# `alpha`: the presence arm sees the BYM2 field at amplitude `sigma`, the cover
# arm at `alpha * sigma`. `control$alpha.grid` is the axis knob for it, declared
# on the positive arm's `field_coef` (`R/cover_hurdle_joint.R`). Every fit below
# pins that axis at `c(0, 0.5, 1.0, 1.5)` with `alpha_true = 1.5`, so the truth
# sits AT its top node -- the configuration the adaptive-grid fix was written
# against (INLAabun Demo 3). A pinned axis is a deliberate pin as far as the
# engine's auto-recenter rescue is concerned, so the placement holds on every
# seed.
#
# A fixed grid carries no cell past its own top node, so at that placement the
# copy axis is truncated exactly where the posterior piles up. The adaptive path
# densifies inside the heaviest boundary cell and adds two outward extension
# points whenever the relative integrand density at the boundary exceeds
# `adaptive_grid_edge_thresh`.
#
# Until the boundary was set on `control$sigma.pos.grid`, an axis of the retired
# (sigma_occ, sigma_pos) parameterization that no cover() route ever read. These
# fits therefore integrated the package default copy axis `c(0, 0.1, 0.234,
# 0.548, 1.282, 3)` with `alpha_true = 1.5` sitting between its 1.282 and 3.0
# nodes, i.e. inside the grid rather than at its edge.
#
# What the boundary costs is the UPPER END of the reported interval, not
# coverage: see the second test for the measurement and for why the coverage
# gap the file used to assert does not open here.

alpha_pin_for_test <- c(0, 0.5, 1.0, 1.5)

simulate_d3_like <- function(seed, alpha_true,
                             N = 300L, n_s = 25L,
                             sigma_true = 0.6, rho_true = 0.7,
                             phi_true = 30,
                             beta_occ = c(-0.3, 0.7),
                             beta_pos = c(0.4, -0.5)) {
  set.seed(seed)
  spatial_idx <- sample.int(n_s, N, replace = TRUE)
  phi_f   <- rnorm(n_s, 0, 1)
  theta_f <- rnorm(n_s, 0, 1)
  w_s     <- sigma_true * (sqrt(rho_true) * phi_f +
                              sqrt(1 - rho_true) * theta_f)
  x <- rnorm(N)
  eta_occ <- beta_occ[1] + beta_occ[2] * x + w_s[spatial_idx]
  occur   <- rbinom(N, 1, plogis(eta_occ))
  eta_pos <- beta_pos[1] + beta_pos[2] * x + alpha_true * w_s[spatial_idx]
  mu_pos  <- plogis(eta_pos)
  y       <- numeric(N)
  is_pos  <- occur == 1L
  y[is_pos] <- rbeta(sum(is_pos), mu_pos[is_pos] * phi_true,
                     (1 - mu_pos[is_pos]) * phi_true)
  y       <- pmin(pmax(y, 0), 1 - 1e-6)
  list(
    data = data.frame(x = x, region = factor(spatial_idx)),
    y    = y,
    truth = list(beta_occ = beta_occ, beta_pos = beta_pos,
                 sigma = sigma_true, rho = rho_true,
                 alpha = alpha_true, phi = phi_true)
  )
}

chain_adj_for_test <- function(n_s) {
  adj <- matrix(0L, n_s, n_s)
  for (s in seq_len(n_s)) {
    for (j in setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L))) adj[s, j] <- 1L
  }
  adj
}

# The copy-axis coordinates the fit actually integrated. `theta_grid` carries
# its axis labels in `theta_names` on this path, so read the column through
# them rather than by position.
alpha_nodes_of <- function(fit) {
  tg <- fit$joint$theta_grid
  cn <- colnames(tg)
  if (is.null(cn)) cn <- fit$joint$theta_names
  sort(unique(as.numeric(tg[, which(cn == "alpha")])))
}

# Posterior weight the fit puts on one copy-axis coordinate.
alpha_node_weight <- function(fit, node) {
  tg <- fit$joint$theta_grid
  cn <- colnames(tg)
  if (is.null(cn)) cn <- fit$joint$theta_names
  a <- as.numeric(tg[, which(cn == "alpha")])
  sum(fit$joint$weights[abs(a - node) < 1e-9])
}

# Axes the adaptive boundary/interior pass fired on. Distinct from the
# var-of-means consistency pass, which reports separately and can place cells
# on the copy axis without the boundary trigger having fired.
adaptive_axes_of <- function(fit) {
  info <- fit$joint$adaptive_grid_info
  if (is.null(info)) return(character(0))
  unlist(strsplit(info$triggered_axes, ","), use.names = FALSE)
}

fit_d3_like <- function(sim, adj, ctrl) {
  tobs(formula = ~ x + bym2(graph = adj, group_var = "region"),
       data = sim$data, family = cover("beta"), y = sim$y,
       method = "nested_laplace", control = ctrl)
}

test_that("adaptive refinement leaves the copy axis at the upper boundary", {
  skip_on_cran()
  skip_if_fast()
  truth_alpha <- 1.5
  n_seeds <- 20L
  n_s     <- 25L
  adj     <- chain_adj_for_test(n_s)

  results <- vapply(seq_len(n_seeds), function(r) {
    sim <- simulate_d3_like(seed = 3400L + r, alpha_true = truth_alpha,
                             n_s = n_s)
    fit <- fit_d3_like(sim, adj, list(
      sigma.grid     = c(0.3, 0.6, 0.9),
      rho.grid       = c(0.5, 0.7, 0.9),
      alpha.grid     = alpha_pin_for_test,
      adaptive.grid  = TRUE
    ))
    expect_s3_class(fit, "cover_fit")
    # The truth is the top node of the pinned axis, so the boundary trigger is
    # what the file is about: the pass must fire ON THE COPY AXIS and place
    # cells strictly beyond the pin. Both hold on 20/20 seeds (tulpa 0.0.163;
    # every seed extends the axis to 3.375). Without them the coverage number
    # below is not evidence about the adaptive path -- the fixed grid covers
    # this truth too, for the reason the next test measures.
    expect_true("alpha" %in% adaptive_axes_of(fit))
    expect_gt(max(alpha_nodes_of(fit)), max(alpha_pin_for_test))
    # Marginalized 95% quantile CI on alpha. Wald mean +/- 1.96*sd over-rejects
    # in this regime because the alpha posterior is right-skewed once the copy
    # amplitude piles against the axis edge. Quantile CI is the engine's
    # documented summary for skew axes.
    alpha_lo <- fit$joint$theta_ci_lo["alpha"]
    alpha_hi <- fit$joint$theta_ci_hi["alpha"]
    alpha_lo <= truth_alpha && truth_alpha <= alpha_hi
  }, logical(1))

  coverage <- mean(results)
  # Nominal 0.95; measured 20/20 across these seeds. 0.85 is a conservative
  # lower bound robust to seed-level noise. Read alongside the two structural
  # assertions above: with the truth ON the top node the upper side of a
  # quantile CI is a weak instrument, so coverage alone does not separate the
  # adaptive path from the fixed one here.
  expect_gte(coverage, 0.85)
})

test_that("the fixed grid's upper CI edge is its own axis geometry, the adaptive one is not", {
  skip_on_cran()
  skip_if_fast()
  truth_alpha <- 1.5
  n_seeds <- 20L
  n_s     <- 25L
  adj     <- chain_adj_for_test(n_s)

  # Pin the 13-point phi axis. On the 7-point package default the adaptive pass
  # also densifies `phi_pos`, so the two arms would differ on two axes at once;
  # at 13 points the pass fires on the copy axis alone (measured: triggered
  # axes are exactly "alpha" on 20/20 adaptive seeds here), which is what makes
  # the paired comparison below a statement about the copy axis.
  ctrl <- list(
    sigma.grid     = c(0.3, 0.6, 0.9),
    rho.grid       = c(0.5, 0.7, 0.9),
    alpha.grid     = alpha_pin_for_test,
    phi.grid       = exp(seq(log(2), log(300), length.out = 13))
  )

  cover_fixed <- logical(n_seeds)
  cover_adapt <- logical(n_seeds)
  hi_fixed    <- numeric(n_seeds)
  hi_adapt    <- numeric(n_seeds)
  for (r in seq_len(n_seeds)) {
    sim <- simulate_d3_like(seed = 3400L + r, alpha_true = truth_alpha,
                             n_s = n_s)
    fit_fix <- fit_d3_like(sim, adj, c(ctrl, list(adaptive.grid = FALSE)))
    fit_ad  <- fit_d3_like(sim, adj, c(ctrl, list(adaptive.grid = TRUE)))

    # The fixed arm integrates the pin and nothing else, and puts most of its
    # posterior weight on the top node: the truncation this fixture exists to
    # create. Measured over these seeds: weight on 1.5 is 0.876 on average,
    # 0.206 at its lowest.
    expect_equal(alpha_nodes_of(fit_fix), alpha_pin_for_test)
    expect_gt(alpha_node_weight(fit_fix, max(alpha_pin_for_test)), 0.10)

    hi_fixed[r] <- fit_fix$joint$theta_ci_hi["alpha"]
    hi_adapt[r] <- fit_ad$joint$theta_ci_hi["alpha"]
    # The adaptive arm reports a strictly higher upper edge on every seed
    # (measured ratio >= 1.05 over 20/20).
    expect_gt(hi_adapt[r], hi_fixed[r])

    cover_fixed[r] <- fit_fix$joint$theta_ci_lo["alpha"] <= truth_alpha &&
      truth_alpha <= fit_fix$joint$theta_ci_hi["alpha"]
    cover_adapt[r] <- fit_ad$joint$theta_ci_lo["alpha"] <= truth_alpha &&
      truth_alpha <= fit_ad$joint$theta_ci_hi["alpha"]
  }

  # What separates the two reads at the boundary is where the upper edge comes
  # from. The fixed arm's 97.5% point is set by the axis's own outer-cell
  # geometry: it lands at essentially the same place on every seed (measured
  # mean 1.733, SD 0.011, full range 1.689-1.738 over 20 seeds). The adaptive
  # arm places real cells past the pin and its upper edge follows the data
  # (mean 2.270, SD 0.373, range 1.825-3.226). Assert the separation in spread
  # rather than a fixed threshold, so the check does not encode this fixture's
  # absolute scale.
  expect_lt(sd(hi_fixed), 0.10)
  expect_gt(sd(hi_adapt), 5 * sd(hi_fixed))

  # Coverage does NOT separate the two arms at this placement, and the file no
  # longer claims it does. Both arms cover 20/20 here. With the truth exactly ON
  # the top node the upper side of the fixed arm's interval cannot miss: since a
  # grid's outer cell contributes its own half-width to the reported support, so
  # the fixed read reaches past 1.5 by construction, and before that change the
  # quantile clamped AT 1.5 and `truth <= ci_hi` held as an equality. The 0.53
  # -> 0.83 gain recorded when the adaptive path landed was measured on a Wald
  # summary, which this file replaced. All that is defended here is
  # non-regression.
  expect_gte(mean(cover_adapt) - mean(cover_fixed), -0.15)
})

test_that("adaptive refinement is a no-op on the copy axis when the integrand has decayed", {
  skip_on_cran()
  skip_if_fast()
  # alpha_true = 0 places the truth at the OPPOSITE edge of the pinned copy
  # axis. The data drives the posterior to concentrate near 0 with quickly
  # decaying tails, so the integrand at the upper boundary (alpha = 1.5) is
  # negligible and the boundary trigger must not fire on that axis. The min
  # side is where the mode sits, but the axis is bounded below at 0 and every
  # proposed extension point is clipped away, so nothing is added there either.
  sim <- simulate_d3_like(seed = 3101L, alpha_true = 0.0)
  n_s <- nlevels(sim$data$region)
  adj <- chain_adj_for_test(n_s)
  fit <- fit_d3_like(sim, adj, list(
    sigma.grid     = c(0.3, 0.6, 0.9),
    rho.grid       = c(0.5, 0.7, 0.9),
    alpha.grid     = alpha_pin_for_test,
    adaptive.grid  = TRUE
  ))

  # The adaptive pass may still fire on `phi_pos` (the 7-point default phi axis
  # is coarse enough to trigger interior densification), and the var-of-means
  # consistency pass -- a separate mechanism, reported separately -- does place
  # cells on the copy axis here. Neither is the boundary trigger this test is
  # about, so assert on the triggered-axes list rather than on the node set.
  expect_false("alpha" %in% adaptive_axes_of(fit))

  # The alpha posterior genuinely concentrates at the truth rather than
  # spreading over the axis, which spans to 1.5. Measured: CI [0, 0.449],
  # median 0. A Wald z-score over-rejects here for the same reason as above.
  expect_lt(fit$joint$theta_ci_lo["alpha"], 0.05)
  expect_lt(fit$joint$theta_ci_hi["alpha"], 0.60)
})

test_that("outer-grid pruning keeps the mode and leaves estimates unchanged", {
  skip_on_cran()
  skip_if_fast()
  # the dense outer tensor concentrates posterior mass on a few cells (ESS ~
  # 1), so most cells run a full-data inner solve for negligible weight. The
  # cheap-pass prune skips them. This test asserts pruning is a no-op on the
  # inference: the modal hyperparameters and the coefficient estimates must
  # match the un-pruned dense-grid fit, and the pruned cells must carry
  # negligible weight (the mode is never pruned). Nothing here depends on where
  # the truth sits relative to a node; the copy axis is pinned so the grid this
  # runs on is stated rather than inherited from the package default.
  n_s <- 25L; adj <- chain_adj_for_test(n_s)
  sim <- simulate_d3_like(seed = 101L, alpha_true = 1.0, N = 400L, n_s = n_s)
  ctrl_grid <- list(
    sigma.grid     = exp(seq(log(0.2), log(1.5), length.out = 5)),
    rho.grid       = c(0.25, 0.5, 0.7, 0.9),
    alpha.grid     = alpha_pin_for_test,
    phi.grid       = exp(seq(log(2), log(300), length.out = 7)),
    adaptive.grid  = FALSE
  )
  run <- function(prune)
    fit_d3_like(sim, adj, c(ctrl_grid, list(prune = prune, prune.tol = 1e-4)))

  f_off <- run(FALSE)
  f_on  <- run(TRUE)

  # Most cells pruned; the prune machinery is engaged.
  expect_true(f_on$joint$prune_n_pruned > 0.5 * f_on$joint$n_grid)

  # The modal hyperparameters (highest-weight cell) agree between fits.
  mode_off <- f_off$joint$theta_grid[which.max(f_off$joint$weights), ]
  mode_on  <- f_on$joint$theta_grid[which.max(f_on$joint$weights), ]
  names(mode_off) <- f_off$joint$theta_names
  names(mode_on)  <- f_on$joint$theta_names
  for (nm in c("sigma", "rho", "alpha"))
    expect_equal(unname(mode_on[nm]), unname(mode_off[nm]), tolerance = 0.05)

  # Coefficient estimates materially unchanged by pruning.
  expect_equal(f_on$beta_occ, f_off$beta_occ, tolerance = 0.02)
  expect_equal(f_on$beta_pos, f_off$beta_pos, tolerance = 0.02)
  expect_equal(unname(f_on$joint$theta_mean["alpha"]),
               unname(f_off$joint$theta_mean["alpha"]), tolerance = 0.05)
})
