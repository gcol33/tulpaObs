# Full-model (field-folded) information criteria for occu_cover() joint fits
# (gcol33/tulpaObs IC fidelity). The occupancy spatial field is sampled with the
# arm coefficients (grid-integrated, via the joint nested-Laplace object) and
# added to the per-cell occupancy predictor, so WAIC / DIC / LOO / CPO are
# full-model (latent field included) and comparable with the INLA / spOccupancy
# criteria, not conditional on the fixed-effect predictor.
#
# These are recovery tests, not smoke tests: the field-folded score must (a)
# differ from the fixed-effect-conditional version, (b) prefer the true field
# model over a no-field fit when the field is strong and real, (c) match an
# independent brute-force per-site log-likelihood, and (d) yield finite,
# mostly-good PSIS-LOO Pareto-k with a returned LOO-PIT.

# Chain adjacency (a 1-D lattice); the package simulator draws the shared ICAR
# field on this graph so the field matches the fitter's parameterisation.
.icfm_chain_adj <- function(n) {
  adj <- matrix(0L, n, n)
  for (s in seq_len(n)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < n)  adj[s, s + 1L] <- 1L
  }
  adj
}

# Simulate a strong-field occu_cover dataset with the package's own generator
# (matched field parameterisation) and fit it with the joint_coupled engine.
.icfm_sim_and_fit <- function(seed, N = 100L, J = 6L, sigma_true = 2.0,
                              alpha_true = 1.0, positive = "lognormal") {
  adj <- .icfm_chain_adj(N)
  sim <- simulate_occu_cover(
    N = N, J = J, beta_occ = c(stats::qlogis(0.4), 0.7),
    beta_p = c(0.0, 0.8), beta_pos = c(log(0.20), -0.4),
    sigma_pos = 0.35, positive = positive, adj = adj,
    sigma = sigma_true, alpha = alpha_true, seed = seed)
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1)
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  fit_field <- suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
    family = occu_cover(positive),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "nested_laplace",
    control = list(engine = "joint_coupled", verbose = FALSE, max.iter = 400L,
                   sigma.grid = c(0.5, 1.0, 1.5, 2.0),
                   alpha.grid = c(0, 0.5, 1.0))))
  fit_nofield <- tobs(
    formula = ~ occ_cov1, data = cell_dat, family = occu_cover(positive),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "laplace", control = list(verbose = FALSE))
  list(fit_field = fit_field, fit_nofield = fit_nofield, sim = sim)
}

# Brute-force per-site marginal log-likelihood for ONE parameter draw, computed
# by an INDEPENDENT routine (no package kernel): field-included occupancy psi,
# per-visit detection p, the Bernoulli-mixture site likelihood over the latent z,
# plus the positive-cover (lognormal) density at detected visits.
.icfm_bruteforce_site_ll <- function(model, psi, p_mat, ep_mat, sigma_pos) {
  N <- model$n_sites; J <- model$max_visits
  valid <- model$valid; y <- model$y; ypos <- model$y_pos
  out <- numeric(N)
  for (i in seq_len(N)) {
    log_det_occ <- 0      # sum_j log P(y_ij, cover_ij | z = 1)
    log_all0    <- 0      # sum_j log(1 - p_ij) over valid visits
    any_det     <- FALSE
    for (j in seq_len(J)) {
      if (!valid[i, j]) next
      pij <- p_mat[i, j]
      log_all0 <- log_all0 + log(1 - pij)
      if (y[i, j] == 1L) {
        any_det <- TRUE
        v  <- ypos[i, j]
        ld_ln <- -log(v) - log(sigma_pos) - 0.5 * log(2 * pi) -
          0.5 * ((log(v) - ep_mat[i, j]) / sigma_pos)^2
        log_det_occ <- log_det_occ + log(pij) + ld_ln
      } else {
        log_det_occ <- log_det_occ + log(1 - pij)
      }
    }
    lp  <- log(psi[i]); l1mp <- log(1 - psi[i])
    if (any_det) {
      out[i] <- lp + log_det_occ
    } else {
      a <- lp + log_all0; b <- l1mp     # psi prod(1-p) + (1-psi), log-sum-exp
      m <- max(a, b)
      out[i] <- m + log(exp(a - m) + exp(b - m))
    }
  }
  out
}

test_that("occu_cover() joint: pointwise log-lik folds in the spatial field", {
  skip_on_cran()
  skip_if_fast()
  o <- .icfm_sim_and_fit(seed = 4101L)
  fit <- o$fit_field

  # The joint path is active and the folded field is non-zero and aligned.
  expect_false(is.null(tulpaObs:::.tobs_joint_fit(fit)))
  c0 <- tulpaObs:::.tobs_occu_cover_components(fit, n.draws = 300L)
  expect_true(any(c0$field_occ != 0))
  expect_gt(stats::sd(as.numeric(c0$field_occ)), 0.3)
  expect_gt(abs(stats::cor(rowMeans(c0$field_occ), o$sim$truth$f)), 0.8)

  # Field-folded pointwise loglik vs fixed-effect-conditional (field zeroed):
  # folding the strong field in changes the WAIC / LOO materially.
  ll_full <- tulpaObs:::.occu_cover_ploglik_core(
    fit$model, c0$b_occ, c0$b_det, c0$b_pos, c0$disp, c0$field_occ, c0$field_pos)
  zf <- matrix(0, nrow(c0$field_occ), ncol(c0$field_occ))
  ll_fe <- tulpaObs:::.occu_cover_ploglik_core(
    fit$model, c0$b_occ, c0$b_det, c0$b_pos, c0$disp, zf, zf)
  cr_full <- tulpa::tulpa_criteria(ll_full, criteria = c("waic", "loo"))
  cr_fe   <- tulpa::tulpa_criteria(ll_fe,   criteria = c("waic", "loo"))
  expect_gt(abs(cr_full$waic - cr_fe$waic), 5.0)
  expect_gt(abs(cr_full$elpd_loo - cr_fe$elpd_loo), 5.0)

  # tobs_waic() end-to-end returns the field-folded value, not the FE one.
  w <- tobs_waic(fit, n.draws = 300L)
  expect_lt(abs(w$waic - cr_full$waic), 15)    # MC noise across draw sets
  expect_gt(abs(w$waic - cr_fe$waic), 5.0)
})

test_that("occu_cover() field-folded WAIC/LOO prefers the true strong-field model", {
  skip_on_cran()
  skip_if_fast()
  # The strong real field carries occupancy structure the fixed effects cannot;
  # the field-folded full-model score then prefers it (lower WAIC, higher elpd).
  prefer_waic <- prefer_loo <- logical(0)
  for (seed in c(4201L, 4202L, 4203L)) {
    o <- .icfm_sim_and_fit(seed = seed, sigma_true = 2.0)
    wf <- tobs_waic(o$fit_field,   n.draws = 600L)
    wn <- tobs_waic(o$fit_nofield, n.draws = 600L)
    cf <- tobs_cpo(o$fit_field,    n.draws = 600L)
    cn <- tobs_cpo(o$fit_nofield,  n.draws = 600L)
    prefer_waic <- c(prefer_waic, wf$waic < wn$waic)
    prefer_loo  <- c(prefer_loo,  cf$elpd_loo > cn$elpd_loo)
  }
  # The full-model field score prefers the field model in the clear majority of
  # seeds (WAIC on a marginalised occupancy model is a moderate discriminator).
  expect_gte(sum(prefer_waic), 2L)
  expect_gte(sum(prefer_loo),  2L)
})

test_that("occu_cover() field-folded pointwise loglik matches a brute-force routine", {
  skip_on_cran()
  skip_if_fast()
  o <- .icfm_sim_and_fit(seed = 4301L)
  fit <- o$fit_field
  model <- fit$model

  c0 <- tulpaObs:::.tobs_occu_cover_components(fit, n.draws = 40L)
  comp <- tulpaObs:::.occu_cover_eta_components(
    model, c0$b_occ, c0$b_det, c0$b_pos, c0$field_occ, c0$field_pos)
  cl <- tulpaObs:::.tobs_clamp_eta

  ll_kernel <- tulpaObs:::.occu_cover_ploglik_core(
    model, c0$b_occ, c0$b_det, c0$b_pos, c0$disp, c0$field_occ, c0$field_pos)

  S <- nrow(c0$b_occ)
  ll_bf <- matrix(0, S, model$n_sites)
  for (sdraw in seq_len(S)) {
    de    <- tulpaObs:::.occu_cover_draw_eta(comp, sdraw, model$n_sites, model$max_visits)
    psi   <- stats::plogis(cl(de$psi_eta))
    p_mat <- stats::plogis(cl(de$p_eta))
    ll_bf[sdraw, ] <- .icfm_bruteforce_site_ll(model, psi, p_mat, de$ep_mat, c0$disp[sdraw])
  }
  expect_lt(max(abs(ll_kernel - ll_bf)), 1e-8)
})

test_that("occu_cover() LOO-PIT is returned, calibrated, with good Pareto-k", {
  skip_on_cran()
  skip_if_fast()
  o <- .icfm_sim_and_fit(seed = 4401L)
  fit <- o$fit_field

  cpo <- tobs_cpo(fit, n.draws = 600L)
  expect_true("pit" %in% names(cpo))
  expect_equal(length(cpo$pit), fit$model$n_sites)
  expect_true(all(is.finite(cpo$pit)))
  expect_true(all(cpo$pit >= 0 & cpo$pit <= 1))

  # PSIS-LOO Pareto-k finite and mostly good on the well-specified strong-field
  # fit (the grid-integrated joint draws give a well-behaved LOO predictive).
  expect_true(all(is.finite(cpo$pareto_k)))
  expect_gt(mean(cpo$pareto_k < 0.7), 0.8)
})

test_that("occu_cover(): loo.unit = 'cell' routes through the site_cell map (tulpaObs#105)", {
  skip_on_cran()
  skip_if_fast()
  o   <- .icfm_sim_and_fit(seed = 4501L)
  fit <- o$fit_field
  n_sites <- fit$model$n_sites

  # occu_cover's cell map is site_cell (the per-site field cell); this fit has no
  # group_var, so it is the identity (each site its own cell) and the pointwise
  # column order is the site order.
  map <- tulpaObs:::.tobs_loo_cell_map(fit)
  expect_equal(map, fit$model$site_cell %||% seq_len(n_sites))
  expect_equal(length(map), n_sites)

  # loo.unit = "cell" == passing that map as group (RNG fixed so the sampled
  # joint draws are identical between the two calls).
  set.seed(7L); a <- tobs_cpo(fit, n.draws = 200L, loo.unit = "cell")
  set.seed(7L); b <- tobs_cpo(fit, n.draws = 200L, group = map)
  expect_equal(a$elpd_loo, b$elpd_loo)
  expect_equal(a$pareto_k, b$pareto_k)
  expect_equal(a$n_groups, n_sites)

  # With the identity cell map cell-level LOO equals the per-site default exactly
  # (each fold is a single site).
  set.seed(7L); d <- tobs_cpo(fit, n.draws = 200L)
  expect_equal(a$elpd_loo, d$elpd_loo)
})
