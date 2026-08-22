# =============================================================================
# test-ms-occu-cover.R - community (multispecies) joint occupancy-detection +
# cover (ms_occu_cover()).
#
# Per-species coefficient RE with Gaussian community covariances across the
# psi / p / cover arms; the latent presence z integrates out in closed form
# (the occu_cover marginal) and the per-species deviations are integrated by a
# Laplace-EM with Louis/Schur community-mean SEs (read at the natural scale).
#
# These cover:
#   - family constructor + capability gates
#   - community-mean point recovery + per-species coefficient recovery
#   - community-mean 95% CI coverage over seeds
#   - the S3 method surface
#   - a beta-arm recovery smoke
#
# The family is status = "working": the community-MEAN 95% CIs are gated at the
# 0.85 working floor of the recovery rubric (measured pooled coverage ~0.92 at 24
# seeds). The community MEANS are at the natural scale and unbiased; the
# binary-detection community VARIANCE component carries the documented Laplace
# small-cluster attenuation -- a lower bound debiased by AGHQ up to the
# re.aghq.maxdim cap, and explicitly tested as a lower bound above it.
# =============================================================================


# Build the visit-level tobs_data() detection design shared across species.
.msoc_visits <- function(N, J, visit_data) {
  long <- data.frame(
    site_id = rep(seq_len(N), each = J),
    visit   = rep(seq_len(J), times = N),
    yy      = 0L,
    det_cov1 = visit_data$det_cov1,
    pos_cov1 = visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "yy", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  od$det.covs
}


test_that("ms_occu_cover() constructor returns a tobs_family", {
  f <- ms_occu_cover()
  expect_s3_class(f, "tobs_family")
  expect_equal(f$name, "ms_occu_cover")
  expect_equal(f$status, "working")
  expect_equal(f$params$positive, "beta")
  expect_equal(ms_occu_cover("lognormal")$params$positive, "lognormal")
})


test_that("ms_occu_cover() enforces its capability gates", {
  sim <- simulate_ms_occu_cover(n_species = 4, N = 30, J = 3,
                                positive = "lognormal", seed = 1)
  vis <- .msoc_visits(30, 3, sim$visit_data)

  # Missing species.
  expect_error(
    tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
         detection = ~ det_cov1, positive = ~ pos_cov1,
         y = sim$y, y_pos = sim$y_pos, visits = vis, method = "laplace"),
    "species"
  )
  # Missing y_pos.
  expect_error(
    tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
         detection = ~ det_cov1, positive = ~ pos_cov1,
         y = sim$y, visits = vis, species = sim$species, method = "laplace"),
    "y_pos"
  )
  # nested_laplace is not offered (community coupled-field engine is upstream-
  # pending).
  expect_error(
    tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
         detection = ~ det_cov1, positive = ~ pos_cov1,
         y = sim$y, y_pos = sim$y_pos, visits = vis, species = sim$species,
         method = "nested_laplace"),
    "not available"
  )
  # A structured term on the detection arm is rejected: shared fields are
  # supported on the occupancy and cover arms only.
  adj <- matrix(0L, 30, 30)
  for (i in seq_len(29)) adj[i, i + 1L] <- adj[i + 1L, i] <- 1L
  expect_error(
    tobs(~ occ_cov1, data = sim$data,
         family = ms_occu_cover("lognormal"),
         detection = ~ det_cov1 + icar(graph = adj),
         positive = ~ pos_cov1, y = sim$y, y_pos = sim$y_pos, visits = vis,
         species = sim$species, method = "laplace"),
    "detection arm must use a plain formula"
  )
})


test_that("ms_occu_cover() recovers community means + per-species coefs (lognormal)", {
  skip_on_cran()
  skip_if_fast()
  set.seed(21)
  sim <- simulate_ms_occu_cover(
    n_species = 14, N = 90, J = 4,
    mu_occ = c(stats::qlogis(0.45), 0.7), mu_p = c(0.2, -0.4),
    mu_pos = c(log(0.12), 0.5), sd_occ = 0.5, sd_p = 0.4, sd_pos = 0.4,
    positive = "lognormal", sigma_pos = 0.4, seed = 21)
  vis <- .msoc_visits(90, 4, sim$visit_data)

  fit <- tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
              detection = ~ det_cov1, positive = ~ pos_cov1,
              y = sim$y, y_pos = sim$y_pos, visits = vis, species = sim$species,
              method = "laplace", control = list(verbose = FALSE))

  expect_true(isTRUE(fit$convergence$converged))

  # Community means: within ~2.5 SE of truth on every coordinate (beta arms).
  beta_means <- fit$means[seq_len(6)]
  truth <- c(sim$truth$mu_occ, sim$truth$mu_p, sim$truth$mu_pos)
  z <- abs(beta_means - truth) / fit$sds[seq_len(6)]
  expect_true(all(z < 2.5))

  # Shared dispersion recovered.
  expect_lt(abs(exp(fit$means[["log_sigma_pos"]]) - 0.4), 0.1)

  # Per-species coefficients: high correlation with simulated truth on the
  # occupancy and cover arms; detection moderate (binary-data RE attenuation).
  # Occupancy BLUPs are detection-filtered and shrunk toward the community
  # mean, so they recover less sharply than the count-informed coefficients of
  # ms_abun (correct Bayesian shrinkage, not bias).
  cm <- fit$ms_community
  expect_gt(min(diag(cor(cm$coef_occ, sim$truth$beta_occ))), 0.78)
  expect_gt(min(diag(cor(cm$coef_pos, sim$truth$beta_pos))), 0.80)
  expect_gt(min(diag(cor(cm$coef_p,   sim$truth$beta_p))),   0.50)

  # Community SDs positive and finite; the occupancy / cover arms track the
  # realized spread, the detection arm is attenuated but non-degenerate.
  expect_true(all(cm$sd_occ > 0 & cm$sd_pos > 0 & cm$sd_p > 0))
  emp_occ <- apply(sim$truth$beta_occ, 2, sd)
  emp_pos <- apply(sim$truth$beta_pos, 2, sd)
  expect_true(all(cm$sd_occ > 0.35 * emp_occ & cm$sd_occ < 1.8 * emp_occ))
  expect_true(all(cm$sd_pos > 0.35 * emp_pos & cm$sd_pos < 1.8 * emp_pos))
})


test_that("ms_occu_cover() community-mean 95% CIs cover near the nominal rate", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 20L
  covered <- logical(0)
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_occu_cover(
      n_species = 12, N = 70, J = 4,
      mu_occ = c(stats::qlogis(0.45), 0.6), mu_p = c(0.2, -0.4),
      mu_pos = c(log(0.12), 0.4), sd_occ = 0.5, sd_p = 0.4, sd_pos = 0.4,
      positive = "lognormal", sigma_pos = 0.4, seed = 300 + s)
    vis <- .msoc_visits(70, 4, sim$visit_data)
    fit <- tryCatch(
      tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
           detection = ~ det_cov1, positive = ~ pos_cov1,
           y = sim$y, y_pos = sim$y_pos, visits = vis, species = sim$species,
           method = "laplace", control = list(verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    bm    <- fit$means[seq_len(6)]
    bsd   <- fit$sds[seq_len(6)]
    truth <- c(sim$truth$mu_occ, sim$truth$mu_p, sim$truth$mu_pos)
    covered <- c(covered, abs(bm - truth) < 1.96 * bsd)
  }
  # Working-family gate: pooled over the six community-mean coefficients x 20
  # seeds at the 0.85 floor. Measured pooled coverage ~0.92.
  expect_gt(mean(covered), 0.85)
})


test_that("ms_occu_cover() S3 methods work", {
  skip_on_cran()
  skip_if_fast()
  set.seed(5)
  sim <- simulate_ms_occu_cover(n_species = 8, N = 45, J = 3,
                                positive = "lognormal", seed = 5)
  vis <- .msoc_visits(45, 3, sim$visit_data)
  fit <- tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
              detection = ~ det_cov1, positive = ~ pos_cov1,
              y = sim$y, y_pos = sim$y_pos, visits = vis, species = sim$species,
              method = "laplace", control = list(verbose = FALSE))

  expect_s3_class(fit, "tobs_fit")
  expect_no_error(print(fit))
  expect_no_error(summary(fit))

  # coef(): per-process list (community means by arm).
  cf <- coef(fit)
  expect_true(is.list(cf))
  expect_setequal(names(cf), c("psi", "p", "pos"))
  expect_setequal(names(cf$psi), c("(Intercept)", "occ_cov1"))

  V <- vcov(fit)
  expect_equal(nrow(V), length(fit$means))
  expect_equal(nrow(confint(fit)), length(fit$means))

  re <- ranef(fit)
  expect_s3_class(re, "data.frame")
  # n_species x (p_occ + p_p + p_pos) = 8 x (2 + 2 + 2)
  expect_equal(nrow(re), 8L * 6L)
  expect_true(all(c("species", "arm", "term", "estimate") %in% names(re)))
  expect_setequal(unique(re$arm), c("psi", "p", "pos"))

  fv <- fitted(fit)
  expect_equal(dim(fv$psi),   c(45L, 8L))
  expect_equal(dim(fv$p),     c(45L, 8L))
  expect_equal(dim(fv$cover), c(45L, 8L))
  expect_true(all(fv$psi > 0 & fv$psi < 1))
  expect_true(all(fv$cover > 0))

  ys <- simulate(fit, nsim = 1)
  expect_equal(dim(ys$y),     c(45L, 3L, 8L))
  expect_equal(dim(ys$y_pos), c(45L, 3L, 8L))

  expect_type(nobs(fit), "integer")
})


test_that("ms_occu_cover() flags community-variance Laplace attenuation", {
  skip_on_cran()
  skip_if_fast()
  set.seed(5)
  sim <- simulate_ms_occu_cover(n_species = 8, N = 45, J = 3,
                                positive = "lognormal", seed = 5)
  vis <- .msoc_visits(45, 3, sim$visit_data)
  fit <- tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
              detection = ~ det_cov1, positive = ~ pos_cov1,
              y = sim$y, y_pos = sim$y_pos, visits = vis, species = sim$species,
              method = "laplace", control = list(verbose = FALSE))

  # Machine-readable marker on the community block: variance attenuated, means not.
  va <- fit$ms_community$var_attenuation
  expect_type(va, "list")
  expect_false(va$means_affected)
  expect_identical(va$source, "laplace_small_cluster")
  expect_identical(va$debias, "none")
  expect_true(all(c("Sigma_occ", "Sigma_p", "Sigma_pos") %in% va$affects))

  # print() surfaces the caveat so the reported community variance is not read
  # as unbiased.
  out <- paste(utils::capture.output(print(fit)), collapse = "\n")
  expect_match(out, "community VARIANCE|Community variance|variance components",
               ignore.case = TRUE)
})


test_that("ms_occu_cover() recovers community means (beta arm, smoke)", {
  skip_on_cran()
  skip_if_fast()
  n_seeds <- 8L
  occ_x <- p_x <- pos_x <- phi_e <- rep(NA_real_, n_seeds)
  occ_dev <- p_dev <- pos_dev <- rep(NA_real_, n_seeds)
  for (s in seq_len(n_seeds)) {
    sim <- simulate_ms_occu_cover(
      n_species = 12, N = 70, J = 4,
      mu_occ = c(stats::qlogis(0.5), 0.6), mu_p = c(0.2, -0.4),
      mu_pos = c(stats::qlogis(0.3), 0.3), sd_occ = 0.5, sd_p = 0.4,
      sd_pos = 0.3, positive = "beta", phi = 25, seed = 500 + s)
    vis <- .msoc_visits(70, 4, sim$visit_data)
    fit <- tryCatch(
      tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("beta"),
           detection = ~ det_cov1, positive = ~ pos_cov1,
           y = sim$y, y_pos = sim$y_pos, visits = vis, species = sim$species,
           method = "laplace", control = list(verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    occ_x[s] <- fit$means[["psi_occ_cov1"]]
    p_x[s]   <- fit$means[["p_det_cov1"]]
    pos_x[s] <- fit$means[["pos_pos_cov1"]]
    phi_e[s] <- exp(fit$means[["log_phi"]])
    occ_dev[s] <- occ_x[s] - colMeans(sim$truth$beta_occ)[2L]
    p_dev[s]   <- p_x[s]   - colMeans(sim$truth$beta_p)[2L]
    pos_dev[s] <- pos_x[s] - colMeans(sim$truth$beta_pos)[2L]
  }
  # Community-mean slopes against the seed's REALIZED mean (colMeans of
  # truth$beta_occ / beta_p / beta_pos), not the nominal 0.6 / -0.4 / 0.3
  # (#155): this loop already averaged over seeds before comparing to the
  # nominal, so calls it statistically valid as it stood; retargeting is a
  # power improvement, not a bug fix. Budget is 5x the SE of an 8-seed mean,
  # from a fresh 16-seed measurement of this exact fixture (occ sd 0.092 ->
  # SE_8 0.033 -> budget 0.17; p sd 0.061 -> SE_8 0.021 -> budget 0.11; pos sd
  # 0.015 -> SE_8 0.0053 -> budget 0.03), all well below the old flat 0.25.
  # The occupancy and detection arms carry real per-seed outliers (occ
  # deviation -0.208 on one seed, p deviation +0.161 on another, of 16
  # measured) consistent with this family's documented binary-data RE
  # attenuation; none of the three means is significantly biased over the 16
  # seeds (largest |mean| / se is occ at 0.66), so the budget is not absorbing
  # a known one-sided shift.
  expect_lt(abs(mean(occ_dev, na.rm = TRUE)), 0.17)
  expect_lt(abs(mean(p_dev,   na.rm = TRUE)), 0.11)
  expect_lt(abs(mean(pos_dev, na.rm = TRUE)), 0.03)
  expect_gt(mean(phi_e, na.rm = TRUE), 8)   # well above the boundary
})

test_that("ms_occu_cover() recovers community means + sigma_pos (gaussian, #127)", {
  skip_on_cran()
  skip_if_fast()
  # Identity-Gaussian positive arm (the delta-normal hurdle): mu = eta, an
  # unbounded magnitude, residual SD sigma_pos. Same community Laplace-EM as the
  # lognormal arm (only the positive density differs), so the SAME recovery
  # config and bar as the lognormal test; the gaussian residual SD recovers to
  # ~0.39 (matched to lognormal at 0.392 in dev_notes/_diag_127_sigma.R). The
  # per-species pos deviation only inflates the residual under a sparser-than-
  # this cover design, exactly as it would on the log scale.
  set.seed(31)
  sim <- simulate_ms_occu_cover(
    n_species = 14, N = 90, J = 4,
    mu_occ = c(stats::qlogis(0.45), 0.7), mu_p = c(0.2, -0.4),
    mu_pos = c(2.0, 0.5), sd_occ = 0.5, sd_p = 0.4, sd_pos = 0.4,
    positive = "gaussian", sigma_pos = 0.4, seed = 31)
  vis <- .msoc_visits(90, 4, sim$visit_data)
  fit <- tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("gaussian"),
              detection = ~ det_cov1, positive = ~ pos_cov1,
              y = sim$y, y_pos = sim$y_pos, visits = vis, species = sim$species,
              method = "laplace", control = list(verbose = FALSE))
  expect_true(isTRUE(fit$convergence$converged))

  # Community means within ~2.5 SE (beta arms), same bar as lognormal.
  beta_means <- fit$means[seq_len(6)]
  truth <- c(sim$truth$mu_occ, sim$truth$mu_p, sim$truth$mu_pos)
  expect_true(all(abs(beta_means - truth) / fit$sds[seq_len(6)] < 2.5))

  # Shared residual SD recovered (identity-Gaussian sigma).
  expect_lt(abs(exp(fit$means[["log_sigma_pos"]]) - 0.4), 0.1)

  # Per-species coefficients correlate with truth. The cover arm (what this
  # gaussian test validates) recovers sharply; the occupancy BLUPs are
  # detection-filtered and shrunk, so a looser, seed-robust floor (as the family
  # note documents -- they recover less sharply than count-informed coefs).
  cm <- fit$ms_community
  expect_gt(min(diag(cor(cm$coef_occ, sim$truth$beta_occ))), 0.70)
  expect_gt(min(diag(cor(cm$coef_pos, sim$truth$beta_pos))), 0.80)
})

test_that("ms_occu_cover() has WAIC / DIC / CPO (per-species cell marginal, #116)", {
  sim <- simulate_ms_occu_cover(n_species = 5, N = 40, J = 3,
                                mu_pos = c(log(0.12), 0.4), positive = "lognormal",
                                sigma_pos = 0.4, seed = 7)
  vis <- .msoc_visits(40, 3, sim$visit_data)
  fit <- tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
              detection = ~ det_cov1, positive = ~ pos_cov1,
              y = sim$y, y_pos = sim$y_pos, visits = vis, species = sim$species,
              method = "laplace", control = list(verbose = FALSE))
  # Pointwise ll is per-(species, cell): [n_draws x (n_species * n_sites)].
  ll <- .tobs_pointwise_loglik(fit, n.draws = 100L)
  expect_equal(dim(ll), c(100L, 5L * 40L))
  expect_true(all(is.finite(ll)))
  w <- waic(fit, n.draws = 100L)
  d <- dic(fit, n.draws = 100L)
  cp <- cpo(fit, n.draws = 100L)
  expect_true(is.finite(w$waic) && w$p_waic > 0)
  expect_true(is.finite(d$dic))
  expect_true(is.finite(cp$lpml))
  # WAIC and the DIC plug-in agree to within a few units on a well-constrained fit.
  expect_lt(abs(w$waic - d$dic), 0.05 * abs(w$waic))
})

test_that("ms_occu_cover(\"gaussian\") WAIC uses the gaussian density (#116/#127)", {
  sim <- simulate_ms_occu_cover(n_species = 5, N = 40, J = 3,
                                mu_pos = c(2.0, 0.4), positive = "gaussian",
                                sigma_pos = 0.4, seed = 9)
  vis <- .msoc_visits(40, 3, sim$visit_data)
  fit <- tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("gaussian"),
              detection = ~ det_cov1, positive = ~ pos_cov1,
              y = sim$y, y_pos = sim$y_pos, visits = vis, species = sim$species,
              method = "laplace", control = list(verbose = FALSE))
  w <- waic(fit, n.draws = 100L)
  expect_true(is.finite(w$waic))
})

test_that("ms_occu_cover(\"gaussian\") simulate() round-trips (#127)", {
  sim <- simulate_ms_occu_cover(n_species = 5, N = 40, J = 3,
                                mu_pos = c(2.0, 0.4), positive = "gaussian",
                                sigma_pos = 0.4, seed = 11)
  expect_true(any(sim$y_pos[!is.na(sim$y_pos)] < 0) ||
              min(sim$y_pos, na.rm = TRUE) < 1)   # gaussian admits low / negative
  expect_true(all(is.na(sim$y_pos[!is.na(sim$y) & sim$y == 0L])))
  expect_equal(sim$truth$positive, "gaussian")
})

test_that("ms_occu_cover() AGHQ debias reduces variance-component attenuation (#56)", {
  skip_on_cran()
  skip_if_fast()

  # Intercept-only per arm -> total RE dim P = 3, so the AGHQ debias is active by
  # default. Small per-species n is the attenuation regime the debias targets.
  sd_occ_t <- 0.7; sd_p_t <- 0.6
  n_seeds <- 6L
  em <- aghq <- matrix(NA_real_, n_seeds, 2L)   # cols: sd_occ, sd_p
  for (s in seq_len(n_seeds)) {
    sim <- simulate_ms_occu_cover(
      n_species = 16, N = 40, J = 3, n_occ_covs = 0, n_det_covs = 0,
      n_pos_covs = 0, sd_occ = sd_occ_t, sd_p = sd_p_t, sd_pos = 0.4,
      positive = "lognormal", sigma_pos = 0.4, seed = 700 + s)
    # Intercept-only arms carry no per-visit covariates, so no `visits` design.
    common <- list(formula = ~ 1, data = sim$data,
                   family = ms_occu_cover("lognormal"),
                   detection = ~ 1, positive = ~ 1, y = sim$y, y_pos = sim$y_pos,
                   species = sim$species, method = "laplace")
    fit_em <- tryCatch(do.call(tobs, c(common, list(
      control = list(verbose = FALSE, re.aghq = FALSE)))), error = function(e) NULL)
    fit_ag <- tryCatch(do.call(tobs, c(common, list(
      control = list(verbose = FALSE)))), error = function(e) NULL)
    if (is.null(fit_em) || is.null(fit_ag)) next
    em[s, ]   <- c(fit_em$ms_community$sd_occ[1], fit_em$ms_community$sd_p[1])
    aghq[s, ] <- c(fit_ag$ms_community$sd_occ[1], fit_ag$ms_community$sd_p[1])
    if (s == 1L) {
      expect_identical(fit_ag$ms_community$var_attenuation$debias, "aghq")
      expect_identical(fit_em$ms_community$var_attenuation$debias, "none")
    }
  }
  ok <- stats::complete.cases(em) & stats::complete.cases(aghq)
  expect_gte(sum(ok), 4L)

  truth <- c(sd_occ_t, sd_p_t)
  em_mean   <- colMeans(em[ok, , drop = FALSE])
  aghq_mean <- colMeans(aghq[ok, , drop = FALSE])
  # The EM (Laplace) variance components are attenuated (biased low); the AGHQ
  # debias inflates them toward the truth. Check on both binary arms.
  expect_true(all(aghq_mean > em_mean),
              info = paste("em:", paste(round(em_mean, 3), collapse = " "),
                           "aghq:", paste(round(aghq_mean, 3), collapse = " ")))
  expect_true(all(abs(aghq_mean - truth) <= abs(em_mean - truth) + 1e-6),
              info = paste("em:", paste(round(em_mean, 3), collapse = " "),
                           "aghq:", paste(round(aghq_mean, 3), collapse = " "),
                           "truth:", paste(truth, collapse = " ")))
})


test_that("ms_occu_cover() variance debias is a hard cap; EM is a tested lower bound (#98)", {
  skip_on_cran()
  skip_if_fast()

  # One covariate per arm -> total RE dim P = 6, above the re.aghq.maxdim cap (4).
  # The tensor AGHQ debias is exponential in P, so it is a hard scope limit: above
  # the cap the EM (Laplace) variance is reported as a documented lower bound,
  # unchanged by re.aghq. This pins both halves of that contract.
  sd_occ_t <- 0.7; sd_p_t <- 0.6; n_seeds <- 5L
  sd_default <- sd_noaghq <- matrix(NA_real_, n_seeds, 2L)   # cols: sd_occ, sd_p
  for (s in seq_len(n_seeds)) {
    sim <- simulate_ms_occu_cover(
      n_species = 16, N = 60, J = 4,
      mu_occ = c(stats::qlogis(0.45), 0.5), mu_p = c(0.2, -0.4),
      mu_pos = c(log(0.12), 0.3), sd_occ = sd_occ_t, sd_p = sd_p_t, sd_pos = 0.4,
      positive = "lognormal", sigma_pos = 0.4, seed = 720 + s)
    vis <- .msoc_visits(60, 4, sim$visit_data)
    common <- list(formula = ~ occ_cov1, data = sim$data,
                   family = ms_occu_cover("lognormal"), detection = ~ det_cov1,
                   positive = ~ pos_cov1, y = sim$y, y_pos = sim$y_pos,
                   visits = vis, species = sim$species, method = "laplace")
    fd <- tryCatch(do.call(tobs, c(common, list(
      control = list(verbose = FALSE)))), error = function(e) NULL)
    fn <- tryCatch(do.call(tobs, c(common, list(
      control = list(verbose = FALSE, re.aghq = FALSE)))), error = function(e) NULL)
    if (is.null(fd) || is.null(fn)) next
    if (s == 1L) {
      # Above the cap the debias is out of scope, flagged "none" even though
      # re.aghq is TRUE by default, and the affected components are named.
      expect_identical(fd$ms_community$var_attenuation$debias, "none")
      expect_true(all(c("Sigma_occ", "Sigma_p", "Sigma_pos") %in%
                        fd$ms_community$var_attenuation$affects))
    }
    sd_default[s, ] <- c(fd$ms_community$sd_occ[1], fd$ms_community$sd_p[1])
    sd_noaghq[s, ]  <- c(fn$ms_community$sd_occ[1], fn$ms_community$sd_p[1])
  }
  ok <- stats::complete.cases(sd_default) & stats::complete.cases(sd_noaghq)
  expect_gte(sum(ok), 3L)

  # Above the cap re.aghq has no effect: the EM variance is exactly what is
  # reported (the debias never runs to inflate it).
  expect_equal(sd_default[ok, , drop = FALSE], sd_noaghq[ok, , drop = FALSE],
               tolerance = 1e-6)

  # The reported community SD is a lower bound on the truth -- attenuated, never
  # over-estimating -- the documented behaviour above the cap.
  truth <- c(sd_occ_t, sd_p_t)
  em_mean <- colMeans(sd_default[ok, , drop = FALSE])
  expect_true(all(em_mean <= truth + 0.05),
              info = paste("em:", paste(round(em_mean, 3), collapse = " "),
                           "truth:", paste(truth, collapse = " ")))
})
