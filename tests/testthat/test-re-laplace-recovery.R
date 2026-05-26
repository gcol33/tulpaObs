# Parameter-recovery tests for formula random effects fit under the
# DETERMINISTIC engine (gcol33/tulpaObs#11). The default Laplace engine fits iid
# intercept RE, uncorrelated random slopes, and correlated slopes via a
# variance-component EM (R/em_laplace_re.R) instead of silently dropping them.
# The raw EM integrates the RE block by Laplace (the glmer nAGQ=1 regime, not
# Breslow-Clayton PQL), which attenuates sigma / the RE correlation for binary
# data; by default an adaptive Gauss-Hermite pass (R/re_aghq.R) debiases them.
# The sigma tolerances here are still generous (single-seed noise + detection
# thinning), with the calibrated check against NUTS and a multi-seed AGHQ-vs-EM
# bias check below.

sim_occu_re_intercept <- function(seed = 101, ng = 30L, per = 25L, J = 6L,
                                  b0 = 0.3, b1 = -0.6, sigma = 0.9, p = 0.45) {
  set.seed(seed)
  N <- ng * per
  g <- rep(seq_len(ng), each = per)
  x <- rnorm(N)
  b_true <- rnorm(ng, 0, sigma)
  z <- rbinom(N, 1, plogis(b0 + b1 * x + b_true[g]))
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, z[i] * p)
  list(y = y, d = data.frame(g = factor(g), x = x), b_true = b_true)
}

# Single-season occupancy with a CORRELATED random intercept + slope:
# (b0_g, b1_g) ~ N(0, Sigma) per group, eta = b0 + b0_g + (b1 + b1_g) * x.
sim_occu_re_corr <- function(seed = 404, ng = 40L, per = 25L, J = 6L,
                             b0 = 0.2, b1 = -0.4, p = 0.5,
                             Sigma = matrix(c(0.6, 0.3, 0.3, 0.4), 2, 2)) {
  set.seed(seed)
  N <- ng * per
  g <- rep(seq_len(ng), each = per)
  x <- rnorm(N)
  U <- matrix(rnorm(ng * 2L), ng, 2L) %*% chol(Sigma)  # rows ~ N(0, Sigma)
  eta <- b0 + U[g, 1] + (b1 + U[g, 2]) * x
  z <- rbinom(N, 1, plogis(eta))
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, z[i] * p)
  list(y = y, d = data.frame(g = factor(g), x = x), U = U, Sigma = Sigma)
}

test_that("iid intercept RE is fit (not dropped) by the default Laplace engine", {
  s <- sim_occu_re_intercept()
  fit <- tobs(~ x + (1 | g), data = s$d, y = s$y, detection = ~ 1,
              family = occu(), method = "laplace",
              control = list(verbose = FALSE))

  expect_identical(fit$method, "laplace")
  # The random effect is present, not silently dropped: a sigma hyperparameter
  # and per-group BLUPs are reported.
  sig_nm <- grep("^sigma_", names(fit$means), value = TRUE)
  expect_length(sig_nm, 1L)
  expect_true(is.finite(fit$means[[sig_nm]]) && fit$means[[sig_nm]] > 0.4 &&
              fit$means[[sig_nm]] < 1.4)

  cf <- coef(fit)$psi
  expect_lt(abs(cf[["(Intercept)"]] - 0.3), 0.3)
  expect_lt(abs(cf[["x"]] - (-0.6)), 0.3)
  expect_lt(abs(plogis(fit$means[["p_(Intercept)"]]) - 0.45), 0.08)

  # Occupancy fixed-effect SE comes from the natural-scale observed info, not
  # the M-inflated M-step Hessian (which would be ~sqrt(M)=~31x too small).
  se_int <- fit$sds[["psi_(Intercept)"]]
  expect_true(is.finite(se_int) && se_int > 0.03 && se_int < 1)

  re <- ranef(fit)
  expect_s3_class(re, "data.frame")
  expect_equal(nrow(re), 30L)
  expect_gt(cor(re$estimate, s$b_true), 0.7)
})

test_that("uncorrelated random slopes (1 + x || g) recover under Laplace", {
  set.seed(202)
  ng <- 30L; per <- 30L; N <- ng * per; J <- 6L
  g <- rep(seq_len(ng), each = per); x <- rnorm(N)
  st <- c(0.7, 0.5)
  b0 <- rnorm(ng, 0, st[1]); b1 <- rnorm(ng, 0, st[2])
  z <- rbinom(N, 1, plogis(0.2 + b0[g] + (-0.4 + b1[g]) * x))
  y <- matrix(0L, N, J); for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, z[i] * 0.5)
  d <- data.frame(g = factor(g), x = x)

  fit <- tobs(~ x + (1 + x || g), data = d, y = y, detection = ~ 1,
              family = occu(), method = "laplace", control = list(verbose = FALSE))

  sig_nm <- grep("^sigma_", names(fit$means), value = TRUE)
  expect_length(sig_nm, 2L)   # intercept + slope sigma
  sig_hat <- fit$means[sig_nm]
  expect_true(all(is.finite(sig_hat) & sig_hat > 0.2 & sig_hat < 1.4))

  re <- ranef(fit)
  expect_equal(sort(unique(re$term)), c("(Intercept)", "x"))
  b0_hat <- re$estimate[re$term == "(Intercept)"]
  b1_hat <- re$estimate[re$term == "x"]
  expect_gt(cor(b0_hat, b0), 0.6)
  expect_gt(cor(b1_hat, b1), 0.4)
})

test_that("deterministic sigma and BLUPs track the NUTS fit on the same data", {
  skip_on_cran()
  skip_if_fast()
  s <- sim_occu_re_intercept(seed = 7, ng = 30L, per = 25L)

  fit_l <- tobs(~ x + (1 | g), data = s$d, y = s$y, detection = ~ 1,
                family = occu(), method = "laplace", control = list(verbose = FALSE))
  fit_n <- tobs(~ x + (1 | g), data = s$d, y = s$y, detection = ~ 1,
                family = occu(), method = "nuts",
                control = list(n.iter = 400, n.warmup = 200, seed = 1, verbose = FALSE))

  sig_l <- fit_l$means[[grep("^sigma_", names(fit_l$means), value = TRUE)]]
  sig_n <- exp(fit_n$means[[grep("^log_sigma_", names(fit_n$means), value = TRUE)]])
  # Deterministic Laplace sigma is in the NUTS ballpark (the Laplace
  # small-cluster bias is modest at this cluster size; both within ~60%).
  expect_lt(abs(sig_l - sig_n) / sig_n, 0.6)

  re_l <- ranef(fit_l); re_n <- ranef(fit_n)
  expect_gt(cor(re_l$estimate, re_n$estimate), 0.85)
})

test_that("correlated random slopes (1 + x | g) recover under Laplace", {
  skip_on_cran()
  skip_if_fast()
  s <- sim_occu_re_corr(seed = 402)

  fit <- tobs(~ x + (1 + x | g), data = s$d, y = s$y, detection = ~ 1,
              family = occu(), method = "laplace", control = list(verbose = FALSE))

  expect_identical(fit$method, "laplace")

  # intercept + slope sigma, plus the off-diagonal reported as a correlation
  # (the block is genuinely full, not forced diagonal like `(1 + x || g)`).
  sig_nm <- grep("^sigma_", names(fit$means), value = TRUE)
  expect_length(sig_nm, 2L)
  sig_hat <- fit$means[sig_nm]
  expect_true(all(is.finite(sig_hat) & sig_hat > 0.2 & sig_hat < 1.6))

  cor_nm <- grep("^cor_", names(fit$means), value = TRUE)
  expect_length(cor_nm, 1L)
  rho_hat <- fit$means[[cor_nm]]
  # Simulated correlation is +0.61; recover the sign / a positive association
  # (the Laplace approximation attenuates the magnitude, so the bound is generous).
  expect_gt(rho_hat, 0.2)
  expect_lt(rho_hat, 0.99)

  cf <- coef(fit)$psi
  expect_lt(abs(cf[["(Intercept)"]] - 0.2), 0.3)
  expect_lt(abs(cf[["x"]] - (-0.4)), 0.3)

  re <- ranef(fit)
  expect_equal(sort(unique(re$term)), c("(Intercept)", "x"))
  b0_hat <- re$estimate[re$term == "(Intercept)"]
  b1_hat <- re$estimate[re$term == "x"]
  expect_gt(cor(b0_hat, s$U[, 1]), 0.6)
  expect_gt(cor(b1_hat, s$U[, 2]), 0.4)
})

test_that("RE forms the deterministic engine cannot fit error toward NUTS", {
  s <- sim_occu_re_intercept(seed = 3, ng = 12L, per = 12L)

  # Non-single model types route RE to NUTS (validator guard, exercised
  # directly to avoid building a full community dataset here).
  re_int <- tulpaObs:::.tobs_term_re(group = s$d$g, type = "intercept")
  stub <- list(model_type = "community", data = s$d, X_det_visit = NULL)
  expect_error(
    tulpaObs:::.validate_re_laplace(re_int, stub, NULL, "gaussian_laplace"),
    "single|nuts")

  # RE on the detection predictor alone is now supported (own RE block), so the
  # validator accepts it.
  re_det <- tulpaObs:::.tobs_term_re(group = s$d$g, type = "intercept")
  re_det$shared <- c(FALSE, TRUE)
  stub2 <- list(model_type = "single", data = s$d, X_det_visit = NULL)
  expect_silent(
    tulpaObs:::.validate_re_laplace(re_det, stub2, NULL, "gaussian_laplace"))

  # A single RE shared across BOTH predictors stays NUTS-only (each arm fits its
  # own block on the deterministic path, not one shared realization).
  re_both <- tulpaObs:::.tobs_term_re(group = s$d$g, type = "intercept")
  re_both$shared <- c(TRUE, TRUE)
  expect_error(
    tulpaObs:::.validate_re_laplace(re_both, stub2, NULL, "gaussian_laplace"),
    "shared|nuts")

  # RE + visit-level detection covariates also stays NUTS-only.
  stub3 <- list(model_type = "single", data = s$d,
                X_det_visit = matrix(0, nrow(s$d), 1L))
  expect_error(
    tulpaObs:::.validate_re_laplace(re_det, stub3, NULL, "gaussian_laplace"),
    "visit|nuts")
})

test_that("AGHQ variance-component debias runs by default and is toggleable", {
  s <- sim_occu_re_intercept(seed = 21, ng = 25L, per = 12L)
  args <- list(formula = ~ x + (1 | g), data = s$d, y = s$y, detection = ~ 1,
               family = occu(), method = "laplace")

  fit_on  <- do.call(tobs, c(args, list(control = list(verbose = FALSE))))
  fit_off <- do.call(tobs, c(args, list(control = list(re.aghq = FALSE,
                                                       verbose = FALSE))))

  # On by default; status surfaced on the fit.
  expect_true(isTRUE(fit_on$aghq$applied))
  expect_identical(fit_on$aghq$n_quad, 9L)
  expect_false(isTRUE(fit_off$aghq$applied))

  # The refine moves the estimate (the two fits are not identical) and keeps a
  # finite, positive sigma. (The quadrature engine itself lives in tulpa and is
  # tested there -- tulpa/tests/testthat/test-re-aghq.R.)
  sig_on  <- fit_on$means[[grep("^sigma_", names(fit_on$means), value = TRUE)]]
  sig_off <- fit_off$means[[grep("^sigma_", names(fit_off$means), value = TRUE)]]
  expect_true(is.finite(sig_on) && sig_on > 0)
  expect_false(isTRUE(all.equal(sig_on, sig_off)))
})

test_that("AGHQ removes the small-cluster sigma attenuation (multi-seed)", {
  skip_on_cran()
  skip_if_fast()
  # True sigma = 0.9; per-group n = 8 is small enough that the Laplace EM
  # attenuates sigma. AGHQ should land closer to truth on the seed average.
  truth <- 0.9
  sig <- function(fit) fit$means[[grep("^sigma_", names(fit$means), value = TRUE)]]
  one <- function(seed, aghq) {
    s <- sim_occu_re_intercept(seed = seed, ng = 30L, per = 8L, sigma = truth)
    fit <- tobs(~ x + (1 | g), data = s$d, y = s$y, detection = ~ 1,
                family = occu(), method = "laplace",
                control = list(re.aghq = aghq, verbose = FALSE))
    sig(fit)
  }
  seeds <- 1:8
  em   <- vapply(seeds, one, numeric(1), aghq = FALSE)
  aghq <- vapply(seeds, one, numeric(1), aghq = TRUE)

  # EM attenuates (mean below truth); AGHQ corrects upward and lands closer.
  expect_lt(mean(em), truth)
  expect_lt(abs(mean(aghq) - truth), abs(mean(em) - truth))
  expect_lt(abs(mean(aghq) - truth), 0.12)
})

test_that("LKJ regularization keeps the RE correlation off the +-1 boundary", {
  skip_on_cran()
  skip_if_fast()
  # Weakly-identified correlated slope (per-group n = 12), true rho = +0.61.
  # Unregularized ML (re.lkj = 1) over-estimates rho and can hit +-1; the
  # default LKJ pulls it off the boundary while staying near-unbiased.
  rho_of <- function(f) f$means[[grep("^cor_", names(f$means), value = TRUE)]]
  one <- function(seed, eta) {
    s <- sim_occu_re_corr(seed = seed, ng = 40L, per = 12L)
    f <- tobs(~ x + (1 + x | g), data = s$d, y = s$y, detection = ~ 1,
              family = occu(), method = "laplace",
              control = list(n.quad = 7L, re.lkj = eta, verbose = FALSE))
    list(rho = rho_of(f), eta = f$aghq$lkj_eta)
  }
  seeds <- 401:410
  off  <- vapply(seeds, function(s) one(s, 1)$rho,   numeric(1))   # ML, no prior
  reg  <- lapply(seeds, function(s) one(s, 1.5))                   # default
  reg_rho <- vapply(reg, `[[`, numeric(1), "rho")

  # Status surfaced; default eta is 1.5.
  expect_equal(reg[[1]]$eta, 1.5)
  # The default keeps every fit strictly inside the boundary; the unregularized
  # ML reaches it on at least one seed in this weak regime.
  expect_true(all(reg_rho < 0.97))
  expect_gt(max(off), max(reg_rho))
  # ... while staying near the truth (0.61) on the seed average.
  expect_lt(abs(mean(reg_rho) - 0.612), 0.12)
})

# Single-season occupancy with a DETECTION random intercept (gcol33/tulpaObs#11
# follow-up): detection p = sigmoid(d0 + b_obs[observer]), b_obs ~ N(0, sigma);
# occupancy psi = sigmoid(b0 + b1 occ_cov). The random effect enters the
# detection predictor, so the AGHQ refine integrates b through p (not psi).
sim_det_re_intercept <- function(seed = 1, N = 400L, J = 6L, ng = 40L,
                                 b0 = 0.4, b1 = -0.7, d0 = 0.2, sigma = 0.8) {
  set.seed(seed)
  occ_cov  <- rnorm(N)
  observer <- sample.int(ng, N, replace = TRUE)
  b_obs    <- rnorm(ng, 0, sigma)
  z <- rbinom(N, 1, plogis(b0 + b1 * occ_cov))
  p <- plogis(d0 + b_obs[observer])
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) if (z[i] == 1L) y[i, ] <- rbinom(J, 1L, p[i])
  list(y = y, d = data.frame(occ_cov = occ_cov, observer = factor(observer)),
       b_obs = b_obs, sigma = sigma, d0 = d0, b1 = b1)
}

test_that("a detection random intercept is fit on its own arm (AGHQ arm = det)", {
  s <- sim_det_re_intercept(seed = 1)
  fit <- tobs(~ occ_cov, detection = ~ (1 | observer), family = occu(),
              data = s$d, y = s$y, method = "laplace",
              control = list(verbose = FALSE))

  expect_identical(fit$method, "laplace")
  # The RE landed on the detection process: the sigma hyperparameter is named
  # for the detection arm (sigma_p<t>), not the occupancy arm (sigma_g<t>).
  sig_nm <- grep("^sigma_", names(fit$means), value = TRUE)
  expect_length(sig_nm, 1L)
  expect_match(sig_nm, "^sigma_p")
  expect_true(is.finite(fit$means[[sig_nm]]) && fit$means[[sig_nm]] > 0)

  # AGHQ ran on the detection arm.
  expect_true(isTRUE(fit$aghq$applied))
  expect_identical(fit$aghq$arm, "det")

  # Per-group detection BLUPs track the simulated observer effects.
  re <- ranef(fit)
  rp <- re[re$group == "p1", ]
  expect_equal(nrow(rp), 40L)
  expect_gt(cor(rp$estimate, s$b_obs), 0.6)

  # Fixed effects recover.
  expect_lt(abs(plogis(fit$means[["p_(Intercept)"]]) - plogis(s$d0)), 0.1)
})

test_that("AGHQ removes the detection-RE sigma attenuation (multi-seed)", {
  skip_on_cran()
  skip_if_fast()
  # Detection RE is only informed by occupied sites, so the raw nAGQ=1 EM
  # attenuates sigma severely; the AGHQ refine on the exact p-marginal restores
  # it. True sigma = 0.8.
  truth <- 0.8
  sig <- function(fit) fit$means[[grep("^sigma_", names(fit$means), value = TRUE)]]
  em <- aghq <- numeric(6L)
  for (k in 1:6) {
    s <- sim_det_re_intercept(seed = 100L + k)
    args <- list(formula = ~ occ_cov, detection = ~ (1 | observer),
                 family = occu(), data = s$d, y = s$y, method = "laplace")
    em[k]   <- sig(do.call(tobs, c(args, list(control = list(re.aghq = FALSE, verbose = FALSE)))))
    aghq[k] <- sig(do.call(tobs, c(args, list(control = list(verbose = FALSE)))))
  }
  # AGHQ is closer to truth than the (heavily attenuated) EM, and near-unbiased.
  expect_lt(abs(mean(aghq) - truth), abs(mean(em) - truth))
  expect_lt(abs(mean(aghq) - truth), 0.15)
})
