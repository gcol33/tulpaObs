# =============================================================================
# test-occu-cover-joint-coupled.R - end-to-end gates for the joint-coupled
# engine wired through the occu_cover_lognormal cell-coupling spec
# (gcol33/tulpa#32 Layer B.2 consumer + R-facing fit wiring).
#
# Routed via `tobs(method = "nested_laplace", control = list(engine =
# "joint_coupled"))`. Compares against the v3 nested-Laplace path that
# profiles z out per outer (alpha, sigma) candidate.
# =============================================================================


test_that("joint_coupled smoke fit runs end-to-end and returns finite betas", {
  N <- 30L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim <- simulate_occu_cover(
    N = N, J = J, positive = "lognormal",
    adj = adj, sigma = 0.8, alpha = 1.0, seed = 12345L
  )
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                   det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  fit <- suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = 500L,
                   engine = "joint_coupled")
  ))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(attr(fit, "tobs_family")$name, "occu_cover")
  expect_true(all(is.finite(fit$means[c("psi_(Intercept)", "p_(Intercept)",
                                          "pos_(Intercept)")])))
  # Law-of-total-variance SDs: var-of-means (across outer grid) + within-cell
  # Laplace via constraint-corrected diag(Q_k^-1). Without the sum-to-zero
  # correction on the ICAR field, the intercept SDs blow up (~50-80 on the
  # smoke fit) because Q has a near-null direction along (intercept,
  # mean(phi)); the constraint pulls them back into a sensible range.
  intercept_sds <- fit$sds[c("psi_(Intercept)", "p_(Intercept)",
                              "pos_(Intercept)")]
  expect_true(all(is.finite(intercept_sds)))
  expect_true(all(intercept_sds > 0))
  expect_true(all(intercept_sds < 5))
  # The joint_fit object is attached for diagnostics; carries the outer
  # (sigma, alpha) grid + per-cell log_marginal.
  expect_s3_class(fit$joint_fit, "tulpa_nested_laplace_joint")
  expect_equal(length(fit$spatial_field), N)
  expect_true(all(is.finite(fit$spatial_field)))
})


test_that("joint_coupled smoke fit runs end-to-end under beta positive arm", {
  N <- 30L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim <- simulate_occu_cover(
    N = N, J = J, positive = "beta", phi = 25,
    adj = adj, sigma = 0.8, alpha = 1.0, seed = 6789L
  )
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                   det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  fit <- suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
    family = occu_cover("beta"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = 500L,
                   engine = "joint_coupled")
  ))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(attr(fit, "tobs_family")$name, "occu_cover")
  expect_true(all(is.finite(fit$means[c("psi_(Intercept)", "p_(Intercept)",
                                          "pos_(Intercept)")])))
  intercept_sds <- fit$sds[c("psi_(Intercept)", "p_(Intercept)",
                              "pos_(Intercept)")]
  expect_true(all(is.finite(intercept_sds)))
  expect_true(all(intercept_sds > 0))
  expect_true(all(intercept_sds < 5))
  expect_s3_class(fit$joint_fit, "tulpa_nested_laplace_joint")
  expect_equal(length(fit$spatial_field), N)
  expect_true(all(is.finite(fit$spatial_field)))
})


test_that("joint_coupled returns the joint betas+field posterior covariance", {
  N <- 30L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim <- simulate_occu_cover(
    N = N, J = J, positive = "lognormal",
    adj = adj, sigma = 0.8, alpha = 1.0, seed = 12345L
  )
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                   det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  fit <- suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = 500L,
                   engine = "joint_coupled")
  ))

  pi_list <- fit$process_info
  p_beta <- length(pi_list[[1L]]$coef_names) +
            length(pi_list[[2L]]$coef_names) +
            length(pi_list[[3L]]$coef_names)
  n_joint <- p_beta + N

  Vj <- fit$joint_vcov
  expect_true(is.matrix(Vj) && is.numeric(Vj))
  expect_equal(dim(Vj), c(n_joint, n_joint))
  expect_true(all(is.finite(Vj)))
  expect_true(isSymmetric(Vj))

  # Materially non-zero beta x field cross-block proves it is not marginal-only.
  cross <- Vj[seq_len(p_beta), p_beta + seq_len(N), drop = FALSE]
  expect_gt(max(abs(cross)), 1e-8)

  # The parameter-surface vcov's beta block is no longer diagonal.
  Vbeta <- fit$vcov[seq_len(p_beta), seq_len(p_beta), drop = FALSE]
  off_beta <- Vbeta[upper.tri(Vbeta)]
  expect_gt(max(abs(off_beta)), 0)

  expect_equal(length(fit$joint_means), n_joint)
  expect_true(all(is.finite(fit$joint_means)))
})


test_that("joint_coupled parameter-surface vcov carries beta x hyper cross-cov (tulpaObs#46)", {
  N <- 40L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim <- simulate_occu_cover(
    N = N, J = J, positive = "lognormal",
    adj = adj, sigma = 0.8, alpha = 1.0, seed = 12345L
  )
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                   det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  fit <- suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = 500L, engine = "joint_coupled")
  ))

  V  <- fit$vcov
  nm <- colnames(V)
  beta_nm  <- grep("^(psi_|p_|pos_)", nm, value = TRUE)
  hyper_nm <- setdiff(nm, beta_nm)
  # sigma + alpha are surfaced as hyperparameters on this single-field fit.
  expect_true(all(c("sigma", "alpha") %in% hyper_nm))

  # The beta-hyper block is no longer dropped to zero: the law-of-total-covariance
  # between term ties betas to the grid coordinates (within-cell zero).
  cross <- V[beta_nm, hyper_nm, drop = FALSE]
  expect_gt(max(abs(cross)), 1e-6)

  # The hyper-hyper block carries the sigma-alpha covariance, not just a diagonal.
  hyhy <- V[hyper_nm, hyper_nm, drop = FALSE]
  expect_gt(abs(hyhy["sigma", "alpha"]), 1e-6)
  # Its diagonal still equals the reported marginal SDs squared.
  expect_equal(diag(hyhy), fit$sds[hyper_nm]^2, tolerance = 1e-8)

  # V stays symmetric and PSD so the posterior draws are well-defined.
  expect_true(isSymmetric(unname(V)))
  expect_true(min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) > -1e-8)
  expect_no_error(chol(V))
})


test_that("joint_coupled errors on non-spatial occu_cover", {
  N <- 30L; J <- 4L
  sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal", seed = 1L)
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                   det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  # No spatial term on the psi formula -> nested_laplace itself errors before
  # the engine pick fires.
  expect_error(
    tobs(formula = ~ occ_cov1, data = cell_dat,
         family = occu_cover("lognormal"),
         detection = ~ det_cov1, positive = ~ pos_cov1,
         y = od$y, y_pos = y_pos, visits = od$det.covs,
         method = "nested_laplace",
         control = list(verbose = FALSE, engine = "joint_coupled")),
    "spatial term"
  )
})


test_that("joint_coupled regularises the cover (pos) intercept by default (tulpaObs#32)", {
  # The pos arm sees the shared field only at detected visits, so its intercept
  # confounds with the field level over those cells -- a direction the
  # sum-to-zero field constraint does not pin when low-occupancy regions carry
  # no cover. Left unpenalised (the engine's flat 1e-4 ridge, sd ~100) the cover
  # intercept floats to a huge posterior SD while the regularised occupancy
  # intercept stays tight, blowing up predict()'s conditional cover via Jensen.
  # The default must hand the pos arm a finite-precision intercept prior, like
  # the psi / p arms get from occu_priors().
  responses <- list(
    psi = list(X = matrix(1, 3, 1, dimnames = list(NULL, "(Intercept)"))),
    p   = list(X = matrix(1, 3, 1, dimnames = list(NULL, "(Intercept)"))),
    pos = list(X = cbind(`(Intercept)` = rep(1, 3), x = rnorm(3)))
  )

  # Default priors -> all three arms regularised.
  ap <- tulpaObs:::.occu_cover_coupled_arm_priors(NULL, responses)
  expect_false(is.null(ap$pos))
  expect_true(all(is.finite(ap$pos$prec)))
  expect_true(all(ap$pos$prec > 1e-4))           # tighter than the flat engine ridge
  expect_false(is.null(ap$psi))
  expect_false(is.null(ap$p))

  # A custom occu_priors() still leaves the cover arm regularised by the
  # cover_priors() default.
  ap2 <- tulpaObs:::.occu_cover_coupled_arm_priors(occu_priors(), responses)
  expect_false(is.null(ap2$pos))
  expect_true(all(is.finite(ap2$pos$prec)))

  # priors = FALSE / "none" disables every arm (back to the flat ridge).
  apN <- tulpaObs:::.occu_cover_coupled_arm_priors(FALSE, responses)
  expect_null(apN$pos); expect_null(apN$psi); expect_null(apN$p)
  apN2 <- tulpaObs:::.occu_cover_coupled_arm_priors("none", responses)
  expect_null(apN2$pos)

  # An explicit cover_priors() narrows the cover arm.
  ap3 <- tulpaObs:::.occu_cover_coupled_arm_priors(
    cover_priors(pos_intercept = list(mean = 0, sd = 0.5)), responses)
  expect_false(is.null(ap3$pos))
  expect_gt(ap3$pos$prec[1], ap$pos$prec[1])      # sd 0.5 -> higher precision than default sd 3
})


test_that("joint_coupled recovers slopes, hypers, field shape (10 seeds)", {
  skip_on_cran()
  skip_if_fast()

  n_seeds <- 18L
  N <- 100L; J <- 6L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }

  beta_occ_truth <- c(stats::qlogis(0.4),  0.7)
  beta_p_truth   <- c(0.0,                  0.8)
  beta_pos_truth <- c(log(0.20),           -0.4)
  sigma_pos_truth <- 0.35
  sigma_truth   <- 1.0
  alpha_truth   <- 1.0

  est_psi_x <- est_p_x <- est_pos_x <- est_alpha <- est_sigma <-
    field_cor <- se_psi_x <- se_p_x <- se_pos_x <- rep(NA_real_, n_seeds)

  for (s in seq_len(n_seeds)) {
    sim <- simulate_occu_cover(
      N = N, J = J, beta_occ = beta_occ_truth, beta_p = beta_p_truth,
      beta_pos = beta_pos_truth, sigma_pos = sigma_pos_truth,
      positive = "lognormal", adj = adj,
      sigma = sigma_truth, alpha = alpha_truth, seed = 5000L + s
    )
    long <- data.frame(
      site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
      y = as.vector(t(sim$y)),
      det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
    )
    od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                     det.covs = c("det_cov1", "pos_cov1"))
    cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
    y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

    fit <- tryCatch(
      suppressWarnings(tobs(
        formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
        family = occu_cover("lognormal"),
        detection = ~ det_cov1, positive = ~ pos_cov1,
        y = od$y, y_pos = y_pos, visits = od$det.covs,
        method = "nested_laplace",
        control = list(verbose = FALSE, max.iter = 80L,
                       engine = "joint_coupled")
      )),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    est_psi_x[s] <- fit$means["psi_occ_cov1"]
    est_p_x[s]   <- fit$means["p_det_cov1"]
    est_pos_x[s] <- fit$means["pos_pos_cov1"]
    se_psi_x[s]  <- fit$sds["psi_occ_cov1"]
    se_p_x[s]    <- fit$sds["p_det_cov1"]
    se_pos_x[s]  <- fit$sds["pos_pos_cov1"]
    est_alpha[s] <- fit$means[["alpha"]]
    est_sigma[s] <- fit$means[["sigma"]]
    field_cor[s] <- abs(cor(fit$spatial_field, sim$truth$f))
  }

  ok <- is.finite(est_psi_x) & is.finite(est_pos_x)
  expect_gte(sum(ok), 7L)

  expect_lt(abs(mean(est_psi_x[ok])     - beta_occ_truth[2L]), 0.35)
  expect_lt(abs(mean(est_p_x[ok])       - beta_p_truth[2L]),   0.35)
  expect_lt(abs(mean(est_pos_x[ok])     - beta_pos_truth[2L]), 0.20)

  expect_lt(abs(mean(est_alpha[ok]) - alpha_truth), 0.50)
  # The joint engine's coarse 5-point sigma_grid + Laplace's known small-N
  # underestimation of variance components together pull sigma's posterior
  # mean below truth. The wider tolerance here matches the v3 nested-Laplace
  # gate's reasoning (truth recovery within ~half a unit on a scale-like
  # axis) -- joint_coupled with a finer grid via control$sigma.grid
  # tightens this.
  expect_lt(abs(mean(est_sigma[ok]) - sigma_truth), 0.90)

  expect_gt(mean(field_cor[ok]), 0.80)

  # 95% Wald CI coverage on the shared-field path (working-family gate,
  # gcol33/tulpaObs#96): pooled over the three slope coefficients x seeds at the
  # 0.85 floor, with a loose 0.55 per-coordinate floor -- the field-coupled psi
  # slope carries the mild Laplace nested under-dispersion, so the pooled measure
  # (not any single coordinate) is the gate. Measured pooled coverage 0.94
  # (lognormal) at 18 seeds.
  cov_cells <- cbind(
    abs(est_psi_x[ok] - beta_occ_truth[2L]) < 1.96 * se_psi_x[ok],
    abs(est_p_x[ok]   - beta_p_truth[2L])   < 1.96 * se_p_x[ok],
    abs(est_pos_x[ok] - beta_pos_truth[2L]) < 1.96 * se_pos_x[ok])
  expect_gte(mean(cov_cells),          0.85)
  expect_gte(min(colMeans(cov_cells)), 0.55)
})


test_that("joint_coupled (beta arm) recovers slopes + field shape (10 seeds)", {
  skip_on_cran()
  skip_if_fast()

  n_seeds <- 18L
  N <- 100L; J <- 6L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }

  beta_occ_truth <- c(stats::qlogis(0.4),  0.7)
  beta_p_truth   <- c(0.0,                  0.8)
  beta_pos_truth <- c(stats::qlogis(0.3),  0.6)
  phi_pos_truth  <- 25
  sigma_truth    <- 1.0
  alpha_truth    <- 1.0

  est_psi_x <- est_p_x <- est_pos_x <- est_alpha <- est_sigma <-
    field_cor <- se_psi_x <- se_p_x <- se_pos_x <- rep(NA_real_, n_seeds)

  for (s in seq_len(n_seeds)) {
    sim <- simulate_occu_cover(
      N = N, J = J, beta_occ = beta_occ_truth, beta_p = beta_p_truth,
      beta_pos = beta_pos_truth, phi = phi_pos_truth,
      positive = "beta", adj = adj,
      sigma = sigma_truth, alpha = alpha_truth, seed = 7000L + s
    )
    long <- data.frame(
      site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
      y = as.vector(t(sim$y)),
      det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
    )
    od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                     det.covs = c("det_cov1", "pos_cov1"))
    cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
    y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

    fit <- tryCatch(
      suppressWarnings(tobs(
        formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
        family = occu_cover("beta"),
        detection = ~ det_cov1, positive = ~ pos_cov1,
        y = od$y, y_pos = y_pos, visits = od$det.covs,
        method = "nested_laplace",
        control = list(verbose = FALSE, max.iter = 80L,
                       engine = "joint_coupled")
      )),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    est_psi_x[s] <- fit$means["psi_occ_cov1"]
    est_p_x[s]   <- fit$means["p_det_cov1"]
    est_pos_x[s] <- fit$means["pos_pos_cov1"]
    se_psi_x[s]  <- fit$sds["psi_occ_cov1"]
    se_p_x[s]    <- fit$sds["p_det_cov1"]
    se_pos_x[s]  <- fit$sds["pos_pos_cov1"]
    est_alpha[s] <- fit$means[["alpha"]]
    est_sigma[s] <- fit$means[["sigma"]]
    field_cor[s] <- abs(cor(fit$spatial_field, sim$truth$f))
  }

  ok <- is.finite(est_psi_x) & is.finite(est_pos_x)
  expect_gte(sum(ok), 7L)

  expect_lt(abs(mean(est_psi_x[ok])     - beta_occ_truth[2L]), 0.35)
  expect_lt(abs(mean(est_p_x[ok])       - beta_p_truth[2L]),   0.35)
  expect_lt(abs(mean(est_pos_x[ok])     - beta_pos_truth[2L]), 0.25)

  expect_lt(abs(mean(est_alpha[ok]) - alpha_truth), 0.50)
  # Same coarse-grid + Laplace small-N underestimation as the lognormal
  # gate; pass control$sigma.grid to tighten.
  expect_lt(abs(mean(est_sigma[ok]) - sigma_truth), 0.90)

  expect_gt(mean(field_cor[ok]), 0.80)

  # 95% Wald CI coverage on the shared-field path (beta arm; working-family gate,
  # gcol33/tulpaObs#96). Measured pooled coverage 0.93 at 18 seeds.
  cov_cells <- cbind(
    abs(est_psi_x[ok] - beta_occ_truth[2L]) < 1.96 * se_psi_x[ok],
    abs(est_p_x[ok]   - beta_p_truth[2L])   < 1.96 * se_p_x[ok],
    abs(est_pos_x[ok] - beta_pos_truth[2L]) < 1.96 * se_pos_x[ok])
  expect_gte(mean(cov_cells),          0.85)
  expect_gte(min(colMeans(cov_cells)), 0.55)
})
